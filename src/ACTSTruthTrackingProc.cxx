#include "ACTSTruthTrackingProc.hxx"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <IMPL/TrackImpl.h>
#include <IMPL/LCRelationImpl.h>

#include <Acts/EventData/MultiTrajectory.hpp>

#include <Acts/Propagator/EigenStepper.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>

using namespace Acts::UnitLiterals;

#include "Helpers.hxx"
#include "MeasurementCalibrator.hxx"
#include "SourceLink.hxx"

using TrackFitterResult =
    Acts::Result<Acts::KalmanFitterResult<ACTSTracking::SourceLink>>;

// sorting by value of R(=x^2+y^2) in global coordinated so the hits are always
// sorted from close to the IP outward
bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2)
{
  double radius1 =
      sqrt((hit1->getPosition()[0]) * (hit1->getPosition()[0]) + (hit1->getPosition()[1]) * (hit1->getPosition()[1]));
  double radius2 =
      sqrt((hit2->getPosition()[0]) * (hit2->getPosition()[0]) + (hit2->getPosition()[1]) * (hit2->getPosition()[1]));
  return radius1 < radius2;
}

ACTSTruthTrackingProc aACTSTruthTrackingProc;

ACTSTruthTrackingProc::ACTSTruthTrackingProc() : ACTSProcBase("ACTSTruthTrackingProc")
{
  // modify processor description
  _description = "Build and fit tracks out of all hits associated to an MC particle." ;

  // Input collections - mc particles, tracker hits and the relationships between them
  registerInputCollections( LCIO::TRACKERHITPLANE,
                            "TrackerHitCollectionNames" ,
                            "Name of the TrackerHit input collections.",
                            _inputTrackerHitCollections ,
                            {} ) ;

  registerInputCollections( LCIO::LCRELATION,
                            "SimTrackerHitRelCollectionNames",
                            "Name of TrackerHit to SimTrackerHit relation collections.",
                            _inputTrackerHitRelationCollections,
                            {} );

  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollectionName",
                           "Name of the MCParticle input collection.",
                           _inputParticleCollection,
                           std::string("MCParticle"));

  // Output collections - tracks and relations
  registerOutputCollection( LCIO::TRACK,
                            "TrackCollectionName",
                            "Name of track output collection.",
                            _outputTrackCollection,
                            std::string("TruthTracks"));
}

void ACTSTruthTrackingProc::init()
{
  ACTSProcBase::init();
	
  // Reset counters
  _runNumber = 0 ;
  _eventNumber = 0 ;
  _fitFails = 0;

  //Initialize CellID encoder
  _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());
}


void ACTSTruthTrackingProc::processRunHeader( LCRunHeader* )
{
  _runNumber++ ;
}

void ACTSTruthTrackingProc::processEvent( LCEvent* evt )
{
  //
  // Caches
  Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();
  Acts::MagneticFieldProvider::Cache magCache = magneticField()->makeCache(magFieldContext);

  //
  // Initialize track finder
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Stepper = Acts::EigenStepper<>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

  // Configurations
  Navigator::Config navigatorCfg{trackingGeometry()};
  navigatorCfg.resolvePassive   = false;
  navigatorCfg.resolveMaterial  = true;
  navigatorCfg.resolveSensitive = true;

  // construct all components for the fitter
  Stepper stepper(magneticField());
  Navigator navigator(navigatorCfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));

  // Get the collection of MC particles
  LCCollection* particleCollection = getCollection(_inputParticleCollection, evt);
  if(particleCollection == nullptr) return;

  // Make objects to hold all of the tracker hit, simulated hit and relation collections
  std::vector<LCCollection*> trackerHitCollections;
  std::vector<LCCollection*> trackerHitRelationCollections;
  std::vector<std::shared_ptr<LCRelationNavigator>> relations;
  
  // Loop over each input collection and get the data
  for(unsigned int collection=0; collection<_inputTrackerHitCollections.size(); collection++)
  {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = getCollection(_inputTrackerHitCollections[collection], evt);
    if(trackerHitCollection == nullptr) continue;
    trackerHitCollections.push_back(trackerHitCollection);

    // Get the collection of tracker hit relations
    LCCollection* trackerHitRelationCollection = getCollection(_inputTrackerHitRelationCollections[collection], evt);
    trackerHitRelationCollections.push_back(trackerHitRelationCollection);

    // Create the relations navigator
    std::shared_ptr<LCRelationNavigator> relation = std::make_shared<LCRelationNavigator>( trackerHitRelationCollection );
    relations.push_back(relation);
  }

  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec( LCIO::TRACK )  ;
	
  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackCollection->setFlag( trkFlag.getFlag()  ) ;
  
  // Make the output particle to track relation collection
  LCCollectionVec* trackRelationCollection = new LCCollectionVec( LCIO::LCRELATION )  ;

  /* 
     Now for each MC particle we want the list of hits belonging to it. The most 
     efficient way is to loop over all hits once, and store the pointers in a 
     map, with the key a pointer to the MC particle. We can then loop over each
     MC particle at the end and get all of the hits, before making a track.
  */
  // Make the container
  std::map<MCParticle*, std::vector<TrackerHit*> > particleHits;

  // Loop over all input collections
  for(unsigned int collection=0; collection<trackerHitCollections.size(); collection++)
  {
    // Loop over tracker hits
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for(int itHit=0;itHit<nHits;itHit++)
    {
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
      const double* globalpos = hit->getPosition();

      // Get the related simulated hit(s)
      const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );

      // Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
      SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));

      // If the hit was produced by a secondary which was not saved to the MCParticle collection
      if(simHit->isProducedBySecondary())
        continue;

      // Get the particle belonging to that hit
      MCParticle* particle = simHit->getMCParticle();

      // Push back the element into the container
      particleHits[particle].push_back(hit);
    }
  }

  // Now loop over all particles and get the list of hits
  int nParticles = particleCollection->getNumberOfElements();
  for(int itP=0;itP<nParticles;itP++)
  {
    // Get the particle
    MCParticle* mcParticle = static_cast<MCParticle*>( particleCollection->getElementAt(itP) ) ;

    // Get the vector of hits from the container
    if(particleHits.count(mcParticle) == 0) continue;
    std::vector<TrackerHit*> trackHits = particleHits[mcParticle];

    // Only make tracks with 3 or more hits
    if(trackHits.size() < 3) continue;

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Remove the hits on the same layers (removing those with higher R)
    EVENT::TrackerHitVec trackFilteredByRHits;
    removeHitsSameLayer(trackHits, trackFilteredByRHits);
    if(trackFilteredByRHits.size() < 3) continue;

    // Make container
    //MeasurementContainer track;
    std::vector<ACTSTracking::SourceLink> trackSourceLinks;
    ACTSTracking::MeasurementContainer track;
    for(EVENT::TrackerHit* hit : trackFilteredByRHits)
    {
      // Convert to Acts hit
      const Acts::Surface* surface=findSurface(hit);

      const double* globalpos=hit->getPosition();
      Acts::Result<Acts::Vector2> lpResult = surface->globalToLocal(geometryContext(),
                                                                    {globalpos[0], globalpos[1], globalpos[2]},
                                                                    {0,0,0},
                                                                    0.5_um);
      if(!lpResult.ok())
        throw std::runtime_error("Global to local transformation did not succeed.");

      Acts::Vector2 loc = lpResult.value();

      Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
      const EVENT::TrackerHitPlane* hitplane=dynamic_cast<const EVENT::TrackerHitPlane*>(hit);
      if(hitplane)
      {
        cov(0, 0) = std::pow(hitplane->getdU()*Acts::UnitConstants::mm, 2);
        cov(1, 1) = std::pow(hitplane->getdV()*Acts::UnitConstants::mm, 2);
      }
      else
      { throw std::runtime_error("Currently only support TrackerHitPlane."); }

      ACTSTracking::SourceLink sourceLink(surface->geometryId(), track.size(), hit);
      ACTSTracking::Measurement meas =
          Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0,
                                Acts::eBoundLoc1);

      track.push_back(meas);
      trackSourceLinks.push_back(sourceLink);
    }

    //
    // Setup tracker
    // Construct a perigee surface as the target surface
    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3{0., 0., 0.});

    // Set the KalmanFitter options
    //std::unique_ptr<const Acts::Logger> logger=Acts::getDefaultLogger("TrackFitting", Acts::Logging::Level::VERBOSE);
    Acts::KalmanFitterOptions<ACTSTracking::MeasurementCalibrator, Acts::VoidOutlierFinder> kfOptions
        =Acts::KalmanFitterOptions<ACTSTracking::MeasurementCalibrator, Acts::VoidOutlierFinder>(
            geometryContext(), magneticFieldContext(), calibrationContext(),
            ACTSTracking::MeasurementCalibrator(track), Acts::VoidOutlierFinder(),
            //Acts::LoggerWrapper{*logger}, Acts::PropagatorPlainOptions(),
            Acts::getDummyLogger(), Acts::PropagatorPlainOptions(),
            &(*perigeeSurface));

    double px=mcParticle->getMomentum()[0];
    double py=mcParticle->getMomentum()[1];
    double pz=mcParticle->getMomentum()[2];
    double pt=sqrt(px*px+py*py);
    double p =sqrt(px*px+py*py+pz*pz);

    Acts::BoundVector params = Acts::BoundVector::Zero();
    // position/time
    params[Acts::eBoundLoc0  ] = 0;
    params[Acts::eBoundLoc1  ] = 0;
    params[Acts::eBoundTime  ] = mcParticle->getTime();
    // direction angles phi,theta
    params[Acts::eBoundPhi   ] = atan2(py,px);
    params[Acts::eBoundTheta ] = atan2(pt,pz);
    // absolute momentum vector
    params[Acts::eBoundQOverP] = mcParticle->getCharge()/p;

    // build the track covariance matrix using the smearing sigmas 
    Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
    cov(Acts::eBoundLoc0  , Acts::eBoundLoc0  ) = std::pow(_initialTrackError_d0              ,2);
    cov(Acts::eBoundLoc1  , Acts::eBoundLoc1  ) = std::pow(_initialTrackError_z0              ,2);
    cov(Acts::eBoundTime  , Acts::eBoundTime  ) = std::pow(1_ns                               ,2);
    cov(Acts::eBoundPhi   , Acts::eBoundPhi   ) = std::pow(_initialTrackError_phi             ,2);
    cov(Acts::eBoundTheta , Acts::eBoundTheta ) = std::pow(_initialTrackError_lambda          ,2);
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = std::pow(_initialTrackError_relP * p /(p*p), 2);

    std::shared_ptr<Acts::PerigeeSurface> particleSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3(mcParticle->getVertex()));

    Acts::BoundTrackParameters initialparams(perigeeSurface, params, mcParticle->getCharge(), cov);
    streamlog_out(DEBUG) << "Initial Paramemeters" << std::endl << initialparams << std::endl;

    TrackFitterResult result=trackFitter.fit(trackSourceLinks, initialparams, kfOptions);

    if (result.ok())
    {
      const Acts::KalmanFitterResult<ACTSTracking::SourceLink>& fitOutput = result.value();
      if (fitOutput.fittedParameters)
      {
        // Make the track object and relations object
        IMPL::LCRelationImpl* relationTrack = new IMPL::LCRelationImpl;

        //
        // Helpful debug output
        Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
            Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates, fitOutput.lastMeasurementIndex);
        streamlog_out(DEBUG) << "Trajectory Summary" << std::endl;
        streamlog_out(DEBUG) << "\tchi2Sum       " << trajState.chi2Sum       << std::endl;
        streamlog_out(DEBUG) << "\tNDF           " << trajState.NDF           << std::endl;
        streamlog_out(DEBUG) << "\tnHoles        " << trajState.nHoles        << std::endl;
        streamlog_out(DEBUG) << "\tnMeasurements " << trajState.nMeasurements << std::endl;
        streamlog_out(DEBUG) << "\tnOutliers     " << trajState.nOutliers     << std::endl;
        streamlog_out(DEBUG) << "\tnStates       " << trajState.nStates       << std::endl;

        const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.value();
        streamlog_out(DEBUG) << "Fitted Paramemeters" << std::endl << params << std::endl;

        // Make track object
        EVENT::Track* track = ACTSTracking::ACTS2Marlin_track(fitOutput, magneticField(), magCache);

        // Save results
        trackCollection->addElement(track);

        // Make the particle to track link
        relationTrack->setFrom(track);
        relationTrack->setTo(mcParticle);
        relationTrack->setWeight(1.0);
        trackRelationCollection->addElement(relationTrack);
      }
      else
      {
        streamlog_out(WARNING) << "No fitted paramemeters for track" << std::endl;
        _fitFails++;
      }
    }
    else
    {
      streamlog_out(WARNING) << "Track fit error: " << result.error() << std::endl;
      _fitFails++;
    }
  }

  // Save the output track collection
  evt->addCollection( trackCollection , _outputTrackCollection ) ;

  // Increment the event number
  _eventNumber++ ;
}

void ACTSTruthTrackingProc::check( LCEvent* )
{
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ACTSTruthTrackingProc::end()
{
  streamlog_out(MESSAGE) << " end()  " << name()
                         << " processed " << _eventNumber << " events in " << _runNumber << " runs "
                         << std::endl ;
}

LCCollection* ACTSTruthTrackingProc::getCollection(const std::string& collectionName, LCEvent* evt)
{
  try
  {
    return evt->getCollection( collectionName );
  }
  catch(DataNotAvailableException &e)
  {
    streamlog_out( DEBUG5 ) << "- cannot get collection. Collection " << collectionName << " is unavailable" << std::endl;
    return nullptr;
  }
}

int ACTSTruthTrackingProc::getSubdetector(const lcio::TrackerHit* hit)
{ _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::subdet()]; }

int ACTSTruthTrackingProc::getLayer(const lcio::TrackerHit* hit)
{ _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::layer ()]; }

void ACTSTruthTrackingProc::removeHitsSameLayer(const std::vector<TrackerHit*> &trackHits, std::vector<TrackerHit*> &trackFilteredHits)
{
  trackFilteredHits.push_back(*(trackHits.begin()));

  for(std::vector<TrackerHit*>::const_iterator it = trackHits.begin()+1; it != trackHits.end(); ++it)
  {
    int subdet = getSubdetector(*it);
    int layer = getLayer(*it);
    if( subdet != getSubdetector(*(it-1)) )
    {
      trackFilteredHits.push_back(*it);
    }
    else if( layer != getLayer(*(it-1)) )
    {
      trackFilteredHits.push_back(*it);
    }
  }
}
