#include "ACTSTruthTrackingProc.hxx"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/LCRelationImpl.h>

#include <Acts/EventData/MultiTrajectory.hpp>

#include <Acts/Propagator/EigenStepper.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>

using namespace Acts::UnitLiterals;

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
  _description = "Build and fit tracks out of all hits associated to an MC particle" ;

  // Input collections - mc particles, tracker hits and the relationships between them
  registerInputCollections( LCIO::TRACKERHITPLANE,
                            "TrackerHitCollectionNames" ,
                            "Name of the TrackerHit input collections",
                            _inputTrackerHitCollections ,
                            {} ) ;

  registerInputCollections( LCIO::LCRELATION,
                            "SimTrackerHitRelCollectionNames",
                            "Name of TrackerHit SimTrackHit relation collections",
                            _inputTrackerHitRelationCollections,
                            {} );

  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollectionName",
                           "Name of the MCParticle input collection",
                           _inputParticleCollection,
                           std::string("MCParticle"));

  // Output collections - tracks and relations
  registerOutputCollection( LCIO::TRACK,
                            "TrackCollectionName",
                            "Name of track output collection",
                            _outputTrackCollection,
                            std::string("TruthTracks"));

  registerOutputCollection( LCIO::LCRELATION,
                            "TrackRelationCollectionName",
                            "Name of track to particle relation output collection",
                            _outputTrackRelationCollection,
                            std::string("TruthTrackRelations"));
}

void ACTSTruthTrackingProc::init()
{
  ACTSProcBase::init();
	
  // Reset counters
  _runNumber = 0 ;
  _eventNumber = 0 ;
  _fitFails = 0;

  /*
  // Set up the track fit factory
  trackFactory =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , nullptr , "" ) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  false) ;
  trackFactory->init() ;

  // Put default values for track fitting
  _initialTrackError_d0 = 1.e6;
  _initialTrackError_phi0 = 1.e2;
  _initialTrackError_omega = 1.e-4;
  _initialTrackError_z0 = 1.e6;
  _initialTrackError_tanL = 1.e2;
  _maxChi2perHit = 1.e2;
  */

  //Initialize CellID encoder
  _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());
}


void ACTSTruthTrackingProc::processRunHeader( LCRunHeader* )
{
  _runNumber++ ;
}

void ACTSTruthTrackingProc::processEvent( LCEvent* evt )
{
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

  // construct all components for the fitter
  Stepper stepper(magneticField());
  Navigator navigator(trackingGeometry());
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
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
    std::vector<const Acts::Surface*> surfSequence;
    ACTSTracking::MeasurementContainer track;
    for(const EVENT::TrackerHit* hit : trackFilteredByRHits)
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
      cov(0, 0) = std::pow(10_um, 2); // hit->getCovMatrix()[0]; // xx
      cov(0, 1) = hit->getCovMatrix()[1]; // yx
      cov(1, 1) = std::pow(10_um, 2); // hit->getCovMatrix()[2]; // yy
      cov(1, 0) = hit->getCovMatrix()[1]; // yx

      ACTSTracking::SourceLink sourceLink(surface->geometryId(), track.size());
      ACTSTracking::Measurement meas =
          Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0,
                                Acts::eBoundLoc1);

      track.push_back(meas);
      trackSourceLinks.push_back(sourceLink);
      surfSequence.push_back(surface);
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
        const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.value();
        streamlog_out(DEBUG) << "Fitted Paramemeters" << std::endl << params << std::endl;

        // Make the track object and relations object
        IMPL::LCRelationImpl* relationTrack = new IMPL::LCRelationImpl;
        IMPL::TrackImpl* track = new IMPL::TrackImpl ;

        Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
            Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates, fitOutput.trackTip);
        track->setChi2(trajState.chi2Sum);
        track->setNdf (trajState.NDF    );

        // TODO: Add hits on track
        
        //
        // AtIP: Overall fit results as fittedParameters

        IMPL::TrackStateImpl* trackStateAtIP = new IMPL::TrackStateImpl();
        trackStateAtIP->setLocation(lcio::TrackState::AtIP);
        track->trackStates().push_back(trackStateAtIP);

        // Fill the parameters
        static const Acts::Vector3 zeropos(0,0,0);

        double d0    =params.parameters()[Acts::eBoundLoc0  ];
        double z0    =params.parameters()[Acts::eBoundLoc1  ];
        double phi   =params.parameters()[Acts::eBoundPhi   ];
        double theta =params.parameters()[Acts::eBoundTheta ];
        double qoverp=params.parameters()[Acts::eBoundQOverP];
        
        double p=1e3/qoverp;
        double Bz=magneticField().getField(zeropos)[2]/Acts::UnitConstants::T;
        double omega=(0.3*Bz)/(p*std::sin(theta));
        double tanlambda=std::tan(theta);

        trackStateAtIP->setPhi      (phi);
        trackStateAtIP->setTanLambda(tanlambda);
        trackStateAtIP->setOmega    (omega);
        trackStateAtIP->setD0       (d0);
        trackStateAtIP->setZ0       (z0);

        // Fill the covariance matrix
        //d0, phi, omega, z0, tan(lambda)
        Acts::BoundTrackParameters::CovarianceMatrix cov=params.covariance().value();

        double var_d0    =cov(Acts::eBoundLoc0  , Acts::eBoundLoc0  );
        double var_z0    =cov(Acts::eBoundLoc1  , Acts::eBoundLoc1  );
        double var_phi   =cov(Acts::eBoundPhi   , Acts::eBoundPhi   );
        double var_theta =cov(Acts::eBoundTheta , Acts::eBoundTheta );
        double var_qoverp=cov(Acts::eBoundQOverP, Acts::eBoundQOverP);

        double var_omega    =
            var_qoverp*std::pow(omega/(qoverp*1e-3)      , 2) +
            var_theta *std::pow(omega/std::tan(var_theta), 2);
        double var_tanlambda=var_theta*std::pow(1/std::cos(theta), 4);
        
        FloatVec lcioCov(15, 0);
        lcioCov[ 0]=var_d0;
        lcioCov[ 2]=var_phi;
        lcioCov[ 5]=var_omega;
        lcioCov[ 9]=var_z0;
        lcioCov[14]=var_tanlambda;
        // TODO: Add off-diagonals

        trackStateAtIP->setCovMatrix(lcioCov);

        //
        // Other track states
        /** Can be used get track states at different layers
        fitOutput.fittedStates.visitBackwards(fitOutput.trackTip, [](Acts::MultiTrajectory<ACTSTracking::SourceLink>::ConstTrackStateProxy state)
        {
          const Acts::TrackStateType& typeFlags = state.typeFlags();
          if(!typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
            return true;

          const Acts::Surface& surface = state.referenceSurface();
          
          const Acts::GeometryIdentifier& geoID = surface.geometryId();
          std::cout << "volume = " << geoID.volume() << std::endl;
          std::cout << "layer = " << geoID.layer() << std::endl;
          std::cout << "sensitive = " << geoID.sensitive() << std::endl;

          const Acts::BoundVector& params=state.smoothed();
          std::cout << params[Acts::eBoundQOverP] << std::endl;
          return true;
        });
        */

        //
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
  // Save the output particle to track relation collection
  evt->addCollection( trackRelationCollection , _outputTrackRelationCollection ) ;

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
