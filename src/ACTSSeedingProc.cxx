#include "ACTSSeedingProc.hxx"

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

#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>

#include <Acts/Seeding/BinFinder.hpp>
#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Seeding/Seedfinder.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>

#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

using namespace Acts::UnitLiterals;

#include "GeometryIdSelector.hxx"
#include "Helpers.hxx"
#include "MeasurementCalibrator.hxx"
#include "SeedSpacePoint.hxx"
#include "SourceLink.hxx"

// Track fitting definitions
using TrackFinderOptions =
    Acts::CombinatorialKalmanFilterOptions<ACTSTracking::MeasurementCalibrator,
                                           Acts::MeasurementSelector>;

using TrackFinderResult =
    Acts::Result<Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>>;

using TrackFinderResultContainer =
    std::vector<TrackFinderResult>;

ACTSSeedingProc aACTSSeedingProc;

ACTSSeedingProc::ACTSSeedingProc() : ACTSProcBase("ACTSSeedingProc")
{
  // modify processor description
  _description = "Build and fit tracks out of all hits associated to an MC particle" ;

  // Settings
  registerProcessorParameter("InitialTrackError_RelP",
                             "Track error estimate, momentum component (relative)",
                             _initialTrackError_relP,
                             0.25);

  registerProcessorParameter("InitialTrackError_Phi",
                             "Track error estimate, phi (radians)",
                             _initialTrackError_phi,
                             1_degree);

  registerProcessorParameter("InitialTrackError_Lambda",
                             "Track error estimate, lambda (radians)",
                             _initialTrackError_lambda,
                             1_degree);

  registerProcessorParameter("InitialTrackError_Pos",
                             "Track error estimate, local position (mm)",
                             _initialTrackError_pos,
                             10_um);

  // Input collections - mc particles, tracker hits and the relationships between them
  registerInputCollections( LCIO::TRACKERHITPLANE,
                            "TrackerHitCollectionNames" ,
                            "Name of the TrackerHit input collections",
                            _inputTrackerHitCollections ,
                            {} ) ;

  // Output collections - tracks and relations
  registerOutputCollection( LCIO::TRACK,
                            "SeedCollectionName",
                            "Name of seed output collection",
                            _outputSeedCollection,
                            std::string("SeedTracks"));

  registerOutputCollection( LCIO::TRACK,
                            "TrackCollectionName",
                            "Name of track output collection",
                            _outputTrackCollection,
                            std::string("SeededCKFTracks"));
}

void ACTSSeedingProc::init()
{
  ACTSProcBase::init();
	
  // Reset counters
  _runNumber = 0 ;
  _eventNumber = 0 ;
  _fitFails = 0;
}


void ACTSSeedingProc::processRunHeader( LCRunHeader* )
{
  _runNumber++ ;
}

void ACTSSeedingProc::processEvent( LCEvent* evt )
{
  //
  // Select layers to use for seeding.
  // Selecting a volume without an explicit layer selects all
  // layers within the volume.
  ACTSTracking::GeometryIdSelector geometrySelection ({
      // vertex negative endcap
      Acts::GeometryIdentifier().setVolume(13).setLayer( 4),
      Acts::GeometryIdentifier().setVolume(13).setLayer( 8),
      Acts::GeometryIdentifier().setVolume(13).setLayer(12),
      Acts::GeometryIdentifier().setVolume(13).setLayer(16),
      // vertex barrel
      Acts::GeometryIdentifier().setVolume(14).setLayer( 2),
      Acts::GeometryIdentifier().setVolume(14).setLayer( 6),
      Acts::GeometryIdentifier().setVolume(14).setLayer(10),
      Acts::GeometryIdentifier().setVolume(14).setLayer(14),
      // vertex positive endcap
      Acts::GeometryIdentifier().setVolume(15).setLayer( 2),
      Acts::GeometryIdentifier().setVolume(15).setLayer( 6),
      Acts::GeometryIdentifier().setVolume(15).setLayer(10),
      Acts::GeometryIdentifier().setVolume(15).setLayer(14)
    });

  //
  // Prepare the output
  // Make the output track collection
  LCCollectionVec* seedCollection  = new LCCollectionVec( LCIO::TRACK )  ;
  LCCollectionVec* trackCollection = new LCCollectionVec( LCIO::TRACK )  ;
  
  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackCollection->setFlag( trkFlag.getFlag()  ) ;

  //
  // Prepare input hits in ACTS format
  ACTSTracking::SourceLinkContainer sourceLinks;
  ACTSTracking::MeasurementContainer measurements;
  ACTSTracking::SeedSpacePointContainer spacePoints;

  for(unsigned int collection=0; collection<_inputTrackerHitCollections.size(); collection++)
  {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = getCollection(_inputTrackerHitCollections[collection], evt);
    if(trackerHitCollection == nullptr) continue;

    for(int itHit=0;itHit<trackerHitCollection->getNumberOfElements();itHit++)
    {
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollection->getElementAt(itHit) ) ;

      // Convert to Acts hit
      const Acts::Surface* surface=findSurface(hit);

      const double* lcioglobalpos = hit->getPosition();
      Acts::Vector3 globalPos={lcioglobalpos[0], lcioglobalpos[1], lcioglobalpos[2]};
      Acts::Result<Acts::Vector2> lpResult = surface->globalToLocal(geometryContext(),
                                                                    globalPos,
                                                                    {0,0,0},
                                                                    0.5_um);
      if(!lpResult.ok())
        throw std::runtime_error("Global to local transformation did not succeed.");

      Acts::Vector2 loc = lpResult.value();

      Acts::SymMatrix2 localCov = Acts::SymMatrix2::Zero();
      const EVENT::TrackerHitPlane* hitplane=dynamic_cast<const EVENT::TrackerHitPlane*>(hit);
      if(hitplane)
      {
        localCov(0, 0) = std::pow(hitplane->getdU()*Acts::UnitConstants::mm, 2);
        localCov(1, 1) = std::pow(hitplane->getdV()*Acts::UnitConstants::mm, 2);
      }
      else
      { throw std::runtime_error("Currently only support TrackerHitPlane."); }

      ACTSTracking::SourceLink sourceLink(surface->geometryId(), measurements.size(), hit);
      ACTSTracking::Measurement meas =
          Acts::makeMeasurement(sourceLink, loc, localCov, Acts::eBoundLoc0,
                                Acts::eBoundLoc1);

      measurements.push_back(meas);
      sourceLinks .push_back(sourceLink);

      //
      // Seed selection and conversion to useful coordinates
      if(geometrySelection.check(surface->geometryId()))
      {
        Acts::RotationMatrix3 rotLocalToGlobal =
            surface->referenceFrame(geometryContext(),
                                    globalPos,
                                    {0,0,0});

        // Convert to a seed space point
        // the space point requires only the variance of the transverse and
        // longitudinal position. reduce computations by transforming the
        // covariance directly from local to rho/z.
        //
        // compute Jacobian from global coordinates to rho/z
        //
        //         rho = sqrt(x² + y²)
        // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
        //             = 2 * {x,y} / r
        //       dz/dz = 1 (duuh!)
        //
        double x = globalPos[Acts::ePos0];
        double y = globalPos[Acts::ePos1];
        double scale = 2 / std::hypot(x, y);
        Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
        jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
        jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
        jacXyzToRhoZ(1, Acts::ePos2) = 1;
        // compute Jacobian from local coordinates to rho/z
        Acts::ActsMatrix<2, 2> jac =
            jacXyzToRhoZ *
            rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
        // compute rho/z variance
        Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

        // Save spacepoint
        spacePoints.push_back(ACTSTracking::SeedSpacePoint(globalPos, var[0], var[1], sourceLink.index()));
      }
    }
  }

  streamlog_out( DEBUG0 )  << "Created " << spacePoints.size() << " space points" << std::endl;

  //
  // Run the seeding algorithm

  // Finder configuration
  static const Acts::Vector3 zeropos(0,0,0);

  Acts::SeedfinderConfig<ACTSTracking::SeedSpacePoint> finderCfg;
  finderCfg.rMax = 120;
  finderCfg.deltaRMin =  5;
  finderCfg.deltaRMax = 30; 
  finderCfg.collisionRegionMin = -250;
  finderCfg.collisionRegionMax = 250;
  finderCfg.zMin = -300;
  finderCfg.zMax =  300;
  finderCfg.maxSeedsPerSpM = 1;
  finderCfg.cotThetaMax = 7.40627;  // 2.7 eta;
  finderCfg.sigmaScattering = 50;
  finderCfg.radLengthPerSeed = 0.1;
  finderCfg.minPt = 500;
  finderCfg.bFieldInZ = magneticField()->getField(zeropos)[2]/Acts::UnitConstants::T*1e-3;
  finderCfg.beamPos = {0,0};
  finderCfg.impactMax = 3;

  Acts::SpacePointGridConfig gridCfg;
  gridCfg.bFieldInZ = finderCfg.bFieldInZ;
  gridCfg.cotThetaMax = finderCfg.cotThetaMax;
  gridCfg.deltaRMax = finderCfg.deltaRMax;
  gridCfg.minPt = finderCfg.minPt;
  gridCfg.rMax = finderCfg.rMax;
  gridCfg.zMax = finderCfg.zMax;
  gridCfg.zMin = finderCfg.zMin;

  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = finderCfg.maxSeedsPerSpM;

  finderCfg.seedFilter = std::make_unique<Acts::SeedFilter<ACTSTracking::SeedSpacePoint>>(
      Acts::SeedFilter<ACTSTracking::SeedSpacePoint>(filterCfg));

  // Create tools
  std::function<Acts::Vector2(const ACTSTracking::SeedSpacePoint&, float, float, float)> extractCovariance
      = [](const ACTSTracking::SeedSpacePoint& sp, float, float, float) -> Acts::Vector2
      { return {sp.varianceR(), sp.varianceZ()}; };


  std::vector<const ACTSTracking::SeedSpacePoint*> spacePointPtrs(spacePoints.size(), nullptr);
  std::transform(spacePoints.begin(), spacePoints.end(), spacePointPtrs.begin(),
                 [](const ACTSTracking::SeedSpacePoint& sp)
                 { return &sp;} );
  
  std::shared_ptr<Acts::BinFinder<ACTSTracking::SeedSpacePoint>> bottomBinFinder
      = std::make_shared<Acts::BinFinder<ACTSTracking::SeedSpacePoint>>();
  std::shared_ptr<Acts::BinFinder<ACTSTracking::SeedSpacePoint>> topBinFinder
      = std::make_shared<Acts::BinFinder<ACTSTracking::SeedSpacePoint>>();

  std::unique_ptr<Acts::SpacePointGrid<ACTSTracking::SeedSpacePoint>> grid
      = Acts::SpacePointGridCreator::createGrid<ACTSTracking::SeedSpacePoint>(gridCfg);
  
  Acts::BinnedSPGroup<ACTSTracking::SeedSpacePoint> spacePointsGrouping (
      spacePointPtrs.begin(), spacePointPtrs.end(), extractCovariance,
      bottomBinFinder, topBinFinder, std::move(grid), finderCfg);

  Acts::Seedfinder<ACTSTracking::SeedSpacePoint> finder(finderCfg);

  std::vector<Acts::Seed<ACTSTracking::SeedSpacePoint>> seeds;
  Acts::BinnedSPGroupIterator<ACTSTracking::SeedSpacePoint> group = spacePointsGrouping.begin();
  Acts::BinnedSPGroupIterator<ACTSTracking::SeedSpacePoint> groupEnd = spacePointsGrouping.end();
  for (; !(group == groupEnd); ++group)
  {
    std::vector<Acts::Seed<ACTSTracking::SeedSpacePoint>> myseeds=finder.createSeedsForGroup(group.bottom(), group.middle(), group.top());

    seeds.insert(seeds.end(), myseeds.begin(), myseeds.end());
  }


  //
  // Loop over seeds and get track parameters
  std::vector<Acts::BoundTrackParameters> paramseeds;
  for(const Acts::Seed<ACTSTracking::SeedSpacePoint>& seed : seeds)
  {
    // Get the bottom space point and its reference surface
    // @todo do we need to sort the sps first
    const ACTSTracking::SeedSpacePoint* bottomSP = seed.sp().front();
    const std::size_t hitIdx = bottomSP->measurementIndex();
    const ACTSTracking::SourceLink& sourceLink = sourceLinks.at(hitIdx);
    const Acts::GeometryIdentifier& geoId = sourceLink.geometryId();

    const Acts::Surface* surface = trackingGeometry()->findSurface(geoId);
    if (surface == nullptr) {
      std::cout << "surface with geoID "
                << geoId << " is not found in the tracking gemetry";
      continue;
    }

    // Get the magnetic field at the bottom space point
    Acts::Vector3 field = magneticField()->getField(
        {bottomSP->x(), bottomSP->y(), bottomSP->z()});

    std::optional<Acts::BoundVector> optParams = Acts::estimateTrackParamsFromSeed(
        geometryContext(), seed.sp().begin(), seed.sp().end(), *surface, field,
        0.1_T);
    if (!optParams.has_value())
    {
      std::cout << "Failed estimation of track parameters for seed." << std::endl;
      continue;
    }

    const Acts::BoundVector& params = optParams.value();

    float charge = std::copysign(1, params[Acts::eBoundQOverP]);
    float p = std::abs(1/params[Acts::eBoundQOverP]);

    // build the track covariance matrix using the smearing sigmas 
    Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
    cov(Acts::eBoundLoc0  , Acts::eBoundLoc0  ) = std::pow(_initialTrackError_pos             ,2);
    cov(Acts::eBoundLoc1  , Acts::eBoundLoc1  ) = std::pow(_initialTrackError_pos             ,2);
    cov(Acts::eBoundTime  , Acts::eBoundTime  ) = std::pow(_initialTrackError_time            ,2);
    cov(Acts::eBoundPhi   , Acts::eBoundPhi   ) = std::pow(_initialTrackError_phi             ,2);
    cov(Acts::eBoundTheta , Acts::eBoundTheta ) = std::pow(_initialTrackError_lambda          ,2);
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = std::pow(_initialTrackError_relP * p /(p*p) ,2);

    Acts::BoundTrackParameters paramseed(surface->getSharedPtr(), params, charge, cov);
    paramseeds.push_back(paramseed);

    //
    // Add seed to LCIO collection
    IMPL::TrackImpl* seedtrack = new IMPL::TrackImpl;
    seedCollection->addElement(seedtrack);

    Acts::Vector3 globalPos = surface->localToGlobal(geometryContext(),
                                                     {params[Acts::eBoundLoc0], params[Acts::eBoundLoc1]},
                                                     {0,0,0});

    EVENT::TrackState* seedTrackState
        = ACTSTracking::ACTS2Marlin_trackState(
            lcio::TrackState::AtFirstHit,
            paramseed,
            magneticField()->getField(globalPos)[2]/Acts::UnitConstants::T
                                               );
    seedtrack->trackStates().push_back(seedTrackState);

    streamlog_out(DEBUG) << "Seed Paramemeters" << std::endl << paramseed << std::endl;
  }

  //
  // Initialize track finder
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Stepper = Acts::EigenStepper<>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using CKF =
      Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother>;

  // construct all components for the fitter
  Stepper stepper(magneticField());
  Navigator navigator(trackingGeometry());
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  Propagator propagator(std::move(stepper), std::move(navigator));
  CKF trackFinder(std::move(propagator));

  // Set the options
  Acts::MeasurementSelector::Config measurementSelectorCfg={{Acts::GeometryIdentifier(), {15,10}}};

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  // Construct a perigee surface as the target surface
  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  //std::unique_ptr<const Acts::Logger> logger=Acts::getDefaultLogger("TrackFitting", Acts::Logging::Level::VERBOSE);

  TrackFinderOptions ckfOptions
      =TrackFinderOptions(
          geometryContext(), magneticFieldContext(), calibrationContext(),
          ACTSTracking::MeasurementCalibrator(std::move(measurements)),
          Acts::MeasurementSelector(measurementSelectorCfg),
          //Acts::LoggerWrapper{*logger}, pOptions,
          Acts::getDummyLogger(), pOptions,
          &(*perigeeSurface));

  //
  // Find the tracks
  TrackFinderResultContainer results=trackFinder.findTracks(sourceLinks, paramseeds, ckfOptions);
  for (TrackFinderResult& result : results)
  {
    if (result.ok())
    {
      const Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>& fitOutput = result.value();
      for(const size_t& trackTip : fitOutput.trackTips)
      {
        if(fitOutput.fittedParameters.count(trackTip)==0)
        {
          streamlog_out(WARNING) << "No fitted track parameters for trajectory with entry index = " << trackTip << std::endl;
          continue;
        }

        //
        // Helpful debug output
        Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
            Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates, trackTip);
        streamlog_out(DEBUG) << "Trajectory Summary" << std::endl;
        streamlog_out(DEBUG) << "\tchi2Sum       " << trajState.chi2Sum       << std::endl;
        streamlog_out(DEBUG) << "\tNDF           " << trajState.NDF           << std::endl;
        streamlog_out(DEBUG) << "\tnHoles        " << trajState.nHoles        << std::endl;
        streamlog_out(DEBUG) << "\tnMeasurements " << trajState.nMeasurements << std::endl;
        streamlog_out(DEBUG) << "\tnOutliers     " << trajState.nOutliers     << std::endl;
        streamlog_out(DEBUG) << "\tnStates       " << trajState.nStates       << std::endl;

        const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.at(trackTip);
        streamlog_out(DEBUG) << "Fitted Paramemeters" << std::endl << params << std::endl;

        // Make track object
        EVENT::Track* track = ACTSTracking::ACTS2Marlin_track(fitOutput, trackTip, magneticField());

        // Save results
        trackCollection->addElement(track);
      }
    }
    else
    {
      streamlog_out(WARNING) << "Track fit error: " << result.error() << std::endl;
      _fitFails++;
    }
  }

  // Save the output seed collection
  evt->addCollection( seedCollection  , _outputSeedCollection  ) ;

  // Save the output track collection
  evt->addCollection( trackCollection , _outputTrackCollection ) ;

  // Increment the event number
  _eventNumber++ ;
}

void ACTSSeedingProc::check( LCEvent* )
{
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ACTSSeedingProc::end()
{
  streamlog_out(MESSAGE) << " end()  " << name()
                         << " processed " << _eventNumber << " events in " << _runNumber << " runs "
                         << std::endl ;
}

LCCollection* ACTSSeedingProc::getCollection(const std::string& collectionName, LCEvent* evt)
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
