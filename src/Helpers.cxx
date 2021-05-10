#include "Helpers.hxx"

#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

namespace ACTSTracking
{

EVENT::Track* ACTS2Marlin_track(const Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>& fitOutput,
                                std::size_t trackTip,
                                std::shared_ptr<Acts::MagneticFieldProvider> magneticField)
{
  IMPL::TrackImpl* track = new IMPL::TrackImpl ;
  Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates, trackTip);

  //
  // Fit state
  track->setChi2(trajState.chi2Sum);
  track->setNdf (trajState.NDF    );

  //
  // Track states

  // Track state: at IP
  static const Acts::Vector3 zeropos(0,0,0);
  const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.at(trackTip);
  EVENT::TrackState* trackStateAtIP
    = ACTSTracking::ACTS2Marlin_trackState(
					   EVENT::TrackState::AtIP,
					   params,
					   magneticField->getField(zeropos)[2]/Acts::UnitConstants::T
					   );
  track->trackStates().push_back(trackStateAtIP);

  // Track state: at first hit
  // placeholder, found when iterating over hits
  EVENT::TrackState* trackStateAtFirstHit=nullptr;

  // Track state: at last hit
  // placeholder, found when iterating over hits
  EVENT::TrackState* trackStateAtLastHit=nullptr;

  //
  // Hits on track
  EVENT::TrackerHitVec hitsOnTrack;
  fitOutput.fittedStates.visitBackwards(trackTip, [&](const Acts::MultiTrajectory<ACTSTracking::SourceLink>::ConstTrackStateProxy& state)
  {
    // No measurement at this state
    if(!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
    { return true; }

    // register all particles that generated this hit
    EVENT::TrackerHit *myHit=state.uncalibrated().lciohit();
    hitsOnTrack.push_back(myHit);

    // Save state at last hit
    // visiting backwards -> first hit seen
    if(trackStateAtLastHit==nullptr)
      {
	const Acts::Vector3 hitpos(myHit->getPosition()[0],myHit->getPosition()[1],myHit->getPosition()[2]);
	trackStateAtLastHit
	  = ACTSTracking::ACTS2Marlin_trackState(
						 EVENT::TrackState::AtLastHit,
						 state.smoothed(), state.smoothedCovariance(),
						 magneticField->getField(hitpos)[2]/Acts::UnitConstants::T
						 );
      }

    const Acts::Vector3 hitpos(myHit->getPosition()[0],myHit->getPosition()[1],myHit->getPosition()[2]);
    trackStateAtFirstHit
      = ACTSTracking::ACTS2Marlin_trackState(
					     EVENT::TrackState::AtFirstHit,
					     state.smoothed(), state.smoothedCovariance(),
					     magneticField->getField(hitpos)[2]/Acts::UnitConstants::T
					     );

    return true;
  });

  // Reverse hits, above creates them backwards
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  for(EVENT::TrackerHit* hit : hitsOnTrack)
    { track->addHit(hit); }

  // Save the track states at hits
  track->trackStates().push_back(trackStateAtFirstHit);
  track->trackStates().push_back(trackStateAtLastHit);

  return track;
}

EVENT::Track* ACTS2Marlin_track(const Acts::KalmanFitterResult<ACTSTracking::SourceLink>& fitOutput,
                                std::shared_ptr<Acts::MagneticFieldProvider> magneticField)
{
  IMPL::TrackImpl* track = new IMPL::TrackImpl ;
  Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates, fitOutput.trackTip);

  //
  // Fit state
  track->setChi2(trajState.chi2Sum);
  track->setNdf (trajState.NDF    );

  //
  // Track states

  // Track state: at IP
  static const Acts::Vector3 zeropos(0,0,0);
  const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.value();
  EVENT::TrackState* trackStateAtIP
    = ACTSTracking::ACTS2Marlin_trackState(
					   EVENT::TrackState::AtIP,
					   params,
					   magneticField->getField(zeropos)[2]/Acts::UnitConstants::T
					   );
  track->trackStates().push_back(trackStateAtIP);

  // Track state: at first hit
  // placeholder, found when iterating over hits
  EVENT::TrackState* trackStateAtFirstHit=nullptr;

  // Track state: at last hit
  // placeholder, found when iterating over hits
  EVENT::TrackState* trackStateAtLastHit=nullptr;

  //
  // Hits on track
  EVENT::TrackerHitVec hitsOnTrack;
  fitOutput.fittedStates.visitBackwards(fitOutput.trackTip, [&](const Acts::MultiTrajectory<ACTSTracking::SourceLink>::ConstTrackStateProxy& state)
  {
    // No measurement at this state
    if(!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
    { return true; }

    // register all particles that generated this hit
    EVENT::TrackerHit *myHit=state.uncalibrated().lciohit();
    hitsOnTrack.push_back(myHit);

    // Save state at last hit
    // visiting backwards -> first hit seen
    if(trackStateAtLastHit==nullptr)
      {
	const Acts::Vector3 hitpos(myHit->getPosition()[0],myHit->getPosition()[1],myHit->getPosition()[2]);
	trackStateAtLastHit
	  = ACTSTracking::ACTS2Marlin_trackState(
						 EVENT::TrackState::AtLastHit,
						 state.smoothed(), state.smoothedCovariance(),
						 magneticField->getField(hitpos)[2]/Acts::UnitConstants::T
						 );
      }

    const Acts::Vector3 hitpos(myHit->getPosition()[0],myHit->getPosition()[1],myHit->getPosition()[2]);
    trackStateAtFirstHit
      = ACTSTracking::ACTS2Marlin_trackState(
					     EVENT::TrackState::AtFirstHit,
					     state.smoothed(), state.smoothedCovariance(),
					     magneticField->getField(hitpos)[2]/Acts::UnitConstants::T
					     );

    return true;
  });

  // Reverse hits, above creates them backwards
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  for(EVENT::TrackerHit* hit : hitsOnTrack)
    { track->addHit(hit); }

  // Save the track states at hits
  track->trackStates().push_back(trackStateAtFirstHit);
  track->trackStates().push_back(trackStateAtLastHit);

  return track;
}

EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                          const Acts::BoundTrackParameters& params,
                                          double Bz)
{
  return ACTS2Marlin_trackState(location, params.parameters(), params.covariance().value(), Bz);
}

EVENT::TrackState* ACTS2Marlin_trackState(int location,
					  const Acts::BoundVector& value, const Acts::BoundMatrix& cov,
                                          double Bz)
{
  // Create new object
  IMPL::TrackStateImpl* trackState = new IMPL::TrackStateImpl();

  // Basic properties
  trackState->setLocation(location);

  //
  // Trajectory parameters

  // Central values
  double d0    =value[Acts::eBoundLoc0  ];
  double z0    =value[Acts::eBoundLoc1  ];
  double phi   =value[Acts::eBoundPhi   ];
  double theta =value[Acts::eBoundTheta ];
  double qoverp=value[Acts::eBoundQOverP];

  double p=1e3/qoverp;
  double omega=(0.3*Bz)/(p*std::sin(theta));
  double lambda=M_PI/2-theta;
  double tanlambda=std::tan(lambda);

  trackState->setPhi      (phi);
  trackState->setTanLambda(tanlambda);
  trackState->setOmega    (omega);
  trackState->setD0       (d0);
  trackState->setZ0       (z0);

  // Uncertainties (covariance matrix)
  Acts::ActsMatrix<6, 6> jac = Acts::ActsMatrix<6, 6>::Zero();

  jac(0, Acts::eBoundLoc0  ) =1;

  jac(1, Acts::eBoundPhi   ) =1;

  jac(2, Acts::eBoundTheta ) =omega/std::tan(theta);
  jac(2, Acts::eBoundQOverP) =omega/qoverp;

  jac(3, Acts::eBoundLoc1  ) =1;

  jac(4, Acts::eBoundTheta ) =std::pow(1/std::cos(theta), 2);
  
  Acts::ActsMatrix<6,6> trcov=(jac * cov * jac.transpose());

  EVENT::FloatVec lcioCov(15, 0);
  lcioCov[ 0]=trcov(0,0);
  lcioCov[ 1]=trcov(0,1);
  lcioCov[ 2]=trcov(1,1);
  lcioCov[ 3]=trcov(0,2);
  lcioCov[ 4]=trcov(1,2);
  lcioCov[ 5]=trcov(2,2);
  lcioCov[ 6]=trcov(0,3);
  lcioCov[ 7]=trcov(1,3);
  lcioCov[ 8]=trcov(2,3);
  lcioCov[ 9]=trcov(3,3);
  lcioCov[10]=trcov(0,4);
  lcioCov[11]=trcov(1,4);
  lcioCov[12]=trcov(2,4);
  lcioCov[13]=trcov(3,4);
  lcioCov[14]=trcov(4,4);

  trackState->setCovMatrix(lcioCov);

  return trackState;
}

EVENT::LCCollection* getCollection(EVENT::LCEvent* evt, const std::string& name)
{
  if( name.size() == 0 )
    return nullptr;

  try
  {
    return evt->getCollection( name );
  }
  catch(const EVENT::DataNotAvailableException& e)
  {
    // TODO: Reenable output
    //streamlog_out( DEBUG2 ) << "getCollection :  DataNotAvailableException : " << name <<  std::endl ;
    return nullptr;
  }
}
}
