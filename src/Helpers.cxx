#include "Helpers.hxx"

#include <IMPL/TrackStateImpl.h>

namespace ACTSTracking
{

EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                          const Acts::BoundTrackParameters& params,
                                          double Bz)
{
  // Create new object
  IMPL::TrackStateImpl* trackState = new IMPL::TrackStateImpl();

  // Basic properties
  trackState->setLocation(location);

  //
  // Trajectory parameters

  // Central values
  const Acts::BoundTrackParameters::ParametersVector& value = params.parameters();
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
  const Acts::BoundTrackParameters::CovarianceMatrix& cov=params.covariance().value();

  double var_d0    =cov(Acts::eBoundLoc0  , Acts::eBoundLoc0  );
  double var_z0    =cov(Acts::eBoundLoc1  , Acts::eBoundLoc1  );
  double var_phi   =cov(Acts::eBoundPhi   , Acts::eBoundPhi   );
  double var_theta =cov(Acts::eBoundTheta , Acts::eBoundTheta );
  double var_qoverp=cov(Acts::eBoundQOverP, Acts::eBoundQOverP);

  double var_omega    =
      var_qoverp*std::pow(omega/(qoverp*1e-3)      , 2) +
      var_theta *std::pow(omega/std::tan(var_theta), 2);
  double var_tanlambda=var_theta*std::pow(1/std::cos(theta), 4);

  EVENT::FloatVec lcioCov(15, 0);
  lcioCov[ 0]=var_d0;
  lcioCov[ 2]=var_phi;
  lcioCov[ 5]=var_omega;
  lcioCov[ 9]=var_z0;
  lcioCov[14]=var_tanlambda;

  trackState->setCovMatrix(lcioCov);

  return trackState;
}
}
