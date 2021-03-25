#pragma once

#include <EVENT/TrackState.h>

#include <Acts/EventData/TrackParameters.hpp>

namespace ACTSTracking
{

//! Convert ACTS track state class to Marlin class
/**
 * \param location Location where the track state is defined (ie: `AtIP`)
 * \param params ACTS track state parameters
 * \params Bz magnetic field at location of track state [Tesla]
 *
 * \return Track state with equivalent parameters of the ACTS track
 */
EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                         const Acts::BoundTrackParameters& params,
                                         double Bz);
}
