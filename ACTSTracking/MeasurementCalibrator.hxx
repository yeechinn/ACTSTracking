#pragma once

#include <Acts/EventData/Measurement.hpp>

#include "SourceLink.hxx"

namespace ACTSTracking {
//! Hit stored as an measurement
using Measurement = Acts::BoundVariantMeasurement<ACTSTracking::SourceLink>;

//! Collection of measurements
using MeasurementContainer = std::vector<Measurement>;

class MeasurementCalibrator {
 public:
  /// Construct an invalid calibrator. Required to allow copying.
  MeasurementCalibrator() = default;
  /// Construct using a user-provided container to chose measurements from.
  MeasurementCalibrator(const MeasurementContainer& measurements)
      : m_measurements(&measurements) {}

  //! Find the measurement corresponding to the source link.
  /**
   * @tparam parameters_t Track parameters type
   * @param sourceLink Input source link
   * @param parameters Input track parameters (unused)
   */
  template <typename parameters_t>
  const Measurement& operator()(const SourceLink& sourceLink,
                                const parameters_t& /* parameters */) const {
    assert(m_measurements and
           "Undefined measurement container in DigitizedCalibrator");
    assert((sourceLink.index() < m_measurements->size()) and
           "Source link index is outside the container bounds");
    return (*m_measurements)[sourceLink.index()];
  }

 private:
  // use pointer so the calibrator is copyable and default constructible.
  const MeasurementContainer* m_measurements = nullptr;
};

}  // namespace ACTSTracking
