#ifndef ACTSProcBase_h
#define ACTSProcBase_h 1

#include <marlin/Processor.h>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>

#include <Acts/Utilities/CalibrationContext.hpp>

#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>

#include "ACTSGeometryIdMappingTool.hxx"

//! \brief Base processor for ACTS tracking
/**
 * Performs tasks common to all ACTS processors
 *  - loading tracking geometry
 *
 * Assumes that the global TGeoManager with the geometry
 * description is already loaded. For example, via the
 * InitializeDD4hep processor.
 *
 * @author Karol Krizka
 * @version $Id$
 */
class ACTSProcBase : public marlin::Processor
{
  using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

 public:

  ACTSProcBase(const ACTSProcBase &) = delete ;
  ACTSProcBase& operator =(const ACTSProcBase &) = delete ;
  ACTSProcBase(const std::string& procname) ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;

 private:
  void buildDetector();

  void buildBfield();

 protected:
  std::string _matFile {};

  std::shared_ptr<ACTSGeometryIdMappingTool> geoIDMappingTool() const;
  
  const Acts::MagneticFieldContext& magneticFieldContext() const;
  const Acts::GeometryContext& geometryContext() const;
  const Acts::CalibrationContext& calibrationContext() const;

  const Acts::ConstantBField& magneticField() const;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

  //! Find surface for hit
  const Acts::Surface* findSurface(const EVENT::TrackerHit* hit) const;

 private:
  std::shared_ptr<ACTSGeometryIdMappingTool> _geoIDMappingTool;
  
  Acts::MagneticFieldContext _magneticFieldContext;
  Acts::ConstantBField _magneticField;

  Acts::GeometryContext _geometryContext;
  DetectorStore _detectorStore;
  std::shared_ptr<const Acts::TrackingGeometry> _trackingGeometry =nullptr;

  Acts::CalibrationContext _calibrationContext;
};

#endif
