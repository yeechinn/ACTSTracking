#ifndef ACTSDetectorLoaderProc_h
#define ACTSDetectorLoaderProc_h 1

#include <marlin/Processor.h>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>

//! \brief Load tracker geometry for ACTS
/**
 * Assumes that the global TGeoManager with the geometry
 * description is already loaded. For example, via the
 * InitializeDD4hep processor.
 *
 * @author Karol Krizka
 * @version $Id$
 */
class ACTSDetectorLoaderProc : public marlin::Processor
{
  using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

 public:

  virtual marlin::Processor* newProcessor() { return new ACTSDetectorLoaderProc ; }

  ACTSDetectorLoaderProc(const ACTSDetectorLoaderProc &) = delete ;
  ACTSDetectorLoaderProc& operator =(const ACTSDetectorLoaderProc &) = delete ;
  ACTSDetectorLoaderProc() ;

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

 protected:
  std::string _matFile {};

 private:
  DetectorStore _detectorStore;
  Acts::GeometryContext _tGeoContext;
  std::shared_ptr<const Acts::TrackingGeometry> _trackingGeometry;
};

#endif
