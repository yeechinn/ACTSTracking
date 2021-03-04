#ifndef ACTSDetectorLoaderProc_h
#define ACTSDetectorLoaderProc_h 1

#include <marlin/Processor.h>

//! \brief Load tracker geometry for ACTS
/**
 * @author Karol Krizka
 * @version $Id$
 */
class ACTSDetectorLoaderProc : public marlin::Processor
{
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

};

#endif
