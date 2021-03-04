#include "ACTSDetectorLoaderProc.hxx"

ACTSDetectorLoaderProc aACTSDetectorLoaderProc ;


ACTSDetectorLoaderProc::ACTSDetectorLoaderProc()
  : Processor("ACTSDetectorLoaderProc")
{
  // modify processor description
  _description = "Load the tracking detector for ACTS." ;
}

void ACTSDetectorLoaderProc::init()
{ }

void ACTSDetectorLoaderProc::processRunHeader( LCRunHeader* run )
{ }

void ACTSDetectorLoaderProc::processEvent( LCEvent * evt )
{ }
  
void ACTSDetectorLoaderProc::check( LCEvent * evt )
{ }

void ACTSDetectorLoaderProc::end()
{ }
