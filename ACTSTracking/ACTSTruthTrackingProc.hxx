#ifndef ACTSTruthTrackingProc_h
#define ACTSTruthTrackingProc_h 1

#include <marlin/Processor.h>

/*
#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
#include <EVENT/LCCollection.h>
#include "MarlinTrk/IMarlinTrkSystem.h"
*/
#include <EVENT/TrackerHit.h>
/*
#include <UTIL/CellIDDecoder.h>
#include "UTIL/LCTrackerConf.h"
#include <AIDA/AIDA.h>
*/

/**
 * This code performs a true pattern recognition by looping over all MC particles and adding all hits
 * associated to them onto a prototrack. This is then fitted and output.
 */
class ACTSTruthTrackingProc : public marlin::Processor
{		
 public:

  virtual marlin::Processor*  newProcessor() { return new ACTSTruthTrackingProc ; }

  ACTSTruthTrackingProc(const ACTSTruthTrackingProc &) = delete ;
  ACTSTruthTrackingProc& operator =(const ACTSTruthTrackingProc &) = delete ;
  ACTSTruthTrackingProc() ;

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
  /** Call to get collections
   */
  LCCollection* getCollection(const std::string&, LCEvent*);	

 protected:
	
  // Encoder
  /*UTIL::BitField64* m_encoder;*/

  // Get the subdetector ID from a hit
  int getSubdetector(const lcio::TrackerHit* hit){ /*m_encoder->setValue(hit->getCellID0()); return (*m_encoder)[lcio::LCTrackerCellID::subdet()];*/ } 

  // Get the layer ID from a hit
  int getLayer(const lcio::TrackerHit* hit){ /*m_encoder->setValue(hit->getCellID0()); return (*m_encoder)[lcio::LCTrackerCellID::layer()];*/ }

  // Remove hits in the same layer of the same subdetector
  void removeHitsSameLayer(const std::vector<lcio::TrackerHit*> &, std::vector<lcio::TrackerHit*> &);


  // Collection names for (in/out)put
  std::vector<std::string> _inputTrackerHitCollections ;
  std::vector<std::string> _inputTrackerHitRelationCollections ;
  std::string _inputParticleCollection ;
  std::string _outputTrackCollection ;
  std::string _outputTrackRelationCollection;

  // Run and event counters
  int m_eventNumber ;
  int m_runNumber ;
  /*
  // Track fit factory
  MarlinTrk::IMarlinTrkSystem* trackFactory;

  // Track fit parameters
  double m_initialTrackError_d0;
  double m_initialTrackError_phi0;
  double m_initialTrackError_omega;
  double m_initialTrackError_z0;
  double m_initialTrackError_tanL;
  double m_maxChi2perHit;
  */
  double m_magneticField;
  int m_fitFails;		
} ;

#endif



