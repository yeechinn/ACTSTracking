#ifndef TrackTruthProc_h
#define TrackTruthProc_h 1

#include <lcio.h>
#include <marlin/Processor.h>
#include <string>
#include <vector>

/**
 * Helper processor that creates LCRelation collections for track to hit
 * associations to be used with LCTuple.
 *
 * @param  TrackCollection                Names of Track input collections
 * @param  Track2HitRelationName          Name of output collection for track to
 * hit relations
 */

class TrackTruthProc : public marlin::Processor {
 public:
  virtual Processor* newProcessor() { return new TrackTruthProc; }

  TrackTruthProc();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

 protected:
  /** Input collection names.
   */
  std::string _inColTrack;
  std::string _inColMCP;
  std::vector<std::string> _inColH2SH;

  /** Output collection names.
   */
  std::string _outColMC2T;
};

#endif
