#include "LUXEDuplicateRemoval.hxx"

#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <algorithm>

namespace ACTSTracking {
inline int tracks_nshared(const EVENT::Track* trk1, const EVENT::Track* trk2) {
  const EVENT::TrackerHitVec& hits1 = trk1->getTrackerHits();
  const EVENT::TrackerHitVec& hits2 = trk2->getTrackerHits();

  // Number of overlapping hist
  uint32_t hitOlap = 0;
  for (const EVENT::TrackerHit* hit1 : hits1) {
    if (std::find(hits2.begin(), hits2.end(), hit1) != hits2.end()) {
      hitOlap++;
    }
  }

  // Smaller track count
  //uint32_t size = std::min(hits1.size(), hits2.size());

  return hitOlap;
}

/**
 * Return true if `trk1` is of higher quality than `trk2`.
 */
bool track_qua_compare(const EVENT::Track* trk1, const EVENT::Track* trk2) {
  // If number of hits are different, then the one with more
  // hits should be chosen first
  if (trk1->getTrackerHits().size() != trk2->getTrackerHits().size())
    return trk1->getTrackerHits().size() > trk2->getTrackerHits().size();

  // Same number of hits means I want smaller chi2
  return trk1->getChi2() < trk2->getChi2();
}

inline bool track_duplicate_compare(const EVENT::Track* trk1,
                                    const EVENT::Track* trk2) {
  return trk1->getTrackState(TrackState::AtIP)->getTanLambda() <
         trk2->getTrackState(TrackState::AtIP)->getTanLambda();
}
}  // namespace ACTSTracking

LUXEDuplicateRemoval aLUXEDuplicateRemoval;

LUXEDuplicateRemoval::LUXEDuplicateRemoval()
    : Processor("LUXEDuplicateRemoval") {
  // Input collections - tracks and relations
  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName",
                          "Name of track input collection",
                          _inputTrackCollection, std::string("TruthTracks"));

  // Output collections - tracks and relations
  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName",
                           "Name of track output collection",
                           _outputTrackCollection,
                           std::string("DedupedTruthTracks"));
}

void LUXEDuplicateRemoval::init() {
  // Print the initial parameters
  printParameters();
}

void LUXEDuplicateRemoval::processRunHeader(LCRunHeader* run) {}

void LUXEDuplicateRemoval::processEvent(LCEvent* evt) {
  LCCollection* incol = evt->getCollection(_inputTrackCollection);

  // Make the output track collection
  LCCollectionVec* outcol = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  outcol->setFlag(trkFlag.getFlag());

  //
  // Sort tracks by tan(lambda)/eta
  EVENT::TrackVec sortedInput;
  for (int i = 0; i < incol->getNumberOfElements(); ++i) {
    EVENT::Track* myTrk = static_cast<EVENT::Track*>(incol->getElementAt(i));

    EVENT::TrackVec::iterator insertion_point =
        std::upper_bound(sortedInput.begin(), sortedInput.end(), myTrk,
                         ACTSTracking::track_duplicate_compare);

    sortedInput.insert(insertion_point, myTrk);
  }

  //
  // Perform overlap removal by checking existing tracks and
  // adding it only if a matching track (50% shared hits)
  // is not found.
  std::vector<EVENT::Track*> finalTracks;
  std::vector<bool> selected(sortedInput.size(), true);
  for (int nsharedhits=8; nsharedhits>1; nsharedhits--){
    for (int itrk1=0; itrk1<sortedInput.size(); itrk1++){
      if (selected[itrk1]==false) continue;
      for (int itrk2=itrk1+1; itrk2<sortedInput.size(); itrk2++){
        if (selected[itrk2]==false) continue;
        if (ACTSTracking::tracks_nshared(sortedInput[itrk1], sortedInput[itrk2]) == nsharedhits){
          if (ACTSTracking::track_qua_compare(sortedInput[itrk1], sortedInput[itrk2])) selected[itrk2]=false; 
          else selected[itrk1]=false;
std::cout << " removing track " << std::endl;
        }
      }
    } 
  }

  for  (int itrk1=0; itrk1<sortedInput.size(); itrk1++){
    if (selected[itrk1]==false) continue;
    finalTracks.push_back(sortedInput[itrk1]);
  } 

  //
  // Create final collection
  for (EVENT::Track* trk : finalTracks) {
    outcol->addElement(trk);
  }

  // Save the output track collection
  outcol->setTransient(false);
  outcol->setSubset(true);

  evt->addCollection(outcol, _outputTrackCollection);
}

void LUXEDuplicateRemoval::check(LCEvent* evt) {}

void LUXEDuplicateRemoval::end() {}
