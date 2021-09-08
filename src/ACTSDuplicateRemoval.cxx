#include "ACTSDuplicateRemoval.hxx"

#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <algorithm>

namespace ACTSTracking
{
  /**
   * Return true if `trk1` and `trk2` share at least 50% of hits.
   */
  inline bool tracks_equal(const EVENT::Track* trk1, const EVENT::Track* trk2)
  {
    const EVENT::TrackerHitVec& hits1=trk1->getTrackerHits();
    const EVENT::TrackerHitVec& hits2=trk2->getTrackerHits();

    // Number of overlapping hist
    uint32_t hitOlap=0;
    for(const EVENT::TrackerHit* hit1 : hits1)
      {
	if(std::find(hits2.begin(), hits2.end(), hit1)!=hits2.end())
	  { hitOlap++; }
      }

    // Smaller track count
    uint32_t size = std::min(hits1.size(), hits2.size());

    return 2*hitOlap>size; // half of smaller track belong to larger track
  }

  /**
   * Return true if `trk1` is of higher quality than `trk2`.
   */
  bool track_quality_compare(const EVENT::Track* trk1, const EVENT::Track* trk2)
  {
    // If number of hits are different, then the one with more
    // hits should be chosen first
    if(trk1->getTrackerHits().size() != trk2->getTrackerHits().size())
      return trk1->getTrackerHits().size() > trk2->getTrackerHits().size();

    // Same number of hits means I want smaller chi2
    return trk1->getChi2()<trk2->getChi2();
  }

  inline bool track_duplicate_compare(const EVENT::Track* trk1, const EVENT::Track* trk2)
  {
    return trk1->getTrackState(TrackState::AtIP)->getTanLambda()
      < trk2->getTrackState(TrackState::AtIP)->getTanLambda();
  }
}

ACTSDuplicateRemoval aACTSDuplicateRemoval;

ACTSDuplicateRemoval::ACTSDuplicateRemoval()
  : Processor("ACTSDuplicateRemoval")
{
  // Input collections - tracks and relations
  registerInputCollection( LCIO::TRACK,
			   "InputTrackCollectionName",
			   "Name of track input collection",
			   _inputTrackCollection,
			   std::string("TruthTracks"));

  // Output collections - tracks and relations
  registerOutputCollection( LCIO::TRACK,
                            "OutputTrackCollectionName",
                            "Name of track output collection",
                            _outputTrackCollection,
                            std::string("DedupedTruthTracks"));
}


void ACTSDuplicateRemoval::init()
{  
  // Print the initial parameters
  printParameters() ;
}

void ACTSDuplicateRemoval::processRunHeader( LCRunHeader* run )
{ }

void ACTSDuplicateRemoval::processEvent( LCEvent * evt )
{
  LCCollection* incol = evt->getCollection(_inputTrackCollection);

  // Make the output track collection
  LCCollectionVec* outcol = new LCCollectionVec( LCIO::TRACK )  ;

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  outcol->setFlag( trkFlag.getFlag()  ) ;

  //
  // Sort tracks by quality (insertion sort)
  EVENT::TrackVec sortedInput;
  for(int i=0; i < incol->getNumberOfElements() ; ++i)
  {
    EVENT::Track* myTrk = static_cast<EVENT::Track*>( incol->getElementAt(i) ) ;

    EVENT::TrackVec::iterator insertion_point =
      std::upper_bound(sortedInput.begin(), sortedInput.end(), myTrk, ACTSTracking::track_duplicate_compare);

    sortedInput.insert(insertion_point, myTrk);
  }

  //
  // Perform overlap removal by checking existing tracks and
  // adding it only if a matching track (50% shared hits)
  // is not found.
  std::vector<EVENT::Track*> finalTracks;
  for(EVENT::Track* myTrk : sortedInput)
    {
      bool foundAnEqual=false;
      for(int i=(finalTracks.size() >= 10) ? finalTracks.size()-10 : 0;
	  i < finalTracks.size();
	  ++i)
	{
	  EVENT::Track* otherTrk = finalTracks[i];

	  // Skip tracks that are not equal
	  if(!ACTSTracking::tracks_equal(myTrk, otherTrk))
	    { continue; }
	  foundAnEqual=true;

	  // Replace if my track is better
	  if(ACTSTracking::track_quality_compare(myTrk, otherTrk))
	    {
	      finalTracks[i]=myTrk;
	      break;
	    }
	}

      if(!foundAnEqual) // Add a new track that does not duplicate
	{ finalTracks.push_back(myTrk);	}
    }

  //
  // Create final collection
  for(EVENT::Track* trk : finalTracks)
    { outcol->addElement(trk); }

  // Save the output track collection
  outcol->setTransient( false ) ;
  outcol->setSubset( true ) ;

  evt->addCollection( outcol , _outputTrackCollection ) ;
}
  
void ACTSDuplicateRemoval::check( LCEvent * evt )
{ }

void ACTSDuplicateRemoval::end()
{ }
