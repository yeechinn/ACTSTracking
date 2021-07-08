#include "ACTSDuplicateRemoval.hxx"

#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <algorithm>

namespace ACTSTracking
{
bool track_duplicate_compare(const EVENT::Track* trk1, const EVENT::Track* trk2)
{
  // If number of hits are different, then the one with more
  // hits should be chosen first
  if(trk1->getTrackerHits().size() != trk2->getTrackerHits().size())
    return trk1->getTrackerHits().size() > trk2->getTrackerHits().size();

  // Same number of hits means I want smaller chi2
  return trk1->getChi2()<trk2->getChi2();
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

  for(EVENT::Track* myTrk : sortedInput)
    {
      bool addme=true;

      const EVENT::TrackerHitVec& myHits=myTrk->getTrackerHits();

      for(int i=0; i < outcol->getNumberOfElements() ; ++i)
	{
	  const EVENT::Track* otherTrk = static_cast<const EVENT::Track*>( outcol->getElementAt(i) ) ;
	  const EVENT::TrackerHitVec& otherHits=otherTrk->getTrackerHits();

	  // Count overlapping hits
	  uint32_t hitOlap=0;
	  for(const EVENT::TrackerHit* myHit : myHits)
	    {
	      if(std::find(otherHits.begin(), otherHits.end(), myHit)!=otherHits.end())
		{ hitOlap++; }
	    }

	  if(2*hitOlap>myHits.size()) // half my hits belong to superior track
	    {
	      addme=false;
	      break;
	    }
	}

      if(addme)
	{ outcol->addElement(myTrk); }
    }

  // Save the output track collection
  outcol->setTransient( false ) ;
  outcol->setSubset( true ) ;

  evt->addCollection( outcol , _outputTrackCollection ) ;
}
  
void ACTSDuplicateRemoval::check( LCEvent * evt )
{ }

void ACTSDuplicateRemoval::end()
{ }
