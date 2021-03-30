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

  //
  // Sort tracks into groups of overlapping tracks
  std::vector<EVENT::TrackVec> trkGroups;
  for(int i=0; i < incol->getNumberOfElements() ; ++i)
  {
    EVENT::Track* myTrk = static_cast<EVENT::Track*>( incol->getElementAt(i) ) ;
    const TrackerHitVec& myHits=myTrk->getTrackerHits();

    bool inGroup=false;
    for(EVENT::TrackVec& trkGroup : trkGroups)
    {
      for(EVENT::Track* otherTrk : trkGroup)
      {
        const TrackerHitVec& otherHits=otherTrk->getTrackerHits();

        uint32_t hitOlap=0;
        for(const EVENT::TrackerHit* myHit : myHits)
        {
          if(std::find(otherHits.begin(), otherHits.end(), myHit)!=otherHits.end())
          { hitOlap++; }
        }

        // 50% matching hits means the tracks overlap
        if(2*hitOlap>myHits.size() || 2*hitOlap>otherHits.size())
        {
          inGroup=true;

          // Insert to matched group, better track candidate first
          EVENT::TrackVec::iterator pos
              = std::lower_bound(trkGroup.begin(), trkGroup.end(), myTrk, ACTSTracking::track_duplicate_compare);
          trkGroup.insert(pos, myTrk);

          break;
        }
      }
    }

    // Not added to any existing group, create a new one
    if(!inGroup)
    { trkGroups.push_back({myTrk}); }
  }

  //
  // Perform overlap removal by the selecting the track with
  // the most hits and best

  // Make the output track collection
  LCCollectionVec* outcol = new LCCollectionVec( LCIO::TRACK )  ;

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  outcol->setFlag( trkFlag.getFlag()  ) ;

  // The overlap check..
  for(const EVENT::TrackVec& trkGroup : trkGroups)
  { outcol->addElement(trkGroup[0]); }

  // Save the output track collection
  outcol->setTransient( false ) ;
  outcol->setSubset( true ) ;

  evt->addCollection( outcol , _outputTrackCollection ) ;
}
  
void ACTSDuplicateRemoval::check( LCEvent * evt )
{ }

void ACTSDuplicateRemoval::end()
{ }
