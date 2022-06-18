#include "TrackTruthProc.hxx"

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>

#include <IMPL/LCCollectionVec.h>

#include <UTIL/LCRelationNavigator.h>

#include "Helpers.hxx"

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

//------------------------------------------------------------------------------------------------

TrackTruthProc aTrackTruthProc;

TrackTruthProc::TrackTruthProc() : Processor("TrackTruthProc") {
  // modify processor description
  _description = "Associate MCParticle to a reconstructed Track";

  registerInputCollection(LCIO::TRACK, "TrackCollection",
                          "Name of reconstructed track input collection",
                          _inColTrack, std::string("Tracks"));

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollectionName",
                          "Name of the MCParticle input collection", _inColMCP,
                          std::string("MCParticle"));

  registerInputCollections(
      LCIO::LCRELATION, "TrackerHit2SimTrackerHitRelationName",
      "Name of TrackerHit to SimTrackHit relation collection", _inColH2SH, {});

  registerOutputCollection(LCIO::LCRELATION, "Particle2TrackRelationName",
                           "Map from MC particle to reconstructed track.",
                           _outColMC2T,
                           std::string("Particle2TrackRelationName"));
}

//============================================================================================================================

void TrackTruthProc::init() {
  // usually a good idea to
  printParameters();
}
//============================================================================================================================

void TrackTruthProc::processRunHeader(LCRunHeader* /*run*/) {}

//============================================================================================================================

void TrackTruthProc::processEvent(LCEvent* evt) {
  //
  // Load relations
  std::vector<std::shared_ptr<LCRelationNavigator>> hit2simhits;
  for (const std::string& name : _inColH2SH) {
    // Get the collection of tracker hit relations
    LCCollection* trackerHitRelationCollection =
        ACTSTracking::getCollection(evt, name);
    std::shared_ptr<LCRelationNavigator> hit2simhit =
        std::make_shared<LCRelationNavigator>(trackerHitRelationCollection);
    hit2simhits.push_back(hit2simhit);
  }

  //
  // MC particles
  LCCollection* particleCollection =
      ACTSTracking::getCollection(evt, _inColMCP);
  if (particleCollection == nullptr) return;
  int nParticles = particleCollection->getNumberOfElements();

  // store best track for each mc particle
  std::map<EVENT::MCParticle*, EVENT::Track*> mcBestMatch_track;
  std::map<EVENT::MCParticle*, float> mcBestMatch_frac;

  //
  // Tracks
  LCCollection* trackCollection = ACTSTracking::getCollection(evt, _inColTrack);
  if (trackCollection == nullptr) return;
  int nTracks = trackCollection->getNumberOfElements();

  for (int itT = 0; itT < nTracks; ++itT) {
    // Get the track
    EVENT::Track* track =
        static_cast<EVENT::Track*>(trackCollection->getElementAt(itT));

    // Loop over all hits in a track and associate it to a MC particle
    std::map<EVENT::MCParticle*, uint32_t> trackhit2mc;
    for (EVENT::TrackerHit* hit : track->getTrackerHits()) {
      // Find the sim hit
      EVENT::SimTrackerHit* simHit = nullptr;
      for (std::shared_ptr<LCRelationNavigator> hit2simhit : hit2simhits) {
        const LCObjectVec& simHitVector = hit2simhit->getRelatedToObjects(hit);
        if (!simHitVector.empty()) {  // Found the sim hit
          simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
          break;
        }
      }

      // Increment MC particle counter
      if (simHit->getMCParticle() != nullptr)
        trackhit2mc[simHit->getMCParticle()]++;
    }

    // Update best matches
    for (const std::pair<EVENT::MCParticle*, uint32_t>& mchit : trackhit2mc) {
      float frac =
          static_cast<float>(mchit.second) / track->getTrackerHits().size();

      bool better =
          mcBestMatch_track.count(mchit.first) == 0 ||  // no best match exists
          mcBestMatch_frac[mchit.first] <
              frac;  // this match is better (more hits on track)

      if (better) {
        mcBestMatch_track[mchit.first] = track;
        mcBestMatch_frac[mchit.first] = frac;
      }
    }
  }

  //
  // Save the best matches
  LCRelationNavigator relMC2T(LCIO::MCPARTICLE, LCIO::TRACK);
  for (const std::pair<EVENT::MCParticle*, EVENT::Track*>& mctrk :
       mcBestMatch_track) {
    relMC2T.addRelation(mctrk.first, mctrk.second,
                        mcBestMatch_frac[mctrk.first]);
  }

  LCCollectionVec* outColMC2T = (LCCollectionVec*)relMC2T.createLCCollection();
  outColMC2T->setTransient(true);

  evt->addCollection(outColMC2T, _outColMC2T);
}

//============================================================================================================================

void TrackTruthProc::check(LCEvent* /*evt*/) {
  // nothing to check here - could be used to fill checkplots in reconstruction
  // processor
}

//============================================================================================================================

void TrackTruthProc::end() {}
