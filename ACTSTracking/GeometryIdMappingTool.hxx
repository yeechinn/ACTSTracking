#pragma once

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

#include <unordered_map>

namespace ACTSTracking {

//! \brief Maps DD4hep cell ID's to ACTS geometry ID's
/**
 * Lots of hardcoded values that match up to mod5 geometry
 * dumped to TGeo.
 *
 * @author Karol Krizka
 * @version $Id$
 */
class GeometryIdMappingTool {
 public:
  /** Create a mapping tool using the provided encoderString to
   * interpret cell ID's.
   */
  GeometryIdMappingTool(const std::string& encoderString);

  /** Decode hit
   */
  uint64_t getGeometryID(const lcio::SimTrackerHit* hit);
  uint64_t getGeometryID(const lcio::TrackerHit* hit);

  uint64_t getGeometryID(uint32_t systemID, uint32_t layerID, int32_t sideID,
                         uint32_t ladderID, uint32_t moduleID);

 private:
  std::string _encoderString;

};

}  // namespace ACTSTracking
