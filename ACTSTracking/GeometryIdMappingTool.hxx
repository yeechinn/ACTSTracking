#pragma once

#include <unordered_map>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

namespace ACTSTracking
{

//! \brief Maps DD4hep cell ID's to ACTS geometry ID's
/**
 * Lots of hardcoded values that match up to mod5 geometry
 * dumped to TGeo.
 * 
 * @author Karol Krizka
 * @version $Id$
 */
class GeometryIdMappingTool
{
public:
  /** Create a mapping tool using the provided encoderString to
   * interpret cell ID's.
   */
  GeometryIdMappingTool(const std::string& encoderString);

  /** Decode hit
   */
  uint64_t getGeometryID(const lcio::SimTrackerHit* hit);
  uint64_t getGeometryID(const lcio::TrackerHit* hit);

  uint64_t getGeometryID(uint32_t systemID, uint32_t layerID, int32_t sideID, uint32_t ladderID, uint32_t moduleID);
  
private:
  std::string _encoderString;

  static const std::unordered_map<int32_t, uint32_t> VolumeMap;

  static const int32_t VertexEndCapNegative;
  static const int32_t VertexBarrel;
  static const int32_t VertexEndCapPositive;
  static const int32_t InnerTrackerEndCapNegative;
  static const int32_t InnerTrackerBarrel;
  static const int32_t InnerTrackerEndCapPositive;
  static const int32_t OuterInnerTrackerEndCapNegative;
  static const int32_t OuterInnerTrackerBarrel;
  static const int32_t OuterInnerTrackerEndCapPositive;
  static const int32_t OuterTrackerEndCapNegative;
  static const int32_t OuterTrackerBarrel;
  static const int32_t OuterTrackerEndCapPositive;

  // Modules in phi ladder per layer
  static const std::unordered_map<uint32_t, uint32_t> NLad_VertexBarrel;
  static const std::unordered_map<uint32_t, uint32_t> NLad_InnerTrackerBarrel;
  static const std::unordered_map<uint32_t, uint32_t> NLad_OuterInnerTrackerBarrel;
  static const std::unordered_map<uint32_t, uint32_t> NLad_OuterTrackerBarrel;

  // Modules in ring per layer
  static const std::unordered_map<uint32_t, uint32_t> NRng_VertexEndCap;
  static const std::unordered_map<uint32_t, uint32_t> NRng_InnerTrackerEndCap;
  static const std::unordered_map<uint32_t, uint32_t> NRng_OuterInnerTrackerEndCap;
  static const std::unordered_map<uint32_t, uint32_t> NRng_OuterTrackerEndCap;

};

}
