#include "GeometryIdMappingTool.hxx"

#include <iomanip>

using namespace ACTSTracking;

GeometryIdMappingTool::GeometryIdMappingTool(const std::string& encoderString)
    : _encoderString(encoderString) {}

uint64_t GeometryIdMappingTool::getGeometryID(const lcio::SimTrackerHit* hit) {
  UTIL::CellIDDecoder<lcio::SimTrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(hit)["system"];
  uint32_t layerID = decoder(hit)["layer"];
  int32_t sideID = decoder(hit)["side"];
  uint32_t ladderID = decoder(hit)["module"];
  uint32_t moduleID = decoder(hit)["sensor"];

  return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

uint64_t GeometryIdMappingTool::getGeometryID(const lcio::TrackerHit* hit) {
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(hit)["system"];
  uint32_t layerID = decoder(hit)["layer"];
  int32_t sideID = decoder(hit)["side"];
  uint32_t ladderID = decoder(hit)["module"];
  uint32_t moduleID = decoder(hit)["sensor"];

  return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

uint64_t GeometryIdMappingTool::getGeometryID(uint32_t systemID,
                                              uint32_t layerID, int32_t sideID,
                                              uint32_t ladderID,
                                              uint32_t moduleID) {
  uint64_t geometry_id = 0;

  uint64_t volume_id=1;
  geometry_id |= volume_id << (14 * 4);

  uint64_t layer_id;
  layer_id = layerID%2==0?layerID+2:layerID+1;
  geometry_id |= layer_id << (9 * 4);

  uint64_t sensitive_id;
  sensitive_id = layerID%2==0?ladderID+1:ladderID+10;
  geometry_id |= sensitive_id << (0 * 4);
  return geometry_id;
}
