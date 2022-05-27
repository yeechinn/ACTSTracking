#pragma once

#include <EVENT/TrackerHit.h>

#include "GeometryContainers.hxx"

namespace ACTSTracking {
//! \brief Link between an ACTS surface and hit index
class SourceLink final {
 public:
  //! \brief Construct from geometry identifier and hit
  SourceLink(Acts::GeometryIdentifier gid, std::size_t index,
             EVENT::TrackerHit* lciohit)
      : m_geometryId(gid), m_index(index), m_lciohit(lciohit) {}

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  SourceLink() = default;
  SourceLink(const SourceLink&) = default;
  SourceLink(SourceLink&&) = default;
  SourceLink& operator=(const SourceLink&) = default;
  SourceLink& operator=(SourceLink&&) = default;

  /// Access the geometry identifier.
  constexpr Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  /// Access the index.
  constexpr std::size_t index() const { return m_index; }
  /// Access the LCIO hit
  constexpr EVENT::TrackerHit* lciohit() const { return m_lciohit; }

 private:
  Acts::GeometryIdentifier m_geometryId;
  std::size_t m_index = -1;
  EVENT::TrackerHit* m_lciohit = nullptr;

  friend constexpr bool operator==(const SourceLink& lhs,
                                   const SourceLink& rhs) {
    return (lhs.m_geometryId == rhs.m_geometryId) and
           (lhs.m_index == rhs.m_index) and (lhs.m_lciohit == rhs.m_lciohit);
  }
  friend constexpr bool operator!=(const SourceLink& lhs,
                                   const SourceLink& rhs) {
    return not(lhs == rhs);
  }
};

/// Container of index source links
using SourceLinkContainer = GeometryIdMultiset<SourceLink>;
// Wrapper for SourceLinkContainer for use with CKF
using SourceLinkAccessor = GeometryIdMultisetAccessor<SourceLink>;
}  // namespace ACTSTracking
