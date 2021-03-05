#pragma once

#include "Defines.hxx"

namespace ACTSTracking
{
//! \brief Link between an ACTS surface and hit index
class SourceLink final
{
 public:
  //! \brief Construct from geometry identifier and hit
  SourceLink(Acts::GeometryIdentifier gid, std::size_t index)
      : m_geometryId(gid), m_index(index) {}

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  SourceLink() = default;
  SourceLink(const SourceLink&) = default;
  SourceLink(SourceLink&&) = default;
  SourceLink& operator=(const SourceLink&) = default;
  SourceLink& operator=(SourceLink&&) = default;

  /// Access the geometry identifier.
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  /// Access the index.
  std::size_t index() const { return m_index; }

 private:
  Acts::GeometryIdentifier m_geometryId;
  std::size_t m_index;

  friend constexpr bool operator==(const SourceLink& lhs,
                                   const SourceLink& rhs) {
    return (lhs.m_geometryId == rhs.m_geometryId) and
           (lhs.m_index == rhs.m_index);
  }
  friend constexpr bool operator!=(const SourceLink& lhs,
                                   const SourceLink& rhs) {
    return not(lhs == rhs);
  }
};

/// Container of index source links.
///
/// Since the source links provide a `.geometryId()` accessor, they can be
/// stored in an ordered geometry container.
using SourceLinkContainer = GeometryIdMultiset<SourceLink>;
}
