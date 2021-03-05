#pragma once

#include <boost/container/flat_set.hpp>

namespace ACTSTracking
{

//! extract the geometry identifier from a variety of types
struct GeometryIdGetter
{
  //! explicit geometry identifier are just forwarded
  constexpr Acts::GeometryIdentifier operator()(
      Acts::GeometryIdentifier geometryId) const
  { return geometryId; }

  //! encoded geometry ids are converted back to geometry identifiers.
  constexpr Acts::GeometryIdentifier operator()(
      Acts::GeometryIdentifier::Value encoded) const
  { return Acts::GeometryIdentifier(encoded); }

  //! support elements in map-like structures.
  template <typename T>
  constexpr Acts::GeometryIdentifier operator()(
      const std::pair<Acts::GeometryIdentifier, T>& mapItem) const
  { return mapItem.first; }

  //! support elements that implement `.geometryId()`.
  template <typename T>
  inline auto operator()(const T& thing) const
      -> decltype(thing.geometryId(), Acts::GeometryIdentifier())
  { return thing.geometryId(); }
};

//! Compare two geometry ID's
struct CompareGeometryId
{
  // indicate that comparisons between keys and full objects are allowed.
  using is_transparent = void;
  // compare two elements using the automatic key extraction.
  template <typename Left, typename Right>
  constexpr bool operator()(Left&& lhs, Right&& rhs) const
  { return GeometryIdGetter()(lhs) < GeometryIdGetter()(rhs);}
};

//! Store elements that know their detector geometry id, e.g. simulation hits.
/**
 * @tparam T type to be stored, must be compatible with `CompareGeometryId`
 *
 * The container stores an arbitrary number of elements for any geometry
 * id. Elements can be retrieved via the geometry id; elements can be selected
 * for a specific geometry id or for a larger range, e.g. a volume or a layer
 * within the geometry hierachy using the helper functions below. Elements can
 * also be accessed by index that uniquely identifies each element regardless
 * of geometry id.
 */
template <typename T>
using GeometryIdMultiset =
    boost::container::flat_multiset<T, CompareGeometryId>;

}
