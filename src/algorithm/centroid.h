// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CENTROID_H_
#define SFCGAL_ALGORITHM_CENTROID_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL {
namespace algorithm {

class WeightedCentroid {
public:
  SFCGAL::Kernel::FT             area;
  CGAL::Vector_3<SFCGAL::Kernel> centroid;
  SFCGAL::Kernel::FT             m;

  WeightedCentroid(
      SFCGAL::Kernel::FT             a    = 0.0,
      CGAL::Vector_3<SFCGAL::Kernel> c    = CGAL::Vector_3<SFCGAL::Kernel>(),
      SFCGAL::Kernel::FT             mTmp = 0.0)
      : area(a), centroid(c), m(mTmp)
  {
  }
};

/**
 * @brief Compute the centroid for a Geometry
 * @ingroup public_api
 * @warning Z component is ignored, there is no projection for 3D geometries
 * @ingroup public_api
 * @pre g is a valid geometry
 */
SFCGAL_API std::unique_ptr<Point>
           centroid(const Geometry &g);

/**
 * @brief Compute the centroid for a Geometry
 * @ingroup public_api
 * @warning Z component is ignored, there is no projection for 3D geometries
 * @ingroup public_api
 * @pre g is a valid geometry
 */
SFCGAL_API std::unique_ptr<Point>
           centroid3D(const Geometry &g);

/**
 * Returns Compute the centroid for a Geometry
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Geometry &g, bool enable3DComputation = false);

/**
 * Returns Compute the centroid for a Triangle
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Triangle &g, bool enable3DComputation = false);

/**
 * Returns Compute the centroid for a Triangle
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Point &a, const Point &b, const Point &c,
                 bool enable3DComputation = false);

/**
 * Returns Compute the centroid for a LineString
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const LineString &g, bool enable3DComputation = false);

/**
 * Returns Compute the centroid for a Polygon
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Polygon &g, bool enable3DComputation = false);

/**
 * Returns the centroid for a GeometryCollection
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const GeometryCollection &g, bool enable3DComputation = false);

/**
 * Returns the centroid for a TriangulatedSurface
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const TriangulatedSurface &g,
                 bool                       enable3DComputation = false);

/**
 * Returns the centroid for a PolyhedralSurface
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const PolyhedralSurface &g, bool enable3DComputation = false);

/**
 * Returns the centroid for a Solid
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Solid &g, bool enable3DComputation = false);

} // namespace algorithm
} // namespace SFCGAL

#endif
