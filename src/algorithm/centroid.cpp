// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/centroid.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

#include "SFCGAL/Exception.h"
#include <boost/format.hpp>

namespace SFCGAL::algorithm {

using Point_2    = CGAL::Point_2<SFCGAL::Kernel>;
using Triangle_2 = CGAL::Triangle_2<SFCGAL::Kernel>;
using Polygon_2  = CGAL::Polygon_2<SFCGAL::Kernel>;
using Vector_2   = CGAL::Vector_2<SFCGAL::Kernel>;

using Point_3    = CGAL::Point_3<SFCGAL::Kernel>;
using Triangle_3 = CGAL::Triangle_3<SFCGAL::Kernel>;
using Plane_3    = CGAL::Plane_3<SFCGAL::Kernel>;
using Vector_3   = CGAL::Vector_3<SFCGAL::Kernel>;

///
///
///
template <typename T>
std::function<WeightedCentroid(const T &)> weightedCentroidLambda =
    [](const T &g) -> WeightedCentroid {
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (typename T::const_iterator ite = g.begin(); ite != g.end(); ite++) {
    // compute geometry area+centroid
    WeightedCentroid wc = weightedCentroid(*ite);

    // update totals
    if (ite == g.begin())
      totalWeightedCentroid = wc.area * wc.centroid;
    else
      totalWeightedCentroid += wc.area * wc.centroid;
    totalM += wc.area * wc.m;
    totalArea += wc.area;
  }

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
};

///
///
///
auto
centroid(const Geometry &g) -> std::unique_ptr<Point>
{
  if (g.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        InappropriateGeometryException("No point in geometry."));
  }

  WeightedCentroid wCent = weightedCentroid(g);

  Point out;
  if (g.is3D())
    out = Point(wCent.centroid.x(), wCent.centroid.y(), wCent.centroid.z());
  else
    out = Point(wCent.centroid.x(), wCent.centroid.y());

  if (g.isMeasured())
    out.setM(CGAL::to_double(wCent.m));

  return std::make_unique<Point>(out);
}

auto
weightedCentroid(const Geometry &g) -> WeightedCentroid
{
  WeightedCentroid wCent;
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    wCent = WeightedCentroid(0.0, g.as<Point>().toVector_3());
    break;
  case TYPE_LINESTRING:
    wCent = weightedCentroid(g.as<LineString>());
    break;

  case TYPE_POLYGON:
    wCent = weightedCentroid(g.as<Polygon>());
    break;

  case TYPE_TRIANGLE:
    wCent = weightedCentroid(g.as<Triangle>());
    break;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    wCent = weightedCentroid(g.as<GeometryCollection>());
    break;

  case TYPE_TRIANGULATEDSURFACE:
    wCent = weightedCentroid(g.as<TriangulatedSurface>());
    break;

  case TYPE_POLYHEDRALSURFACE:
    // wCent = weightedCentroid(g.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        "Centroid of solid gemotry is not implemented."));
  }

  return wCent;
  BOOST_THROW_EXCEPTION(
      Exception((boost::format("Unexpected geometry type (%s) in "
                               "SFCGAL::algorithm::weightedCentroid") %
                 g.geometryType())
                    .str()));
}

///
///
///
auto
weightedCentroid(const Triangle &g) -> WeightedCentroid
{
  g.toTriangle_2().area();
  return weightedCentroid(g.vertex(0), g.vertex(1), g.vertex(2));
}

///
///
///
auto
weightedCentroid(const Point &a, const Point &b, const Point &c)
    -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(Triangle)\n";
  // compute triangle area
  SFCGAL::Kernel::FT area = SFCGAL::Kernel().compute_area_2_object()(
      a.toPoint_2(), b.toPoint_2(), c.toPoint_2());

  // compute triangle centroid
  Vector_3 out = c.toVector_3();
  out += b.toVector_3();
  out += a.toVector_3();
  out /= 3.0;

  SFCGAL::Kernel::FT m = 0.0;
  if (a.isMeasured() && b.isMeasured() && c.isMeasured()) {
    m = (a.m() + b.m() + c.m()) / 3.0;
  }

  // std::cout << "Triangle weightedCentroid: " << out << "\n";
  // std::cout << "Triangle area: " << area << "\n";

  return {area, out, m};
}

///
///
///
auto
weightedCentroid(const LineString &g) -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(LineString)\n";
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  if (g.isClosed()) { // ie. a polygon
    for (size_t i = 1; i < g.numPoints() - 2; i++) {
      // compute triangle area+centroid
      WeightedCentroid wc =
          weightedCentroid(g.pointN(0), g.pointN(i), g.pointN(i + 1));

      // update totals
      if (i == 1)
        totalWeightedCentroid = wc.area * wc.centroid;
      else
        totalWeightedCentroid += wc.area * wc.centroid;
      totalM += wc.area * wc.m;
      totalArea += wc.area;
      // std::cout << "Polygon weightedCentroid: " << wc.area * wc.centroid <<
      // "\n"; std::cout << "Polygon area: " << wc.area << "\n"; std::cout <<
      // "Polygon totalweightedCentroid: " << totalWeightedCentroid << "\n";
      // std::cout << "Polygon totalArea: " << totalArea << "\n";
    }

  } else { // ie. a linestring
    for (size_t i = 0; i < g.numPoints() - 1; i++) {
      Point    a  = g.pointN(i);
      Point    b  = g.pointN(i + 1);
      Vector_3 ba = b.toPoint_3() - a.toPoint_3();

      // TODO: 3D length! should be 2D?
      SFCGAL::Kernel::FT len = CGAL::sqrt(CGAL::to_double(ba.squared_length()));
      if (i == 0)
        totalWeightedCentroid = (a.toVector_3() + (ba / 2.0)) * len;
      else
        totalWeightedCentroid += (a.toVector_3() + (ba / 2.0)) * len;
      if (a.isMeasured() && b.isMeasured())
        totalM += (a.m() + b.m()) * 0.5 * len;
      totalArea += len;
      // std::cout << "LineString weightedCentroid: " << (a.toVector_3() + (ba
      // / 2.0)) * len << "\n"; std::cout << "LineString len: " << len << "\n";
      // std::cout << "LineString totalweightedCentroid: " <<
      // totalWeightedCentroid << "\n"; std::cout << "LineString totalArea: " <<
      // totalArea << "\n";
    }
  }

  if (totalArea == 0.0)
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(
        "SFCGAL::algorithm::Centroid of LineString without 2D area is not "
        "valid."));

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
}

///
///
///
auto
weightedCentroid(const Polygon &g) -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(Polygon)\n";
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (size_t i = 0; i < g.numRings(); i++) {
    try {
      WeightedCentroid ringCentroid = weightedCentroid(g.ringN(i));
      // Kernel::FT const ringCentroid = CGAL::abs( signedCentroid( g.ringN( i )
      // ) );

      if (i == 0) {
        // exterior ring
        totalArea             = ringCentroid.area;
        totalM                = ringCentroid.area * ringCentroid.m;
        totalWeightedCentroid = ringCentroid.area * ringCentroid.centroid;
      } else {
        // interior ring
        totalArea -= ringCentroid.area;
        totalM -= ringCentroid.area * ringCentroid.m;
        totalWeightedCentroid -= ringCentroid.area * ringCentroid.centroid;
      }
      // std::cout << "sub Polygon weightedCentroid: " << ringCentroid.area *
      // ringCentroid.centroid
      //           << "\n";
      // std::cout << "sub Polygon area: " << ringCentroid.area << "\n";
      // std::cout << "sub Polygon totalweightedCentroid: " <<
      // totalWeightedCentroid << "\n"; std::cout << "sub Polygon totalArea: "
      // << totalArea << "\n";

    } catch (InappropriateGeometryException &e) {
      BOOST_THROW_EXCEPTION(InappropriateGeometryException(
          (boost::format("SFCGAL::algorithm::Centroid of Polygon failed at "
                         "ring %d with cause: \n%s") %
           i % e.what())
              .str()));
    }
  }

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
}

///
///
///
auto
weightedCentroid(const GeometryCollection &g) -> WeightedCentroid
{
  return weightedCentroidLambda<GeometryCollection>(g);
}

///
///
///
auto
weightedCentroid(const TriangulatedSurface &g) -> WeightedCentroid
{
  return weightedCentroidLambda<TriangulatedSurface>(g);
}

///
///
///
auto
weightedCentroid(const PolyhedralSurface &g) -> WeightedCentroid
{
  return weightedCentroidLambda<PolyhedralSurface>(g);
}

///----------------------------------------------------------------------------------
/// -- centroid3D
///----------------------------------------------------------------------------------

///
///
///
// auto
//   centroid3D( const Geometry &g, NoValidityCheck /*unused*/ ) ->
//   std::unique_ptr<Point>
// {
//   switch ( g.geometryTypeId() )
//   {
//     case TYPE_POINT:
//     case TYPE_LINESTRING:
//       return 0;

//     case TYPE_POLYGON:
//       return centroid3D( g.as<Polygon>() );

//     case TYPE_TRIANGLE:
//       return centroid3D( g.as<Triangle>() );

//     case TYPE_MULTIPOINT:
//     case TYPE_MULTILINESTRING:
//     case TYPE_MULTIPOLYGON:
//     case TYPE_GEOMETRYCOLLECTION:
//       return centroid3D( g.as<GeometryCollection>() );

//     case TYPE_TRIANGULATEDSURFACE:
//       return centroid3D( g.as<TriangulatedSurface>() );

//     case TYPE_POLYHEDRALSURFACE:
//       return centroid3D( g.as<PolyhedralSurface>() );

//     case TYPE_SOLID:
//     case TYPE_MULTISOLID:
//       return 0;
//   }

//   BOOST_THROW_EXCEPTION( Exception( "missing case in
//   SFCGAL::algorithm::centroid3D" ) );
// }

// auto
//   centroid3D( const Geometry &g ) -> std::unique_ptr<Point>
// {
//   SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D( g );
//   return centroid3D( g, NoValidityCheck() );
// }
// ///
// ///
// ///
// auto
//   centroid3D( const Polygon &g ) -> std::unique_ptr<Point>
// {
//   double result = 0.0;

//   if ( g.isEmpty() )
//   {
//     return result;
//   }

//   Point_3 a;
//   Point_3 b;
//   Point_3 c;
//   algorithm::plane3D<Kernel>( g, a, b, c );

//   /*
//    * compute polygon basis (CGAL doesn't build an orthonormal basis so that
//    * computing the 2D centroid in this basis would lead to scale effects) ux
//    = bc uz
//    * = bc^ba uy = uz^ux
//    *
//    * Note that the basis is rounded to double (CGAL::sqrt)
//    */
//   Point_3 ux = c - b;
//   Point_3 uz = CGAL::cross_product( ux, a - b );
//   ux = ux / CGAL::sqrt( CGAL::to_double( ux.squared_length() ) );
//   uz = uz / CGAL::sqrt( CGAL::to_double( uz.squared_length() ) );
//   Point_3 const uy = CGAL::cross_product( uz, ux );

//   /*
//    * compute the centroid for each ring in the local basis
//    */
//   for ( size_t i = 0; i < g.numRings(); i++ )
//   {
//     const LineString &ring = g.ringN( i );

//     CGAL::Polygon_2<Kernel> projectedPolygon;

//     for ( size_t j = 0; j < ring.numPoints() - 1; j++ )
//     {
//       Point_3 const point = ring.pointN( j ).toPoint_3();
//       CGAL::Point_2<Kernel> const projectedPoint( ( point - b ) * ux, ( point
//       - b ) * uy ); projectedPolygon.push_back( projectedPoint );
//     }

//     if ( i == 0 )
//     {
//       // exterior ring
//       result += CGAL::to_double( CGAL::abs( projectedPolygon.centroid() ) );
//     }
//     else
//     {
//       // interior ring
//       result -= CGAL::to_double( CGAL::abs( projectedPolygon.centroid() ) );
//     }
//   }

//   return result;
// }

// ///
// ///
// ///
// auto
//   centroid3D( const Triangle &g ) -> std::unique_ptr<Point>
// {
//   CGAL::Triangle_3<Kernel> const triangle( g.vertex( 0 ).toPoint_3(),
//   g.vertex( 1 ).toPoint_3(), g.vertex( 2 ).toPoint_3() ); return sqrt(
//   CGAL::to_double( triangle.squared_centroid() ) );
// }

// ///
// ///
// ///
// auto
//   centroid3D( const GeometryCollection &g ) -> std::unique_ptr<Point>
// {
//   double result = 0.0;

//   for ( size_t i = 0; i < g.numGeometries(); i++ )
//   {
//     result += centroid3D( g.geometryN( i ) );
//   }

//   return result;
// }

// ///
// ///
// ///
// auto
//   centroid3D( const PolyhedralSurface &g ) -> std::unique_ptr<Point>
// {
//   double centroid = 0.0;

//   for ( size_t i = 0; i < g.numPolygons(); i++ )
//   {
//     centroid += centroid3D( g.polygonN( i ) );
//   }

//   return centroid;
// }

// ///
// ///
// ///
// auto
//   centroid3D( const TriangulatedSurface &g ) -> std::unique_ptr<Point>
// {
//   double result = 0.0;

//   for ( size_t i = 0; i < g.numGeometries(); i++ )
//   {
//     result += centroid3D( g.geometryN( i ) );
//   }

//   return result;
// }

} // namespace SFCGAL::algorithm
