/*
 * Copyright © 2021 Enovea (fabien.meurisse@enovea.net)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package dilivia.s2.region

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.coords.S2Coords
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2ShapeIndexRegion.Companion.makeS2ShapeIndexRegion
import org.apache.commons.math3.util.FastMath.min

// This class provides a way to expand an arbitrary collection of geometry by
// a fixed radius (an operation variously known as "buffering", "offsetting",
// or "Minkowski sum with a disc") in order to compute an S2CellId covering
// (see S2RegionCoverer).  The resulting covering contains all points within
// the given radius of any point in the original geometry.
//
// This class does not actually buffer the geometry; instead it implements the
// S2Region API by computing the distance from candidate S2CellIds to the
// original geometry.  If this distance is below the given radius then the
// S2CellId intersects the buffered geometry.  For example, if the original
// geometry consists of a single S2Point then the buffered geometry is exactly
// equivalent to an S2Cap with the given radius.  (Note that the region is not
// approximated as a polygonal loop.)
//
// Example usage:
//
// S2CellUnion GetBufferedCovering(const S2ShapeIndex& index, S1Angle radius) {
//   S2RegionCoverer coverer;
//   coverer.mutable_options()->set_max_cells(20);
//   S2CellUnion covering;
//   S2ShapeIndexBufferedRegion region(&index, radius);
//   coverer.GetCovering(region, &covering);
//   return covering;
// }
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
class S2ShapeIndexBufferedRegion : S2Region {

    private lateinit var radius: S1ChordAngle

    // In order to handle (radius_ == 0) corectly, we need to test whether
    // distances are less than or equal to "radius_".  This is done by testing
    // whether distances are less than radius_.Successor().
    private lateinit var radiusSuccessor: S1ChordAngle

    private lateinit var query: S2ClosestEdgeQuery  // This class is not thread-safe!

  // Default constructor; requires Init() to be called.
  constructor()

  // Constructs a region representing all points within the given radius of
  // any point in the given S2ShapeIndex.
  constructor(index: S2ShapeIndex, radius: S1ChordAngle) {
      this.radius = radius
      this.radiusSuccessor = radius.successor()
      val options = S2ClosestEdgeQuery.Options()
      options.includeInteriors = true
      this.query = S2ClosestEdgeQuery(index, options)
  }

  // Convenience constructor that accepts an S1Angle for the radius.
  // REQUIRES: radius >= S1Angle::Zero()
  constructor(index: S2ShapeIndex, radius: S1Angle): this(index, S1ChordAngle(radius))

  // Equivalent to the constructor above.
  fun init(index: S2ShapeIndex, radius: S1ChordAngle) {
      this.radius = radius
      this.radiusSuccessor = radius.successor()
      query.init(index)
      query.options.includeInteriors = true
  }

  fun index(): S2ShapeIndex = query.index()
  fun radius(): S1ChordAngle = radius

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  // Clone() returns a *shallow* copy; it does not make a copy of the
  // underlying S2ShapeIndex.

    override fun clone(): S2ShapeIndexBufferedRegion = S2ShapeIndexBufferedRegion(index(), radius)

    override val capBound: S2Cap
        get() {
            val origCap = makeS2ShapeIndexRegion(index()).capBound
            return S2Cap(origCap.center, origCap.radius + radius)
        }

    override val rectBound: S2LatLngRect
        get() {
            val origRect = makeS2ShapeIndexRegion(index()).rectBound
            return origRect.expandedByDistance(radius.toAngle())
        }

  // This method returns a small non-optimal covering that may include
  // duplicate or overlapping cells.  It should not be used directly.
  // Instead, use S2RegionCoverer::GetCovering or GetFastCovering.
    override fun getCellUnionBound(cellIds: MutableList<S2CellId>) {
      // We start with a covering of the original S2ShapeIndex, and then expand it
      // by replacing each cell with a block of 4 cells whose union contains the
      // original cell buffered by the given radius.
      //
      // This increases the number of cells in the covering by a factor of 4 and
      // increases the covered area by a factor of 16, so it is not a very good
      // covering, but it is much better than always returning the 6 face cells.
      val orig_cellids = mutableListOf<S2CellId>()
      makeS2ShapeIndexRegion(index()).getCellUnionBound(orig_cellids)

      val radians = radius.toAngle().radians
      val max_level = S2Coords.projection.kMinWidth.getLevelForMinValue(radians) - 1
      if (max_level < 0) {
          return S2Cap.full.getCellUnionBound(cellIds)
      }
      cellIds.clear()
      for (id in orig_cellids) {
          if (id.isFace) {
              return S2Cap.full.getCellUnionBound(cellIds)
          }
          id.appendVertexNeighbors(min(max_level, id.level() - 1), cellIds)
      }
    }

  // The implementation is approximate but conservative; it always returns
  // "false" if the cell is not contained by the buffered region, but it may
  // also return false in some cases where "cell" is in fact contained.
    override fun contains(cell: S2Cell): Boolean {
      // To implement this method perfectly would require computing the directed
      // Hausdorff distance, which is expensive (and not currently implemented).
      // However the following heuristic is almost as good in practice and much
      // cheaper to compute.

      // Return true if the unbuffered region contains this cell.
      if (makeS2ShapeIndexRegion(index()).contains(cell)) return true

      // Otherwise approximate the cell by its bounding cap.
      //
      // NOTE(ericv): It would be slightly more accurate to first find the closest
      // point in the indexed geometry to the cell, and then measure the actual
      // maximum distance from that point to the cell (a poor man's Hausdorff
      // distance).  But based on actual tests this is not worthwhile.
      val cap = cell.capBound
      if (radius < cap.radius) return false

      // Return true if the distance to the cell center plus the radius of the
      // cell's bounding cap is less than or equal to "radius_".
      val target = S2ClosestEdgeQuery.PointTarget(cell.getCenter())
      return query.isDistanceLess(target, radiusSuccessor - cap.radius)
    }

  // Returns true if any buffered shape intersects "cell" (to within a very
  // small error margin).
    override fun mayIntersect(cell: S2Cell): Boolean {
      // Return true if the distance is less than or equal to "radius_".
      val target = S2ClosestEdgeQuery.CellTarget(cell)
      return query.isDistanceLess(target, radiusSuccessor)
    }

  // Returns true if the given point is contained by the buffered region,
  // i.e. if it is within the given radius of any original shape.
    override fun contains(p: S2Point): Boolean {
      val target = S2ClosestEdgeQuery.PointTarget(p)
      // Return true if the distance is less than or equal to "radius_".
      return query.isDistanceLess(target, radiusSuccessor)
    }

}
