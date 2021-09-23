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
package dilivia.s2.index

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.cell.S2CellIndex
import dilivia.s2.index.cell.S2ClosestCellQuery
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.index.shape.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2ShapeIndexRegion.Companion.makeS2ShapeIndexRegion
import org.apache.commons.math3.util.FastMath.sqrt

// This file defines a collection of classes that are useful for computing
// minimum distances on the sphere.  Their purpose is to allow code to be
// shared among the various query classes that find nearby geometry, such as
// S2ClosestEdgeQuery, S2ClosestPointQuery, and S2ClosestCellQuery.


// S2MinDistance is a thin wrapper around S1ChordAngle that is used by classes
// such as S2ClosestEdgeQuery to compute minimum distances on the sphere (as
// opposed to maximum distances, ellipsoidal distances, etc).
//
// It implements the Distance concept defined by S2DistanceTarget (see
// s2distance_target.h for details).
class S2MinDistance constructor(distance: S1ChordAngle) : Distance<S2MinDistance> {

    val value: S1ChordAngle = distance

    constructor() : this(S1ChordAngle.zero())
    constructor(x: S1Angle) : this(S1ChordAngle(S1ChordAngle(x).length2))

    override fun compareTo(other: S2MinDistance): Int = value.compareTo(other.value)
    override operator fun minus(delta: Delta): S2MinDistance = S2MinDistance(value - delta)
    override fun getChordAngleBound(): S1ChordAngle = value.plusError(value.getS1AngleConstructorMaxError())

    // If (dist < *this), updates *this and returns true (used internally).
    fun updateMin(dist: S2MinDistance): Boolean {
        if (dist.value < value) {
            value.length2 = dist.value.length2
            return true
        }
        return false
    }

    override fun toString(): String {
        return value.toString() // "(${value.toAngle().degrees()} degrees)"
    }

    override fun clone(): S2MinDistance = S2MinDistance(value.clone())

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S2MinDistance) return false

        if (value != other.value) return false

        return true
    }

    override fun hashCode(): Int {
        return value.hashCode()
    }


}

object S2MinDistanceFactory : DistanceFactory<S2MinDistance> {

    override fun distance(length2: Double): S2MinDistance = S2MinDistance(S1ChordAngle.fromLength2(length2))

    override fun distance(other: S2MinDistance): S2MinDistance = S2MinDistance(other.value.clone())

    override fun zero(): S2MinDistance = S2MinDistance(S1ChordAngle.zero())

    override fun infinity(): S2MinDistance = S2MinDistance(S1ChordAngle.infinity())

    override fun negative(): S2MinDistance = S2MinDistance(S1ChordAngle.negative())

}

// S2MinDistanceTarget represents a geometric object to which distances are
// measured.  Specifically, it is used to compute minimum distances on the
// sphere (as opposed to maximum distances, ellipsoidal distances, etc).
//
// Subtypes are defined below for measuring the distance to a point, an edge,
// an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
typealias S2MinDistanceTarget = S2DistanceTarget<S2MinDistance>

// An S2DistanceTarget subtype for computing the minimum distance to a point.
open class S2MinDistancePointTarget(val point: S2Point) : S2MinDistanceTarget {

    override fun getCapBound(): S2Cap = S2Cap(point, S1ChordAngle.zero())

    override fun distance(p: S2Point): S2MinDistance = S2MinDistance(S1ChordAngle.between(p, point))

    override fun updateMinDistance(p: S2Point, minDist: S2MinDistance): Boolean = minDist.updateMin(distance(p))

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MinDistance): Boolean =
        S2EdgeDistances.updateMinDistance(point, v0, v1, minDist.value)

    override fun updateMinDistance(cell: S2Cell, minDist: S2MinDistance): Boolean =
        minDist.updateMin(S2MinDistance(cell.getDistance(point)))

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        return makeS2ContainsPointQuery(queryIndex).visitContainingShapes(point) { shape ->
            visitor.visit(
                shape,
                point
            )
        }
    }


    override fun toString(): String {
        return "${this::class.java.simpleName}(point = ${S2PointUtil.toDegreesString(point)})"
    }

}

// An S2DistanceTarget subtype for computing the minimum distance to a edge.
open class S2MinDistanceEdgeTarget(val a: S2Point, val b: S2Point) : S2MinDistanceTarget {

    override fun getCapBound(): S2Cap {
        // The following computes a radius equal to half the edge length in an
        // efficient and numerically stable way.
        val d2 = S1ChordAngle.between(a, b).length2
        val r2 = (0.5 * d2) / (1 + sqrt(1 - 0.25 * d2))
        return S2Cap((a + b).normalize(), S1ChordAngle.fromLength2(r2))
    }

    override fun distance(p: S2Point): S2MinDistance = S2MinDistance(S2EdgeDistances.getDistance(p, a, b))

    override fun updateMinDistance(p: S2Point, minDist: S2MinDistance): Boolean =
        S2EdgeDistances.updateMinDistance(p, a, b, minDist.value)

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MinDistance): Boolean =
        S2EdgeDistances.updateEdgePairMinDistance(a, b, v0, v1, minDist.value)

    override fun updateMinDistance(cell: S2Cell, minDist: S2MinDistance): Boolean =
        minDist.updateMin(S2MinDistance(cell.getDistance(a, b)))

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // We test the center of the edge in order to ensure that edge targets AB
        // and BA yield identical results (which is not guaranteed by the API but
        // users might expect).  Other options would be to test both endpoints, or
        // return different results for AB and BA in some cases.
        val target = S2MinDistancePointTarget((a + b).normalize())
        return target.visitContainingShapes(queryIndex, visitor)
    }

    override fun toString(): String {
        return "${this::class.java.simpleName}(a=$a, b=$b)"
    }


}

// An S2DistanceTarget subtype for computing the minimum distance to an S2Cell
// (including the interior of the cell).
open class S2MinDistanceCellTarget(val cellTarget: S2Cell) : S2MinDistanceTarget {

    override fun getCapBound(): S2Cap = cellTarget.capBound

    override fun distance(p: S2Point): S2MinDistance = S2MinDistance(cellTarget.getDistance(p))

    override fun updateMinDistance(p: S2Point, minDist: S2MinDistance): Boolean = minDist.updateMin(distance(p))

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MinDistance): Boolean =
        minDist.updateMin(S2MinDistance(cellTarget.getDistance(v0, v1)))

    override fun updateMinDistance(cell: S2Cell, minDist: S2MinDistance): Boolean =
        minDist.updateMin(S2MinDistance(this.cellTarget.getDistance(cell)))

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // The simplest approach is simply to return the polygons that contain the
        // cell center.  Alternatively, if the index cell is smaller than the target
        // cell then we could return all polygons that are present in the
        // S2ShapeIndexCell, but since the index is built conservatively this may
        // include some polygons that don't quite intersect the cell.  So we would
        // either need to recheck for intersection more accurately, or weaken the
        // VisitContainingShapes contract so that it only guarantees approximate
        // intersection, neither of which seems like a good tradeoff.
        val target = S2MinDistancePointTarget(cellTarget.getCenter())
        return target.visitContainingShapes(queryIndex, visitor)
    }

    override fun toString(): String {
        return "${this::class.java.simpleName}(cell=$cellTarget)"
    }


}

// An S2DistanceTarget subtype for computing the minimum distance to an
// S2CellUnion (including the interior of all cells).
open class S2MinDistanceCellUnionTarget(val cellUnion: S2CellUnion) : S2MinDistanceTarget {

    private val index: S2CellIndex = S2CellIndex()
    private val query: S2ClosestCellQuery

    init {
        cellUnion.forEach { cellId -> index.add(cellId, 0) }
        index.build()
        query = S2ClosestCellQuery(index)
    }

    // Specifies that the distances should be computed by examining every cell
    // in the S2CellIndex (for testing and debugging purposes).
    //
    // DEFAULT: false
    fun useBruteForce(): Boolean = query.options.useBruteForce
    fun setUseBruteForce(useBruteForce: Boolean) {
        query.options.useBruteForce = useBruteForce
    }

    // Note that set_max_error() should not be called directly by clients; it is
    // used internally by the S2Closest*Query implementations.

    override fun setMaxError(max_error: Delta): Boolean {
        query.options.maxError = max_error
        return true
    }

    override fun getCapBound(): S2Cap = cellUnion.capBound

    override fun distance(p: S2Point): S2MinDistance {
        TODO("Not yet implemented")
    }

    override fun updateMinDistance(p: S2Point, minDist: S2MinDistance): Boolean {
        val target = S2ClosestCellQuery.PointTarget(p)
        return updateMinDistance(target, minDist)
    }

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MinDistance): Boolean {
        val target = S2ClosestCellQuery.EdgeTarget(v0, v1)
        return updateMinDistance(target, minDist)
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MinDistance): Boolean {
        val target = S2ClosestCellQuery.CellTarget(cell)
        return updateMinDistance(target, minDist)
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        for (cell_id in cellUnion) {
            val target = S2MinDistancePointTarget(cell_id.toPoint())
            if (!target.visitContainingShapes(queryIndex, visitor)) {
                return false;
            }
        }
        return true;
    }

    private fun updateMinDistance(target: S2MinDistanceTarget, min_dist: S2MinDistance): Boolean {
        query.options.maxDistance = min_dist
        val r = query.findClosestCell(target)
        return if (r.isEmpty()) {
            false
        } else {
            min_dist.value.length2 = r.distance.value.length2
            true
        }
    }

    override fun toString(): String {
        return "${this::class.java.simpleName}(cellUnion=$cellUnion)"
    }


}

// An S2DistanceTarget subtype for computing the minimum distance to an
// S2ShapeIndex (a collection of points, polylines, and/or polygons).
//
// Note that ShapeIndexTarget has its own options:
//
//   include_interiors()
//     - specifies that distances are measured to the boundary and interior
//       of polygons in the S2ShapeIndex.  (If set to false, distance is
//       measured to the polygon boundary only.)
//       DEFAULT: true.
//
//   brute_force()
//     - specifies that the distances should be computed by examining every
//       edge in the S2ShapeIndex (for testing and debugging purposes).
//       DEFAULT: false.
//
// These options are specified independently of the corresponding
// S2ClosestEdgeQuery options.  For example, if include_interiors is true for
// a ShapeIndexTarget but false for the S2ClosestEdgeQuery where the target
// is used, then distances will be measured from the boundary of one
// S2ShapeIndex to the boundary and interior of the other.
//
// Note that when the distance to a ShapeIndexTarget is zero because the
// target intersects the interior of the query index, you can find a point
// that achieves this zero distance by calling the VisitContainingShapes()
// method directly.  For example:
//
//   S2ClosestEdgeQuery::ShapeIndexTarget target(&target_index);
//   target.VisitContainingShapes(
//       query_index, [](S2Shape* containing_shape,
//                       const S2Point& target_point) {
//         ... do something with "target_point" ...
//         return false;  // Terminate search
//       }));
open class S2MinDistanceShapeIndexTarget(val index: S2ShapeIndex) : S2MinDistanceTarget {

    private val query: S2ClosestEdgeQuery = S2ClosestEdgeQuery(index)

    // Specifies that distance will be measured to the boundary and interior
    // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
    //
    // DEFAULT: true
    fun includeInteriors(): Boolean = query.options.includeInteriors
    fun setIncludeInteriors(include_interiors: Boolean): Unit {
        query.options.includeInteriors = include_interiors
    }

    // Specifies that the distances should be computed by examining every edge
    // in the S2ShapeIndex (for testing and debugging purposes).
    //
    // DEFAULT: false
    fun useBruteForce(): Boolean = query.options.useBruteForce
    fun setUseBruteForce(useBruteForce: Boolean): Unit {
        query.options.useBruteForce = useBruteForce
    }

    // Note that set_max_error() should not be called directly by clients; it is
    // used internally by the S2Closest*Query implementations.
    override fun setMaxError(max_error: Delta): Boolean {
        query.options.maxError = max_error
        return true  // Indicates that we may return suboptimal results.
    }

    override fun getCapBound(): S2Cap = makeS2ShapeIndexRegion(index).capBound

    override fun distance(p: S2Point): S2MinDistance {
        TODO("Not yet implemented")
    }

    override fun updateMinDistance(p: S2Point, min_dist: S2MinDistance): Boolean {
        val target = S2ClosestEdgeQuery.PointTarget(p)
        return updateMinDistance(target, min_dist)
    }

    override fun updateMinDistance(v0: S2Point, v1: S2Point, min_dist: S2MinDistance): Boolean {
        val target = S2ClosestEdgeQuery.EdgeTarget(v0, v1)
        return updateMinDistance(target, min_dist)
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MinDistance): Boolean {
        val target = S2ClosestEdgeQuery.CellTarget(cell)
        return updateMinDistance(target, minDist)
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // It is sufficient to find the set of chain starts in the target index
        // (i.e., one vertex per connected component of edges) that are contained by
        // the query index, except for one special case to handle full polygons.
        //
        // TODO(ericv): Do this by merge-joining the two S2ShapeIndexes, and share
        // the code with S2BooleanOperation.

        for (shape in index) {
            if (shape == null) continue;
            val num_chains = shape.numChains
            // Shapes that don't have any edges require a special case (below).
            var tested_point = false
            for (c in 0 until num_chains) {
                val chain = shape.chain(c)
                if (chain.length == 0) continue
                tested_point = true;
                val v0 = shape.chainEdge(c, 0).v0
                val target = S2MinDistancePointTarget(v0)
                if (!target.visitContainingShapes(queryIndex, visitor)) {
                    return false;
                }
            }
            if (!tested_point) {
                // Special case to handle full polygons.
                val ref = shape.getReferencePoint()
                if (!ref.contained) continue
                val target = S2MinDistancePointTarget(ref.point)
                if (!target.visitContainingShapes(queryIndex, visitor)) {
                    return false
                }
            }
        }
        return true
    }


    private fun updateMinDistance(target: S2MinDistanceTarget, min_dist: S2MinDistance): Boolean {
        query.options.maxDistance = min_dist
        val r = query.findClosestEdge(target)
        if (r.isEmpty()) return false
        min_dist.value.length2 = r.distance.value.length2
        return true
    }

    override fun toString(): String {
        return "${this::class.java.simpleName}()"
    }
}

