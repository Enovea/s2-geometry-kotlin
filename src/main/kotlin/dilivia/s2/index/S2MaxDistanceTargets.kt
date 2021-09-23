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
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import org.apache.commons.math3.util.FastMath.sqrt


class S2FurthestEdgeQuery(val index: S2ShapeIndex)

// S2MaxDistance is a class that allows maximum distances to be computed using
// a minimum distance algorithm.  Specifically, S2MaxDistance(x) represents the
// supplementary distance (Pi - x).  This has the effect of inverting the sort
// order, i.e.
//
//  (S2MaxDistance(x) < S2MaxDistance(y))  <=>  (Pi - x < Pi - y)  <=>  (x > y)
//
// All other operations are implemented similarly (using the supplementary
// distance Pi - x).  For example, S2MaxDistance(x) - S2MaxDistance(y) ==
// S2MaxDistance(x + y).
class S2MaxDistance constructor(distance: S1ChordAngle) : Distance<S2MaxDistance> {

    val value: S1ChordAngle = distance.clone()

    constructor() : this(S1ChordAngle.zero())
    constructor(x: S1Angle) : this(S1ChordAngle(S1ChordAngle(x).length2))

    override fun compareTo(other: S2MaxDistance): Int = -value.compareTo(other.value)
    override operator fun minus(delta: Delta): S2MaxDistance = S2MaxDistance(value + delta)
    override fun getChordAngleBound(): S1ChordAngle = S1ChordAngle.straight() - value

    // If (dist < *this), updates *this and returns true (used internally).
    fun updateMin(dist: S2MaxDistance): Boolean {
        if (dist < this) {
            this.value.length2 = dist.value.length2
            return true
        }
        return false
    }

    override fun clone(): S2MaxDistance = S2MaxDistance(value)

}

object S2MaxDistanceFactory : DistanceFactory<S2MaxDistance> {

    override fun distance(length2: Double): S2MaxDistance = S2MaxDistance(S1ChordAngle.fromLength2(length2))

    override fun distance(other: S2MaxDistance): S2MaxDistance = S2MaxDistance(other.value)

    override fun zero(): S2MaxDistance = S2MaxDistance(S1ChordAngle.straight())

    override fun infinity(): S2MaxDistance = S2MaxDistance(S1ChordAngle.negative())

    override fun negative(): S2MaxDistance = S2MaxDistance(S1ChordAngle.infinity())

}

// S2MaxDistanceTarget represents a geometric object to which maximum distances
// on the sphere are measured.
//
// Subtypes are defined below for measuring the distance to a point, an edge,
// an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
typealias S2MaxDistanceTarget = S2DistanceTarget<S2MaxDistance>

// An S2DistanceTarget subtype for computing the maximum distance to a point.
open class S2MaxDistancePointTarget(val point: S2Point) : S2MaxDistanceTarget {

    // This method returns an S2Cap that bounds the antipode of the target.  (This
    // is the set of points whose S2MaxDistance to the target is
    // S2MaxDistance::Zero().)
    override fun getCapBound(): S2Cap = S2Cap(-point, S1ChordAngle.zero())

    override fun distance(p: S2Point): S2MaxDistance = S2MaxDistance(S1ChordAngle.between(p, point))

    override fun updateMinDistance(p: S2Point, minDist: S2MaxDistance): Boolean = minDist.updateMin(distance(p))

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MaxDistance): Boolean {
        val dist = S1ChordAngle(minDist.value)
        if (S2EdgeDistances.updateMaxDistance(point, v0, v1, dist)) {
            minDist.updateMin(S2MaxDistance(dist))
            return true
        }
        return false
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MaxDistance): Boolean {
        return minDist.updateMin(S2MaxDistance(cell.getMaxDistance(point)))
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // For furthest points, we visit the polygons whose interior contains the
        // antipode of the target point.  (These are the polygons whose
        // S2MaxDistance to the target is S2MaxDistance::Zero().)
        TODO()
//        return makeS2ContainsPointQuery(query_index).visitContainingShapes(-point, [this, &visitor](S2Shape* shape) {
//            return visitor(shape, point_);
//        });
    }

}

// An S2DistanceTarget subtype for computing the maximum distance to an edge.
class S2MaxDistanceEdgeTarget(a: S2Point, b: S2Point) : S2MaxDistanceTarget {

    val a: S2Point = a.normalized()
    val b: S2Point = b.normalized()

    override fun getCapBound(): S2Cap {
        // The following computes a radius equal to half the edge length in an
        // efficient and numerically stable way.
        val d2 = S1ChordAngle.between(a, b).length2
        val r2 = (0.5 * d2) / (1 + sqrt(1 - 0.25 * d2))
        return S2Cap(-(a + b).normalize(), S1ChordAngle.fromLength2(r2))
    }

    override fun distance(p: S2Point): S2MaxDistance = S2MaxDistance(S2EdgeDistances.getDistance(p, a, b))

    override fun updateMinDistance(p: S2Point, minDist: S2MaxDistance): Boolean {
        val dist = S1ChordAngle(minDist.value)
        if (S2EdgeDistances.updateMaxDistance(p, a, b, dist)) {
            minDist.updateMin(S2MaxDistance(dist))
            return true
        }
        return false
    }

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MaxDistance): Boolean {
        val dist = S1ChordAngle(minDist.value)
        if (S2EdgeDistances.updateEdgePairMaxDistance(a, b, v0, v1, dist)) {
            minDist.updateMin(S2MaxDistance(dist))
            return true
        }
        return false
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MaxDistance): Boolean {
        return minDist.updateMin(S2MaxDistance(cell.getMaxDistance(a, b)))
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // We only need to test one edge point.  That is because the method *must*
        // visit a polygon if it fully contains the target, and *is allowed* to
        // visit a polygon if it intersects the target.  If the tested vertex is not
        // contained, we know the full edge is not contained; if the tested vertex is
        // contained, then the edge either is fully contained (must be visited) or it
        // intersects (is allowed to be visited).  We visit the center of the edge so
        // that edge AB gives identical results to BA.
        val target = S2MaxDistancePointTarget((a + b).normalize())
        return target.visitContainingShapes(queryIndex, visitor)
    }

}

// An S2DistanceTarget subtype for computing the maximum distance to an S2Cell
// (including the interior of the cell).
class S2MaxDistanceCellTarget(val cell: S2Cell) : S2MaxDistanceTarget {

    override fun getCapBound(): S2Cap {
        val cap = cell.capBound
        return S2Cap(-cap.center, cap.radius)
    }

    override fun distance(p: S2Point): S2MaxDistance = S2MaxDistance(cell.getMaxDistance(p))

    override fun updateMinDistance(p: S2Point, minDist: S2MaxDistance): Boolean {
        return minDist.updateMin(S2MaxDistance(cell.getMaxDistance(p)))
    }

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MaxDistance): Boolean {
        return minDist.updateMin(S2MaxDistance(cell.getMaxDistance(v0, v1)))
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MaxDistance): Boolean {
        return minDist.updateMin(S2MaxDistance(cell.getMaxDistance(cell)))
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        // We only need to check one point here - cell center is simplest.
        // See comment at S2MaxDistanceEdgeTarget::VisitContainingShapes.
        val target = S2MaxDistancePointTarget(cell.getCenter())
        return target.visitContainingShapes(queryIndex, visitor)
    }

}

// An S2DistanceTarget subtype for computing the maximum distance to an
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
// S2FurthestEdgeQuery options.  For example, if include_interiors is true for
// a ShapeIndexTarget but false for the S2FurthestEdgeQuery where the target
// is used, then distances will be measured from the boundary of one
// S2ShapeIndex to the boundary and interior of the other.
//
class S2MaxDistanceShapeIndexTarget(val index: S2ShapeIndex) : S2MaxDistanceTarget {

    private val query: S2FurthestEdgeQuery = S2FurthestEdgeQuery(index)


    // Specifies that distance will be measured to the boundary and interior
    // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
    //
    // DEFAULT: true
    fun includeInteriors(): Boolean = TODO()
    fun setIncludeInteriors(includeInteriors: Boolean): Unit = TODO()

    // Specifies that the distances should be computed by examining every edge
    // in the S2ShapeIndex (for testing and debugging purposes).
    //
    // DEFAULT: false
    fun useBruteForce(): Boolean = TODO()
    fun setUseBruteForce(useBruteForce: Boolean): Unit = TODO()

    override fun setMaxError(max_error: Delta): Boolean {
        TODO()
    }

    // This method returns an S2Cap that bounds the antipode of the target.  (This
    // is the set of points whose S2MaxDistance to the target is
    // S2MaxDistance::Zero().)
    override fun getCapBound(): S2Cap {
        TODO()
//        val cap = makeS2ShapeIndexRegion(index).capBound
//        return S2Cap(-cap.center, cap.radius)
    }

    override fun distance(p: S2Point): S2MaxDistance {
        TODO("Not yet implemented")
    }

    override fun updateMinDistance(p: S2Point, minDist: S2MaxDistance): Boolean {
        TODO("Not yet implemented")
    }

    override fun updateMinDistance(v0: S2Point, v1: S2Point, minDist: S2MaxDistance): Boolean {
        TODO("Not yet implemented")
    }

    override fun updateMinDistance(cell: S2Cell, minDist: S2MaxDistance): Boolean {
        TODO("Not yet implemented")
    }

    override fun visitContainingShapes(queryIndex: S2ShapeIndex, visitor: S2DistanceTarget.ShapeVisitor): Boolean {
        TODO("Not yet implemented")
    }

}
