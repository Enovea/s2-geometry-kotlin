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

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireGE
import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.R1Interval
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.radians
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S1ChordAngle.Companion.sin
import dilivia.s2.S1Interval
import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeDistances
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.asin
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.max


/**
 * S2Cap represents a disc-shaped region defined by a center and radius. Technically this shape is called a
 * "spherical cap" (rather than disc) because it is not planar; the cap represents a portion of the sphere that has
 * been cut off by a plane.  The boundary of the cap is the circle defined by the intersection of the sphere and the
 * plane. For containment purposes, the cap is a closed set, i.e. it contains its boundary.
 *
 * For the most part, you can use a spherical cap wherever you would use a disc in planar geometry.  The radius of the
 * cap is measured along the surface of the sphere (rather than the straight-line distance through the interior).
 * Thus a cap of radius Pi/2 is a hemisphere, and a cap of radius Pi covers the entire sphere.
 *
 * A cap can also be defined by its center point and height. The height is simply the distance from the center point
 * to the cutoff plane. There is also support for empty and full caps, which contain no points and all points
 * respectively.
 *
 * Here are some useful relationships between the cap height (h), the cap opening angle (theta), the maximum chord
 * length from the cap's center (d), and the radius of cap's base (a). All formulas assume a unit radius.
 *
 * h = 1 - cos(theta) = 2 sin^2(theta/2) d^2 = 2 h = a^2 + h^2
 *
 * This class is a port of the S2Cap class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @property center The center of the cap (a point on the unit sphere).
 * @property radius The radius of the cap given as a chord of the sphere.
 * @constructor Constructs a cap where the angle is expressed as an S1ChordAngle.  This constructor is more efficient
 * than the one with S1Angle.
 * @param center The center of the cap.
 * @param radius The radius of the cap.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2Cap(val center: S2Point, val radius: S1ChordAngle) : S2Region {

    init {
        checkState { isValid }
    }

    /** The default constructor returns an empty S2Cap. */
    constructor() : this(S2Point(1.0, 0.0, 0.0), S1ChordAngle.negative())

    /**
     * Constructs a cap with the given center and radius. A negative radius yields an empty cap; a radius of 1
     * 80 degrees or more yields a full cap (containing the entire sphere).  "center" should be unit length.
     * The radius is computed as the min of radius and S2.M_PI.
     *
     * @param center The center of the cap (must be on the unit sphere.
     * @param radius The radius of the cap.
     */
    private constructor(center: S2Point, radius: S1Angle) : this(center, S1ChordAngle(minOf(radius, radians(M_PI))))

    /** the height of the cap, i.e. the distance from the center point to the cutoff plane. */
    val height: Double
        get() = 0.5 * radius.length2

    /**
     * Get the cap radius as an S1Angle. (Note that the cap angle is stored internally as an S1ChordAngle, so this
     * method requires a trigonometric operation and may yield a slightly different result than the value passed
     * to the (S2Point, S1Angle) constructor.)
     *
     * @return The radius as S1Angle.
     */
    fun radius(): S1Angle {
        return radius.toAngle()
    }

    /** the area of the cap. */
    val area: Double
        get() = 2 * M_PI * max(0.0, height)

    /**
     * The true centroid of the cap multiplied by its surface area (see S2Centroids for details on centroids).
     * The result lies on the ray from the origin through the cap's center, but it is not unit length. Note that if
     * you just want the "surface centroid", i.e. the normalized result, then it is much simpler just to call center.
     *
     * The reason for multiplying the result by the cap area is to make it easier to compute the centroid of more
     * complicated shapes.  The centroid of a union of disjoint regions can be computed simply by adding their
     * centroids. Caveat: for caps that contain a single point (i.e., zero radius), this method always returns the
     * origin (0, 0, 0).
     *
     * This is because shapes with no area don't affect the centroid of a union whose total area is positive.
     */
    val centroid: S2Point
        get() {
            // From symmetry, the centroid of the cap must be somewhere on the line from the origin to the center of
            // the cap on the surface of the sphere.
            // When a sphere is divided into slices of constant thickness by a set of parallel planes, all slices have
            // the same surface area. This implies that the radial component of the centroid is simply the midpoint of
            // the range of radial distances spanned by the cap. That is easily computed from the cap height.
            if (isEmpty) return S2Point()
            val r: Double = 1.0 - 0.5 * height
            return (r * area) * center
        }

    /**
     * Indicates if the cap is valid.
     * We allow negative radius (to represent empty caps) but not greater than 4.0.
     */
    val isValid: Boolean
        get() = isUnitLength(center) && radius.length2 <= 4.0

    /** Indicates if the cap is empty. true if the cap is empty, i.e. it contains no points.  */
    val isEmpty: Boolean
        get() = radius.isNegative()

    /** Indicates if the cap is full. true if the cap is full, i.e. it contains all points.  */
    val isFull: Boolean
        get() = radius.length2 == 4.0

    /**
     * The complement of the interior of the cap. A cap and its complement have the same boundary but do not share
     * any interior points. The complement operator is not a bijection, since the complement of a singleton cap
     * (containing a single point) is the same as the complement of an empty cap.
     */
    val complement: S2Cap
        get() =
        // The complement of a full cap is an empty cap, not a singleton.
                // Also make sure that the complement of an empty cap has height 2.
            when {
                isFull -> empty
                isEmpty -> full
                else -> S2Cap(-center, S1ChordAngle.fromLength2(4.0 - radius.length2))
            }

    /**
     * Indicates if this cap contains another one.
     *
     * @param other A cap.
     * @return true if and only if this cap contains the given other cap (in a set containment sense, e.g. every cap
     * contains the empty cap).
     */
    operator fun contains(other: S2Cap): Boolean =
            if (isFull || other.isEmpty) true
            else radius >= S1ChordAngle.between(center, other.center) + other.radius

    /**
     * Indicates if this cap intersects another one.
     *
     * @param other a cap.
     * @return true if and only if this cap intersects the given other cap, i.e. whether they have any points in common.
     */
    fun intersects(other: S2Cap): Boolean =
            if (isEmpty || other.isEmpty) false
            else radius + other.radius >= S1ChordAngle.between(center, other.center)

    /**
     * Indicates if the interior of this cap intersects an other cap.
     *
     * @param other A cap.
     * @return true if and only if the interior of this cap intersects the given other cap. (This relationship is not
     * symmetric, since only the interior of this cap is used.)
     */
    fun interiorIntersects(other: S2Cap): Boolean =
            // Make sure this cap has an interior and the other cap is non-empty.
            if (radius.length2 <= 0 || other.isEmpty) false
            else radius + other.radius > S1ChordAngle.between(center, other.center)

    override operator fun contains(p: S2Point): Boolean {
        requireArgument { isUnitLength(p) }
        return S1ChordAngle.between(center, p) <= radius
    }

    /**
     * Return true if and only if the given point is contained in the interior of
     * the region (i.e. the region excluding its boundary). 'p' should be a
     * unit-length vector.
     */
    fun interiorContains(p: S2Point): Boolean {
        requireArgument { isUnitLength(p) }
        return isFull || S1ChordAngle.between(center, p) < radius
    }

    /**
     * Increase the cap if necessary to include the given point. If the cap is empty the center is set to the given
     * point, but otherwise it is left unchanged. 'p' should be a unit-length vector.
     *
     * @param p A point
     * @return The new cap instance that includes this cap and the given point.
     */
    fun addPoint(p: S2Point): S2Cap {
        // Compute the squared chord length, then convert it into a height.
        requireArgument { isUnitLength(p) }
        return if (isEmpty) {
            S2Cap(p, S1ChordAngle.zero())
        } else {
            // After calling cap.addPoint(p), cap.contains(p) must be true.  However
            // we don't need to do anything special to achieve this because Contains()
            // does exactly the same distance calculation that we do here.
            S2Cap(center, maxOf(radius, S1ChordAngle.between(center, p)))
        }
    }

    /**
     * Increase the cap if necessary to include "other". If the current cap is empty it is set to the given other cap.
     *
     * @param other A cap.
     * @return The new cap instance that includes this cap and the other one.
     */
    fun addCap(other: S2Cap): S2Cap {
        return if (isEmpty) {
            other
        } else if (!other.isEmpty) {
            // We round up the distance to ensure that the cap is actually contained.
            // TODO(ericv): Do some error analysis in order to guarantee this.
            val dist = S1ChordAngle.between(center, other.center) + other.radius
            val newRadius = maxOf(radius, dist.plusError(DoubleType.epsilon * dist.length2))
            return S2Cap(center, newRadius)
        } else this
    }

    /**
     * Build an instance of this cap expanded by the given distance.
     *
     * @param distance The expansion amount.
     * @return a cap that contains all points within a given distance of this cap.  Note that any expansion of the
     * empty cap is still empty.
     */
    fun expanded(distance: S1Angle): S2Cap {
        requireGE(distance.radians, 0.0)
        return if (isEmpty) empty
        else S2Cap(center, radius + S1ChordAngle(distance))
    }

    /**
     * Build the union of this cap and another one.
     *
     * @param other A cap.
     * @return the smallest cap which encloses this cap and "other".
     */
    fun union(other: S2Cap): S2Cap {
        if (radius < other.radius) {
            return other.union(this)
        } else if (isFull || other.isEmpty) {
            return this
        } else {
            // This calculation would be more efficient using S1ChordAngles.
            val thisRadius = radius()
            val otherRadius = other.radius()
            val distance = S1Angle(center, other.center)
            if (thisRadius >= distance + otherRadius) {
                return this;
            } else {
                val resultRadius = 0.5 * (distance + thisRadius + otherRadius)
                val resultCenter = S2EdgeDistances.interpolateAtDistance(
                        0.5 * (distance - thisRadius + otherRadius),
                        center,
                        other.center)
                return S2Cap(resultCenter, resultRadius)
            }
        }
    }

    override fun equals(other: Any?): Boolean {
        if (other !is S2Cap) {
            return false
        }
        return (center == other.center && radius == other.radius || isEmpty && other.isEmpty || isFull && other.isFull)
    }

    override fun hashCode(): Int {
        var result = center.hashCode()
        result = 31 * result + radius.hashCode()
        return result
    }

    /**
     * Return true if the cap axis and height differ by at most "max_error" from
     * the given cap "other".
     */
    @JvmOverloads
    fun approxEquals(other: S2Cap, maxErrorAngle: S1Angle = S1Angle.radians(1e-14)): Boolean {
        val maxError = maxErrorAngle.radians
        val r2 = radius.length2
        val otherR2 = other.radius.length2

        return (S2PointUtil.approxEquals(center, other.center, maxErrorAngle) && abs(r2 - otherR2) <= maxError)
                || (isEmpty && otherR2 <= maxError)
                || (other.isEmpty && r2 <= maxError)
                || (isFull && otherR2 >= 2 - maxError)
                || (other.isFull && r2 >= 2 - maxError)
    }

    override fun toString(): String {
        return "[Center = $center, Radius = ${radius()}, Chord = $radius]"
    }

    // -- S2Region interface --//

    override fun clone(): S2Cap {
        return S2Cap(center, radius)
    }

    override val capBound: S2Cap
        get() = this

    override val rectBound: S2LatLngRect
        get() {
            if (isEmpty) return S2LatLngRect.empty()

            // Convert the center to a (lat,lng) pair, and compute the cap angle.
            val centerLl = S2LatLng.fromPoint(center)
            val capAngle = radius().radians

            var allLongitudes = false
            val lat = DoubleArray(2) { 0.0 }
            val lng = DoubleArray(2) { 0.0 }
            lng[0] = -M_PI
            lng[1] = M_PI

            // Check whether cap includes the south pole.
            lat[0] = centerLl.lat().radians - capAngle
            if (lat[0] <= -M_PI_2) {
                lat[0] = -M_PI_2;
                allLongitudes = true
            }
            // Check whether cap includes the north pole.
            lat[1] = centerLl.lat().radians + capAngle
            if (lat[1] >= M_PI_2) {
                lat[1] = M_PI_2
                allLongitudes = true;
            }
            if (!allLongitudes) {
                // Compute the range of longitudes covered by the cap.  We use the law
                // of sines for spherical triangles.  Consider the triangle ABC where
                // A is the north pole, B is the center of the cap, and C is the point
                // of tangency between the cap boundary and a line of longitude.  Then
                // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
                // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
                // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
                // minus the latitude).  This formula also works for negative latitudes.
                //
                // The formula for sin(a) follows from the relationship h = 1 - cos(a).
                val sinA = sin(radius)
                val sinC = cos(centerLl.lat().radians)
                if (sinA <= sinC) {
                    val angleA = asin(sinA / sinC)
                    lng[0] = IEEEremainder(centerLl.lng().radians - angleA, 2 * M_PI)
                    lng[1] = IEEEremainder(centerLl.lng().radians + angleA, 2 * M_PI)
                }
            }
            return S2LatLngRect(
                    R1Interval(lat[0], lat[1]),
                    S1Interval(lng[0], lng[1])
            )
        }

    // Computes a covering of the S2Cap.  In general the covering consists of at
    // most 4 cells except for very large caps, which may need up to 6 cells.
    // The output is not sorted.
    override fun getCellUnionBound(cellIds: MutableList<S2CellId>) {
        // TODO(ericv): The covering could be made quite a bit tighter by mapping
        // the cap to a rectangle in (i,j)-space and finding a covering for that.
        cellIds.clear()

        // Find the maximum level such that the cap contains at most one cell vertex
        // and such that S2CellId::AppendVertexNeighbors() can be called.
        val level = S2Coords.projection.kMinWidth.getLevelForMinValue(radius().radians) - 1// S2::kMinWidth.GetLevelForMinValue(GetRadius().radians()) - 1;

        // If level < 0, then more than three face cells are required.
        if (level < 0) {
            for (face in 0..5) {
                cellIds.add(S2CellId.fromFace(face))
            }
        } else {
            // The covering consists of the 4 cells at the given level that share the
            // cell vertex that is closest to the cap center.
            S2CellId.fromPoint(center).appendVertexNeighbors(level, cellIds)
        }
    }

    /**
     * Return true if the cap intersects 'cell', given that the cap vertices have
     * alrady been checked.
     */
    fun intersects(cell: S2Cell, vertices: Array<S2Point>): Boolean {
        // Return true if this cap intersects any point of 'cell' excluding its  vertices (which are assumed to already
        // have been checked).

        // If the cap is a hemisphere or larger, the cell and the complement of the
        // cap are both convex. Therefore since no vertex of the cell is contained,
        // no other interior point of the cell is contained either.
        if (radius >= S1ChordAngle.right()) {
            return false
        }

        // We need to check for empty caps due to the axis check just below.
        if (isEmpty) {
            return false
        }

        // Optimization: return true if the cell contains the cap axis. (This
        // allows half of the edge checks below to be skipped.)
        if (cell.contains(center)) {
            return true
        }

        // At this point we know that the cell does not contain the cap axis,
        // and the cap does not contain any cell vertex. The only way that they
        // can intersect is if the cap intersects the interior of some edge.
        val sin2Angle = radius.sin2() // sin^2(capAngle)
        for (k in 0..3) {
            val edge = cell.getEdgeRaw(k)
            val dot = center.dotProd(edge)
            if (dot > 0) {
                // The axis is in the interior half-space defined by the edge. We don't
                // need to consider these edges, since if the cap intersects this edge
                // then it also intersects the edge on the opposite side of the cell
                // (because we know the axis is not contained with the cell).
                continue
            }
            // The norm2() factor is necessary because "edge" is not normalized.
            if (dot * dot > sin2Angle * edge.norm2()) {
                return false // Entire cap is on the exterior side of this edge.
            }
            // Otherwise, the great circle containing this edge intersects
            // the interior of the cap. We just need to check whether the point
            // of closest approach occurs between the two edge endpoints.
            val dir = edge.crossProd(center)
            if (dir.dotProd(vertices[k]) < 0 && dir.dotProd(vertices[k + 1 and 3]) > 0) {
                return true
            }
        }
        return false
    }

    override fun contains(cell: S2Cell): Boolean {
        // If the cap does not contain all cell vertices, return false.
        // We check the vertices before taking the Complement() because we can't
        // accurately represent the complement of a very small cap (a height
        // of 2-epsilon is rounded off to 2).
        val vertices = Array(4) { k ->
            val vertex = cell.getVertex(k)
            if (!contains(vertex))
                return false
            vertex
        }

        // Otherwise, return true if the complement of the cap does not intersect
        // the cell. (This test is slightly conservative, because technically we
        // want Complement().InteriorIntersects() here.)
        return !complement.intersects(cell, vertices)
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        // If the cap contains any cell vertex, return true.
        val vertices = Array(4) { k ->
            val vertex = cell.getVertex(k)!!
            if (contains(vertex)) return true
            vertex
        }

        return intersects(cell, vertices)
    }

    companion object {

        /**
         * Convenience function that creates a cap containing a single point. This method is more efficient that the
         * S2Cap(center, radius) constructor.
         *
         * @param center A point.
         * @return A cap instance that only include the given point.
         */
        @JvmStatic
        fun fromPoint(center: S2Point): S2Cap = S2Cap(center, S1ChordAngle.zero())

        /**
         * Create a cap with the given center and height. A negative height yields an empty cap; a height of 2 or more
         * yields a full cap.  "center" should be unit length.
         *
         * @param center The center of the cap.
         * @param height The height of the cap.
         * @return The new cap instance.
         */
        @JvmStatic
        fun fromCenterHeight(center: S2Point, height: Double): S2Cap {
            requireArgument { isUnitLength(center) }
            return S2Cap(center, S1ChordAngle.fromLength2(2.0 * height))
        }

        /**
         * Create a cap with the given center and surface area.  Note that the area can also be interpreted as the
         * solid angle subtended by the cap (because the sphere has unit radius).  A negative area yields an empty cap;
         * an area of 4*Pi or more yields a full cap.  "center" should be unit length.
         *
         * @param center The center of the cap.
         * @param area The area of the cap.
         * @return The new cap instance.
         */
        @JvmStatic
        fun fromCenterArea(center: S2Point, area: Double): S2Cap {
            requireArgument { isUnitLength(center) }
            return S2Cap(center, S1ChordAngle.fromLength2(area / M_PI))
        }

        /** Empty cap, i.e. a cap that contains no points.  */
        @JvmField
        val empty: S2Cap = S2Cap()

        /** Full cap, i.e. a cap that contains all points. */
        @JvmField
        val full: S2Cap = S2Cap(S2Point(1, 0, 0), S1ChordAngle.straight())

        /**
         * Create a cap given its center and the cap opening angle, i.e. maximum angle between the center and a point on
         * the cap. 'center' should be a unit-length vector, and 'angle' should be between 0 and 180 degrees.
         *
         * @param center The center point of the cap.
         * @param angle the maximum angle of the cap.
         * @return The new cap instance.
         */
        @JvmStatic
        fun fromCenterAngle(center: S2Point, angle: S1Angle): S2Cap {
            requireArgument { isUnitLength(center) }
            return S2Cap(center, angle)
        }

        /**
         * Multiply a positive number by this constant to ensure that the result of a
         * floating point operation is at least as large as the true
         * infinite-precision result.
         */
        //private const val ROUND_UP = 1.0 + 1.0 / (1L shl 52)

    }
}
