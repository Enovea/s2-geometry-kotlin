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
import dilivia.PreConditions.requireLE
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.R1Interval
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.cos
import dilivia.s2.S1Angle.Companion.sin
import dilivia.s2.S1ChordAngle
import dilivia.s2.S1Interval
import dilivia.s2.S2LatLng
import dilivia.s2.S2LatLng.Companion.times
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2PointUtil.robustCrossProd
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.edge.S2EdgeDistances
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.acos
import org.apache.commons.math3.util.FastMath.asin
import org.apache.commons.math3.util.FastMath.atan2
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sin
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * An S2LatLngRect represents a closed latitude-longitude rectangle. It is capable of representing the empty and full
 * rectangles as well as single points.  Note that the latitude-longitude space is considered to have a *cylindrical*
 * topology rather than a spherical one, i.e. the poles have multiple lat/lng representations.
 * An S2LatLngRect may be defined so that includes some representations of a pole but not others.
 * Use the polarClosure() method if you want to expand a rectangle so that it contains all possible representations of
 * any contained poles.
 *
 * Because S2LatLngRect uses S1Interval to store the longitude range, longitudes of -180 degrees are treated specially.
 * Except for empty and full longitude spans, -180 degree longitudes will turn into +180 degrees. This sign flip causes
 * lngLo() to be greater than lngHi(), indicating that the rectangle will wrap around through -180 instead of through
 * +179. Thus the math is consistent within the library, but the sign flip can be surprising, especially when working
 * with map projections where -180 and +180 are at opposite ends of the flattened map.  See the comments on S1Interval
 * for more details.
 *
 * This class is a port of the S2LatLngRect class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @property lat The latitude interval.
 * @property lng The longitude interval.
 * @constructor Construct a rectangle from latitude and longitude intervals. The lat/lng ranges must either be both
 * empty or both non-empty.
 * @param lat The latitude interval.
 * @param lng The longitude interval.
 */
class S2LatLngRect(val lat: R1Interval, val lng: S1Interval) : S2Region {

    init {
        checkState { isValid }
    }

    /**
     * Construct a rectangle from minimum and maximum latitudes and longitudes. If lo.lng() > hi.lng(), the rectangle
     * spans the 180 degree longitude line. Both points must be normalized, with lo.lat() <= hi.lat().
     * The rectangle contains all the points p such that 'lo' <= p <= 'hi', where '<=' is defined in the obvious way.
     *
     * @param lo The lower-left corner.
     * @param hi The upper-right corner.
     */
    constructor(lo: S2LatLng, hi: S2LatLng) : this(
            lat = R1Interval(lo.lat().radians, hi.lat().radians),
            lng = S1Interval(lo.lng().radians, hi.lng().radians)
    )

    /** The default constructor creates an empty S2LatLngRect. */
    constructor() : this(
            lat = R1Interval.empty(),
            lng = S1Interval.empty()
    )

    fun init(other: S2LatLngRect) {
        lat.set(other.lat)
        lng.set(other.lng)
    }

    fun setEmpty() {
        lat.setEmpty()
        lng.setEmpty()
    }

    fun setFull() {
        lat.set(fullLat())
        lng.set(fullLng())
    }

    // Accessor methods.
    fun latLo(): S1Angle = S1Angle.radians(lat.lo)
    fun latHi(): S1Angle = S1Angle.radians(lat.hi)
    fun lngLo(): S1Angle = S1Angle.radians(lng.lo)
    fun lngHi(): S1Angle = S1Angle.radians(lng.hi)
    fun lo(): S2LatLng = S2LatLng.fromLatLng(latLo(), lngLo())
    fun hi(): S2LatLng = S2LatLng.fromLatLng(latHi(), lngHi())

    /**
     * Return true if the rectangle is valid, which essentially just means that the latitude bounds do not exceed Pi/2
     * in absolute value and the longitude bounds do not exceed Pi in absolute value. Also, if either the latitude or
     * longitude bound is empty then both must be.
     */
    // The lat/lng ranges must either be both empty or both non-empty.
    val isValid: Boolean
        get() = abs(lat.lo) <= M_PI_2 && abs(lat.hi) <= M_PI_2 && lng.isValid && lat.isEmpty == lng.isEmpty

    /**
     * Return true if the rectangle is empty, i.e. it contains no points at all.
     */
    val isEmpty: Boolean
        get() = lat.isEmpty

    // Return true if the rectangle is full, i.e. it contains all points.
    val isFull: Boolean
        get() = lat == fullLat() && lng.isFull

    // Returns true if the rectangle is a point, i.e. lo() == hi()
    val isPoint: Boolean
        get() = lat.lo == lat.hi && lng.lo == lng.hi

    /**
     * Return true if lng_.getLo() > lng_.getHi(), i.e. the rectangle crosses the 180
     * degree latitude line.
     */
    val isInverted: Boolean
        get() = lng.isInverted

    /** Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.  */
    fun getVertex(k: Int): S2LatLng {
        // Twiddle bits to return the points in CCW order (lower left, lower right,
        // upper right, upper left).
        val i = (k shr 1) and 1
        return S2LatLng.fromRadians(lat[i], lng[i xor (k and 1)])
    }

    fun getLowerLeftVertex() = getVertex(0)
    fun getLowerRightVertex() = getVertex(1)
    fun getUpperRightVertex() = getVertex(2)
    fun getUpperLeftVertex() = getVertex(3)

    /**
     * Return the center of the rectangle in latitude-longitude space (in general
     * this is not the center of the region on the sphere).
     */
    val center: S2LatLng
        get() = S2LatLng.fromRadians(lat.center, lng.center)

    /**
     * Return the width and height of this rectangle in latitude-longitude space.
     * Empty rectangles have a negative width and height.
     */
    val size: S2LatLng
        get() = S2LatLng.fromRadians(lat.length, lng.length)

    /** Return the surface area of this rectangle on the unit sphere.  */
    val area: Double
        get() {
            return if (isEmpty) 0.0
            else lng.length * (sin(latHi().radians) - sin(latLo().radians))
            // This is the size difference of the two spherical caps, multiplied by the longitude ratio.
        }

    // Returns the true centroid of the rectangle multiplied by its surface area
    // (see S2Centroids for details on centroids).  The result is not unit
    // length, so you may want to normalize it.  Note that in general the
    // centroid is *not* at the center of the rectangle, and in fact it may not
    // even be contained by the rectangle.  (It is the "center of mass" of the
    // rectangle viewed as subset of the unit sphere, i.e. it is the point in
    // space about which this curved shape would rotate.)
    //
    // The reason for multiplying the result by the rectangle area is to make it
    // easier to compute the centroid of more complicated shapes.  The centroid
    // of a union of disjoint regions can be computed simply by adding their
    // GetCentroid() results.
    val centroid: S2Point
        get() {
            // When a sphere is divided into slices of constant thickness by a set of
            // parallel planes, all slices have the same surface area.  This implies
            // that the z-component of the centroid is simply the midpoint of the
            // z-interval spanned by the S2LatLngRect.
            //
            // Similarly, it is easy to see that the (x,y) of the centroid lies in the
            // plane through the midpoint of the rectangle's longitude interval.  We
            // only need to determine the distance "d" of this point from the z-axis.
            //
            // Let's restrict our attention to a particular z-value.  In this z-plane,
            // the S2LatLngRect is a circular arc.  The centroid of this arc lies on a
            // radial line through the midpoint of the arc, and at a distance from the
            // z-axis of
            //
            //     r * (sin(alpha) / alpha)
            //
            // where r = sqrt(1-z^2) is the radius of the arc, and "alpha" is half of
            // the arc length (i.e., the arc covers longitudes [-alpha, alpha]).
            //
            // To find the centroid distance from the z-axis for the entire rectangle,
            // we just need to integrate over the z-interval.  This gives
            //
            //    d = Integrate[sqrt(1-z^2)*sin(alpha)/alpha, z1..z2] / (z2 - z1)
            //
            // where [z1, z2] is the range of z-values covered by the rectangle.  This
            // simplifies to
            //
            //    d = sin(alpha)/(2*alpha*(z2-z1))*(z2*r2 - z1*r1 + theta2 - theta1)
            //
            // where [theta1, theta2] is the latitude interval, z1=sin(theta1),
            // z2=sin(theta2), r1=cos(theta1), and r2=cos(theta2).
            //
            // Finally, we want to return not the centroid itself, but the centroid
            // scaled by the area of the rectangle.  The area of the rectangle is
            //
            //    A = 2 * alpha * (z2 - z1)
            //
            // which fortunately appears in the denominator of "d".

            if (isEmpty) return S2Point()
            val z1 = sin(latLo())
            val z2 = sin(latHi())
            val r1 = cos(latLo())
            val r2 = cos(latHi())
            val alpha = 0.5 * lng.length
            val r = sin(alpha) * (r2 * z2 - r1 * z1 + lat.length)
            val lng = lng.center
            val z = alpha * (z2 + z1) * (z2 - z1)  // scaled by the area
            return S2Point(r * cos(lng), r * sin(lng), z)
        }

    /**
     * More efficient version of Contains() that accepts a S2LatLng rather than an S2Point. The argument must be
     * normalized.
     */
    operator fun contains(ll: S2LatLng): Boolean {
        requireArgument { ll.isValid }
        return lat.contains(ll.lat().radians) && lng.contains(ll.lng().radians)
    }

    /**
     * Return true if and only if the given point is contained in the interior of
     * the region (i.e. the region excluding its boundary). The point 'p' does not
     * need to be normalized.
     */
    fun interiorContains(p: S2Point): Boolean {
        return interiorContains(S2LatLng.fromPoint(p))
    }

    /**
     * More efficient version of interiorContains() that accepts a S2LatLng rather than an S2Point.
     * The argument must be normalized.
     */
    fun interiorContains(ll: S2LatLng): Boolean {
        requireArgument { ll.isValid }
        return lat.interiorContains(ll.lat().radians) && lng.interiorContains(ll.lng().radians)
    }

    /**
     * Return true if and only if the rectangle contains the given other rectangle.
     */
    operator fun contains(other: S2LatLngRect): Boolean {
        return lat.contains(other.lat) && lng.contains(other.lng)
    }

    /**
     * Return true if and only if the interior of this rectangle contains all
     * points of the given other rectangle (including its boundary).
     */
    fun interiorContains(other: S2LatLngRect): Boolean {
        return lat.interiorContains(other.lat) && lng.interiorContains(other.lng)
    }

    /** Return true if this rectangle and the given other rectangle have any
     * points in common.  */
    fun intersects(other: S2LatLngRect): Boolean {
        return lat.intersects(other.lat) && lng.intersects(other.lng)
    }

    /**
     * Return true if and only if the interior of this rectangle intersects any
     * point (including the boundary) of the given other rectangle.
     */
    fun interiorIntersects(other: S2LatLngRect): Boolean {
        return lat.interiorIntersects(other.lat) && lng.interiorIntersects(other.lng)
    }

    // Returns true if the boundary of this rectangle intersects the given
    // geodesic edge (v0, v1).
    fun boundaryIntersects(v0: S2Point, v1: S2Point): Boolean {
        if (isEmpty) return false
        if (!lng.isFull) {
            if (intersectsLngEdge(v0, v1, lat, lng.lo)) return true
            if (intersectsLngEdge(v0, v1, lat, lng.hi)) return true
        }
        if (lat.lo != -M_PI_2 && intersectsLatEdge(v0, v1, lat.lo, lng)) {
            return true
        }
        if (lat.hi != M_PI_2 && intersectsLatEdge(v0, v1, lat.hi, lng)) {
            return true
        }
        return false
    }

    fun addPoint(p: S2Point): S2LatLngRect {
        return addPoint(S2LatLng.fromPoint(p))
    }

    // Increase the size of the bounding rectangle to include the given point.
    // The rectangle is expanded by the minimum amount possible.
    fun addPoint(ll: S2LatLng): S2LatLngRect {
        requireArgument { ll.isValid }
        val newLat = lat.clone().addPoint(ll.lat().radians)
        val newLng = lng.clone().addPoint(ll.lng().radians)
        return S2LatLngRect(newLat, newLng)
    }

    // Returns a rectangle that has been expanded by margin.lat() on each side in
    // the latitude direction, and by margin.lng() on each side in the longitude
    // direction.  If either margin is negative, then shrinks the rectangle on
    // the corresponding sides instead.  The resulting rectangle may be empty.
    //
    // As noted above, the latitude-longitude space has the topology of a
    // cylinder.  Longitudes "wrap around" at +/-180 degrees, while latitudes
    // are clamped to range [-90, 90].  This means that any expansion (positive
    // or negative) of the full longitude range remains full (since the
    // "rectangle" is actually a continuous band around the cylinder), while
    // expansion of the full latitude range remains full only if the margin is
    // positive.
    //
    // If either the latitude or longitude interval becomes empty after
    // expansion by a negative margin, the result is empty.
    //
    // Note that if an expanded rectangle contains a pole, it may not contain
    // all possible lat/lng representations of that pole (see header above).
    // Use the PolarClosure() method if you do not want this behavior.
    //
    // If you are trying to grow a rectangle by a certain *distance* on the
    // sphere (e.g. 5km), use the ExpandedByDistance() method instead.
    fun expanded(margin: S2LatLng): S2LatLngRect {
        // assert (margin.lat().getRadians() >= 0 && margin.lng().getRadians() >= 0);
        val latitudes = lat.expanded(margin.lat().radians)
        val longitudes = lng.expanded(margin.lng().radians)
        return if (latitudes.isEmpty || longitudes.isEmpty) empty()
        else S2LatLngRect(latitudes.intersection(fullLat()), longitudes)
    }

    // If the rectangle does not include either pole, returns it unmodified.
    // Otherwise expands the longitude range to Full() so that the rectangle
    // contains all possible representations of the contained pole(s).
    fun polarClosure(): S2LatLngRect {
        if (lat.lo == -M_PI_2 || lat.hi == M_PI_2) {
            return S2LatLngRect(lat, S1Interval.full())
        }
        return this
    }

    /**
     * Return the smallest rectangle containing the union of this rectangle and
     * the given rectangle.
     */
    fun union(other: S2LatLngRect): S2LatLngRect {
        return S2LatLngRect(lat.union(other.lat), lng.union(other.lng))
    }

    /**
     * Return the smallest rectangle containing the intersection of this rectangle
     * and the given rectangle. Note that the region of intersection may consist
     * of two disjoint rectangles, in which case a single rectangle spanning both
     * of them is returned.
     */
    fun intersection(other: S2LatLngRect): S2LatLngRect {
        val intersectLat = lat.intersection(other.lat)
        val intersectLng = lng.intersection(other.lng)
        return if (intersectLat.isEmpty || intersectLng.isEmpty) {
            // The lat/lng ranges must either be both empty or both non-empty.
            empty()
        } else S2LatLngRect(intersectLat, intersectLng)
    }

    // Expands this rectangle so that it contains all points within the given
    // distance of the boundary, and return the smallest such rectangle.  If the
    // distance is negative, then instead shrinks this rectangle so that it
    // excludes all points within the given absolute distance of the boundary,
    // and returns the largest such rectangle.
    //
    // Unlike Expanded(), this method treats the rectangle as a set of points on
    // the sphere, and measures distances on the sphere.  For example, you can
    // use this method to find a rectangle that contains all points within 5km
    // of a given rectangle.  Because this method uses the topology of the
    // sphere, note the following:
    //
    //  - The full and empty rectangles have no boundary on the sphere.  Any
    //    expansion (positive or negative) of these rectangles leaves them
    //    unchanged.
    //
    //  - Any rectangle that covers the full longitude range does not have an
    //    east or west boundary, therefore no expansion (positive or negative)
    //    will occur in that direction.
    //
    //  - Any rectangle that covers the full longitude range and also includes
    //    a pole will not be expanded or contracted at that pole, because it
    //    does not have a boundary there.
    //
    //  - If a rectangle is within the given distance of a pole, the result will
    //    include the full longitude range (because all longitudes are present
    //    at the poles).
    //
    // Expansion and contraction are defined such that they are inverses whenver
    // possible, i.e.
    //
    //   rect.ExpandedByDistance(x).ExpandedByDistance(-x) == rect
    //
    // (approximately), so long as the first operation does not cause a
    // rectangle boundary to disappear (i.e., the longitude range newly becomes
    // full or empty, or the latitude range expands to include a pole).
    fun expandedByDistance(distance: S1Angle): S2LatLngRect {
        if (distance >= S1Angle.zero()) {
            // The most straightforward approach is to build a cap centered on each
            // vertex and take the union of all the bounding rectangles (including the
            // original rectangle; this is necessary for very large rectangles).

            // TODO(ericv): Update this code to use an algorithm like the one below.
            val radius = S1ChordAngle(distance)
            var r = this
            for (k in 0..3) {
                r = r.union(S2Cap(getVertex(k).toPoint(), radius).rectBound)
            }
            return r
        } else {
            // Shrink the latitude interval unless the latitude interval contains a pole
            // and the longitude interval is full, in which case the rectangle has no
            // boundary at that pole.
            val latResult = R1Interval(
                    if (lat.lo <= fullLat().lo && lng.isFull) fullLat().lo else lat.lo - distance.radians,
                    if (lat.hi >= fullLat().hi && lng.isFull) fullLat().hi else lat.hi + distance.radians
            )
            if (latResult.isEmpty) {
                return empty()
            }

            // Maximum absolute value of a latitude in lat_result. At this latitude,
            // the cap occupies the largest longitude interval.
            val maxAbsLat = max(-latResult.lo, latResult.hi)

            // Compute the largest longitude interval that the cap occupies. We use the
            // law of sines for spherical triangles. For the details, see the comment in
            // S2Cap::GetRectBound().
            //
            // When sin_a >= sin_c, the cap covers all the latitude.
            val sinA = sin(-distance.radians)
            val sinC = cos(maxAbsLat)
            val maxLngMargin = if (sinA < sinC) asin(sinA / sinC) else M_PI_2

            val lngResult = lng.expanded(-maxLngMargin)
            if (lngResult.isEmpty) {
                return empty()
            }
            return S2LatLngRect(latResult, lngResult)
        }
    }


    /**
     * Returns true if this rectangle intersects the given cell. (This is an exact
     * test and may be fairly expensive, see also MayIntersect below.)
     */
    fun intersects(cell: S2Cell): Boolean {
        // First we eliminate the cases where one region completely contains the
        // other. Once these are disposed of, then the regions will intersect
        // if and only if their boundaries intersect.
        if (isEmpty) {
            return false
        }
        if (contains(cell.getCenter())) {
            return true
        }
        if (cell.contains(center.toPoint())) {
            return true
        }

        // Quick rejection test (not required for correctness).
        if (!intersects(cell.rectBound)) {
            return false
        }

        // Now check whether the boundaries intersect. Unfortunately, a
        // latitude-longitude rectangle does not have straight edges -- two edges
        // are curved, and at least one of them is concave.

        // Precompute the cell vertices as points and latitude-longitudes.
        val cellV = arrayOfNulls<S2Point>(4)
        val cellLl = arrayOfNulls<S2LatLng>(4)
        for (i in 0..3) {
            cellV[i] = cell.getVertex(i) // Must be normalized.
            cellLl[i] = S2LatLng.fromPoint(cellV[i]!!)
            if (contains(cellLl[i]!!)) {
                return true // Quick acceptance test.
            }
        }
        for (i in 0..3) {
            val edgeLng = S1Interval.fromPointPair(
                    cellLl[i]!!.lng().radians, cellLl[i + 1 and 3]!!.lng().radians)
            if (!lng.intersects(edgeLng)) {
                continue
            }
            val a = cellV[i]!!
            val b = cellV[i + 1 and 3]!!
            if (edgeLng.contains(lng.lo)) {
                if (intersectsLngEdge(a, b, lat, lng.lo)) {
                    return true
                }
            }
            if (edgeLng.contains(lng.hi)) {
                if (intersectsLngEdge(a, b, lat, lng.hi)) {
                    return true
                }
            }
            if (intersectsLatEdge(a, b, lat.lo, lng)) {
                return true
            }
            if (intersectsLatEdge(a, b, lat.hi, lng)) {
                return true
            }
        }
        return false
    }


    /**
     * Return the minimum distance (measured along the surface of the sphere) to
     * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
     */
    fun getDistance(other: S2LatLngRect): S1Angle {
        check(!isEmpty)
        require(!other.isEmpty)
        val a = this
        val b = other

        // First, handle the trivial cases where the longitude intervals overlap.
        if (a.lng.intersects(b.lng)) {
            if (a.lat.intersects(b.lat)) {
                return S1Angle.radians(0.0) // Intersection between a and b.
            }

            // We found an overlap in the longitude interval, but not in the latitude
            // interval. This means the shortest path travels along some line of
            // longitude connecting the high-latitude of the lower rect with the
            // low-latitude of the higher rect.
            val lo: S1Angle
            val hi: S1Angle
            if (a.lat.lo > other.lat.hi) {
                lo = b.latHi()
                hi = a.latLo()
            } else {
                lo = a.latHi()
                hi = b.latLo()
            }
            return hi - lo
        }

        // The longitude intervals don't overlap. In this case, the closest points
        // occur somewhere on the pair of longitudinal edges which are nearest in
        // longitude-space.
        val aLng: S1Angle
        val bLng: S1Angle
        val loHi = S1Interval.fromPointPair(a.lng.lo, b.lng.hi)
        val hiLo = S1Interval.fromPointPair(a.lng.hi, b.lng.lo)
        if (loHi.length < hiLo.length) {
            aLng = a.lngLo()
            bLng = b.lngHi()
        } else {
            aLng = a.lngHi()
            bLng = b.lngLo()
        }

        // The shortest distance between the two longitudinal segments will include
        // at least one segment endpoint. We could probably narrow this down further
        // to a single point-edge distance by comparing the relative latitudes of the
        // endpoints, but for the sake of clarity, we'll do all four point-edge
        // distance tests.
        val aLo = S2LatLng.fromLatLng(a.latLo(), aLng).toPoint()
        val aHi = S2LatLng.fromLatLng(a.latHi(), aLng).toPoint()
        val bLo = S2LatLng.fromLatLng(b.latLo(), bLng).toPoint()
        val bHi = S2LatLng.fromLatLng(b.latHi(), bLng).toPoint()
        return minOf(
            S2EdgeDistances.getDistance(aLo, bLo, bHi),
                minOf(S2EdgeDistances.getDistance(aHi, bLo, bHi),
                        minOf(
                                S2EdgeDistances.getDistance(bLo, aLo, aHi),
                                S2EdgeDistances.getDistance(bHi, aLo, aHi)
                        )
                )
        )
    }

    /**
     * Return the minimum distance (measured along the surface of the sphere)
     * from a given point to the rectangle (both its boundary and its interior).
     * The latLng must be valid.
     */
    fun getDistance(p: S2LatLng): S1Angle {
        check(!isEmpty)
        require(p.isValid)
        // The algorithm here is the same as in getDistance(S2LagLngRect), only
        // with simplified calculations.
        val a = this
        if (a.lng.contains(p.lng().radians)) {
            return S1Angle.radians(max(0.0, max(p.lat().radians - a.lat.hi, a.lat.lo - p.lat().radians)))
        }
        val interval = S1Interval(a.lng.hi, a.lng.complementCenter)
        var aLng = a.lng.lo
        if (interval.contains(p.lng().radians)) {
            aLng = a.lng.hi
        }
        val lo = S2LatLng.fromRadians(a.lat.lo, aLng).toPoint()
        val hi = S2LatLng.fromRadians(a.lat.hi, aLng).toPoint()
        return S2EdgeDistances.getDistance(p.toPoint(), lo, hi)
    }

    // Returns the (directed or undirected) Hausdorff distance (measured along the
    // surface of the sphere) to the given S2LatLngRect. The directed Hausdorff
    // distance from rectangle A to rectangle B is given by
    //     h(A, B) = max_{p in A} min_{q in B} d(p, q).
    // The Hausdorff distance between rectangle A and rectangle B is given by
    //     H(A, B) = max{h(A, B), h(B, A)}.
    fun getHausdorffDistance(other: S2LatLngRect): S1Angle {
        return maxOf(getDirectedHausdorffDistance(other), other.getDirectedHausdorffDistance(this))
    }

    fun getDirectedHausdorffDistance(other: S2LatLngRect): S1Angle {
        if (isEmpty) {
            return S1Angle.radians(0)
        }
        if (other.isEmpty) {
            return S1Angle.radians(M_PI)  // maximum possible distance on S2
        }

        val lng_distance = lng.getDirectedHausdorffDistance(other.lng)
        assert(lng_distance >= 0)
        return getDirectedHausdorffDistance(lng_distance, lat, other.lat)
    }

    /** Return true if two rectangles contains the same set of points.  */
    override fun equals(other: Any?): Boolean {
        if (other !is S2LatLngRect) {
            return false
        }
        return lat == other.lat && lng == other.lng
    }

    override fun hashCode(): Int {
        var value = 17
        value = 37 * value + lat.hashCode()
        return 37 * value + lng.hashCode()
    }

    /**
     * Return true if the latitude and longitude intervals of the two rectangles
     * are the same up to the given tolerance (see r1interval.h and s1interval.h
     * for details).
     */
    @JvmOverloads
    fun approxEquals(other: S2LatLngRect, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        return (lat.approxEquals(other.lat, maxError.radians) && lng.approxEquals(other.lng, maxError.radians))
    }

    fun approxEquals(other: S2LatLngRect, maxError: S2LatLng): Boolean {
        return (lat.approxEquals(other.lat, maxError.lat().radians) && lng.approxEquals(other.lng, maxError.lng().radians))
    }

    /**
     * Return a rectangle that contains the convolution of this rectangle with a
     * cap of the given angle. This expands the rectangle by a fixed distance (as
     * opposed to growing the rectangle in latitude-longitude space). The returned
     * rectangle includes all points whose minimum distance to the original
     * rectangle is at most the given angle.
     */
    fun convolveWithCap(angle: S1Angle): S2LatLngRect {
        // The most straightforward approach is to build a cap centered on each
        // vertex and take the union of all the bounding rectangles (including the
        // original rectangle; this is necessary for very large rectangles).

        // Optimization: convert the angle to a height exactly once.
        val cap = S2Cap.fromCenterAngle(S2Point(1, 0, 0), angle)
        var r = this
        for (k in 0..3) {
            val vertexCap = S2Cap.fromCenterHeight(getVertex(k).toPoint(), cap.height)
            r = r.union(vertexCap.rectBound)
        }
        return r
    }

    override fun toString(): String {
        return "[Lo=" + lo() + ", Hi=" + hi() + "]"
    }

    // -- S2Region interface --//

    public override fun clone(): S2LatLngRect {
        return S2LatLngRect(lo(), hi())
    }

    // For bounding rectangles that span 180 degrees or less in longitude, the
    // maximum cap size is achieved at one of the rectangle vertices. For
    // rectangles that are larger than 180 degrees, we punt and always return a
    // bounding cap centered at one of the two poles.
    // We consider two possible bounding caps, one whose axis passes
    // through the center of the lat-long rectangle and one whose axis
    // is the north or south pole. We return the smaller of the two caps.
    override val capBound: S2Cap
        get() {
            // We consider two possible bounding caps, one whose axis passes
            // through the center of the lat-long rectangle and one whose axis
            // is the north or south pole. We return the smaller of the two caps.
            if (isEmpty) {
                return S2Cap.empty
            }
            val poleZ: Double
            val poleAngle: Double
            if (lat.lo + lat.hi < 0) {
                // South pole axis yields smaller cap.
                poleZ = -1.0
                poleAngle = M_PI_2 + lat.hi
            } else {
                poleZ = 1.0
                poleAngle = M_PI_2 - lat.lo
            }
            val poleCap = S2Cap.fromCenterAngle(S2Point(0.0, 0.0, poleZ), S1Angle.radians(poleAngle))

            // For bounding rectangles that span 180 degrees or less in longitude, the
            // maximum cap size is achieved at one of the rectangle vertices. For
            // rectangles that are larger than 180 degrees, we punt and always return a
            // bounding cap centered at one of the two poles.
            val lngSpan = lng.hi - lng.lo
            if (IEEEremainder(lngSpan, 2 * M_PI) >= 0 && lngSpan < 2 * M_PI_2) {
                var midCap = S2Cap.fromCenterAngle(center.toPoint(), S1Angle.radians(0.0))
                for (k in 0..3) {
                    midCap = midCap.addPoint(getVertex(k).toPoint())
                }
                if (midCap.height < poleCap.height) {
                    return midCap
                }
            }
            return poleCap
        }

    override val rectBound: S2LatLngRect
        get() = this

    override fun contains(cell: S2Cell): Boolean {
        // A latitude-longitude rectangle contains a cell if and only if it contains
        // the cell's bounding rectangle.  This test is exact from a mathematical
        // point of view, assuming that the bounds returned by S2Cell::GetRectBound()
        // are tight.  However, note that there can be a loss of precision when
        // converting between representations -- for example, if an S2Cell is
        // converted to a polygon, the polygon's bounding rectangle may not contain
        // the cell's bounding rectangle.  This has some slightly unexpected side
        // effects; for instance, if one creates an S2Polygon from an S2Cell, the
        // polygon will contain the cell, but the polygon's bounding box will not.
        return contains(cell.rectBound)
    }

    /**
     * This test is cheap but is NOT exact. Use Intersects() if you want a more
     * accurate and more expensive test. Note that when this method is used by an
     * S2RegionCoverer, the accuracy isn't all that important since if a cell may
     * intersect the region then it is subdivided, and the accuracy of this method
     * goes up as the cells get smaller.
     */
    override fun mayIntersect(cell: S2Cell): Boolean {
        // This test is cheap but is NOT exact (see s2latlngrect.h).
        return intersects(cell.rectBound)
    }

    /** The point 'p' does not need to be normalized.  */
    override operator fun contains(p: S2Point): Boolean {
        return contains(S2LatLng.fromPoint(p))
    }

    companion object {

        /**
         * Construct a rectangle from a center point (in lat-lng space) and size in each dimension. If size.lng() is
         * greater than 360 degrees it is clamped, and latitudes greater than +/- 90 degrees are also clamped.
         * So for example, fromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
         *
         * @param center The center of the rectangle.
         * @param size The size of the rectangle.
         * @return The rectangle instance.
         */
        @JvmStatic
        fun fromCenterSize(center: S2LatLng, size: S2LatLng): S2LatLngRect = fromPoint(center).expanded(0.5 * size)

        /**
         * Convenience method to construct a rectangle containing a single (normalized) point.
         *
         * @param p A point.
         * @return The rectangle instance.
         */
        @JvmStatic
        fun fromPoint(p: S2LatLng): S2LatLngRect {
            requireArgument { p.isValid }
            return S2LatLngRect(p, p)
        }

        /**
         * Convenience method to construct the minimal bounding rectangle containing the two given points. This is
         * equivalent to starting with an empty rectangle and calling addPoint() twice. Note that it is different than
         * the S2LatLngRect(lo, hi) constructor, where the first point is always used as the lower-left corner of the
         * resulting rectangle.
         *
         * @param p1 The first point
         * @param p2 The second point
         * @return The rectangle instance.
         */
        @JvmStatic
        fun fromPointPair(p1: S2LatLng, p2: S2LatLng): S2LatLngRect {
            requireArgument { p1.isValid }
            requireArgument { p2.isValid }
            return S2LatLngRect(
                    R1Interval.fromPointPair(p1.lat().radians, p2.lat().radians),
                    S1Interval.fromPointPair(p1.lng().radians, p2.lng().radians)
            )
        }

        /**
         * The canonical empty rectangle.
         * Empty: lat_lo=1, lat_hi=0, lng_lo=Pi, lng_hi=-Pi (radians)
         */
        fun empty(): S2LatLngRect = S2LatLngRect(R1Interval.empty(), S1Interval.empty())

        /**
         * The canonical full rectangle.
         * Full: lat_lo=-Pi/2, lat_hi=Pi/2, lng_lo=-Pi, lng_hi=Pi (radians)
         */
        fun full(): S2LatLngRect = S2LatLngRect(fullLat(), fullLng())

        /** The full allowable range of latitudes.  */
        @JvmStatic
        fun fullLat(): R1Interval {
            return R1Interval(-M_PI_2, M_PI_2)
        }

        /**
         * The full allowable range of longitudes.
         */
        @JvmStatic
        fun fullLng(): S1Interval {
            return S1Interval.full()
        }

        /**
         * Return a latitude-longitude rectangle that contains the edge from "a" to
         * "b". Both points must be unit-length. Note that the bounding rectangle of
         * an edge can be larger than the bounding rectangle of its endpoints.
         */
        @JvmStatic
        fun fromEdge(a: S2Point, b: S2Point): S2LatLngRect {
            // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
            val r = fromPointPair(S2LatLng.fromPoint(a), S2LatLng.fromPoint(b))

            // Check whether the min/max latitude occurs in the edge interior.
            // We find the normal to the plane containing AB, and then a vector "dir" in
            // this plane that also passes through the equator. We use RobustCrossProd
            // to ensure that the edge normal is accurate even when the two points are
            // very close together.
            val ab = robustCrossProd(a, b)
            val dir = ab.crossProd(S2Point(0, 0, 1))
            val da = dir.dotProd(a)
            val db = dir.dotProd(b)
            if (da * db >= 0) {
                // Minimum and maximum latitude are attained at the vertices.
                return r
            }
            // Minimum/maximum latitude occurs in the edge interior. This affects the
            // latitude bounds but not the longitude bounds.
            val absLat = acos(abs(ab.z / ab.norm()))
            return if (da < 0) {
                S2LatLngRect(R1Interval(r.lat.lo, absLat), r.lng)
            } else {
                S2LatLngRect(R1Interval(-absLat, r.lat.hi), r.lng)
            }
        }

        /**
         * Return true if the edge AB intersects the given edge of constant longitude.
         */
        fun intersectsLngEdge(a: S2Point, b: S2Point, lat: R1Interval, lng: Double): Boolean {
            // Return true if the segment AB intersects the given edge of constant
            // longitude. The nice thing about edges of constant longitude is that
            // they are straight lines on the sphere (geodesics).

            // TODO replace with return S2::CrossingSign(
            //      a, b, S2LatLng::FromRadians(lat.lo(), lng).ToPoint(),
            //      S2LatLng::FromRadians(lat.hi(), lng).ToPoint()) > 0;

            return S2EdgeCrossings.crossingSign(a, b, S2LatLng.fromRadians(lat.lo, lng).toPoint(), S2LatLng.fromRadians(lat.hi, lng).toPoint()) > 0

            //return S2.simpleCrossing(a, b, S2LatLng.fromRadians(lat.lo, lng).toPoint(), S2LatLng.fromRadians(lat.hi, lng).toPoint())
        }

        /**
         * Return true if the edge AB intersects the given edge of constant latitude.
         */
        fun intersectsLatEdge(a: S2Point, b: S2Point, lat: Double, lng: S1Interval): Boolean {
            // Return true if the segment AB intersects the given edge of constant
            // latitude. Unfortunately, lines of constant latitude are curves on
            // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
            requireArgument { isUnitLength(a) }
            requireArgument { isUnitLength(b) }

            // First, compute the normal to the plane AB that points vaguely north.
            var z = robustCrossProd(a, b).normalize()
            if (z[2] < 0) {
                z = -z
            }

            // Extend this to an orthonormal frame (x,y,z) where x is the direction
            // where the great circle through AB achieves its maximium latitude.
            val y = robustCrossProd(z, S2Point(0, 0, 1)).normalize()
            val x = y.crossProd(z)
            checkState { isUnitLength(x) }
            assert(x[2] >= 0)

            // Compute the angle "theta" from the x-axis (in the x-y plane defined
            // above) where the great circle intersects the given line of latitude.
            val sinLat = sin(lat)
            if (abs(sinLat) >= x[2]) {
                return false // The great circle does not reach the given latitude.
            }
            assert(x[2] > 0);
            val cosTheta = sinLat / x[2]
            val sinTheta = sqrt(1 - cosTheta * cosTheta)
            val theta = atan2(sinTheta, cosTheta)

            // The candidate intersection points are located +/- theta in the x-y
            // plane. For an intersection to be valid, we need to check that the
            // intersection point is contained in the interior of the edge AB and
            // also that it is contained within the given longitude interval "lng".

            // Compute the range of theta values spanned by the edge AB.
            val abTheta = S1Interval.fromPointPair(
                    atan2(a.dotProd(y), a.dotProd(x)),
                    atan2(b.dotProd(y), b.dotProd(x))
            )
            if (abTheta.contains(theta)) {
                // Check if the intersection point is also in the given "lng" interval.
                val isect = x * cosTheta + y * sinTheta
                if (lng.contains(atan2(isect.y, isect.x))) {
                    return true
                }
            }
            if (abTheta.contains(-theta)) {
                // Check if the intersection point is also in the given "lng" interval.
                val isect = x * cosTheta - y * sinTheta
                if (lng.contains(atan2(isect.y, isect.x))) {
                    return true
                }
            }
            return false
        }


        // Return the directed Hausdorff distance from one longitudinal edge spanning
        // latitude range 'a_lat' to the other longitudinal edge spanning latitude
        // range 'b_lat', with their longitudinal difference given by 'lng_diff'.
        @JvmStatic
        fun getDirectedHausdorffDistance(lngDiff: Double, a: R1Interval, b: R1Interval): S1Angle {
            // By symmetry, we can assume a's longtitude is 0 and b's longtitude is
            // lng_diff. Call b's two endpoints b_lo and b_hi. Let H be the hemisphere
            // containing a and delimited by the longitude line of b. The Voronoi diagram
            // of b on H has three edges (portions of great circles) all orthogonal to b
            // and meeting at b_lo cross b_hi.
            // E1: (b_lo, b_lo cross b_hi)
            // E2: (b_hi, b_lo cross b_hi)
            // E3: (-b_mid, b_lo cross b_hi), where b_mid is the midpoint of b
            //
            // They subdivide H into three Voronoi regions. Depending on how longitude 0
            // (which contains edge a) intersects these regions, we distinguish two cases:
            // Case 1: it intersects three regions. This occurs when lng_diff <= M_PI_2.
            // Case 2: it intersects only two regions. This occurs when lng_diff > M_PI_2.
            //
            // In the first case, the directed Hausdorff distance to edge b can only be
            // realized by the following points on a:
            // A1: two endpoints of a.
            // A2: intersection of a with the equator, if b also intersects the equator.
            //
            // In the second case, the directed Hausdorff distance to edge b can only be
            // realized by the following points on a:
            // B1: two endpoints of a.
            // B2: intersection of a with E3
            // B3: farthest point from b_lo to the interior of D, and farthest point from
            //     b_hi to the interior of U, if any, where D (resp. U) is the portion
            //     of edge a below (resp. above) the intersection point from B2.

            requireGE(lngDiff, 0.0)
            requireLE(lngDiff, M_PI)

            if (lngDiff == 0.0) {
                return S1Angle.radians(a.directedHausdorffDistance(b))
            }

            // Assumed longtitude of b.
            val bLng = lngDiff
            // Two endpoints of b.
            val bLo = S2LatLng.fromRadians(b.lo, bLng).toPoint()
            val bHi = S2LatLng.fromRadians(b.hi, bLng).toPoint()

            // Handling of each case outlined at the top of the function starts here.
            // This is initialized a few lines below.
            var maxDistance: S1Angle

            // Cases A1 and B1.
            val aLo = S2LatLng.fromRadians(a.lo, 0.0).toPoint()
            val aHi = S2LatLng.fromRadians(a.hi, 0.0).toPoint()
            maxDistance = S2EdgeDistances.getDistance(aLo, bLo, bHi)
            maxDistance = maxOf(maxDistance, S2EdgeDistances.getDistance(aHi, bLo, bHi))

            if (lngDiff <= M_PI_2) {
                // Case A2.
                if (a.contains(0.0) && b.contains(0.0)) {
                    maxDistance = maxOf(maxDistance, S1Angle.radians(lngDiff))
                }
            } else {
                // Case B2.
                val p = getBisectorIntersection(b, bLng)
                val p_lat = S2LatLng.latitude(p).radians
                if (a.contains(p_lat)) {
                    maxDistance = maxOf(maxDistance, S1Angle(p, bLo))
                }

                // Case B3.
                if (p_lat > a.lo) {
                    maxDistance = maxOf(maxDistance, getInteriorMaxDistance(R1Interval(a.lo, min(p_lat, a.hi)), bLo))
                }
                if (p_lat < a.hi) {
                    maxDistance = maxOf(maxDistance, getInteriorMaxDistance(R1Interval(max(p_lat, a.lo), a.hi), bHi))
                }
            }

            return maxDistance;
        }

        // Return the intersection of longitude 0 with the bisector of an edge
        // on longitude 'lng' and spanning latitude range 'lat'.
        @JvmStatic
        fun getBisectorIntersection(lat: R1Interval, lng: Double): S2Point {
            val lngAbs = abs(lng)
            val latCenter = lat.center
            // A vector orthogonal to the bisector of the given longitudinal edge.
            val orthoBisector: S2LatLng = if (latCenter >= 0) {
                S2LatLng.fromRadians(latCenter - M_PI_2, lngAbs)
            } else {
                S2LatLng.fromRadians(-latCenter - M_PI_2, lngAbs - M_PI)
            }
            // A vector orthogonal to longitude 0.
            val orthoLng = S2Point(0, -1, 0)
            return robustCrossProd(orthoLng, orthoBisector.toPoint())
        }

        // Return max distance from a point b to the segment spanning latitude range
        // a_lat on longitude 0, if the max occurs in the interior of a_lat. Otherwise
        // return -1.
        @JvmStatic
        fun getInteriorMaxDistance(a_lat: R1Interval, b: S2Point): S1Angle {
            // Longitude 0 is in the y=0 plane. b.x() >= 0 implies that the maximum
            // does not occur in the interior of a_lat.
            if (a_lat.isEmpty || b.x >= 0) return S1Angle.radians(-1)

            // Project b to the y=0 plane. The antipodal of the normalized projection is
            // the point at which the maxium distance from b occurs, if it is contained
            // in a_lat.
            val intersectionPoint = S2Point(-b.x, 0.0, -b.z).normalize()
            if (a_lat.interiorContains(S2LatLng.latitude(intersectionPoint).radians)) {
                return S1Angle(b, intersectionPoint)
            } else {
                return S1Angle.radians(-1)
            }
        }

    }
}
