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
package dilivia.s2.edge

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireNE
import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.math.M_SQRT3
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2PointUtil.robustCrossProd
import dilivia.s2.S2PointUtil.simpleCCW
import dilivia.s2.S2Predicates
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sin
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * Defines a collection of functions for computing the distance to an edge, interpolating along an edge, projecting
 * points onto edges, etc.
 *
 * This class is a port of the s2edge_distances methods of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2EdgeDistances {

    private val logger = KotlinLogging.logger {  }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////                        (point, edge) functions                        ///////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Gets the minimum distance from X to any point on the edge AB. All arguments should be unit length. The result is
     * very accurate for small distances but may have some numerical error if the distance is large
     * (approximately Pi/2 or greater).  The case A == B is handled correctly.
     *
     * If you want to compare a distance against a fixed threshold, e.g.
     *       if (S2EdgeDistances.getDistance(x, a, b) < limit)
     * then it is significantly faster to use UpdateMinDistance() below.
     *
     * @param x A point.
     * @param a The begin point of the edge.
     * @param b The end point of the edge.
     * @return The minimum distance from X to the edge AB.
     */
    fun getDistance(x: S2Point, a: S2Point, b: S2Point): S1Angle {
        logger.trace { "getDistance(x = $x, a = $a, b = $b)" }
        val minDist = S1ChordAngle(Double.MAX_VALUE)
        alwaysUpdateMinDistance(x, a, b, minDist, true)
        logger.trace { "getDistance(x = $x, a = $a, b = $b) = $minDist" }
        return minDist.toAngle()
    }

    /**
     * Indicates if the distance from X to the edge AB is less than "limit". (Specify limit.successor() for "less than
     * or equal to".) This method is significantly faster than GetDistance(). If you want to compare against a fixed
     * S1Angle, you should convert it to an S1ChordAngle once and save the value, since this step is relatively
     * expensive.
     *
     * @param x A point.
     * @param a The begin point of the edge.
     * @param b The end point of the edge.
     * @param limit The limit distance.
     * @see S2Predicates.compareDistances for an exact version of this predicate.
     */
    fun isDistanceLess(x: S2Point, a: S2Point, b: S2Point, limit: S1ChordAngle): Boolean = updateMinDistance(x, a, b, limit.clone())

    /**
     * If the distance from X to the edge AB is less than "minDist", this method updates "minDist" and returns true.
     * Otherwise it returns false. The case A == B is handled correctly.
     *
     * Use this method when you want to compute many distances and keep track of the minimum. It is significantly faster
     * than using getDistance(), because (1) using S1ChordAngle is much faster than S1Angle, and (2) it can save a lot
     * of work by not actually computing the distance when it is obviously larger than the current minimum.
     *
     * @param x A point.
     * @param a The begin point of the edge.
     * @param b The end point of the edge.
     * @param minDist The minimum distance to update if the distance from X to AB is less than minDist.
     * @return true if the distance has been updated and false otherwise.
     */
    fun updateMinDistance(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle): Boolean = alwaysUpdateMinDistance(x, a, b, minDist, false)

    /**
     * If the maximum distance from X to the edge AB is greater than "max_dist", this method updates "max_dist" and
     * returns true. Otherwise it returns false.
     * The case A == B is handled correctly.
     *
     * @param x A point.
     * @param a The begin point of the edge.
     * @param b The end point of the edge.
     * @param maxDist The maximum distance to update if the distance from X to AB is greater than maxDist.
     * @return true if the distance has been updated and false otherwise.
     */
    fun updateMaxDistance(x: S2Point, a: S2Point, b: S2Point, maxDist: S1ChordAngle): Boolean {
        val dist = maxOf(S1ChordAngle.between(x, a), S1ChordAngle.between(x, b))
        if (dist > S1ChordAngle.right()) {
            alwaysUpdateMinDistance(-x, a, b, dist, true)
            dist.length2 = (S1ChordAngle.straight() - dist).length2
        }
        if (maxDist < dist) {
            maxDist.length2 = dist.length2
            return true
        }
        return false
    }

    /**
     * Gets the maximum error in the result of updateMinDistance (and associated functions such as
     * updateMinInteriorDistance, isDistanceLess, etc), assuming that all input points are normalized to within
     * the bounds guaranteed by S2Point.normalize(). The error can be added or subtracted from an S1ChordAngle "x"
     * using x.plusError(error).
     *
     * Note that accuracy goes down as the distance approaches 0 degrees or 180 degrees (for different reasons).
     * Near 0 degrees the error is acceptable for all practical purposes (about 1.2e-15 radians ~= 8 nanometers). For
     * exactly antipodal points the maximum error is quite high (0.5 meters), but this error drops rapidly as the points
     * move away from antipodality (approximately 1 millimeter for points that are 50 meters from antipodal,
     * and 1 micrometer for points that are 50km from antipodal).
     *
     * TODO(ericv): Currently the error bound does not hold for edges whose endpoints are antipodal to within about
     * 1e-15 radians (less than 1 micron). This could be fixed by extending S2PointUtil.robustCrossProd to use higher
     * precision when necessary.
     *
     * @param dist
     * @return The maximum error.
     */
    fun getUpdateMinDistanceMaxError(dist: S1ChordAngle): Double {
        // There are two cases for the maximum error in UpdateMinDistance(),
        // depending on whether the closest point is interior to the edge.
        return max(getUpdateMinInteriorDistanceMaxError(dist), dist.getS2PointConstructorMaxError())
    }

    // Returns the maximum error in the result of UpdateMinInteriorDistance,
    // assuming that all input points are normalized to within the bounds
    // guaranteed by S2Point::Normalize().  The error can be added or subtracted
    // from an S1ChordAngle "x" using x.PlusError(error).
    private fun getUpdateMinInteriorDistanceMaxError(dist: S1ChordAngle): Double {
        // If a point is more than 90 degrees from an edge, then the minimum
        // distance is always to one of the endpoints, not to the edge interior.
        if (dist >= S1ChordAngle.right()) return 0.0

        // This bound includes all source of error, assuming that the input points
        // are normalized to within the bounds guaranteed to S2Point::Normalize().
        // "a" and "b" are components of chord length that are perpendicular and
        // parallel to the plane containing the edge respectively.
        val b = min(1.0, 0.5 * dist.length2)
        val a = sqrt(b * (2 - b))
        return ((2.5 + 2 * M_SQRT3 + 8.5 * a) * a +
                (2 + 2 * M_SQRT3 / 3 + 6.5 * (1 - b)) * b +
                (23 + 16 / M_SQRT3) * DoubleType.epsilon) * DoubleType.epsilon
    }

    /**
     * Indicates if the minimum distance from X to the edge AB is attained at an interior point of AB (i.e., not an
     * endpoint), and that distance is less than "limit". (Specify limit.successor() for "less than or equal to".)
     */
    fun isInteriorDistanceLess(x: S2Point, a: S2Point, b: S2Point, limit: S1ChordAngle): Boolean = updateMinInteriorDistance(x, a, b, limit.clone())

    /**
     * If the minimum distance from X to AB is attained at an interior point of AB (i.e., not an endpoint), and that
     * distance is less than "min_dist", then this method updates "min_dist" and returns true.  Otherwise returns false.
     */
    fun updateMinInteriorDistance(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle): Boolean {
        val xa2 = (x - a).norm2()
        val xb2 = (x - b).norm2()
        return alwaysUpdateMinInteriorDistance(x, a, b, xa2, xb2, minDist, false)
    }

    /**
     * Returns the point along the edge AB that is closest to the point X. The fractional distance of this point along
     * the edge AB can be obtained using getDistanceFraction() above.  Requires that all vectors have
     * unit length.
     */
    fun project(x: S2Point, a: S2Point, b: S2Point): S2Point = project(x, a, b, robustCrossProd(a, b))

    /**
     * A slightly more efficient version of project() where the cross product of the two endpoints has been precomputed.
     * The cross product does not need to be normalized, but should be computed using S2::RobustCrossProd() for the
     * most accurate results.  Requires that x, a, and b have unit length.
     */
    fun project(x: S2Point, a: S2Point, b: S2Point, aCrossB: S2Point): S2Point {
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        requireArgument { isUnitLength(x) }

        // Find the closest point to X along the great circle through AB.
        val p = x - (x.dotProd(aCrossB) / aCrossB.norm2()) * aCrossB

        // If this point is on the edge AB, then it's the closest point.
        if (simpleCCW(aCrossB, a, p) && simpleCCW(p, b, aCrossB)) {
            return p.normalize()
        }
        // Otherwise, the closest point is either A or B.
        return if ((x - a).norm2() <= (x - b).norm2()) a else b
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////                     (point along edge) functions                      ///////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Given a point X and an edge AB, returns the distance ratio AX / (AX + BX). If X happens to be on the line segment
     * AB, this is the fraction "t" such that X == Interpolate(t, A, B).  Requires that A and B are distinct.
     */
    fun getDistanceFraction(x: S2Point, a: S2Point, b: S2Point): Double {
        assert(a <= b)
        val d0 = x.angle(a)
        val d1 = x.angle(b)
        return d0 / (d0 + d1)
    }

    // Returns the point X along the line segment AB whose distance from A is the
    // given fraction "t" of the distance AB.  Does NOT require that "t" be
    // between 0 and 1.  Note that all distances are measured on the surface of
    // the sphere, so this is more complicated than just computing (1-t)*a + t*b
    // and normalizing the result.
    fun interpolate(t: Double, a: S2Point, b: S2Point): S2Point {
        if (t == 0.0) return a
        if (t == 1.0) return b
        val ab = S1Angle(a, b)
        return interpolateAtDistance(t * ab, a, b)
    }

    // Like Interpolate(), except that the parameter "ax" represents the desired
    // distance from A to the result X rather than a fraction between 0 and 1.
    fun interpolateAtDistance(axAngle: S1Angle, a: S2Point, b: S2Point): S2Point {
        val ax = axAngle.radians

        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }

        // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
        // result is always perpendicular to A, even if A=B or A=-B, but it is not
        // necessarily unit length.  (We effectively normalize it below.)
        val normal = robustCrossProd(a, b)
        val tangent = normal.crossProd(a)
        requireNE(tangent, S2Point(0, 0, 0))

        // Now compute the appropriate linear combination of A and "tangent".  With
        // infinite precision the result would always be unit length, but we
        // normalize it anyway to ensure that the error is within acceptable bounds.
        // (Otherwise errors can build up when the result of one interpolation is
        // fed into another interpolation.)
        return (cos(ax) * a + (sin(ax) / tangent.norm()) * tangent).normalize()
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////                        (edge, edge) functions                         ///////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    // Like UpdateMinDistance(), but computes the minimum distance between the
    // given pair of edges.  (If the two edges cross, the distance is zero.)
    // The cases a0 == a1 and b0 == b1 are handled correctly.
    fun updateEdgePairMinDistance(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point, minDist: S1ChordAngle): Boolean {
        if (minDist == S1ChordAngle.zero()) {
            return false
        }
        if (S2EdgeCrossings.crossingSign(a0, a1, b0, b1) > 0) {
            minDist.length2 = S1ChordAngle.zero().length2
            return true
        }
        // Otherwise, the minimum distance is achieved at an endpoint of at least
        // one of the two edges.  We use "|" rather than "||" below to ensure that
        // all four possibilities are always checked.
        //
        // The calculation below computes each of the six vertex-vertex distances
        // twice (this could be optimized).
        return (updateMinDistance(a0, b0, b1, minDist) or
                updateMinDistance(a1, b0, b1, minDist) or
                updateMinDistance(b0, a0, a1, minDist) or
                updateMinDistance(b1, a0, a1, minDist))
    }

    // As above, but for maximum distances.  If one edge crosses the antipodal
    // reflection of the other, the distance is Pi.
    fun updateEdgePairMaxDistance(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point, maxDist: S1ChordAngle): Boolean {
        if (maxDist == S1ChordAngle.straight()) {
            return false
        }
        if (S2EdgeCrossings.crossingSign(a0, a1, -b0, -b1) > 0) {
            maxDist.length2 = S1ChordAngle.straight().length2
            return true
        }
        // Otherwise, the maximum distance is achieved at an endpoint of at least
        // one of the two edges.  We use "|" rather than "||" below to ensure that
        // all four possibilities are always checked.
        //
        // The calculation below computes each of the six vertex-vertex distances
        // twice (this could be optimized).
        return (updateMaxDistance(a0, b0, b1, maxDist) or
                updateMaxDistance(a1, b0, b1, maxDist) or
                updateMaxDistance(b0, a0, a1, maxDist) or
                updateMaxDistance(b1, a0, a1, maxDist))
    }

    // Returns the pair of points (a, b) that achieves the minimum distance
    // between edges a0a1 and b0b1, where "a" is a point on a0a1 and "b" is a
    // point on b0b1.  If the two edges intersect, "a" and "b" are both equal to
    // the intersection point.  Handles a0 == a1 and b0 == b1 correctly.
    fun getEdgePairClosestPoints(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): Pair<S2Point, S2Point> {
        if (S2EdgeCrossings.crossingSign(a0, a1, b0, b1) > 0) {
            val x = S2EdgeCrossings.getIntersection(a0, a1, b0, b1)
            return x to x
        }
        // We save some work by first determining which vertex/edge pair achieves
        // the minimum distance, and then computing the closest point on that edge.
        val minDist = S1ChordAngle(0.0)
        alwaysUpdateMinDistance(a0, b0, b1, minDist, true)
        var closestVertex = Vertex.A0
        if (updateMinDistance(a1, b0, b1, minDist)) { closestVertex = Vertex.A1; }
        if (updateMinDistance(b0, a0, a1, minDist)) { closestVertex = Vertex.B0; }
        if (updateMinDistance(b1, a0, a1, minDist)) { closestVertex = Vertex.B1; }
        return when (closestVertex) {
            Vertex.A0 -> Pair(a0, project(a0, b0, b1))
            Vertex.A1 -> Pair(a1, project(a1, b0, b1))
            Vertex.B0 -> Pair(project(b0, a0, a1), b0)
            Vertex.B1 -> Pair(project(b1, a0, a1), b1)
        }
    }
    enum class Vertex { A0, A1, B0, B1 }

    // Returns true if every point on edge B=b0b1 is no further than "tolerance"
    // from some point on edge A=a0a1.  Equivalently, returns true if the directed
    // Hausdorff distance from B to A is no more than "tolerance".
    // Requires that tolerance is less than 90 degrees.
    fun isEdgeBNearEdgeA(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point, tolerance: S1Angle): Boolean {
        assert(tolerance.radians <= M_PI / 2)
        assert(tolerance.radians >= 0)
        // The point on edge B=b0b1 furthest from edge A=a0a1 is either b0, b1, or
        // some interior point on B.  If it is an interior point on B, then it must be
        // one of the two points where the great circle containing B (circ(B)) is
        // furthest from the great circle containing A (circ(A)).  At these points,
        // the distance between circ(B) and circ(A) is the angle between the planes
        // containing them.

        var aOrtho = robustCrossProd(a0, a1).normalize()
        val aNearestB0 = project(b0, a0, a1, aOrtho)
        val aNearestB1 = project(b1, a0, a1, aOrtho)
        // If a_nearest_b0 and a_nearest_b1 have opposite orientation from a0 and a1,
        // we invert a_ortho so that it points in the same direction as a_nearest_b0 x
        // a_nearest_b1.  This helps us handle the case where A and B are oppositely
        // oriented but otherwise might be near each other.  We check orientation and
        // invert rather than computing a_nearest_b0 x a_nearest_b1 because those two
        // points might be equal, and have an unhelpful cross product.
        if (S2Predicates.sign(aOrtho, aNearestB0, aNearestB1) < 0) aOrtho = -aOrtho

        // To check if all points on B are within tolerance of A, we first check to
        // see if the endpoints of B are near A.  If they are not, B is not near A.
        val b0Distance = S1Angle(b0, aNearestB0)
        val b1Distance = S1Angle(b1, aNearestB1)
        if (b0Distance > tolerance || b1Distance > tolerance)
            return false

        // If b0 and b1 are both within tolerance of A, we check to see if the angle
        // between the planes containing B and A is greater than tolerance.  If it is
        // not, no point on B can be further than tolerance from A (recall that we
        // already know that b0 and b1 are close to A, and S2Edges are all shorter
        // than 180 degrees).  The angle between the planes containing circ(A) and
        // circ(B) is the angle between their normal vectors.
        val bOrtho = robustCrossProd(b0, b1).normalize()
        val planarAngle = S1Angle(aOrtho, bOrtho)
        if (planarAngle <= tolerance)
            return true

        // As planar_angle approaches M_PI, the projection of a_ortho onto the plane
        // of B approaches the null vector, and normalizing it is numerically
        // unstable.  This makes it unreliable or impossible to identify pairs of
        // points where circ(A) is furthest from circ(B).  At this point in the
        // algorithm, this can only occur for two reasons:
        //
        //  1.) b0 and b1 are closest to A at distinct endpoints of A, in which case
        //      the opposite orientation of a_ortho and b_ortho means that A and B are
        //      in opposite hemispheres and hence not close to each other.
        //
        //  2.) b0 and b1 are closest to A at the same endpoint of A, in which case
        //      the orientation of a_ortho was chosen arbitrarily to be that of a0
        //      cross a1.  B must be shorter than 2*tolerance and all points in B are
        //      close to one endpoint of A, and hence to A.
        //
        // The logic applies when planar_angle is robustly greater than M_PI/2, but
        // may be more computationally expensive than the logic beyond, so we choose a
        // value close to M_PI.
        if (planarAngle >= S1Angle.radians(M_PI - 0.01)) {
            return (S1Angle(b0, a0) < S1Angle(b0, a1)) == (S1Angle(b1, a0) < S1Angle(b1, a1))
        }

        // Finally, if either of the two points on circ(B) where circ(B) is furthest
        // from circ(A) lie on edge B, edge B is not near edge A.
        //
        // The normalized projection of a_ortho onto the plane of circ(B) is one of
        // the two points along circ(B) where it is furthest from circ(A).  The other
        // is -1 times the normalized projection.
        val furthest = (aOrtho - aOrtho.dotProd(bOrtho) * bOrtho).normalize()
        checkState { isUnitLength(furthest) }
        val furthestInv = -furthest

        // A point p lies on B if you can proceed from b_ortho to b0 to p to b1 and
        // back to b_ortho without ever turning right.  We test this for furthest and
        // furthest_inv, and return true if neither point lies on B.
        return !((S2Predicates.sign(bOrtho, b0, furthest) > 0 && S2Predicates.sign(furthest, b1, bOrtho) > 0) ||
                (S2Predicates.sign(bOrtho, b0, furthestInv) > 0 && S2Predicates.sign(furthestInv, b1, bOrtho) > 0))
    }

    //////////////////////////////////////////////////////////////////////////////
    ////////////////              internal functions              ///////////////
    /////////////////////////////////////////////////////////////////////////////

    // This function computes the distance from a point X to a line segment AB.
    // If the distance is less than "min_dist" or "always_update" is true, it
    // updates "min_dist" and returns true.  Otherwise it returns false.
    //
    // The "Always" in the function name refers to the template argument, i.e.
    // AlwaysUpdateMinDistance<true> always updates the given distance, while
    // AlwaysUpdateMinDistance<false> does not.  This optimization increases the
    // speed of GetDistance() by about 10% without creating code duplication.
    private fun alwaysUpdateMinDistance(x: S2Point, a: S2Point, b: S2Point, min_dist: S1ChordAngle, always_update: Boolean = true): Boolean {
        requireArgument { isUnitLength(x) }
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }

        val xa2 = (x - a).norm2()
        val xb2 = (x - b).norm2()
        if (alwaysUpdateMinInteriorDistance(x, a, b, xa2, xb2, min_dist, always_update)) {
            return true  // Minimum distance is attained along the edge interior.
        }
        // Otherwise the minimum distance is to one of the endpoints.
        val dist2 = min(xa2, xb2)
        if (!always_update && dist2 >= min_dist.length2) {
            return false
        }
        min_dist.length2 = dist2
        return true
    }

    // If the minimum distance from X to AB is attained at an interior point of AB
    // (i.e., not an endpoint), and that distance is less than "min_dist" or
    // "always_update" is true, then update "min_dist" and return true.  Otherwise
    // return false.
    //
    // The "Always" in the function name refers to the template argument, i.e.
    // AlwaysUpdateMinInteriorDistance<true> always updates the given distance,
    // while AlwaysUpdateMinInteriorDistance<false> does not.  This optimization
    // increases the speed of GetDistance() by about 10% without creating code
    // duplication.
    private fun alwaysUpdateMinInteriorDistance(x: S2Point, a: S2Point, b: S2Point, xa2: Double, xb2: Double, minDist: S1ChordAngle, alwaysUpdate: Boolean = true): Boolean {
        requireArgument { isUnitLength(x) }
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        requireArgument { xa2 == (x - a).norm2() }
        requireArgument { xb2 == (x - b).norm2() }

        // The closest point on AB could either be one of the two vertices (the
        // "vertex case") or in the interior (the "interior case").  Let C = A x B.
        // If X is in the spherical wedge extending from A to B around the axis
        // through C, then we are in the interior case.  Otherwise we are in the
        // vertex case.
        //
        // Check whether we might be in the interior case.  For this to be true, XAB
        // and XBA must both be acute angles.  Checking this condition exactly is
        // expensive, so instead we consider the planar triangle ABX (which passes
        // through the sphere's interior).  The planar angles XAB and XBA are always
        // less than the corresponding spherical angles, so if we are in the
        // interior case then both of these angles must be acute.
        //
        // We check this by computing the squared edge lengths of the planar
        // triangle ABX, and testing acuteness using the law of cosines:
        //
        //             max(XA^2, XB^2) < min(XA^2, XB^2) + AB^2
        //
        if (max(xa2, xb2) >= min(xa2, xb2) + (a - b).norm2()) {
            return false
        }
        // The minimum distance might be to a point on the edge interior.  Let R
        // be closest point to X that lies on the great circle through AB.  Rather
        // than computing the geodesic distance along the surface of the sphere,
        // instead we compute the "chord length" through the sphere's interior.
        // If the squared chord length exceeds min_dist.length2() then we can
        // return "false" immediately.
        //
        // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
        // is the point X projected onto the plane through the great circle AB.
        // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
        // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
        // is faster and the corresponding distance on the Earth's surface is
        // accurate to within 1% for distances up to about 1800km.
        val c = robustCrossProd(a, b)
        val c2 = c.norm2()
        val xDotC = x.dotProd(c)
        val xDotC2 = xDotC * xDotC
        if (!alwaysUpdate && xDotC2 > c2 * minDist.length2) {
            // The closest point on the great circle AB is too far away.  We need to
            // test this using ">" rather than ">=" because the actual minimum bound
            // on the distance is (x_dot_c2 / c2), which can be rounded differently
            // than the (more efficient) multiplicative test above.
            return false
        }
        // Otherwise we do the exact, more expensive test for the interior case.
        // This test is very likely to succeed because of the conservative planar
        // test we did initially.
        val cx = c.crossProd(x)
        if (a.dotProd(cx) >= 0 || b.dotProd(cx) <= 0) {
            return false
        }
        // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
        // This calculation has good accuracy for all chord lengths since it
        // is based on both the dot product and cross product (rather than
        // deriving one from the other).  However, note that the chord length
        // representation itself loses accuracy as the angle approaches Pi.
        val qr = 1 - sqrt(cx.norm2() / c2)
        val dist2 = (xDotC2 / c2) + (qr * qr)
        if (!alwaysUpdate && dist2 >= minDist.length2) {
            return false
        }
        minDist.length2 = dist2
        return true
    }

}
