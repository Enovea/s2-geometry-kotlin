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

import ch.obermuhlner.math.big.BigDecimalMath
import dilivia.PreConditions.checkLE
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.math.DoubleType
import dilivia.math.ExactFloatType
import dilivia.math.M_SQRT3
import dilivia.math.vectors.R3VectorDouble
import dilivia.math.vectors.R3VectorExactFloat
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2PointUtil.robustCrossProd
import dilivia.s2.S2Predicates
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.sqrt
import java.math.BigDecimal

/**
 * Defines functions related to determining whether two geodesic edges cross and for computing intersection points.
 *
 * The predicates crossingSign(), vertexCrossing(), and edgeOrVertexCrossing() are robust, meaning that they produce
 * correct, consistent results even in pathological cases.  See S2Predicates for additional robust predicates.
 *
 * See also S2EdgeCrosser (which efficiently tests an edge against a sequence of other edges) and S2CrossingEdgeQuery
 * (which uses an index to speed up the process).
 *
 * This class is a port of the s2edge_crossing methods of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2EdgeCrossings {

    private val logger = KotlinLogging.logger { }

    /**
     * kIntersectionError is an upper bound on the distance from the intersection point returned by GetIntersection()
     * to the true intersection point.
     */
    val kIntersectionError: S1Angle = S1Angle.radians(8.0 * DoubleType.roundingEpsilon)

    /**
     * This value can be used as the S2Builder snapRadius() to ensure that edges that have been displaced by up to
     * kIntersectionError are merged back together again. For example this can happen when geometry is intersected
     * with a set of tiles and then unioned.  It is equal to twice the intersection error because input edges might have
     * been displaced in opposite directions.
     */
    val kIntersectionMergeRadius: S1Angle = 2.0 * kIntersectionError

    /** The list of intersection methods implemented by getIntersection(). */
    enum class IntersectionMethod(val methodName: String) {
        SIMPLE("Simple"),
        SIMPLE_LD("Simple_ld"),
        STABLE("Stable"),
        STABLE_LD("Stable_ld"),
        EXACT("Exact"),
    }

    var kUseSimpleMethod = false
    var kUseLongDouble = false

    val intersectionMethodCounter = mutableMapOf<IntersectionMethod, Int>()

    /**
     * This function determines whether the edge AB intersects the edge CD.
     *
     * Note that if an edge is degenerate (A == B or C == D), the return value is 0, if two vertices from different
     * edges are the same and -1 otherwise.
     *
     * Properties of CrossingSign:
     *
     *  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
     *  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
     *  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
     *  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
     *
     * This function implements an exact, consistent perturbation model such that no three points are ever considered to
     * be collinear. This means that even if you have 4 points A, B, C, D that lie exactly in a line (say, around the
     * equator), C and D will be treated as being slightly to one side or the other of AB.  This is done in a way such
     * that the results are always consistent (see S2Predicates.sign).
     *
     * Note that if you want to check an edge against a collection of other edges, it is much more efficient to use an
     * S2EdgeCrosser.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param c The first endpoint of the edge CD.
     * @param d The second endpoint of the edge CD.
     * @return +1 if AB crosses CD at a point that is interior to both edges, 0 if any two vertices from different edges
     * are the same and -1 otherwise.
     */
    fun crossingSign(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Int {
        val crosser = S2EdgeCrosser(a, b, c)
        return crosser.crossingSign(d)
    }

    /**
     * Given two edges AB and CD where at least two vertices are identical (i.e. crossingSign(a,b,c,d) == 0), this
     * function defines whether the two edges "cross" in a such a way that point-in-polygon containment tests can be
     * implemented by counting the number of edge crossings.  The basic rule is that a "crossing" occurs if AB is
     * encountered after CD during a CCW sweep around the shared vertex starting from a fixed reference point.
     *
     * Note that according to this rule, if AB crosses CD then in general CD does not cross AB. However, this leads to
     * the correct result when counting polygon edge crossings.  For example, suppose that A,B,C are three consecutive
     * vertices of a CCW polygon.  If we now consider the edge crossings of a segment BP as P sweeps around B, the
     * crossing number changes parity exactly when BP crosses BA or BC.
     *
     * Useful properties of vertexCrossing (VC):
     *
     *  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
     *  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
     *  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
     *  (3) If exactly one of a,b equals one of c,d, then exactly one of VC(a,b,c,d) and VC(c,d,a,b) is true
     *
     * It is an error to call this method with 4 distinct vertices.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param c The first endpoint of the edge CD.
     * @param d The second endpoint of the edge CD.
     * @return true if AB crosses CD and false otherwise.
     */
    fun vertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Boolean {
        logger.trace { "vertexCrossing(a = $a, b = $b, c = $c, d = $d)" }
        // If A == B or C == D there is no intersection.  We need to check this
        // case first in case 3 or more input points are identical.
        if (a == b) {
            logger.trace { "vertexCrossing, A == B: No intersection." }
            return false
        }
        if (c == d) {
            logger.trace { "vertexCrossing, C == D: No intersection." }
            return false
        }

        // If any other pair of vertices is equal, there is a crossing if and only
        // if OrderedCCW() indicates that the edge AB is further CCW around the
        // shared vertex O (either A or B) than the edge CD, starting from an
        // arbitrary fixed reference point.
        //
        // Optimization: if AB=CD or AB=DC, we can avoid most of the calculations.
        when {
            a == c -> {
                logger.trace { "vertexCrossing, A == C" }
                if (b == d) {
                    logger.trace { "vertexCrossing, B == D" }; return true
                }
                val aOrtho = S2PointUtil.ortho(a)
                val orderedCCW = S2Predicates.orderedCCW(aOrtho, d, b, a)
                logger.trace { "vertexCrossing, a ortho = $aOrtho, S2Predicates.orderedCCW(aOrtho, d, b, a) = $orderedCCW" }
                return orderedCCW
            }
            b == d -> {
                logger.trace { "vertexCrossing, B == D" }
                val bOrtho = S2PointUtil.ortho(b)
                val orderedCCW = S2Predicates.orderedCCW(bOrtho, c, a, b)
                logger.trace { "vertexCrossing, b ortho = $bOrtho,  S2Predicates.orderedCCW(bOrtho, c, a, b) = $orderedCCW" }
                return orderedCCW
            }
            a == d -> {
                logger.trace { "vertexCrossing, A == D" }
                if (b == c) {
                    logger.trace { "vertexCrossing, B == C" }; return true
                }
                val aOrtho = S2PointUtil.ortho(a)
                val orderedCCW = S2Predicates.orderedCCW(aOrtho, c, b, a)
                logger.trace { "vertexCrossing, a ortho = $aOrtho, S2Predicates.orderedCCW(aOrtho, c, b, a) = $orderedCCW" }
                return orderedCCW
            }
            b == c -> {
                logger.trace { "vertexCrossing, B == C" }
                val bOrtho = S2PointUtil.ortho(b);
                val orderedCCW = S2Predicates.orderedCCW(bOrtho, d, a, b)
                logger.trace { "vertexCrossing, b ortho = $bOrtho, S2Predicates.orderedCCW(bOrtho, d, a, b) = $orderedCCW" }
                return orderedCCW
            }
            else -> {
                logger.error("VertexCrossing called with 4 distinct vertices")
                return false
            }
        }

    }

    /**
     * A convenience function that calls crossingSign() to handle cases where all four vertices are distinct, and
     * vertexCrossing() to handle cases where two or more vertices are the same. This defines a crossing function such
     * that point-in-polygon containment tests can be implemented by simply counting edge crossings.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param c The first endpoint of the edge CD.
     * @param d The second endpoint of the edge CD.
     * @return true if AB crosses CD and false otherwise.
     */
    fun edgeOrVertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Boolean {
        val crossing = crossingSign(a, b, c, d)
        return crossing > 0 || (crossing == 0 && vertexCrossing(a, b, c, d))
    }

    /**
     * Given two edges AB and CD such that crossingSign(A, B, C, D) > 0, returns their intersection point.
     * Useful properties of getIntersection (GI):
     *
     *  (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
     *  (2) GI(c,d,a,b) == GI(a,b,c,d)
     *
     * The returned intersection point X is guaranteed to be very close to the true intersection point of AB and CD,
     * even if the edges intersect at a very small angle.  See "kIntersectionError" below for details.
     *
     * REQUIRES: crossingSign(A, B, C, D) > 0
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param c The first endpoint of the edge CD.
     * @param d The second endpoint of the edge CD.
     */
    fun getIntersection(a: S2Point, b: S2Point, c: S2Point, d: S2Point): S2Point {
//        requireArgument { crossingSign(a, b, c, d) > 0 }

        // It is difficult to compute the intersection point of two edges accurately
        // when the angle between the edges is very small.  Previously we handled
        // this by only guaranteeing that the returned intersection point is within
        // kIntersectionError of each edge.  However, this means that when the edges
        // cross at a very small angle, the computed result may be very far from the
        // true intersection point.
        //
        // Instead this function now guarantees that the result is always within
        // kIntersectionError of the true intersection.  This requires using more
        // sophisticated techniques and in some cases extended precision.
        //
        // Three different techniques are implemented, but only two are used:
        //
        //  - getIntersectionSimple() computes the intersection point using
        //    numerically stable cross products in "long double" precision.
        //
        //  - getIntersectionStable() computes the intersection point using
        //    projection and interpolation, taking care to minimize cancellation
        //    error.  This method exists in "double" and "long double" versions.
        //
        //  - getIntersectionExact() computes the intersection point using exact
        //    arithmetic and converts the final result back to an S2Point.
        //
        // We don't actually use the first method (GetIntersectionSimple) because it
        // turns out that GetIntersectionStable() is twice as fast and also much
        // more accurate (even in double precision).  The "long double" version
        // (only available on Intel platforms) uses 80-bit precision and is about
        // twice as slow.  The exact arithmetic version is about 100x slower.
        //
        // So our strategy is to first call GetIntersectionStable() in double
        // precision; if that doesn't work and this platform supports "long double",
        // then we try again in "long double"; if that doesn't work then we fall
        // back to exact arithmetic.

        var result: S2Point? = null
        var method: IntersectionMethod = IntersectionMethod.STABLE

        if (kUseSimpleMethod) {
            result = getIntersectionSimple(a, b, c, d)
            method = IntersectionMethod.SIMPLE
        }

        if (result == null) {
            result = getIntersectionStable(a, b, c, d)
            method = IntersectionMethod.STABLE
        }

        if (result == null) {
            result = getIntersectionExact(a, b, c, d)
            method = IntersectionMethod.EXACT
        }

        intersectionMethodCounter.compute(method) { _, value -> if (value == null) 1 else value + 1 }

        // Make sure the intersection point is on the correct side of the sphere.
        // Since all vertices are unit length, and edges are less than 180 degrees,
        // (a0 + a1) and (b0 + b1) both have positive dot product with the
        // intersection point.  We use the sum of all vertices to make sure that the
        // result is unchanged when the edges are swapped or reversed.
        if (result.dotProd((a + b) + (c + d)) < 0) result = -result;

        // Make sure that the intersection point lies on both edges.
        checkState({ approximatelyOrdered(a, result, b, kIntersectionError.radians) }, { "!approximatelyOrdered(a0 = $a, result = $result, a1 = $b, kIntersectionError.radians = ${kIntersectionError.radians})" })
        checkState { approximatelyOrdered(c, result, d, kIntersectionError.radians) }

        return result
    }

    //////////////////////////// Internal methods //////////////////////////////////////

    // If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
    // to within an error of at most kIntersectionError by this function, then set
    // "result" to the intersection point and return true.
    //
    // The intersection point is not guaranteed to have the correct sign
    // (i.e., it may be either "result" or "-result").
    private fun getIntersectionSimple(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): S2Point? {
        // The code below computes the intersection point as
        //
        //    (a0.CrossProd(a1)).CrossProd(b0.CrossProd(b1))
        //
        // except that it has better numerical stability and also computes a
        // guaranteed error bound.
        //
        // Each cross product is computed as (X-Y).CrossProd(X+Y) using unit-length
        // input vectors, which eliminates most of the cancellation error.  However
        // the error in the direction of the cross product can still become large if
        // the two points are extremely close together.  We can show that as long as
        // the length of the cross product is at least (16 * sqrt(3) + 24) * DBL_ERR
        // (about 6e-15), then the directional error is at most 5 * T_ERR (about
        // 3e-19 when T == "long double").  (DBL_ERR appears in the first formula
        // because the inputs are assumed to be normalized in double precision
        // rather than in the given type T.)
        //
        // The third cross product is different because its inputs already have some
        // error.  Letting "result_len" be the length of the cross product, it can
        // be shown that the error is at most
        //
        //   (2 + 2 * sqrt(3) + 12 / result_len) * T_ERR
        //
        // We want this error to be at most kIntersectionError, which is true as
        // long as "result_len" is at least kMinResultLen defined below.

        val kMinNormalLength = (16 * M_SQRT3 + 24) * DoubleType.epsilon
        val kMinResultLen = 12 / (kIntersectionError.radians / DoubleType.roundingEpsilon - (2 + 2 * M_SQRT3))

        // On some platforms "long double" is the same as "double", and on these
        // platforms this method always returns false (e.g. ARM, Win32).  Rather
        // than testing this directly, instead we look at kMinResultLen since this
        // is a direct measure of whether "long double" has sufficient accuracy to
        // be useful.  If kMinResultLen > 0.5, it means that this method will fail
        // even for edges that meet at an angle of 30 degrees.  (On Intel platforms
        // kMinResultLen corresponds to an intersection angle of about 0.04
        // degrees.)
        checkLE(kMinResultLen, 0.5)

        val (a_norm, a_length) = robustNormalWithLength(a0, a1)
        if (a_length >= kMinNormalLength) {
            val (b_norm, b_length) = robustNormalWithLength(b0, b1)
            if (b_length >= kMinNormalLength) {
                val (result, result_length) = robustNormalWithLength(a_norm, b_norm)
                if (result_length >= kMinResultLen) {
                    return result
                }
            }
        }
        return null
    }

    // If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
    // to within an error of at most kIntersectionError by this function, then set
    // "result" to the intersection point and return true.
    //
    // The intersection point is not guaranteed to have the correct sign
    // (i.e., it may be either "result" or "-result").
    private fun getIntersectionStable(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): S2Point? {
        // Sort the two edges so that (a0,a1) is longer, breaking ties in a
        // deterministic way that does not depend on the ordering of the endpoints.
        // This is desirable for two reasons:
        //  - So that the result doesn't change when edges are swapped or reversed.
        //  - It reduces error, since the first edge is used to compute the edge
        //    normal (where a longer edge means less error), and the second edge
        //    is used for interpolation (where a shorter edge means less error).
        val aLen2 = (a1 - a0).norm2()
        val bLen2 = (b1 - b0).norm2()
        if (aLen2 < bLen2 || (aLen2 == bLen2 && compareEdges(a0, a1, b0, b1))) {
            return getIntersectionStableSorted(b0, b1, a0, a1)
        } else {
            return getIntersectionStableSorted(a0, a1, b0, b1)
        }
    }

    // Returns whether (a0,a1) is less than (b0,b1) with respect to a total
    // ordering on edges that is invariant under edge reversals.
    fun compareEdges(a0: R3VectorDouble, a1: R3VectorDouble, b0: R3VectorDouble, b1: R3VectorDouble): Boolean {
        var pa0 = a0
        var pa1 = a1
        var pb0 = b0
        var pb1 = b1
        if (pa0 >= pa1) {
            val temp = pa0
            pa0 = pa1
            pa1 = temp
        }
        if (pb0 >= pb1) {
            val temp = pb0
            pb0 = pb1
            pb1 = temp
        }
        return pa0 < pb0 || (pa0 == pb0 && pb0 < pb1)
    }


    // Helper function for GetIntersectionStable().  It expects that the edges
    // (a0,a1) and (b0,b1) have been sorted so that the first edge is longer.
    fun getIntersectionStableSorted(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): S2Point? {
        requireArgument { ((a1 - a0).norm2() >= (b1 - b0).norm2()) }

        // Compute the normal of the plane through (a0, a1) in a stable way.
        val aNorm = (a0 - a1).crossProd(a0 + a1)
        val aNormLen = aNorm.norm()
        val bLen = (b1 - b0).norm()

        // Compute the projection (i.e., signed distance) of b0 and b1 onto the
        // plane through (a0, a1).  Distances are scaled by the length of a_norm.
        ;
        val (b0_dist, b0_error) = getProjection(b0, aNorm, aNormLen, a0, a1)
        val (b1_dist, b1_error) = getProjection(b1, aNorm, aNormLen, a0, a1)

        // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
        // |b0_dist - b1_dist|.  Note that b0_dist and b1_dist generally have
        // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
        // code below finds the intersection point by interpolating along the edge
        // (b0, b1) to a fractional distance of b0_dist / (b0_dist - b1_dist).
        //
        // It can be shown that the maximum error in the interpolation fraction is
        //
        //     (b0_dist * b1_error - b1_dist * b0_error) /
        //        (dist_sum * (dist_sum - error_sum))
        //
        // We save ourselves some work by scaling the result and the error bound by
        // "dist_sum", since the result is normalized to be unit length anyway.
        val distSum = abs(b0_dist - b1_dist)
        val errorSum = b0_error + b1_error
        if (distSum <= errorSum) {
            return null  // Error is unbounded in this case.
        }
        val x = b0_dist * b1 - b1_dist * b0
        val err = DoubleType.roundingEpsilon
        val error = bLen * abs(b0_dist * b1_error - b1_dist * b0_error) / (distSum - errorSum) + 2 * err * distSum

        // Finally we normalize the result, compute the corresponding error, and
        // check whether the total error is acceptable.
        val xLen2 = x.norm2()
        if (xLen2 < Double.MIN_VALUE) { // TODO(fmeurisse) Is it possible ?
            // If x.Norm2() is less than the minimum normalized value of T, x_len might
            // lose precision and the result might fail to satisfy S2::IsUnitLength().
            // TODO(ericv): Implement S2::RobustNormalize().
            return null
        }
        val xLen = sqrt(xLen2)
        val kMaxError = kIntersectionError.radians
        if (error > (kMaxError - err) * xLen) {
            return null
        }
        return (1 / xLen) * x
    }

    // Given a point X and a vector "a_norm" (not necessarily unit length),
    // compute x.DotProd(a_norm) and return a bound on the error in the result.
    // The remaining parameters allow this dot product to be computed more
    // accurately and efficiently.  They include the length of "a_norm"
    // ("a_norm_len") and the edge endpoints "a0" and "a1".
    private fun getProjection(x: S2Point, aNorm: S2Point, aNormLen: Double, a0: S2Point, a1: S2Point): Pair<Double, Double> {
        // The error in the dot product is proportional to the lengths of the input
        // vectors, so rather than using "x" itself (a unit-length vector) we use
        // the vectors from "x" to the closer of the two edge endpoints.  This
        // typically reduces the error by a huge factor.
        val x0 = x - a0
        val x1 = x - a1
        val x0Dist2 = x0.norm2()
        val x1Dist2 = x1.norm2()

        // If both distances are the same, we need to be careful to choose one
        // endpoint deterministically so that the result does not change if the
        // order of the endpoints is reversed.
        val dist: Double
        val result: Double
        if (x0Dist2 < x1Dist2 || (x0Dist2 == x1Dist2 && x0 < x1)) {
            dist = sqrt(x0Dist2)
            result = x0.dotProd(aNorm)
        } else {
            dist = sqrt(x1Dist2);
            result = x1.dotProd(aNorm)
        }
        // This calculation bounds the error from all sources: the computation of
        // the normal, the subtraction of one endpoint, and the dot product itself.
        // (DBL_ERR appears because the input points are assumed to be normalized in
        // double precision rather than in the given type T.)
        //
        // For reference, the bounds that went into this calculation are:
        // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * DBL_ERR) * T_ERR
        // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * T_ERR
        // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * T_ERR
        val error = (((3.5 + 2 * M_SQRT3) * aNormLen + 32 * M_SQRT3 * DoubleType.epsilon)
                * dist + 1.5 * abs(result)) * DoubleType.roundingEpsilon
        return result to error
    }

    /**
     * Compute the intersection point of (a0, a1) and (b0, b1) using exact arithmetic (BigDecima). Note that the result
     * is not exact because it is rounded to double precision. Also, the intersection point is not guaranteed to have
     * the correct sign (i.e., the return value may need to be negated).
     *
     * @param a0 First point of the edge a
     * @param a1 Second point of the edge a
     * @param b0 First point of the edge b
     * @param b1 Second point of the edge b
     * @return The intersection point
     */
    internal fun getIntersectionExact(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): S2Point {
        logger.trace { "getIntersectionExact(a0 = $a0, a1 = $a1, b0 = $b0, b1 = $b1)" }
        // Since we are using exact arithmetic, we don't need to worry about
        // numerical stability.
        val a0Xf = a0.toExactFloat()
        val a1Xf = a1.toExactFloat()
        val b0Xf = b0.toExactFloat()
        val b1Xf = b1.toExactFloat()
        val aNormXf = a0Xf.crossProd(a1Xf)
        val bNormXf = b0Xf.crossProd(b1Xf)
        val xXf = aNormXf.crossProd(bNormXf);

        logger.trace { """
            |
            |a0XF: ${S2Point(a0Xf)}
            |a1XF: ${S2Point(a1Xf)}
            |b0XF: ${S2Point(b0Xf)}
            |b1XF: ${S2Point(b1Xf)}
            |aNormXF: ${S2Point(aNormXf)}
            |bNormXF: ${S2Point(bNormXf)}
            |xXF: ${S2Point(xXf)}
        """.trimMargin() }

        // The final Normalize() call is done in double precision, which creates a
        // directional error of up to 2 * DBL_ERR.  (ToDouble() and Normalize() each
        // contribute up to DBL_ERR of directional error.)
        var x = convertS2PointFromExact(xXf)

        if (x == S2Point(0, 0, 0)) {
            // The two edges are exactly collinear, but we still consider them to be
            // "crossing" because of simulation of simplicity.  Out of the four
            // endpoints, exactly two lie in the interior of the other edge.  Of
            // those two we return the one that is lexicographically smallest.
            x = S2Point(10, 10, 10);  // Greater than any valid S2Point
            val aNorm = convertS2PointFromExact(aNormXf)
            val bNorm = convertS2PointFromExact(bNormXf)
            if (aNorm == S2Point(0, 0, 0) || bNorm == S2Point(0, 0, 0)) {
                // TODO(ericv): To support antipodal edges properly, we would need to
                // add an s2pred::CrossProd() function that computes the cross product
                // using simulation of simplicity and rounds the result to the nearest
                // floating-point representation.
                logger.error { "Exactly antipodal edges not supported by GetIntersection" }
            }
            if (S2Predicates.orderedCCW(b0, a0, b1, bNorm) && a0 < x) {
                logger.trace { "orderedCCW(b0, a0, b1, b_norm) && a0 < x => x = a0" }
                x = a0
            }
            if (S2Predicates.orderedCCW(b0, a1, b1, bNorm) && a1 < x) {
                logger.trace { "orderedCCW(b0, a1, b1, b_norm) && a1 < x => x = a1" }
                x = a1
            }
            if (S2Predicates.orderedCCW(a0, b0, a1, aNorm) && b0 < x) {
                logger.trace { "orderedCCW(a0, b0, a1, a_norm) && b0 < x => x = b0" }
                x = b0
            }
            if (S2Predicates.orderedCCW(a0, b1, a1, aNorm) && b1 < x)  {
                logger.trace { "orderedCCW(a0, b1, a1, a_norm) && b1 < x => x = b1" }
                x = b1
            }
        }
        logger.trace { "getIntersectionExact(a0 = $a0, a1 = $a1, b0 = $b0, b1 = $b1) = $x" }
        checkState { isUnitLength(x) }
        return x
    }

    // Computes the cross product of "x" and "y", normalizes it to be unit length,
    // and stores the result in "result".  Also returns the length of the cross
    // product before normalization, which is useful for estimating the amount of
    // error in the result.  For numerical stability, "x" and "y" should both be
    // approximately unit length.
    private fun robustNormalWithLength(x: S2Point, y: S2Point): Pair<S2Point, Double> {
        // This computes 2 * (x.CrossProd(y)), but has much better numerical
        // stability when "x" and "y" are unit length.
        val tmp: R3VectorDouble = (x - y).crossProd(x + y)
        val length = tmp.norm()
        val result = if (length != 0.0) {
            (1 / length) * tmp
        } else S2Point()
        return result to (0.5 * length)  // Since tmp == 2 * (x.CrossProd(y))
    }

    // Given three points "a", "x", "b", returns true if these three points occur
    // in the given order along the edge (a,b) to within the given tolerance.
    // More precisely, either "x" must be within "tolerance" of "a" or "b", or
    // when "x" is projected onto the great circle through "a" and "b" it must lie
    // along the edge (a,b) (i.e., the shortest path from "a" to "b").
    private fun approximatelyOrdered(a: S2Point, x: S2Point, b: S2Point, tolerance: Double): Boolean {
        if ((x - a).norm2() <= tolerance * tolerance) return true
        if ((x - b).norm2() <= tolerance * tolerance) return true
        return S2Predicates.orderedCCW(a, x, b, robustCrossProd(a, b).normalize())
    }

    private fun convertS2PointFromExact(xf: R3VectorExactFloat): S2Point {
        logger.trace { "Convert $xf to S2Point" }
        // If all components of "x" have absolute value less than about 1e-154,
        // then x.Norm2() is zero in double precision due to underflow.  Therefore
        // we need to scale "x" by an appropriate power of 2 before the conversion.
        val x = S2Point(xf[0].toDouble(), xf[1].toDouble(), xf[2].toDouble())
        if (x.norm2() > 0) {
            logger.trace { "Norm of $x > 0 => $x" }
            return x.normalize()
        }

        logger.trace { "Norm of $x == 0" }

        // Scale so that the largest component magnitude is in the range [0.5, 1).
        var exp = -181;
        for (i in 0..2) {
            if (BigDecimalMath.mantissa(xf[i]) != BigDecimal.ZERO) {
                exp = max(BigDecimalMath.exponent(xf[i]), exp)
            }
            //if (xf[i].is_normal()) exp = std::max(exp, xf[i].exp());
        }
        if (exp < -180) {
            return S2Point(0, 0, 0);
        }
        return S2Point(
                (xf[0] * 10.0.toBigDecimal(ExactFloatType.mathContext).pow(-exp)).toDouble(),
                (xf[1] * 10.0.toBigDecimal(ExactFloatType.mathContext).pow(-exp)).toDouble(),
                (xf[2] * 10.0.toBigDecimal(ExactFloatType.mathContext).pow(-exp)).toDouble()
        ).normalize()
    }

}