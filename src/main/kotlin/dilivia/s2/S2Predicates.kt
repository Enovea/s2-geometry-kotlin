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
package dilivia.s2

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.math.DoubleType
import dilivia.math.ExactFloatType
import dilivia.math.M_SQRT1_2
import dilivia.math.M_SQRT2
import dilivia.math.M_SQRT3
import dilivia.math.vectors.R3VectorExactFloat
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sqrt
import java.math.BigDecimal

/**
 * This object contains various predicates that are guaranteed to produce correct, consistent results.
 * They are also relatively efficient.  This is achieved by computing conservative error bounds and falling back
 * to high precision arithmetic when the result is uncertain.
 * Such predicates are useful in implementing robust algorithms.
 *
 * See also S2EdgeCrosser, which implements various exact edge-crossing predicates more efficiently than can be done
 * here.
 *
 * This class is a port of the s2predicates namespace of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2Predicates {

    /** Logger of the class. */
    private val logger = KotlinLogging.logger { }

    /** A predefined S1ChordAngle representing (approximately) 45 degrees. */
    private val k45Degrees = S1ChordAngle.fromLength2(2 - M_SQRT2)

    // ---------------------------------------------------------------------------------
    //  Sign
    // ---------------------------------------------------------------------------------

    /**
     * Indicates if the points A, B, C are counterclockwise.
     *
     * This function is essentially like taking the sign of the determinant of ABC, except that it has
     * additional logic to make sure that the above properties hold even when the three points are coplanar, and to
     * deal with the limitations of floating-point arithmetic.
     *
     * Sign satisfies the following conditions:
     *
     * (1) sign(a,b,c) == 0 if and only if a == b, b == c, or c == a
     * (2) sign(b,c,a) == sign(a,b,c) for all a,b,c
     * (3) sign(c,b,a) == -sign(a,b,c) for all a,b,c
     *
     * In other words:
     *
     *   (1) The result is zero if and only if two points are the same.
     *   (2) Rotating the order of the arguments does not affect the result.
     *   (3) Exchanging any two arguments inverts the result.
     *
     * On the other hand, note that it is not true in general that sign(-a,b,c) == -sign(a,b,c), or any similar
     * identities involving antipodal points.
     *
     * @param a Point a
     * @param b Point b
     * @param c Point c
     * @param aCrossB Pre-computed a cross b. Computed if not specified.
     *
     * @return +1 if the points A, B, C are counterclockwise, -1 if the points are clockwise, and 0 if any two points
     * are the same.
     */
    fun sign(a: S2Point, b: S2Point, c: S2Point, aCrossB: S2Point = a.crossProd(b)): Int {
        // We don't need RobustCrossProd() here because Sign() does its own
        // error estimation and calls ExpensiveSign() if there is any uncertainty
        // about the result.
        logger.trace { ">> Sign(a = $a, b = $b, c = $c)" }
        var sign = triageSign(a, b, c, aCrossB)
        if (sign == 0) sign = expensiveSign(a, b, c)
        logger.trace { "<< Sign(a = $a, b = $b, c = $c) = $sign" }
        return sign
    }

    // ---- Low-Level Methods ----
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.

    /**
     * This version of Sign returns +1 if the points are definitely CCW, -1 if they are definitely CW, and 0 if two
     * points are identical or the result is uncertain. Uncertain cases can be resolved, if desired, by calling
     * expensiveSign.
     *
     * The purpose of this method is to allow additional cheap tests to be done, where possible, in order to avoid
     * calling expensiveSign unnecessarily.
     *
     * @param a Point a
     * @param b Point b
     * @param c Point c
     * @param aCrossB Pre-computed a cross b.
     *
     * @return +1 if the points are definitely CCW, -1 if they are definitely CW, and 0 if two
     * points are identical or the result is uncertain.
     */
    fun triageSign(a: S2Point, b: S2Point, c: S2Point, aCrossB: S2Point): Int {
        // kMaxDetError is the maximum error in computing (AxB).C where all vectors
        // are unit length.  Using standard inequalities, it can be shown that
        //
        //  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
        //
        // where "fl()" denotes a calculation done in floating-point arithmetic,
        // |x| denotes either absolute value or the L2-norm as appropriate, and
        // e = 0.5*DBL_EPSILON.  Similarly,
        //
        //  fl(B.C) = B.C + d where |d| <= (1.5*|B.C| + 1.5*|B|*|C|) * e .
        //
        // Applying these bounds to the unit-length vectors A,B,C and neglecting
        // relative error (which does not affect the sign of the result), we get
        //
        //  fl((AxB).C) = (AxB).C + d where |d| <= (2.5 + 2/sqrt(3)) * e
        //
        // which is about 3.6548 * e, or 1.8274 * DBL_EPSILON.
        val kMaxDetError = 1.8274 * DoubleType.epsilon
//        checkState { S2PointUtil.isUnitLength(a) && S2PointUtil.isUnitLength(b) && S2PointUtil.isUnitLength(c) }
        val det = aCrossB.dotProd(c)

        // Double-check borderline cases in debug mode.
        checkState {
            !debug || (
                    abs(det) <= kMaxDetError
                            || abs(det) >= 100 * kMaxDetError
                            || det * expensiveSign(a, b, c) > 0.0
                    )
        }

        val sign = if (det > kMaxDetError) 1 else if (det < -kMaxDetError) -1 else 0
        logger.trace { "triageSign(a = $a, b = $b, c = $c, a^b = $aCrossB) = $sign (det = $det, kMaxDetError = $kMaxDetError)" }
        return sign
    }

    // This function is invoked by Sign() if the sign of the determinant is
    // uncertain.  It always returns a non-zero result unless two of the input
    // points are the same.  It uses a combination of multiple-precision
    // arithmetic and symbolic perturbations to ensure that its results are
    // always self-consistent (cf. Simulation of Simplicity, Edelsbrunner and
    // Muecke).  The basic idea is to assign an infinitesimal symbolic
    // perturbation to every possible S2Point such that no three S2Points are
    // collinear and no four S2Points are coplanar.  These perturbations are so
    // small that they do not affect the sign of any determinant that was
    // non-zero before the perturbations.  If "perturb" is false, then instead
    // the exact sign of the unperturbed input points is returned, which can be
    // zero even when all three points are distinct.
    //
    // Unlike Sign(), this method does not require the input points to be
    // normalized.

    fun expensiveSign(a: S2Point, b: S2Point, c: S2Point, perturb: Boolean = true): Int {
        // Return zero if and only if two points are the same.  This ensures (1).
        if (a == b || b == c || c == a) return 0

        // Next we try recomputing the determinant still using floating-point
        // arithmetic but in a more precise way.  This is more expensive than the
        // simple calculation done by TriageSign(), but it is still *much* cheaper
        // than using arbitrary-precision arithmetic.  This optimization is able to
        // compute the correct determinant sign in virtually all cases except when
        // the three points are truly collinear (e.g., three points on the equator).
        val detSign = stableSign(a, b, c)
        if (detSign != 0) return detSign

        // TODO(ericv): Create a templated version of StableSign so that we can
        // retry in "long double" precision before falling back to ExactFloat.

        // TODO(ericv): Optimize ExactFloat so that it stores up to 32 bytes of
        // mantissa inline (without requiring memory allocation).

        // Otherwise fall back to exact arithmetic and symbolic permutations.
        return exactSign(a, b, c, perturb)
    }

    // Compute the determinant in a numerically stable way.  Unlike TriageSign(),
    // this method can usually compute the correct determinant sign even when all
    // three points are as collinear as possible.  For example if three points are
    // spaced 1km apart along a random line on the Earth's surface using the
    // nearest representable points, there is only a 0.4% chance that this method
    // will not be able to find the determinant sign.  The probability of failure
    // decreases as the points get closer together; if the collinear points are
    // 1 meter apart, the failure rate drops to 0.0004%.
    //
    // This method could be extended to also handle nearly-antipodal points (and
    // in fact an earlier version of this code did exactly that), but antipodal
    // points are rare in practice so it seems better to simply fall back to
    // exact arithmetic in that case.

    fun stableSign(a: S2Point, b: S2Point, c: S2Point): Int {
        val ab = b - a
        val bc = c - b
        val ca = a - c
        val ab2 = ab.norm2()
        val bc2 = bc.norm2()
        val ca2 = ca.norm2()

        // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
        // cyclically permuted if necessary so that AB is the longest edge.  (This
        // minimizes the magnitude of cross product.)  At the same time we also
        // compute the maximum error in the determinant.  Using a similar technique
        // to the one used for kMaxDetError, the error is at most
        //
        //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
        //
        // where e = 0.5 * DBL_EPSILON.  If the determinant magnitude is larger than
        // this value then we know its sign with certainty.
        val kDetErrorMultiplier = 3.2321 * DoubleType.epsilon  // see above
        val det: Double
        val maxError: Double
        when {
            // AB is the longest edge, so compute (A-C)x(B-C).C.
            ab2 >= bc2 && ab2 >= ca2 -> {
                det = -(ca.crossProd(bc).dotProd(c))
                maxError = kDetErrorMultiplier * sqrt(ca2 * bc2)
            }
            // BC is the longest edge, so compute (B-A)x(C-A).A.
            bc2 >= ca2 -> {
                det = -(ab.crossProd(ca).dotProd(a))
                maxError = kDetErrorMultiplier * sqrt(ab2 * ca2)
            }
            // CA is the longest edge, so compute (C-B)x(A-B).B.
            else -> {
                det = -(bc.crossProd(ab).dotProd(b))
                maxError = kDetErrorMultiplier * sqrt(bc2 * ab2)
            }
        }
        val sign = if (abs(det) <= maxError) 0 else if (det > 0) 1 else -1
        logger.trace { "stableSign(a = $a, b = $b, c = $c) = $sign (det = $det, maxError = $maxError" }
        return sign
    }

    // Compute the determinant using exact arithmetic and/or symbolic
    // permutations.  Requires that the three points are distinct.

    fun exactSign(a: S2Point, b: S2Point, c: S2Point, perturb: Boolean): Int {
        requireArgument { a != b && b != c && c != a }

        // Sort the three points in lexicographic order, keeping track of the sign
        // of the permutation.  (Each exchange inverts the sign of the determinant.)
        var permSign = 1
        var pa: S2Point = a
        var pb: S2Point = b
        var pc: S2Point = c
        if (pa > pb) {
            logger.trace { "exactSign(a: $a, b: $b, c: $c, perturb: $perturb): ($pa, $pb, $pc) pa > pb => permut pa and pb" }
            val temp = pa; pa = pb; pb = temp; permSign = -permSign
            logger.trace { "($pa, $pb, $pc), perm sign = $permSign" }
        }
        if (pb > pc) {
            logger.trace { "exactSign(a: $a, b: $b, c: $c, perturb: $perturb): ($pa, $pb, $pc) pb > pc => permut pb and pc" }
            val temp = pb; pb = pc; pc = temp; permSign = -permSign
            logger.trace { "($pa, $pb, $pc), perm sign = $permSign" }
        }
        if (pa > pb) {
            logger.trace { "exactSign(a: $a, b: $b, c: $c, perturb: $perturb): ($pa, $pb, $pc) pa > pb => permut pa and pb" }
            val temp = pa; pa = pb; pb = temp; permSign = -permSign
            logger.trace { "($pa, $pb, $pc), perm sign = $permSign" }
        }
        checkState { pa < pb && pb < pc }

        // Construct multiple-precision versions of the sorted points and compute
        // their exact 3x3 determinant.
        val xa = R3VectorExactFloat(pa.x, pa.y, pa.z)
        val xb = R3VectorExactFloat(pb.x, pb.y, pb.z)
        val xc = R3VectorExactFloat(pc.x, pc.y, pc.z)
        val xbCrossXc: R3VectorExactFloat = xb.crossProd(xc)
        val det = xa.dotProd(xbCrossXc)

        // If the exact determinant is non-zero, we're done.
        var detSign = det.signum()
        logger.trace { "exactSign(a = $a, b = $b, c = $c, perturb = $perturb) [xa = $xa, xb = $xb, xc = $xc, xbCrossXc = $xbCrossXc, det = $det, precision = ${det.precision()}, signum = $detSign]" }
        if (detSign == 0 && perturb) {
            // Otherwise, we need to resort to symbolic perturbations to resolve the
            // sign of the determinant.
            detSign = symbolicallyPerturbedSign(xa, xb, xc, xbCrossXc)
            checkState { 0 != detSign }
        }
        logger.trace { "permSign (=$permSign) * detSign (=$detSign) = ${permSign * detSign}" }
        return permSign * detSign
    }

    // The following function returns the sign of the determinant of three points
    // A, B, C under a model where every possible S2Point is slightly perturbed by
    // a unique infinitesmal amount such that no three perturbed points are
    // collinear and no four points are coplanar.  The perturbations are so small
    // that they do not change the sign of any determinant that was non-zero
    // before the perturbations, and therefore can be safely ignored unless the
    // determinant of three points is exactly zero (using multiple-precision
    // arithmetic).
    //
    // Since the symbolic perturbation of a given point is fixed (i.e., the
    // perturbation is the same for all calls to this method and does not depend
    // on the other two arguments), the results of this method are always
    // self-consistent.  It will never return results that would correspond to an
    // "impossible" configuration of non-degenerate points.
    //
    // Requirements:
    //   The 3x3 determinant of A, B, C must be exactly zero.
    //   The points must be distinct, with A < B < C in lexicographic order.
    //
    // Returns:
    //   +1 or -1 according to the sign of the determinant after the symbolic
    // perturbations are taken into account.
    //
    // Reference:
    //   "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on
    //   Graphics, 1990).
    //

    private fun symbolicallyPerturbedSign(
        a: R3VectorExactFloat, b: R3VectorExactFloat, c: R3VectorExactFloat,
        bCrossC: R3VectorExactFloat
    ): Int {
        logger.trace { "symbolicallyPerturbedSign(a: $a, b: $b, c: $c, b_cross_c: $bCrossC)" }
        // This method requires that the points are sorted in lexicographically
        // increasing order.  This is because every possible S2Point has its own
        // symbolic perturbation such that if A < B then the symbolic perturbation
        // for A is much larger than the perturbation for B.
        //
        // Alternatively, we could sort the points in this method and keep track of
        // the sign of the permutation, but it is more efficient to do this before
        // converting the inputs to the multi-precision representation, and this
        // also lets us re-use the result of the cross product B x C.
        checkState({ a < b && b < c }, { "a < b && b < c is false: a = $a, b = $b, c = $c" })

        // Every input coordinate x[i] is assigned a symbolic perturbation dx[i].
        // We then compute the sign of the determinant of the perturbed points,
        // i.e.
        //               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
        //               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
        //               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
        //
        // The perturbations are chosen such that
        //
        //   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
        //
        // where each perturbation is so much smaller than the previous one that we
        // don't even need to consider it unless the coefficients of all previous
        // perturbations are zero.  In fact, it is so small that we don't need to
        // consider it unless the coefficient of all products of the previous
        // perturbations are zero.  For example, we don't need to consider the
        // coefficient of db[1] unless the coefficient of db[2]*da[0] is zero.
        //
        // The follow code simply enumerates the coefficients of the perturbations
        // (and products of perturbations) that appear in the determinant above, in
        // order of decreasing perturbation magnitude.  The first non-zero
        // coefficient determines the sign of the result.  The easiest way to
        // enumerate the coefficients in the correct order is to pretend that each
        // perturbation is some tiny value "eps" raised to a power of two:
        //
        // eps**    1      2      4      8     16     32     64     128    256
        //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
        //
        // Essentially we can then just count in binary and test the corresponding
        // subset of perturbations at each step.  So for example, we must test the
        // coefficient of db[2]*da[0] before db[1] because eps**12 > eps**16.
        //
        // Of course, not all products of these perturbations appear in the
        // determinant above, since the determinant only contains the products of
        // elements in distinct rows and columns.  Thus we don't need to consider
        // da[2]*da[1], db[1]*da[1], etc.  Furthermore, sometimes different pairs of
        // perturbations have the same coefficient in the determinant; for example,
        // da[1]*db[0] and db[1]*da[0] have the same coefficient (c[2]).  Therefore
        // we only need to test this coefficient the first time we encounter it in
        // the binary order above (which will be db[1]*da[0]).
        //
        // The sequence of tests below also appears in Table 4-ii of the paper
        // referenced above, if you just want to look it up, with the following
        // translations: [a,b,c] -> [i,j,k] and [0,1,2] -> [1,2,3].  Also note that
        // some of the signs are different because the opposite cross product is
        // used (e.g., B x C rather than C x B).

        var detSign = bCrossC[2].signum()                         // da[2]
        if (detSign != 0) {
            logger.trace { "sign(bCrossC[2]) = ${bCrossC[2].signum()} != 0" }
            return detSign
        }
        detSign = bCrossC[1].signum()                             // da[1]
        if (detSign != 0) return detSign
        detSign = bCrossC[0].signum()                             // da[0]
        if (detSign != 0) return detSign

        detSign = (c[0] * a[1] - c[1] * a[0]).signum()            // db[2]
        if (detSign != 0) return detSign
        detSign = c[0].signum()                                   // db[2] * da[1]
        if (detSign != 0) return detSign
        detSign = -(c[1].signum())                                // db[2] * da[0]
        if (detSign != 0) return detSign
        detSign = (c[2] * a[0] - c[0] * a[2]).signum()            // db[1]
        if (detSign != 0) return detSign
        detSign = c[2].signum()                                  // db[1] * da[0]
        if (detSign != 0) return detSign
        // The following test is listed in the paper, but it is redundant because
        // the previous tests guarantee that C == (0, 0, 0).
        checkState { 0 == (c[1] * a[2] - c[2] * a[1]).signum() }  // db[0]

        detSign = (a[0] * b[1] - a[1] * b[0]).signum()            // dc[2]
        if (detSign != 0) return detSign
        detSign = -(b[0].signum())                                // dc[2] * da[1]
        if (detSign != 0) return detSign
        detSign = b[1].signum()                                   // dc[2] * da[0]
        if (detSign != 0) return detSign
        detSign = a[0].signum()                                   // dc[2] * db[1]
        if (detSign != 0) return detSign
        return 1                                                  // dc[2] * db[1] * da[0]
    }

    // ----------------------------------------------------
    //  Ordered CCW
    // ----------------------------------------------------

    // Given 4 points on the unit sphere, return true if the edges OA, OB, and
    // OC are encountered in that order while sweeping CCW around the point O.
    // You can think of this as testing whether A <= B <= C with respect to the
    // CCW ordering around O that starts at A, or equivalently, whether B is
    // contained in the range of angles (inclusive) that starts at A and extends
    // CCW to C.  Properties:
    //
    //  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
    //  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
    //  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
    //  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
    //  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
    fun orderedCCW(a: S2Point, b: S2Point, c: S2Point, o: S2Point): Boolean {
        // The last inequality below is ">" rather than ">=" so that we return true
        // if A == B or B == C, and otherwise false if A == C.  Recall that
        // Sign(x,y,z) == -Sign(z,y,x) for all x,y,z.

        var sum = 0
        if (sign(b, o, a) >= 0) ++sum
        if (sign(c, o, b) >= 0) ++sum
        if (sign(a, o, c) > 0) ++sum

        logger.trace { "orderedCCW(a = $a, b = $b, c = $c, o = $o) = $sum >= 2" }
        return sum >= 2
    }

    // Returns -1, 0, or +1 according to whether AX < BX, A == B, or AX > BX
    // respectively.  Distances are measured with respect to the positions of X,
    // A, and B as though they were reprojected to lie exactly on the surface of
    // the unit sphere.  Furthermore, this method uses symbolic perturbations to
    // ensure that the result is non-zero whenever A != B, even when AX == BX
    // exactly, or even when A and B project to the same point on the sphere.
    // Such results are guaranteed to be self-consistent, i.e. if AB < BC and
    // BC < AC, then AB < AC.
    fun compareDistances(x: S2Point, a: S2Point, b: S2Point): Int {
        logger.trace {
            """
            |--------------------------------------------
            | CompareDistances(x = $x, a = $a, b = $b)
            |--------------------------------------------
        """.trimMargin()
        }
        // We start by comparing distances using dot products (i.e., cosine of the
        // angle), because (1) this is the cheapest technique, and (2) it is valid
        // over the entire range of possible angles.  (We can only use the sin^2
        // technique if both angles are less than 90 degrees or both angles are
        // greater than 90 degrees.)
        var sign = triageCompareCosDistances(x, a, b)
        logger.trace { "triageCompareCosDistances(x, a, b) = $sign" }
        if (sign != 0) return sign

        // Optimization for (a == b) to avoid falling back to exact arithmetic.
        if (a == b) {
            logger.trace { "a == b => sign = 0" }
            return 0
        }

        // It is much better numerically to compare distances using cos(angle) if
        // the distances are near 90 degrees and sin^2(angle) if the distances are
        // near 0 or 180 degrees.  We only need to check one of the two angles when
        // making this decision because the fact that the test above failed means
        // that angles "a" and "b" are very close together.
        val cosAx = a.dotProd(x)
        if (cosAx > M_SQRT1_2) {
            // Angles < 45 degrees.
            sign = compareSin2Distances(x, a, b)
            logger.trace { "cosAx > M_SQRT1_2: sign = compareSin2Distances(x, a, b) = $sign" }
        } else if (cosAx < -M_SQRT1_2) {
            // Angles > 135 degrees.  sin^2(angle) is decreasing in this range.
            sign = -compareSin2Distances(x, a, b)
            logger.trace { "cosAx < -M_SQRT1_2: sign = -compareSin2Distances(x, a, b) = $sign" }
        }
//        else {
//            // We've already tried double precision, so continue with "long double".
//            sign = triageCompareCosDistances(x.toDecimal128(), a.toDecimal128(), b.toDecimal128())
//            logger.trace { "triageCompareCosDistances(x.toLongDouble(), a.toLongDouble(), b.toLongDouble()) = $sign" }
//        }
        if (sign != 0) return sign
        sign = exactCompareDistances(x.toExactFloat(), a.toExactFloat(), b.toExactFloat())
        logger.trace { "exactCompareDistances(x.toExactFloat(), a.toExactFloat(), b.toExactFloat()) = $sign" }
        if (sign != 0) return sign

        sign = symbolicCompareDistances(x, a, b)
        logger.trace { "symbolicCompareDistances(x, a, b) = $sign" }
        return sign
    }

    // Returns -1, 0, or +1 according to whether the distance XY is less than,
    // equal to, or greater than "r" respectively.  Distances are measured with
    // respect the positions of all points as though they are projected to lie
    // exactly on the surface of the unit sphere.
    fun compareDistance(x: S2Point, y: S2Point, r: S1ChordAngle): Int {
        // As with CompareDistances(), we start by comparing dot products because
        // the sin^2 method is only valid when the distance XY and the limit "r" are
        // both less than 90 degrees.
        var sign = triageCompareCosDistance(x, y, r.length2)
        if (sign != 0) return sign

        // Unlike with CompareDistances(), it's not worth using the sin^2 method
        // when the distance limit is near 180 degrees because the S1ChordAngle
        // representation itself has has a rounding error of up to 2e-8 radians for
        // distances near 180 degrees.
        if (r < k45Degrees) {
            sign = triageCompareSin2Distance(x, y, r.length2)
            if (sign != 0) return sign
            //sign = triageCompareSin2Distance(x.toLongDouble(), y.toLongDouble(), LongDoubleType.fromDouble(r.length2))
        }
        /*else {
            sign = triageCompareCosDistance(x.toLongDouble(), y.toLongDouble(), LongDoubleType.fromDouble(r.length2))
        }*/
        if (sign != 0) return sign
        return exactCompareDistance(x.toExactFloat(), y.toExactFloat(), ExactFloatType.cast(r.length2))
    }

    // Returns -1, 0, or +1 according to whether the distance from the point X to
    // the edge A is less than, equal to, or greater than "r" respectively.
    // Distances are measured with respect the positions of all points as though
    // they were projected to lie exactly on the surface of the unit sphere.
    //
    // REQUIRES: A0 and A1 do not project to antipodal points (e.g., A0 == -A1).
    //           This requires that (A0 != C * A1) for any constant C < 0.
    //
    // NOTE(ericv): All of the predicates defined here could be extended to handle
    // edges consisting of antipodal points by implementing additional symbolic
    // perturbation logic (similar to Sign) in order to rigorously define the
    // direction of such edges.
    fun compareEdgeDistance(x: S2Point, a0: S2Point, a1: S2Point, r: S1ChordAngle): Int {
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactCompareEdgeDistance.)
        assert(a0 != -a1)

        val sign = triageCompareEdgeDistance(x, a0, a1, r.length2)
        if (sign != 0) return sign

        // Optimization for the case where the edge is degenerate.
        if (a0 == a1) return compareDistance(x, a0, r)

        //sign = triageCompareEdgeDistance(x.toLongDouble(), a0.toLongDouble(), a1.toLongDouble(), ExactFloatType.fromDouble(r.length2))
        if (sign != 0) return sign
        return exactCompareEdgeDistance(x, a0, a1, r)
    }

    // Returns -1, 0, or +1 according to whether the normal of edge A has
    // negative, zero, or positive dot product with the normal of edge B.  This
    // essentially measures whether the edges A and B are closer to proceeding in
    // the same direction or in opposite directions around the sphere.
    //
    // This method returns an exact result, i.e. the result is zero if and only if
    // the two edges are exactly perpendicular or at least one edge is degenerate.
    // (i.e., both edge endpoints project to the same point on the sphere).
    //
    // CAVEAT: This method does not use symbolic perturbations.  Therefore it can
    // return zero even when A0 != A1 and B0 != B1, e.g. if (A0 == C * A1) exactly
    // for some constant C > 0 (which is possible even when both points are
    // considered "normalized").
    //
    // REQUIRES: Neither edge can consist of antipodal points (e.g., A0 == -A1)
    //           (see comments in CompareEdgeDistance).
    fun compareEdgeDirections(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): Int {
        // Check that no edge consists of antipodal points.  (This catches the most
        // common case -- the full test is in ExactCompareEdgeDirections.)
        assert(a0 != -a1)
        assert(b0 != -b1)

        var sign = triageCompareEdgeDirections(a0, a1, b0, b1)
        if (sign != 0) return sign

        // Optimization for the case where either edge is degenerate.
        if (a0 == a1 || b0 == b1) return 0

        //sign = triageCompareEdgeDirections(a0.toLongDouble(), a1.toLongDouble(), b0.toLongDouble(), b1.toLongDouble())
        if (sign != 0) return sign
        return exactCompareEdgeDirections(a0.toExactFloat(), a1.toExactFloat(), b0.toExactFloat(), b1.toExactFloat())
    }

    // Returns Sign(X0, X1, Z) where Z is the circumcenter of triangle ABC.
    // The return value is +1 if Z is to the left of edge X, and -1 if Z is to the
    // right of edge X.  The return value is zero if A == B, B == C, or C == A
    // (exactly), and also if X0 and X1 project to identical points on the sphere
    // (e.g., X0 == X1).
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when Z lies exactly on edge X.
    //
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    fun edgeCircumcenterSign(x0: S2Point, x1: S2Point, a: S2Point, b: S2Point, c: S2Point): Int {
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactEdgeCircumcenterSign.)
        assert(x0 != -x1)

        val abc_sign = sign(a, b, c)
        var sign = triageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign)
        if (sign != 0) return sign

        // Optimization for the cases that are going to return zero anyway, in order
        // to avoid falling back to exact arithmetic.
        if (x0 == x1 || a == b || b == c || c == a) return 0

        //sign = triageEdgeCircumcenterSign(x0.toLongDouble(), x1.toLongDouble(), a.toLongDouble(), b.toLongDouble(), c.toLongDouble(), abc_sign)
        if (sign != 0) return sign
        sign = exactEdgeCircumcenterSign(
            x0.toExactFloat(),
            x1.toExactFloat(),
            a.toExactFloat(),
            b.toExactFloat(),
            c.toExactFloat(),
            abc_sign
        )
        if (sign != 0) return sign

        // Unlike the other methods, SymbolicEdgeCircumcenterSign does not depend
        // on the sign of triangle ABC.
        return symbolicEdgeCircumcenterSign(x0, x1, a, b, c)
    }

    // This is a specialized method that is used to compute the intersection of an
    // edge X with the Voronoi diagram of a set of points, where each Voronoi
    // region is intersected with a disc of fixed radius "r".
    //
    // Given two sites A and B and an edge (X0, X1) such that d(A,X0) < d(B,X0)
    // and both sites are within the given distance "r" of edge X, this method
    // intersects the Voronoi region of each site with a disc of radius r and
    // determines whether either region has an empty intersection with edge X.  It
    // returns FIRST if site A has an empty intersection, SECOND if site B has an
    // empty intersection, NEITHER if neither site has an empty intersection, or
    // UNCERTAIN if A == B exactly.  Note that it is not possible for both
    // intersections to be empty because of the requirement that both sites are
    // within distance r of edge X.  (For example, the only reason that Voronoi
    // region A can have an empty intersection with X is that site B is closer to
    // all points on X that are within radius r of site A.)
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when A and B lie on opposite sides of X
    // such that the Voronoi edge between them exactly coincides with edge X, or
    // when A and B are distinct but project to the same point on the sphere
    // (i.e., they are linearly dependent).
    //
    // REQUIRES: r < S1ChordAngle::Right() (90 degrees)
    // REQUIRES: s2pred::CompareDistances(x0, a, b) < 0
    // REQUIRES: s2pred::CompareEdgeDistance(a, x0, x1, r) <= 0
    // REQUIRES: s2pred::CompareEdgeDistance(b, x0, x1, r) <= 0
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    enum class Excluded { FIRST, SECOND, NEITHER, UNCERTAIN }

    fun getVoronoiSiteExclusion(a: S2Point, b: S2Point, x0: S2Point, x1: S2Point, r: S1ChordAngle): Excluded {
        requireArgument { r < S1ChordAngle.right() }
        //Assertions.assert { compareDistances(x0, a, b) < 0 }  // (implies a != b)
        //Assertions.assert { compareEdgeDistance(a, x0, x1, r) <= 0 }
        //Assertions.assert { compareEdgeDistance(b, x0, x1, r) <= 0 }
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactVoronoiSiteExclusion.)
        requireArgument { x0 != -x1 }

        // If one site is closer than the other to both endpoints of X, then it is
        // closer to every point on X.  Note that this also handles the case where A
        // and B are equidistant from every point on X (i.e., X is the perpendicular
        // bisector of AB), because CompareDistances uses symbolic perturbations to
        // ensure that either A or B is considered closer (in a consistent way).
        // This also ensures that the choice of A or B does not depend on the
        // direction of X.
        if (compareDistances(x1, a, b) < 0) {
            return Excluded.SECOND;  // Site A is closer to every point on X.
        }

        val result = triageVoronoiSiteExclusion(a, b, x0, x1, r.length2)
        if (result != Excluded.UNCERTAIN) return result

        return exactVoronoiSiteExclusion(
            a.toExactFloat(),
            b.toExactFloat(),
            x0.toExactFloat(),
            x1.toExactFloat(),
            ExactFloatType.cast(r.length2)
        )
    }

    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.


    // Returns cos(XY), and sets "error" to the maximum error in the result.
    // REQUIRES: "x" and "y" satisfy S2::IsNormalized().
    // A high precision "long double" version of the function above.

    internal fun getCosDistance(x: S2Point, y: S2Point): ValueAndError<Double> {
        val c = x.dotProd(y)
        val error = 9.5 * DoubleType.roundingEpsilon * abs(c) + 1.5 * DoubleType.roundingEpsilon
        return ValueAndError(c, error)
    }

    // A high precision "long double" version of the function above.

//    internal fun getCosDistance(x: R3VectorDecimal128, y: R3VectorDecimal128): ValueAndError<BigDecimal> {
//        // With "long double" precision it is worthwhile to compensate for length
//        // errors in "x" and "y", since they are only unit length to within the
//        // precision of "double".  (This would also reduce the error constant
//        // slightly in the method above but is not worth the additional effort.)
//        val c = x.dotProd(y) / Decimal128Type.sqrt(x.norm2() * y.norm2())
//        val error =
//            9.5.toDecimal128() * Decimal128Type.roundingEpsilon * Decimal128Type.abs(c) + 1.5.toDecimal128() * Decimal128Type.roundingEpsilon
//        return ValueAndError(c, error)
//    }

    // Returns sin**2(XY), where XY is the angle between X and Y, and sets "error"
    // to the maximum error in the result.
    //
    // REQUIRES: "x" and "y" satisfy S2::IsNormalized().

    fun getSin2Distance(x: S2Point, y: S2Point): ValueAndError<Double> {
        // The (x-y).CrossProd(x+y) trick eliminates almost all of error due to "x"
        // and "y" being not quite unit length.  This method is extremely accurate
        // for small distances; the *relative* error in the result is O(DBL_ERR) for
        // distances as small as DBL_ERR.
        val n = (x - y).crossProd(x + y)
        val d2 = 0.25 * n.norm2()
        val error = ((21 + 4 * M_SQRT3) * DoubleType.roundingEpsilon * d2 +
                32 * M_SQRT3 * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * sqrt(d2) +
                768 * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon)
        return ValueAndError(d2, error)
    }

    // A high precision "long double" version of the function above.

//    fun getSin2Distance(x: R3VectorDecimal128, y: R3VectorDecimal128): ValueAndError<BigDecimal> {
//        // In "long double" precision it is worthwhile to compensate for length
//        // errors in "x" and "y", since they are only unit length to within the
//        // precision of "double".  Otherwise the "d2" error coefficient below would
//        // be (16 * DBL_ERR + (5 + 4 * sqrt(3)) * LD_ERR), which is much larger.
//        // (Dividing by the squared norms of "x" and "y" would also reduce the error
//        // constant slightly in the double-precision version, but this is not worth
//        // the additional effort.)
//        val n = (x - y).crossProd(x + y)
//        val d2 = 0.25.toDecimal128() * n.norm2() / (x.norm2() * y.norm2())
//        val roundingEpsilonD128 = Decimal128Type.epsilon * BigDecimal.TEN * BigDecimal.TEN
//        val roundingEpsilonDouble = DoubleType.epsilon.toDecimal128() * BigDecimal.TEN * BigDecimal.TEN
//        val error = (13 + 4 * M_SQRT3).toDecimal128() * roundingEpsilonD128 * d2 +
//                32.toDecimal128() * M_SQRT3.toDecimal128() * roundingEpsilonDouble * roundingEpsilonD128 * Decimal128Type.sqrt(
//            d2
//        ) +
//                (768.toDecimal128() * roundingEpsilonDouble * roundingEpsilonDouble * roundingEpsilonD128 * roundingEpsilonD128)
//        val valueAndError = ValueAndError(d2, error)
//        logger.debug { valueAndError.toString() }
//        return valueAndError
//    }


    fun triageCompareCosDistances(x: S2Point, a: S2Point, b: S2Point): Int {
        val (cos_ax, cos_ax_error) = getCosDistance(a, x)
        val (cos_bx, cos_bx_error) = getCosDistance(b, x)
        val diff = cos_ax - cos_bx
        val error = cos_ax_error + cos_bx_error
        return if (diff > error) -1 else if (diff < -error) 1 else 0
    }


//    fun triageCompareCosDistances(x: R3VectorDecimal128, a: R3VectorDecimal128, b: R3VectorDecimal128): Int {
//        val (cos_ax, cos_ax_error) = getCosDistance(a, x)
//        val (cos_bx, cos_bx_error) = getCosDistance(b, x)
//        val diff = cos_ax - cos_bx
//        val error = cos_ax_error + cos_bx_error
//        return if (diff > error) -1 else if (diff < -error) 1 else 0
//    }


    fun triageCompareSin2Distances(x: S2Point, a: S2Point, b: S2Point): Int {
        val (sin2_ax, sin2_ax_error) = getSin2Distance(a, x)
        val (sin2_bx, sin2_bx_error) = getSin2Distance(b, x)
        val diff = sin2_ax - sin2_bx
        val error = sin2_ax_error + sin2_bx_error
        return if (diff > error) 1 else if (diff < -error) -1 else 0
    }


//    fun triageCompareSin2Distances(x: R3VectorDecimal128, a: R3VectorDecimal128, b: R3VectorDecimal128): Int {
//        val (sin2_ax, sin2_ax_error) = getSin2Distance(a, x)
//        logger.trace { " getSin2Distance(a, x) = ($sin2_ax, $sin2_ax_error)" }
//        val (sin2_bx, sin2_bx_error) = getSin2Distance(b, x)
//        logger.trace { " getSin2Distance(b, x) = ($sin2_bx, $sin2_bx_error)" }
//        val diff = sin2_ax - sin2_bx
//        val error = sin2_ax_error + sin2_bx_error
//        logger.trace { " diff = $diff, error = $error " }
//        return if (diff > error) 1 else if (diff < -error) -1 else 0
//    }


    fun exactCompareDistances(x: R3VectorExactFloat, a: R3VectorExactFloat, b: R3VectorExactFloat): Int {
        // This code produces the same result as though all points were reprojected
        // to lie exactly on the surface of the unit sphere.  It is based on testing
        // whether x.DotProd(a.Normalize()) < x.DotProd(b.Normalize()), reformulated
        // so that it can be evaluated using exact arithmetic.
        val cos_ax = x.dotProd(a)
        val cos_bx = x.dotProd(b)
        // If the two values have different signs, we need to handle that case now
        // before squaring them below.
        val a_sign = cos_ax.signum()
        val b_sign = cos_bx.signum()
        if (a_sign != b_sign) {
            return if (a_sign > b_sign) -1 else 1  // If cos(AX) > cos(BX), then AX < BX.
        }
        val cmp = ExactFloatType.minus(
            cos_bx.multiply(cos_bx, ExactFloatType.mathContext).multiply(a.norm2(), ExactFloatType.mathContext),
            cos_ax.multiply(cos_ax, ExactFloatType.mathContext).multiply(b.norm2(), ExactFloatType.mathContext)
        )
        return a_sign * cmp.signum()
    }

    // Given three points such that AX == BX (exactly), returns -1, 0, or +1
    // according whether AX < BX, AX == BX, or AX > BX after symbolic
    // perturbations are taken into account.

    internal fun symbolicCompareDistances(x: S2Point, a: S2Point, b: S2Point): Int {
        // Our symbolic perturbation strategy is based on the following model.
        // Similar to "simulation of simplicity", we assign a perturbation to every
        // point such that if A < B, then the symbolic perturbation for A is much,
        // much larger than the symbolic perturbation for B.  We imagine that
        // rather than projecting every point to lie exactly on the unit sphere,
        // instead each point is positioned on its own tiny pedestal that raises it
        // just off the surface of the unit sphere.  This means that the distance AX
        // is actually the true distance AX plus the (symbolic) heights of the
        // pedestals for A and X.  The pedestals are infinitesmally thin, so they do
        // not affect distance measurements except at the two endpoints.  If several
        // points project to exactly the same point on the unit sphere, we imagine
        // that they are placed on separate pedestals placed close together, where
        // the distance between pedestals is much, much less than the height of any
        // pedestal.  (There are a finite number of S2Points, and therefore a finite
        // number of pedestals, so this is possible.)
        //
        // If A < B, then A is on a higher pedestal than B, and therefore AX > BX.
        return if (a < b) 1 else if (a > b) -1 else 0
    }


    fun compareSin2Distances(x: S2Point, a: S2Point, b: S2Point): Int {
        val sign = triageCompareSin2Distances(x, a, b)
        logger.trace { "triageCompareSin2Distances(x, a, b) = $sign" }
//        if (sign != 0)
        return sign
//        sign = triageCompareSin2Distances(x.toDecimal128(), a.toDecimal128(), b.toDecimal128())
//        logger.trace { "triageCompareSin2Distances(x.toDecimal128(), a.toDecimal128(), b.toDecimal128()) = $sign" }
//        return sign
    }


    internal fun triageCompareCosDistance(x: S2Point, y: S2Point, r2: Double): Int {
        val (cos_xy, cos_xy_error) = getCosDistance(x, y)
        val cos_r = 1 - 0.5 * r2
        val cos_r_error = 2 * DoubleType.roundingEpsilon * cos_r
        val diff = cos_xy - cos_r
        val error = cos_xy_error + cos_r_error
        return if (diff > error) -1 else if (diff < -error) 1 else 0
    }


    internal fun triageCompareSin2Distance(x: S2Point, y: S2Point, r2: Double): Int {
        assert(r2 < 2.0)  // Only valid for distance limits < 90 degrees.

        val (sin2_xy, sin2_xy_error) = getSin2Distance(x, y)
        val sin2_r = r2 * (1 - 0.25 * r2)
        val sin2_r_error = 3 * DoubleType.roundingEpsilon * sin2_r
        val diff = sin2_xy - sin2_r
        val error = sin2_xy_error + sin2_r_error
        return if (diff > error) 1 else if (diff < -error) -1 else 0
    }


    internal fun exactCompareDistance(x: R3VectorExactFloat, y: R3VectorExactFloat, r2: BigDecimal): Int {
        // This code produces the same result as though all points were reprojected
        // to lie exactly on the surface of the unit sphere.  It is based on
        // comparing the cosine of the angle XY (when both points are projected to
        // lie exactly on the sphere) to the given threshold.
        val cosXy = x.dotProd(y)
        val cosR = ExactFloatType.cast(1.0) - ExactFloatType.cast(0.5) * r2
        // If the two values have different signs, we need to handle that case now
        // before squaring them below.
        val xySign = cosXy.signum()
        val rSign = cosR.signum()
        if (xySign != rSign) {
            return if (xySign > rSign) -1 else 1  // If cos(XY) > cos(r), then XY < r.
        }
        val cmp = cosR * cosR * x.norm2() * y.norm2() - cosXy * cosXy
        return xySign * cmp.signum()
    }

    // Helper function that compares the distance XY against the squared chord
    // distance "r2" using the given precision "T".

    private fun triageCompareDistance(x: S2Point, y: S2Point, r2: Double): Int {
        // The Sin2 method is much more accurate for small distances, but it is only
        // valid when the actual distance and the distance limit are both less than
        // 90 degrees.  So we always start with the Cos method.
        var sign = triageCompareCosDistance(x, y, r2)
        if (sign == 0 && r2 < k45Degrees.length2) {
            sign = triageCompareSin2Distance(x, y, r2)
        }
        return sign
    }

    // Helper function that returns "a0" or "a1", whichever is closer to "x".
    // Also returns the squared distance from the returned point to "x" in "ax2".

    private fun getClosestVertex(x: S2Point, a0: S2Point, a1: S2Point): Pair<S2Point, Double> {
        val a0x2 = (a0 - x).norm2()
        val a1x2 = (a1 - x).norm2()
        return if (a0x2 < a1x2 || (a0x2 == a1x2 && a0 < a1)) {
            a0 to a0x2
        } else {
            a1 to a1x2
        }
    }

    // Helper function that returns -1, 0, or +1 according to whether the distance
    // from "x" to the great circle through (a0, a1) is less than, equal to, or
    // greater than the given squared chord length "r2".  This method computes the
    // squared sines of the distances involved, which is more accurate when the
    // distances are small (less than 45 degrees).
    //
    // The remaining parameters are functions of (a0, a1) and are passed in
    // because they have already been computed: n = (a0 - a1) x (a0 + a1),
    // n1 = n.Norm(), and n2 = n.Norm2().

    private fun triageCompareLineSin2Distance(
        x: S2Point,
        a0: S2Point,
        a1: S2Point,
        r2: Double,
        n: S2Point,
        n1: Double,
        n2: Double
    ): Int {
        //constexpr T T_ERR = rounding_epsilon<T>()

        // The minimum distance is to a point on the edge interior.  Since the true
        // distance to the edge is always less than 90 degrees, we can return
        // immediately if the limit is 90 degrees or larger.
        if (r2 >= 2.0) return -1  // distance < limit

        // Otherwise we compute sin^2(distance to edge) to get the best accuracy
        // when the distance limit is small (e.g., S2::kIntersectionError).
        val n2sin2_r = n2 * r2 * (1 - 0.25 * r2)
        var n2sin2_r_error = 6 * DoubleType.roundingEpsilon * n2sin2_r
        val (cv, ax2) = getClosestVertex(x, a0, a1)
        val xDn = (x - cv).dotProd(n)
        val xDn2 = xDn * xDn
        val c1 =
            (((3.5 + 2 * M_SQRT3) * n1 + 32 * M_SQRT3 * DoubleType.roundingEpsilon) * DoubleType.roundingEpsilon * sqrt(
                ax2
            ))
        val xDn2_error = 4 * DoubleType.roundingEpsilon * xDn2 + (2 * abs(xDn) + c1) * c1

        n2sin2_r_error += 8 * DoubleType.roundingEpsilon * n2sin2_r
        val diff = xDn2 - n2sin2_r
        val error = xDn2_error + n2sin2_r_error
        return if (diff > error) 1 else if (diff < -error) -1 else 0
    }

    // Like TriageCompareLineSin2Distance, but this method computes the squared
    // cosines of the distances involved.  It is more accurate when the distances
    // are large (greater than 45 degrees).

    private fun triageCompareLineCos2Distance(
        x: S2Point,
        a0: S2Point,
        a1: S2Point,
        r2: Double,
        n: S2Point,
        n1: Double,
        n2: Double
    ): Int {
        // The minimum distance is to a point on the edge interior.  Since the true
        // distance to the edge is always less than 90 degrees, we can return
        // immediately if the limit is 90 degrees or larger.
        if (r2 >= 2.0) return -1  // distance < limit

        // Otherwise we compute cos^2(distance to edge).
        val cos_r = 1 - 0.5 * r2
        val n2cos2_r = n2 * cos_r * cos_r
        var n2cos2_r_error = 7 * DoubleType.roundingEpsilon * n2cos2_r

        // The length of M = X.CrossProd(N) is the cosine of the distance.
        val m2 = x.crossProd(n).norm2()
        val m1 = sqrt(m2)
        val m1_error = ((1 + 8 / M_SQRT3) * n1 + 32 * M_SQRT3 * DoubleType.roundingEpsilon) * DoubleType.roundingEpsilon
        val m2_error = 3 * DoubleType.roundingEpsilon * m2 + (2 * m1 + m1_error) * m1_error

        n2cos2_r_error += 8 * DoubleType.roundingEpsilon * n2cos2_r

        val diff = m2 - n2cos2_r
        val error = m2_error + n2cos2_r_error
        return if (diff > error) -1 else if (diff < -error) 1 else 0
    }


    private fun triageCompareLineDistance(
        x: S2Point,
        a0: S2Point,
        a1: S2Point,
        r2: Double,
        n: S2Point,
        n1: Double,
        n2: Double
    ): Int =
        if (r2 < k45Degrees.length2) {
            triageCompareLineSin2Distance(x, a0, a1, r2, n, n1, n2)
        } else {
            triageCompareLineCos2Distance(x, a0, a1, r2, n, n1, n2)
        }


    internal fun triageCompareEdgeDistance(x: S2Point, a0: S2Point, a1: S2Point, r2: Double): Int {
        // First we need to decide whether the closest point is an edge endpoint or
        // somewhere in the interior.  To determine this we compute a plane
        // perpendicular to (a0, a1) that passes through X.  Letting M be the normal
        // to this plane, the closest point is in the edge interior if and only if
        // a0.M < 0 and a1.M > 0.  Note that we can use "<" rather than "<=" because
        // if a0.M or a1.M is zero exactly then it doesn't matter which code path we
        // follow (since the distance to an endpoint and the distance to the edge
        // interior are exactly the same in this case).
        val n = (a0 - a1).crossProd(a0 + a1)
        val m = n.crossProd(x)
        // For better accuracy when the edge (a0,a1) is very short, we subtract "x"
        // before computing the dot products with M.
        val a0_dir = a0 - x
        val a1_dir = a1 - x
        val a0_sign = a0_dir.dotProd(m)
        val a1_sign = a1_dir.dotProd(m)
        val n2 = n.norm2()
        val n1 = sqrt(n2)
        val n1_error =
            ((3.5 + 8 / M_SQRT3) * n1 + 32 * M_SQRT3 * DoubleType.roundingEpsilon) * DoubleType.roundingEpsilon
        val a0_sign_error = n1_error * a0_dir.norm()
        val a1_sign_error = n1_error * a1_dir.norm()
        return if (abs(a0_sign) < a0_sign_error || abs(a1_sign) < a1_sign_error) {
            // It is uncertain whether minimum distance is to an edge vertex or to the
            // edge interior.  We handle this by computing both distances and checking
            // whether they yield the same result.
            val vertex_sign = min(triageCompareDistance(x, a0, r2), triageCompareDistance(x, a1, r2))
            val line_sign = triageCompareLineDistance(x, a0, a1, r2, n, n1, n2)
            if (vertex_sign == line_sign) line_sign else 0
        } else if (a0_sign >= 0 || a1_sign <= 0) {
            // The minimum distance is to an edge endpoint.
            min(triageCompareDistance(x, a0, r2), triageCompareDistance(x, a1, r2))
        } else {
            // The minimum distance is to the edge interior.
            triageCompareLineDistance(x, a0, a1, r2, n, n1, n2)
        }
    }

    // REQUIRES: the closest point to "x" is in the interior of edge (a0, a1).

    fun exactCompareLineDistance(
        x: R3VectorExactFloat,
        a0: R3VectorExactFloat,
        a1: R3VectorExactFloat,
        r2: BigDecimal
    ): Int {
        // Since we are given that the closest point is in the edge interior, the
        // true distance is always less than 90 degrees (which corresponds to a
        // squared chord length of 2.0).
        if (r2 >= ExactFloatType.cast(2.0)) return -1  // distance < limit

        // Otherwise compute the edge normal
        val n = a0.crossProd(a1)
        val sinD = x.dotProd(n)
        val sin2R = r2 * (ExactFloatType.cast(1.0) - ExactFloatType.cast(0.25) * r2)
        val cmp = sinD * sinD - sin2R * x.norm2() * n.norm2()
        return cmp.signum()
    }


    fun exactCompareEdgeDistance(x: S2Point, a0: S2Point, a1: S2Point, r: S1ChordAngle): Int {
        // Even if previous calculations were uncertain, we might not need to do
        // *all* the calculations in exact arithmetic here.  For example it may be
        // easy to determine whether "x" is closer to an endpoint or the edge
        // interior.  The only calculation where we always use exact arithmetic is
        // when measuring the distance to the extended line (great circle) through
        // "a0" and "a1", since it is virtually certain that the previous floating
        // point calculations failed in that case.
        //
        // CompareEdgeDirections also checks that no edge has antipodal endpoints.
        if (compareEdgeDirections(a0, a1, a0, x) > 0 && compareEdgeDirections(a0, a1, x, a1) > 0) {
            // The closest point to "x" is along the interior of the edge.
            return exactCompareLineDistance(
                x.toExactFloat(),
                a0.toExactFloat(),
                a1.toExactFloat(),
                ExactFloatType.cast(r.length2)
            )
        } else {
            // The closest point to "x" is one of the edge endpoints.
            return min(compareDistance(x, a0, r), compareDistance(x, a1, r))
        }
    }


    internal fun triageCompareEdgeDirections(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): Int {
        val na = (a0 - a1).crossProd(a0 + a1)
        val nb = (b0 - b1).crossProd(b0 + b1)
        val na_len = na.norm()
        val nb_len = nb.norm()
        val cos_ab = na.dotProd(nb)
        val cos_ab_error =
            ((5 + 4 * M_SQRT3) * na_len * nb_len + 32 * M_SQRT3 * DoubleType.roundingEpsilon * (na_len + nb_len)) * DoubleType.roundingEpsilon
        return if (cos_ab > cos_ab_error) 1 else if (cos_ab < -cos_ab_error) -1 else 0
    }


    fun arePointsLinearlyDependent(x: R3VectorExactFloat, y: R3VectorExactFloat): Boolean {
        val n = x.crossProd(y)
        return n[0].signum() == 0 && n[1].signum() == 0 && n[2].signum() == 0
    }


    fun arePointsAntipodal(x: R3VectorExactFloat, y: R3VectorExactFloat): Boolean {
        return arePointsLinearlyDependent(x, y) && x.dotProd(y).signum() < 0
    }


    internal fun exactCompareEdgeDirections(
        a0: R3VectorExactFloat,
        a1: R3VectorExactFloat,
        b0: R3VectorExactFloat,
        b1: R3VectorExactFloat
    ): Int {
        assert(!arePointsAntipodal(a0, a1))
        assert(!arePointsAntipodal(b0, b1))
        return a0.crossProd(a1).dotProd(b0.crossProd(b1)).signum()
    }

    // If triangle ABC has positive sign, returns its circumcenter.  If ABC has
    // negative sign, returns the negated circumcenter.

    private fun getCircumcenter(a: S2Point, b: S2Point, c: S2Point): Pair<S2Point, Double> {
        // We compute the circumcenter using the intersection of the perpendicular
        // bisectors of AB and BC.  The formula is essentially
        //
        //    Z = ((A x B) x (A + B)) x ((B x C) x (B + C)),
        //
        // except that we compute the cross product (A x B) as (A - B) x (A + B)
        // (and similarly for B x C) since this is much more stable when the inputs
        // are unit vectors.
        val ab_diff = a - b
        val ab_sum = a + b
        val bc_diff = b - c
        val bc_sum = b + c
        val nab = ab_diff.crossProd(ab_sum)
        val nab_len = nab.norm()
        val ab_len = ab_diff.norm()
        val nbc = bc_diff.crossProd(bc_sum)
        val nbc_len = nbc.norm()
        val bc_len = bc_diff.norm()
        val mab = nab.crossProd(ab_sum)
        val mbc = nbc.crossProd(bc_sum)
        val error = (((16 + 24 * M_SQRT3) * DoubleType.roundingEpsilon +
                8 * DoubleType.roundingEpsilon * (ab_len + bc_len)) * nab_len * nbc_len +
                128 * M_SQRT3 * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * (nab_len + nbc_len) +
                3 * 4096 * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon * DoubleType.roundingEpsilon)
        return mab.crossProd(mbc) to error
    }


    internal fun triageEdgeCircumcenterSign(
        x0: S2Point,
        x1: S2Point,
        a: S2Point,
        b: S2Point,
        c: S2Point,
        abc_sign: Int
    ): Int {
        val T_ERR = DoubleType.roundingEpsilon

        // Compute the circumcenter Z of triangle ABC, and then test which side of
        // edge X it lies on.
        val (z, z_error) = getCircumcenter(a, b, c)
        val nx = (x0 - x1).crossProd(x0 + x1)
        // If the sign of triangle ABC is negative, then we have computed -Z and the
        // result should be negated.
        val result = abc_sign * nx.dotProd(z)

        val z_len = z.norm()
        val nx_len = nx.norm()
        val nx_error = ((1 + 2 * M_SQRT3) * nx_len + 32 * M_SQRT3 * DoubleType.roundingEpsilon) * T_ERR
        val result_error = ((3 * T_ERR * nx_len + nx_error) * z_len + z_error * nx_len)
        return if (result > result_error) 1 else if (result < -result_error) -1 else 0
    }


    fun exactEdgeCircumcenterSign(
        x0: R3VectorExactFloat,
        x1: R3VectorExactFloat,
        a: R3VectorExactFloat,
        b: R3VectorExactFloat,
        c: R3VectorExactFloat,
        abc_sign: Int
    ): Int {
        // Return zero if the edge X is degenerate.  (Also see the comments in
        // SymbolicEdgeCircumcenterSign.)
        if (arePointsLinearlyDependent(x0, x1)) {
            assert(x0.dotProd(x1) > BigDecimal.ZERO);  // Antipodal edges not allowed.
            return 0
        }
        // The simplest predicate for testing whether the sign is positive is
        //
        // (1)  (X0 x X1) . (|C|(A x B) + |A|(B x C) + |B|(C x A)) > 0
        //
        // where |A| denotes A.Norm() and the expression after the "." represents
        // the circumcenter of triangle ABC.  (This predicate is terrible from a
        // numerical accuracy point of view, but that doesn't matter since we are
        // going to use exact arithmetic.)  This predicate also assumes that
        // triangle ABC is CCW (positive sign); we correct for that below.
        //
        // The only problem with evaluating this inequality is that computing |A|,
        // |B| and |C| requires square roots.  To avoid this problem we use the
        // standard technique of rearranging the inequality to isolate at least one
        // square root and then squaring both sides.  We need to repeat this process
        // twice in order to eliminate all the square roots, which leads to a
        // polynomial predicate of degree 20 in the input arguments.
        //
        // Rearranging (1) we get
        //
        //      (X0 x X1) . (|C|(A x B) + |A|(B x C)) > |B|(X0 x X1) . (A x C)
        //
        // Before squaring we need to check the sign of each side.  If the signs are
        // different then we know the result without squaring, and if the signs are
        // both negative then after squaring both sides we need to invert the
        // result.  Define
        //
        //      dAB = (X0 x X1) . (A x B)
        //      dBC = (X0 x X1) . (B x C)
        //      dCA = (X0 x X1) . (C x A)
        //
        // Then we can now write the inequality above as
        //
        // (2)  |C| dAB + |A| dBC > -|B| dCA
        //
        // The RHS of (2) is positive if dCA < 0, and the LHS of (2) is positive if
        // (|C| dAB + |A| dBC) > 0.  Since the LHS has square roots, we need to
        // eliminate them using the same process.  Rewriting the LHS as
        //
        // (3)  |C| dAB > -|A| dBC
        //
        // we again need to check the signs of both sides.  Let's start with that.
        // We also precompute the following values because they are used repeatedly
        // when squaring various expressions below:
        //
        //     abc2 = |A|^2 dBC^2
        //     bca2 = |B|^2 dCA^2
        //     cab2 = |C|^2 dAB^2
        val nx = x0.crossProd(x1)
        val dab = nx.dotProd(a.crossProd(b))
        val dbc = nx.dotProd(b.crossProd(c))
        val dca = nx.dotProd(c.crossProd(a))
        val abc2 = a.norm2() * (dbc * dbc)
        val bca2 = b.norm2() * (dca * dca)
        val cab2 = c.norm2() * (dab * dab)

        // If the two sides of (3) have different signs (including the case where
        // one side is zero) then we know the result.  Also, if both sides are zero
        // then we know the result.  The following logic encodes this.
        val lhs3_sgn = dab.signum()
        val rhs3_sgn = -dbc.signum()
        var lhs2_sgn = max(-1, min(1, lhs3_sgn - rhs3_sgn))
        if (lhs2_sgn == 0 && lhs3_sgn != 0) {
            // Both sides of (3) have the same non-zero sign, so square both sides.
            // If both sides were negative then invert the result.
            lhs2_sgn = (cab2 - abc2).signum() * lhs3_sgn
        }
        // Now if the two sides of (2) have different signs then we know the result
        // of this entire function.
        val rhs2_sgn = -dca.signum()
        var result = max(-1, min(1, lhs2_sgn - rhs2_sgn))
        if (result == 0 && lhs2_sgn != 0) {
            // Both sides of (2) have the same non-zero sign, so square both sides.
            // (If both sides were negative then we invert the result below.)
            // This gives
            //
            //        |C|^2 dAB^2 + |A|^2 dBC^2 + 2 |A| |C| dAB dBC > |B|^2 dCA^2
            //
            // This expression still has square roots (|A| and |C|), so we rewrite as
            //
            // (4)    2 |A| |C| dAB dBC > |B|^2 dCA^2 - |C|^2 dAB^2 - |A|^2 dBC^2 .
            //
            // Again, if the two sides have different signs then we know the result.
            val lhs4_sgn = dab.signum() * dbc.signum()
            val rhs4 = bca2 - cab2 - abc2
            result = max(-1, min(1, lhs4_sgn - rhs4.signum()))
            if (result == 0 && lhs4_sgn != 0) {
                // Both sides of (4) have the same non-zero sign, so square both sides.
                // If both sides were negative then invert the result.
                result = (ExactFloatType.cast(4.0) * abc2 * cab2 - rhs4 * rhs4).signum() * lhs4_sgn
            }
            // Correct the sign if both sides of (2) were negative.
            result *= lhs2_sgn
        }
        // If the sign of triangle ABC is negative, then we have computed -Z and the
        // result should be negated.
        return abc_sign * result
    }

    // Like Sign, except this method does not use symbolic perturbations when
    // the input points are exactly coplanar with the origin (i.e., linearly
    // dependent).  Clients should never use this method, but it is useful here in
    // order to implement the combined pedestal/axis-aligned perturbation scheme
    // used by some methods (such as EdgeCircumcenterSign).

    fun unperturbedSign(a: S2Point, b: S2Point, c: S2Point): Int {
        var sign = triageSign(a, b, c, a.crossProd(b))
        if (sign == 0) sign = expensiveSign(a, b, c, false /*perturb*/)
        return sign
    }

    // Given arguments such that ExactEdgeCircumcenterSign(x0, x1, a, b, c) == 0,
    // returns the value of Sign(X0, X1, Z) (where Z is the circumcenter of
    // triangle ABC) after symbolic perturbations are taken into account.  The
    // result is zero only if X0 == X1, A == B, B == C, or C == A.  (It is nonzero
    // if these pairs are exactly proportional to each other but not equal.)

    fun symbolicEdgeCircumcenterSign(x0: S2Point, x1: S2Point, a_arg: S2Point, b_arg: S2Point, c_arg: S2Point): Int {
        // We use the same perturbation strategy as SymbolicCompareDistances.  Note
        // that pedestal perturbations of X0 and X1 do not affect the result,
        // because Sign(X0, X1, Z) does not change when its arguments are scaled
        // by a positive factor.  Therefore we only need to consider A, B, C.
        // Suppose that A is the smallest lexicographically and therefore has the
        // largest perturbation.  This has the effect of perturbing the circumcenter
        // of ABC slightly towards A, and since the circumcenter Z was previously
        // exactly collinear with edge X, this implies that after the perturbation
        // Sign(X0, X1, Z) == UnperturbedSign(X0, X1, A).  (We want the result
        // to be zero if X0, X1, and A are linearly dependent, rather than using
        // symbolic perturbations, because these perturbations are defined to be
        // much, much smaller than the pedestal perturbation of B and C that are
        // considered below.)
        //
        // If A is also exactly collinear with edge X, then we move on to the next
        // smallest point lexicographically out of {B, C}.  It is easy to see that
        // as long as A, B, C are all distinct, one of these three Sign calls
        // will be nonzero, because if A, B, C are all distinct and collinear with
        // edge X then their circumcenter Z coincides with the normal of X, and
        // therefore Sign(X0, X1, Z) is nonzero.
        //
        // This function could be extended to handle the case where X0 and X1 are
        // linearly dependent as follows.  First, suppose that every point has both
        // a pedestal peturbation as described above, and also the three
        // axis-aligned perturbations described in the "Simulation of Simplicity"
        // paper, where all pedestal perturbations are defined to be much, much
        // larger than any axis-aligned perturbation.  Note that since pedestal
        // perturbations have no effect on Sign, we can use this model for *all*
        // the S2 predicates, which ensures that all the various predicates are
        // fully consistent with each other.
        //
        // With this model, the strategy described above yields the correct result
        // unless X0 and X1 are exactly linearly dependent.  When that happens, then
        // no perturbation (pedestal or axis-aligned) of A,B,C affects the result,
        // and no pedestal perturbation of X0 or X1 affects the result, therefore we
        // need to consider the smallest axis-aligned perturbation of X0 or X1.  The
        // first perturbation that makes X0 and X1 linearly independent yields the
        // result.  Supposing that X0 < X1, this is the perturbation of X0[2] unless
        // both points are multiples of [0, 0, 1], in which case it is the
        // perturbation of X0[1].  The sign test can be implemented by computing the
        // perturbed cross product of X0 and X1 and taking the dot product with the
        // exact value of Z.  For example if X0[2] is perturbed, the perturbed cross
        // product is proportional to (0, 0, 1) x X1 = (-X1[1], x1[0], 0).  Note
        // that if the dot product with Z is exactly zero, then it is still
        // necessary to fall back to pedestal perturbations of A, B, C, but one of
        // these perturbations is now guaranteed to succeed.

        // If any two triangle vertices are equal, the result is zero.
        if (a_arg == b_arg || b_arg == c_arg || c_arg == a_arg) return 0

        // Sort A, B, C in lexicographic order.
        var a = a_arg
        var b = b_arg
        var c = c_arg
        if (b < a) {
            val temp = a; a = b; b = temp; }
        if (c < b) {
            val temp = b; b = c; c = temp; }
        if (b < a) {
            val temp = a; a = b; b = temp; }

        // Now consider the perturbations in decreasing order of size.
        var sign = unperturbedSign(x0, x1, a)
        if (sign != 0) return sign
        sign = unperturbedSign(x0, x1, b)
        if (sign != 0) return sign
        return unperturbedSign(x0, x1, c)
    }


    fun triageVoronoiSiteExclusion(a: S2Point, b: S2Point, x0: S2Point, x1: S2Point, r2: Double): Excluded {
        val T_ERR = DoubleType.roundingEpsilon

        // Define the "coverage disc" of a site S to be the disc centered at S with
        // radius r (i.e., squared chord angle length r2).  Similarly, define the
        // "coverage interval" of S along an edge X to be the intersection of X with
        // the coverage disc of S.  The coverage interval can be represented as the
        // point at the center of the interval and an angle that measures the
        // semi-width or "radius" of the interval.
        //
        // To test whether site A excludes site B along the input edge X, we test
        // whether the coverage interval of A contains the coverage interval of B.
        // Let "ra" and "rb" be the radii (semi-widths) of the two intervals, and
        // let "d" be the angle between their center points.  Then "a" properly
        // contains "b" if (ra - rb > d), and "b" contains "a" if (rb - ra > d).
        // Note that only one of these conditions can be true.  Therefore we can
        // determine whether one site excludes the other by checking whether
        //
        // (1)   |rb - ra| > d
        //
        // and use the sign of (rb - ra) to determine which site is excluded.
        //
        // The actual code is based on the following.  Let A1 and B1 be the unit
        // vectors A and B scaled by cos(r) (these points are inside the sphere).
        // The planes perpendicular to OA1 and OA2 cut off two discs of radius r
        // around A and B.  Now consider the two lines (inside the sphere) where
        // these planes intersect the plane containing the input edge X, and let A2
        // and B2 be the points on these lines that are closest to A and B.  The
        // coverage intervals of A and B can be represented as an interval along
        // each of these lines, centered at A2 and B2.  Let P1 and P2 be the
        // endpoints of the coverage interval for A, and let Q1 and Q2 be the
        // endpoints of the coverage interval for B.  We can view each coverage
        // interval as either a chord through the sphere's interior, or as a segment
        // of the original edge X (by projecting the chord onto the sphere's
        // surface).
        //
        // To check whether B's interval is contained by A's interval, we test
        // whether both endpoints of B's interval (Q1 and Q2) are contained by A's
        // interval.  E.g., we could test whether Qi.DotProd(A2) > A2.Norm2().
        //
        // However rather than constructing the actual points A1, A2, and so on, it
        // turns out to be more efficient to compute the sines and cosines
        // ("components") of the various angles and then use trigonometric
        // identities.  Predicate (1) can be expressed as
        //
        //      |sin(rb - ra)| > sin(d)
        //
        // provided that |d| <= Pi/2 (which must be checked), and then expanded to
        //
        // (2)  |sin(rb) cos(ra) - sin(ra) cos(rb)| > sin(d) .
        //
        // The components of the various angles can be expressed using dot and cross
        // products based on the construction above:
        //
        //   sin(ra) = sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2) / |aXn|
        //   cos(ra) = cos(r) |a| |n| / |aXn|
        //   sin(rb) = sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) / |bXn|
        //   cos(rb) = cos(r) |b| |n| / |bXn|
        //   sin(d)  = (aXb).n |n| / (|aXn| |bXn|)
        //   cos(d)  = (aXn).(bXn) / (|aXn| |bXn|)
        //
        // Also, the squared chord length r2 is equal to 4 * sin^2(r / 2), which
        // yields the following relationships:
        //
        //   sin(r)  = sqrt(r2 (1 - r2 / 4))
        //   cos(r)  = 1 - r2 / 2
        //
        // We then scale both sides of (2) by |aXn| |bXn| / |n| (in order to
        // minimize the number of calculations and to avoid divisions), which gives:
        //
        //    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
        //            |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
        //
        // Furthermore we can substitute |a| = |b| = 1 (as long as this is taken
        // into account in the error bounds), yielding
        //
        // (3)   cos(r) |sqrt(sin^2(r) |n|^2 - |b.n|^2) -
        //               sqrt(sin^2(r) |n|^2 - |a.n|^2)| > (aXb).n
        //
        // The code below is more complicated than this because many expressions
        // have been modified for better numerical stability.  For example, dot
        // products between unit vectors are computed using (x - y).DotProd(x + y),
        // and the dot product between a point P and the normal N of an edge X is
        // measured using (P - Xi).DotProd(N) where Xi is the endpoint of X that is
        // closer to P.

        val n = (x0 - x1).crossProd(x0 + x1)  // 2 * x0.CrossProd(x1)
        val n2 = n.norm2()
        val n1 = sqrt(n2)
        // This factor is used in the error terms of dot products with "n" below.
        val Dn_error = ((3.5 + 2 * M_SQRT3) * n1 + 32 * M_SQRT3 * DoubleType.roundingEpsilon) * T_ERR

        val cos_r = 1 - 0.5 * r2
        val sin2_r = r2 * (1 - 0.25 * r2)
        val n2sin2_r = n2 * sin2_r

        // "ra" and "rb" denote sin(ra) and sin(rb) after the scaling above.
        val (acv, ax2) = getClosestVertex(a, x0, x1)
        val aDn = (a - acv).dotProd(n)
        val aDn2 = aDn * aDn
        val aDn_error = Dn_error * sqrt(ax2)
        val ra2 = n2sin2_r - aDn2
        val ra2_error =
            (8 * DoubleType.roundingEpsilon + 4 * T_ERR) * aDn2 + (2 * abs(aDn) + aDn_error) * aDn_error + 6 * T_ERR * n2sin2_r
        // This is the minimum possible value of ra2, which is used to bound the
        // derivative of sqrt(ra2) in computing ra_error below.
        val min_ra2 = ra2 - ra2_error
        if (min_ra2 < 0) return Excluded.UNCERTAIN
        val ra = sqrt(ra2)
        // Includes the ra2 subtraction error above.
        val ra_error = 1.5 * T_ERR * ra + 0.5 * ra2_error / sqrt(min_ra2)

        val (bcv, bx2) = getClosestVertex(b, x0, x1)
        val bDn = (b - bcv).dotProd(n)
        val bDn2 = bDn * bDn
        val bDn_error = Dn_error * sqrt(bx2)
        val rb2 = n2sin2_r - bDn2
        val rb2_error =
            (8 * DoubleType.roundingEpsilon + 4 * T_ERR) * bDn2 + (2 * abs(bDn) + bDn_error) * bDn_error + 6 * T_ERR * n2sin2_r
        val min_rb2 = rb2 - rb2_error
        if (min_rb2 < 0) return Excluded.UNCERTAIN
        val rb = sqrt(rb2)
        // Includes the rb2 subtraction error above.
        val rb_error = 1.5 * T_ERR * rb + 0.5 * rb2_error / sqrt(min_rb2)

        // The sign of LHS(3) determines which site may be excluded by the other.
        val lhs3 = cos_r * (rb - ra)
        val abs_lhs3 = abs(lhs3)
        val lhs3_error = cos_r * (ra_error + rb_error) + 3 * T_ERR * abs_lhs3

        // Now we evaluate the RHS of (3), which is proportional to sin(d).
        val aXb = (a - b).crossProd(a + b)  // 2 * a.CrossProd(b)
        val aXb1 = aXb.norm()
        val sin_d = 0.5 * aXb.dotProd(n)
        val sin_d_error = (4 * DoubleType.roundingEpsilon + (2.5 + 2 * M_SQRT3) * T_ERR) * aXb1 * n1 +
                16 * M_SQRT3 * DoubleType.roundingEpsilon * T_ERR * (aXb1 + n1)

        // If LHS(3) is definitely less than RHS(3), neither site excludes the other.
        val result = abs_lhs3 - sin_d
        val result_error = lhs3_error + sin_d_error
        if (result < -result_error) return Excluded.NEITHER

        // Otherwise, before proceeding further we need to check that |d| <= Pi/2.
        // In fact, |d| < Pi/2 is enough because of the requirement that r < Pi/2.
        // The following expression represents cos(d) after scaling; it is
        // equivalent to (aXn).(bXn) but has about 30% less error.
        val cos_d = a.dotProd(b) * n2 - aDn * bDn
        val cos_d_error = ((8 * DoubleType.roundingEpsilon + 5 * T_ERR) * abs(aDn) + aDn_error) * abs(bDn) +
                (abs(aDn) + aDn_error) * bDn_error + (8 * DoubleType.roundingEpsilon + 8 * T_ERR) * n2
        if (cos_d <= -cos_d_error) return Excluded.NEITHER

        // Potential optimization: if the sign of cos(d) is uncertain, then instead
        // we could check whether cos(d) >= cos(r).  Unfortunately this is fairly
        // expensive since it requires computing denominator |aXn||bXn| of cos(d)
        // and the associated error bounds.  In any case this case is relatively
        // rare so it seems better to punt.
        if (cos_d < cos_d_error) return Excluded.UNCERTAIN

        // Normally we have d > 0 because the sites are sorted so that A is closer
        // to X0 and B is closer to X1.  However if the edge X is longer than Pi/2,
        // and the sites A and B are beyond its endpoints, then AB can wrap around
        // the sphere in the opposite direction from X.  In this situation d < 0 but
        // each site is closest to one endpoint of X, so neither excludes the other.
        //
        // It turns out that this can happen only when the site that is further away
        // from edge X is less than 90 degrees away from whichever endpoint of X it
        // is closer to.  It is provable that if this distance is less than 90
        // degrees, then it is also less than r2, and therefore the Voronoi regions
        // of both sites intersect the edge.
        if (sin_d < -sin_d_error) {
            val r90 = S1ChordAngle.right().length2
            // "ca" is negative if Voronoi region A definitely intersects edge X.
            val ca = if (lhs3 < -lhs3_error) -1 else triageCompareCosDistance(a, x0, r90)
            val cb = if (lhs3 > lhs3_error) -1 else triageCompareCosDistance(b, x1, r90)
            if (ca < 0 && cb < 0) return Excluded.NEITHER
            if (ca <= 0 && cb <= 0) return Excluded.UNCERTAIN
            if (abs_lhs3 <= lhs3_error) return Excluded.UNCERTAIN
        } else if (sin_d <= sin_d_error) {
            return Excluded.UNCERTAIN
        }
        // Now we can finish checking the results of predicate (3).
        if (result <= result_error) return Excluded.UNCERTAIN
        assert(abs_lhs3 > lhs3_error)
        return if (lhs3 > 0) Excluded.FIRST else Excluded.SECOND
    }


    fun exactVoronoiSiteExclusion(
        a: R3VectorExactFloat,
        b: R3VectorExactFloat,
        x0: R3VectorExactFloat,
        x1: R3VectorExactFloat,
        r2: BigDecimal
    ): Excluded {
        assert(!arePointsAntipodal(x0, x1))

        // Recall that one site excludes the other if
        //
        // (1)  |sin(rb - ra)| > sin(d)
        //
        // and that the sign of (rb - ra) determines which site is excluded (see the
        // comments in TriageVoronoiSiteExclusion).  To evaluate this using exact
        // arithmetic, we expand this to the same predicate as before:
        //
        // (2)    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
        //                |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
        //
        // We also need to verify that d <= Pi/2, which is implemented by checking
        // that sin(d) >= 0 and cos(d) >= 0.
        //
        // To eliminate the square roots we use the standard technique of
        // rearranging the inequality to isolate at least one square root and then
        // squaring both sides.  We need to repeat this process twice in order to
        // eliminate all the square roots, which leads to a polynomial predicate of
        // degree 20 in the input arguments (i.e., degree 4 in each of "a", "b",
        // "x0", "x1", and "r2").
        //
        // Before squaring we need to check the sign of each side.  We also check
        // the condition that cos(d) >= 0.  Given what else we need to compute, it
        // is cheaper use the identity (aXn).(bXn) = (a.b) |n|^2 - (a.n)(b.n) .
        val n = x0.crossProd(x1)
        val n2 = n.norm2()
        val aDn = a.dotProd(n)
        val bDn = b.dotProd(n)
        val cos_d = a.dotProd(b) * n2 - aDn * bDn
        if (cos_d.signum() < 0) return Excluded.NEITHER

        // Otherwise we continue evaluating the LHS of (2), defining
        //    sa = |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)
        //    sb = |a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) .
        // The sign of the LHS of (2) (before taking the absolute value) determines
        // which coverage interval is larger and therefore which site is potentially
        // being excluded.
        val a2 = a.norm2()
        val b2 = b.norm2()
        val n2sin2_r = r2 * (ExactFloatType.cast(1.0) - ExactFloatType.cast(0.25) * r2) * n2
        val sa2 = b2 * (n2sin2_r * a2 - aDn * aDn)
        val sb2 = a2 * (n2sin2_r * b2 - bDn * bDn)
        val lhs2_sgn = (sb2 - sa2).signum()

        // If the RHS of (2) is negative (corresponding to sin(d) < 0), then we need
        // to consider the possibility that the edge AB wraps around the sphere in
        // the opposite direction from edge X, with the result that neither site
        // excludes the other (see TriageVoronoiSiteExclusion).
        val rhs2 = a.crossProd(b).dotProd(n)
        val rhs2_sgn = rhs2.signum()
        if (rhs2_sgn < 0) {
            val r90 = ExactFloatType.cast(S1ChordAngle.right().length2)
            val ca = if (lhs2_sgn < 0) -1 else exactCompareDistance(a, x0, r90)
            val cb = if (lhs2_sgn > 0) -1 else exactCompareDistance(b, x1, r90)
            if (ca <= 0 && cb <= 0) return Excluded.NEITHER
            assert(ca != 1 || cb != 1)
            return if (ca == 1) Excluded.FIRST else Excluded.SECOND
        }
        if (lhs2_sgn == 0) {
            // If the RHS of (2) is zero as well (i.e., d == 0) then both sites are
            // equidistant from every point on edge X.  This case requires symbolic
            // perturbations, but it should already have been handled in
            // GetVoronoiSiteExclusion() (see the call to CompareDistances).
// TODO            assert(rhs2_sgn > 0)
            return Excluded.NEITHER
        }
        // Next we square both sides of (2), yielding
        //
        //      cos^2(r) (sb^2 + sa^2 - 2 sa sb) > (aXb.n)^2
        //
        // which can be rearranged to give
        //
        // (3)  cos^2(r) (sb^2 + sa^2) - (aXb.n)^2 > 2 cos^2(r) sa sb .
        //
        // The RHS of (3) is always non-negative, but we still need to check the
        // sign of the LHS.
        val cos_r = ExactFloatType.cast(1.0) - ExactFloatType.cast(0.5) * r2
        val cos2_r = cos_r * cos_r
        val lhs3 = cos2_r * (sa2 + sb2) - rhs2 * rhs2
        if (lhs3.signum() < 0) return Excluded.NEITHER

        // Otherwise we square both sides of (3) to obtain:
        //
        // (4)  LHS(3)^2  >  4 cos^4(r) sa^2 sb^2
        val lhs4 = lhs3 * lhs3
        val rhs4 = ExactFloatType.cast(4.0) * cos2_r * cos2_r * sa2 * sb2
        val result = (lhs4 - rhs4).signum()
        if (result < 0) return Excluded.NEITHER
        if (result == 0) {
            // We have |rb - ra| = d and d > 0.  This implies that one coverage
            // interval contains the other, but not properly: the two intervals share
            // a common endpoint.  The distance from each site to that point is
            // exactly "r", therefore we need to use symbolic perturbations.  Recall
            // that site A is considered closer to an equidistant point if and only if
            // A > B.  Therefore if (rb > ra && A > B) or (ra > rb && B > A) then each
            // site is closer to at least one point and neither site is excluded.
            //
            // Ideally this logic would be in a separate SymbolicVoronoiSiteExclusion
            // method for better testing, but this is not convenient because it needs
            // lhs_sgn (which requires exact arithmetic to compute).
            if ((lhs2_sgn > 0) == (a > b)) return Excluded.NEITHER
        }
        // At this point we know that one of the two sites is excluded.  The sign of
        // the LHS of (2) (before the absolute value) determines which one.
        return if (lhs2_sgn > 0) Excluded.FIRST else Excluded.SECOND
    }

}


data class ValueAndError<T>(val value: T, val error: T)
