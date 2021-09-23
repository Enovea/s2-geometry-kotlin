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

import dilivia.math.DoubleType
import dilivia.math.ExactFloatType
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.toExactFloat
import dilivia.math.vectors.times
import dilivia.s2.S1Angle.Companion.tan
import dilivia.s2.S2Predicates.Excluded.FIRST
import dilivia.s2.S2Predicates.Excluded.NEITHER
import dilivia.s2.S2Predicates.Excluded.SECOND
import dilivia.s2.S2Predicates.Excluded.UNCERTAIN
import dilivia.s2.S2PredicatesUnitTest.PrecisionStats.Precision.DOUBLE
import dilivia.s2.S2PredicatesUnitTest.PrecisionStats.Precision.EXACT
import dilivia.s2.S2PredicatesUnitTest.PrecisionStats.Precision.LONG_DOUBLE
import dilivia.s2.S2PredicatesUnitTest.SignTest.testGreatCircle
import dilivia.s2.S2Random.oneIn
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomPoint
import dilivia.s2.edge.S2EdgeDistances
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Disabled
import org.junit.jupiter.api.Test
import tech.units.indriya.quantity.Quantities
import tech.units.indriya.unit.Units.METRE
import java.lang.Math.pow
import kotlin.math.nextTowards
import kotlin.math.pow


class S2PredicatesUnitTest {

    private val logger = KotlinLogging.logger { }
    private val kConsistencyIters: Int = 5000

    @Test
    fun signCollinearPoints() {
        // The following points happen to be *exactly collinear* along a line that it
        // approximate tangent to the surface of the unit sphere.  In fact, C is the
        // exact midpoint of the line segment AB.  All of these points are close
        // enough to unit length to satisfy .isUnitLength().
        val a = S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852)
        val b = S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738)
        val c = S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795)
        assertThat(c - a).isEqualTo(b - c)
        assertThat(S2Predicates.sign(a, b, c)).isNotEqualTo(0)
        assertThat(S2Predicates.sign(a, b, c)).isEqualTo(S2Predicates.sign(b, c, a))
        assertThat(S2Predicates.sign(a, b, c)).isEqualTo(-S2Predicates.sign(c, b, a))

        // The points "x1" and "x2" are exactly proportional, i.e. they both lie
        // on a common line through the origin.  Both points are considered to be
        // normalized, and in fact they both satisfy (x == x.normalize()).
        // Therefore the triangle (x1, x2, -x1) consists of three distinct points
        // that all lie on a common line through the origin.
        val x1 = S2Point(0.99999999999999989, 1.4901161193847655e-08, 0.0)
        val x2 = S2Point(1.0, 1.4901161193847656e-08, 0.0)
        assertThat(x1.normalized()).isEqualTo(x1)
        assertThat(x2.normalized()).isEqualTo(x2)
        assertThat(S2Predicates.sign(x1, x2, -x1)).isNotEqualTo(0)
        assertThat(S2Predicates.sign(x1, x2, -x1)).isEqualTo(S2Predicates.sign(x2, -x1, x1))
        assertThat(S2Predicates.sign(x1, x2, -x1)).isEqualTo(-S2Predicates.sign(-x1, x2, x1))

        // Here are two more points that are distinct, exactly proportional, and
        // that satisfy (x == x.normalize()).
        val x3 = S2Point(1, 1, 1).normalize()
        val x4 = 0.99999999999999989 * x3
        assertThat(x3.normalize()).isEqualTo(x3)
        assertThat(x4.normalize()).isEqualTo(x4)
        assertThat(x3).isNotEqualTo(x4)
        assertThat(S2Predicates.sign(x3, x4, -x3)).isNotEqualTo(0)

        // The following two points demonstrate that Normalize() is not idempotent,
        // i.e. y0.normalize() != y0.normalize().normalize().  Both points satisfy
        // S2::IsNormalized(), though, and the two points are exactly proportional.
        val y0 = S2Point(1, 1, 0)
        val y1 = y0.normalized()
        val y2 = y1.normalized()
        assertThat(y1).isNotEqualTo(y2)
        assertThat(y2.normalized()).isEqualTo(y2)
        assertThat(S2Predicates.sign(y1, y2, -y1)).isNotEqualTo(0)
        assertThat(S2Predicates.sign(y1, y2, -y1)).isEqualTo(S2Predicates.sign(y2, -y1, y1))
        assertThat(S2Predicates.sign(y1, y2, -y1)).isEqualTo(-S2Predicates.sign(-y1, y2, y1))
    }

    // This test repeatedly constructs some number of points that are on or nearly
    // on a given great circle.  Then it chooses one of these points as the
    // "origin" and sorts the other points in CCW order around it.  Of course,
    // since the origin is on the same great circle as the points being sorted,
    // nearly all of these tests are degenerate.  It then does various consistency
    // checks to verify that the points are indeed sorted in CCW order.
    //
    // It is easier to think about what this test is doing if you imagine that the
    // points are in general position rather than on a great circle.
    object SignTest {
        private val logger = KotlinLogging.logger { }

        // The following method is used to sort a collection of points in CCW order
        // around a given origin.  It returns true if A comes before B in the CCW
        // ordering (starting at an arbitrary fixed direction).
        class LessCCW(private val origin: S2Point, private val start: S2Point) : Comparator<S2Point> {

            override fun compare(a: S2Point, b: S2Point): Int {
                // OrderedCCW() acts like "<=", so we need to invert the comparison.
                return if (S2Predicates.orderedCCW(start, b, a, origin)) 1 else -1
            }

        }

        // Given a set of points with no duplicates, first remove "origin" from
        // "points" (if it exists) and then sort the remaining points in CCW order
        // around "origin" putting the result in "sorted".
        fun sortCCW(points: List<S2Point>, origin: S2Point): List<S2Point> {
            // Make a copy of the points with "origin" removed.
            val sorted = points.filter { it != origin }.toMutableList()

            // Sort the points CCW around the origin starting at (*sorted)[0].
            sorted.sortWith(LessCCW(origin, sorted[0]))
            return sorted
        }

        // Given a set of points sorted circularly CCW around "origin", and the
        // index "start" of a point A, count the number of CCW triangles OAB over
        // all sorted points B not equal to A.  Also check that the results of the
        // CCW tests are consistent with the hypothesis that the points are sorted.
        fun countCCW(sorted: List<S2Point>, origin: S2Point, start: Int): Int {
            logger.trace {
                """
                |sorted=$sorted
                |origin=$origin
                |start=$start
            """.trimMargin()
            }
            var numCcw = 0
            var lastSign = 1
            val n = sorted.size
            for (j in 1 until n) {
                val sign = S2Predicates.sign(origin, sorted[start], sorted[(start + j) % n])
                assertThat(sign).isNotEqualTo(0)
                if (sign > 0) ++numCcw

                // Since the points are sorted around the origin, we expect to see a
                // (possibly empty) sequence of CCW triangles followed by a (possibly
                // empty) sequence of CW triangles.
                assertThat(sign > 0 && lastSign < 0).isFalse()
                lastSign = sign
            }
            return numCcw
        }

        // Test exhaustively whether the points in "sorted" are sorted circularly
        // CCW around "origin".
        fun testCCW(sorted: List<S2Point>, origin: S2Point) {
            val n = sorted.size
            var totalNumCcw = 0
            var lastNumCcw = countCCW(sorted, origin, n - 1)
            for (start in 0 until n) {
                val numCcw = countCCW(sorted, origin, start)
                // Each iteration we increase the start index by 1, therefore the number
                // of CCW triangles should decrease by at most 1.
                assertThat(numCcw).isGreaterThanOrEqualTo(lastNumCcw - 1)
                totalNumCcw += numCcw
                lastNumCcw = numCcw
            }
            // We have tested all triangles of the form OAB.  Exactly half of these
            // should be CCW.
            assertThat(totalNumCcw).isEqualTo(n * (n - 1) / 2)
        }

        fun addNormalized(a: S2Point, points: MutableList<S2Point>) {
            val element = a.normalized()
            for (i in 0..2) if (element[i] == -0.0) a[i] = 0.0
            points.add(element)
        }

        // Add two points A1 and A2 that are slightly offset from A along the
        // tangent toward B, and such that A, A1, and A2 are exactly collinear
        // (i.e. even with infinite-precision arithmetic).
        fun addTangentPoints(a: S2Point, b: S2Point, points: MutableList<S2Point>) {
            val dir = S2PointUtil.robustCrossProd(a, b).crossProd(a).normalize()
            if (dir == S2Point(0, 0, 0)) return
            while (true) {
                val delta = 1e-15 * S2Random.randomDouble() * dir
                if ((a + delta) != a && (a + delta) - a == a - (a - delta) &&
                    S2PointUtil.isUnitLength(a + delta) && S2PointUtil.isUnitLength(a - delta)
                ) {
                    points.add(a + delta)
                    points.add(a - delta)
                    return
                }
            }
        }

        // Add zero or more (but usually one) point that is likely to trigger
        // Sign() degeneracies among the given points.
        fun addDegeneracy(points: MutableList<S2Point>) {
            val a = points[S2Random.randomInt(points.size)].clone()
            val b = points[S2Random.randomInt(points.size)].clone()
            val coord = S2Random.randomInt(3)
            when (S2Random.randomInt(8)) {
                // Add a random point (not uniformly distributed) along the great circle AB.
                0 -> {
                    val c = S2Random.randomDouble(-1.0, 1.0) * a + S2Random.randomDouble(-1.0, 1.0) * b
                    if (c[0] == -0.0) c[0] = 0.0
                    if (c[1] == -0.0) c[1] = 0.0
                    if (c[2] == -0.0) c[2] = 0.0
                    addNormalized(c, points)
                }
                // Perturb one coordinate by the minimum amount possible.
                1 -> {
                    var perturbed = false
                    var i = 0
                    while (!perturbed && i < 3) {
                        if (a[i] != 0.0 && a[i] != -0.0) {
                            a[i] = a[i].nextTowards(if (oneIn(2)) 2.0 else -2.0)
                            addNormalized(a, points)
                            perturbed = true
                        }
                        ++i
                    }
                }
                // Perturb one coordinate by up to 1e-15.
                2 -> {
                    a[coord] += 1e-15 * S2Random.randomDouble(-1.0, 1.0)
                    addNormalized(a, points)
                }
                // Scale a point just enough so that it is different while still being
                // considered normalized.
                3 -> {
                    a *= if (oneIn(2)) (1 + 2e-16) else (1 - 1e-16)
                    if (S2PointUtil.isUnitLength(a)) points.add(a)
                }
                // Add the intersection point of AB with X=0, Y=0, or Z=0.
                4 -> {
                    val dir = S2Point(0, 0, 0)
                    dir[coord] = if (oneIn(2)) 1.0 else -1.0
                    val norm = S2PointUtil.robustCrossProd(a, b).normalize()
                    if (norm.norm2() > 0) {
                        addNormalized(S2PointUtil.robustCrossProd(dir, norm), points)
                    }
                }
                // Add two closely spaced points along the tangent at A to the great
                // circle through AB.
                5 -> {
                    addTangentPoints(a, b, points)
                }
                // Add two closely spaced points along the tangent at A to the great
                // circle through A and the X-axis.
                6 -> {
                    addTangentPoints(a, S2Point(1, 0, 0), points)
                }
                // Add the negative of a point.
                7 -> {
                    points.add(-a)
                }
            }
        }

        // Sort the points around the given origin, and then do some consistency
        // checks to verify that they are actually sorted.
        fun sortAndTest(points: List<S2Point>, origin: S2Point) {
            val sorted = sortCCW(points, origin)

            logger.debug {
                """
                | Points: origin = $origin
                | --------------------------------------
                | ${
                    sorted.mapIndexed { i, p ->
                        "$p: ${
                            if (i > 0) S2Predicates.orderedCCW(
                                sorted[0],
                                p,
                                sorted[i - 1],
                                origin
                            ) else ""
                        }"
                    }.joinToString("\n")
                }
                | --------------------------------------
            """.trimMargin()
            }

            testCCW(sorted, origin)
        }

        // Construct approximately "n" points near the great circle through A and B,
        // then sort them and test whether they are sorted.
        fun testGreatCircle(a: S2Point, b: S2Point, n: Int) {
            logger.debug { "testGreatCircle a = $a ; b = $b : n = $n" }
            val points = mutableListOf(a.normalized(), b.normalized())
            while (points.size < n) {
                addDegeneracy(points)
            }
            // Remove any (0, 0, 0) points that were accidentically created, then sort
            // the points and remove duplicates.
            val pointsSet = points.filter { it != S2Point(0, 0, 0) }.toSortedSet()
            points.clear()
            points.addAll(pointsSet)
//            assertThat(points.size).isGreaterThanOrEqualTo(n / 2)

            sortAndTest(points, a)
            sortAndTest(points, b)
            for (origin in points) {
                sortAndTest(points, origin)
            }
        }
    }

    @Test
    @Disabled("to investigate")
    fun stressTest() {
        // The run time of this test is *cubic* in the parameter below.
        val kNumPointsPerCircle = 20

        // This test is randomized, so it is beneficial to run it several times.
        repeat(3) { iter ->

            logger.debug { "\niter = $iter\n" }

            // The most difficult great circles are the ones in the X-Y, Y-Z, and X-Z
            // planes, for two reasons.  First, when one or more coordinates are close
            // to zero then the perturbations can be much smaller, since floating
            // point numbers are spaced much more closely together near zero.  (This
            // tests the handling of things like underflow.)  The second reason is
            // that most of the cases of SymbolicallyPerturbedSign() can only be
            // reached when one or more input point coordinates are zero.
            testGreatCircle(S2Point(1, 0, 0), S2Point(0, 1, 0), kNumPointsPerCircle)
            testGreatCircle(S2Point(1, 0, 0), S2Point(0, 0, 1), kNumPointsPerCircle)
            testGreatCircle(S2Point(0, -1, 0), S2Point(0, 0, 1), kNumPointsPerCircle)

            // This tests a great circle where at least some points have X, Y, and Z
            // coordinates with exactly the same mantissa.  One useful property of
            // such points is that when they are scaled (e.g. multiplying by 1+eps),
            // all such points are exactly collinear with the origin.
            testGreatCircle(S2Point(1 shl 25, 1, -8), S2Point(-4, -(1 shl 20), 1), kNumPointsPerCircle)
        }
    }


    inner class StableSignTest {
        val logger = KotlinLogging.logger { }
        val kIters = 1000

        // Estimate the probability that S2::StableSign() will not be able to compute
        // the determinant sign of a triangle A, B, C consisting of three points
        // that are as collinear as possible and spaced the given distance apart.
        fun getFailureRate(km: Double): Double {
            var failure_count = 0
            val m = tan(S2Earth.toAngle(Quantities.getQuantity(km, METRE)))
            for (iter in 0 until kIters) {
                val (a, x, _) = S2Random.randomFrame()
                val b = (a - m * x).normalize()
                val c = (a + m * x).normalize()
                val sign = S2Predicates.stableSign(a, b, c)
                if (sign != 0) {
                    assertThat(sign).isEqualTo(S2Predicates.exactSign(a, b, c, true))
                } else {
                    ++failure_count
                }
            }
            val rate = failure_count.toDouble() / kIters
            logger.info { "StableSign failure rate for $km km = $rate" }
            return rate
        }
    }

    @Test
    fun stableSignTestFailureRate() {
        // Verify that StableSign() is able to handle most cases where the three
        // points are as collinear as possible.  (For reference, TriageSign() fails
        // virtually 100% of the time on this test.)
        //
        // Note that the failure rate *decreases* as the points get closer together,
        // and the decrease is approximately linear.  For example, the failure rate
        // is 0.4% for collinear points spaced 1km apart, but only 0.0004% for
        // collinear points spaced 1 meter apart.
        val stableSignTest = StableSignTest()
        assertThat(stableSignTest.getFailureRate(1.0) < 0.01).isTrue()
        //  1km spacing: <  1% (actual 0.4%)
        assertThat(stableSignTest.getFailureRate(10.0) < 0.1).isTrue()
// 10km spacing: < 10% (actual 4%)
    }

    // Given 3 points A, B, C that are exactly coplanar with the origin and where
// A < B < C in lexicographic order, verify that ABC is counterclockwise (if
// expected == 1) or clockwise (if expected == -1) using ExpensiveSign().
//
// This method is intended specifically for checking the cases where
// symbolic perturbations are needed to break ties.
    fun checkSymbolicSign(expected: Int, a: S2Point, b: S2Point, c: S2Point) {
        assert(a < b)
        assert(b < c)
        assert(0.0 == a.dotProd(b.crossProd(c)))

        // Use ASSERT rather than EXPECT to suppress spurious error messages.
        assertThat(S2Predicates.expensiveSign(a, b, c)).isEqualTo(expected)
        assertThat(S2Predicates.expensiveSign(b, c, a)).isEqualTo(expected)
        assertThat(S2Predicates.expensiveSign(c, a, b)).isEqualTo(expected)
        assertThat(S2Predicates.expensiveSign(c, b, a)).isEqualTo(-expected)
        assertThat(S2Predicates.expensiveSign(b, a, c)).isEqualTo(-expected)
        assertThat(S2Predicates.expensiveSign(a, c, b)).isEqualTo(-expected)
    }

    @Test
    fun signSymbolicPerturbationCodeCoverage() {
        // The purpose of this test is simply to get code coverage of
        // SymbolicallyPerturbedSign().  Let M_1, M_2, ... be the sequence of
        // submatrices whose determinant sign is tested by that function.  Then the
        // i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
        //
        //    det(M) = 0
        //    det(M_j) = 0 for j < i
        //    det(M_i) != 0
        //    A < B < C in lexicographic order.
        //
        // I checked that reversing the sign of any of the "return" statements in
        // SymbolicallyPerturbedSign() will cause this test to fail.

        // det(M_1) = b0*c1 - b1*c0
        checkSymbolicSign(1, S2Point(-3, -1, 0), S2Point(-2, 1, 0), S2Point(1, -2, 0))

        // det(M_2) = b2*c0 - b0*c2
        checkSymbolicSign(1, S2Point(-6, 3, 3), S2Point(-4, 2, -1), S2Point(-2, 1, 4))

        // det(M_3) = b1*c2 - b2*c1
        checkSymbolicSign(1, S2Point(0, -1, -1), S2Point(0, 1, -2), S2Point(0, 2, 1))
        // From this point onward, B or C must be zero, or B is proportional to C.

        // det(M_4) = c0*a1 - c1*a0
        checkSymbolicSign(1, S2Point(-1, 2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8))

        // det(M_5) = c0
        checkSymbolicSign(1, S2Point(-4, -2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8))

        // det(M_6) = -c1
        checkSymbolicSign(1, S2Point(0, -5, 7), S2Point(0, -4, 8), S2Point(0, -2, 4))

        // det(M_7) = c2*a0 - c0*a2
        checkSymbolicSign(1, S2Point(-5, -2, 7), S2Point(0, 0, -2), S2Point(0, 0, -1))

        // det(M_8) = c2
        checkSymbolicSign(1, S2Point(0, -2, 7), S2Point(0, 0, 1), S2Point(0, 0, 2))
        // From this point onward, C must be zero.

        // det(M_9) = a0*b1 - a1*b0
        checkSymbolicSign(1, S2Point(-3, 1, 7), S2Point(-1, -4, 1), S2Point(0, 0, 0))

        // det(M_10) = -b0
        checkSymbolicSign(1, S2Point(-6, -4, 7), S2Point(-3, -2, 1), S2Point(0, 0, 0))

        // det(M_11) = b1
        checkSymbolicSign(-1, S2Point(0, -4, 7), S2Point(0, -2, 1), S2Point(0, 0, 0))

        // det(M_12) = a0
        checkSymbolicSign(-1, S2Point(-1, -4, 5), S2Point(0, 0, -3), S2Point(0, 0, 0))

        // det(M_13) = 1
        checkSymbolicSign(1, S2Point(0, -4, 5), S2Point(0, 0, -5), S2Point(0, 0, 0))
    }


    // A helper class that keeps track of how often each precision was used and
// generates a string for logging purposes.
    class PrecisionStats {
        enum class Precision { DOUBLE, LONG_DOUBLE, EXACT, SYMBOLIC }

        val kPrecisionNames = arrayOf("double", "long double", "exact", "symbolic")
        val counts: MutableMap<Precision, Int> = mutableMapOf()

        fun tally(precision: Precision) {
            counts.putIfAbsent(precision, 0)
            counts[precision] = counts[precision]!! + 1
        }

        override fun toString(): String {
            return """Precision Stats
               |---
               |${kPrecisionNames[DOUBLE.ordinal]}: ${counts[DOUBLE]}
               |${kPrecisionNames[LONG_DOUBLE.ordinal]}: ${counts[LONG_DOUBLE]}
               |${kPrecisionNames[EXACT.ordinal]}: ${counts[EXACT]}
               |${kPrecisionNames[Precision.SYMBOLIC.ordinal]}: ${counts[Precision.SYMBOLIC]}
           """.trimMargin()
        }
    }


    // Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
    fun choosePoint(): S2Point {
        val x = randomPoint().clone()
        for (i in 0..2) {
            if (oneIn(3)) {
                x[i] = x[i].times(1e-50.pow(randomDouble()))
            }
        }
        return x.normalize()
    }

    // The following helper classes allow us to test the various distance
// calculation methods using a common test framework.
    interface CompareDistancesWrapper {
        fun triage(x: S2Point, a: S2Point, b: S2Point): Int
        // fun triage(x: R3VectorLongDouble, a: R3VectorLongDouble, b: R3VectorLongDouble): Int
    }

    class Sin2Distances : CompareDistancesWrapper {
        override fun triage(x: S2Point, a: S2Point, b: S2Point): Int {
            return S2Predicates.triageCompareSin2Distances(x, a, b)
        }

//    override fun triage(x: R3VectorLongDouble, a: R3VectorLongDouble, b: R3VectorLongDouble): Int {
//        return S2Predicates.triageCompareSin2Distances(x, a, b)
//    }
    }

    class CosDistances : CompareDistancesWrapper {
        override fun triage(x: S2Point, a: S2Point, b: S2Point): Int {
            return S2Predicates.triageCompareCosDistances(x, a, b)
        }

//    override fun triage(x: R3VectorLongDouble, a: R3VectorLongDouble, b: R3VectorLongDouble): Int {
//        return S2Predicates.triageCompareCosDistances(x, a, b)
//    }
    }

    // Compares distances greater than 90 degrees using sin^2(distance).
    class MinusSin2Distances : CompareDistancesWrapper {
        override fun triage(x: S2Point, a: S2Point, b: S2Point): Int {
            return -S2Predicates.triageCompareSin2Distances(-x, a, b)
        }

//    override fun triage(x: R3VectorLongDouble, a: R3VectorLongDouble, b: R3VectorLongDouble): Int {
//        return -S2Predicates.triageCompareSin2Distances(-x, a, b)
//    }
    }

    // Verifies that CompareDistances(x, a, b) == expected_sign, and furthermore
// checks that the minimum required precision is "expected_prec" when the
// distance calculation method defined by CompareDistancesWrapper is used.
    fun testCompareDistances(
        x: S2Point,
        a: S2Point,
        b: S2Point,
        expected_sign: Int,
        expected_prec: PrecisionStats.Precision,
        wrapper: CompareDistancesWrapper
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        x.normalize()
        a.normalize()
        b.normalize()

        val dbl_sign = wrapper.triage(x, a, b)
        //val ld_sign = wrapper.triage(x.toLongDouble(), a.toLongDouble(), b.toLongDouble())
        val exact_sign = S2Predicates.exactCompareDistances(x.toExactFloat(), a.toExactFloat(), b.toExactFloat())
        val actual_sign = if (exact_sign != 0) exact_sign else S2Predicates.symbolicCompareDistances(x, a, b)

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        assertThat(actual_sign).isEqualTo(expected_sign)
        if (exact_sign != 0) assertThat(actual_sign).isEqualTo(exact_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)

        val actual_prec = if (dbl_sign != 0) DOUBLE
        //else if (ld_sign != 0) LONG_DOUBLE
        else if (exact_sign != 0) EXACT
        else PrecisionStats.Precision.SYMBOLIC
        assertThat(actual_prec).isEqualTo(expected_prec)

        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.compareDistances(x, a, b)).isEqualTo(expected_sign)

        // Check that reversing the arguments negates the result.
        assertThat(S2Predicates.compareDistances(x, b, a)).isEqualTo(-expected_sign)
    }

    @Test
    fun compareDistancesCoverage() {
        // This test attempts to exercise all the code paths in all precisions.
        // Test TriageCompareSin2Distances.
        val sin2Distances = Sin2Distances()
        testCompareDistances(
            S2Point(1, 1, 1),
            S2Point(1.0, 1 - 1e-15, 1.0),
            S2Point(1.0, 1.0, 1 + 2e-15),
            -1,
            DOUBLE,
            sin2Distances
        )
        testCompareDistances(
            S2Point(1, 1, 0),
            S2Point(1.0, 1 - 1e-15, 1e-21),
            S2Point(1.0, 1 - 1e-15, 0.0),
            1,
            DOUBLE,
            sin2Distances
        )
        testCompareDistances(
            S2Point(2, 0, 0),
            S2Point(2, -1, 0),
            S2Point(2.0, 1.0, 1e-8),
            -1,
            EXACT,
            sin2Distances
        )
        testCompareDistances(S2Point(2, 0, 0), S2Point(2, -1, 0), S2Point(2.0, 1.0, 1e-100), -1, EXACT, sin2Distances)
        testCompareDistances(
            S2Point(1, 0, 0),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            1,
            PrecisionStats.Precision.SYMBOLIC,
            sin2Distances
        )
        testCompareDistances(
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            0,
            PrecisionStats.Precision.SYMBOLIC,
            sin2Distances
        )

        // Test TriageCompareCosDistances.
        val cosDistances = CosDistances()
        testCompareDistances(S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1.0, 1.0, 3e-15), 1, DOUBLE, cosDistances)
        testCompareDistances(
            S2Point(1, 0, 0),
            S2Point(1.0, 1e-30, 0.0),
            S2Point(-1.0, 1e-40, 0.0),
            -1,
            DOUBLE,
            cosDistances
        )
        testCompareDistances(
            S2Point(1, 1, 1),
            S2Point(1, -1, 0),
            S2Point(-1.0, 1.0, 3e-18),
            1,
            EXACT,
            cosDistances
        )
        testCompareDistances(S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1.0, 1.0, 1e-100), 1, EXACT, cosDistances)
        testCompareDistances(
            S2Point(1, 1, 1),
            S2Point(1, -1, 0),
            S2Point(-1, 1, 0),
            -1,
            PrecisionStats.Precision.SYMBOLIC,
            cosDistances
        )
        testCompareDistances(
            S2Point(1, 1, 1),
            S2Point(1, -1, 0),
            S2Point(1, -1, 0),
            0,
            PrecisionStats.Precision.SYMBOLIC,
            cosDistances
        )

        // Test TriageCompareSin2Distances using distances greater than 90 degrees.
        val minSin2Distances = MinusSin2Distances()
        testCompareDistances(
            S2Point(1, 1, 0),
            S2Point(-1.0, -1 + 1e-15, 0.0),
            S2Point(-1, -1, 0),
            -1,
            DOUBLE,
            minSin2Distances
        )
        testCompareDistances(
            S2Point(-1, -1, 0),
            S2Point(1.0, 1 - 1e-15, 0.0),
            S2Point(1.0, 1 - 1e-15, 1e-21),
            1,
            DOUBLE,
            minSin2Distances
        )
        testCompareDistances(
            S2Point(-1, -1, 0),
            S2Point(2, 1, 0),
            S2Point(2.0, 1.0, 1e-8),
            1,
            EXACT,
            minSin2Distances
        )
        testCompareDistances(S2Point(-1, -1, 0), S2Point(2, 1, 0), S2Point(2.0, 1.0, 1e-30), 1, EXACT, minSin2Distances)
        testCompareDistances(
            S2Point(-1, -1, 0),
            S2Point(2, 1, 0),
            S2Point(1, 2, 0),
            -1,
            PrecisionStats.Precision.SYMBOLIC,
            minSin2Distances
        )
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testCompareDistancesConsistency(
        x: S2Point,
        a: S2Point,
        b: S2Point,
        wrapper: CompareDistancesWrapper
    ): PrecisionStats.Precision {
        val dbl_sign = wrapper.triage(x, a, b)
        //val ld_sign = wrapper.triage(x.toLongDouble(), a.toLongDouble(), b.toLongDouble())
        val exact_sign = S2Predicates.exactCompareDistances(x.toExactFloat(), a.toExactFloat(), b.toExactFloat())
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        if (exact_sign != 0) {
            val compareDistances = S2Predicates.compareDistances(x, a, b)
            assertThat(compareDistances)
                .withFailMessage(
                    "$exact_sign != S2Predicates.compareDistances(x = $x, a = $a, b = $b) (=${S2Predicates.compareDistances(x, a, b)})"
                ).isEqualTo(exact_sign)
            return if (dbl_sign == 0) EXACT else DOUBLE
            //return if (ld_sign == 0) PrecisionStats.Precision.EXACT else if (dbl_sign == 0) PrecisionStats.Precision.LONG_DOUBLE else PrecisionStats.Precision.DOUBLE
        } else {
            // Unlike the other methods, SymbolicCompareDistances has the
            // precondition that the exact sign must be zero.
            val symbolic_sign = S2Predicates.symbolicCompareDistances(x, a, b)
            assertThat(S2Predicates.compareDistances(x, a, b)).isEqualTo(symbolic_sign)
            return PrecisionStats.Precision.SYMBOLIC
        }
    }

    @Test
    fun compareDistancesConsistency() {
        // This test chooses random point pairs that are nearly equidistant from a
        // target point, and then checks that the answer given by a method at one
        // level of precision is consistent with the answer given at the next higher
        // level of precision.
        //
        // The way the .cc file is structured, we can only do comparisons using a
        // specific precision if we also choose the specific distance calculation
        // method.  The code below checks that the Cos, Sin2, and MinusSin2 methods
        // are consistent across their entire valid range of inputs, and also
        // simulates the logic in CompareDistance that chooses which method to use
        // in order to gather statistics about how often each precision is needed.
        // (These statistics are only useful for coverage purposes, not benchmarks,
        // since the input points are chosen to be pathological worst cases.)
        val cosDistances = CosDistances()
        val sin2Distances = Sin2Distances()
        val minusSin2Distances = MinusSin2Distances()
        testCompareDistancesConsistency(S2Point(1, 0, 0), S2Point(0, -1, 0), S2Point(0, 1, 0), cosDistances)
        val sin2_stats = PrecisionStats()
        val cos_stats = PrecisionStats()
        val minus_sin2_stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val x = choosePoint()
            val dir = choosePoint()
            var r = S1Angle.radians(M_PI_2 * pow(1e-30, S2Random.randomDouble()))
            if (S2Random.oneIn(2)) r = S1Angle.radians(M_PI_2) - r
            if (S2Random.oneIn(2)) r = S1Angle.radians(M_PI_2) + r
            val a = S2EdgeDistances.interpolateAtDistance(r, x, dir)
            val b = S2EdgeDistances.interpolateAtDistance(r, x, -dir)
            var prec = testCompareDistancesConsistency(x, a, b, cosDistances)
            if (r.degrees() in 45.0..135.0) cos_stats.tally(prec)
            // The Sin2 method is only valid if both distances are less than 90
            // degrees, and similarly for the MinusSin2 method.  (In the actual
            // implementation these methods are only used if both distances are less
            // than 45 degrees or greater than 135 degrees respectively.)
            if (r.radians < M_PI_2 - 1e-14) {
                prec = testCompareDistancesConsistency(x, a, b, sin2Distances)
                if (r.degrees() < 45) {
                    // Don't skew the statistics by recording degenerate inputs.
                    if (a == b) {
                        assertThat(prec).isEqualTo(PrecisionStats.Precision.SYMBOLIC)
                    } else {
                        sin2_stats.tally(prec)
                    }
                }
            } else if (r.radians > M_PI_2 + 1e-14) {
                prec = testCompareDistancesConsistency(x, a, b, minusSin2Distances)
                if (r.degrees() > 135) minus_sin2_stats.tally(prec)
            }
        }

        logger.error { "\nsin2:\n  $sin2_stats\ncos:\n   $cos_stats\n-sin2:\n $minus_sin2_stats" }
    }

    interface CompareDistanceWrapper {
        fun triage(x: S2Point, y: S2Point, r: S1ChordAngle): Int
        // fun triage(x: R3VectorLongDouble, y: R3VectorLongDouble, r: S1ChordAngle): Int
    }

    // Helper classes for testing the various distance calculation methods.
    class Sin2Distance : CompareDistanceWrapper {
        override fun triage(x: S2Point, y: S2Point, r: S1ChordAngle): Int =
            S2Predicates.triageCompareSin2Distance(x, y, r.length2)
        //override fun triage(x: R3VectorLongDouble, y: R3VectorLongDouble, r: S1ChordAngle): Int = S2Predicates.triageCompareSin2Distance(x, y, LongDoubleType.fromDouble(r.length2))
    }

    class CosDistance : CompareDistanceWrapper {
        override fun triage(x: S2Point, y: S2Point, r: S1ChordAngle): Int =
            S2Predicates.triageCompareCosDistance(x, y, r.length2)
        //override fun triage(x: R3VectorLongDouble, y: R3VectorLongDouble, r: S1ChordAngle): Int = S2Predicates.triageCompareCosDistance(x, y, LongDoubleType.fromDouble(r.length2))

    }

    // Verifies that CompareDistance(x, y, r) == expected_sign, and furthermore
// checks that the minimum required precision is "expected_prec" when the
// distance calculation method defined by CompareDistanceWrapper is used.
    fun testCompareDistance(
        x: S2Point,
        y: S2Point,
        r: S1ChordAngle,
        expected_sign: Int,
        expected_prec: PrecisionStats.Precision,
        wrapper: CompareDistanceWrapper
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!S2PointUtil.isUnitLength(x)) x.normalize()
        if (!S2PointUtil.isUnitLength(y)) y.normalize()

        val dbl_sign = wrapper.triage(x, y, r)
        //val ld_sign = wrapper.triage(x.toLongDouble(), y.toLongDouble(), r)
        val exact_sign =
            S2Predicates.exactCompareDistance(x.toExactFloat(), y.toExactFloat(), ExactFloatType.cast(r.length2))

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        assertThat(exact_sign).isEqualTo(expected_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)

        //val actual_prec = if (dbl_sign != 0) DOUBLE else if (ld_sign != 0) LONG_DOUBLE else EXACT
        val actual_prec = if (dbl_sign != 0) PrecisionStats.Precision.DOUBLE else PrecisionStats.Precision.EXACT
        assertThat(actual_prec).isEqualTo(expected_prec)

        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.compareDistance(x, y, r)).isEqualTo(expected_sign)

        // Mathematically, if d(X, Y) < r then d(-X, Y) > (Pi - r).  Unfortunately
        // there can be rounding errors when computing the supplementary distance,
        // so to ensure the two distances are exactly supplementary we need to do
        // the following.
        val r_supp = S1ChordAngle.straight() - r
        val r = S1ChordAngle.straight() - r_supp
        assertThat(S2Predicates.compareDistance(-x, y, r_supp)).isEqualTo(-S2Predicates.compareDistance(x, y, r))
    }

    @Test
    fun compareDistanceCoverage() {
        val sin2Distance = Sin2Distance()
        val cosDistance = CosDistance()
        // Test TriageCompareSin2Distance.
        testCompareDistance(
            S2Point(1, 1, 1),
            S2Point(1.0, 1 - 1e-15, 1.0),
            S1ChordAngle.radians(1e-15),
            -1,
            DOUBLE,
            sin2Distance
        )
        testCompareDistance(
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(M_PI_4),
            -1,
            EXACT,
            sin2Distance
        )
        testCompareDistance(
            S2Point(1.0, 1e-40, 0.0),
            S2Point(1.0 + DoubleType.epsilon, 1e-40, 0.0),
            S1ChordAngle.radians(0.9 * DoubleType.epsilon * 1e-40),
            1,
            EXACT,
            sin2Distance
        )
        testCompareDistance(
            S2Point(1.0, 1e-40, 0.0),
            S2Point(1.0 + DoubleType.epsilon, 1e-40, 0.0),
            S1ChordAngle.radians(1.1 * DoubleType.epsilon * 1e-40),
            -1,
            EXACT,
            sin2Distance
        )
        testCompareDistance(
            S2Point(1, 0, 0),
            S2Point(1 + DoubleType.epsilon, 0.0, 0.0),
            S1ChordAngle.zero(),
            0,
            EXACT,
            sin2Distance
        )

        // Test TriageCompareCosDistance.
        testCompareDistance(
            S2Point(1, 0, 0),
            S2Point(1.0, 1e-8, 0.0),
            S1ChordAngle.radians(1e-7),
            -1,
            DOUBLE,
            cosDistance
        )
        testCompareDistance(
            S2Point(1, 0, 0),
            S2Point(-1.0, 1e-8, 0.0),
            S1ChordAngle.radians(M_PI - 1e-7),
            1,
            DOUBLE,
            cosDistance
        )
        testCompareDistance(
            S2Point(1, 1, 0),
            S2Point(1.0, -1 - 2 * DoubleType.epsilon, 0.0),
            S1ChordAngle.right(),
            1,
            DOUBLE,
            cosDistance
        )
        testCompareDistance(
            S2Point(1, 1, 0),
            S2Point(1.0, -1 - DoubleType.epsilon, 0.0),
            S1ChordAngle.right(),
            1,
            DOUBLE,
            cosDistance
        )
        testCompareDistance(S2Point(1, 1, 0), S2Point(1.0, -1.0, 1e-30), S1ChordAngle.right(), 0, EXACT, cosDistance)
        // The angle between these two points is exactly 60 degrees.
        testCompareDistance(S2Point(1, 1, 0), S2Point(0, 1, 1), S1ChordAngle.fromLength2(1), 0, EXACT, cosDistance)
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testCompareDistanceConsistency(
        x: S2Point,
        y: S2Point,
        r: S1ChordAngle,
        wrapper: CompareDistanceWrapper
    ): PrecisionStats.Precision {
        val dbl_sign = wrapper.triage(x, y, r)
        //val ld_sign = wrapper.triage(x.toLongDouble(), y.toLongDouble(), r)
        val exact_sign =
            S2Predicates.exactCompareDistance(x.toExactFloat(), y.toExactFloat(), ExactFloatType.cast(r.length2))
        assertThat(S2Predicates.compareDistance(x, y, r)).isEqualTo(exact_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        return if (dbl_sign == 0) EXACT else DOUBLE
        //return if(ld_sign == 0) PrecisionStats.Precision.EXACT else if(dbl_sign == 0) PrecisionStats.Precision.LONG_DOUBLE else PrecisionStats.Precision.DOUBLE
    }

    @Test
    fun compareDistanceConsistency() {
        // This test chooses random inputs such that the distance between points X
        // and Y is very close to the threshold distance "r".  It then checks that
        // the answer given by a method at one level of precision is consistent with
        // the answer given at the next higher level of precision.  See also the
        // comments in the CompareDistances consistency test.
        val sin2_stats = PrecisionStats()
        val cos_stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val x = choosePoint()
            val dir = choosePoint()
            var r = S1Angle.radians(M_PI_2 * pow(1e-30, randomDouble()))
            if (oneIn(2)) r = S1Angle.radians(M_PI_2) - r
            if (oneIn(5)) r = S1Angle.radians(M_PI_2) + r
            val y = S2EdgeDistances.interpolateAtDistance(r, x, dir)
            var prec = testCompareDistanceConsistency(x, y, S1ChordAngle(r), CosDistance())
            if (r.degrees() >= 45) cos_stats.tally(prec)
            if (r.radians < M_PI_2 - 1e-14) {
                prec = testCompareDistanceConsistency(x, y, S1ChordAngle(r), Sin2Distance())
                if (r.degrees() < 45) sin2_stats.tally(prec)
            }
        }

        logger.error { "\nsin2:\n  $sin2_stats\ncos:\n   $cos_stats" }
    }

    // Verifies that CompareEdgeDistance(x, a0, a1, r) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
    fun testCompareEdgeDistance(
        x: S2Point,
        a0: S2Point,
        a1: S2Point,
        r: S1ChordAngle,
        expected_sign: Int,
        expected_prec: PrecisionStats.Precision
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!S2PointUtil.isUnitLength(x)) x.normalize()
        if (!S2PointUtil.isUnitLength(a0)) a0.normalize()
        if (!S2PointUtil.isUnitLength(a1)) a1.normalize()

        val dbl_sign = S2Predicates.triageCompareEdgeDistance(x, a0, a1, r.length2)
        //val ld_sign = S2Predicates.triageCompareEdgeDistance(x.toLongDouble(), a0.toLongDouble(), a1.toLongDouble(), r.length2.toLD())
        val exact_sign = S2Predicates.exactCompareEdgeDistance(x, a0, a1, r)

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        assertThat(exact_sign).isEqualTo(expected_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)

        //val actual_prec = if (dbl_sign != 0) DOUBLE else if (ld_sign != 0) LONG_DOUBLE else EXACT
        val actual_prec = if (dbl_sign != 0) DOUBLE else EXACT
        assertThat(actual_prec).isEqualTo(expected_prec)

        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.compareEdgeDistance(x, a0, a1, r)).isEqualTo(expected_sign)
    }

    @Test
    fun compareEdgeDistanceCoverage() {
        // Test TriageCompareLineSin2Distance.
        testCompareEdgeDistance(
            S2Point(1.0, 1e-10, 1e-15),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.radians(1e-15 + DoubleType.epsilon),
            -1,
            DOUBLE
        )
        testCompareEdgeDistance(
            S2Point(1.0, 1.0, 1e-15),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.radians(1e-15 + DoubleType.epsilon),
            -1,
            EXACT
        )
        testCompareEdgeDistance(
            S2Point(1.0, 1.0, 1e-40),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.radians(1e-40),
            -1,
            EXACT
        )
        testCompareEdgeDistance(S2Point(1, 1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), S1ChordAngle.zero(), 0, EXACT)

        // Test TriageCompareLineCos2Distance.
        testCompareEdgeDistance(
            S2Point(1e-15, 0.0, 1.0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.radians(M_PI_2 - 1e-15 - 3 * DoubleType.epsilon),
            1,
            DOUBLE
        )
        testCompareEdgeDistance(
            S2Point(1e-15, 0.0, 1.0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.radians(M_PI_2 - 1e-15 - DoubleType.epsilon),
            1,
            EXACT
        )
        testCompareEdgeDistance(
            S2Point(1e-40, 0.0, 1.0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S1ChordAngle.right(),
            -1,
            EXACT
        )
        testCompareEdgeDistance(S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), S1ChordAngle.right(), 0, EXACT)

        // Test cases where the closest point is an edge endpoint.
        testCompareEdgeDistance(
            S2Point(1e-15, -1.0, 0.0),
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.right(),
            -1,
            DOUBLE
        )
        testCompareEdgeDistance(
            S2Point(1e-18, -1.0, 0.0),
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.right(),
            -1,
            EXACT
        )
        testCompareEdgeDistance(
            S2Point(1e-100, -1.0, 0.0),
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.right(),
            -1,
            EXACT
        )
        testCompareEdgeDistance(S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0), S1ChordAngle.right(), 0, EXACT)
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testCompareEdgeDistanceConsistency(
        x: S2Point,
        a0: S2Point,
        a1: S2Point,
        r: S1ChordAngle
    ): PrecisionStats.Precision {
        val dbl_sign = S2Predicates.triageCompareEdgeDistance(x, a0, a1, r.length2)
        //val ld_sign = S2Predicates.triageCompareEdgeDistance(x.toLongDouble(), a0.toLongDouble(), a1.toLongDouble(), r.length2.toLD())
        val exact_sign = S2Predicates.exactCompareEdgeDistance(x, a0, a1, r)
        assertThat(S2Predicates.compareEdgeDistance(x, a0, a1, r)).isEqualTo(exact_sign)

        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)
        return if (exact_sign == 0) EXACT else DOUBLE

        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        //return if(ld_sign == 0) PrecisionStats.Precision.EXACT else if(dbl_sign == 0) PrecisionStats.Precision.LONG_DOUBLE else PrecisionStats.Precision.DOUBLE
    }

    @Test
    fun compareEdgeDistanceConsistency() {
        // This test chooses random inputs such that the distance between "x" and
        // the line (a0, a1) is very close to the threshold distance "r".  It then
        // checks that the answer given by a method at one level of precision is
        // consistent with the answer given at the next higher level of precision.
        // See also the comments in the CompareDistances consistency test.
        val stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val a0 = choosePoint()
            val len = S1Angle.radians(M_PI * pow(1e-20, randomDouble()))
            var a1 = S2EdgeDistances.interpolateAtDistance(len, a0, choosePoint())
            if (oneIn(2)) a1 = -a1
            if (a0 == -a1) continue;  // Not allowed by API.
            val n = S2PointUtil.robustCrossProd(a0, a1).normalize()
            val f = pow(1e-20, randomDouble())
            val a = ((1 - f) * a0 + f * a1).normalize()
            var r = S1Angle.radians(M_PI_2 * pow(1e-20, randomDouble()))
            if (oneIn(2)) r = S1Angle.radians(M_PI_2) - r
            var x = S2EdgeDistances.interpolateAtDistance(r, a, n)
            if (oneIn(5)) {
                // Replace "x" with a random point that is closest to an edge endpoint.
                do {
                    x = choosePoint()
                } while (S2Predicates.compareEdgeDirections(a0, x, a0, a1) > 0 && S2Predicates.compareEdgeDirections(
                        x,
                        a1,
                        a0,
                        a1
                    ) > 0
                )
                r = minOf(S1Angle(x, a0), S1Angle(x, a1))
            }
            val prec = testCompareEdgeDistanceConsistency(x, a0, a1, S1ChordAngle(r))
            stats.tally(prec)
        }
        logger.info { stats.toString() }
    }

    // Verifies that S2Predicates.compareEdgeDirections(a0, a1, b0, b1) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
    fun testCompareEdgeDirections(
        a0: S2Point,
        a1: S2Point,
        b0: S2Point,
        b1: S2Point,
        expected_sign: Int,
        expected_prec: PrecisionStats.Precision
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!S2PointUtil.isUnitLength(a0)) a0.normalize()
        if (!S2PointUtil.isUnitLength(a1)) a1.normalize()
        if (!S2PointUtil.isUnitLength(b0)) b0.normalize()
        if (!S2PointUtil.isUnitLength(b1)) b1.normalize()

        val dbl_sign = S2Predicates.triageCompareEdgeDirections(a0, a1, b0, b1)
//        val ld_sign = S2Predicates.triageCompareEdgeDirections(
//            a0.toLongDouble(),
//            a1.toLongDouble(),
//            b0.toLongDouble(),
//            b1.toLongDouble()
//        )
        val exact_sign = S2Predicates.exactCompareEdgeDirections(
            a0.toExactFloat(),
            a1.toExactFloat(),
            b0.toExactFloat(),
            b1.toExactFloat()
        )

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        assertThat(exact_sign).isEqualTo(expected_sign)
        //  if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        // if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)

        //val actual_prec = if (dbl_sign != 0) DOUBLE else if (ld_sign != 0) LONG_DOUBLE else EXACT
        val actual_prec = if (dbl_sign != 0) DOUBLE else EXACT
        assertThat(actual_prec).isEqualTo(expected_prec)

        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, b0, b1)).isEqualTo(expected_sign)

        // Check various identities involving swapping or negating arguments.
        assertThat(S2Predicates.compareEdgeDirections(b0, b1, a0, a1)).isEqualTo(expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(-a0, -a1, b0, b1)).isEqualTo(expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, -b0, -b1)).isEqualTo(expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a1, a0, b0, b1)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, b1, b0)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(-a0, a1, b0, b1)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a0, -a1, b0, b1)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, -b0, b1)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, b0, -b1)).isEqualTo(-expected_sign)
    }

    @Test
    fun compareEdgeDirectionsCoverage() {
        testCompareEdgeDirections(S2Point(1, 0, 0), S2Point(1, 1, 0), S2Point(1, -1, 0), S2Point(1, 0, 0), 1, DOUBLE)
        testCompareEdgeDirections(
            S2Point(1.0, 0.0, 1.5e-15),
            S2Point(1, 1, 0),
            S2Point(0, -1, 0),
            S2Point(0, 0, 1),
            1,
            DOUBLE
        )
        testCompareEdgeDirections(
            S2Point(1.0, 0.0, 1e-18),
            S2Point(1, 1, 0),
            S2Point(0, -1, 0),
            S2Point(0, 0, 1),
            1,
            EXACT
        )
        testCompareEdgeDirections(
            S2Point(1.0, 0.0, 1e-50),
            S2Point(1, 1, 0),
            S2Point(0, -1, 0),
            S2Point(0, 0, 1),
            1,
            EXACT
        )
        testCompareEdgeDirections(S2Point(1, 0, 0), S2Point(1, 1, 0), S2Point(0, -1, 0), S2Point(0, 0, 1), 0, EXACT)
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testCompareEdgeDirectionsConsistency(
        a0: S2Point,
        a1: S2Point,
        b0: S2Point,
        b1: S2Point
    ): PrecisionStats.Precision {
        val dbl_sign = S2Predicates.triageCompareEdgeDirections(a0, a1, b0, b1)
//        val ld_sign = S2Predicates.triageCompareEdgeDirections(
//            a0.toLongDouble(),
//            a1.toLongDouble(),
//            b0.toLongDouble(),
//            b1.toLongDouble()
//        )
        val exact_sign = S2Predicates.exactCompareEdgeDirections(
            a0.toExactFloat(),
            a1.toExactFloat(),
            b0.toExactFloat(),
            b1.toExactFloat()
        )
        assertThat(S2Predicates.compareEdgeDirections(a0, a1, b0, b1)).isEqualTo(exact_sign)
//        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
//        if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
//        return if (ld_sign == 0) EXACT else if (dbl_sign == 0) LONG_DOUBLE else DOUBLE
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)
        //  if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        return if (dbl_sign == 0) EXACT else DOUBLE
    }

    @Test
    fun compareEdgeDirectionsConsistency() {
        // This test chooses random pairs of edges that are nearly perpendicular,
        // then checks that the answer given by a method at one level of precision
        // is consistent with the answer given at the next higher level of
        // precision.  See also the comments in the CompareDistances test.
        val stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val a0 = choosePoint()
            val a_len = S1Angle.radians(M_PI * 1e-20.pow(randomDouble()))
            val a1 = S2EdgeDistances.interpolateAtDistance(a_len, a0, choosePoint())
            val a_norm = S2PointUtil.robustCrossProd(a0, a1).normalize()
            val b0 = choosePoint()
            val b_len = S1Angle.radians(M_PI * 1e-20.pow(randomDouble()))
            val b1 = S2EdgeDistances.interpolateAtDistance(b_len, b0, a_norm)
            if (a0 == -a1 || b0 == -b1) continue;  // Not allowed by API.
            val prec = testCompareEdgeDirectionsConsistency(a0, a1, b0, b1)
            // Don't skew the statistics by recording degenerate inputs.
            if (a0 == a1 || b0 == b1) {
                assertThat(prec).isEqualTo(EXACT)
            } else {
                stats.tally(prec)
            }
        }
        logger.info { stats.toString() }
    }

    // Verifies that S2Predicates.edgeCircumcenterSign(x0, x1, a, b, c) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
    fun testEdgeCircumcenterSign(
        x0: S2Point,
        x1: S2Point,
        a: S2Point,
        b: S2Point,
        c: S2Point,
        expected_sign: Int,
        expected_prec: PrecisionStats.Precision
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!S2PointUtil.isUnitLength(x0)) x0.normalize()
        if (!S2PointUtil.isUnitLength(x1)) x1.normalize()
        if (!S2PointUtil.isUnitLength(a)) a.normalize()
        if (!S2PointUtil.isUnitLength(b)) b.normalize()
        if (!S2PointUtil.isUnitLength(c)) c.normalize()

        val abc_sign = S2Predicates.sign(a, b, c)
        val dbl_sign = S2Predicates.triageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign)
//        val ld_sign = S2Predicates.triageEdgeCircumcenterSign(
//            x0.toLongDouble(),
//            x1.toLongDouble(),
//            a.toLongDouble(),
//            b.toLongDouble(),
//            c.toLongDouble(),
//            abc_sign
//        )
        val exact_sign = S2Predicates.exactEdgeCircumcenterSign(
            x0.toExactFloat(),
            x1.toExactFloat(),
            a.toExactFloat(),
            b.toExactFloat(),
            c.toExactFloat(),
            abc_sign
        )
        val actual_sign =
            if (exact_sign != 0) exact_sign else S2Predicates.symbolicEdgeCircumcenterSign(x0, x1, a, b, c)

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        assertThat(actual_sign).isEqualTo(expected_sign)
        if (exact_sign != 0) assertThat(actual_sign).isEqualTo(exact_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        //if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)

        //val actual_prec = if (dbl_sign != 0) DOUBLE else if (ld_sign != 0) LONG_DOUBLE else if (exact_sign != 0) EXACT else PrecisionStats.Precision.SYMBOLIC
        val actual_prec =
            if (dbl_sign != 0) DOUBLE else if (exact_sign != 0) EXACT else PrecisionStats.Precision.SYMBOLIC
        assertThat(actual_prec).isEqualTo(expected_prec)

        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, a, b, c)).isEqualTo(expected_sign)

        // Check various identities involving swapping or negating arguments.
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, a, c, b)).isEqualTo(expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, b, a, c)).isEqualTo(expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, b, c, a)).isEqualTo(expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, c, a, b)).isEqualTo(expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, c, b, a)).isEqualTo(expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(x1, x0, a, b, c)).isEqualTo(-expected_sign)
        assertThat(S2Predicates.edgeCircumcenterSign(-x0, -x1, a, b, c)).isEqualTo(expected_sign)
        if (actual_sign == exact_sign) {
            // Negating the input points may not preserve the result when symbolic
            // perturbations are used, since -X is not an exact multiple of X.
            assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, -a, -b, -c)).isEqualTo(-expected_sign)
        }
    }

    @Test
    fun edgeCircumcenterSignCoverage() {
        testEdgeCircumcenterSign(
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S2Point(0, 0, 1),
            S2Point(1, 0, 1),
            S2Point(0, 1, 1),
            1,
            DOUBLE
        )
        testEdgeCircumcenterSign(
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S2Point(0, 0, -1),
            S2Point(1, 0, -1),
            S2Point(0, 1, -1),
            -1,
            DOUBLE
        )
        testEdgeCircumcenterSign(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1.0, -1e-5, 1.0),
            S2Point(1.0, 1e-5, -1.0),
            S2Point(1.0, 1 - 1e-5, 1e-5),
            -1,
            DOUBLE
        )
        testEdgeCircumcenterSign(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1.0, -1e-5, 1.0),
            S2Point(1.0, 1e-5, -1.0),
            S2Point(1.0, 1 - 1e-9, 1e-5),
            -1,
            EXACT
        )
        testEdgeCircumcenterSign(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1.0, -1e-5, 1.0),
            S2Point(1.0, 1e-5, -1.0),
            S2Point(1.0, 1 - 1e-15, 1e-5),
            -1,
            EXACT
        )
        testEdgeCircumcenterSign(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1.0, -1e-5, 1.0),
            S2Point(1.0, 1e-5, -1.0),
            S2Point(1.0, 1.0, 1e-5),
            1,
            PrecisionStats.Precision.SYMBOLIC
        )

        // This test falls back to the second symbolic perturbation:
        testEdgeCircumcenterSign(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(0, -1, 0),
            S2Point(0, 0, -1),
            S2Point(0, 0, 1),
            -1,
            PrecisionStats.Precision.SYMBOLIC
        )

        // This test falls back to the third symbolic perturbation:
        testEdgeCircumcenterSign(
            S2Point(0, -1, 1),
            S2Point(0, 1, 1),
            S2Point(0, 1, 0),
            S2Point(0, -1, 0),
            S2Point(1, 0, 0),
            -1,
            PrecisionStats.Precision.SYMBOLIC
        )
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testEdgeCircumcenterSignConsistency(
        x0: S2Point,
        x1: S2Point,
        a: S2Point,
        b: S2Point,
        c: S2Point
    ): PrecisionStats.Precision {
        val abc_sign = S2Predicates.sign(a, b, c)
        val dbl_sign = S2Predicates.triageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign)
//        val ld_sign = S2Predicates.triageEdgeCircumcenterSign(
//            x0.toLongDouble(),
//            x1.toLongDouble(),
//            a.toLongDouble(),
//            b.toLongDouble(),
//            c.toLongDouble(),
//            abc_sign
//        )
        val exact_sign = S2Predicates.exactEdgeCircumcenterSign(
            x0.toExactFloat(),
            x1.toExactFloat(),
            a.toExactFloat(),
            b.toExactFloat(),
            c.toExactFloat(),
            abc_sign
        )
        // if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(ld_sign)
        if (dbl_sign != 0) assertThat(dbl_sign).isEqualTo(exact_sign)
        //if (ld_sign != 0) assertThat(ld_sign).isEqualTo(exact_sign)
        if (exact_sign != 0) {
            assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, a, b, c)).isEqualTo(exact_sign)
            //  return if (ld_sign == 0) EXACT else if (dbl_sign == 0) LONG_DOUBLE else DOUBLE
            return if (dbl_sign == 0) EXACT else DOUBLE
        } else {
            // Unlike the other methods, SymbolicEdgeCircumcenterSign has the
            // precondition that the exact sign must be zero.
            val symbolic_sign = S2Predicates.symbolicEdgeCircumcenterSign(x0, x1, a, b, c)
            assertThat(S2Predicates.edgeCircumcenterSign(x0, x1, a, b, c)).isEqualTo(symbolic_sign)
            return PrecisionStats.Precision.SYMBOLIC
        }
    }

    @Test
    fun edgeCircumcenterSignConsistency() {
        // This test chooses random a random edge X, then chooses a random point Z
        // on the great circle through X, and finally choose three points A, B, C
        // that are nearly equidistant from X.  It then checks that the answer given
        // by a method at one level of precision is consistent with the answer given
        // at the next higher level of precision.
        val stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val x0 = choosePoint()
            val x1 = choosePoint()
            if (x0 == -x1) continue;  // Not allowed by API.
            val c0 = (if (oneIn(2)) -1 else 1) * pow(1e-20, randomDouble())
            val c1 = (if (oneIn(2)) -1 else 1) * pow(1e-20, randomDouble())
            val z = (c0 * x0 + c1 * x1).normalize()
            val r = S1Angle.radians(M_PI * pow(1e-30, randomDouble()))
            val a = S2EdgeDistances.interpolateAtDistance(r, z, choosePoint())
            val b = S2EdgeDistances.interpolateAtDistance(r, z, choosePoint())
            val c = S2EdgeDistances.interpolateAtDistance(r, z, choosePoint())
            val prec = testEdgeCircumcenterSignConsistency(x0, x1, a, b, c)
            // Don't skew the statistics by recording degenerate inputs.
            if (x0 == x1) {
                // This precision would be SYMBOLIC if we handled this degeneracy.
                assertThat(prec).isEqualTo(EXACT)
            } else if (a == b || b == c || c == a) {
                assertThat(prec).isEqualTo(PrecisionStats.Precision.SYMBOLIC)
            } else {
                stats.tally(prec)
            }
        }
        logger.info { stats.toString() }
    }

    // Verifies that VoronoiSiteExclusion(a, b, x0, x1, r) == expected_result, and
// furthermore checks that the minimum required precision is "expected_prec".
    private fun testVoronoiSiteExclusion(
        a: S2Point,
        b: S2Point,
        x0: S2Point,
        x1: S2Point,
        r: S1ChordAngle,
        expected_result: S2Predicates.Excluded,
        expected_prec: PrecisionStats.Precision
    ) {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!S2PointUtil.isUnitLength(a)) a.normalize()
        if (!S2PointUtil.isUnitLength(b)) b.normalize()
        if (!S2PointUtil.isUnitLength(x0)) x0.normalize()
        if (!S2PointUtil.isUnitLength(x1)) x1.normalize()

        // The internal methods (Triage, Exact, etc) require that site A is closer
        // to X0 and site B is closer to X1.  GetVoronoiSiteExclusion has special
        // code to handle the case where this is not true.  We need to duplicate
        // that code here.  Essentially, since the API requires site A to be closer
        // than site B to X0, then if site A is also closer to X1 then site B must
        // be excluded.
        if (S2Predicates.compareDistances(x1, a, b) < 0) {
            assertThat(SECOND).isEqualTo(expected_result)
            // We don't know what precision was used by CompareDistances(), but we
            // arbitrarily require the test to specify it as DOUBLE.
            assertThat(DOUBLE).isEqualTo(expected_prec)
        } else {
            val dbl_result = S2Predicates.triageVoronoiSiteExclusion(a, b, x0, x1, r.length2)
//            val ld_result = S2Predicates.triageVoronoiSiteExclusion(
//                a.toLongDouble(),
//                b.toLongDouble(),
//                x0.toLongDouble(),
//                x1.toLongDouble(),
//                r.length2.toLD()
//            )
            val exact_result = S2Predicates.exactVoronoiSiteExclusion(
                a.toExactFloat(),
                b.toExactFloat(),
                x0.toExactFloat(),
                x1.toExactFloat(),
                r.length2.toExactFloat()
            )

            // Check that the results are correct (if not UNCERTAIN), and also that if
            // dbl_result is not UNCERTAIN then so is ld_result, etc.
            assertThat(exact_result).isEqualTo(expected_result)
            //if (ld_result != UNCERTAIN) assertThat(ld_result).isEqualTo(exact_result)
            //if (dbl_result != UNCERTAIN) assertThat(dbl_result).isEqualTo(ld_result)
            if (dbl_result != UNCERTAIN) assertThat(dbl_result).isEqualTo(exact_result)

            //val actual_prec = if (dbl_result != UNCERTAIN) DOUBLE else if (ld_result != UNCERTAIN) LONG_DOUBLE else EXACT
            val actual_prec = if (dbl_result != UNCERTAIN) DOUBLE else EXACT
            assertThat(actual_prec).isEqualTo(expected_prec)
        }
        // Make sure that the top-level function returns the expected result.
        assertThat(S2Predicates.getVoronoiSiteExclusion(a, b, x0, x1, r)).isEqualTo(expected_result)

        // If site B is closer to X1, then the same site should be excluded (if any)
        // when we swap the sites and the edge direction.
        val swapped_result = if (expected_result == FIRST) SECOND
        else if (expected_result == SECOND) FIRST
        else expected_result
        if (S2Predicates.compareDistances(x1, b, a) < 0) {
            assertThat(S2Predicates.getVoronoiSiteExclusion(b, a, x1, x0, r)).isEqualTo(swapped_result)
        }
    }

    @Test
    fun voronoiSiteExclusionCoverage() {
        // Both sites are closest to edge endpoint X0.
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-5, 0.0),
            S2Point(1.0, -2e-5, 0.0),
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-3),
            SECOND,
            DOUBLE
        )

        // Both sites are closest to edge endpoint X1.
        testVoronoiSiteExclusion(
            S2Point(1.0, 1.0, 1e-30),
            S2Point(1.0, 1.0, -1e-20),
            S2Point(1, 0, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-10),
            SECOND,
            DOUBLE
        )

        // Test cases where neither site is excluded.
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-10, 1e-5),
            S2Point(1.0, 1e-10, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-4),
            NEITHER,
            DOUBLE
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-10, 1e-5),
            S2Point(1.0, 1e-10, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-5),
            NEITHER,
            EXACT
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-17, 1e-5),
            S2Point(1.0, 1e-17, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-4),
            NEITHER,
            EXACT
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-20, 1e-5),
            S2Point(1.0, 1e-20, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1e-5),
            NEITHER,
            EXACT
        )

        // Test cases where the first site is excluded.  (Tests where the second
        // site is excluded are constructed by TestVoronoiSiteExclusion.)
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-6, 1.0049999999e-5),
            S2Point(1.0, 0.0, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1.005e-5),
            FIRST,
            DOUBLE
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1.00105e-6, 1.0049999999e-5),
            S2Point(1.0, 0.0, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1.005e-5),
            FIRST,
            EXACT
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-6, 1.005e-5),
            S2Point(1.0, 0.0, -1e-5),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1.005e-5),
            FIRST,
            EXACT
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-31, 1.005e-30),
            S2Point(1.0, 0.0, -1e-30),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1.005e-30),
            FIRST,
            EXACT
        )
        testVoronoiSiteExclusion(
            S2Point(1.0, -1e-31, 1.005e-30),
            S2Point(1.0, 0.0, -1e-30),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S1ChordAngle.radians(1.005e-30),
            FIRST,
            EXACT
        )

        // These two sites are exactly 60 degrees away from the point (1, 1, 0),
        // which is the midpoint of edge X.  This case requires symbolic
        // perturbations to resolve correctly.  Site A is closer to every point in
        // its coverage interval except for (1, 1, 0), but site B is considered
        // closer to that point symbolically.
        testVoronoiSiteExclusion(
            S2Point(0, 1, 1),
            S2Point(1, 0, 1),
            S2Point(0, 1, 1),
            S2Point(1, 0, -1),
            S1ChordAngle.fromLength2(1),
            NEITHER,
            EXACT
        )

        // This test is similar except that site A is considered closer to the
        // equidistant point (-1, 1, 0), and therefore site B is excluded.
        testVoronoiSiteExclusion(
            S2Point(0, 1, 1),
            S2Point(-1, 0, 1),
            S2Point(0, 1, 1),
            S2Point(-1, 0, -1),
            S1ChordAngle.fromLength2(1),
            SECOND,
            EXACT
        )
    }

    // Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
    fun testVoronoiSiteExclusionConsistency(
        a: S2Point,
        b: S2Point,
        x0: S2Point,
        x1: S2Point,
        r: S1ChordAngle
    ): PrecisionStats.Precision {
        // The internal methods require this (see TestVoronoiSiteExclusion).
        if (S2Predicates.compareDistances(x1, a, b) < 0) return DOUBLE

        val dbl_result = S2Predicates.triageVoronoiSiteExclusion(a, b, x0, x1, r.length2)
//        val ld_result = S2Predicates.triageVoronoiSiteExclusion(
//            a.toLongDouble(),
//            b.toLongDouble(),
//            x0.toLongDouble(),
//            x1.toLongDouble(),
//            r.length2.toLD()
//        )
        val exact_result = S2Predicates.exactVoronoiSiteExclusion(
            a.toExactFloat(),
            b.toExactFloat(),
            x0.toExactFloat(),
            x1.toExactFloat(),
            r.length2.toExactFloat()
        )
        assertThat(S2Predicates.getVoronoiSiteExclusion(a, b, x0, x1, r)).isEqualTo(exact_result)

        assertThat(UNCERTAIN != exact_result).isTrue()
//        if (ld_result == UNCERTAIN) {
//            assertThat(dbl_result).isEqualTo(UNCERTAIN)
//            return EXACT
//        }
//        assertThat(ld_result).isEqualTo(exact_result)
        if (dbl_result == UNCERTAIN) {
            return EXACT
        }
        assertThat(dbl_result).isEqualTo(exact_result)
        return DOUBLE
    }

    @Test
    fun voronoiSiteExclusionConsistency() {
        // This test chooses random a random edge X, a random point P on that edge,
        // and a random threshold distance "r".  It then choose two sites A and B
        // whose distance to P is almost exactly "r".  This ensures that the
        // coverage intervals for A and B will (almost) share a common endpoint.  It
        // then checks that the answer given by a method at one level of precision
        // is consistent with the answer given at higher levels of precision.
        val stats = PrecisionStats()
        for (iter in 0 until kConsistencyIters) {
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val x0 = choosePoint()
            val x1 = choosePoint()
            if (x0 == -x1) continue;  // Not allowed by API.
            val f = pow(1e-20, randomDouble())
            val p = ((1 - f) * x0 + f * x1).normalize()
            val r1 = S1Angle.radians(M_PI_2 * 1e-20.pow(randomDouble()))
            val a = S2EdgeDistances.interpolateAtDistance(r1, p, choosePoint())
            val b = S2EdgeDistances.interpolateAtDistance(r1, p, choosePoint())
            // Check that the other API requirements are met.
            val r = S1ChordAngle(r1)
            if (S2Predicates.compareEdgeDistance(a, x0, x1, r) > 0) continue
            if (S2Predicates.compareEdgeDistance(b, x0, x1, r) > 0) continue
            if (S2Predicates.compareDistances(x0, a, b) > 0) {
                var temp = a; }
            if (a == b) continue

            val prec = testVoronoiSiteExclusionConsistency(a, b, x0, x1, r)
            // Don't skew the statistics by recording degenerate inputs.
            if (x0 == x1) {
                assertThat(prec).isEqualTo(DOUBLE)
            } else {
                stats.tally(prec)
            }
        }
        logger.info { stats.toString() }
    }

    @Test
    fun sign() {
        var a = S2Point(1.0, 0.0, 0.0)
        var b = S2Point(0.0, 1.0, 0.0)
        val c = S2Point(0.0, 0.9999999999999999, 0.0)
        assertThat(S2Predicates.sign(a, b, c)).isOne()

        val deltaBCrossC = S2Point(
            7.8615562344198168526508285956e-6,
            6.24999561121962149754537261329e-5,
            -0.000143766135723317032302845730007
        ) -
                S2Point(
                    0.000007861556234467316471419675022,
                    0.000062499956112179674477218761648,
                    -0.00014376613572334472465439473944
                )
        assertThat(
            S2Predicates.exactSign(
                a = S2Point(0.4903891754856072, 0.7891211432339005, 0.36987332678603396),
                b = S2Point(0.4905257354088808, 0.7890477243764918, 0.3698488766480931),
                c = S2Point(0.49025260348073774, 0.7891945426499283, 0.3698977678114979),
                perturb = true
            )
        ).isOne()

        val x = S2Point(-2.368580570649774E-13, -1.0, 1.9413907697072796E-11)
        a = S2Point(-2.368580570649762E-13, -1.0, 1.9413907697072796E-11)
        b = S2Point(-2.3685805706497854E-13, -1.0, 1.9413907697072796E-11)
        assertThat(S2Predicates.compareSin2Distances(x, a, b)).isZero()
        assertThat(S2Predicates.compareDistances(x, a, b)).isOne()
    }

    @Test
    fun compDistance() {
        assertThat(
            S2Predicates.compareDistances(
                x = S2Point(-0.3889376284948148, -0.921264088706832, 4.3101474383497625E-9),
                a = S2Point(-0.3889376284948148, -0.921264088706832, 4.3101474382150685E-9),
                b = S2Point(-0.3889376284948148, -0.921264088706832, 4.3101474384844565E-9)
            )
        ).isOne()
    }
}
