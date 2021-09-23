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

import dilivia.math.DoubleType
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random
import dilivia.s2.edge.S2EdgeCrossings.intersectionMethodCounter
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.lang.Math.pow
import kotlin.math.pow


class S2EdgeCrossingsUnitTest {

    private val logger = KotlinLogging.logger { }

    fun printIntersectionStats() {
        var total = 0
        val totals = IntArray(S2EdgeCrossings.IntersectionMethod.values().size)
        for (i in S2EdgeCrossings.IntersectionMethod.values()) {
            total += intersectionMethodCounter[i] ?: 0
            totals[i.ordinal] = total
        }
        println("%10s %16s %16s  %6s".format("Method", "Successes", "Attempts", "Rate"))
        for (i in S2EdgeCrossings.IntersectionMethod.values()) {
            if (intersectionMethodCounter[i] == null || intersectionMethodCounter[i] == 0) continue
            println(
                "%10s %9d %5.1f%% %9d %5.1f%%  %5.1f%%".format(
                    i.methodName,
                    intersectionMethodCounter[i],
                    100.0 * intersectionMethodCounter[i]!! / total,
                    totals[i.ordinal],
                    100.0 * totals[i.ordinal] / total,
                    100.0 * intersectionMethodCounter[i]!! / totals[i.ordinal]
                )
            )
        }
        for (i in S2EdgeCrossings.IntersectionMethod.values()) intersectionMethodCounter[i] = 0
    }

    // This returns the true intersection point of two line segments (a0,a1) and
    // (b0,b1), with a relative error of at most DoubleType.epsilon in each coordinate
    // (i.e., one ulp, or twice the double precision rounding error).
    fun getIntersectionExact(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): S2Point {
        var x = S2EdgeCrossings.getIntersectionExact(a0, a1, b0, b1)
        if (x.dotProd((a0 + a1) + (b0 + b1)) < 0) x = -x
        return x
    }

    // The approximate maximum error in GetDistance() for small distances.
    val kGetDistanceAbsError = S1Angle.radians(3 * DoubleType.epsilon)

    @Test
    fun intersectionError() {
        // We repeatedly construct two edges that cross near a random point "p", and
        // measure the distance from the actual intersection point "x" to the
        // exact intersection point and also to the edges.

        var max_point_dist: S1Angle = S1Angle.zero()
        var max_edge_dist: S1Angle = S1Angle.zero()
        for (iter in 0 until 5000) {
            // We construct two edges AB and CD that intersect near "p".  The angle
            // between AB and CD (expressed as a slope) is chosen randomly between
            // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
            // Similarly, two edge lengths approximately between 1e-15 and 1 are
            // chosen.  The edge endpoints are chosen such that they are often very
            // close to the other edge (i.e., barely crossing).  Taken together this
            // ensures that we test both long and very short edges that intersect at
            // both large and very small angles.
            //
            // Sometimes the edges we generate will not actually cross, in which case
            // we simply try again.
            var (p, d1, d2) = S2Random.randomFrame()
            val slope = 1e-15 * 1e30.pow(S2Random.randomDouble())
            d2 = (d1 + slope * d2).normalize()
            var a: S2Point
            var b: S2Point
            var c: S2Point
            var d: S2Point
            var expected: S2Point
            do {
                val abLen = 1e-15.pow(S2Random.randomDouble())
                val cdLen = 1e-15.pow(S2Random.randomDouble())
                var aFraction = 1e-5.pow(S2Random.randomDouble())
                if (S2Random.oneIn(2)) aFraction = 1 - aFraction
                var cFraction = 1e-5.pow(S2Random.randomDouble())
                if (S2Random.oneIn(2)) cFraction = 1 - cFraction
                a = (p - aFraction * abLen * d1).normalize()
                b = (p + (1 - aFraction) * abLen * d1).normalize()
                c = (p - cFraction * cdLen * d2).normalize()
                d = (p + (1 - cFraction) * cdLen * d2).normalize()
                expected = S2EdgeCrossings.getIntersectionExact(a, b, c, d)
            } while (S2EdgeCrossings.crossingSign(a, b, c, d) <= 0
                || S2EdgeDistances.getDistance(
                    expected,
                    a,
                    b
                ) > S1Angle.radians(3 * DoubleType.epsilon) + kGetDistanceAbsError
                || S2EdgeDistances.getDistance(
                    expected,
                    c,
                    d
                ) > S1Angle.radians(3 * DoubleType.epsilon) + kGetDistanceAbsError
            )

            // Each constructed edge should be at most 1.5 * DoubleType.epsilon away from the
            // original point P.
            assertThat(
                S2EdgeDistances.getDistance(
                    p,
                    a,
                    b
                )
            ).isLessThanOrEqualTo(S1Angle.radians(1.5 * DoubleType.epsilon) + kGetDistanceAbsError)
            assertThat(
                S2EdgeDistances.getDistance(
                    p,
                    c,
                    d
                )
            ).isLessThanOrEqualTo(S1Angle.radians(1.5 * DoubleType.epsilon) + kGetDistanceAbsError)

            // Verify that the expected intersection point is close to both edges and
            // also close to the original point P.  (It might not be very close to P
            // if the angle between the edges is very small.)
            assertThat(
                S2EdgeDistances.getDistance(
                    expected,
                    a,
                    b
                )
            ).isLessThanOrEqualTo(S1Angle.radians(3 * DoubleType.epsilon) + kGetDistanceAbsError)
            assertThat(
                S2EdgeDistances.getDistance(
                    expected,
                    c,
                    d
                )
            ).isLessThanOrEqualTo(S1Angle.radians(3 * DoubleType.epsilon) + kGetDistanceAbsError)
            assertThat(
                S1Angle(
                    expected,
                    p
                )
            ).isLessThanOrEqualTo(S1Angle.radians(3 * DoubleType.epsilon / slope) + S2EdgeCrossings.kIntersectionError)

            // Now we actually test the GetIntersection() method.
            val actual = S2EdgeCrossings.getIntersection(a, b, c, d)
            val dist_ab = S2EdgeDistances.getDistance(actual, a, b)
            val dist_cd = S2EdgeDistances.getDistance(actual, c, d)
            assertThat(dist_ab).isLessThanOrEqualTo(S2EdgeCrossings.kIntersectionError + kGetDistanceAbsError)
            assertThat(dist_cd).isLessThanOrEqualTo(S2EdgeCrossings.kIntersectionError + kGetDistanceAbsError)
            max_edge_dist = maxOf(max_edge_dist, maxOf(dist_ab, dist_cd))
            val point_dist = S1Angle(expected, actual)
            assertThat(point_dist).isLessThanOrEqualTo(S2EdgeCrossings.kIntersectionError)
            max_point_dist = maxOf(max_point_dist, point_dist)
        }
        printIntersectionStats()
        logger.info { "Max distance to either edge being intersected: ${max_edge_dist.radians}" }
        logger.info { "Maximum distance to expected intersection point: ${max_point_dist.radians}" }
    }

    // Chooses a point in the XY plane that is separated from X by at least 1e-15
    // (to avoid choosing too many duplicate points) and by at most Pi/2 - 1e-3
    // (to avoid nearly-diametric edges, since the test below is not sophisticated
    // enough to test such edges).
    fun chooseSemicirclePoint(x: S2Point, y: S2Point): S2Point {
        val sign = (2 * S2Random.randomInt(2)) - 1
        return (x + sign * 1e3 * pow(1e-18, S2Random.randomDouble()) * y).normalize()
    }

    @Test
    fun grazingIntersections() {
        // This test choose 5 points along a great circle (i.e., as collinear as
        // possible), and uses them to construct an edge AB and a triangle CDE such
        // that CD and CE both cross AB.  It then checks that the intersection
        // points returned by GetIntersection() have the correct relative ordering
        // along AB (to within kIntersectionError).
        for (iter in 0 until 1000) {
            val (x, y, _) = S2Random.randomFrame()
            var a: S2Point
            var b: S2Point
            var c: S2Point
            var d: S2Point
            var e: S2Point
            var ab: S2Point
            do {
                a = chooseSemicirclePoint(x, y)
                b = chooseSemicirclePoint(x, y)
                c = chooseSemicirclePoint(x, y)
                d = chooseSemicirclePoint(x, y)
                e = chooseSemicirclePoint(x, y)
                ab = (a - b).crossProd(a + b)
            } while (ab.norm() < 50 * DoubleType.epsilon ||
                S2EdgeCrossings.crossingSign(a, b, c, d) <= 0 ||
                S2EdgeCrossings.crossingSign(a, b, c, e) <= 0
            )
            val xcd = S2EdgeCrossings.getIntersection(a, b, c, d)
            val xce = S2EdgeCrossings.getIntersection(a, b, c, e)
            // Essentially this says that if CDE and CAB have the same orientation,
            // then CD and CE should intersect along AB in that order.
            ab = ab.normalize()
            if (S1Angle(xcd, xce) > 2 * S2EdgeCrossings.kIntersectionError) {
                assertThat(S2Predicates.sign(c, d, e) == S2Predicates.sign(c, a, b)).isEqualTo(
                    S2Predicates.sign(
                        ab,
                        xcd,
                        xce
                    ) > 0
                )
            }
        }
        printIntersectionStats()
    }

    @Test
    fun exactIntersectionUnderflow() {
        // Tests that a correct intersection is computed even when two edges are
        // exactly collinear and the normals of both edges underflow in double
        // precision when normalized (see S2PointFromExact function for details).
        val a0 = S2Point(1, 0, 0)
        val a1 = S2Point(1.0, 2e-300, 0.0)
        val b0 = S2Point(1.0, 1e-300, 0.0)
        val b1 = S2Point(1.0, 3e-300, 0.0)
        assertThat(S2EdgeCrossings.getIntersection(a0, a1, b0, b1)).isEqualTo(S2Point(1.0, 1e-300, 0.0))
    }

    @Test
    fun getIntersectionInvariants() {
        // Test that the result of GetIntersection does not change when the edges
        // are swapped and/or reversed.  The number of iterations is high because it
        // is difficult to generate test cases that show that CompareEdges() is
        // necessary and correct, for example.
        val kIters = 5000 // 50000
        for (iter in 0 until kIters) {
            var a: S2Point
            var b: S2Point
            var c: S2Point
            var d: S2Point
            do {
                // GetIntersectionStable() sorts the two edges by length, so construct
                // edges (a,b) and (c,d) that cross and have exactly the same length.
                // This can be done by swapping the "x" and "y" coordinates.
                // [Swapping other coordinate pairs doesn't work because it changes the
                // order of addition in Norm2() == (x**2 + y**2) + z**2.]
                a = S2Random.randomPoint()
                b = S2Random.randomPoint()
                c = S2Point(a[1], a[0], a[2])
                d = S2Point(b[1], b[0], b[2])
            } while (S2EdgeCrossings.crossingSign(a, b, c, d) <= 0)
            assertThat((a - b).norm2()).isEqualTo((c - d).norm2())

            // Now verify that GetIntersection returns exactly the same result when
            // the edges are swapped and/or reversed.
            val result = S2EdgeCrossings.getIntersection(a, b, c, d)
            if (S2Random.oneIn(2)) {
                val temp = a
                a = b
                b = temp
            }
            if (S2Random.oneIn(2)) {
                val temp = c
                c = d
                d = temp
            }
            if (S2Random.oneIn(2)) {
                var temp = a
                a = c
                c = temp
                temp = b
                b = d
                d = temp
            }
            assertThat(result).isEqualTo(S2EdgeCrossings.getIntersection(a, b, c, d))
        }
    }

}




