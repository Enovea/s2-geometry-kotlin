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

import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.M_SQRT2
import dilivia.math.M_SQRT3
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random.oneIn
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomPoint
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import java.lang.Math.pow
import kotlin.math.acos
import kotlin.math.asin
import kotlin.math.pow
import kotlin.math.sqrt


class S2EdgeDistancesUnitTest {

    // Checks that the error returned by S2::GetUpdateMinDistanceMaxError() for
    // the distance "input" (measured in radians) corresponds to a distance error
    // of less than "max_error" (measured in radians).
    //
    // The reason for the awkward phraseology above is that the value returned by
    // GetUpdateMinDistanceMaxError() is not a distance; it represents an error in
    // the *squared* distance.
    fun checkUpdateMinDistanceMaxError(actual: Double, maxError: Double) {
        val ca = S1ChordAngle(S1Angle.radians(actual))
        val bound = ca.plusError(S2EdgeDistances.getUpdateMinDistanceMaxError(ca)).toAngle()
        assertThat(bound.radians - actual).isLessThanOrEqualTo(maxError)
    }

    @Test
    fun getUpdateMinDistanceMaxError() {
        // Verify that the error is "reasonable" for a sampling of distances.
        checkUpdateMinDistanceMaxError(0.0, 1.5e-15)
        checkUpdateMinDistanceMaxError(1e-8, 1e-15)
        checkUpdateMinDistanceMaxError(1e-5, 1e-15)
        checkUpdateMinDistanceMaxError(0.05, 1e-15)
        checkUpdateMinDistanceMaxError(M_PI_2 - 1e-8, 2e-15)
        checkUpdateMinDistanceMaxError(M_PI_2, 2e-15)
        checkUpdateMinDistanceMaxError(M_PI_2 + 1e-8, 2e-15)
        checkUpdateMinDistanceMaxError(M_PI - 1e-5, 2e-10)
        checkUpdateMinDistanceMaxError(M_PI, 0.0)
    }

    @Test
    fun getUpdateMinInteriorDistanceMaxError() {
        // Check that the error bound returned by
        // GetUpdateMinInteriorDistanceMaxError() is large enough.
        for (iter in 0 until 10000) {
            val a0 = randomPoint()
            val len = S1Angle.radians(M_PI * 1e-20.pow(randomDouble()))
            val a1 = S2EdgeDistances.interpolateAtDistance(len, a0, randomPoint())
            // TODO(ericv): If s2pred::RobustCrossProd() is implemented, then we can
            // also test nearly-antipodal points here.  In theory the error bound can
            // be exceeded when the edge endpoints are antipodal to within 0.8e-13
            // radians, but the only examples found in testing require the endpoints
            // to be nearly-antipodal to within 1e-16 radians.
            val n = S2PointUtil.robustCrossProd(a0, a1).normalize()
            val f = pow(1e-20, randomDouble())
            val a = ((1 - f) * a0 + f * a1).normalize()
            var r = S1Angle.radians(M_PI_2 * pow(1e-20, randomDouble()))
            if (oneIn(2)) r = S1Angle.radians(M_PI_2) - r
            val x = S2EdgeDistances.interpolateAtDistance(r, a, n)
            val minDist = S1ChordAngle.infinity().clone()
            if (!S2EdgeDistances.updateMinInteriorDistance(x, a0, a1, minDist)) {
                continue
            }
            val error = S2EdgeDistances.getUpdateMinDistanceMaxError(minDist)
            assertThat(S2Predicates.compareEdgeDistance(x, a0, a1, minDist.plusError(error))).isLessThanOrEqualTo(0)
            assertThat(S2Predicates.compareEdgeDistance(x, a0, a1, minDist.plusError(-error))).isGreaterThanOrEqualTo(0)
        }
    }

    // Given a point X and an edge AB, check that the distance from X to AB is
    // "distance_radians" and the closest point on AB is "expected_closest".
    fun checkDistance(x: S2Point, a: S2Point, b: S2Point, distance_radians: Double, expected_closest: S2Point) {
        x.normalize()
        a.normalize()
        b.normalize()
        val expectedClosest = if (expected_closest.norm2() != 0.0) expected_closest.normalize() else expected_closest
        assertThat(S2EdgeDistances.getDistance(x, a, b).radians).isCloseTo(distance_radians, Offset.offset(1e-15))
        val closest = S2EdgeDistances.project(x, a, b)
        if (expectedClosest == S2Point(0, 0, 0)) {
            // This special value says that the result should be A or B.
            assertThat(closest == a || closest == b).isTrue()
        } else {
            assertThat(S2PointUtil.approxEquals(closest, expectedClosest)).isTrue()
        }
        val minDistance = S1ChordAngle.zero().clone()
        assertThat(S2EdgeDistances.updateMinDistance(x, a, b, minDistance)).isFalse()
        minDistance.length2 = S1ChordAngle.infinity().length2
        assertThat(S2EdgeDistances.updateMinDistance(x, a, b, minDistance)).isTrue()
        assertThat(minDistance.toAngle().radians).isCloseTo(distance_radians, Offset.offset(1e-15))
    }

    @Test
    fun distance() {
        checkDistance(S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), 0.0, S2Point(1, 0, 0))
        checkDistance(S2Point(0, 1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), 0.0, S2Point(0, 1, 0))
        checkDistance(S2Point(1, 3, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), 0.0, S2Point(1, 3, 0))
        checkDistance(S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2, S2Point(1, 0, 0))
        checkDistance(S2Point(0, 0, -1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2, S2Point(1, 0, 0))
        checkDistance(S2Point(-1, -1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), 0.75 * M_PI, S2Point(0, 0, 0))

        checkDistance(S2Point(0, 1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0), M_PI_4, S2Point(1, 1, 0))
        checkDistance(S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0), M_PI_2, S2Point(1, 0, 0))

        checkDistance(S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(-1, 1, 0), M_PI_2, S2Point(1, 0, 0))
        checkDistance(S2Point(-1, -1, 0), S2Point(1, 0, 0), S2Point(-1, 1, 0), M_PI_2, S2Point(-1, 1, 0))

        checkDistance(S2Point(1, 1, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), asin(sqrt(1.0 / 3.0)), S2Point(1, 1, 0))
        checkDistance(S2Point(1, 1, -1), S2Point(1, 0, 0), S2Point(0, 1, 0), asin(sqrt(1.0 / 3.0)), S2Point(1, 1, 0))

        checkDistance(S2Point(-1, 0, 0), S2Point(1, 1, 0), S2Point(1, 1, 0), 0.75 * M_PI, S2Point(1, 1, 0))
        checkDistance(S2Point(0, 0, -1), S2Point(1, 1, 0), S2Point(1, 1, 0), M_PI_2, S2Point(1, 1, 0))
        checkDistance(S2Point(-1, 0, 0), S2Point(1, 0, 0), S2Point(1, 0, 0), M_PI, S2Point(1, 0, 0))
    }

    @Test
    fun distanceOptimizationIsConservative() {
        // Verifies that AlwaysUpdateMinInteriorDistance() computes the lower bound
        // on the true distance conservatively.  (This test used to fail.)
        val x = S2Point(-0.017952729194524016, -0.30232422079175203, 0.95303607751077712)
        val a = S2Point(-0.017894725505830295, -0.30229974986194175, 0.95304493075220664)
        val b = S2Point(-0.017986591360900289, -0.30233851195954353, 0.95303090543659963)
        var minDistance = S1ChordAngle.infinity().clone()
        assertThat(S2EdgeDistances.updateMinDistance(x, a, b, minDistance)).isTrue()
        minDistance = minDistance.successor().clone()
        assertThat(S2EdgeDistances.updateMinDistance(x, a, b, minDistance)).isTrue()
    }

    fun checkMaxDistance(x: S2Point, a: S2Point, b: S2Point, distance_radians: Double) {
        x.normalize()
        a.normalize()
        b.normalize()
        var maxDistance = S1ChordAngle.straight().clone()
        assertThat(S2EdgeDistances.updateMaxDistance(x, a, b, maxDistance)).isFalse()
        maxDistance = S1ChordAngle.negative().clone()
        assertThat(S2EdgeDistances.updateMaxDistance(x, a, b, maxDistance)).isTrue()

        assertThat(maxDistance.radians()).isCloseTo(distance_radians, Offset.offset(1e-15))
    }

    @Test
    fun maxDistance() {
        checkMaxDistance(S2Point(1, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2)
        checkMaxDistance(S2Point(1, 0, -1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2)
        checkMaxDistance(S2Point(0, 1, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2)
        checkMaxDistance(S2Point(0, 1, -1), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2)

        checkMaxDistance(S2Point(1, 1, 1), S2Point(1, 0, 0), S2Point(0, 1, 0), asin(sqrt(2.0 / 3.0)))
        checkMaxDistance(S2Point(1, 1, -1), S2Point(1, 0, 0), S2Point(0, 1, 0), asin(sqrt(2.0 / 3.0)))

        checkMaxDistance(S2Point(1, 0, 0), S2Point(1, 1, 0), S2Point(1, -1, 0), M_PI_4)
        checkMaxDistance(S2Point(0, 1, 0), S2Point(1, 1, 0), S2Point(-1, 1, 0), M_PI_4)
        checkMaxDistance(S2Point(0, 0, 1), S2Point(0, 1, 1), S2Point(0, -1, 1), M_PI_4)

        checkMaxDistance(S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(1, 0, -1), 3 * M_PI_4)
        checkMaxDistance(S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(1.0, 1.0, -M_SQRT2), 3 * M_PI_4)

        checkMaxDistance(S2Point(0, 0, 1), S2Point(0, 0, -1), S2Point(0, 0, -1), M_PI)
    }

    fun checkInterpolate(t: Double, a: S2Point, b: S2Point, expected: S2Point) {
        a.normalize()
        b.normalize()
        val actual = S2EdgeDistances.interpolate(t, a, b)

        // We allow a bit more than the usual 1e-15 error tolerance because
        // Interpolate() uses trig functions.
        assertThat(S2PointUtil.approxEquals(expected, actual, S1Angle.radians(3e-15))).isTrue()
    }

    @Test
    fun interpolate() {
        // Choose test points designed to expose floating-point errors.
        val p1 = S2Point(0.1, 1e-30, 0.3).normalize()
        val p2 = S2Point(-0.7, -0.55, -1e30).normalize()

        // A zero-length edge.
        checkInterpolate(0.0, p1, p1, p1)
        checkInterpolate(1.0, p1, p1, p1)

        // Start, end, and middle of a medium-length edge.
        checkInterpolate(0.0, p1, p2, p1)
        checkInterpolate(1.0, p1, p2, p2)
        checkInterpolate(0.5, p1, p2, 0.5 * (p1 + p2))

        // Test that interpolation is done using distances on the sphere rather than
        // linear distances.
        checkInterpolate(1.0 / 3.0, S2Point(1, 0, 0), S2Point(0, 1, 0), S2Point(M_SQRT3, 1.0, 0.0))
        checkInterpolate(2.0 / 3.0, S2Point(1, 0, 0), S2Point(0, 1, 0), S2Point(1.0, M_SQRT3, 0.0))

        // Test that interpolation is accurate on a long edge (but not so long that
        // the definition of the edge itself becomes too unstable).

        val kLng = M_PI - 1e-2
        val a = S2LatLng.fromRadians(0.0, 0.0).toPoint()
        val b = S2LatLng.fromRadians(0.0, kLng).toPoint()
        var f = 0.4
        while (f > 1e-15) {
            checkInterpolate(f, a, b, S2LatLng.fromRadians(0.0, f * kLng).toPoint())
            checkInterpolate(1 - f, a, b, S2LatLng.fromRadians(0.0, (1 - f) * kLng).toPoint())
            f *= 0.1
        }

        // Test that interpolation on a 180 degree edge (antipodal endpoints) yields
        // a result with the correct distance from each endpoint.
        var t = 0.0
        while (t <= 1) {
            val actual = S2EdgeDistances.interpolate(t, p1, -p1)
            assertThat(S1Angle(actual, p1).radians).isCloseTo(t * M_PI, Offset.offset(3e-15))
            t += 0.125
        }
    }

    @Test
    fun interpolateCanExtrapolate() {
        val i = S2Point(1, 0, 0)
        val j = S2Point(0, 1, 0)
        // Initial vectors at 90 degrees.
        checkInterpolate(0.0, i, j, S2Point(1, 0, 0))
        checkInterpolate(1.0, i, j, S2Point(0, 1, 0))
        checkInterpolate(1.5, i, j, S2Point(-1, 1, 0))
        checkInterpolate(2.0, i, j, S2Point(-1, 0, 0))
        checkInterpolate(3.0, i, j, S2Point(0, -1, 0))
        checkInterpolate(4.0, i, j, S2Point(1, 0, 0))

        // Negative values of t.
        checkInterpolate(-1.0, i, j, S2Point(0, -1, 0))
        checkInterpolate(-2.0, i, j, S2Point(-1, 0, 0))
        checkInterpolate(-3.0, i, j, S2Point(0, 1, 0))
        checkInterpolate(-4.0, i, j, S2Point(1, 0, 0))

        // Initial vectors at 45 degrees.
        checkInterpolate(2.0, i, S2Point(1, 1, 0), S2Point(0, 1, 0))
        checkInterpolate(3.0, i, S2Point(1, 1, 0), S2Point(-1, 1, 0))
        checkInterpolate(4.0, i, S2Point(1, 1, 0), S2Point(-1, 0, 0))

        // Initial vectors at 135 degrees.
        checkInterpolate(2.0, i, S2Point(-1, 1, 0), S2Point(0, -1, 0))

        // Take a small fraction along the curve.
        val p = S2EdgeDistances.interpolate(0.001, i, j)
        // We should get back where we started.
        checkInterpolate(1000.0, i, p, j)
    }

    @Test
    fun repeatedInterpolation() {
        // Check that points do not drift away from unit length when repeated
        // interpolations are done.
        for (i in 0 until 100) {
            var a = randomPoint()
            val b = randomPoint()
            for (j in 0 until 1000) {
                a = S2EdgeDistances.interpolate(0.01, a, b)
            }
            assertThat(isUnitLength(a)).isTrue()
        }
    }

    // Given two edges a0a1 and b0b1, check that the minimum distance between them
    // is "distance_radians", and that GetEdgePairClosestPoints() returns
    // "expected_a" and "expected_b" as the points that achieve this distance.
    // S2Point(0, 0, 0) may be passed for "expected_a" or "expected_b" to indicate
    // that both endpoints of the corresponding edge are equally distant, and
    // therefore either one might be returned.
    //
    // Parameters are passed by value so that this function can normalize them.
    fun checkEdgePairMinDistance(
        a0: S2Point,
        a1: S2Point,
        b0: S2Point,
        b1: S2Point,
        distance_radians: Double,
        expectedA: S2Point,
        expectedB: S2Point
    ) {
        a0.normalize()
        a1.normalize()
        b0.normalize()
        b1.normalize()
        if (expectedA.norm2() != 0.0) expectedA.normalize()
        if (expectedB.norm2() != 0.0) expectedB.normalize()
        val closest = S2EdgeDistances.getEdgePairClosestPoints(a0, a1, b0, b1)
        val actualA = closest.first
        val actualB = closest.second
        if (expectedA == S2Point(0, 0, 0)) {
            // This special value says that the result should be a0 or a1.
            assertThat(actualA == a0 || actualA == a1).isTrue()
        } else {
            assertThat(S2PointUtil.approxEquals(expectedA, actualA)).isTrue()
        }
        if (expectedB == S2Point(0, 0, 0)) {
            // This special value says that the result should be b0 or b1.
            assertThat(actualB == b0 || actualB == b1).isTrue()
        } else {
            assertThat(S2PointUtil.approxEquals(expectedB, actualB)).isTrue()
        }
        var minDistance = S1ChordAngle.zero().clone()
        assertThat(S2EdgeDistances.updateEdgePairMinDistance(a0, a1, b0, b1, minDistance)).isFalse()
        minDistance = S1ChordAngle.infinity().clone()
        assertThat(S2EdgeDistances.updateEdgePairMinDistance(a0, a1, b0, b1, minDistance)).isTrue()
        assertThat(minDistance.radians()).isCloseTo(distance_radians, Offset.offset(1e-15))
    }

    @Test
    fun edgePairMinDistance() {
        // One edge is degenerate.
        checkEdgePairMinDistance(
            S2Point(1, 0, 1),
            S2Point(1, 0, 1),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            M_PI_4,
            S2Point(1, 0, 1),
            S2Point(1, 0, 0)
        )
        checkEdgePairMinDistance(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1, 0, 1),
            S2Point(1, 0, 1),
            M_PI_4,
            S2Point(1, 0, 0),
            S2Point(1, 0, 1)
        )

        // Both edges are degenerate.
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(0, 1, 0),
            M_PI_2,
            S2Point(1, 0, 0),
            S2Point(0, 1, 0)
        )

        // Both edges are degenerate and antipodal.
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            S2Point(-1, 0, 0),
            S2Point(-1, 0, 0),
            M_PI,
            S2Point(1, 0, 0),
            S2Point(-1, 0, 0)
        )

        // Two identical edges.
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            0.0,
            S2Point(0, 0, 0),
            S2Point(0, 0, 0)
        )

        // Both edges are degenerate and identical.
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            S2Point(1, 0, 0),
            0.0,
            S2Point(1, 0, 0),
            S2Point(1, 0, 0)
        )

        // Edges that share exactly one vertex (all 4 possibilities).
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(0, 1, 0),
            S2Point(0, 1, 1),
            0.0,
            S2Point(0, 1, 0),
            S2Point(0, 1, 0)
        )
        checkEdgePairMinDistance(
            S2Point(0, 1, 0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(0, 1, 1),
            0.0,
            S2Point(0, 1, 0),
            S2Point(0, 1, 0)
        )
        checkEdgePairMinDistance(
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(0, 1, 1),
            S2Point(0, 1, 0),
            0.0,
            S2Point(0, 1, 0),
            S2Point(0, 1, 0)
        )
        checkEdgePairMinDistance(
            S2Point(0, 1, 0),
            S2Point(1, 0, 0),
            S2Point(0, 1, 1),
            S2Point(0, 1, 0),
            0.0,
            S2Point(0, 1, 0),
            S2Point(0, 1, 0)
        )

        // Two edges whose interiors cross.
        checkEdgePairMinDistance(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(1, 0, -1),
            S2Point(1, 0, 1),
            0.0,
            S2Point(1, 0, 0),
            S2Point(1, 0, 0)
        )

        // The closest distance occurs between two edge endpoints, but more than one
        // endpoint pair is equally distant.
        checkEdgePairMinDistance(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(-1, 0, 0),
            S2Point(-1, 0, 1),
            acos(-0.5),
            S2Point(0, 0, 0),
            S2Point(-1, 0, 1)
        )
        checkEdgePairMinDistance(
            S2Point(-1, 0, 0),
            S2Point(-1, 0, 1),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            acos(-0.5),
            S2Point(-1, 0, 1),
            S2Point(0, 0, 0)
        )
        checkEdgePairMinDistance(
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            S2Point(-1, 0, -1),
            S2Point(-1, 0, 1),
            acos(-0.5),
            S2Point(0, 0, 0),
            S2Point(0, 0, 0)
        )
    }

    // Given two edges a0a1 and b0b1, check that the maximum distance between them
    // is "distance_radians".  Parameters are passed by value so that this
    // function can normalize them.
    fun checkEdgePairMaxDistance(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point, distanceRadians: Double) {
        a0.normalize()
        a1.normalize()
        b0.normalize()
        b1.normalize()
        var maxDistance = S1ChordAngle.straight().clone()
        assertThat(S2EdgeDistances.updateEdgePairMaxDistance(a0, a1, b0, b1, maxDistance)).isFalse()
        maxDistance = S1ChordAngle.negative().clone()
        assertThat(S2EdgeDistances.updateEdgePairMaxDistance(a0, a1, b0, b1, maxDistance)).isTrue()
        assertThat(maxDistance.radians()).isCloseTo(distanceRadians, Offset.offset(1e-15))
    }

    @Test
    fun edgePairMaxDistance() {
        // Standard situation.  Same hemisphere, not degenerate.
        checkEdgePairMaxDistance(
            S2Point(1, 0, 0),
            S2Point(0, 1, 0),
            S2Point(1, 1, 0),
            S2Point(1, 1, 1),
            acos(1 / M_SQRT3)
        )

        // One edge is degenerate.
        checkEdgePairMaxDistance(S2Point(1, 0, 1), S2Point(1, 0, 1), S2Point(1, -1, 0), S2Point(1, 1, 0), acos(0.5))
        checkEdgePairMaxDistance(S2Point(1, -1, 0), S2Point(1, 1, 0), S2Point(1, 0, 1), S2Point(1, 0, 1), acos(0.5))

        // Both edges are degenerate.
        checkEdgePairMaxDistance(S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), S2Point(0, 1, 0), M_PI_2)

        // Both edges are degenerate and antipodal.
        checkEdgePairMaxDistance(S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(-1, 0, 0), S2Point(-1, 0, 0), M_PI)

        // Two identical edges.
        checkEdgePairMaxDistance(S2Point(1, 0, 0), S2Point(0, 1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0), M_PI_2)

        // Both edges are degenerate and identical.
        checkEdgePairMaxDistance(S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(1, 0, 0), 0.0)

        // Antipodal reflection of one edge crosses the other edge.
        checkEdgePairMaxDistance(S2Point(1, 0, 1), S2Point(1, 0, -1), S2Point(-1, -1, 0), S2Point(-1, 1, 0), M_PI)

        // One vertex of one edge touches the interior of the antipodal reflection
        // of the other edge.
        checkEdgePairMaxDistance(S2Point(1, 0, 1), S2Point(1, 0, 0), S2Point(-1, -1, 0), S2Point(-1, 1, 0), M_PI)
    }

    fun isEdgeBNearEdgeA(a_str: String, b_str: String, max_error_degrees: Double): Boolean {
        val a = makePolyline(a_str)
        assertThat(a.numVertices).isEqualTo(2)
        val b = makePolyline(b_str)
        assertThat(b.numVertices).isEqualTo(2)
        return S2EdgeDistances.isEdgeBNearEdgeA(
            a.vertex(0),
            a.vertex(1),
            b.vertex(0),
            b.vertex(1),
            S1Angle.degrees(max_error_degrees)
        )
    }

    @Test
    fun edgeBNearEdgeA() {
        // Edge is near itself.
        assertThat(isEdgeBNearEdgeA("5:5, 10:-5", "5:5, 10:-5", 1e-6)).isTrue()

        // Edge is near its reverse
        assertThat(isEdgeBNearEdgeA("5:5, 10:-5", "10:-5, 5:5", 1e-6)).isTrue()

        // Short edge is near long edge.
        assertThat(isEdgeBNearEdgeA("10:0, -10:0", "2:1, -2:1", 1.0)).isTrue()

        // Long edges cannot be near shorter edges.
        assertThat(isEdgeBNearEdgeA("2:1, -2:1", "10:0, -10:0", 1.0)).isFalse()

        // Orthogonal crossing edges are not near each other...
        assertThat(isEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 1.0)).isFalse()

        // ... unless all points on B are within tolerance of A.
        assertThat(isEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 2.0)).isTrue()

        // Very long edges whose endpoints are close may have interior points that are
        // far apart.  An implementation that only considers the vertices of polylines
        // will incorrectly consider such edges as "close" when they are not.
        // Consider, for example, two consecutive lines of longitude.  As they
        // approach the poles, they become arbitrarily close together, but along the
        // equator they bow apart.
        assertThat(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 0.5)).isFalse()
        assertThat(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 1.5)).isTrue()

        // The two arcs here are nearly as long as S2 edges can be (just shy of 180
        // degrees), and their endpoints are less than 1 degree apart.  Their
        // midpoints, however, are at opposite ends of the sphere along its equator.
        assertThat(
            isEdgeBNearEdgeA(
                "0:-179.75, 0:-0.25", "0:179.75, 0:0.25", 1.0
            )
        ).isFalse()

        // At the equator, the second arc here is 9.75 degrees from the first, and
        // closer at all other points.  However, the southern point of the second arc
        // (-1, 9.75) is too far from the first arc for the short-circuiting logic in
        // isEdgeBNearEdgeA to apply.
        assertThat(isEdgeBNearEdgeA("40:0, -5:0", "39:0.975, -1:0.975", 1.0)).isTrue()

        // Same as above, but B's orientation is reversed, causing the angle between
        // the normal vectors of circ(B) and circ(A) to be (180-9.75) = 170.5 degrees,
        // not 9.75 degrees.  The greatest separation between the planes is still 9.75
        // degrees.
        assertThat(isEdgeBNearEdgeA("10:0, -10:0", "-.4:0.975, 0.4:0.975", 1.0)).isTrue()

        // A and B are on the same great circle, A and B partially overlap, but the
        // only part of B that does not overlap A is shorter than tolerance.
        assertThat(isEdgeBNearEdgeA("0:0, 1:0", "0.9:0, 1.1:0", 0.25)).isTrue()

        // A and B are on the same great circle, all points on B are close to A at its
        // second endpoint, (1,0).
        assertThat(isEdgeBNearEdgeA("0:0, 1:0", "1.1:0, 1.2:0", 0.25)).isTrue()

        // Same as above, but B's orientation is reversed.  This case is special
        // because the projection of the normal defining A onto the plane containing B
        // is the null vector, and must be handled by a special case.
        assertThat(isEdgeBNearEdgeA("0:0, 1:0", "1.2:0, 1.1:0", 0.25)).isTrue()
    }

}
