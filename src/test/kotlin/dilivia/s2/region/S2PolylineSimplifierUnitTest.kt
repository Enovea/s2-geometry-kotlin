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

import dilivia.math.DoubleType
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.edge.S2EdgeDistances
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2PolylineSimplifierUnitTest {

    fun checkSimplify(
        src: String,
        dst: String,
        target: String,
        avoid: String,
        disc_on_left: List<Boolean>,
        radiusDegrees: Double,
        expectedResult: Boolean
    ) {
        val radius = S1ChordAngle(S1Angle.degrees(radiusDegrees))
        val s = S2PolylineSimplifier()
        s.init(S2TextParser.makePoint(src))
        for (p in S2TextParser.parsePoints(target)) {
            s.targetDisc(p, radius)
        }
        var i = 0
        for (p in S2TextParser.parsePoints(avoid)) {
            s.avoidDisc(p, radius, disc_on_left[i++])
        }
        assertThat(s.extend(S2TextParser.makePoint(dst)))
            .withFailMessage("\nsrc = $src\ndst = $dst\ntarget = $target\navoid = $avoid")
            .isEqualTo(expectedResult)

    }

    @Test
    fun reuse() {
        // Check that Init() can be called more than once.
        val s = S2PolylineSimplifier()
        val radius = S1ChordAngle(S1Angle.degrees(10))
        s.init(S2Point(1, 0, 0))
        assertThat(s.targetDisc(S2Point(1, 1, 0).normalize(), radius)).isTrue()
        assertThat(s.targetDisc(S2Point(1.0, 1.0, 0.1).normalize(), radius)).isTrue()
        assertThat(s.extend(S2Point(1.0, 1.0, 0.4).normalize())).isFalse()

        s.init(S2Point(0, 1, 0));
        assertThat(s.targetDisc(S2Point(1.0, 1.0, 0.3).normalize(), radius)).isTrue()
        assertThat(s.targetDisc(S2Point(1.0, 1.0, 0.2).normalize(), radius)).isTrue()
        assertThat(s.extend(S2Point(1.0, 1.0, 0.0).normalize())).isFalse()
    }

    @Test
    fun noConstraints() {
        // No constraints, dst == src.
        checkSimplify("0:1", "0:1", "", "", listOf(), 0.0, true);

        // No constraints, dst != src.
        checkSimplify("0:1", "1:0", "", "", listOf(), 0.0, true);

        // No constraints, (src, dst) longer than 90 degrees (not supported).
        checkSimplify("0:0", "0:91", "", "", listOf(), 0.0, false);
    }

    @Test
    fun targetOnePoint() {
        // Three points on a straight line.  In theory zero tolerance should work,
        // but in practice there are floating point errors.
        checkSimplify("0:0", "0:2", "0:1", "", listOf(), 1e-10, true);

        // Three points where the middle point is too far away.
        checkSimplify("0:0", "0:2", "1:1", "", listOf(), 0.9, false);

        // A target disc that contains the source vertex.
        checkSimplify("0:0", "0:2", "0:0.1", "", listOf(), 1.0, true);

        // A target disc that contains the destination vertex.
        checkSimplify("0:0", "0:2", "0:2.1", "", listOf(), 1.0, true);
    }

    @Test
    fun avoidOnePoint() {
        // Three points on a straight line, attempting to avoid the middle point.
        checkSimplify("0:0", "0:2", "", "0:1", listOf(true), 1e-10, false);

        // Three points where the middle point can be successfully avoided.
        checkSimplify("0:0", "0:2", "", "1:1", listOf(true), 0.9, true);

        // Three points where the middle point is on the left, but where the client
        // requires the point to be on the right of the edge.
        checkSimplify("0:0", "0:2", "", "1:1", listOf(false), 1e-10, false);
    }

    @Test
    fun targetAndAvoid() {
        // Target several points that are separated from the proposed edge by about
        // 0.7 degrees, and avoid several points that are separated from the
        // proposed edge by about 1.4 degrees.
        checkSimplify(
            "0:0", "10:10", "2:3, 4:3, 7:8",
            "4:2, 7:5, 7:9", listOf(true, true, false), 1.0, true
        );

        // The same example, but one point to be targeted is 1.4 degrees away.
        checkSimplify(
            "0:0", "10:10", "2:3, 4:6, 7:8",
            "4:2, 7:5, 7:9", listOf(true, true, false), 1.0, false
        );

        // The same example, but one point to be avoided is 0.7 degrees away.
        checkSimplify(
            "0:0", "10:10", "2:3, 4:3, 7:8",
            "4:2, 6:5, 7:9", listOf(true, true, false), 1.0, false
        );
    }

    @Test
    fun precision() {
        // This is a rough upper bound on both the error in constructing the disc
        // locations (i.e., S2::InterpolateAtDistance, etc), and also on the
        // padding that S2PolylineSimplifier uses to ensure that its results are
        // conservative (i.e., the error calculated by GetSemiwidth).
        val kMaxError = S1Angle.radians(25 * DoubleType.epsilon)

        // We repeatedly generate a random edge.  We then target several discs that
        // barely overlap the edge, and avoid several discs that barely miss the
        // edge.  About half the time, we choose one disc and make it slightly too
        // large or too small so that targeting fails.
        val kIters = 1000;  // Passes with 1 million iterations.
        val simplifier = S2PolylineSimplifier()
        repeat(kIters) { iter ->
            S2Random.reset(iter + 1);  // Easier to reproduce a specific case.
            val src = S2Random.randomPoint()
            simplifier.init(src)
            val dst = S2EdgeDistances.interpolateAtDistance(
                S1Angle.radians(S2Random.randomDouble()),
                src,
                S2Random.randomPoint()
            )
            val n = S2PointUtil.robustCrossProd(src, dst).normalize()

            // If bad_disc >= 0, then we make targeting fail for that disc.
            val kNumDiscs = 5
            val bad_disc = S2Random.randomInt(2 * kNumDiscs) - kNumDiscs
            repeat(kNumDiscs) { i ->
                val f = S2Random.randomDouble()
                val a = (src * (1 - f) + dst * f).normalize()
                val r = S1Angle.radians(S2Random.randomDouble())
                val on_left = S2Random.oneIn(2)
                val x = S2EdgeDistances.interpolateAtDistance(r, a, if (on_left) n else -n)
                // We grow the radius slightly if we want to target the disc and shrink
                // it otherwise, *unless* we want targeting to fail for this disc, in
                // which case these actions are reversed.
                val avoid = S2Random.oneIn(2)
                val grow_radius = (avoid == (i == bad_disc))
                val radius = S1ChordAngle(if (grow_radius) (r + kMaxError) else (r - kMaxError))
                if (avoid) {
                    simplifier.avoidDisc(x, radius, on_left)
                } else {
                    simplifier.targetDisc(x, radius);
                }
            }
            // The result is true iff all the disc constraints were satisfiable.
            assertThat(simplifier.extend(dst)).isEqualTo(bad_disc < 0)
        }
    }

}
