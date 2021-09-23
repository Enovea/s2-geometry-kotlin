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
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.s2.S1ChordAngle.Companion.cos
import dilivia.s2.S1ChordAngle.Companion.sin
import dilivia.s2.S1ChordAngle.Companion.tan
import dilivia.s2.edge.S2EdgeDistances
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.atan
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.tan

class S1ChordAngleUnitTest {

    @Test
    fun defaultConstructor() {
        // Check that the default constructor returns an angle of 0.
        val a = S1ChordAngle()
        assertThat(a).isEqualTo(S1ChordAngle.zero())
        assertThat(a.length2).isEqualTo(0.0)
    }

    @Test
    fun cloneAngle() {
        val a = S1ChordAngle.fromLength2(0.23)
        val b = a.clone()
        b.length2 = 0.43
        assertThat(a.length2).isEqualTo(0.23)
        assertThat(b).isNotEqualTo(a)
    }

    @Test
    fun twoPointConstructor() {
        for (iter in 0 until 100) {
            val (x, y, z) = S2Random.randomFrame()
            assertThat(S1ChordAngle.between(z, z).toAngle()).isEqualTo(S1Angle.zero())
            assertThat(S1ChordAngle.between(z.unaryMinus(), z).radians())
                .withFailMessage("[x=$x;y=$y;z=$z] angle between z and -z is not close to M_PI.")
                .isCloseTo(
                    M_PI,
                    Offset.offset(1e-7)
                )
            assertThat(S1ChordAngle.between(x, z).radians()).isCloseTo(M_PI_2, Offset.offset(1e-15))
            val w = (y + z).normalize()
            assertThat(S1ChordAngle.between(w, z).radians()).isCloseTo(M_PI_4, Offset.offset(1e-15))
        }
    }

    @Test
    fun fromLength2() {
        assertThat(S1ChordAngle.fromLength2(0).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat(S1ChordAngle.fromLength2(1).degrees()).isCloseTo(60.0, Offset.offset(8e-15))
        assertThat(S1ChordAngle.fromLength2(2).degrees()).isCloseTo(90.0, Offset.offset(2e-14))
        assertThat(S1ChordAngle.fromLength2(4).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat(S1ChordAngle.fromLength2(5).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
    }

    @Test
    fun zero() {
        assertThat(S1ChordAngle.zero().toAngle()).isEqualTo(S1Angle.zero())
    }

    @Test
    fun right() {
        assertThat(S1ChordAngle.right().degrees()).isCloseTo(90.0, Offset.offset(2e-14))
    }

    @Test
    fun straight() {
        assertThat(S1ChordAngle.straight().toAngle()).isEqualTo(S1Angle.degrees(180.0))
    }

    @Test
    fun infinity() {
        assertThat(S1ChordAngle.straight() < S1ChordAngle.infinity()).isTrue()
        assertThat(S1ChordAngle.infinity().length2).isEqualTo(Double.MAX_VALUE)
        assertThat(S1ChordAngle.infinity().toAngle()).isEqualTo(S1Angle.infinity())
    }

    @Test
    fun negative() {
        assertThat(S1ChordAngle.negative() < S1ChordAngle.zero()).isTrue()
        assertThat(S1ChordAngle.negative().length2).isEqualTo(-1.0)
        assertThat(S1ChordAngle.negative().toAngle() < S1Angle.zero()).isTrue()
    }

    @Test
    fun e5e6e7Representations() {
        // Check that E5/E6/E7 representations work as expected.
        assertThat(S1ChordAngle.e5(-4500000).length2).isCloseTo(S1ChordAngle.degrees(-45).length2, Offset.offset(1e-15))
        assertThat(S1ChordAngle.e6(-60000000).length2).isCloseTo(S1ChordAngle.degrees(-60).length2, Offset.offset(1e-15))
        assertThat(S1ChordAngle.e7(750000000).length2).isCloseTo(S1ChordAngle.degrees(75).length2, Offset.offset(1e-15))
        assertThat(S1ChordAngle.degrees(12.3456789).e5()).isEqualTo(1234568)
        assertThat(S1ChordAngle.degrees(12.3456789).e6()).isEqualTo(12345679)
        assertThat(S1ChordAngle.degrees(12.3456789).e7()).isEqualTo(123456789)
    }

    @Test
    fun predicates() {
        assertThat(S1ChordAngle.zero().isZero()).isTrue()
        assertThat(S1ChordAngle.zero().isNegative()).isFalse()
        assertThat(S1ChordAngle.zero().isSpecial()).isFalse()
        assertThat(S1ChordAngle.straight().isSpecial()).isFalse()
        assertThat(S1ChordAngle.negative().isNegative()).isTrue()
        assertThat(S1ChordAngle.negative().isSpecial()).isTrue()
        assertThat(S1ChordAngle.infinity().isInfinity()).isTrue()
        assertThat(S1ChordAngle.infinity().isSpecial()).isTrue()
    }

    @Test
    fun toFromS1Angle() {
        assertThat(S1ChordAngle(S1Angle.zero()).radians()).isEqualTo(0.0)
        assertThat(S1ChordAngle(S1Angle.radians(M_PI)).length2).isEqualTo(4.0)
        assertThat(S1ChordAngle(S1Angle.radians(M_PI)).radians()).isEqualTo(M_PI)
        assertThat(S1ChordAngle(S1Angle.infinity()).toAngle()).isEqualTo(S1Angle.infinity())
        assertThat(S1ChordAngle(S1Angle.radians(-1)).radians() < 0).isTrue()
        assertThat(S1ChordAngle(S1Angle.radians(1.0)).radians()).isEqualTo(1.0)
    }

    @Test
    fun successor() {
        assertThat(S1ChordAngle.negative().successor()).isEqualTo(S1ChordAngle.zero())
        assertThat(S1ChordAngle.straight().successor()).isEqualTo(S1ChordAngle.infinity())
        assertThat(S1ChordAngle.infinity().successor()).isEqualTo(S1ChordAngle.infinity())
        var x = S1ChordAngle.negative()
        for (i in 0 until 10) {
            assertThat(x < x.successor()).isTrue()
            x = x.successor()
        }
    }

    @Test
    fun predecessor() {
        assertThat(S1ChordAngle.infinity().predecessor()).isEqualTo(S1ChordAngle.straight())
        assertThat(S1ChordAngle.zero().predecessor()).isEqualTo(S1ChordAngle.negative())
        assertThat(S1ChordAngle.negative().predecessor()).isEqualTo(S1ChordAngle.negative())
        var x = S1ChordAngle.infinity()
        for (i in 0 until 10) {
            assertThat(x > x.predecessor()).isTrue()
            x = x.predecessor()
        }
    }

    @Test
    fun arithmetic() {
        val zero = S1ChordAngle.zero()
        val degree30 = S1ChordAngle.degrees(30)
        val degree60 = S1ChordAngle.degrees(60)
        val degree90 = S1ChordAngle.degrees(90)
        val degree120 = S1ChordAngle.degrees(120)
        val degree180 = S1ChordAngle.straight()
        assertThat((zero + zero).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((zero - zero).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree60 - degree60).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree180 - degree180).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((zero - degree60).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree30 - degree90).degrees()).isCloseTo(0.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree60 + zero).degrees()).isCloseTo(60.0, Offset.offset(8e-15))
        assertThat((degree60 - zero).degrees()).isCloseTo(60.0, Offset.offset(8e-15))
        assertThat((zero + degree60).degrees()).isCloseTo(60.0, Offset.offset(8e-15))
        assertThat((degree30 + degree60).degrees()).isCloseTo(90.0, Offset.offset(2e-14))
        assertThat((degree60 + degree30).degrees()).isCloseTo(90.0, Offset.offset(2e-14))
        assertThat((degree90 - degree30).degrees()).isCloseTo(60.0, Offset.offset(8e-15))
        assertThat((degree90 - degree60).degrees()).isCloseTo(30.0, Offset.offset(8e-15))
        assertThat((degree180 + zero).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree180 - zero).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree90 + degree90).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree120 + degree90).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree120 + degree120).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree30 + degree180).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
        assertThat((degree180 + degree180).degrees()).isCloseTo(180.0, Offset.offset(DoubleType.epsilon))
    }

    @Test
    fun trigonometry() {
        val kIters = 20
        for (iter in 0..kIters) {
            val radians = M_PI * iter / kIters
            val angle = S1ChordAngle(S1Angle.radians(radians))
            assertThat(sin(angle)).isCloseTo(sin(radians), Offset.offset(1e-15))
            assertThat(cos(angle)).isCloseTo(cos(radians), Offset.offset(1e-15))
            // Since the tan(x) is unbounded near Pi/4, we map the result back to an
            // angle before comparing.  (The assertion is that the result is equal to
            // the tangent of a nearby angle.)
            assertThat(atan(tan(angle))).isCloseTo(atan(tan(radians)), Offset.offset(1e-15))
        }

        // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
        val angle90 = S1ChordAngle.fromLength2(2)
        val angle180 = S1ChordAngle.fromLength2(4)
        assertThat(sin(angle90)).isEqualTo(1.0)
        assertThat(cos(angle90)).isEqualTo(0.0)
        assertThat(tan(angle90)).isEqualTo(Double.POSITIVE_INFINITY)
        assertThat(sin(angle180)).isEqualTo(0.0)
        assertThat(cos(angle180)).isEqualTo(-1.0)
        assertThat(tan(angle180)).isEqualTo(0.0)
    }

    @Test
    fun plusError() {
        assertThat(S1ChordAngle.negative().plusError(5)).isEqualTo(S1ChordAngle.negative())
        assertThat(S1ChordAngle.infinity().plusError(-5)).isEqualTo(S1ChordAngle.infinity())
        assertThat(S1ChordAngle.straight().plusError(5)).isEqualTo(S1ChordAngle.straight())
        assertThat(S1ChordAngle.zero().plusError(-5)).isEqualTo(S1ChordAngle.zero())
        assertThat(S1ChordAngle.fromLength2(1).plusError(0.25)).isEqualTo(S1ChordAngle.fromLength2(1.25))
        assertThat(S1ChordAngle.fromLength2(1).plusError(-0.25)).isEqualTo(S1ChordAngle.fromLength2(0.75))
    }

    @Test
    fun getS2PointConstructorMaxError() {
        // Check that the error bound returned by GetS2PointConstructorMaxError() is
        // large enough.
        for (iter in 0 until 100000) {
            S2Random.reset(iter)     // Easier to reproduce a specific case.
            val x = S2Random.randomPoint()
            var y = S2Random.randomPoint()
            if (S2Random.randomInt(0, 10) == 0) {
                // Occasionally test a point pair that is nearly identical or antipodal.
                val r = S1Angle.radians(1e-15 * S2Random.randomDouble())
                y = S2EdgeDistances.interpolateAtDistance(r, x, y)
                if (S2Random.oneIn(2)) y = -y
            }
            val dist = S1ChordAngle(x, y)
            val error = dist.getS2PointConstructorMaxError()
            assertThat(S2Predicates.compareDistance(x, y, dist.plusError(error)) <= 0.0).isTrue()
            assertThat(S2Predicates.compareDistance(x, y, dist.plusError(-error)) >= 0.0).isTrue()
        }
    }


}
