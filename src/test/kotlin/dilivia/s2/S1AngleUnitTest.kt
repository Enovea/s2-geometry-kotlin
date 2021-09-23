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
import dilivia.s2.S1Angle.Companion.cos
import dilivia.s2.S1Angle.Companion.degrees
import dilivia.s2.S1Angle.Companion.e5
import dilivia.s2.S1Angle.Companion.e6
import dilivia.s2.S1Angle.Companion.e7
import dilivia.s2.S1Angle.Companion.radians
import dilivia.s2.S1Angle.Companion.sin
import dilivia.s2.S1Angle.Companion.tan
import dilivia.s2.S1Angle.Companion.times
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test

class S1AngleUnitTest {

    val epsilon = Offset.offset(DoubleType.epsilon)

    @Test
    fun basic() {
        // Check that the conversion between Pi radians and 180 degrees is exact.
        assertThat(radians(Math.PI).radians).isEqualTo(M_PI)
        assertThat(radians(Math.PI).degrees()).isEqualTo(180.0)
        assertThat(degrees(180.0).radians).isEqualTo(Math.PI)
        assertThat(degrees(180.0).degrees()).isEqualTo(180.0)
        assertThat(radians(Math.PI / 2).degrees()).isEqualTo(90.0)

        // Check negative angles.
        assertThat(radians(-Math.PI / 2).degrees()).isEqualTo(-90.0)
        assertThat(degrees(-45.0).radians).isEqualTo(-Math.PI / 4)

        // Check that E5/E6/E7 representations work as expected.
        assertThat(e5(2000000)).isEqualTo(degrees(20.0))
        assertThat(e6(-60000000)).isEqualTo(degrees(-60.0))
        assertThat(e7(750000000)).isEqualTo(degrees(75.0))
        assertThat(degrees(12.34567).e5()).isEqualTo(1234567)
        assertThat(degrees(12.345678).e6()).isEqualTo(12345678)
        assertThat(degrees(-12.3456789).e7()).isEqualTo(-123456789)
    }

    @Test
    fun defaultConstructor() {
        // Check that the default constructor returns an angle of 0.
        assertThat(S1Angle().radians).isZero()
    }

    @Test
    fun infinity() {
        assertThat(radians(1e30) < S1Angle.infinity()).isTrue()
        assertThat(-S1Angle.infinity() < S1Angle.zero()).isTrue()
        assertThat(S1Angle.infinity()).isEqualTo(S1Angle.infinity())
    }

    @Test
    fun zero() {
        assertThat(radians(0)).isEqualTo(S1Angle.zero())
        assertThat(radians(0.0)).isEqualTo(S1Angle.zero())
    }

    @Test
    fun piRadiansExactly180Degrees() {
        // Check that the conversion between Pi radians and 180 degrees is exact.
        assertThat(radians(M_PI).radians).isEqualTo(M_PI)
        assertThat(radians(M_PI).degrees()).isEqualTo(180.0)
        assertThat(degrees(180).radians).isEqualTo(M_PI)
        assertThat(degrees(180).degrees()).isEqualTo(180.0)

        assertThat(radians(M_PI_2).degrees()).isEqualTo(90.0)

        // Check negative angles.
        assertThat(radians(-M_PI_2).degrees()).isEqualTo(-90.0)
        assertThat(degrees(-45).radians).isEqualTo(-M_PI_4)
    }

    @Test
    fun e5e6e7Representations() {
        // Check that E5/E6/E7 representations work as expected.
        assertThat(e5(-4500000).radians).isCloseTo(degrees(-45).radians, epsilon)
        assertThat(e6(-60000000).radians).isCloseTo(degrees(-60).radians, epsilon)
        assertThat(e7(750000000).radians).isCloseTo(degrees(75).radians, epsilon)
        assertThat(degrees(-172.56123).e5()).isEqualTo(-17256123)
        assertThat(degrees(12.345678).e6()).isEqualTo(12345678)
        assertThat(degrees(-12.3456789).e7()).isEqualTo(-123456789)
    }

    @Test
    fun normalizeCorrectlyCanonicalizesAngles() {
        assertThat(degrees(360.0).normalized().degrees()).isEqualTo(0.0)
        assertThat(degrees(-90.0).normalized().degrees()).isEqualTo(-90.0)
        assertThat(degrees(-180.0).normalized().degrees()).isEqualTo(180.0)
        assertThat(degrees(180.0).normalized().degrees()).isEqualTo(180.0)
        assertThat(degrees(540.0).normalized().degrees()).isEqualTo(180.0)
        assertThat(degrees(-270.0).normalized().degrees()).isEqualTo(90.0)

        assertThat(degrees(360.0).normalize().degrees()).isEqualTo(0.0)
        assertThat(degrees(-90.0).normalize().degrees()).isEqualTo(-90.0)
        assertThat(degrees(-180.0).normalize().degrees()).isEqualTo(180.0)
        assertThat(degrees(180.0).normalize().degrees()).isEqualTo(180.0)
        assertThat(degrees(540.0).normalize().degrees()).isEqualTo(180.0)
        assertThat(degrees(-270.0).normalize().degrees()).isEqualTo(90.0)
    }

    @Test
    fun arithmeticOperationsOnAngles() {
        assertThat(radians(-0.3).abs().radians).isCloseTo(0.3, epsilon)
        assertThat((-radians(0.1)).radians).isCloseTo(-0.1, epsilon)
        assertThat((radians(0.1) + radians(0.3)).radians).isCloseTo(0.4, epsilon)
        assertThat((radians(0.1) - radians(0.3)).radians).isCloseTo(-0.2, epsilon)
        assertThat((radians(2) * radians(0.3)).radians).isCloseTo(0.6, epsilon)
        assertThat((2.0 * radians(0.3)).radians).isCloseTo(0.6, epsilon)
        assertThat((radians(0.3) * 2.0).radians).isCloseTo(0.6, epsilon)
        assertThat((radians(0.3) / 2.0).radians).isCloseTo(0.15, epsilon)
        assertThat((radians(0.3) / radians(0.6)).radians).isCloseTo(0.5, epsilon)

        val tmp = radians(1.0)
        tmp += radians(0.5)
        assertThat(tmp.radians).isEqualTo(1.5)
        tmp -= radians(1.0)
        assertThat(tmp.radians).isEqualTo(0.5)
        tmp *= 5.0
        assertThat(tmp.radians).isEqualTo(2.5)
        tmp /= 2.0
        assertThat(tmp.radians).isEqualTo(1.25)
    }

    @Test
    fun trigonometry() {
        // Spot check a few angles to ensure that the correct function is called.
        assertThat(cos(degrees(0))).isCloseTo(1.0, epsilon)
        assertThat(sin(degrees(90))).isCloseTo(1.0, epsilon)
        assertThat(tan(degrees(45))).isCloseTo(1.0, epsilon)
    }

    @Test
    fun constructorsThatMeasureAngles() {
        assertThat(S1Angle(S2Point(1, 0, 0), S2Point(0, 0, 2)).radians).isEqualTo(M_PI_2)
        assertThat(S1Angle(S2Point(1, 0, 0), S2Point(1, 0, 0)).radians).isEqualTo(0.0)
        assertThat(S1Angle(S2LatLng.fromDegrees(20, 20), S2LatLng.fromDegrees(70, 20)).degrees()).isCloseTo(
            50.0,
            Offset.offset(1.5e-14)
        )
    }

    @Test
    fun textFormatting() {
        assertThat(degrees(180.0).toString()).isEqualTo("180.0000000d")
    }

    @Test
    fun degreesVsE6() {
        //  The current implementation guarantees exact conversions between
        //  Degrees() and E6() when the Degrees() argument is an integer.
        for (i in 0..180) {
            assertThat(e6(1000000 * i)).isEqualTo(degrees(i))
        }
    }

    @Test
    fun degreesVsE7() {
        // The current implementation guarantees exact conversions between
        // Degrees() and E7() when the Degrees() argument is an integer.
        for (i in 0..180) {
            assertThat(e7(10000000 * i)).isEqualTo(degrees(i))
        }
    }

    @Test
    fun e6VsE7() {
        // The current implementation guarantees exact conversions between
        // E6() and E7() when the E6() argument is an integer.
        for (iter in 0..1000) {
            val i = S2Random.randomInt(0, 180000000)
            assertThat(e6(i)).isEqualTo(e7(10 * i))
        }
    }

    @Test
    fun degreesVsRadians() {
        // The current implementation guarantees certain exact conversions between
        // degrees and radians (see the header file for details).
        for (k in -8..8) {
            assertThat(degrees(45 * k)).isEqualTo(radians(k * M_PI / 4))
            assertThat(45.0 * k).isEqualTo(degrees(45 * k).degrees())
        }
        for (k in 0..30) {
            val n = 1 shl k
            assertThat(degrees(180.0 / n)).isEqualTo(radians(M_PI / n))
            assertThat(degrees(60.0 / n)).isEqualTo(radians(M_PI / (3.0 * n)))
            assertThat(degrees(36.0 / n)).isEqualTo(radians(M_PI / (5.0 * n)))
            assertThat(degrees(20.0 / n)).isEqualTo(radians(M_PI / (9.0 * n)))
            assertThat(degrees(4.0 / n)).isEqualTo(radians(M_PI / (45.0 * n)))
        }
        // We also spot check a couple of non-identities.
        assertThat(degrees(3) != radians(M_PI / 60)).isTrue()
        assertThat(60.0 != degrees(60).degrees()).isTrue()
    }

}
