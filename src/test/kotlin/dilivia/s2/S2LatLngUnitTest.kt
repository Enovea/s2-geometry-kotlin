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

import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.s2.S2LatLng.Companion.times
import dilivia.s2.S2Random.randomPoint
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.abs

class S2LatLngUnitTest {

    @Test
    fun testBasic() {
        val llRad = S2LatLng.fromRadians(M_PI_4, M_PI_2)
        assertThat(llRad.lat().radians).isEqualTo(M_PI_4)
        assertThat(llRad.lng().radians).isEqualTo(M_PI_2)
        assertThat(llRad.isValid).isTrue()
        val llDeg = S2LatLng.fromDegrees(45, 90)
        assertThat(llDeg).isEqualTo(llRad)
        assertThat(llDeg.isValid).isTrue()
        assertThat(S2LatLng.fromDegrees(-91, 0).isValid).isFalse()
        assertThat(S2LatLng.fromDegrees(0, 181).isValid).isFalse()

        var bad = S2LatLng.fromDegrees(120, 200)
        assertThat(bad.isValid).isFalse()
        var better = bad.normalized()
        assertThat(better.isValid).isTrue()
        assertThat(better.lat()).isEqualTo(S1Angle.degrees(90))
        assertThat(better.lng().radians).isEqualTo(S1Angle.degrees(-160).radians)

        bad = S2LatLng.fromDegrees(-100, -360)
        assertThat(bad.isValid).isFalse()
        better = bad.normalized()
        assertThat(better.isValid).isTrue()
        assertThat(better.lat()).isEqualTo(S1Angle.degrees(-90))
        assertThat(better.lng().radians).isCloseTo(0.0, Offset.offset(1e-15))

        assertThat(
            (S2LatLng.fromDegrees(10, 20) + S2LatLng.fromDegrees(20, 30)).approxEquals(
                S2LatLng.fromDegrees(
                    30,
                    50
                )
            )
        ).isTrue()
        assertThat(
            (S2LatLng.fromDegrees(10, 20) - S2LatLng.fromDegrees(20, 30)).approxEquals(
                S2LatLng.fromDegrees(
                    -10,
                    -10
                )
            )
        ).isTrue()
        assertThat((0.5 * S2LatLng.fromDegrees(10, 20)).approxEquals(S2LatLng.fromDegrees(5, 10))).isTrue()

        // Check that Invalid() returns an invalid point.
        assertThat(S2LatLng.invalid.isValid).isFalse()

        // Check that the default constructor sets latitude and longitude to 0.
        val center = S2LatLng.center
        assertThat(center.isValid).isTrue()
        assertThat(center.lat().radians).isEqualTo(0.0)
        assertThat(center.lng().radians).isEqualTo(0.0)
    }

    @Test
    fun conversion() {
        // Test special cases: poles, "date line"
        assertThat(S2LatLng.fromPoint(S2LatLng.fromDegrees(90.0, 65.0).toPoint()).lat().degrees()).isEqualTo(90.0)
        assertThat(S2LatLng.fromPoint(S2LatLng.fromRadians(-M_PI_2, 1.0).toPoint()).lat().radians).isEqualTo(-M_PI_2)
        assertThat(
            abs(
                S2LatLng.fromPoint(S2LatLng.fromDegrees(12.2, 180.0).toPoint()).lng().degrees()
            )
        ).isEqualTo(180.0)
        assertThat(abs(S2LatLng.fromPoint(S2LatLng.fromRadians(0.1, -M_PI).toPoint()).lng().radians)).isEqualTo(M_PI)

        // Test a bunch of random points.
        for (i in 0 until 100000) {
            val p = randomPoint()
            assertThat(p.approxEquals(S2LatLng.fromPoint(p).toPoint())).isTrue()
        }
    }

    @Test
    fun distance() {
        assertThat(S2LatLng.fromDegrees(90, 0).getDistance(S2LatLng.fromDegrees(90, 0)).radians).isEqualTo(
            0.0
        )
        assertThat(S2LatLng.fromDegrees(-37, 25).getDistance(S2LatLng.fromDegrees(-66, -155)).degrees()).isCloseTo(
            77.0,
            Offset.offset(1e-13)
        )
        assertThat(S2LatLng.fromDegrees(0, 165).getDistance(S2LatLng.fromDegrees(0, -80)).degrees()).isCloseTo(
            115.0,
            Offset.offset(1e-13)
        )
        assertThat(S2LatLng.fromDegrees(47, -127).getDistance(S2LatLng.fromDegrees(-47, 53)).degrees()).isCloseTo(
            180.0,
            Offset.offset(2e-6)
        )
    }

}
