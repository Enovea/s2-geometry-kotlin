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
import org.assertj.core.api.Assertions
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.EnumSource
import tech.units.indriya.unit.Units
import kotlin.math.pow

class S2EarthUnitTest {

    @Test
    fun angleConversion() {
        assertThat(S2Earth.toAngle(S2Earth.radius).radians).isCloseTo(1.0, Offset.offset(1e-15))
        assertThat(S2Earth.toChordAngle(S2Earth.radius).radians()).isCloseTo(1.0, Offset.offset(1e-15))
        assertThat(S2Earth.toDistance(S1Angle.radians(2)).to(Units.METRE).value.toDouble()).isCloseTo(
            2 * S2Earth.radiusMeters,
            Offset.offset(1e-7)
        )
        assertThat(
            S2Earth.toDistance(S1ChordAngle.radians(2.0)).to(Units.METRE).value.toDouble()
        ).isCloseTo(2 * S2Earth.radiusMeters, Offset.offset(1e-7))
        assertThat(S2Earth.toRadians(S2Earth.radius)).isCloseTo(1.0, Offset.offset(1e-15))
        assertThat(S2Earth.toMeters(S1Angle.degrees(180))).isCloseTo(S2Earth.radiusMeters * M_PI, Offset.offset(1e-15))
        assertThat(S2Earth.toMeters(S1ChordAngle.degrees(180))).isCloseTo(
            S2Earth.radiusMeters * M_PI,
            Offset.offset(1e-15)
        )
        assertThat(S2Earth.toKm(S1Angle.radians(0.5))).isCloseTo(0.5 * S2Earth.radiusKm, Offset.offset(1e-15))
        assertThat(S2Earth.toKm(S1ChordAngle.radians(0.5))).isCloseTo(0.5 * S2Earth.radiusKm, Offset.offset(1e-15))
        assertThat(S2Earth.kmToRadians(S2Earth.radiusMeters / 1000)).isCloseTo(1.0, Offset.offset(1e-15))
        assertThat(S2Earth.radiansToKm(0.5)).isCloseTo(0.5 * S2Earth.radiusKm, Offset.offset(1e-15))
        assertThat(S2Earth.metersToRadians(S2Earth.radiansToKm(0.3) * 1000)).isCloseTo(0.3, Offset.offset(1e-15))
        assertThat(S2Earth.radiansToMeters(S2Earth.kmToRadians(2.5))).isCloseTo(2500.0, Offset.offset(1e-15))
    }

    @Test
    fun solidAngleConversion() {
        assertThat(S2Earth.squareKmToSteradians((S2Earth.radiusMeters / 1000).pow(2.0))).isCloseTo(
            1.0,
            Offset.offset(1e-15)
        )
        assertThat(S2Earth.steradiansToSquareKm(0.5.pow(2.0))).isCloseTo(
            (0.5 * S2Earth.radiusKm).pow(2.0),
            Offset.offset(1e-15)
        )
        assertThat(S2Earth.squareMetersToSteradians((S2Earth.radiansToKm(0.3) * 1000).pow(2.0))).isCloseTo(
            0.3.pow(2.0),
            Offset.offset(1e-15)
        )
        assertThat(S2Earth.steradiansToSquareMeters(S2Earth.kmToRadians(2.5).pow(2.0))).isCloseTo(
            2500.0.pow(2.0),
            Offset.offset(1e-9)
        )
    }

    @Test
    fun toLongitudeRadians() {
        // At the equator, ToLongitudeRadians behaves exactly like ToRadians.
        assertThat(S2Earth.toLongitudeRadians(S2Earth.radius, 0.0)).isCloseTo(1.0, Offset.offset(1e-15))

        // The closer we get to the poles, the more radians we need to go the same
        // distance.
        assertThat(S2Earth.toLongitudeRadians(S2Earth.radius, 0.5)).isGreaterThan(S2Earth.toLongitudeRadians(S2Earth.radius, 0.4))

        // At the poles, we should return 2PI radians instead of dividing by 0.
        assertThat(S2Earth.toLongitudeRadians(S2Earth.radius, M_PI_2)).isCloseTo(M_PI * 2, Offset.offset(1e-15))

        // Within epsilon of the poles, we should still return 2PI radians instead
        // of directing the caller to take thousands of radians around.
        assertThat(S2Earth.toLongitudeRadians(S2Earth.radius, M_PI_2 - 1e-4)).isCloseTo(M_PI * 2, Offset.offset(1e-15))
    }

    enum class TestConfig(
        val description: String,
        val a: S2LatLng,
        val b: S2LatLng,
        val bearing: S1Angle,
    ) {
        WESTWARD_IN_EQUATOR(
            "Westward on equator",
            S2LatLng.fromDegrees(0, 50),
            S2LatLng.fromDegrees(0, 100),
            S1Angle.degrees(90)
        ),
        EASTWARD_ON_EQUATOR(
            "Eastward on equator",
            S2LatLng.fromDegrees(0, 50),
            S2LatLng.fromDegrees(0, 0),
            S1Angle.degrees(-90)
        ),
        NORTHWARD_ON_MERIDIAN(
            "Northward on meridian",
            S2LatLng.fromDegrees(16, 28),
            S2LatLng.fromDegrees(81, 28),
            S1Angle.degrees(0)
        ),
        SOUTHWARD_ON_MERIDIAN(
            "Southward on meridian",
            S2LatLng.fromDegrees(24, 64),
            S2LatLng.fromDegrees(-27, 64),
            S1Angle.degrees(180)
        ),
        TOWARDS_NORTH_POLE(
            "Towards north pole",
            S2LatLng.fromDegrees(12, 76),
            S2LatLng.fromDegrees(90, 50),
            S1Angle.degrees(0)
        ),
        TOWARDS_SOUTH_POLE(
            "Towards south pole",
            S2LatLng.fromDegrees(-35, 105),
            S2LatLng.fromDegrees(-90, -120),
            S1Angle.degrees(180)
        ),
        SPAIN_TO_JAPAN(
            "Spain to Japan",
            S2LatLng.fromDegrees(40.4379332, -3.749576),
            S2LatLng.fromDegrees(35.6733227, 139.6403486),
            S1Angle.degrees(29.2)
        ),
        JAPAN_TO_SPAIN(
            "Japan to Spain",
            S2LatLng.fromDegrees(35.6733227, 139.6403486),
            S2LatLng.fromDegrees(40.4379332, -3.749576),
            S1Angle.degrees(-27.2)
        );
    }

    @ParameterizedTest
    @EnumSource(value = TestConfig::class)
    fun getInitialBearing(config: TestConfig) {
        val bearing = S2Earth.getInitialBearing(config.a, config.b)
        val angleDiff = (bearing - config.bearing).normalize().abs()
        Assertions.assertThat(angleDiff.degrees())
            .withFailMessage("getInitialBearing() test failed on: ${config.description}. Expected ${config.bearing}, got $bearing")
            .isLessThanOrEqualTo(1e-2)
    }

    @Test
    fun getDistance() {
        val north = S2Point(0, 0, 1)
        val south = S2Point(0, 0, -1)
        val west = S2Point(0, -1, 0)

        assertThat(
            S2Earth.getDistance(north, south).to(Units.METRE).value.toDouble()
        ).isCloseTo(M_PI * S2Earth.radiusMeters, Offset.offset(1e-7))
        assertThat(S2Earth.getDistanceKm(west, west)).isCloseTo(0.0, Offset.offset(1e-15))
        assertThat(S2Earth.getDistanceMeters(north, west)).isCloseTo(
            M_PI_2 * S2Earth.radiusMeters,
            Offset.offset(1e-15)
        )

        assertThat(
            S2Earth.getDistance(S2LatLng.fromDegrees(0, -90), S2LatLng.fromDegrees(-90, -38))
                .to(Units.METRE).value.toDouble()
        ).isCloseTo(S2Earth.getDistance(west, south).to(Units.METRE).value.toDouble(), Offset.offset(1e-7))

        assertThat(S2Earth.getDistanceKm(S2LatLng.fromRadians(0.0, 0.6), S2LatLng.fromRadians(0.0, -0.4))).isCloseTo(
            S2Earth.radiusKm,
            Offset.offset(1e-15)
        )

        assertThat(S2Earth.getDistanceMeters(S2LatLng.fromDegrees(80, 27), S2LatLng.fromDegrees(55, -153))).isCloseTo(
            1000 * S2Earth.radiusKm * M_PI / 4,
            Offset.offset(1e-9)
        )
    }

}
