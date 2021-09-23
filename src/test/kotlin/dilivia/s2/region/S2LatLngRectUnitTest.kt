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

import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.R1Interval
import dilivia.math.vectors.R2VectorDouble
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.min
import dilivia.s2.S1Angle.Companion.zero
import dilivia.s2.S1Interval
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random.randomPoint
import dilivia.s2.edge.S2EdgeDistances
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.IEEErem
import kotlin.math.abs
import kotlin.math.sin
import kotlin.random.Random

@Suppress("UsePropertyAccessSyntax")
class S2LatLngRectUnitTest {

    @Test
    fun emptyAndFull() {
        // Test basic properties of empty and full rectangles.
        val empty = S2LatLngRect.empty()
        val full = S2LatLngRect.full()
        assertThat(empty.isValid).isTrue()
        assertThat(empty.isEmpty).isTrue()
        assertThat(empty.isPoint).isFalse()
        assertThat(full.isValid).isTrue()
        assertThat(full.isFull).isTrue()
        assertThat(full.isPoint).isFalse()
        // Check that the default S2LatLngRect is identical to Empty().
        val defaultEmpty = S2LatLngRect()
        assertThat(defaultEmpty.isValid).isTrue()
        assertThat(defaultEmpty.isEmpty).isTrue()
        assertThat(empty.lat.bounds).isEqualTo(defaultEmpty.lat.bounds)
        assertThat(empty.lng.bounds).isEqualTo(defaultEmpty.lng.bounds)
    }

    @Test
    fun accessors() {
        // Check various accessor methods.
        val d1 = rectFromDegrees(-90.0, 0.0, -45.0, 180.0)
        assertThat(d1.latLo().degrees()).isEqualTo(-90.0)
        assertThat(d1.latHi().degrees()).isEqualTo(-45.0)
        assertThat(d1.lngLo().degrees()).isEqualTo(0.0)
        assertThat(d1.lngHi().degrees()).isEqualTo(180.0)
        assertThat(d1.lat).isEqualTo(R1Interval(-M_PI_2, -M_PI_4))
        assertThat(d1.lng).isEqualTo(S1Interval(0.0, M_PI))
    }

    @Test
    fun approxEquals() {
        // S1Interval and R1Interval have additional testing.

        assertThat(S2LatLngRect.empty().approxEquals(rectFromDegrees(1.0, 5.0, 1.0, 5.0))).isTrue()
        assertThat(rectFromDegrees(1.0, 5.0, 1.0, 5.0).approxEquals(S2LatLngRect.empty())).isTrue()
        assertThat(rectFromDegrees(1.0, 5.0, 1.0, 5.0).approxEquals(rectFromDegrees(2.0, 7.0, 2.0, 7.0))).isFalse()

        // Test the max_error (double) parameter.
        assertThat(rectFromDegrees(10.0, 10.0, 20.0, 20.0).approxEquals(rectFromDegrees(11.0, 11.0, 19.0, 19.0), S1Angle.degrees(1.001))).isTrue()
        assertThat(rectFromDegrees(10.0, 10.0, 20.0, 20.0).approxEquals(rectFromDegrees(11.0, 11.0, 19.0, 19.0), S1Angle.degrees(0.999))).isFalse()

        // Test the max_error (S2LatLng) parameter.
        assertThat(rectFromDegrees(0.0, 10.0, 20.0, 30.0).approxEquals(rectFromDegrees(-1.0, 8.0, 21.0, 32.0), S2LatLng.fromDegrees(1.001, 2.001))).isTrue()
        assertThat(rectFromDegrees(0.0, 10.0, 20.0, 30.0).approxEquals(rectFromDegrees(-1.0, 8.0, 21.0, 32.0), S2LatLng.fromDegrees(0.999, 1.999))).isFalse()
    }

    @Test
    fun fromCenterSize() {
        assertThat(S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(80, 170), S2LatLng.fromDegrees(40, 60)).approxEquals(rectFromDegrees(60.0, 140.0, 90.0, -160.0))).isTrue()
        assertThat(S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(10, 40), S2LatLng.fromDegrees(210, 400)).isFull).isTrue()
        assertThat(S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(-90, 180), S2LatLng.fromDegrees(20, 50)).approxEquals(rectFromDegrees(-90.0, 155.0, -80.0, -155.0))).isTrue()
    }

    @Test
    fun fromPoint() {
        val p = S2LatLng.fromDegrees(23, 47)
        assertThat(S2LatLngRect.fromPoint(p)).isEqualTo(S2LatLngRect(p, p))
        assertThat(S2LatLngRect.fromPoint(p).isPoint).isTrue()
    }

    @Test
    fun fromPointPair() {
        assertThat(rectFromDegrees(-35.0, 155.0, 15.0, -140.0)).isEqualTo(S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(-35, -140), S2LatLng.fromDegrees(15, 155)))
        assertThat(rectFromDegrees(-90.0, -70.0, 25.0, 80.0)  ).isEqualTo(S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(25, -70), S2LatLng.fromDegrees(-90, 80))  )
    }

    @Test
    fun getCenterSize() {
        val r1 = S2LatLngRect(R1Interval(0.0, M_PI_2), S1Interval(-M_PI, 0.0))
        assertThat(r1.center).isEqualTo(S2LatLng.fromRadians(M_PI_4, -M_PI_2))
        assertThat(r1.size  ).isEqualTo(S2LatLng.fromRadians(M_PI_2, M_PI))
        assertThat(S2LatLngRect.empty().size.lat().radians < 0).isTrue()
        assertThat(S2LatLngRect.empty().size.lng().radians < 0).isTrue()
    }

    @Test
    fun getVertex() {
        val r1 = S2LatLngRect(R1Interval(0.0, M_PI_2), S1Interval(-M_PI, 0.0))
        assertThat(r1.getVertex(0)).isEqualTo(S2LatLng.fromRadians(0.0, M_PI))
        assertThat(r1.getVertex(1)).isEqualTo(S2LatLng.fromRadians(0.0, 0.0))
        assertThat(r1.getVertex(2)).isEqualTo(S2LatLng.fromRadians(M_PI_2, 0.0))
        assertThat(r1.getVertex(3)).isEqualTo(S2LatLng.fromRadians(M_PI_2, M_PI))

        // Make sure that GetVertex() returns vertices in CCW order.
        for (i in 0..3) {
            val lat = M_PI_4 * (i - 2)
            val lng = M_PI_2 * (i - 2) + 0.2
            val r = S2LatLngRect(R1Interval(lat, lat + M_PI_4), S1Interval(lng.IEEErem(2 * M_PI), (lng + M_PI_2).IEEErem(2 * M_PI)))
            for (k in 0..3) {
                assertThat(S2PointUtil.simpleCCW(r.getVertex(k - 1).toPoint(), r.getVertex(k).toPoint(), r.getVertex(k + 1).toPoint())).isTrue()
            }
        }
    }

    @Test
    fun contains() {
        // Contains(S2LatLng), InteriorContains(S2LatLng), Contains()
        val eqM180 = S2LatLng.fromRadians(0.0, -M_PI)
        val northPole = S2LatLng.fromRadians(M_PI_2, 0.0)
        val r1 = S2LatLngRect(eqM180, northPole)

        assertThat(r1.contains(S2LatLng.fromDegrees(30, -45))).isTrue()
        assertThat(r1.interiorContains(S2LatLng.fromDegrees(30, -45))).isTrue()
        assertThat(r1.contains(S2LatLng.fromDegrees(30, 45))).isFalse()
        assertThat(r1.interiorContains(S2LatLng.fromDegrees(30, 45))).isFalse()
        assertThat(r1.contains(eqM180)).isTrue()
        assertThat(r1.interiorContains(eqM180)).isFalse()
        assertThat(r1.contains(northPole)).isTrue()
        assertThat(r1.interiorContains(northPole)).isFalse()
        assertThat(r1.contains(S2Point(0.5, -0.3, 0.1))).isTrue()
        assertThat(r1.contains(S2Point(0.5, 0.2, 0.1))).isFalse()
    }

    @Test
    fun intervalOps() {
        // Contains(S2LatLngRect), InteriorContains(S2LatLngRect),
        // Intersects(), InteriorIntersects(), Union(), Intersection().
        //
        // Much more testing of these methods is done in s1interval_test
        // and r1interval_test.

        // Rectangle "r1" covers one-quarter of the sphere.
        val r1 = rectFromDegrees(0.0, -180.0, 90.0, 0.0)

        // Test operations where one rectangle consists of a single point.
        val r1Mid = rectFromDegrees(45.0, -90.0, 45.0, -90.0)
        testIntervalOps(r1, r1Mid, "TTTT", r1, r1Mid)

        val reqM180 = rectFromDegrees(0.0, -180.0, 0.0, -180.0)
        testIntervalOps(r1, reqM180, "TFTF", r1, reqM180)

        val rnorthPole = rectFromDegrees(90.0, 0.0, 90.0, 0.0)
        testIntervalOps(r1, rnorthPole, "TFTF", r1, rnorthPole)

        testIntervalOps(r1, rectFromDegrees(-10, -1, 1, 20), "FFTT",
                rectFromDegrees(-10, 180, 90, 20),
                rectFromDegrees(0, -1, 1, 0))
        testIntervalOps(r1, rectFromDegrees(-10, -1, 0, 20), "FFTF",
                rectFromDegrees(-10, 180, 90, 20),
                rectFromDegrees(0, -1, 0, 0))
        testIntervalOps(r1, rectFromDegrees(-10, 0, 1, 20), "FFTF",
                rectFromDegrees(-10, 180, 90, 20),
                rectFromDegrees(0, 0, 1, 0))

        testIntervalOps(rectFromDegrees(-15, -160, -15, -150),
                rectFromDegrees(20, 145, 25, 155), "FFFF",
                rectFromDegrees(-15, 145, 25, -150),
                S2LatLngRect.empty())
        testIntervalOps(rectFromDegrees(70, -10, 90, -140),
                rectFromDegrees(60, 175, 80, 5), "FFTT",
                rectFromDegrees(60, -180, 90, 180),
                rectFromDegrees(70, 175, 80, 5))

        // Check that the intersection of two rectangles that overlap in latitude
        // but not longitude is valid, and vice versa.
        testIntervalOps(rectFromDegrees(12, 30, 60, 60),
                rectFromDegrees(0, 0, 30, 18), "FFFF",
                rectFromDegrees(0, 0, 60, 60), S2LatLngRect.empty())
        testIntervalOps(rectFromDegrees(0, 0, 18, 42),
                rectFromDegrees(30, 12, 42, 60), "FFFF",
                rectFromDegrees(0, 0, 42, 60), S2LatLngRect.empty())
    }

    @Test
    fun boundaryIntersectsEmptyRectangle() {
        val rect = S2LatLngRect.empty()
        val lo = rect.lo().toPoint()
        val hi = rect.hi().toPoint()
        assertThat(rect.boundaryIntersects(lo, lo)).isFalse()
        assertThat(rect.boundaryIntersects(lo, hi)).isFalse()
    }

    @Test
    fun boundaryIntersectsFullRectangle() {
        val rect = S2LatLngRect.full()
        val lo = rect.lo().toPoint()
        val hi = rect.hi().toPoint()
        assertThat(rect.boundaryIntersects(lo, lo)).isFalse()
        assertThat(rect.boundaryIntersects(lo, hi)).isFalse()
    }

    @Test
    fun boundaryIntersectsSphericalLune() {
        // This rectangle only has two non-degenerate sides.
        val rect = rectFromDegrees(-90, 100, 90, 120)
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(60, 60).toPoint(), S2LatLng.fromDegrees(90, 60).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-60, 110).toPoint(), S2LatLng.fromDegrees(60, 110).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-60, 95).toPoint(), S2LatLng.fromDegrees(60, 110).toPoint())).isTrue()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(60, 115).toPoint(), S2LatLng.fromDegrees(80, 125).toPoint())).isTrue()
    }

    @Test
    fun boundaryIntersectsNorthHemisphere() {
        // This rectangle only has only one non-degenerate side.
        val rect = rectFromDegrees(0, -180, 90, 180)
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(60, -180).toPoint(), S2LatLng.fromDegrees(90, -180).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(60, -170).toPoint(), S2LatLng.fromDegrees(60, 170).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-10, -180).toPoint(), S2LatLng.fromDegrees(10, -180).toPoint())).isTrue()
    }

    @Test
    fun boundaryIntersectsSouthHemisphere() {
        // This rectangle only has only one non-degenerate side.
        val rect = rectFromDegrees(-90, -180, 0, 180)
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-90, -180).toPoint(), S2LatLng.fromDegrees(-60, -180).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-60, -170).toPoint(), S2LatLng.fromDegrees(-60, 170).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(-10, -180).toPoint(), S2LatLng.fromDegrees(10, -180).toPoint())).isTrue()
    }

    @Test
    fun boundaryIntersectsRectCrossingAntiMeridian() {
        val rect = rectFromDegrees(20, 170, 40, -170)
        assertThat(rect.contains(S2LatLng.fromDegrees(30, 180))).isTrue()

        // Check that crossings of all four sides are detected.
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(25, 160).toPoint(), S2LatLng.fromDegrees(25, 180).toPoint())).isTrue()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(25, -160).toPoint(), S2LatLng.fromDegrees(25, -180).toPoint())).isTrue()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(15, 175).toPoint(), S2LatLng.fromDegrees(30, 175).toPoint())).isTrue()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(45, 175).toPoint(), S2LatLng.fromDegrees(30, 175).toPoint())).isTrue()

        // Check that the edges on the opposite side of the sphere but at the same
        // latitude do not intersect the rectangle boundary.
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(25, -20).toPoint(), S2LatLng.fromDegrees(25, 0).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(25, 20).toPoint(), S2LatLng.fromDegrees(25, 0).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(15, -5).toPoint(), S2LatLng.fromDegrees(30, -5).toPoint())).isFalse()
        assertThat(rect.boundaryIntersects(S2LatLng.fromDegrees(45, -5).toPoint(), S2LatLng.fromDegrees(30, -5).toPoint())).isFalse()
    }

    @Test
    fun addPoint() {
        var p = S2LatLngRect.empty()
        p = p.addPoint(S2LatLng.fromDegrees(0, 0))
        assertThat(p.isPoint).isTrue()
        p = p.addPoint(S2LatLng.fromRadians(0.0, -M_PI_2))
        assertThat(p.isPoint).isFalse()
        p = p.addPoint(S2LatLng.fromRadians(M_PI_4, -M_PI))
        p = p.addPoint(S2Point(0, 0, 1))
        assertThat(p).isEqualTo(rectFromDegrees(0, -180, 90, 0))
    }

    @Test
    fun expanded() {
        assertThat(rectFromDegrees(70, 150, 80, 170).expanded(S2LatLng.fromDegrees(20, 30)).approxEquals(rectFromDegrees(50, 120, 90, -160))).isTrue()
        assertThat(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(20, 30)).isEmpty).isTrue()
        assertThat(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(500, 500)).isFull).isTrue()
        assertThat(rectFromDegrees(-90, 170, 10, 20).expanded(S2LatLng.fromDegrees(30, 80)).approxEquals(rectFromDegrees(-90, -180, 40, 180))).isTrue()

        // Negative margins.
        assertThat(rectFromDegrees(10, -50, 60, 70).expanded(S2LatLng.fromDegrees(-10, -10)).approxEquals(rectFromDegrees(20, -40, 50, 60))).isTrue()
        assertThat(rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(-10, -10)).approxEquals(rectFromDegrees(-10, -180, 10, 180))).isTrue()
        assertThat(rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(-30, -30)).isEmpty).isTrue()
        assertThat(rectFromDegrees(-90, 10, 90, 11).expanded(S2LatLng.fromDegrees(-10, -10)).isEmpty).isTrue()
        assertThat(rectFromDegrees(-90, 10, 90, 100).expanded(S2LatLng.fromDegrees(-10, -10)).approxEquals(rectFromDegrees(-80, 20, 80, 90))).isTrue()
        assertThat(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(-50, -500)).isEmpty).isTrue()
        assertThat(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(-50, -50)).approxEquals(rectFromDegrees(-40, -180, 40, 180))).isTrue()

        // Mixed margins.
        assertThat(rectFromDegrees(10, -50, 60, 70).expanded(S2LatLng.fromDegrees(-10, 30)).approxEquals(rectFromDegrees(20, -80, 50, 100))).isTrue()
        assertThat(rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(10, -500)).approxEquals(rectFromDegrees(-30, -180, 30, 180))).isTrue()
        assertThat(rectFromDegrees(-90, -180, 80, 180).expanded(S2LatLng.fromDegrees(-30, 500)).approxEquals(rectFromDegrees(-60, -180, 50, 180))).isTrue()
        assertThat(rectFromDegrees(-80, -100, 80, 150).expanded(S2LatLng.fromDegrees(30, -50)).approxEquals(rectFromDegrees(-90, -50, 90, 100))).isTrue()
        assertThat(rectFromDegrees(0, -180, 50, 180).expanded(S2LatLng.fromDegrees(-30, 500)).isEmpty).isTrue()
        assertThat(rectFromDegrees(-80, 10, 70, 20).expanded(S2LatLng.fromDegrees(30, -200)).isEmpty).isTrue()
        assertThat(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(100, -100)).isEmpty).isTrue()
        assertThat(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(100, -100)).isFull).isTrue()

    }

    @Test
    fun polarClosure() {
        assertThat(rectFromDegrees(-89, 0, 89, 1).polarClosure()).isEqualTo(rectFromDegrees(-89, 0, 89, 1))
        assertThat(rectFromDegrees(-90, -30, -45, 100).polarClosure()).isEqualTo(rectFromDegrees(-90, -180, -45, 180))
        assertThat(rectFromDegrees(89, 145, 90, 146).polarClosure()).isEqualTo(rectFromDegrees(89, -180, 90, 180))
        assertThat(rectFromDegrees(-90, -145, 90, -144).polarClosure()).isEqualTo(S2LatLngRect.full())
    }

    @Test
    fun expandedByDistancePositiveDistance() {
        assertThat(rectFromDegrees(0, 170, 0, -170).expandedByDistance(S1Angle.degrees(15)).approxEquals(rectFromDegrees(-15, 155, 15, -155))).isTrue()
        assertThat(rectFromDegrees(60, 150, 80, 10).expandedByDistance(S1Angle.degrees(15)).approxEquals(rectFromDegrees(45, -180, 90, 180))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceNorthEast() {
        val inRect = rectFromDegrees(0.0, 0.0, 30.0, 90.0)
        val distance = S1Angle.degrees(5.0)
        val outRect = inRect.expandedByDistance(distance).expandedByDistance(-distance)
        assertThat(outRect.approxEquals(inRect)).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceSouthWest() {
        val inRect = rectFromDegrees(-30.0, -90.0, 0.0, 0.0)
        val distance = S1Angle.degrees(5.0)
        val outRect = inRect.expandedByDistance(distance).expandedByDistance(-distance)
        assertThat(outRect.approxEquals(inRect)).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLatWithNorthPole() {
        val rect = rectFromDegrees(0.0, -90.0, 90.0, 180.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.approxEquals(rectFromDegrees(5.0, 0.0, 85.0, 90.0))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLatWithNorthPoleAndLngFull() {
        val rect = rectFromDegrees(0.0, -180.0, 90.0, 180.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.approxEquals(rectFromDegrees(5.0, -180.0, 90.0, 180.0))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLatWithSouthPole() {
        val rect = rectFromDegrees(-90.0, -90.0, 0.0, 180.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.approxEquals(rectFromDegrees(-85.0, 0.0, -5.0, 90.0))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLatWithSouthPoleAndLngFull() {
        val rect = rectFromDegrees(-90.0, -180.0, 0.0, 180.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.approxEquals(rectFromDegrees(-90.0, -180.0, -5.0, 180.0))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLngFull() {
        val rect = rectFromDegrees(0.0, -180.0, 30.0, 180.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.approxEquals(rectFromDegrees(5.0, -180.0, 25.0, 180.0))).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLatResultEmpty() {
        val rect = rectFromDegrees(0.0, 0.0, 9.9, 90.0).expandedByDistance(-S1Angle.degrees(5.0))
        assertThat(rect.isEmpty).isTrue()
    }

    @Test
    fun expandedByDistanceNegativeDistanceLngResultEmpty() {
        val rect = rectFromDegrees(0.0, 0.0, 30.0, 11.0).expandedByDistance(-S1Angle.degrees(5.0))

        // The cap center is at latitude 30 - 5 = 25 degrees. The length of the
        // latitude 25 degree line is 0.906 times the length of the equator. Thus the
        // cap whose radius is 5 degrees covers the rectangle whose latitude interval
        // is 11 degrees.
        assertThat(rect.isEmpty).isTrue()
    }

    @Test
    fun getCapBound() {
        // Bounding cap at center is smaller:
        assertThat(rectFromDegrees(-45, -45, 45, 45).capBound.approxEquals(S2Cap.fromCenterHeight(S2Point(1, 0, 0), 0.5))).isTrue()

        // Bounding cap at north pole is smaller:
        assertThat(rectFromDegrees(88, -80, 89, 80).capBound.approxEquals(S2Cap.fromCenterAngle(S2Point(0, 0, 1), S1Angle.degrees(2)))).isTrue()

        // Longitude span > 180 degrees:
        assertThat(rectFromDegrees(-30, -150, -10, 50).capBound.approxEquals(S2Cap.fromCenterAngle(S2Point(0, 0, -1), S1Angle.degrees(80)))).isTrue()
    }

    @Test
    fun cellOps() {
        // Contains(S2Cell), MayIntersect(S2Cell), Intersects(S2Cell)

        // Special cases.
        testCellOps(S2LatLngRect.empty(), S2Cell.fromFacePosLevel(3, 0, 0), 0)
        testCellOps(S2LatLngRect.full(), S2Cell.fromFacePosLevel(2, 0, 0), 4)
        testCellOps(S2LatLngRect.full(), S2Cell.fromFacePosLevel(5, 0, 25), 4)

        // This rectangle includes the first quadrant of face 0.  It's expanded
        // slightly because cell bounding rectangles are slightly conservative.
        val r4 = rectFromDegrees(-45.1, -45.1, 0.1, 0.1)
        testCellOps(r4, S2Cell.fromFacePosLevel(0, 0, 0), 3)
        testCellOps(r4, S2Cell.fromFacePosLevel(0, 0, 1), 4)
        testCellOps(r4, S2Cell.fromFacePosLevel(1, 0, 1), 0)

        // This rectangle intersects the first quadrant of face 0.
        val r5 = rectFromDegrees(-10, -45, 10, 0)
        testCellOps(r5, S2Cell.fromFacePosLevel(0, 0, 0), 3)
        testCellOps(r5, S2Cell.fromFacePosLevel(0, 0, 1), 3)
        testCellOps(r5, S2Cell.fromFacePosLevel(1, 0, 1), 0)

        // Rectangle consisting of a single point.
        testCellOps(rectFromDegrees(4, 4, 4, 4), S2Cell.fromFace(0), 3)

        // Rectangles that intersect the bounding rectangle of a face
        // but not the face itself.
        testCellOps(rectFromDegrees(41, -87, 42, -79), S2Cell.fromFace(2), 1)
        testCellOps(rectFromDegrees(-41, 160, -40, -160), S2Cell.fromFace(5), 1)

        // This is the leaf cell at the top right hand corner of face 0.
        // It has two angles of 60 degrees and two of 120 degrees.
        val cell0tr = S2Cell(S2Point(1 + 1e-12, 1.0, 1.0))
        val bound0tr = cell0tr.rectBound
        val v0 = S2LatLng.fromPoint(cell0tr.getVertexRaw(0))
        testCellOps(rectFromDegrees(v0.lat().degrees() - 1e-8,
                v0.lng().degrees() - 1e-8,
                v0.lat().degrees() - 2e-10,
                v0.lng().degrees() + 1e-10),
                cell0tr, 1)

        // Rectangles that intersect a face but where no vertex of one region
        // is contained by the other region.  The first one passes through
        // a corner of one of the face cells.
        testCellOps(rectFromDegrees(-37, -70, -36, -20), S2Cell.fromFace(5), 2)

        // These two intersect like a diamond and a square.
        val cell202 = S2Cell.fromFacePosLevel(2, 0, 2)
        val bound202 = cell202.rectBound
        testCellOps(rectFromDegrees(bound202.lo().lat().degrees() + 3,
                bound202.lo().lng().degrees() + 3,
                bound202.hi().lat().degrees() - 3,
                bound202.hi().lng().degrees() - 3),
                cell202, 2)
    }

    @Test
    fun area() {
        assertThat(S2LatLngRect.empty().area).isEqualTo(0.0)
        assertThat(S2LatLngRect.full().area).isEqualTo(4 * M_PI)
        assertThat(rectFromDegrees(0, 0, 90, 90).area).isEqualTo(M_PI_2)
    }

    @Test
    fun getCentroid() {
        // Empty and full rectangles.
        assertThat(S2LatLngRect.empty().centroid).isEqualTo(S2Point())
        assertThat(S2LatLngRect.full().centroid.norm() <= 1e-15).isTrue()

        // Rectangles that cover the full longitude range.
        repeat (100) {
            val lat1 = Random.nextDouble(-M_PI_2, M_PI_2)
            val lat2 = Random.nextDouble(-M_PI_2, M_PI_2)
            val r = S2LatLngRect(R1Interval.fromPointPair(lat1, lat2), S1Interval.full())
            val centroid = r.centroid
            assertThat(centroid.z).isCloseTo(0.5 * (sin(lat1) + sin(lat2)) * r.area, Offset.offset(1e-15))
            assertThat(R2VectorDouble(centroid.x, centroid.y).norm() <= 1e-15).isTrue()
        }

        // Rectangles that cover the full latitude range.
        repeat(100) {
            val lng1 = Random.nextDouble(-M_PI, M_PI)
            val lng2 = Random.nextDouble(-M_PI, M_PI)
            val r = S2LatLngRect(S2LatLngRect.fullLat(), S1Interval.fromPointPair(lng1, lng2))
            val centroid = r.centroid
            assertThat(abs(centroid.z) <= 1e-15).isTrue()
            assertThat(S2LatLng.fromPoint(centroid).lng().radians).isCloseTo(r.lng.center, Offset.offset(1e-15))
            val alpha = 0.5 * r.lng.length
            assertThat(R2VectorDouble(centroid.x, centroid.y).norm()).isCloseTo(0.25 * M_PI * sin(alpha) / alpha * r.area, Offset.offset(1e-15))
        }

        // Finally, verify that when a rectangle is recursively split into pieces,
        // the centroids of the pieces add to give the centroid of their parent.
        // To make the code simpler we avoid rectangles that cross the 180 degree
        // line of longitude.
        testCentroidSplitting(S2LatLngRect(S2LatLngRect.fullLat(), S1Interval(-3.14, 3.14)), 10 /*splits_left*/)
    }

    @Test
    fun getDistanceOverlapping() {
        // Check pairs of rectangles that overlap: (should all return 0):
        val a = rectFromDegrees(0, 0, 2, 2)
        val b = pointrectFromDegrees(0.0, 0.0)
        assertThat(a.getDistance(a)).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(b)).isEqualTo(S1Angle.radians(0))
        assertThat(b.getDistance(b)).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(S2LatLng.fromDegrees(0, 0))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(0, 1, 2, 3))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(0, 2, 2, 4))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(1, 0, 3, 2))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(2, 0, 4, 2))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(1, 1, 3, 3))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDistance(rectFromDegrees(2, 2, 4, 4))).isEqualTo(S1Angle.radians(0))
    }

    @Test
    fun getDistanceRectVsPoint() {
        // Rect that spans 180.
        val a = rectFromDegrees(-1, -1, 2, 1)
        verifyGetDistance(a, pointrectFromDegrees(-2, -1))
        verifyGetDistance(a, pointrectFromDegrees(1, 2))

        verifyGetDistance(pointrectFromDegrees(-2, -1), a)
        verifyGetDistance(pointrectFromDegrees(1, 2), a)

        verifyGetRectPointDistance(a, S2LatLng.fromDegrees(-2, -1))
        verifyGetRectPointDistance(a, S2LatLng.fromDegrees(1, 2))

        // Tests near the north pole.
        val b = rectFromDegrees(86, 0, 88, 2)
        verifyGetDistance(b, pointrectFromDegrees(87, 3))
        verifyGetDistance(b, pointrectFromDegrees(87, -1))
        verifyGetDistance(b, pointrectFromDegrees(89, 1))
        verifyGetDistance(b, pointrectFromDegrees(89, 181))
        verifyGetDistance(b, pointrectFromDegrees(85, 1))
        verifyGetDistance(b, pointrectFromDegrees(85, 181))
        verifyGetDistance(b, pointrectFromDegrees(90, 0))

        verifyGetDistance(pointrectFromDegrees(87, 3), b)
        verifyGetDistance(pointrectFromDegrees(87, -1), b)
        verifyGetDistance(pointrectFromDegrees(89, 1), b)
        verifyGetDistance(pointrectFromDegrees(89, 181), b)
        verifyGetDistance(pointrectFromDegrees(85, 1), b)
        verifyGetDistance(pointrectFromDegrees(85, 181), b)
        verifyGetDistance(pointrectFromDegrees(90, 0), b)

        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, 3))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, -1))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 1))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 181))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 1))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 181))
        verifyGetRectPointDistance(b, S2LatLng.fromDegrees(90, 0))

        // Rect that touches the north pole.
        val c = rectFromDegrees(88, 0, 90, 2)
        verifyGetDistance(c, pointrectFromDegrees(89, 3))
        verifyGetDistance(c, pointrectFromDegrees(89, 90))
        verifyGetDistance(c, pointrectFromDegrees(89, 181))
        verifyGetDistance(pointrectFromDegrees(89, 3), c)
        verifyGetDistance(pointrectFromDegrees(89, 90), c)
        verifyGetDistance(pointrectFromDegrees(89, 181), c)
    }

    @Test
    fun getDistanceRectVsRect() {
        // Rect that spans 180.
        val a = rectFromDegrees(-1, -1, 2, 1)
        verifyGetDistance(a, rectFromDegrees(0, 2, 1, 3))
        verifyGetDistance(a, rectFromDegrees(-2, -3, -1, -2))

        // Tests near the south pole.
        val b = rectFromDegrees(-87, 0, -85, 3)
        verifyGetDistance(b, rectFromDegrees(-89, 1, -88, 2))
        verifyGetDistance(b, rectFromDegrees(-84, 1, -83, 2))
        verifyGetDistance(b, rectFromDegrees(-88, 90, -86, 91))
        verifyGetDistance(b, rectFromDegrees(-84, -91, -83, -90))
        verifyGetDistance(b, rectFromDegrees(-90, 181, -89, 182))
        verifyGetDistance(b, rectFromDegrees(-84, 181, -83, 182))
    }

    @Test
    fun getDistanceRandomPairs() {
        // Test random pairs.
        for (i in 0..9999) {
            val a = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(randomPoint()), S2LatLng.fromPoint(randomPoint()))
            val b = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(randomPoint()), S2LatLng.fromPoint(randomPoint()))
            verifyGetDistance(a, b)

            val c = S2LatLng.fromPoint(randomPoint())
            verifyGetRectPointDistance(a, c)
            verifyGetRectPointDistance(b, c)
        }
    }

    @Test
    fun getDirectedHausdorffDistanceRandomPairs() {
        // Test random pairs.
        val kIters = 1000
        for (i in 1..kIters) {
            val a = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(randomPoint()), S2LatLng.fromPoint(randomPoint()))
            val b = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(randomPoint()), S2LatLng.fromPoint(randomPoint()))
            // a and b are *minimum* bounding rectangles of two random points, in
            // particular, their Voronoi diagrams are always of the same topology. We
            // take the "complements" of a and b for more thorough testing.
            val a2 = S2LatLngRect(a.lat, a.lng.complement)
            val b2 = S2LatLngRect(b.lat, b.lng.complement)

            // Note that "a" and "b" come from the same distribution, so there is no
            // need to test pairs such as (b, a), (b, a2), etc.
            verifyGetDirectedHausdorffDistance(a, b)
            verifyGetDirectedHausdorffDistance(a, b2)
            verifyGetDirectedHausdorffDistance(a2, b)
            verifyGetDirectedHausdorffDistance(a2, b2)
        }
    }

    @Test
    fun getDirectedHausdorffDistanceContained() {
        // Caller rect is contained in callee rect. Should return 0.
        val a = rectFromDegrees(-10, 20, -5, 90)
        assertThat(a.getDirectedHausdorffDistance(rectFromDegrees(-10, 20, -5, 90))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDirectedHausdorffDistance(rectFromDegrees(-10, 19, -5, 91))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDirectedHausdorffDistance(rectFromDegrees(-11, 20, -4, 90))).isEqualTo(S1Angle.radians(0))
        assertThat(a.getDirectedHausdorffDistance(rectFromDegrees(-11, 19, -4, 91))).isEqualTo(S1Angle.radians(0))
    }

    @Test
    fun getDirectHausdorffDistancePointToRect() {
        // The Hausdorff distance from a point to a rect should be the same as its
        // distance to the rect.
        val a1 = pointrectFromDegrees (5, 8)
        val a2 = pointrectFromDegrees (90, 10)  // north pole

        var b = rectFromDegrees (-85, -50, -80, 10)
        assertThat(a1.getDirectedHausdorffDistance(b).radians).isEqualTo(a1.getDistance(b).radians)
        assertThat(a2.getDirectedHausdorffDistance(b).radians).isEqualTo(a2.getDistance(b).radians)

        b = rectFromDegrees(4, -10, 80, 10)
        assertThat(a1.getDirectedHausdorffDistance(b).radians).isEqualTo(a1.getDistance(b).radians)
        assertThat(a2.getDirectedHausdorffDistance(b).radians).isEqualTo(a2.getDistance(b).radians)

        b = rectFromDegrees(70, 170, 80, -170)
        assertThat(a1.getDirectedHausdorffDistance(b).radians).isCloseTo(a1.getDistance(b).radians, Offset.offset(1e-15))
        assertThat(a2.getDirectedHausdorffDistance(b).radians).isCloseTo(a2.getDistance(b).radians, Offset.offset(1e-15))
    }

    @Test
    fun getDirectedHausdorffDistanceRectToPoint() {
        val a = rectFromDegrees (1, -8, 10, 20)
        verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(5, 8))
        verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(-6, -100))
        // south pole
        verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(-90, -20))
        // north pole
        verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(90, 0))
    }

    @Test
    fun getDirectedHausdorffDistanceRectToRectNearPole() {
        // Tests near south pole.
        val a = rectFromDegrees (-87, 0, -85, 3)
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-89, 1, -88, 2))
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 1, -83, 2))
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-88, 90, -86, 91))
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, -91, -83, -90))
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-90, 181, -89, 182))
        verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 181, -83, 182))
    }

    @Test
    fun getDirectedHausdorffDistanceRectToRectDegenerateCases() {
        // Rectangles that contain poles.
        verifyGetDirectedHausdorffDistance(rectFromDegrees(0, 10, 90, 20), rectFromDegrees(-4, -10, 4, 0))
        verifyGetDirectedHausdorffDistance(rectFromDegrees(-4, -10, 4, 0), rectFromDegrees(0, 10, 90, 20))

        // Two rectangles share same or complement longitudinal intervals.
        val a = rectFromDegrees (-50, -10, 50, 10)
        val b = rectFromDegrees (30, -10, 60, 10)
        verifyGetDirectedHausdorffDistance(a, b)
        val c = S2LatLngRect(a.lat, a.lng.complement)
        verifyGetDirectedHausdorffDistance(c, b)

        // rectangle a touches b_opposite_lng.
        verifyGetDirectedHausdorffDistance(rectFromDegrees(10, 170, 30, 180), rectFromDegrees(-50, -10, 50, 10))
        verifyGetDirectedHausdorffDistance(rectFromDegrees(10, -180, 30, -170), rectFromDegrees(-50, -10, 50, 10))

        // rectangle b's Voronoi diagram is degenerate (lng interval spans 180
        // degrees), and a touches the degenerate Voronoi vertex.
        verifyGetDirectedHausdorffDistance(rectFromDegrees(-30, 170, 30, 180), rectFromDegrees(-10, -90, 10, 90))
        verifyGetDirectedHausdorffDistance(rectFromDegrees(-30, -180, 30, -170), rectFromDegrees(-10, -90, 10, 90))

        // rectangle a touches a voronoi vertex of rectangle b.
        verifyGetDirectedHausdorffDistance(rectFromDegrees(-20, 105, 20, 110), rectFromDegrees(-30, 5, 30, 15))
        verifyGetDirectedHausdorffDistance(rectFromDegrees(-20, 95, 20, 105), rectFromDegrees(-30, 5, 30, 15))
    }


    companion object {

        fun testIntervalOps(x: S2LatLngRect, y: S2LatLngRect, expected_relation: String, expected_union: S2LatLngRect, expected_intersection: S2LatLngRect) {
            // Test all of the interval operations on the given pair of intervals.
            // "expected_relation" is a sequence of "T" and "F" characters corresponding
            // to the expected results of Contains(), InteriorContains(), Intersects(),
            // and InteriorIntersects() respectively.

            assertThat(x.contains(y)).isEqualTo(expected_relation[0] == 'T')
            assertThat(x.interiorContains(y)).isEqualTo(expected_relation[1] == 'T')
            assertThat(x.intersects(y)).isEqualTo(expected_relation[2] == 'T')
            assertThat(x.interiorIntersects(y)).isEqualTo(expected_relation[3] == 'T')

            assertThat(x.contains(y)).isEqualTo(x.union(y) == x)
            assertThat(x.intersects(y)).isEqualTo(!x.intersection(y).isEmpty)

            assertThat(x.union(y)).isEqualTo(expected_union)
            assertThat(x.intersection(y)).isEqualTo(expected_intersection)

            if (y.size == S2LatLng.fromRadians(0.0, 0.0)) {
                var r = x
                r = r.addPoint(y.lo())
                assertThat(r).isEqualTo(expected_union)
            }
        }

        fun testCellOps(r: S2LatLngRect, cell: S2Cell, level: Int) {
            // Test the relationship between the given rectangle and cell:
            // 0 == no intersection, 1 == MayIntersect, 2 == Intersects,
            // 3 == Vertex Containment, 4 == Contains

            var vertex_contained = false
            for (i in 0..3) {
                if (r.contains(cell.getVertexRaw(i)) || (!r.isEmpty && cell.contains(r.getVertex(i).toPoint())))
                    vertex_contained = true
            }
            assertThat(r.mayIntersect(cell)).isEqualTo(level >= 1)
            assertThat(r.intersects(cell)).isEqualTo(level >= 2)
            assertThat(vertex_contained).isEqualTo(level >= 3)
            assertThat(r.contains(cell)).isEqualTo(level >= 4)
        }


        private fun rectFromDegrees(latLo: Double, lngLo: Double, latHi: Double, lngHi: Double): S2LatLngRect {
            // Convenience method to construct a rectangle. This method is
            // intentionally *not* in the S2LatLngRect interface because the
            // argument order is ambiguous, but hopefully it's not too confusing
            // within the context of this unit test.
            return S2LatLngRect(S2LatLng.fromDegrees(latLo, lngLo).normalized(), S2LatLng.fromDegrees(latHi, lngHi).normalized())
        }

        private fun rectFromDegrees(latLo: Int, lngLo: Int, latHi: Int, lngHi: Int): S2LatLngRect =
                rectFromDegrees(latLo.toDouble(), lngLo.toDouble(), latHi.toDouble(), lngHi.toDouble())

        // Recursively verify that when a rectangle is split into two pieces, the
        // centroids of the children sum to give the centroid of the parent.
        fun testCentroidSplitting(r: S2LatLngRect, splits_left: Int) {
            val child0: S2LatLngRect
            val child1: S2LatLngRect

            if (Random.nextBoolean()) {
                val lat = Random.nextDouble(r.lat.lo, r.lat.hi)
                child0 = S2LatLngRect(R1Interval(r.lat.lo, lat), r.lng)
                child1 = S2LatLngRect(R1Interval(lat, r.lat.hi), r.lng)
            } else {
                assertThat(r.lng.lo <= r.lng.hi).isTrue()
                val lng = Random.nextDouble(r.lng.lo, r.lng.hi)
                child0 = S2LatLngRect(r.lat, S1Interval(r.lng.lo, lng))
                child1 = S2LatLngRect(r.lat, S1Interval(lng, r.lng.hi))
            }
            assertThat((r.centroid - child0.centroid - child1.centroid).norm()).isCloseTo(0.0, Offset.offset(1e-15))
            if (splits_left > 0) {
                testCentroidSplitting(child0, splits_left - 1)
                testCentroidSplitting(child1, splits_left - 1)
            }
        }

        // Returns the minimum distance from X to the latitude line segment defined by
        // the given latitude and longitude interval.
        fun GetDistance(x: S2LatLng, lat: S1Angle, interval: S1Interval): S1Angle {
            assertThat(x.isValid).isTrue()
            assertThat(interval.isValid).isTrue()

            // Is X inside the longitude interval?
            if (interval.contains(x.lng().radians))
                return (x.lat() - lat).abs()

            // Return the distance to the closer endpoint.
            return min(
                    x.getDistance(S2LatLng.fromLatLng(lat, S1Angle.radians(interval.lo))),
                    x.getDistance(S2LatLng.fromLatLng(lat, S1Angle.radians(interval.hi)))
            )
        }

        fun bruteForceDistance(a: S2LatLngRect, b: S2LatLngRect): S1Angle {
            if (a.intersects(b))
                return S1Angle.radians(0)

            // Compare every point in 'a' against every latitude edge and longitude edge
            // in 'b', and vice-versa, for a total of 16 point-vs-latitude-edge tests and
            // 16 point-vs-longitude-edge tests.
            val pnt_a = Array<S2LatLng>(4) { return@Array S2LatLng() }
            val pnt_b = Array<S2LatLng>(4) { S2LatLng() }
            pnt_a[0] = S2LatLng.fromLatLng(a.latLo(), a.lngLo())
            pnt_a[1] = S2LatLng.fromLatLng(a.latLo(), a.lngHi())
            pnt_a[2] = S2LatLng.fromLatLng(a.latHi(), a.lngHi())
            pnt_a[3] = S2LatLng.fromLatLng(a.latHi(), a.lngLo())
            pnt_b[0] = S2LatLng.fromLatLng(b.latLo(), b.lngLo())
            pnt_b[1] = S2LatLng.fromLatLng(b.latLo(), b.lngHi())
            pnt_b[2] = S2LatLng.fromLatLng(b.latHi(), b.lngHi())
            pnt_b[3] = S2LatLng.fromLatLng(b.latHi(), b.lngLo())

            // Make arrays containing the lo/hi latitudes and the lo/hi longitude edges.
            val lat_a = arrayOf(a.latLo(), a.latHi())
            val lat_b = arrayOf(b.latLo(), b.latHi())
            val lng_edge_a = arrayOf(
                    arrayOf(pnt_a[0].toPoint(), pnt_a[3].toPoint()),
                    arrayOf(pnt_a[1].toPoint(), pnt_a[2].toPoint())
            )
            val lng_edge_b = arrayOf(
                    arrayOf(pnt_b[0].toPoint(), pnt_b[3].toPoint()),
                    arrayOf(pnt_b[1].toPoint(), pnt_b[2].toPoint())
            )

            var min_distance = S1Angle.degrees(180.0)
            for (i in 0..3) {
                // For each point in a and b.
                val current_a = pnt_a[i]
                val current_b = pnt_b[i]

                for (j in 0..1) {
                    // Get distances to latitude and longitude edges.
                    val a_to_lat = GetDistance(current_a, lat_b[j], b.lng)
                    val b_to_lat = GetDistance(current_b, lat_a[j], a.lng)
                    val a_to_lng = S2EdgeDistances.getDistance(current_a.toPoint(), lng_edge_b[j][0], lng_edge_b[j][1])
                    val b_to_lng = S2EdgeDistances.getDistance(current_b.toPoint(), lng_edge_a[j][0], lng_edge_a[j][1])

                    min_distance = minOf(min_distance, minOf(a_to_lat, minOf(b_to_lat, minOf(a_to_lng, b_to_lng))))
                }
            }
            return min_distance
        }

        fun BruteForceRectPointDistance(a: S2LatLngRect, b: S2LatLng): S1Angle {
            if (a.contains(b)) {
                return S1Angle.radians(0)
            }

            val b_to_lo_lat = GetDistance(b, a.latLo(), a.lng)
            val b_to_hi_lat = GetDistance(b, a.latHi(), a.lng)
            val b_to_lo_lng = S2EdgeDistances.getDistance(b.toPoint(), S2LatLng.fromLatLng(a.latLo(), a.lngLo()).toPoint(), S2LatLng.fromLatLng(a.latHi(), a.lngLo()).toPoint())
            val b_to_hi_lng = S2EdgeDistances.getDistance(b.toPoint(), S2LatLng.fromLatLng(a.latLo(), a.lngHi()).toPoint(), S2LatLng.fromLatLng(a.latHi(), a.lngHi()).toPoint())
            return min(b_to_lo_lat, min(b_to_hi_lat, min(b_to_lo_lng, b_to_hi_lng)))
        }

        // This method verifies a.getDistance(b) by comparing its result against a
        // brute-force implementation. The correctness of the brute-force version is
        // much easier to verify by inspection.
        fun verifyGetDistance(a: S2LatLngRect, b: S2LatLngRect) {
            val distance1 = bruteForceDistance(a, b)
            val distance2 = a.getDistance(b)
            assertThat(distance1.radians - distance2.radians).isCloseTo(0.0, Offset.offset(1e-10))
        }

        fun pointrectFromDegrees(lat: Double, lng: Double): S2LatLngRect {
            return S2LatLngRect.fromPoint(S2LatLng.fromDegrees(lat, lng).normalized())
        }

        fun pointrectFromDegrees(lat: Int, lng: Int): S2LatLngRect = pointrectFromDegrees(lat.toDouble(), lng.toDouble())

        // This method verifies a.getDistance(b), where b is a S2LatLng, by comparing
        // its result against a.getDistance(c), c being the point rectangle created
        // from b.
        fun verifyGetRectPointDistance(a: S2LatLngRect, p: S2LatLng) {
            val distance1 = BruteForceRectPointDistance(a, p.normalized())
            val distance2 = a.getDistance(p.normalized())
            assertThat(abs(distance1.radians - distance2.radians)).isCloseTo( 0.0, Offset.offset(1e-10))
        }

        // This function assumes that GetDirectedHausdorffDistance() always returns
        // a distance from some point in a to b. So the function mainly tests whether
        // the returned distance is large enough, and only does a weak test on whether
        // it is small enough.
        fun verifyGetDirectedHausdorffDistance(a: S2LatLngRect, b: S2LatLngRect) {
            val hausdorff_distance = a.getDirectedHausdorffDistance(b)

            val kResolution = 0.1
            // Record the max sample distance as well as the sample point realizing the
            // max for easier debugging.
            var max_distance: S1Angle = zero()
            var lat_max: Double
            var lng_max: Double

            val sample_size_on_lat = ((a.lat.length / kResolution) + 1).toInt()
            val sample_size_on_lng = ((a.lng.length / kResolution) + 1).toInt()
            val delta_on_lat = a.lat.length / sample_size_on_lat
            val delta_on_lng = a.lng.length / sample_size_on_lng

            var lng = a.lng.lo
            for (i in 0..sample_size_on_lng) {
                var lat = a.lat.lo
                for (j in 0..sample_size_on_lat) {
                    val latlng = S2LatLng.fromRadians(lat, lng).normalized()
                    val distance_to_b = b.getDistance(latlng)

                    if (distance_to_b >= max_distance) {
                        max_distance = distance_to_b
                        lat_max = lat
                        lng_max = lng
                    }
                    lat += delta_on_lat
                }
                lng += delta_on_lng
            }

            assertThat(max_distance.radians <= hausdorff_distance.radians + 1e-10).isTrue()
            assertThat(max_distance.radians >= hausdorff_distance.radians - kResolution).isTrue()
        }
    }
}
