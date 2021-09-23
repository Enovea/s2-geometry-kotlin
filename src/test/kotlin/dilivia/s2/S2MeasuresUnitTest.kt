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
import dilivia.math.vectors.times
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test

class S2MeasuresUnitTest {

    private val logger = KotlinLogging.logger { }

    @Test
    fun angleMethods() {
        val pz = S2Point(0, 0, 1)
        val p000 = S2Point(1, 0, 0)
        val p045 = S2Point(1, 1, 0).normalize()
        val p090 = S2Point(0, 1, 0)
        val p180 = S2Point(-1, 0, 0)

        assertThat(S2Measures.angle(p000, pz, p045)).isEqualTo(M_PI_4)
        assertThat(S2Measures.turnAngle(p000, pz, p045)).isEqualTo(-3 * M_PI_4)

        assertThat(S2Measures.angle(p045, pz, p180)).isEqualTo(3 * M_PI_4)
        assertThat(S2Measures.turnAngle(p045, pz, p180)).isEqualTo(-M_PI_4)

        assertThat(S2Measures.angle(p000, pz, p180)).isEqualTo(M_PI)
        assertThat(S2Measures.turnAngle(p000, pz, p180)).isEqualTo(-0.0)

        assertThat(S2Measures.angle(pz, p000, p045)).isEqualTo(M_PI_2)
        assertThat(S2Measures.turnAngle(pz, p000, p045)).isEqualTo(M_PI_2)

        assertThat(S2Measures.angle(pz, p000, pz)).isEqualTo(0.0)
        assertThat(abs(S2Measures.turnAngle(pz, p000, pz))).isEqualTo(M_PI)
    }

    @Test
    fun areaMethods() {
        val pz = S2Point(0, 0, 1)
        val p000 = S2Point(1, 0, 0)
        val p045 = S2Point(1, 1, 0).normalize()
        val p090 = S2Point(0, 1, 0)
        val p180 = S2Point(-1, 0, 0)

        assertThat(S2Measures.area(p000, p090, pz)).isEqualTo(M_PI_2)
        assertThat(S2Measures.area(p045, pz, p180)).isEqualTo(3 * M_PI_4)

        // Make sure that Area() has good *relative* accuracy even for
        // very small areas.
        val eps = 1e-10
        val pepsx = S2Point(eps, 0.0, 1.0).normalize()
        val pepsy = S2Point(0.0, eps, 1.0).normalize()
        val expected1 = 0.5 * eps * eps
        assertThat(S2Measures.area(pepsx, pepsy, pz)).isCloseTo(expected1, Offset.offset(1e-14 * expected1))

        // Make sure that it can handle degenerate triangles.
        val pr = S2Point(0.257, -0.5723, 0.112).normalize()
        val pq = S2Point(-0.747, 0.401, 0.2235).normalize()
        assertThat(S2Measures.area(pr, pr, pr)).isEqualTo(0.0)
        // The following test is not exact due to rounding error.
        assertThat(S2Measures.area(pr, pq, pr)).isCloseTo(0.0, Offset.offset(1e-15))
        assertThat(S2Measures.area(p000, p045, p090)).isEqualTo(0.0)

        var maxGirard = 0.0
        repeat(10000) {
            val p0 = S2Random.randomPoint()
            val d1 = S2Random.randomPoint()
            val d2 = S2Random.randomPoint()
            val p1 = (p0 + 1e-15 * d1).normalize()
            val p2 = (p0 + 1e-15 * d2).normalize()
            // The actual displacement can be as much as 1.2e-15 due to roundoff.
            // This yields a maximum triangle area of about 0.7e-30.
            assertThat(S2Measures.area(p0, p1, p2)).isLessThanOrEqualTo(0.7e-30)
            maxGirard = max(maxGirard, S2Measures.girardArea(p0, p1, p2))
        }
        // This check only passes if GirardArea() uses RobustCrossProd().
        logger.info { "Worst case Girard for triangle area 1e-30: $maxGirard" }
        assertThat(maxGirard).isLessThanOrEqualTo(1e-14)

        // Try a very long and skinny triangle.
        val p045eps = S2Point(1.0, 1.0, eps).normalize()
        val expected2 = 5.8578643762690495119753e-11;  // Mathematica.
        assertThat(S2Measures.area(p000, p045eps, p090)).isCloseTo(expected2, Offset.offset(1e-9 * expected2))

        // Triangles with near-180 degree edges that sum to a quarter-sphere.
        val eps2 = 1e-14
        val p000eps2 = S2Point(1.0, 0.1 * eps2, eps2).normalize()
        val quarter_area1 = S2Measures.area(p000eps2, p000, p045) +
                S2Measures.area(p000eps2, p045, p180) +
                S2Measures.area(p000eps2, p180, pz) +
                S2Measures.area(p000eps2, pz, p000)
        assertThat(quarter_area1).isCloseTo(M_PI, Offset.offset(1e-15))

        // Four other triangles that sum to a quarter-sphere.
        val p045eps2 = S2Point(1.0, 1.0, eps2).normalize()
        val quarter_area2 = S2Measures.area(p045eps2, p000, p045) +
                S2Measures.area(p045eps2, p045, p180) +
                S2Measures.area(p045eps2, p180, pz) +
                S2Measures.area(p045eps2, pz, p000)
        assertThat(quarter_area2).isEqualTo(M_PI)

        // Compute the area of a hemisphere using four triangles with one near-180
        // degree edge and one near-degenerate edge.
        repeat(100) {
            val lng = 2 * M_PI * S2Random.randomDouble()
            val p0 = S2LatLng.fromRadians(1e-20, lng).normalized().toPoint()
            val p1 = S2LatLng.fromRadians(0.0, lng).normalized().toPoint()
            val p2_lng = lng + S2Random.randomDouble()
            val p2 = S2LatLng.fromRadians(0.0, p2_lng).normalized().toPoint()
            val p3 = S2LatLng.fromRadians(0.0, lng + M_PI).normalized().toPoint()
            val p4 = S2LatLng.fromRadians(0.0, lng + 5.0).normalized().toPoint()
            val area = (S2Measures.area(p0, p1, p2) + S2Measures.area(p0, p2, p3) + S2Measures.area(
                p0,
                p3,
                p4
            ) + S2Measures.area(p0, p4, p1))
            assertThat(area).isCloseTo(2 * M_PI, Offset.offset(2e-15))
        }

        // This tests a case where the triangle has zero area, but S2Measures.area()
        // computes (dmin > 0) due to rounding errors.
        assertThat(
            S2Measures.area(
                S2LatLng.fromDegrees(-45, -170).toPoint(),
                S2LatLng.fromDegrees(45, -170).toPoint(),
                S2LatLng.fromDegrees(0, -170).toPoint()
            )
        ).isEqualTo(0.0)
    }

}
