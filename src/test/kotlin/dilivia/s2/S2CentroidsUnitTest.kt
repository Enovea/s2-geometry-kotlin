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
import dilivia.math.vectors.times
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomFrame
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test
import kotlin.math.cos
import kotlin.math.pow
import kotlin.math.sin

class S2CentroidsUnitTest {


    @Test
    fun triangleTrueCentroidSmallTriangles() {
        // Test TrueCentroid() with very small triangles.  This test assumes that
        // the triangle is small enough so that it is nearly planar.
        for (iter in 0 until 100) {
            val (p, x, y) = randomFrame()
            val d = 1e-4 * 1e-4.pow(randomDouble())
            val p0 = (p - d * x).normalize()
            val p1 = (p + d * x).normalize()
            val p2 = (p + 3 * d * y).normalize()
            val centroid = S2Centroids.trueCentroid(p0, p1, p2).normalize()

            // The centroid of a planar triangle is at the intersection of its
            // medians, which is two-thirds of the way along each median.
            val expectedCentroid = (p + d * y).normalize()
            assertThat(centroid.angle(expectedCentroid)).isLessThanOrEqualTo(2e-8)
            assertThat(centroid.angle(S2Centroids.planarCentroid(p0, p1, p2).normalize())).isLessThanOrEqualTo(2e-8)
        }
    }

    @Test
    fun edgeTrueCentroidSemiEquator() {
        // Test the centroid of polyline ABC that follows the equator and consists
        // of two 90 degree edges (i.e., C = -A).  The centroid (multiplied by
        // length) should point toward B and have a norm of 2.0.  (The centroid
        // itself has a norm of 2/Pi, and the total edge length is Pi.)
        val a = S2Point(0, -1, 0)
        val b = S2Point(1, 0, 0)
        val c = S2Point(0, 1, 0)
        val centroid = S2Centroids.trueCentroid(a, b) + S2Centroids.trueCentroid(b, c)
        println(S2Centroids.trueCentroid(a, b))
        println(S2Centroids.trueCentroid(b, c))
        assertThat(S2PointUtil.approxEquals(b, centroid.normalized())).isTrue()
        println(centroid)
        assertThat(centroid.norm()).isCloseTo(2.0, Offset.offset(1e-15))
    }

    @Test
    fun edgeTrueCentroidGreatCircles() {
        // Construct random great circles and divide them randomly into segments.
        // Then make sure that the centroid is approximately at the center of the
        // sphere.  Note that because of the way the centroid is computed, it does
        // not matter how we split the great circle into segments.
        //
        // Note that this is a direct test of the properties that the centroid
        // should have, rather than a test that it matches a particular formula.

        for (iter in 0 until 100) {
            var centroid = S2Point()
            // Choose a coordinate frame for the great circle.
            val (x, y, _) = randomFrame()

            var v0 = x
            var theta = 0.0
            while (theta < 2 * M_PI) {
                val v1 = cos(theta) * x + sin(theta) * y
                centroid = centroid + S2Centroids.trueCentroid(v0, v1)
                v0 = v1;
                theta += randomDouble().pow(10.0)
            }
            // Close the circle.
            centroid = centroid + S2Centroids.trueCentroid(v0, x)
            assertTrue(centroid.norm() <= 2e-14)
        }
    }

}
