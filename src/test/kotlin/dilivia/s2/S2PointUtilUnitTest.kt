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
import dilivia.math.matrix.Matrix3x3Double
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.region.S2Cell
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.IEEErem
import kotlin.math.acos
import kotlin.math.asin
import kotlin.math.pow

internal class S2PointUtilUnitTest {

    private val logger = KotlinLogging.logger {  }

    @Test
    fun testFrames() {
        val m = Matrix3x3Double()
        val z = S2Point(0.2, 0.5, -3.3).normalize()
        S2PointUtil.getFrame(z, m)
        assertThat(S2PointUtil.approxEquals(m.col(2), z)).isTrue()
        assertThat(S2PointUtil.isUnitLength(m.col(0))).isTrue()
        assertThat(S2PointUtil.isUnitLength(m.col(1))).isTrue()
        assertThat(m.det()).isEqualTo(1.0, Offset.offset(1e-15))

        assertThat(S2PointUtil.approxEquals(S2PointUtil.toFrame(m, m.col(0)), S2Point(1, 0, 0))).isTrue()
        assertThat(S2PointUtil.approxEquals(S2PointUtil.toFrame(m, m.col(1)), S2Point(0, 1, 0))).isTrue()
        assertThat(S2PointUtil.approxEquals(S2PointUtil.toFrame(m, m.col(2)), S2Point(0, 0, 1))).isTrue()

        assertThat(S2PointUtil.approxEquals(S2PointUtil.fromFrame(m, S2Point(1, 0, 0)), m.col(0))).isTrue()
        assertThat(S2PointUtil.approxEquals(S2PointUtil.fromFrame(m, S2Point(0, 1, 0)), m.col(1))).isTrue()
        assertThat(S2PointUtil.approxEquals(S2PointUtil.fromFrame(m, S2Point(0, 0, 1)), m.col(2))).isTrue()
    }

    fun testRotate(p: S2Point, axis: S2Point, angle: S1Angle) {
        val result = S2PointUtil.rotate(p, axis, angle)

        // "result" should be unit length.
        assertThat(S2PointUtil.isUnitLength(result)).isTrue()

        // "result" and "p" should be the same distance from "axis".
        val kMaxPositionError = 1e-15
        assertThat((S1Angle(result, axis) - S1Angle(p, axis)).abs().radians).isLessThanOrEqualTo(kMaxPositionError)

        // Check that the rotation angle is correct.  We allow a fixed error in the
        // *position* of the result, so we need to convert this into a rotation
        // angle.  The allowable error can be very large as "p" approaches "axis".
        val axisDistance = p.crossProd(axis).norm()
        val maxRotationError: Double = if (axisDistance < kMaxPositionError) {
            2 * M_PI
        } else {
            asin(kMaxPositionError / axisDistance)
        }
        val actualRotation = S2Measures.turnAngle(p, axis, result) + M_PI
        val rotationError = (angle.radians - actualRotation).IEEErem(2 * M_PI)
        assertThat(rotationError).isLessThanOrEqualTo(maxRotationError)
    }

    @Test
    fun testRotate() {
        repeat(1000) {
            val axis = S2Random.randomPoint()
            val target = S2Random.randomPoint()
            // Choose a distance whose logarithm is uniformly distributed.
            var distance = M_PI * 1e-15.pow(S2Random.randomDouble())
            // Sometimes choose points near the far side of the axis.
            if (S2Random.oneIn(5)) distance = M_PI - distance
            val p = S2EdgeDistances.interpolateAtDistance(S1Angle.radians(distance), axis, target)
            // Choose the rotation angle.
            var angle = 2 * M_PI * 1e-15.pow(S2Random.randomDouble())
            if (S2Random.oneIn(3)) angle = -angle
            if (S2Random.oneIn(10)) angle = 0.0
            testRotate(p, axis, S1Angle.radians(angle))
        }
    }

    // Given a point P, return the minimum level at which an edge of some S2Cell
    // parent of P is nearly collinear with S2::Origin().  This is the minimum
    // level for which Sign() may need to resort to expensive calculations in
    // order to determine which side of an edge the origin lies on.
    fun getMinExpensiveLevel(p: S2Point): Int {
        val id = S2CellId.fromPoint(p)
        for (level in 0..S2CellId.kMaxLevel) {
            val cell = S2Cell(id.parent(level))
            repeat(4) { k ->
                val a = cell.getVertex(k)
                val b = cell.getVertex(k + 1)
                if (S2Predicates.triageSign(a, b, S2PointUtil.origin(), a.crossProd(b)) == 0) {
                    return level
                }
            }
        }
        return S2CellId.kMaxLevel + 1
    }

    @Test
    fun testOrigin() {
        // To minimize the number of expensive Sign() calculations,
        // S2::Origin() should not be nearly collinear with any commonly used edges.
        // Two important categories of such edges are:
        //
        //  - edges along a line of longitude (reasonably common geographically)
        //  - S2Cell edges (used extensively when computing S2Cell coverings)
        //
        // This implies that the origin:
        //
        //  - should not be too close to either pole (since all lines of longitude
        //    converge at the poles)
        //  - should not be colinear with edges of any S2Cell except for very small
        //    ones (which are used less frequently)
        //
        // The point chosen below is about 66km from the north pole towards the East
        // Siberian Sea.  The purpose of the STtoUV(2/3) calculation is to keep the
        // origin as far away as possible from the longitudinal edges of large
        // S2Cells.  (The line of longitude through the chosen point is always 1/3
        // or 2/3 of the way across any S2Cell with longitudinal edges that it
        // passes through.)

        assertThat(S2Point(-0.01, 0.01 * S2Coords.stToUv(2.0/ 3.0 ), 1.0).normalize()).isEqualTo(S2PointUtil.origin())

        // Check that the origin is not too close to either pole.  (We don't use
        // S2Earth because we don't want to depend on that package.)
        val distanceKm = acos (S2PointUtil.origin().z) * S2Earth.radiusKm
        assertThat(distanceKm).isGreaterThanOrEqualTo(50.0)
        logger.info { "\nS2PointUtil.origin() coordinates: ${ S2LatLng.fromPoint(S2PointUtil.origin()) }, distance from pole: $distanceKm  km" }

        // Check that S2::Origin() is not collinear with the edges of any large
        // S2Cell.  We do this is two parts.  For S2Cells that belong to either
        // polar face, we simply need to check that S2::Origin() is not nearly
        // collinear with any edge of any cell that contains it (except for small
        // cells < 3 meters across).
        assertThat(getMinExpensiveLevel(S2PointUtil.origin())).isGreaterThanOrEqualTo(22)

        // For S2Cells that belong to the four non-polar faces, only longitudinal
        // edges can possibly be colinear with S2::Origin().  We check these edges
        // by projecting S2::Origin() onto the equator, and then testing all S2Cells
        // that contain this point to make sure that none of their edges are nearly
        // colinear with S2::Origin() (except for small cells < 3 meters across).
        val equatorPoint = S2Point(S2PointUtil.origin().x, S2PointUtil.origin().y, 0.0)
        assertThat(getMinExpensiveLevel(equatorPoint)).isGreaterThanOrEqualTo(22)
    }

}
