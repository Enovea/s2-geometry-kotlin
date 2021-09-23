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
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2Debug
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.pow
import kotlin.math.sin
import kotlin.math.tan

class S2PolylineUnitTest {

    private val logger = KotlinLogging.logger {}
    
    @Test
    fun basic() {
        val vertices = emptyList<S2Point>().toMutableList()
        var empty = S2Polyline(vertices)
        assertThat(empty.rectBound).isEqualTo(S2LatLngRect.empty())
        empty = empty.reversed()
        assertThat(empty.numVertices).isEqualTo(0)

        val latlngs = listOf(
                S2LatLng.fromDegrees(0, 0),
                S2LatLng.fromDegrees(0, 90),
                S2LatLng.fromDegrees(0, 180)
        )
        var semiEquator = S2Polyline.fromLatLng(latlngs)
        assertThat(S2PointUtil.approxEquals(semiEquator.interpolate(0.5), S2Point(0, 1, 0))).isTrue()
        semiEquator = semiEquator.reversed()
        assertThat(semiEquator.vertex(2)).isEqualTo(S2Point(1, 0, 0))
    }

    @Test
    fun getLengthAndCentroid() {
        // Construct random great circles and divide them randomly into segments.
        // Then make sure that the length and centroid are correct.  Note that
        // because of the way the centroid is computed, it does not matter how
        // we split the great circle into segments.

        repeat(100) {
            // Choose a coordinate frame for the great circle.
            val (x, y, _) = S2Random.randomFrame()

            val vertices = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                val p = cos(theta) * x + sin(theta) * y
                if (vertices.isEmpty() || p != vertices.last())
                    vertices.add(p)
                theta += S2Random.randomDouble().pow(10.0)
            }
            // Close the circle.
            vertices.add(vertices[0])
            val line = S2Polyline(vertices)
            val length = line.length
            assertThat(abs(length.radians - 2 * M_PI) <= 2e-14).isTrue()
            val centroid = line.centroid
            assertThat(centroid.norm() <= 2e-14).isTrue()
        }
    }

    @Test
    fun mayIntersect() {
        val vertices = mutableListOf(
                S2Point(1.0, -1.1, 0.8).normalize(),
                S2Point(1.0, -0.8, 1.1).normalize()
        )
        val line = S2Polyline(vertices)
        for (face in 0..5) {
            val cell = S2Cell.fromFace(face)
            assertThat(line.mayIntersect(cell)).isEqualTo((face and 1) == 0)
        }
    }

    @Test
    fun interpolate() {
        var vertices = mutableListOf(
                S2Point(1, 0, 0),
                S2Point(0, 1, 0),
                S2Point(0, 1, 1).normalize(),
                S2Point(0, 0, 1)
        )
        val line = S2Polyline(vertices)
        assertThat(line.interpolate(-0.1)).isEqualTo(vertices[0])
        assertThat(S2PointUtil.approxEquals(line.interpolate(0.1), S2Point(1.0, tan(0.2 * M_PI / 2), 0.0).normalize())).isTrue()
        assertThat(S2PointUtil.approxEquals(line.interpolate(0.25), S2Point(1, 1, 0).normalize())).isTrue()
        assertThat(line.interpolate(0.5)).isEqualTo(vertices[1])
        assertThat(S2PointUtil.approxEquals(vertices[2], line.interpolate(0.75))).isTrue()
        val (vertex0, next_vertex0) = line.getSuffix(-0.1)
        assertThat(vertex0).isEqualTo(vertices[0])
        assertThat(next_vertex0).isEqualTo(1)
        val (vertex2, next_vertex2) = line.getSuffix(0.75)
        assertThat(S2PointUtil.approxEquals(vertices[2], vertex2)).isTrue()
        assertThat(next_vertex2).isEqualTo(3)
        val (vertex3, next_vertex3) = line.getSuffix(1.1)
        assertThat(vertex3).isEqualTo(vertices[3])
        assertThat(next_vertex3).isEqualTo(4)

        // Check the case where the interpolation fraction is so close to 1 that
        // the interpolated point is identical to the last vertex.
        vertices = mutableListOf(
            S2Point(1, 1, 1).normalize(),
            S2Point(1.0, 1.0, 1 + 1e-15).normalize(),
            S2Point(1.0, 1.0, 1 + 2e-15).normalize()
        )
        val shortLine = S2Polyline(vertices)
        val (vertex, next_vertex) = shortLine.getSuffix(1.0 - 2e-16)
        assertThat(vertex).isEqualTo(vertices[2])
        assertThat(next_vertex).isEqualTo(3)
    }

    @Test
    fun unInterpolate() {
        val vertices = mutableListOf( S2Point(1, 0, 0) )
        val pointLine = S2Polyline(vertices)
        assertThat(pointLine.unInterpolate(S2Point(0, 1, 0), 1)).isZero()

        vertices.add(S2Point(0, 1, 0))
        vertices.add(S2Point(0, 1, 1).normalize())
        vertices.add(S2Point(0, 0, 1))
        val line = S2Polyline(vertices)

        var interpolated: Pair<S2Point, Int>
        interpolated = line.getSuffix(-0.1)
        assertThat(line.unInterpolate(interpolated.first, interpolated.second)).isEqualTo(0.0)
        interpolated = line.getSuffix(0.0)
        assertThat(line.unInterpolate(interpolated.first, interpolated.second)).isEqualTo(0.0)
        interpolated = line.getSuffix(0.5)
        assertThat(line.unInterpolate(interpolated.first, interpolated.second)).isEqualTo(0.5)
        interpolated = line.getSuffix(0.75)
        assertThat(line.unInterpolate(interpolated.first, interpolated.second)).isEqualTo(0.75)
        interpolated = line.getSuffix(1.1)
        assertThat(line.unInterpolate(interpolated.first, interpolated.second)).isEqualTo(1.0)

        // Check that the return value is clamped to 1.0.
        assertThat(line.unInterpolate(S2Point(0, 1, 0), vertices.size)).isEqualTo(1.0)
    }

    @Test
    fun project() {
        val latlngs = listOf(
                S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
                S2LatLng.fromDegrees(0, 2), S2LatLng.fromDegrees(1, 2)
        )
        val line = S2Polyline.fromLatLng(latlngs)


        var projected = line.project(S2LatLng.fromDegrees(0.5, -0.5).toPoint())
        assertThat(S2PointUtil.approxEquals(projected.first, S2LatLng.fromDegrees(0, 0).toPoint())).isTrue()
        assertThat(projected.second).isEqualTo(1)
        projected = line.project(S2LatLng.fromDegrees(0.5, 0.5).toPoint())
        assertThat(S2PointUtil.approxEquals(projected.first, S2LatLng.fromDegrees(0.0, 0.5).toPoint())).isTrue()
        assertThat(projected.second).isEqualTo(1)
        projected = line.project(S2LatLng.fromDegrees(0.5, 1.0).toPoint())
        assertThat(S2PointUtil.approxEquals(projected.first, S2LatLng.fromDegrees(0, 1).toPoint())).isTrue()
        assertThat(projected.second).isEqualTo(2)
        projected = line.project(S2LatLng.fromDegrees(-0.5, 2.5).toPoint())
        assertThat(S2PointUtil.approxEquals(projected.first, S2LatLng.fromDegrees(0, 2).toPoint())).isTrue()
        assertThat(projected.second).isEqualTo(3)
        projected = line.project(S2LatLng.fromDegrees(2, 2).toPoint())
        assertThat(S2PointUtil.approxEquals(projected.first, S2LatLng.fromDegrees(1, 2).toPoint())).isTrue()
        assertThat(projected.second).isEqualTo(4)
    }

    @Test
    fun isOnRight() {
        var latlngs = listOf(
                S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
                S2LatLng.fromDegrees(0, 2), S2LatLng.fromDegrees(1, 2)
        )
        val line = S2Polyline.fromLatLng(latlngs)

        assertThat(line.isOnRight(S2LatLng.fromDegrees(-0.5, 0.5).toPoint())).isTrue()
        assertThat(line.isOnRight(S2LatLng.fromDegrees(0.5, -0.5).toPoint())).isFalse()
        assertThat(line.isOnRight(S2LatLng.fromDegrees(0.5, 0.5).toPoint())).isFalse()
        assertThat(line.isOnRight(S2LatLng.fromDegrees(0.5, 1.0).toPoint())).isFalse()
        assertThat(line.isOnRight(S2LatLng.fromDegrees(-0.5, 2.5).toPoint())).isTrue()
        assertThat(line.isOnRight(S2LatLng.fromDegrees(1.5, 2.5).toPoint())).isTrue()

        // Explicitly test the case where the closest point is an interior vertex.
        latlngs = listOf(
                S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
                S2LatLng.fromDegrees(-1, 0)
        )
        val line2 = S2Polyline.fromLatLng(latlngs)

        // The points are chosen such that they are on different sides of the two
        // edges that the interior vertex is on.
        assertThat(line2.isOnRight(S2LatLng.fromDegrees(-0.5, 5.0).toPoint())).isFalse()
        assertThat(line2.isOnRight(S2LatLng.fromDegrees(5.5, 5.0).toPoint())).isFalse()
    }

    @Test
    fun intersectsEmptyPolyline() {
        val line1 = makePolyline("1:1, 4:4")
        val emptyPolyline = S2Polyline()
        assertThat(emptyPolyline.intersects(line1)).isFalse()
    }

    @Test
    fun intersectsOnePointPolyline() {
        val line1 = makePolyline("1:1, 4:4")
        val line2 = makePolyline("1:1")
        assertThat(line1.intersects(line2)).isFalse()
    }

    @Test
    fun intersects() {
        val line1  = makePolyline("1:1, 4:4")
        val smallCrossing  = makePolyline("1:2, 2:1")
        val smallNoncrossing  = makePolyline("1:2, 2:3")
        val bigCrossing  = makePolyline("1:2, 2:3, 4:3")

        assertThat(line1.intersects(smallCrossing)).isTrue()
        assertThat(line1.intersects(smallNoncrossing)).isFalse()
        assertThat(line1.intersects(bigCrossing)).isTrue()
    }

    @Test
    fun intersectsAtVertex() {
        val line1  = makePolyline("1:1, 4:4, 4:6")
        val line2  = makePolyline("1:1, 1:2")
        val line3  = makePolyline("5:1, 4:4, 2:2")
        assertThat(line1.intersects(line2)).isTrue()
        assertThat(line1.intersects(line3)).isTrue()
    }

    @Test
    fun intersectsVertexOnEdge() {
        val horizontalLeftToRight  = makePolyline("0:1, 0:3")
        val verticalBottomToTop  = makePolyline("-1:2, 0:2, 1:2")
        val horizontalRightToLeft  = makePolyline("0:3, 0:1")
        val verticalTopToBottom  = makePolyline("1:2, 0:2, -1:2")
        assertThat(horizontalLeftToRight.intersects(verticalBottomToTop)).isTrue()
        assertThat(horizontalLeftToRight.intersects(verticalTopToBottom)).isTrue()
        assertThat(horizontalRightToLeft.intersects(verticalBottomToTop)).isTrue()
        assertThat(horizontalRightToLeft.intersects(verticalTopToBottom)).isTrue()
    }

    fun checkSubsample(polyline_str: String, tolerance_degrees: Double, expected_str: String, check: S2Debug = S2Debug.ALLOW) {
        logger.info { """"$polyline_str, tolerance $tolerance_degrees""" }
        val polyline = makePolyline(polyline_str, check)
        val indices = polyline.subsampleVertices(S1Angle.degrees(tolerance_degrees))
        assertThat(indices.joinToString(",")).isEqualTo(expected_str)
    }

    @Test
    fun subsampleVerticesTrivialInputs() {
        // No vertices.
        checkSubsample("", 1.0, "")
        // One vertex.
        checkSubsample("0:1", 1.0, "0")
        // Two vertices.
        checkSubsample("10:10, 11:11", 5.0, "0,1")
        // Three points on a straight line.
        // In theory, zero tolerance should work, but in practice there are floating
        // point errors.
        checkSubsample("-1:0, 0:0, 1:0", 1e-15, "0,2")
        // Zero tolerance on a non-straight line.
        checkSubsample("-1:0, 0:0, 1:1", 0.0, "0,1,2")
        // Negative tolerance should return all vertices.
        checkSubsample("-1:0, 0:0, 1:1", -1.0, "0,1,2")
        // Non-zero tolerance with a straight line.
        checkSubsample("0:1, 0:2, 0:3, 0:4, 0:5", 1.0, "0,4")

        // And finally, verify that we still do something reasonable if the client
        // passes in an invalid polyline with two or more adjacent vertices.
        checkSubsample("0:1, 0:1, 0:1, 0:2", 0.0, "0,3", check = S2Debug.DISABLE)
    }

    @Test
    fun subsampleVerticesSimpleExample() {
        val polyStr = "0:0, 0:1, -1:2, 0:3, 0:4, 1:4, 2:4.5, 3:4, 3.5:4, 4:4"
        checkSubsample(polyStr, 3.0, "0,9")
        checkSubsample(polyStr, 2.0, "0,6,9")
        checkSubsample(polyStr, 0.9, "0,2,6,9")
        checkSubsample(polyStr, 0.4, "0,1,2,3,4,6,9")
        checkSubsample(polyStr, 0.0, "0,1,2,3,4,5,6,7,8,9")
    }

    @Test
    fun subsampleVerticesGuarantees() {
        // Check that duplicate vertices are never generated.
        checkSubsample("10:10, 12:12, 10:10", 5.0, "0")
        checkSubsample("0:0, 1:1, 0:0, 0:120, 0:130", 5.0, "0,3,4")

        // Check that points are not collapsed if they would create a line segment
        // longer than 90 degrees, and also that the code handles original polyline
        // segments longer than 90 degrees.
        checkSubsample("90:0, 50:180, 20:180, -20:180, -50:180, -90:0, 30:0, 90:0", 5.0, "0,2,4,5,6,7")

        // Check that the output polyline is parametrically equivalent and not just
        // geometrically equivalent, i.e. that backtracking is preserved.  The
        // algorithm achieves this by requiring that the points must be encountered
        // in increasing order of distance along each output segment, except for
        // points that are within "tolerance" of the first vertex of each segment.
        checkSubsample("10:10, 10:20, 10:30, 10:15, 10:40", 5.0, "0,2,3,4")
        checkSubsample("10:10, 10:20, 10:30, 10:10, 10:30, 10:40", 5.0, "0,2,3,5")
        checkSubsample("10:10, 12:12, 9:9, 10:20, 10:30", 5.0, "0,4")
    }


    fun testEquals(a_str: String, b_str: String, max_error: S1Angle): Boolean {
        val a = makePolyline(a_str)
        val b = makePolyline(b_str)
        return a.approxEquals(b, max_error)
    }

    @Test
    fun approxEquals() {
        val degree = S1Angle.degrees(1)

        // Close lines, differences within max_error.
        assertThat(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.5 * degree)).isTrue()

        // Close lines, differences outside max_error.
        assertThat(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.01 * degree)).isFalse()

        // Same line, but different number of vertices.
        assertThat(testEquals("0:0, 0:10, 0:20", "0:0, 0:20", 0.1 * degree)).isFalse()

        // Same vertices, in different order.
        assertThat(testEquals("0:0, 5:5, 0:10", "5:5, 0:10, 0:0", 0.1 * degree)).isFalse()
    }

    @Test
    fun polylineShapeBasic() {
        val polyline  = makePolyline("0:0, 1:0, 1:1, 2:1")
        val shape = S2Polyline.Shape(polyline = polyline)
        assertThat(shape.polyline).isEqualTo(polyline)
        assertThat(shape.numEdges).isEqualTo(3)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.chain(0).start).isEqualTo(0)
        assertThat(shape.chain(0).length).isEqualTo(3)
        val edge2 = shape . edge (2)
        assertThat(edge2.v0).isEqualTo(S2LatLng.fromDegrees(1, 1).toPoint())
        assertThat(edge2.v1).isEqualTo(S2LatLng.fromDegrees(2, 1).toPoint())
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun polylineShapeEmptyPolyline() {
        val polyline = S2Polyline()
        val shape = S2Polyline.Shape(polyline = polyline)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    fun testNearlyCovers(a_str: String, b_str: String, max_error_degrees: Double, expect_b_covers_a: Boolean, expect_a_covers_b: Boolean, check: S2Debug = S2Debug.ALLOW) {
        logger.info { "a=\"$a_str\", b=\"$b_str\", max error=$max_error_degrees" }
        val a = makePolyline(a_str, check)
        val b = makePolyline(b_str, check)
        val maxError = S1Angle.degrees(max_error_degrees)
        assertThat(b.nearlyCovers(a, maxError)).isEqualTo(expect_b_covers_a)
        assertThat(a.nearlyCovers(b, maxError)).isEqualTo(expect_a_covers_b)
    }

    @Test
    fun polylineCoveringTestPolylineOverlapsSelf() {
        val pline = "1:1, 2:2, -1:10"
        testNearlyCovers(pline, pline, 1e-10, true, true)
    }

    @Test
    fun polylineCoveringTestPolylineDoesNotOverlapReverse() {
        testNearlyCovers("1:1, 2:2, -1:10", "-1:10, 2:2, 1:1", 1e-10, false, false)
    }

    @Test
    fun polylineCoveringTestPolylineOverlapsEquivalent() {
        // These two polylines trace the exact same polyline, but the second one uses
        // three points instead of two.
        testNearlyCovers("1:1, 2:1", "1:1, 1.5:1, 2:1", 1e-10, true, true)
    }

    @Test
    fun polylineCoveringTestShortCoveredByLong() {
        // The second polyline is always within 0.001 degrees of the first polyline,
        // but the first polyline is too long to be covered by the second.
        testNearlyCovers(
                "-5:1, 10:1, 10:5, 5:10", "9:1, 9.9995:1, 10.0005:5", 1e-3, false, true)
    }

    @Test
    fun polylineCoveringTestPartialOverlapOnly() {
        // These two polylines partially overlap each other, but neither fully
        // overlaps the other.
        testNearlyCovers("-5:1, 10:1", "0:1, 20:1", 1.0, false, false)
    }

    @Test
    fun polylineCoveringTestShortBacktracking() {
        // Two lines that backtrack a bit (less than 1.5 degrees) on different edges.
        // A simple greedy matching algorithm would fail on this example.
        val t1 = "0:0, 0:2, 0:1, 0:4, 0:5"
        val t2 = "0:0, 0:2, 0:4, 0:3, 0:5"
        testNearlyCovers(t1, t2, 1.5, true, true)
        testNearlyCovers(t1, t2, 0.5, false, false)
    }

    @Test
    fun polylineCoveringTestLongBacktracking() {
        // Two arcs with opposite direction do not overlap if the shorter arc is
        // longer than max_error, but do if the shorter arc is shorter than max-error.
        testNearlyCovers("5:1, -5:1", "1:1, 3:1", 1.0, false, false)
        testNearlyCovers("5:1, -5:1", "1:1, 3:1", 2.5, false, true)
    }

    @Test
    fun polylineCoveringTestIsResilientToDuplicatePoints() {
        // S2Polyines are not generally supposed to contain adjacent, identical
        // points, but it happens in practice.  We also set S2Debug::DISABLE so
        // debug binaries won't abort on such polylines.
        testNearlyCovers("0:1, 0:2, 0:2, 0:3", "0:1, 0:1, 0:1, 0:3", 1e-10, true, true, S2Debug.DISABLE)
    }

    @Test
    fun polylineCoveringTestCanChooseBetweenTwoPotentialStartingPoints() {
        // Can handle two possible starting points, only one of which leads to finding
        // a correct path.  In the first polyline, the edge from 0:1.1 to 0:0 and the
        // edge from 0:0.9 to 0:2 might be lucrative starting states for covering the
        // second polyline, because both edges are with the max_error of 1.5 degrees
        // from 0:10.  However, only the latter is actually effective.
        testNearlyCovers("0:11, 0:0, 0:9, 0:20", "0:10, 0:15", 1.5, false, true)
    }

    @Test
    fun polylineCoveringTestStraightAndWigglyPolylinesCoverEachOther() {
        testNearlyCovers("40:1, 20:1",
                "39.9:0.9, 40:1.1, 30:1.15, 29:0.95, 28:1.1, 27:1.15, 26:1.05, 25:0.85, 24:1.1, 23:0.9, 20:0.99",
                0.2, true, true)
    }

    @Test
    fun polylineCoveringTestMatchStartsAtLastVertex() {
        // The first polyline covers the second, but the matching segment starts at
        // the last vertex of the first polyline.
        testNearlyCovers(
                "0:0, 0:2", "0:2, 0:3", 1.5, false, true)
    }

    @Test
    fun polylineCoveringTestMatchStartsAtDuplicatedLastVertex() {
        testNearlyCovers(
                "0:0, 0:2, 0:2, 0:2", "0:2, 0:3", 1.5, false, true, S2Debug.DISABLE)
    }

    @Test
    fun polylineCoveringTestEmptyPolylines() {
        // We expect:
        //    anything.covers(empty) = true
        //    empty.covers(nonempty) = false
        testNearlyCovers("0:1, 0:2", "", 0.0, false, true)
        testNearlyCovers("", "", 0.0, true, true)
    }

}
