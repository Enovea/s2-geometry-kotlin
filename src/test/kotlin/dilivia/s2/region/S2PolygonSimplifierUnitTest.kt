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
@file:Suppress("UsePropertyAccessSyntax")

package dilivia.s2.region

import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S2TextParser
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.edge.S2EdgeDistances
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.cos

fun makeRegularPolygon(center: String, num_points: Int, radius_in_degrees: Double): S2Polygon {
    val radius = S1Angle.degrees(radius_in_degrees)
    return S2Polygon(S2Loop.makeRegularLoop(S2TextParser.makePoint(center), radius, num_points))
}

class S2PolygonSimplifierUnitTest {

    private lateinit var simplified: S2Polygon
    private lateinit var original: S2Polygon

    protected fun setInput(poly: S2Polygon, tolerance_in_degrees: Double) {
        original = poly

        simplified = S2Polygon()
        simplified.initToSimplified(
            original, IdentitySnapFunction(S1Angle.degrees(tolerance_in_degrees))
        )
    }

    protected fun setInput(poly: String, tolerance_in_degrees: Double) {
        setInput(S2TextParser.makePolygon(poly), tolerance_in_degrees)
    }

    // Returns the diameter of a loop (maximum distance between any two
    // points in the loop).
    fun loopDiameter(loop: S2Loop): S1Angle {
        var diameter = S1Angle()
                for (i in 0 until loop.numVertices) {
            val test_point = loop.vertex(i)
            for (j in (i + 1) until loop.numVertices) {
            diameter = maxOf(diameter, S2EdgeDistances.getDistance(test_point, loop.vertex(j), loop.vertex(j + 1)))
        }
        }
        return diameter
    }

    // Returns the maximum distance from any vertex of poly_a to poly_b, that is,
    // the directed Haussdorf distance of the set of vertices of poly_a to the
    // boundary of poly_b.
    //
    // Doesn't consider loops from poly_a that have diameter less than min_diameter
    // in degrees.
    fun maximumDistanceInDegrees(polyA: S2Polygon, polyB: S2Polygon, minDiameterInDegrees: Double): Double {
        var minDistance = 360.0
        var hasBigLoops = false
        for (l in 0 until polyA.numLoops()) {
            val aLoop = polyA.loop(l)
            if (loopDiameter(aLoop).degrees() <= minDiameterInDegrees) {
                continue
            }
            hasBigLoops = true
            for (v in 0 until aLoop.numVertices) {
            val distance = polyB.getDistance(aLoop.vertex(v)).degrees()
            if (distance < minDistance) {
                minDistance = distance
            }
        }
        }
        if (hasBigLoops) {
            return minDistance
        } else {
            return 0.0  // As if the first polygon were empty.
        }
    }

    @Test
    fun noSimplification() {
        setInput("0:0, 0:20, 20:20, 20:0", 1.0)
        assertThat(simplified.numVertices()).isEqualTo(4)

        assertThat(maximumDistanceInDegrees(simplified, original, 0.0)).isEqualTo(0.0)
        assertThat(maximumDistanceInDegrees(original, simplified, 0.0)).isEqualTo(0.0)
    }

    // Here, 10:-2 will be removed and  0:0-20:0 will intersect two edges.
// (The resulting polygon will in fact probably have more edges.)
    @Test
    fun simplifiedLoopSelfIntersects() {
        setInput("0:0, 0:20, 10:-0.1, 20:20, 20:0, 10:-0.2", 0.22)

        // The simplified polygon has the same number of vertices but it should now
        // consists of two loops rather than one.
        assertThat(simplified.numLoops()).isEqualTo(2)
        assertThat(maximumDistanceInDegrees(simplified, original, 0.0)).isLessThanOrEqualTo(0.22)
        assertThat(maximumDistanceInDegrees(original, simplified, 0.22)).isLessThanOrEqualTo(0.22)
    }

    @Test
    fun noSimplificationManyLoops() {
        setInput(
            "0:0,    0:1,   1:0;   0:20, 0:21, 1:20; " +
            "20:20, 20:21, 21:20; 20:0, 20:1, 21:0", 0.01
        )
        assertThat(maximumDistanceInDegrees(simplified, original, 0.0)).isEqualTo(0.0)
        assertThat(maximumDistanceInDegrees(original, simplified, 0.0)).isEqualTo(0.0)
    }

    @Test
    fun tinyLoopDisappears() {
        setInput("0:0, 0:1, 1:1, 1:0", 1.1)
        assertThat(simplified.isEmpty()).isTrue()
    }

    @Test
    fun straightLinesAreSimplified() {
        setInput(
            "0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0," +
            "6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1", 0.01
        )
        assertThat(simplified.numVertices()).isEqualTo(4)
    }

    @Test
    fun edgeSplitInManyPieces() {
        // near_square's right four-point side will be simplified to a vertical
        // line at lng=7.9, that will cut the 9 teeth of the saw (the edge will
        // therefore be broken into 19 pieces).
        val saw: String = "1:1, 1:8, 2:2, 2:8, 3:2, 3:8, 4:2, 4:8, 5:2, 5:8," +
        "6:2, 6:8, 7:2, 7:8, 8:2, 8:8, 9:2, 9:8, 10:1"
        val near_square = "0:0, 0:7.9, 1:8.1, 10:8.1, 11:7.9, 11:0"
        setInput("$saw;$near_square", 0.21)

        assertThat(simplified.isValid()).isTrue()
        assertThat(maximumDistanceInDegrees(simplified, original, 0.0)).isLessThanOrEqualTo(0.11)
        assertThat(maximumDistanceInDegrees(original, simplified, 0.0)).isLessThanOrEqualTo(0.11)
        // The resulting polygon's 9 little teeth are very small and disappear
        // due to the vertex_merge_radius of the polygon builder.  There remains
        // nine loops.
        assertThat(simplified.numLoops()).isEqualTo(9)
    }

    @Test
    fun edgesOverlap() {
        // Two loops, One edge of the second one ([0:1 - 0:2]) is part of an
        // edge of the first one..
        setInput("0:0, 0:3, 1:0; 0:1, -1:1, 0:2", 0.01)
        val true_poly = S2TextParser.makePolygon("0:3, 1:0, 0:0, 0:1, -1:1, 0:2")
        assertThat(simplified.boundaryApproxEquals(true_poly, S1Angle.radians(1e-15))).isTrue()
    }

    // Tests that a regular polygon with many points gets simplified
    // enough.
    @Test fun largeRegularPolygon() {
        val kRadius = 2.0  // in degrees
        val num_initial_points = 1000
        val num_desired_points = 250
        val tolerance = 1.05 * kRadius * (1 - cos(M_PI / num_desired_points))

        setInput(makeRegularPolygon("0:0", num_initial_points, kRadius), tolerance)

        assertThat(maximumDistanceInDegrees(simplified, original, 0.0)).isLessThanOrEqualTo(tolerance)
        assertThat(maximumDistanceInDegrees(original, simplified, 0.0)).isLessThanOrEqualTo(tolerance)
        assertThat(simplified.numVertices()).isBetween(200, 250)
    }

}
