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

import dilivia.PreConditions
import dilivia.collections.resizeWith
import dilivia.collections.swap
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Factory
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import org.assertj.core.api.Assertions
import org.junit.jupiter.api.AfterEach
import org.junit.jupiter.api.Test
import kotlin.math.pow

class S2PolygonIsValidUnitTest() {
    private var init_oriented = false
    private var modifyPolygonHook: (S2Polygon) -> Unit = {}
    private val vloops: MutableList<MutableList<S2Point>> = mutableListOf()

    init {
        S2Random.reset(1)
    }

    fun addLoop(): MutableList<S2Point> {
        vloops.add(mutableListOf<S2Point>())
        return vloops.last()
    }

    // Create "num_loops" nested regular loops around a common center point.
    // All loops have the same number of vertices (at least "min_vertices").
    // Furthermore, the vertices at the same index position are collinear with
    // the common center point of all the loops.  The loop radii decrease
    // exponentially in order to prevent accidental loop crossings when one of
    // the loops is modified.
    fun addConcentricLoops(num_loops: Int, min_vertices: Int) {
        PreConditions.requireLE(num_loops, 10)  // Because radii decrease exponentially.
        val center = S2Random.randomPoint()
        val num_vertices = min_vertices + S2Random.randomInt(10)
        for (i in 0 until num_loops) {
            val radius = S1Angle.degrees(80 * Math.pow(0.1, i.toDouble()))
            addLoop().addAll(S2Factory.makeRegularPoints(center, radius, num_vertices))
        }
    }

    @AfterEach
    fun reset() {
        vloops.clear()
    }

    fun checkInvalid(snippet: String) {
        val loops = mutableListOf<S2Loop>()
        for (vloop in vloops) {
            loops.add(S2Loop(vloop, debugOverride = S2Debug.DISABLE))
        }

        loops.shuffle(S2Random.random)
        val polygon = S2Polygon()
        polygon.debugOverride = S2Debug.DISABLE
        if (init_oriented) {
            polygon.initOriented(loops)
        } else {
            polygon.initNested(loops)
        }
        modifyPolygonHook.invoke(polygon)

        val error = S2Error()
        Assertions.assertThat(polygon.findValidationError(error)).isTrue()
        Assertions.assertThat(error.text)
            .withFailMessage("\nActual error: $error\nExpected substring: $snippet")
            .contains(snippet)
        reset()
    }

    @Test
    fun unitLength() {
        // This test can only be run in optimized builds because there are
        // S2_DCHECK(IsUnitLength()) calls scattered throughout the S2 code.
        if (PreConditions.enabled) return
        repeat(kIters) {
            addConcentricLoops(1 + S2Random.randomInt(6), 3 /*min_vertices*/)
            val vloop = vloops[S2Random.randomInt(vloops.size)]
            val p = vloop[S2Random.randomInt(vloop.size)]
            when (S2Random.randomInt(3)) {
                0 -> p *= 0.0
                1 -> p *= 1e-30 * 1e60.pow(S2Random.randomDouble())
                2 -> p *= Double.NaN
            }
            checkInvalid("unit length")
        }
    }

    @Test
    fun vertexCount() {
        repeat(kIters) {
            val vloop = addLoop()
            if (S2Random.oneIn(2)) {
                vloop.add(S2Random.randomPoint())
                vloop.add(S2Random.randomPoint())
            }
            checkInvalid("at least 3 vertices")
        }
    }

    @Test
    fun duplicateVertex() {
        repeat(kIters) {
            addConcentricLoops(1, 3 /*min_vertices*/)
            val vloop = vloops[0]
            val n = vloop.size
            val i = S2Random.randomInt(n)
            val j = S2Random.randomInt(n - 1)
            vloop[i] = vloop[j + if (j >= i) 1 else 0]
            checkInvalid("duplicate vertex")
        }
    }

    @Test
    fun selfIntersection() {
        repeat(kIters) {
            // Use multiple loops so that we can test both holes and shells.  We need
            // at least 5 vertices so that the modified edges don't intersect any
            // nested loops.
            addConcentricLoops(1 + S2Random.randomInt(6), 5 /*min_vertices*/)
            val vloop = vloops[S2Random.randomInt(vloops.size)]
            val n = vloop.size
            val i = S2Random.randomInt(n)
            vloop.swap(i, (i + 1) % n)
            checkInvalid("crosses edge")
        }
    }

    @Test
    fun emptyLoop() {
        repeat(kIters) {
            addConcentricLoops(S2Random.randomInt(5), 3 /*min_vertices*/)
            addLoop().add(S2Loop.kEmptyVertex())
            checkInvalid("empty loop")
        }
    }

    @Test
    fun fullLoop() {
        repeat(kIters) {
            // This is only an error if there is at least one other loop.
            addConcentricLoops(1 + S2Random.randomInt(5), 3 /*min_vertices*/)
            addLoop().add(S2Loop.kFullVertex())
            checkInvalid("full loop")
        }
    }

    @Test
    fun loopsCrossing() {
        repeat(kIters) {
            addConcentricLoops(2, 4 /*min_vertices*/)
            // Both loops have the same number of vertices, and vertices at the same
            // index position are collinear with the center point, so we can create a
            // crossing by simply exchanging two vertices at the same index position.
            val n = vloops[0].size
            val i = S2Random.randomInt(n)
            val tmp = vloops[0][i]
            vloops[0][i] = vloops[1][i]
            vloops[1][i] = tmp
            if (S2Random.oneIn(2)) {
                // By copy the two adjacent vertices from one loop to the other, we can
                // ensure that the crossings happen at vertices rather than edges.
                vloops[0][(i + 1) % n] = vloops[1][(i + 1) % n]
                vloops[0][(i + n - 1) % n] = vloops[1][(i + n - 1) % n]
            }
            checkInvalid("crosses loop")
        }
    }

    @Test
    fun duplicateEdge() {
        repeat(kIters) {
            addConcentricLoops(2, 4 /*min_vertices*/)
            val n = vloops[0].size
            if (S2Random.oneIn(2)) {
                // Create a shared edge (same direction in both loops).
                val i = S2Random.randomInt(n)
                vloops[0][i] = vloops[1][i]
                vloops[0][(i + 1) % n] = vloops[1][(i + 1) % n]
            } else {
                // Create a reversed edge (opposite direction in each loop) by cutting
                // loop 0 into two halves along one of its diagonals and replacing both
                // loops with the result.
                val split = 2 + S2Random.randomInt(n - 3)
                vloops[1].clear()
                vloops[1].add(vloops[0][0])
                for (i in split until n) {
                    vloops[1].add(vloops[0][i])
                }
                vloops[0].resizeWith(split + 1) { S2Point() }
            }
            checkInvalid("has duplicate")
        }
    }

    @Test
    fun inconsistentOrientations() {
        repeat(kIters) { iter ->
            addConcentricLoops(2 + S2Random.randomInt(5), 3 /*min_vertices*/)
            init_oriented = true
            checkInvalid("Inconsistent loop orientations")
        }
    }

    @Test
    fun loopDepthNegative() {
        modifyPolygonHook = { polygon: S2Polygon ->
            val i = S2Random.randomInt(polygon.numLoops())
            if (i == 0 || S2Random.oneIn(3)) {
                polygon.loop(i).depth = -1
            } else {
                polygon.loop(i).depth = polygon.loop(i - 1).depth + 2
            }
        }
        repeat(kIters) {
            addConcentricLoops(1 + S2Random.randomInt(4), 3 /*min_vertices*/)
            checkInvalid("invalid loop depth")
        }
    }

    @Test
    fun loopNestingInvalid() {
        modifyPolygonHook = { polygon: S2Polygon ->
            val i = S2Random.randomInt(polygon.numLoops())
            polygon.loop(i).invert()
        }
        repeat(kIters) {
            addConcentricLoops(2 + S2Random.randomInt(4), 3 /*min_vertices*/)
            // Randomly invert all the loops in order to generate cases where the
            // outer loop encompasses almost the entire sphere.  This tests different
            // code paths because bounding box checks are not as useful.
            if (S2Random.oneIn(2)) {
                for (loop in vloops) {
                    loop.reverse()
                }
            }
            checkInvalid("Invalid nesting")
        }
    }

    @Test
    fun fuzzTest() {
        // Check that the S2Loop/S2Polygon constructors and IsValid() don't crash
        // when they receive arbitrary invalid input.  (We don't test large inputs
        // it is assumed that the client enforces their own size limits before even
        // attempting to construct geometric objects.)
        if (PreConditions.enabled)
            return;  // Requires unit length vertices.
        repeat(kIters) {
            val num_loops = 1 + S2Random.randomInt(10)
            for (i in 0 until num_loops) {
                val num_vertices = S2Random.randomInt(10)
                val vloop = addLoop()
                while (vloop.size < num_vertices) {
                    // Since the number of vertices is random, we automatically test empty
                    // loops, full loops, and invalid vertex counts.  Also since most
                    // vertices are random, we automatically get self-intersections and
                    // loop crossings.  That leaves zero and NaN vertices, duplicate
                    // vertices, and duplicate edges to be created explicitly.
                    if (S2Random.oneIn(10)) {
                        // Zero vertex.
                        vloop.add(S2Point(0, 0, 0))
                    } else if (S2Random.oneIn(10)) {
                        // NaN vertex.
                        vloop += Double.NaN * S2Point()
                    } else if (S2Random.oneIn(10) && vloop.isNotEmpty()) {
                        // Duplicate vertex.
                        vloop.add(vloop[S2Random.randomInt(vloop.size)])
                    } else if (S2Random.oneIn(10) && vloop.size + 2 <= num_vertices) {
                        // Try to copy an edge from a random loop.
                        val other = vloops[S2Random.randomInt(vloops.size)]
                        val n = other.size
                        if (n >= 2) {
                            var k0 = S2Random.randomInt(n)
                            var k1 = (k0 + 1) % n
                            if (S2Random.oneIn(2)) k0 = k1.also { k1 = k0 } // swap(k0, k1);  // Copy reversed edge.
                            vloop.add(other[k0])
                            vloop.add(other[k1])
                        }
                    } else {
                        // Random non-unit-length point.
                        val p = S2Random.randomPoint()
                        vloop.add(1e-30 * 1e60.pow(S2Random.randomDouble()) * p)
                    }
                }
            }
            checkInvalid("");  // We could get any error message.
        }
    }

    companion object {
        val kIters = 100
    }

}
