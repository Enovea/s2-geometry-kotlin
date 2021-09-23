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
package dilivia.s2.index.shape

import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.index.CrossingType
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.shape.S2EdgeVectorShape
import dilivia.s2.shape.ShapeEdge
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import org.assertj.core.api.Assertions
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

// A set of edge pairs within an S2ShapeIndex.
typealias EdgePairVector = Set<Pair<ShapeEdgeId, ShapeEdgeId>>

class S2CrossingEdgePairsScannerTest {

    fun getCrossings(index: S2ShapeIndex, type: CrossingType): EdgePairVector {
        val edge_pairs = HashSet<Pair<ShapeEdgeId, ShapeEdgeId>>()
        S2CrossingEdgePairsScanner.visitCrossingEdgePairs(
            index, type, object : EdgePairVisitor {
                override fun visit(a: ShapeEdge, b: ShapeEdge, isInterior: Boolean): Boolean {
                    edge_pairs.add(Pair(a.id, b.id))
                    return true  // Continue visiting.
                }
            })
        return edge_pairs
    }

    fun getCrossingEdgePairsBruteForce(index: S2ShapeIndex, type: CrossingType): EdgePairVector {
        val result = HashSet<Pair<ShapeEdgeId, ShapeEdgeId>>()
        val min_sign = if (type == CrossingType.ALL) 0 else 1
        val a_iter = S2ShapeIndex.EdgeIterator(index)
        while (!a_iter.done()) {
            val a = a_iter.edge()
            val b_iter = S2ShapeIndex.EdgeIterator(a_iter)
            b_iter.next()
            while (!b_iter.done()) {
                val b = b_iter.edge()
                if (S2EdgeCrossings.crossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
                    result.add(Pair(a_iter.shapeEdgeId(), b_iter.shapeEdgeId()))
                }
                b_iter.next()
            }
            a_iter.next()
        }
        return result
    }

    fun testGetCrossingEdgePairs(index: S2ShapeIndex, type: CrossingType) {
        val expected = getCrossingEdgePairsBruteForce(index, type)
        val actual = getCrossings(index, type)
        if (actual != expected) {
            var message = """
            |Unexpected edge pairs; see details below.
            |Expected number of edge pairs: ${expected.size}
            |Actual number of edge pairs: ${actual.size}
            |
          """.trimMargin()
            expected.filter { edgePair -> !actual.contains(edgePair) }
                .forEach { edgePair -> message += "Missing value: $edgePair\n" }
            actual.filter { edgePair -> !expected.contains(edgePair) }
                .forEach { edgePair -> message += "Extra value: $edgePair\n" }
            Assertions.fail<Void>(message)
        }
    }

    @Test
    fun getCrossingEdgePairsNoIntersections() {
        val index = MutableS2ShapeIndex()
        testGetCrossingEdgePairs(index, CrossingType.ALL)
        testGetCrossingEdgePairs(index, CrossingType.INTERIOR)
    }

    @Test
    fun getCrossingEdgePairsEdgeGrid() {
        val kGridSize = 10  // (kGridSize + 1) * (kGridSize + 1) crossings
        val index = MutableS2ShapeIndex()
        val shape = S2EdgeVectorShape()
        for (i in 0..kGridSize) {
            shape.add(S2LatLng.fromDegrees(0, i).toPoint(), S2LatLng.fromDegrees(kGridSize, i).toPoint())
            shape.add(S2LatLng.fromDegrees(i, 0).toPoint(), S2LatLng.fromDegrees(i, kGridSize).toPoint())
        }
        index.add(shape)
        testGetCrossingEdgePairs(index, CrossingType.ALL)
        testGetCrossingEdgePairs(index, CrossingType.INTERIOR)
    }

    // This function recursively verifies that HasCrossing returns the given
    // result for all possible cyclic permutations of the loop vertices for the
    // given set of loops.
    fun testHasCrossingPermutations(loops: MutableList<S2Loop>, i: Int, has_crossing: Boolean) {
        if (i == loops.size) {
            val index = MutableS2ShapeIndex()
            val polygon = S2Polygon(loops, debugOverride = S2Debug.DISABLE)
            index.add(S2Polygon.Shape(polygon = polygon))
            assertThat(has_crossing).isEqualTo(hasSelfIntersection(index))
        } else {
            val orig_loop = loops[i]
            for (j in 0 until orig_loop.numVertices) {
                val vertices = mutableListOf<S2Point>()
                for (k in 0 until orig_loop.numVertices) {
                    vertices.add(orig_loop.vertex(j + k))
                }
                loops[i] = S2Loop(vertices, debugOverride = S2Debug.DISABLE)
                testHasCrossingPermutations(loops, i + 1, has_crossing)
            }
            loops[i] = orig_loop
        }
    }

    // Given a string reprsenting a polygon, and a boolean indicating whether this
    // polygon has any self-intersections or loop crossings, verify that all
    // HasSelfIntersection returns the expected result for all possible cyclic
    // permutations of the loop vertices.
    fun testHasCrossing(polygon_str: String, has_crossing: Boolean) {
        // Set S2Debug::DISABLE to allow invalid polygons.
        val polygon = S2TextParser.makePolygon(polygon_str, S2Debug.DISABLE)
        val loops = polygon.loops()
        testHasCrossingPermutations(loops.toMutableList(), 0, has_crossing)
    }

    @Test
    fun findSelfIntersectionBasic() {
        // Coordinates are (lat,lng), which can be visualized as (y,x).
        testHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false)
        testHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true)  // duplicate vertex
        testHasCrossing("0:0, 0:1, 1:0, 1:1", true)  // edge crossing
        testHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true)  // duplicate edge
        testHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true)  // reversed edge
        testHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0", true)  // vertex crossing
    }

    companion object {

        private val logger = KotlinLogging.logger { }

        // Return true if any loop crosses any other loop (including vertex crossings
        // and duplicate edges), or any loop has a self-intersection (including
        // duplicate vertices).
        fun hasSelfIntersection(index: MutableS2ShapeIndex): Boolean {
            val error = S2Error()
            if (S2CrossingEdgePairsScanner.findSelfIntersection(index, error)) {
                logger.error { error }
                return true
            }
            return false
        }

    }
}
