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

package dilivia.s2.builder.layers

import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Polyline
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2PolylineLayerUnitTest {

    private val logger = KotlinLogging.logger { }

    fun testS2Polyline(
        input_strs: List<String>,
        expected_str: String,
        edge_type: EdgeType,
        options: S2Builder.Options = S2Builder.Options()
    ) {
        logger.trace { edge_type }
        val builder = S2Builder(options)
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output, options = S2PolylineLayer.Options(edge_type)))
        for (input_str in input_strs) {
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output)).isEqualTo(expected_str)
    }

    // Convenience function that tests both directed and undirected edges.
    fun testS2Polyline(
        input_strs: List<String>,
        expected_str: String,
        options: S2Builder.Options = S2Builder.Options()
    ) {
        testS2Polyline(input_strs, expected_str, EdgeType.DIRECTED, options)
        testS2Polyline(input_strs, expected_str, EdgeType.UNDIRECTED, options)
    }

    fun testS2PolylineUnchanged(input_str: String) {
        testS2Polyline(listOf(input_str), input_str)
    }

    @Test
    fun noEdges() {
        testS2Polyline(emptyList(), "")
    }

    @Test
    fun oneEdge() {
        // Even with undirected edges, S2PolylineLayer prefers to reconstruct edges
        // in their original direction.
        testS2PolylineUnchanged("3:4, 1:1")
        testS2PolylineUnchanged("1:1, 3:4")
    }

    @Test
    fun straightLineWithBacktracking() {
        testS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0")
    }

    @Test
    fun earlyWalkTerminationWithEndLoop1() {
        // Test that the "early walk termination" code (which is needed by
        // S2PolylineVectorLayer in order to implement idempotency) does not create
        // two polylines when it is possible to assemble the edges into one.
        //
        // This example tests a code path where the early walk termination code
        // should not be triggered at all (but was at one point due to a bug).
        val options = S2Builder.Options()
        options.snapFunction = IntLatLngSnapFunction(2)
        testS2Polyline(listOf("0:0, 0:2, 0:1"), "0:0, 0:1, 0:2, 0:1", options)
    }

    @Test
    fun earlyWalkTerminationWithEndLoop2() {
        // This tests a different code path where the walk is terminated early
        // (yield a polyline with one edge), and then the walk is "maximimzed" by
        // appending a two-edge loop to the end.
        testS2Polyline(
            listOf("0:0, 0:1", "0:2, 0:1", "0:1, 0:2"),
            "0:0, 0:1, 0:2, 0:1"
        )
    }

    @Test
    fun simpleLoop() {
        testS2PolylineUnchanged("0:0, 0:5, 5:5, 5:0, 0:0")
    }

    @Test
    fun manyLoops() {
        // This polyline consists of many overlapping loops that keep returning to
        // the same starting vertex (2:2).  This tests whether the implementation is
        // able to assemble the polyline in the original order.
        testS2PolylineUnchanged(
            "0:0, 2:2, 2:4, 2:2, 2:4, 4:4, 4:2, 2:2, 4:4, 4:2, 2:2, 2:0, 2:2, " +
                    "2:0, 4:0, 2:2, 4:2, 2:2, 0:2, 0:4, 2:2, 0:4, 0:2, 2:2, 0:4, 2:2, " +
                    "0:2, 2:2, 0:0, 0:2, 2:2, 0:0"
        )
    }

    @Test
    fun unorderedLoops() {
        // This test consists of 5 squares that touch diagonally, similar to the 5
        // white squares of a 3x3 chessboard.  The edges of these squares need to be
        // reordered to assemble them into a single unbroken polyline.
        testS2Polyline(
            listOf(
                "3:3, 3:2, 2:2, 2:3, 3:3",
                "1:0, 0:0, 0:1, 1:1, 1:0",
                "3:1, 3:0, 2:0, 2:1, 3:1",
                "1:3, 1:2, 0:2, 0:1, 1:3",
                "1:1, 1:2, 2:2, 2:1, 1:1",  // Central square
            ),
            "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, 0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3"
        )
    }

    @Test
    fun splitEdges() {
        // Test reconstruction of a polyline where two edges have been split into
        // many pieces by crossing edges.  This example is particularly challenging
        // because (1) the edges form a loop, and (2) the first and last edges are
        // identical (but reversed).  This is designed to test the heuristics that
        // attempt to find the first edge of the input polyline.
        val options = S2Builder.Options()
        options.splitCrossingEdges = true
        options.snapFunction = IntLatLngSnapFunction(7)
        testS2Polyline(
            listOf("0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10"),
            "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, -1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
            options
        )
    }

    @Test
    fun simpleEdgeLabels() {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polyline()
        val labelSetIds: S2PointVectorLabelSetIds = mutableListOf()
        val labelSetLexicon = IdSetLexicon()
        builder.startLayer(
            S2PolylineLayer(
                output,
                labelSetIds,
                labelSetLexicon,
                S2PolylineLayer.Options(EdgeType.UNDIRECTED)
            )
        )
        builder.setLabel(5)
        builder.addPolyline(makePolyline("0:0, 0:1, 0:2"))
        builder.pushLabel(7)
        builder.addPolyline(makePolyline("0:3, 0:2"))
        builder.clearLabels()
        builder.addPolyline(makePolyline("0:3, 0:4, 0:5"))
        builder.setLabel(11)
        builder.addPolyline(makePolyline("0:6, 0:5"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = listOf(listOf(5), listOf(5), listOf(5, 7), emptyList(), emptyList(), listOf(11))
        assertThat(labelSetIds.size).isEqualTo(expected.size)
        for (i in expected.indices) {
            assertThat(labelSetLexicon.idSet(labelSetIds[i]).size()).isEqualTo(expected[i].size)
            var j = 0
            for (label in labelSetLexicon.idSet(labelSetIds[i])) {
                assertThat(label).isEqualTo(expected[i][j++])
            }
        }
    }

    @Test
    fun invalidPolyline() {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polyline()
        val options = S2PolylineLayer.Options()
        options.setValidate(true)
        builder.startLayer(S2PolylineLayer(output, options = options))
        val vertices = mutableListOf<S2Point>()
        vertices.add(S2Point(1, 0, 0))
        vertices.add(S2Point(-1, 0, 0))
        val input = S2Polyline(vertices, S2Debug.DISABLE)
        builder.addPolyline(input)
        val error = S2Error()
        assertThat(builder.build(error)).isFalse()
        assertThat(error.code).isEqualTo(S2Error.ANTIPODAL_VERTICES)
    }

    @Test
    fun indexedS2PolylineLayerAddsShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PolylineLayer(index))
        val polylineStr = "0:0, 0:10"
        builder.addPolyline(makePolyline(polylineStr))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(1)
        val polyline = (index.shape(0) as S2Polyline.Shape).polyline
        assertThat(S2TextParser.toString(polyline)).isEqualTo(polylineStr)
    }

    @Test
    fun indexedS2PolylineLayerAddsEmptyShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PolylineLayer(index))
        val line = S2Polyline()
        builder.addPolyline(line)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(0)
    }

}
