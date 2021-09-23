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

import dilivia.s2.S2Error
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2TextParser
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph.PolylineType
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Polyline
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test


class S2PolylineVectorLayerUnitTest {

    private fun testS2PolylineVector(
        inputStrs: List<String>,
        expectedStrs: List<String>,
        edgeType: EdgeType,
        layerOptions: S2PolylineVectorLayer.Options = S2PolylineVectorLayer.Options(),
        builderOptions: S2Builder.Options = S2Builder.Options()
    ) {
        layerOptions.edgeType = edgeType
        val builder = S2Builder(builderOptions)
        val output = mutableListOf<S2Polyline>()
        builder.startLayer(S2PolylineVectorLayer(output, options = layerOptions))
        for (input_str in inputStrs) {
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val outputStrs = mutableListOf<String>()
        for (polyline in output) {
            outputStrs.add(S2TextParser.toString(polyline))
        }
        assertThat(outputStrs.joinToString("; ")).isEqualTo(expectedStrs.joinToString("; "))
    }

    // Convenience function that tests both directed and undirected edges.
    private fun testS2PolylineVector(
        inputStrs: List<String>,
        expectedStrs: List<String>,
        layerOptions: S2PolylineVectorLayer.Options = S2PolylineVectorLayer.Options(),
        builderOptions: S2Builder.Options = S2Builder.Options()
    ) {
        testS2PolylineVector(inputStrs, expectedStrs, EdgeType.DIRECTED, layerOptions, builderOptions)
        testS2PolylineVector(inputStrs, expectedStrs, EdgeType.UNDIRECTED, layerOptions, builderOptions)
    }

    private fun testS2PolylineVectorUnchanged(inputStrs: List<String>) {
        testS2PolylineVector(inputStrs, inputStrs)
    }

    @Test
    fun noEdges() {
        testS2PolylineVectorUnchanged(emptyList())
    }

    @Test
    fun twoPolylines() {
        testS2PolylineVectorUnchanged(listOf("0:0, 1:1, 2:2", "4:4, 3:3"))
    }

    @Test
    fun joiningPolylines() {
        // Check that polylines are joined together when possible, even if they were
        // not adjacent in the input.  For undirected edges, the polyline direction
        // should be chosen such that the first edge of the polyline was added to
        // S2Builder before the last edge of the polyline.
        testS2PolylineVector(
            listOf("1:1, 2:2", "3:3, 2:2", "0:0, 1:1"),
            listOf("3:3, 2:2", "0:0, 1:1, 2:2"),
            EdgeType.DIRECTED
        )
        testS2PolylineVector(
            listOf("1:1, 2:2", "3:3, 2:2", "0:0, 1:1"),
            listOf("3:3, 2:2, 1:1, 0:0"),
            EdgeType.UNDIRECTED
        )
    }

    @Test
    fun segmentNetwork() {
        // Test a complex network of polylines that meet at shared vertices.
        testS2PolylineVectorUnchanged(
            listOf(
                "0:0, 1:1, 2:2",
                "2:2, 2:3, 2:4",
                "2:4, 3:4, 4:4",
                "2:2, 3:2, 4:2",
                "4:2, 4:3, 4:4",
                "1:0, 2:2",
                "0:1, 2:2",
                "5:4, 4:4",
                "4:5, 4:4",
                "2:4, 2:5, 1:5, 1:4, 2:4",
                "4:2, 6:1, 5:0",  // Two nested loops
                "4:2, 7:0, 6:-1",
                "11:1, 11:0, 10:0, 10:1, 11:1"  // Isolated loop
            )
        )
    }

    @Test
    fun multipleIntersectingWalks() {
        // This checks idempotency for directed edges in the case of several
        // polylines that share edges (and that even share loops).  The test
        // happens to pass for undirected edges as well.
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.polylineType = PolylineType.WALK
        val input = listOf(
            "5:5, 5:6, 6:5, 5:5, 5:4, 5:3",
            "4:4, 5:5, 6:5, 5:6, 5:5, 5:6, 6:5, 5:5, 4:5",
            "3:5, 5:5, 5:6, 6:5, 5:5, 5:6, 6:6, 7:7",
        )
        testS2PolylineVector(input, input, layerOptions)
    }

    @Test
    fun earlyWalkTermination() {
        // This checks idempotency for cases where earlier polylines in the input
        // happen to terminate in the middle of later polylines.  This requires
        // building non-maximal polylines.
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.polylineType = PolylineType.WALK
        val input = listOf(
            "0:1, 1:1",
            "1:0, 1:1, 1:2",
            "0:2, 1:2, 2:2",
            "2:1, 2:2, 2:3"
        )
        testS2PolylineVector(input, input, layerOptions)
    }

    @Test
    fun inputEdgeStartsMultipleLoops() {
        // A single input edge is split into several segments by removing portions
        // of it, and then each of those segments becomes one edge of a loop.
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.polylineType = PolylineType.WALK
        layerOptions.setSiblingPairs(SiblingPairs.DISCARD)
        val builderOptions = S2Builder.Options()
        builderOptions.snapFunction = IntLatLngSnapFunction(7)
        val input = listOf(
            "0:10, 0:0",
            "0:6, 1:6, 1:7, 0:7, 0:8",
            "0:8, 1:8, 1:9, 0:9, 0:10",
            "0:2, 1:2, 1:3, 0:3, 0:4",
            "0:0, 1:0, 1:1, 0:1, 0:2",
            "0:4, 1:4, 1:5, 0:5, 0:6",
        )
        val expected = listOf(
            "0:1, 0:0, 1:0, 1:1, 0:1",
            "0:3, 0:2, 1:2, 1:3, 0:3",
            "0:5, 0:4, 1:4, 1:5, 0:5",
            "0:7, 0:6, 1:6, 1:7, 0:7",
            "0:9, 0:8, 1:8, 1:9, 0:9",
        )
        testS2PolylineVector(input, expected, layerOptions, builderOptions)
    }

    @Test
    fun simpleEdgeLabels() {
        val builder = S2Builder(S2Builder.Options())
        val output = mutableListOf<S2Polyline>()
        val labelSetIds: PolylineVectorLabelSetIds = mutableListOf()
        val labelSetLexicon = IdSetLexicon()
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.edgeType = EdgeType.UNDIRECTED
        layerOptions.duplicateEdges = DuplicateEdges.MERGE
        builder.startLayer(S2PolylineVectorLayer(output, labelSetIds, labelSetLexicon, layerOptions))
        builder.setLabel(1)
        builder.addPolyline(makePolyline("0:0, 0:1, 0:2"))
        builder.setLabel(2)
        builder.addPolyline(makePolyline("0:3, 0:2, 0:1"))
        builder.clearLabels()
        builder.addPolyline(makePolyline("0:4, 0:5"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected: List<List<List<Int>>> = listOf(listOf(listOf(1), listOf(1, 2), listOf(2)), listOf(listOf()))
        assertThat(labelSetIds.size).isEqualTo(expected.size)
        for (i in expected.indices) {
            assertThat(labelSetIds[i].size).isEqualTo(expected[i].size)
            for (j in expected[i].indices) {
                assertThat(labelSetLexicon.idSet(labelSetIds[i][j]).size()).isEqualTo(expected[i][j].size)
                for ((k, label) in labelSetLexicon.idSet(labelSetIds[i][j]).withIndex()) {
                    assertThat(label).isEqualTo(expected[i][j][k])
                }
            }
        }
    }

    @Test
    fun indexedS2PolylineVectorLayerAddsShapes() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PolylineVectorLayer(index))
        val polyline0Str = "0:0, 1:1"
        val polyline1Str = "2:2, 3:3"
        builder.addPolyline(makePolyline(polyline0Str))
        builder.addPolyline(makePolyline(polyline1Str))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(2)
        val polyline0 = (index.shape(0) as S2Polyline.Shape).polyline
        val polyline1 = (index.shape(1) as S2Polyline.Shape).polyline
        assertThat(S2TextParser.toString(polyline0)).isEqualTo(polyline0Str)
        assertThat(S2TextParser.toString(polyline1)).isEqualTo(polyline1Str)
    }

}
