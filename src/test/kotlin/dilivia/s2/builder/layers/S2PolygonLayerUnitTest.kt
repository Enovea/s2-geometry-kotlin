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
import dilivia.s2.S2Factory.makeLoop
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.IsFullPolygon
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.S2Builder
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

// Since we don't expect to have any crossing edges, the key for each edge is
// simply the sum of its endpoints.  This key has the advantage of being
// unchanged when the endpoints of an edge are swapped.
typealias S2PolygonLayerEdgeLabelMap = MutableMap<S2Point, MutableSet<Int>>


class S2PolygonLayerUnitTest {

    private fun testS2Polygon(inputStrs: List<String>, expectedStr: String, edgeType: EdgeType) {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output, options = S2PolygonLayer.Options(edgeType)))
        var isFull = false
        for (inputStr in inputStrs) {
            if (inputStr == "full") isFull = true
            builder.addPolygon(S2TextParser.makeVerbatimPolygon(inputStr))
        }
        builder.addIsFullPolygonPredicate(IsFullPolygon(isFull))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        // The input strings in tests may not be in normalized form, so we build an
        // S2Polygon and convert it back to a string.
        val expected = S2TextParser.makePolygon(expectedStr)
        assertThat(S2TextParser.toString(output)).isEqualTo(S2TextParser.toString(expected))
    }

    private fun testS2Polygon(inputStrs: List<String>, expectedStr: String) {
        testS2Polygon(inputStrs, expectedStr, EdgeType.DIRECTED)
        testS2Polygon(inputStrs, expectedStr, EdgeType.UNDIRECTED)
    }

    private fun testS2PolygonUnchanged(inputStr: String) {
        testS2Polygon(listOf(inputStr), inputStr)
    }

    // Unlike the methods above, the input consists of a set of *polylines*.
    private fun testS2PolygonError(inputStrs: List<String>, expectedError: Int, edgeType: EdgeType) {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        val options = S2PolygonLayer.Options(edgeType)
        options.validate = true
        builder.startLayer(S2PolygonLayer(output, options = options))
        for (input_str in inputStrs) {
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isFalse()
        assertThat(error.code).isEqualTo(expectedError)
    }

    private fun testS2PolygonError(inputStrs: List<String>, expectedError: Int) {
        testS2PolygonError(inputStrs, expectedError, EdgeType.DIRECTED)
        testS2PolygonError(inputStrs, expectedError, EdgeType.UNDIRECTED)
    }

    @Test
    fun empty() {
        testS2PolygonUnchanged("")
    }

    @Test
    fun full() {
        testS2PolygonUnchanged("full")
    }

    @Test
    fun smallLoop() {
        testS2PolygonUnchanged("0:0, 0:1, 1:1")
    }

    @Test
    fun threeLoops() {
        // The second two loops are nested.
        testS2PolygonUnchanged("0:1, 1:1, 0:0; 3:3, 3:6, 6:6, 6:3; 4:4, 4:5, 5:5, 5:4")
    }

    @Test
    fun partialLoop() {
        testS2PolygonError(listOf("0:1, 2:3, 4:5"), S2Error.BUILDER_EDGES_DO_NOT_FORM_LOOPS)
    }

    @Test
    fun invalidPolygon() {
        testS2PolygonError(listOf("0:0, 0:10, 10:0, 10:10, 0:0"), S2Error.LOOP_SELF_INTERSECTION)
    }

    @Test
    fun duplicateInputEdges() {
        // Check that S2PolygonLayer can assemble polygons even when there are
        // duplicate edges (after sibling pairs are removed), and then report the
        // duplicate edges as an error.
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        val options = S2PolygonLayer.Options()
        options.validate = true
        builder.startLayer(S2PolygonLayer(output, options = options))
        builder.addPolyline(makePolyline("0:0, 0:2, 2:2, 1:1, 0:2, 2:2, 2:0, 0:0"))
        val error = S2Error()
        assertThat(builder.build(error)).isFalse()
        assertThat(error.code).isEqualTo(S2Error.POLYGON_LOOPS_SHARE_EDGE)
        assertThat(output.numLoops()).isEqualTo(2)
        val loop0 = makeLoop("0:0, 0:2, 2:2, 2:0")
        val loop1 = makeLoop("0:2, 2:2, 1:1")
        assertThat(loop0 == output.loop(0)).isTrue()
        assertThat(loop1 == output.loop(1)).isTrue()
    }
    
    private fun addPolylineWithLabels(
        polyline: S2Polyline,
        edgeType: EdgeType,
        labelBegin: Int,
        builder: S2Builder,
        edgeLabelMap: S2PolygonLayerEdgeLabelMap
    ) {
        var i = 0
        while (i + 1 < polyline.numVertices) {
            val label = labelBegin + i
            builder.setLabel(label)
            // With undirected edges, reverse the direction of every other input edge.
            val dir = if (edgeType == EdgeType.DIRECTED) 1 else (i and 1)
            builder.addEdge(polyline.vertex(i + (1 - dir)), polyline.vertex(i + dir))
            val key = polyline.vertex(i) + polyline.vertex(i + 1)
            edgeLabelMap[key] = edgeLabelMap.getOrDefault(key, mutableSetOf()).let { it.add(label); it }
            ++i
        }
    }

    private fun testEdgeLabels(edge_type: EdgeType) {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        val labelSetIds: MutableList<MutableList<LabelSetId>> = mutableListOf()
        val labelSetLexicon = IdSetLexicon()
        builder.startLayer(S2PolygonLayer(output, labelSetIds, labelSetLexicon, S2PolygonLayer.Options(edge_type)))

        // We use a polygon consisting of 3 loops.  The loops are reordered and
        // some of the loops are inverted during S2Polygon construction.
        val edgeLabelMap: S2PolygonLayerEdgeLabelMap = mutableMapOf()
        addPolylineWithLabels(
            makePolyline("0:0, 9:1, 1:9, 0:0, 2:8, 8:2, 0:0, 0:10, 10:10, 10:0, 0:0"),
            edge_type,
            0,
            builder,
            edgeLabelMap
        )
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expectedLoopSizes = listOf(4, 3, 3)
        assertThat(labelSetIds.size).isEqualTo(expectedLoopSizes.size)
        for (i in expectedLoopSizes.indices) {
            assertThat(labelSetIds[i].size).isEqualTo(expectedLoopSizes[i])
            for (j in 0 until labelSetIds[i].size) {
                val key = output.loop(i).vertex(j) + output.loop(i).vertex(j + 1)
                val expectedLabels = edgeLabelMap.getOrDefault(key, mutableSetOf())
                assertThat(labelSetLexicon.idSet(labelSetIds[i][j]).size()).isEqualTo(expectedLabels.size)
                assertThat(labelSetLexicon.idSet(labelSetIds[i][j])).containsAll(expectedLabels)
            }
        }
    }

    @Test
    fun directedEdgeLabels() {
        testEdgeLabels(EdgeType.DIRECTED)
    }

    @Test
    fun undirectedEdgeLabels() {
        testEdgeLabels(EdgeType.UNDIRECTED)
    }

    @Test
    fun threeLoopsIntoOne() {
        // Three loops (two shells and one hole) that combine into one.
        testS2Polygon(
            listOf("10:0, 0:0, 0:10, 5:10, 10:10, 10:5", "0:10, 0:15, 5:15, 5:10", "10:10, 5:10, 5:5, 10:5"),
            "10:5, 10:0, 0:0, 0:10, 0:15, 5:15, 5:10, 5:5"
        )
    }

    @Test
    fun trianglePyramid() {
        // A big CCW triangle containing 3 CW triangular holes.  The whole thing
        // looks like a pyramid of nine triangles.  The output consists of 6
        // positive triangles with no holes.
        testS2Polygon(
            listOf(
                "0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1",
                "0:2, 1:1, 1:3",
                "0:4, 1:3, 1:5",
                "1:3, 2:2, 2:4"
            ),
            "0:4, 0:6, 1:5; 2:4, 3:3, 2:2; 2:2, 1:1, 1:3; 1:1, 0:0, 0:2; 1:3, 0:2, 0:4; 1:3, 1:5, 2:4"
        )
    }

    @Test
    fun complexNesting() {
        // A complex set of nested polygons, with the loops in random order and the
        // vertices in random cyclic order within each loop.  This test checks that
        // the order (after S2Polygon.InitNested is called) is preserved exactly,
        // whether directed or undirected edges are used.
        testS2PolygonUnchanged(
            "47:15, 47:5, 5:5, 5:15; " +
                    "35:12, 35:7, 27:7, 27:12; " +
                    "1:50, 50:50, 50:1, 1:1; " +
                    "42:22, 10:22, 10:25, 42:25; " +
                    "47:30, 47:17, 5:17, 5:30; " +
                    "7:27, 45:27, 45:20, 7:20; " +
                    "37:7, 37:12, 45:12, 45:7; " +
                    "47:47, 47:32, 5:32, 5:47; " +
                    "50:60, 50:55, 1:55, 1:60; " +
                    "25:7, 17:7, 17:12, 25:12; " +
                    "7:7, 7:12, 15:12, 15:7"
        )
    }

    @Test
    fun fiveLoopsTouchingAtOneCommonPoint() {
        // Five nested loops that touch at one common point.
        testS2PolygonUnchanged("0:0, 0:10, 10:10, 10:0; 0:0, 1:9, 9:9, 9:1; 0:0, 2:8, 8:8, 8:2; 0:0, 3:7, 7:7, 7:3; 0:0, 4:6, 6:6, 6:4")
    }

    @Test
    fun fourNestedDiamondsTouchingAtTwoPointsPerPair() {
        // Four diamonds nested inside each other, where each diamond shares two
        // vertices with the diamond inside it and shares its other two vertices
        // with the diamond that contains it.  The resulting shape looks vaguely
        // like an eye made out of chevrons.
        testS2Polygon(
            listOf(
                "0:10, -10:0, 0:-10, 10:0",
                "0:-20, -10:0, 0:20, 10:0",
                "0:-10, -5:0, 0:10, 5:0",
                "0:5, -5:0, 0:-5, 5:0"
            ),
            "10:0, 0:10, -10:0, 0:20; 0:-20, -10:0, 0:-10, 10:0; 5:0, 0:-10, -5:0, 0:-5; 0:5, -5:0, 0:10, 5:0"
        )
    }

    @Test
    fun sevenDiamondsTouchingAtOnePointPerPair() {
        // Seven diamonds nested within each other touching at one
        // point between each nested pair.
        testS2PolygonUnchanged(
            "0:-70, -70:0, 0:70, 70:0; 0:-70, -60:0, 0:60, 60:0; " +
                    "0:-50, -60:0, 0:50, 50:0; 0:-40, -40:0, 0:50, 40:0; 0:-30, -30:0, 0:30, 40:0; " +
                    "0:-20, -20:0, 0:30, 20:0; 0:-10, -20:0, 0:10, 10:0"
        )
    }

    @Test
    fun indexedS2PolygonLayerAddsShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PolygonLayer(index))
        val polygonStr = "0:0, 0:10, 10:0"
        builder.addPolygon(S2TextParser.makePolygon(polygonStr))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(1)
        val polygon = (index.shape(0) as S2Polygon.Shape).polygon
        assertThat(S2TextParser.toString(polygon)).isEqualTo(polygonStr)
    }

    @Test
    fun indexedS2PolygonLayerIgnoresEmptyShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PolygonLayer(index))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isZero()
    }

}
