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
package dilivia.s2.builder.layers

import dilivia.s2.S2Error
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2TextParser
import dilivia.s2.S2TextParser.makeLaxPolygon
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.IsFullPolygon
import dilivia.s2.builder.IsFullPolygonUnspecified
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.layers.LaxPolygonLayer.DegenerateBoundaries
import dilivia.s2.builder.layers.LaxPolygonLayer.DegenerateBoundaries.DISCARD
import dilivia.s2.builder.layers.LaxPolygonLayer.DegenerateBoundaries.DISCARD_HOLES
import dilivia.s2.builder.layers.LaxPolygonLayer.DegenerateBoundaries.DISCARD_SHELLS
import dilivia.s2.builder.layers.LaxPolygonLayer.DegenerateBoundaries.KEEP
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2LaxPolygonShape
import dilivia.s2.shape.S2Shape
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2LaxPolygonLayerUnitTest {

    private val logger = KotlinLogging.logger {  }

    fun testLaxPolygon(
        inputStr: String,
        expectedStr: String,
        edge_type: EdgeType,
        degenerate_boundaries: DegenerateBoundaries
    ) {
        //SCOPED_TRACE(edge_type == EdgeType.DIRECTED ? "DIRECTED" : "UNDIRECTED");
        //SCOPED_TRACE(ToString(degenerate_boundaries));
        val builder = S2Builder(options = S2Builder.Options())
        val output = S2LaxPolygonShape()
        val options = LaxPolygonLayer.Options()
        options.edgeType = edge_type
        options.degenerateBoundaries = degenerate_boundaries
        builder.startLayer(LaxPolygonLayer(output, options = options))

        val polygon = makeLaxPolygon(inputStr)
        builder.addShape(polygon)

        // In order to construct polygons that are full except possibly for a
        // collection of degenerate holes, we must supply S2Builder with a predicate
        // that distinguishes empty polygons from full ones (modulo degeneracies).
        var has_full_loop = false
        for (i in 0 until polygon.numLoops()) {
            if (polygon.numLoopVertices(i) == 0) {
                has_full_loop = true
            }
        }
        builder.addIsFullPolygonPredicate(IsFullPolygon(has_full_loop))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val actual_str = S2TextParser.toString(output, "; ")
        assertThat(actual_str).isEqualTo(expectedStr)
    }

    private fun testLaxPolygon(input_str: String, expectedStr: String, degenerate_boundaries: DegenerateBoundaries) {
        testLaxPolygon(input_str, expectedStr, EdgeType.DIRECTED, degenerate_boundaries)
        // TODO(ericv): Implement.
        //testLaxPolygon(input_str, expected_str, EdgeType.UNDIRECTED, degenerate_boundaries)
    }

    fun testLaxPolygonUnchanged(input_str: String, degenerate_boundaries: DegenerateBoundaries) {
        testLaxPolygon(input_str, input_str, degenerate_boundaries);
    }

    fun allDegenerateBoundaries(): List<DegenerateBoundaries> {
        return listOf(DISCARD, DISCARD_HOLES, DISCARD_SHELLS, KEEP)
    }

    @Test
    fun empty() {
        for (db in allDegenerateBoundaries()) {
            logger.info { "Test empty with DegenerateBoundaries = $db" }
            testLaxPolygonUnchanged("", db)
        }
    }

    @Test
    fun full() {
        for (db in allDegenerateBoundaries()) {
            testLaxPolygonUnchanged("full", db)
        }
    }

    @Test
    fun oneNormalShell() {
        for (db in allDegenerateBoundaries()) {
            testLaxPolygonUnchanged("0:0, 0:1, 1:1", db);
        }
    }

    @Test
    fun isFullPolygonPredicateNotCalled() {
        // Test that the IsFullPolygonPredicate is not called when at least one
        // non-degenerate loop is present.
        for (degenerate_boundaries in allDegenerateBoundaries()) {
            val builder = S2Builder(S2Builder.Options())
            val output = S2LaxPolygonShape()
            val options = LaxPolygonLayer.Options()
            options.edgeType = EdgeType.DIRECTED
            options.degenerateBoundaries = degenerate_boundaries
            builder.startLayer(LaxPolygonLayer(output, options = options))
            val polygon = makeLaxPolygon("0:0, 0:1, 1:1")
            builder.addShape(polygon)
            // If the predicate is called, it will return an error.
            builder.addIsFullPolygonPredicate(IsFullPolygonUnspecified())
            val error = S2Error()
            assertThat(builder.build(error)).isTrue()
        }
    }

    @Test
    fun twoNormalShellsOneNormalHole() {
        // The second two loops are nested.  Note that S2LaxPolygon and S2Polygon
        // require opposite vertex orderings for holes.
        for (db in allDegenerateBoundaries()) {
            testLaxPolygonUnchanged("0:1, 1:1, 0:0; 3:3, 3:6, 6:6, 6:3; 4:4, 5:4, 5:5, 4:5", db)
        }
    }

    @Test
    fun allDegenerateShells() {
        for (db in listOf(KEEP, DISCARD_HOLES)) {
            testLaxPolygonUnchanged("1:1; 2:2, 3:3", db)
        }
        for (db in listOf(DISCARD, DISCARD_SHELLS)) {
            testLaxPolygon("1:1; 2:2, 3:3", "", db)
        }
    }

    @Test
    fun allDegenerateHoles() {
        for (db in listOf(KEEP, DISCARD_SHELLS)) {
            testLaxPolygonUnchanged("full; 1:1; 2:2, 3:3", db)
        }
        for (db in listOf(DISCARD, DISCARD_HOLES)) {
            testLaxPolygon("full; 1:1; 2:2, 3:3", "full", db)
        }
    }

    @Test
    fun someDegenerateShells() {
        val kNormal = "0:0, 0:9, 9:0; 1:1, 7:1, 1:7"
        val kInput = kNormal + "; 3:2; 2:2, 2:3"
        testLaxPolygonUnchanged(kInput, KEEP);
        testLaxPolygonUnchanged(kInput, DISCARD_HOLES);
        testLaxPolygon(kInput, kNormal, DISCARD);
        testLaxPolygon(kInput, kNormal, DISCARD_SHELLS);
    }

    @Test
    fun someDegenerateHoles() {
        for (db in listOf(KEEP, DISCARD_SHELLS)) {
            testLaxPolygonUnchanged("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", db)
        }
        for (db in listOf(DISCARD, DISCARD_HOLES)) {
            testLaxPolygon("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", "0:0, 0:9, 9:0", db)
        }
    }

    @Test
    fun normalAndDegenerateShellsAndHoles() {
        // We start with two normal shells and one normal hole.
        val kNormal = "0:0, 0:9, 9:9, 9:0; 0:10, 0:19, 9:19, 9:10; 1:11, 8:11, 8:18, 1:18"
        // These are the same loops augmented with degenerate interior filaments
        // (holes).  Note that one filament connects the second shell and hole
        // above, transforming them into a single loop.
        val kNormalWithDegenHoles =
            "0:0, 0:9, 1:8, 1:7, 1:8, 0:9, 9:9, 9:0; 0:10, 0:19, 9:19, 9:10, 0:10, 1:11, 8:11, 8:18, 1:18, 1:11"
        // Then we add other degenerate shells and holes, including a sibling pair
        // that connects the two shells above.
        val kDegenShells = "0:9, 0:10; 2:12; 3:13, 3:14; 20:20; 10:0, 10:1"
        val kDegenHoles = "2:5; 3:6, 3:7; 8:8"
        val kInput = "$kNormalWithDegenHoles; $kDegenShells; $kDegenHoles"
        testLaxPolygon(kInput, kNormal, DISCARD);
        testLaxPolygon(kInput, "$kNormal; $kDegenShells", DISCARD_HOLES)
        testLaxPolygon(kInput, "$kNormalWithDegenHoles; $kDegenHoles", DISCARD_SHELLS)
        testLaxPolygon(kInput, kInput, KEEP)
    }

    @Test
    fun partialLoop() {
        val builder = S2Builder(S2Builder.Options())
        val output = S2LaxPolygonShape()
        builder.startLayer(LaxPolygonLayer(output))
        builder.addPolyline(makePolyline("0:1, 2:3, 4:5"))
        val error = S2Error()
        assertThat(builder.build(error)).isFalse()
        assertThat(error.code).isEqualTo(S2Error.BUILDER_EDGES_DO_NOT_FORM_LOOPS)
        assertThat(output.isEmpty()).isTrue()
    }

// TODO(ericv): Implement validation of S2LaxPolygonShape.
    /*
  @Test fun InvalidPolygon() {
    val builder = S2Builder(S2Builder.Options())
    val output = S2LaxPolygonShape()
    val options = LaxPolygonLayer.Options()
    options.set_validate(true);
    builder.startLayer(LaxPolygonLayer(output, options));
    builder.addPolyline(makePolyline("0:0, 0:10, 10:0, 10:10, 0:0"));
    val error = S2Error()
    assertThat(builder.build(error)).isFalse()
    assertThat(error.code).isEqualTo(S2Error::LOOP_SELF_INTERSECTION)
  }
     */

    @Test
    fun duplicateInputEdges() {
        // Check that LaxPolygonLayer removes duplicate edges in such a way that
        // degeneracies are not lost.
        val builder = S2Builder(S2Builder.Options())
        val output = S2LaxPolygonShape()
        val options = LaxPolygonLayer.Options()
        options.degenerateBoundaries = KEEP
        builder.startLayer(LaxPolygonLayer(output, options = options))
        builder.addShape(makeLaxPolygon("0:0, 0:5, 5:5, 5:0"))
        builder.addPoint(makePoint("0:0"))
        builder.addPoint(makePoint("1:1"))
        builder.addPoint(makePoint("1:1"))
        builder.addShape(makeLaxPolygon("2:2, 2:3"))
        builder.addShape(makeLaxPolygon("2:2, 2:3"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output, "; ")).isEqualTo("0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 2:3")
    }


    fun getKey(edge: Edge, edge_type: EdgeType): Edge {
        var key = edge
        // For undirected edges, sort the vertices in lexicographic order.
        if (edge_type == EdgeType.UNDIRECTED && edge.v0 > edge.v1) {
            key = edge.copy(v0 = edge.v1, v1 = edge.v0)
        }
        return key
    }

    fun addShapeWithLabels(shape: S2Shape, edge_type: EdgeType, builder: S2Builder, edge_label_map: EdgeLabelMap) {
        val kLabelBegin = 1234  // Arbitrary.
        for (e in 0 until shape.numEdges) {
            val label = kLabelBegin + e
            builder.setLabel(label)
            // For undirected edges, reverse the direction of every other input edge.
            var edge = shape.edge(e)
            if (edge_type == EdgeType.UNDIRECTED && (e and 1 != 0)) {
                edge = edge.reversed()
            }
            builder.addEdge(edge.v0, edge.v1)
            edge_label_map.compute(getKey(edge, edge_type)) { _, set ->
                if (set == null) mutableSetOf(label) else {
                    set.add(label); set
                }
            }
        }
    }

    // Converts "input_str" to an S2LaxPolygonShape, assigns labels to its edges,
// then uses LaxPolygonLayer with the given arguments to build a new
// S2LaxPolygonShape and verifies that all edges have the expected labels.
// (This function does not test whether the output edges are correct.)
    fun testEdgeLabels(input_str: String, edge_type: EdgeType, degenerate_boundaries: DegenerateBoundaries) {
        val builder = S2Builder(S2Builder.Options())
        val output = S2LaxPolygonShape()
        val label_set_ids = mutableListOf(mutableListOf<LabelSetId>())
        val label_set_lexicon = IdSetLexicon()
        val options = LaxPolygonLayer.Options()
        options.edgeType = edge_type
        options.degenerateBoundaries = degenerate_boundaries
        builder.startLayer(LaxPolygonLayer(output, label_set_ids, label_set_lexicon, options))

        val edge_label_map: EdgeLabelMap = mutableMapOf()
        addShapeWithLabels(makeLaxPolygon(input_str), edge_type, builder, edge_label_map)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        for (i in 0 until output.numChains) {
            for (j in 0 until output.chain(i).length) {
                val edge = output.chainEdge(i, j)
                val expected_labels = edge_label_map[getKey(edge, edge_type)]!!
                assertThat(label_set_lexicon.idSet(label_set_ids[i][j]).size()).isEqualTo(expected_labels.size)
                assertThat(label_set_lexicon.idSet(label_set_ids[i][j])).containsExactly(*expected_labels.toTypedArray())
            }
        }
    }

    @Test
    fun edgeLabels() {
        // TODO(ericv): Implement EdgeType.UNDIRECTED.
        for (edge_type in listOf(EdgeType.DIRECTED)) {
            for (db in allDegenerateBoundaries()) {
                // Test a polygon with normal and degenerate shells and holes.  Note
                // that this S2LaxPolygonShape has duplicate edges and is therefore not
                // valid in most contexts.
                testEdgeLabels(
                    "1:1, 1:2; 0:0, 0:9, 9:9, 9:0; 1:2, 1:1; 3:3, 8:3, 8:8, 3:8; 4:4; 4:5, 5:5; 4:4",
                    edge_type,
                    db
                )
            }
        }
    }

    @Test
    fun indexedLaxPolygonLayerAddsShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedLaxPolygonLayer(index))
        val polygon_str = "0:0, 0:10, 10:0";
        builder.addPolygon(S2TextParser.makePolygon(polygon_str))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(1)
        val polygon: S2LaxPolygonShape = index.shape(0) as S2LaxPolygonShape
        assertThat(S2TextParser.toString(polygon)).isEqualTo(polygon_str)
    }

    @Test
    fun indexedLaxPolygonLayerIgnoresEmptyShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedLaxPolygonLayer(index))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(0)
    }

}  // namespace
