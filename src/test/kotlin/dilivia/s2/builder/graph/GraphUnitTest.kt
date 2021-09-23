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

package dilivia.s2.builder.graph

import dilivia.collections.get
import dilivia.s2.S2Error
import dilivia.s2.S2TextParser.makeLaxPolyline
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType.DIRECTED
import dilivia.s2.builder.EdgeType.UNDIRECTED
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.InputEdgeIdSetId
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.graph.DegenerateEdges.DISCARD_EXCESS
import dilivia.s2.builder.graph.DuplicateEdges.KEEP
import dilivia.s2.builder.graph.DuplicateEdges.MERGE
import dilivia.s2.builder.graph.Graph.PolylineType
import dilivia.s2.builder.graph.SiblingPairs.CREATE
import dilivia.s2.builder.layers.Layer
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import dilivia.s2.builder.graph.SiblingPairs.KEEP as KEEP1

/**
 * Most of S2Builder.Graph is tested by the S2Builder.Layer implementations
 */

class GraphUnitTest {


    // A layer type that copies an S2Builder.Graph into a GraphClone object
// (which owns the underlying data, unlike S2Builder.Graph itself).
    class NoOpLayer(val graphOptions: GraphOptions) : Layer() {

        lateinit var g: Graph

        override fun graphOptions(): GraphOptions = graphOptions

        override fun build(g: Graph, error: S2Error) {
            // Does nothing
            this.g = g
        }

    }

    @Test
    fun getDirectedLoopsDegenerateEdges() {
        val builder = S2Builder(S2Builder.Options())
        val graphOptions = GraphOptions(
            edgeType = DIRECTED,
            degenerateEdges = DISCARD_EXCESS,
            duplicateEdges = KEEP,
            siblingPairs = KEEP1
        )
        val layer = NoOpLayer(graphOptions)
        builder.startLayer(layer)
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        builder.addShape(makeLaxPolyline("0:0, 0:2, 2:2, 2:0, 0:0"))
        builder.addShape(makeLaxPolyline("0:3, 3:3, 0:3"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val g = layer.g
        val loops = mutableListOf<List<EdgeId>>()
        assertThat(g.getDirectedLoops(Graph.LoopType.SIMPLE, loops, error)).isTrue()
        assertThat(loops.size).isEqualTo(3)
        assertThat(loops[0].size).isEqualTo(1)
        assertThat(loops[1].size).isEqualTo(4)
        assertThat(loops[2].size).isEqualTo(2)
    }

    @Test
    fun getDirectedComponentsDegenerateEdges() {
        val builder = S2Builder(S2Builder.Options())
        val graphOptions = GraphOptions(DIRECTED, DISCARD_EXCESS, MERGE, CREATE)
        val layer = NoOpLayer(graphOptions)
        builder.startLayer(layer)
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        builder.addShape(makeLaxPolyline("0:0, 0:2, 2:2, 2:0, 0:0"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val g = layer.g
        val components = mutableListOf<DirectedComponent>()
        assertThat(g.getDirectedComponents(Graph.DegenerateBoundaries.KEEP, components, error)).isTrue()
        assertThat(components.size).isEqualTo(2)
        assertThat(components[0].size).isEqualTo(1)
        assertThat(components[0][0].size).isEqualTo(1)
        assertThat(components[1].size).isEqualTo(2)
        assertThat(components[1][0].size).isEqualTo(4)
        assertThat(components[1][1].size).isEqualTo(4)
    }

    @Test
    fun getUndirectedComponentsDegenerateEdges() {
        val builder = S2Builder(S2Builder.Options())
        val graphOptions = GraphOptions(UNDIRECTED, DISCARD_EXCESS, KEEP, SiblingPairs.DISCARD_EXCESS)
        val layer = NoOpLayer(graphOptions)
        builder.startLayer(layer)
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        builder.addShape(makeLaxPolyline("0:0, 0:2, 2:2, 2:0, 0:0"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val g = layer.g
        val components = mutableListOf<UndirectedComponent>()
        assertThat(g.getUndirectedComponents(Graph.LoopType.CIRCUIT, components, error)).isTrue()
        // The result consists of two components, each with two complements.  Each
        // complement in this example has exactly one loop.  The loops in both
        // complements of the first component have 1 vertex, while the loops in both
        // complements of the second component have 4 vertices.
        assertThat(components.size).isEqualTo(2)
        assertThat(components[0][0].size).isEqualTo(1)
        assertThat(components[0][0][0].size).isEqualTo(1)
        assertThat(components[0][1].size).isEqualTo(1)
        assertThat(components[0][1][0].size).isEqualTo(1)
        assertThat(components[1][0].size).isEqualTo(1)
        assertThat(components[1][0][0].size).isEqualTo(4)
        assertThat(components[1][1].size).isEqualTo(1)
        assertThat(components[1][1][0].size).isEqualTo(4)
    }

    @Test
    fun getPolylinesUndirectedDegeneratePaths() {
        val builder = S2Builder(S2Builder.Options())
        val graphOptions = GraphOptions(UNDIRECTED, DegenerateEdges.KEEP, KEEP, SiblingPairs.KEEP)
        val layer = NoOpLayer(graphOptions)
        builder.startLayer(layer)
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        builder.addShape(makeLaxPolyline("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"))
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val g = layer.g
        val polylines = g.getPolylines(PolylineType.PATH)
        assertThat(polylines.size).isEqualTo(7)
    }

    @Test
    fun getPolylinesUndirectedDegenerateWalks() {
        val builder = S2Builder(S2Builder.Options())
        val graphOptions = GraphOptions(UNDIRECTED, DegenerateEdges.KEEP, KEEP, SiblingPairs.KEEP)
        val layer = NoOpLayer(graphOptions)
        builder.startLayer(layer)
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        builder.addShape(makeLaxPolyline("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"))
        builder.addShape(makeLaxPolyline("1:1, 1:1"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val g = layer.g
        val polylines = g.getPolylines(PolylineType.WALK)
        assertThat(polylines.size).isEqualTo(2)
        assertThat(polylines[0].size).isEqualTo(2)
        assertThat(polylines[1].size).isEqualTo(5)
    }

    fun testProcessEdges(input: List<TestEdge>, expected: List<TestEdge>, options: GraphOptions, expected_code: Int = S2Error.OK) {
        val edges = ArrayList<Edge>()
        val inputIdSetIds = ArrayList<InputEdgeIdSetId>()
        val idSetLexicon = IdSetLexicon()
        for (e in input) {
            edges.add(Edge(e.first, e.second))
            inputIdSetIds.add(idSetLexicon.add(e.input_ids))
        }
        val error = S2Error()
        Graph.processEdges(options, edges, inputIdSetIds, idSetLexicon, error)
        assertThat(expected_code).isEqualTo(error.code)
        assertThat(edges.size).isEqualTo(inputIdSetIds.size)
        for (i in expected.indices) {
            val e = expected[i]
            assertThat(edges.size)
                .withFailMessage("Not enough output edges")
                .isGreaterThan(i)
            assertThat(edges[i])
                .withFailMessage("(edge $i)")
                .isEqualTo(Edge(e.first, e.second))
            val idSet = idSetLexicon.idSet(inputIdSetIds[i])
            assertThat(idSet)
                .withFailMessage("(edge $i)")
                .hasSize(e.input_ids.size)
                .containsAll(e.input_ids)
        }
        assertThat(edges.size)
            .withFailMessage("Too many output edges")
            .isEqualTo(expected.size)
    }
    
    @Test
    fun processEdgesDiscardDegenerateEdges() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.KEEP)
        testProcessEdges(listOf(TestEdge(0, 0), TestEdge(0, 0) ), emptyList(), options)
    }

    @Test
    fun processEdgesKeepDuplicateDegenerateEdges() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.KEEP, KEEP, SiblingPairs.KEEP)
        testProcessEdges(
            listOf(TestEdge(0, 0), TestEdge(0, 0)), 
            listOf(TestEdge(0, 0), TestEdge(0, 0)), 
            options
        )
    }

    @Test
    fun processEdgesMergeDuplicateDegenerateEdges() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.KEEP, MERGE, SiblingPairs.KEEP)
        testProcessEdges(
            listOf(TestEdge( 0, 0, listOf(1)), TestEdge(0, 0, listOf(2))),
            listOf(TestEdge(0, 0, listOf(1, 2))), 
            options
        )
    }

    @Test
    fun processEdgesMergeUndirectedDuplicateDegenerateEdges() {
        // Edge count should be reduced to 2 (i.e., one undirected edge), and all
        // labels should be merged.
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.KEEP, MERGE, SiblingPairs.KEEP)
        testProcessEdges(
            listOf(TestEdge( 0, 0, listOf(1)), TestEdge(0, 0), TestEdge(0, 0),TestEdge(0, 0, listOf(2))),
            listOf(TestEdge(0, 0, listOf(1, 2)), TestEdge(0, 0, listOf(1, 2))),
            options
        )
    }

    @Test
    fun processEdgesConvertedUndirectedDegenerateEdges() {
        // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
        // merges any edge labels.
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.KEEP, KEEP, SiblingPairs.REQUIRE)
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0), TestEdge(0, 0), TestEdge(0, 0,listOf(2))),
            listOf(TestEdge(0, 0, listOf(1, 2)), TestEdge(0, 0, listOf(1, 2))),
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)
    }

    @Test
    fun processEdgesMergeConvertedUndirectedDuplicateDegenerateEdges() {
        // Like the test above, except that we also merge duplicates.
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.KEEP, MERGE, SiblingPairs.REQUIRE)
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0), TestEdge(0, 0), TestEdge(0, 0,listOf(2))),
            listOf(TestEdge(0, 0,listOf(1, 2))),
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)
    }

    @Test
    fun processEdgesDiscardExcessConnectedDegenerateEdges() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD_EXCESS,
        KEEP, SiblingPairs.KEEP);
        // Test that degenerate edges are discarded if they are connnected to any
        // non-degenerate edges (whether they are incoming or outgoing, and whether
        // they are lexicographically before or after the degenerate edge).
        testProcessEdges(listOf(TestEdge(0, 0), TestEdge(0, 1)), listOf(TestEdge(0, 1)), options)
        testProcessEdges(listOf(TestEdge(0, 0), TestEdge(1, 0)), listOf(TestEdge(1, 0)), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 1)), listOf(TestEdge(0, 1)), options)
        testProcessEdges(listOf(TestEdge(1, 0), TestEdge(1, 1)), listOf(TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesDiscardExcessIsolatedDegenerateEdges() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD_EXCESS, KEEP, SiblingPairs.KEEP)
        // Test that DISCARD_EXCESS does not merge any remaining duplicate
        // degenerate edges together.
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0,listOf(2))),
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0,listOf(2))),
            options
        )
    }

    @Test
    fun processEdgesDiscardExcessUndirectedIsolatedDegenerateEdges() {
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.DISCARD_EXCESS, KEEP, SiblingPairs.KEEP)
        // Test that DISCARD_EXCESS does not merge any remaining duplicate
        // undirected degenerate edges together.
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0), TestEdge(0, 0,listOf(2)), TestEdge(0, 0)),
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0), TestEdge(0, 0,listOf(2)), TestEdge(0, 0)),
            options
        )
    }

    @Test
    fun processEdgesDiscardExcessConvertedUndirectedIsolatedDegenerateEdges() {
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.DISCARD_EXCESS, KEEP, SiblingPairs.REQUIRE)
        // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
        // merges edge labels.
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0,listOf(2)), TestEdge(0, 0,listOf(3)), TestEdge(0, 0)),
            listOf(TestEdge(0, 0,listOf(1, 2, 3)), TestEdge(0, 0,listOf(1, 2, 3))),
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)
    }

    @Test
    fun processEdgesSiblingPairsDiscardMergesDegenerateEdgeLabels() {
        // Test that when SiblingPairs.DISCARD or SiblingPairs.DISCARD_EXCESS
        // is specified, the edge labels of degenerate edges are merged together
        // (for consistency, since these options merge the labels of all
        // non-degenerate edges as well).
        val options = GraphOptions(DIRECTED, DegenerateEdges.KEEP, KEEP, SiblingPairs.DISCARD)
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0,listOf(2)), TestEdge(0, 0,listOf(3))),
            listOf(TestEdge(0, 0,listOf(1, 2, 3)), TestEdge(0, 0,listOf(1, 2, 3)), TestEdge(0, 0,listOf(1, 2, 3))),
            options
        )
        options.siblingPairs = SiblingPairs.DISCARD_EXCESS
        testProcessEdges(
            listOf(TestEdge(0, 0,listOf(1)), TestEdge(0, 0,listOf(2)), TestEdge(0, 0,listOf(3))),
            listOf(TestEdge(0, 0,listOf(1, 2, 3)), TestEdge(0, 0,listOf(1, 2, 3)), TestEdge(0, 0,listOf(1, 2, 3))),
            options
        )
    }

    @Test
    fun processEdgesKeepSiblingPairs() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.KEEP)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesMergeDuplicateSiblingPairs() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, MERGE, SiblingPairs.KEEP)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesDiscardSiblingPairs() {
        // Check that matched pairs are discarded, leaving behind any excess edges.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.DISCARD)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), emptyList(), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)), emptyList(), options)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(0, 1)), 
            options
        )
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(1, 0), TestEdge(1, 0)),
            options
        )
    }

    @Test
    fun processEdgesDiscardSiblingPairsMergeDuplicates() {
        // Check that matched pairs are discarded, and then any remaining edges
        // are merged.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, MERGE, SiblingPairs.DISCARD)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)), emptyList(), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1)), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)), listOf(TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesDiscardUndirectedSiblingPairs() {
        // An undirected sibling pair consists of four edges, two in each direction
        // (see s2builder.h).  Since undirected edges always come in pairs, this
        // means that the result always consists of either 0 or 2 edges.
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.DISCARD,
        KEEP, SiblingPairs.DISCARD);
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)), emptyList(), options)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options
        )
    }

    @Test
    fun processEdgesDiscardExcessSiblingPairs() {
        // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
        // result would otherwise be empty.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.DISCARD_EXCESS)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)),
            options
        )
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(0, 1)),
            options
        )
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(1, 0), TestEdge(1, 0)),
            options
        )
    }

    @Test
    fun processEdgesDiscardExcessSiblingPairsMergeDuplicates() {
        // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
        // result would otherwise be empty.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, MERGE, SiblingPairs.DISCARD_EXCESS)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)),
            options
        )
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1)), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)), listOf(TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesDiscardExcessUndirectedSiblingPairs() {
        // Like SiblingPairs.DISCARD, except that one undirected sibling pair
        // (4 edges) is kept if the result would otherwise be empty.
        val options = GraphOptions(UNDIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.DISCARD_EXCESS)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)), 
            options
        )
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)),
            options
        )
    }

    @Test
    fun processEdgesCreateSiblingPairs() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, CREATE)
        testProcessEdges(listOf(TestEdge(0, 1)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1)),
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)),
            options
        )
    }

    @Test
    fun processEdgesRequireSiblingPairs() {
        // Like SiblingPairs.CREATE, but generates an error.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, SiblingPairs.REQUIRE)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(1, 0)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(
            listOf(TestEdge(0, 1)), 
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options, 
            S2Error.BUILDER_MISSING_EXPECTED_SIBLING_EDGES
        )
    }

    @Test
    fun processEdgesCreateUndirectedSiblingPairs() {
        // An undirected sibling pair consists of 4 edges, but SiblingPairs.CREATE
        // also converts the graph to EdgeType.DIRECTED and cuts the number of
        // edges in half.
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, KEEP, CREATE)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)

        options.edgeType = UNDIRECTED
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)

        options.edgeType = UNDIRECTED
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)), 
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)
    }

    @Test
    fun processEdgesCreateSiblingPairsMergeDuplicates() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, MERGE, CREATE)
        testProcessEdges(listOf(TestEdge(0, 1)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
        testProcessEdges(listOf(TestEdge(0, 1), TestEdge(0, 1)), listOf(TestEdge(0, 1), TestEdge(1, 0)), options)
    }

    @Test
    fun processEdgesCreateUndirectedSiblingPairsMergeDuplicates() {
        val options = GraphOptions(DIRECTED, DegenerateEdges.DISCARD, MERGE, CREATE)
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)

        options.edgeType = UNDIRECTED
        testProcessEdges(
            listOf(TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)),
            listOf(TestEdge(0, 1), TestEdge(1, 0)), 
            options
        )
        assertThat(options.edgeType).isEqualTo(DIRECTED)
    }


    data class TestEdge(val first: VertexId, val second: VertexId, val input_ids: List<InputEdgeId> = emptyList())

}

