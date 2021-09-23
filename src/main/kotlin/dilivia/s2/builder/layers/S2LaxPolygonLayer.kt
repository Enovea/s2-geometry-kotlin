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

import dilivia.PreConditions.requireArgument
import dilivia.collections.isSorted
import dilivia.s2.S2Error
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeIdSetId
import dilivia.s2.builder.Label
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.S2FindPolygonDegeneracies.findPolygonDegeneracies
import dilivia.s2.builder.S2FindPolygonDegeneracies.isFullyDegenerate
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.EdgeLoop
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.shape.S2LaxPolygonShape
import mu.KotlinLogging

//
// Note that there are two supported output types for polygons: S2Polygon and
// S2LaxPolygonShape.  Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use S2LaxPolygonShape if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from S2LaxPolygonShape to S2Polygon you will need to use S2Builder
// again.
//
// Similarly, there are two supported output formats for polygon meshes:
// S2PolygonMesh and S2LaxPolygonShapeVector.  Use S2PolygonMesh if you need
// to be able to determine which polygons are adjacent to each edge or vertex;
// otherwise use S2LaxPolygonShapeVector, which uses less memory and is faster
// to construct.

// A layer type that assembles edges (directed or undirected) into an
// S2LaxPolygonShape.  Returns an error if the edges cannot be assembled into
// loops.
//
// If the input edges are directed, they must be oriented such that the
// polygon interior is to the left of all edges.  Directed edges are always
// preferred (see S2Builder::EdgeType).
//
// LaxPolygonLayer is implemented such that if the input to S2Builder is a
// polygon and is not modified, then the output has the same cyclic ordering
// of loop vertices and the same loop ordering as the input polygon.
//
// If the given edge graph is degenerate (i.e., it consists entirely of
// degenerate edges and sibling pairs), then the IsFullPolygonPredicate
// associated with the edge graph is called to determine whether the output
// polygon should be empty (possibly with degenerate shells) or full (possibly
// with degenerate holes).  This predicate can be specified as part of the
// S2Builder input geometry.

// Specifies that a polygon should be constructed using the given options,
// and that any labels attached to the input edges should be returned in
// "label_set_ids" and "label_set_lexicion".
//
// The labels associated with the edge "polygon.chain_edge(i, j)"
// can be retrieved as follows:
//
//   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
class LaxPolygonLayer(
        val polygon: S2LaxPolygonShape,
        val labelSetIds: MutableList<MutableList<LabelSetId>>? = null,
        val labelSetLexicon: IdSetLexicon? = null,
        val options: Options = Options()
) : Layer() {

    init {
        require((labelSetIds == null) == (labelSetLexicon == null))
    }

    enum class DegenerateBoundaries {
        DISCARD, DISCARD_HOLES, DISCARD_SHELLS, KEEP
    }

    data class Options(

            // Indicates whether the input edges provided to S2Builder are directed or
            // undirected.  Directed edges should be used whenever possible (see
            // S2Builder::EdgeType for details).
            //
            // If the input edges are directed, they should be oriented so that the
            // polygon interior is to the left of all edges.  This means that for a
            // polygon with holes, the outer loops ("shells") should be directed
            // counter-clockwise while the inner loops ("holes") should be directed
            // clockwise.  Note that S2Builder::AddPolygon() does this automatically.
            //
            // DEFAULT: S2Builder::EdgeType::DIRECTED
            var edgeType: EdgeType = EdgeType.DIRECTED,

            // Specifies whether degenerate boundaries should be discarded or kept.
            // (A degenerate boundary consists of either a sibling edge pair or an
            // edge from a vertex to itself.)  Optionally, degenerate boundaries may
            // be kept only if they represent shells, or only if they represent holes.
            //
            // This option is useful for normalizing polygons with various boundary
            // conditions.  For example, DISCARD_HOLES can be used to normalize closed
            // polygons (those that include their boundary), since degenerate holes do
            // not affect the set of points contained by such polygons.  Similarly,
            // DISCARD_SHELLS can be used to normalize polygons with open boundaries.
            // DISCARD is used to normalize polygons with semi-open boundaries (since
            // degenerate loops do not affect point containment in that case), and
            // finally KEEP is useful for working with any type of polygon where
            // degeneracies are assumed to contain an infinitesmal interior.  (This
            // last model is the most useful for working with simplified geometry,
            // since it maintains the closest fidelity to the original geometry.)
            //
            // DEFAULT: DegenerateBoundaries::KEEP
            var degenerateBoundaries: DegenerateBoundaries = DegenerateBoundaries.KEEP

    )

    override fun graphOptions(): GraphOptions =
            if (options.degenerateBoundaries == DegenerateBoundaries.DISCARD) {
                // There should not be any duplicate edges, but if there are then we keep
                // them since this yields more comprehensible error messages.
                GraphOptions(options.edgeType, DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.DISCARD)
            } else {
                // Keep at most one copy of each sibling pair and each isolated vertex.
                // We need DuplicateEdges::MERGE because DegenerateEdges::DISCARD_EXCESS
                // can still keep multiple copies (it only discards degenerate edges that
                // are connected to non-degenerate edges).
                GraphOptions(options.edgeType, DegenerateEdges.DISCARD_EXCESS, DuplicateEdges.MERGE, SiblingPairs.DISCARD_EXCESS)
            }

    override fun build(g: Graph, error: S2Error) {
        labelSetIds?.clear()
        if (g.options.edgeType == EdgeType.DIRECTED) {
            buildDirected(g, error)
        } else {
            error.init(S2Error.UNIMPLEMENTED, "Undirected edges not supported yet")
        }
    }

    private fun appendPolygonLoops(g: Graph, edge_loops: List<EdgeLoop>, loops: MutableList<MutableList<S2Point>>) {
        for (edge_loop in edge_loops) {
            val vertices = ArrayList<S2Point>(edge_loop.size)
            for (edge_id in edge_loop) {
                vertices.add(g.vertex(g.edge(edge_id).first))
            }
            loops.add(vertices)
        }
    }

    private fun appendEdgeLabels(g: Graph, edgeLoops: List<EdgeLoop>) {
        if (labelSetIds == null) return

        val labels = mutableListOf<Label>()  // Temporary storage for labels.
        val fetcher = Graph.LabelFetcher(g, options.edgeType)
        for (edge_loop in edgeLoops) {
            val loopLabelSetIds = ArrayList<LabelSetId>(edge_loop.size)
            for (edge_id in edge_loop) {
                fetcher.fetch(edge_id, labels)
                loopLabelSetIds.add(labelSetLexicon!!.add(labels))
            }
            labelSetIds.add(loopLabelSetIds)
        }
    }

    private fun buildDirected(graph: Graph, error: S2Error) {
        var g = graph
        // Some cases are implemented by constructing a new graph with certain
        // degenerate edges removed (overwriting "g").  "new_edges" is where the
        // edges for the new graph are stored.
        val newEdges = mutableListOf<Edge>()
        val newInputEdgeIdSetIds = mutableListOf<InputEdgeIdSetId>()
        val loops = mutableListOf<MutableList<S2Point>>()
        val degenerateBoundaries = options.degenerateBoundaries
        if (degenerateBoundaries == DegenerateBoundaries.DISCARD) {
            // This is the easiest case, since there are no degeneracies.
            if (g.numEdges == 0) maybeAddFullLoop(g, loops, error)
        } else if (degenerateBoundaries == DegenerateBoundaries.KEEP) {
            // S2LaxPolygonShape doesn't need to distinguish degenerate shells from
            // holes except when the entire graph is degenerate, in which case we need
            // to decide whether it represents an empty polygons possibly with
            // degenerate shells, or a full polygon possibly with degenerate holes.
            if (isFullyDegenerate(g)) {
                maybeAddFullLoop(g, loops, error)
            }
        } else {
            // For DISCARD_SHELLS and DISCARD_HOLES we first determine whether any
            // degenerate loops of the given type exist, and if so we construct a new
            // graph with those edges removed (overwriting "g").
            val discardHoles = (degenerateBoundaries == DegenerateBoundaries.DISCARD_HOLES)
            val degeneracies = findPolygonDegeneracies(g, error)
            if (!error.isOk()) return
            if (degeneracies.size == g.numEdges) {
                if (degeneracies.isEmpty()) {
                    maybeAddFullLoop(g, loops, error)
                } else if (degeneracies[0].isHole) {
                    loops.add(mutableListOf())  // Full loop.
                }
            }
            val edgesToDiscard = mutableListOf<EdgeId>()
            for (degeneracy in degeneracies) {
                if (degeneracy.isHole == discardHoles) {
                    edgesToDiscard.add(degeneracy.edgeId)
                }
            }
            if (edgesToDiscard.isNotEmpty()) {
                // Construct a new graph that discards the unwanted edges.
                edgesToDiscard.sort()
                discardEdges(g, edgesToDiscard, newEdges, newInputEdgeIdSetIds)
                g = Graph(g.options, g.vertices, newEdges, newInputEdgeIdSetIds,
                        g.inputEdgeIdSetLexicon, g.labelSetIds, g.labelSetLexicon, g.isFullPolygonPredicate)
            }
        }
        val edgeLoops = mutableListOf<EdgeLoop>()
        if (!g.getDirectedLoops(Graph.LoopType.CIRCUIT, edgeLoops, error)) {
            return
        }
        appendPolygonLoops(g, edgeLoops, loops)
        logger.trace { "Loops after polygon loops appending: loops = ${loops.map { l -> l.map { S2LatLng.fromPoint(it)} }} | edgeloops = $edgeLoops" }
        appendEdgeLabels(g, edgeLoops)
        polygon.init(loops)
    }

    companion object {

        private val logger = KotlinLogging.logger(LaxPolygonLayer::class.java.name)

        // Returns all edges of "g" except for those identified by "edges_to_discard".
        fun discardEdges(g: Graph, edgesToDiscard: List<EdgeId>, newEdges: MutableList<Edge>, newInputEdgeIdSetIds: MutableList<InputEdgeIdSetId>) {
            requireArgument { edgesToDiscard.isSorted() }
            newEdges.clear()
            newInputEdgeIdSetIds.clear()
            if (newEdges is ArrayList) newEdges.ensureCapacity(g.numEdges)
            if (newInputEdgeIdSetIds is ArrayList) newInputEdgeIdSetIds.ensureCapacity(g.numEdges)
            var idx = 0
            for (e in 0 until g.numEdges) {
                if (idx != edgesToDiscard.size && e == edgesToDiscard[idx]) {
                    ++idx
                } else {
                    newEdges.add(g.edge(e))
                    newInputEdgeIdSetIds.add(g.inputEdgeIdSetId(e))
                }
            }
            assert(idx == edgesToDiscard.size)
        }

        fun maybeAddFullLoop(g: Graph, loops: MutableList<MutableList<S2Point>>, error: S2Error) {
            if (g.isFullPolygon(error)) {
                loops.add(mutableListOf())  // Full loop.
            }
        }
    }

}

