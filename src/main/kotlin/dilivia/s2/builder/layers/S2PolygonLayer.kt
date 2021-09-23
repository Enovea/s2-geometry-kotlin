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

import dilivia.PreConditions.requireEQ
import dilivia.collections.assignWith
import dilivia.collections.get
import dilivia.collections.reverse
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.Label
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.EdgeLoop
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.graph.UndirectedComponent
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon

//
// Note that there are two supported output types for polygons: S2Polygon and
// S2LaxPolygonShape.  Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use S2LaxPolygonShape if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from a S2LaxPolygonShape to an S2Polygon you will need to use
// S2Builder again.
//
// Similarly, there are two supported output formats for polygon meshes:
// S2LaxPolygonShapeVector and S2PolygonMesh.  Use S2PolygonMesh if you need
// to be able to determine which polygons are adjacent to each edge or vertex;
// otherwise use S2LaxPolygonShapeVector, which uses less memory and is faster
// to construct.

typealias LoopMap = MutableMap<S2Loop, Pair<Int, Boolean>>

// A layer type that assembles edges (directed or undirected) into an
// S2Polygon.  Returns an error if the edges cannot be assembled into loops.
//
// If the input edges are directed, they must be oriented such that the
// polygon interior is to the left of all edges.  Directed edges are always
// preferred (see S2Builder::EdgeType).
//
// Before the edges are assembled into loops, "sibling pairs" consisting of an
// edge and its reverse edge are automatically removed.  Such edge pairs
// represent zero-area degenerate regions, which S2Polygon does not allow.
// (If you need to build polygons with degeneracies, use LaxPolygonLayer
// instead.)
//
// S2PolygonLayer is implemented such that if the input to S2Builder is a
// polygon and is not modified, then the output has the same cyclic ordering
// of loop vertices and the same loop ordering as the input polygon.
//
// If the polygon has no edges, then the graph's IsFullPolygonPredicate is
// called to determine whether the output polygon should be empty (containing
// no points) or full (containing all points).  This predicate can be
// specified as part of the S2Builder input geometry.

// Specifies that a polygon should be constructed using the given options,
// and that any labels attached to the input edges should be returned in
// "label_set_ids" and "label_set_lexicion".
//
// The labels associated with the edge "polygon.loop(i).vertex({j, j+1})"
// can be retrieved as follows:
//
//   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
class S2PolygonLayer(
    private val polygon: S2Polygon,
    private val labelSetIds: MutableList<MutableList<LabelSetId>>? = null,
    private val labelSetLexicon: IdSetLexicon? = null,
    private val options: Options = Options()
) : Layer() {

    init {
        requireEQ(labelSetIds == null, labelSetLexicon == null)

        if (options.validate) {
            polygon.debugOverride = S2Debug.DISABLE
        }
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

        // If true, calls FindValidationError() on the output polygon.  If any
        // error is found, it will be returned by S2Builder::Build().
        //
        // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
        // order to turn off the default error checking in debug builds.
        //
        // DEFAULT: false
        var validate: Boolean = false
    )

    // Layer interface:

    override fun graphOptions(): GraphOptions {
        // Prevent degenerate edges and sibling edge pairs.  There should not be any
        // duplicate edges if the input is valid, but if there are then we keep them
        // since this tends to produce more comprehensible errors.
        return GraphOptions(options.edgeType, DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.DISCARD)
    }

    override fun build(g: Graph, error: S2Error) {
        labelSetIds?.clear()

        // It's tricky to compute the edge labels for S2Polygons because the
        // S2Polygon::Init methods can reorder and/or invert the loops.  We handle
        // this by remembering the original vector index of each loop and whether or
        // not the loop contained S2::Origin().  By comparing this with the final
        // S2Polygon loops we can fix up the edge labels appropriately.
        val loopMap = mutableMapOf<S2Loop, Pair<Int, Boolean>>()
        if (g.numEdges == 0) {
            // The polygon is either full or empty.
            if (g.isFullPolygon(error)) {
                polygon.init(S2Loop(S2Loop.kFull))
            } else {
                polygon.initNested(mutableListOf())
            }
        } else if (g.options.edgeType == EdgeType.DIRECTED) {
            val edgeLoops = mutableListOf<EdgeLoop>()
            if (!g.getDirectedLoops(Graph.LoopType.SIMPLE, edgeLoops, error)) {
                return;
            }
            val loops = mutableListOf<S2Loop>()
            appendS2Loops(g, edgeLoops, loops)
            appendEdgeLabels(g, edgeLoops)
            edgeLoops.clear() // Release memory
            initLoopMap(loops, loopMap)
            polygon.initOriented(loops)
        } else {
            val components = mutableListOf<UndirectedComponent>()
            if (!g.getUndirectedComponents(Graph.LoopType.SIMPLE, components, error)) {
                return
            }
            // It doesn't really matter which complement of each component we use,
            // since below we normalize all the loops so that they enclose at most
            // half of the sphere (to ensure that the loops can always be nested).
            //
            // The only reason to prefer one over the other is that when there are
            // multiple loops that touch, only one of the two complements matches the
            // structure of the input loops.  GetUndirectedComponents() tries to
            // ensure that this is always complement 0 of each component.
            val loops = mutableListOf<S2Loop>()
            for (component in components) {
                appendS2Loops(g, component[0], loops)
                appendEdgeLabels(g, component[0])
            }
            components.clear()  // Release memory
            initLoopMap(loops, loopMap);
            for (loop in loops) loop.normalize()
            polygon.initNested(loops)
        }
        reorderEdgeLabels(loopMap);


        if (options.validate) {
            polygon.findValidationError(error)
        }
    }


    private fun appendS2Loops(g: Graph, edge_loops: List<EdgeLoop>, loops: MutableList<S2Loop>) {
        val vertices = ArrayList<S2Point>()
        for (edge_loop in edge_loops) {
            vertices.ensureCapacity(edge_loop.size)
            for (edge_id in edge_loop) {
                vertices.add(g.vertex(g.edge(edge_id).first));
            }
            loops.add(S2Loop(vertices.toList(), debugOverride = S2Debug.DISABLE))
            vertices.clear()
        }
    }

    private fun appendEdgeLabels(g: Graph, edgeLoops: List<EdgeLoop>) {
        if (labelSetIds == null || labelSetLexicon == null) return

        val labels = ArrayList<Label>();  // Temporary storage for labels.
        val fetcher = Graph.LabelFetcher(g, options.edgeType)
        for (edgeLoop in edgeLoops) {
            val loopLabelSetIds = ArrayList<LabelSetId>()
            loopLabelSetIds.ensureCapacity(edgeLoop.size)
            for (edgeId in edgeLoop) {
                fetcher.fetch(edgeId, labels)
                loopLabelSetIds.add(labelSetLexicon.add(labels))
            }
            labelSetIds.add(loopLabelSetIds)
        }
    }

    private fun initLoopMap(loops: List<S2Loop>, loop_map: LoopMap) {
        if (labelSetIds == null || labelSetLexicon == null) return
        loops.forEachIndexed { index, loop ->
            loop_map[loop] = Pair(index, loop.containsOrigin())
        }
    }

    private fun reorderEdgeLabels(loop_map: LoopMap) {
        if (labelSetIds == null || labelSetLexicon == null) return

        val newIds = ArrayList<MutableList<LabelSetId>>(labelSetIds.size)
        newIds.assignWith(labelSetIds.size) { mutableListOf() }
        for (i in 0 until polygon.numLoops()) {
            val loop = polygon.loop(i)
            val old = loop_map.getValue(loop)
            val temp = newIds[i]
            newIds[i] = labelSetIds[old.first]
            labelSetIds[old.first] = temp
            if (loop.containsOrigin() != old.second) {
                // S2Loop::Invert() reverses the order of the vertices, which leaves
                // the last edge unchanged.  For example, the loop ABCD (with edges
                // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
                newIds[i].reverse(0, newIds[i].lastIndex)
            }
        }
        labelSetIds.clear()
        labelSetIds.addAll(newIds)
    }

};

