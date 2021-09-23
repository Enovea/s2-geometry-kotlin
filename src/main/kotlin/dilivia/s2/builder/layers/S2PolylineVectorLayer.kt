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
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.Label
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.region.S2Polyline
import mu.KotlinLogging

typealias PolylineVectorLabelSetIds = MutableList<MutableList<LabelSetId>>

// A layer type that assembles edges (directed or undirected) into multiple
// S2Polylines.  Returns an error if S2Builder found any problem with the
// input edges; this layer type does not generate any errors of its own.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of the input edges whenever
// possible, so that if the input is a polyline and it is not modified by
// S2Builder, then the output will be the same polyline even if the polyline
// forms a loop.  However, note that this is not guaranteed when undirected
// edges are used: for example, if the input consists of a single undirected
// edge, then either directed edge may be returned.
class S2PolylineVectorLayer(
        private val polylines: MutableList<S2Polyline>,
        private val labelSetIds: PolylineVectorLabelSetIds? = null,
        private val labelSetLexicon: IdSetLexicon? = null,
        private val options: Options = Options()
) : Layer() {

    init {
        requireEQ(labelSetIds == null, labelSetLexicon == null)
    }

    data class Options(

            // Indicates whether the input edges provided to S2Builder are directed or
            // undirected.
            //
            // Directed edges should be used whenever possible to avoid ambiguity.
            // The implementation attempts to preserve the structure of directed input
            // edges whenever possible, so that if the input is a vector of disjoint
            // polylines and none of them need to be modified then the output will be
            // the same polylines in the same order.  With undirected edges, there are
            // no such guarantees.
            //
            // DEFAULT: S2Builder::EdgeType::DIRECTED
            var edgeType: EdgeType = EdgeType.DIRECTED,

            // Indicates whether polylines should be "paths" (which don't allow
            // duplicate vertices, except possibly the first and last vertex) or
            // "walks" (which allow duplicate vertices and edges).
            //
            // If your input consists of polylines, and you want to split them into
            // separate pieces whenever they self-intersect or cross each other, then
            // use PolylineType::PATH (and probably use split_crossing_edges()).  If
            // you don't mind if your polylines backtrack or contain loops, then use
            // PolylineType::WALK.
            //
            // DEFAULT: PolylineType::PATH
            var polylineType: Graph.PolylineType = Graph.PolylineType.PATH,

            // Indicates whether duplicate edges in the input should be kept (KEEP) or
            // merged together (MERGE).  Note you can use edge labels to determine
            // which input edges were merged into a given output edge.
            //
            // DEFAULT: DuplicateEdges::KEEP
            var duplicateEdges: DuplicateEdges = DuplicateEdges.KEEP,

            // Indicates whether sibling edge pairs (i.e., pairs consisting of an edge
            // and its reverse edge) should be kept (KEEP) or discarded (DISCARD).
            // For example, if a polyline backtracks on itself, the DISCARD option
            // would cause this section of the polyline to be removed.  Note that this
            // option may cause a single polyline to split into several pieces (e.g.,
            // if a polyline has a "lollipop" shape).
            //
            // REQUIRES: sibling_pairs == { DISCARD, KEEP }
            //           (the CREATE and REQUIRE options are not allowed)
            //
            // DEFAULT: SiblingPairs::KEEP
            private var siblingPairs: SiblingPairs = SiblingPairs.KEEP,

            // If true, calls FindValidationError() on each output polyline.  If any
            // error is found, it will be returned by S2Builder::Build().
            //
            // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
            // order to turn off the default error checking in debug builds.
            //
            // DEFAULT: false
            private var validate: Boolean = false,

            // This method can turn off the automatic validity checks triggered by the
            // --s2debug flag (which is on by default in debug builds).  The main
            // reason to do this is if your code already does its own error checking,
            // or if you need to work with invalid geometry for some reason.
            //
            // In any case, polylines have very few restrictions so they are unlikely
            // to have errors.  Errors include vertices that aren't unit length (which
            // can only happen if they are present in the input data), or adjacent
            // vertices that are at antipodal points on the sphere (unlikely with real
            // data).  The other possible error is adjacent identical vertices, but
            // this can't happen because S2Builder does not generate such polylines.
            //
            // DEFAULT: S2Debug::ALLOW
            var debug: S2Debug = S2Debug.ALLOW
    ) {

        fun getSiblingPairs() = siblingPairs

        fun setSiblingPairs(siblingPairs: SiblingPairs) {
            check(siblingPairs == SiblingPairs.KEEP || siblingPairs == SiblingPairs.DISCARD)
            this.siblingPairs = siblingPairs
        }

        fun getValidate(): Boolean = validate

        fun setValidate(validate: Boolean) {
            this.validate = validate
            debug = S2Debug.DISABLE
        }

    }


    // Layer interface:

    override fun graphOptions(): GraphOptions {
        return GraphOptions(options.edgeType, DegenerateEdges.DISCARD, options.duplicateEdges, options.getSiblingPairs())
    }

    override fun build(g: Graph, error: S2Error) {
        val edgePolylines = g.getPolylines(options.polylineType)
        logger.trace { "build |\ng = ${g.toDebugString()}\nedgePolylines = $edgePolylines" }
        if (polylines is ArrayList) {
            polylines.ensureCapacity(edgePolylines.size)
        }
        if (labelSetIds != null && labelSetIds is ArrayList) {
            labelSetIds.ensureCapacity(edgePolylines.size)
        }
        val vertices = mutableListOf<S2Point>()  // Temporary storage for vertices.
        val labels = mutableListOf<Label>()  // Temporary storage for labels.
        for (edge_polyline in edgePolylines) {
            vertices.add(g.vertex(g.edge(edge_polyline[0]).first))
            for (e in edge_polyline) {
            vertices.add(g.vertex(g.edge(e).second))
        }
            val polyline = S2Polyline(vertices, options.debug)
            vertices.clear()
            if (options.getValidate()) {
                polyline.findValidationError(error)
            }
            polylines.add(polyline)

            if (labelSetIds != null && labelSetLexicon != null) {
                val fetcher = Graph.LabelFetcher(g, options.edgeType)
                val polyline_labels = ArrayList<LabelSetId>()
                polyline_labels.ensureCapacity(edge_polyline.size)
                for (e in edge_polyline) {
                    fetcher.fetch(e, labels)
                    polyline_labels.add(labelSetLexicon.add(labels))
                }
                labelSetIds.add(polyline_labels)
            }

        }
    }

    companion object {
        val logger = KotlinLogging.logger(S2PolylineVectorLayer::class.java.name)
    }
}


