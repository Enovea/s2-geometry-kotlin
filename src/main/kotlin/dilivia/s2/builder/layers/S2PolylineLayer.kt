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
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.Label
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.region.S2Polyline
import mu.KotlinLogging


// A layer type that assembles edges (directed or undirected) into an
// S2Polyline.  Returns an error if the edges cannot be assembled into a
// single unbroken polyline.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of directed input edges
// whenever possible, so that if the input is a polyline and it is not
// modified by S2Builder, then the output will be the same polyline (even if
// the polyline backtracks on itself or forms a loop).  With undirected edges,
// there are no such guarantees; for example, even if the input consists of a
// single undirected edge, then either directed edge may be returned.
//
// S2PolylineLayer does not support options such as discarding sibling pairs
// or merging duplicate edges because these options can split the polyline
// into several pieces.  Use S2PolylineVectorLayer if you need these features.
class S2PolylineLayer(
    private val polyline: S2Polyline,
    private val labelSetIds: S2PointVectorLabelSetIds? = null,
    private val labelSetLexicon: IdSetLexicon? = null,
    private val options: Options = Options()
) : Layer() {

    init {
        requireEQ(labelSetIds == null, labelSetLexicon == null)
    }

    data class Options(

        // Indicates whether the input edges provided to S2Builder are directed or
        // undirected.  Directed edges should be used whenever possible to avoid
        // ambiguity.
        //
        // DEFAULT: S2Builder::EdgeType::DIRECTED
        var edgeType: EdgeType = EdgeType.DIRECTED,

        // If true, calls FindValidationError() on the output polyline.  If any
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
        fun getValidate(): Boolean = validate

        fun setValidate(validate: Boolean) {
            this.validate = validate
            debug = S2Debug.DISABLE
        }
    }


    // Layer interface:

    override fun graphOptions(): GraphOptions {
        // Remove edges that collapse to a single vertex, but keep duplicate and
        // sibling edges, since merging duplicates or discarding siblings can make
        // it impossible to assemble the edges into a single polyline.
        return GraphOptions(options.edgeType, DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.KEEP)
    }

    override fun build(g: Graph, error: S2Error) {
        logger.trace { """
                |
                |Build layer
                |--------------
                |options = ${g.options}
                |vertices: ${g.vertices.map { S2LatLng.fromPoint(it) }}
                |edges: ${g.edges}
                |--------------
            """.trimMargin() }
        if (g.numEdges == 0) {
            polyline.init(mutableListOf())
            return
        }
        val edge_polylines = g.getPolylines(Graph.PolylineType.WALK)
        logger.trace { "Walk Edge polylines = $edge_polylines" }

        if (edge_polylines.size != 1) {
            error.init(S2Error.BUILDER_EDGES_DO_NOT_FORM_POLYLINE, "Input edges cannot be assembled into polyline")
            return
        }
        val edge_polyline = edge_polylines[0]
        val vertices = ArrayList<S2Point>(edge_polyline.size)  // Temporary storage for vertices.
        vertices.add(g.vertex(g.edge(edge_polyline[0]).first))
        for (e in edge_polyline) {
            vertices.add(g.vertex(g.edge(e).second))
        }
        logger.trace { "vertices = $vertices" }
        if (labelSetIds != null && labelSetLexicon != null) {
            val fetcher = Graph.LabelFetcher(g, options.edgeType)
            val labels = mutableListOf<Label>();  // Temporary storage for labels.
            if (labelSetIds is ArrayList) {
                labelSetIds.ensureCapacity(edge_polyline.size)
            }
            for (e in edge_polyline) {
                fetcher.fetch(e, labels)
                labelSetIds.add(labelSetLexicon.add(labels))
            }
        }
        polyline.init(vertices, options.debug)
        if (options.getValidate()) {
            polyline.findValidationError(error)
        }
    }

    companion object {
        val logger = KotlinLogging.logger(S2PolylineLayer::class.java.name)
    }
}


