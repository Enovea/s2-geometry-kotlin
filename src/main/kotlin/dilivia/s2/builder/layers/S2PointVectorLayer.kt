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


typealias S2PointVectorLabelSetIds = MutableList<LabelSetId>

// A layer type that collects degenerate edges as points.
// This layer expects all edges to be degenerate. In case of finding
// non-degenerate edges it sets S2Error but it still generates the
// output with degenerate edges.
class S2PointVectorLayer(
    private val points: MutableList<S2Point>,
    private val labelSetIds: S2PointVectorLabelSetIds? = null,
    private val labelSetLexicon: IdSetLexicon? = null,
    private val options: Options = Options()
) : Layer() {

  data class Options(var duplicateEdges: DuplicateEdges = DuplicateEdges.MERGE)

  // Layer interface:

    override fun graphOptions(): GraphOptions = GraphOptions(
            edgeType = EdgeType.DIRECTED,
            degenerateEdges = DegenerateEdges.KEEP,
            duplicateEdges = options.duplicateEdges,
            siblingPairs = SiblingPairs.KEEP
    )

    override fun build(g: Graph, error: S2Error) {
        val fetcher = Graph.LabelFetcher(g, EdgeType.DIRECTED)

        val labels = mutableListOf<Label>()  // Temporary storage for labels.
        for (edge_id in g.edges.indices) {
            val edge = g.edge(edge_id)
            if (edge.first != edge.second) {
                error.init(S2Error.INVALID_ARGUMENT, "Found non-degenerate edges");
                continue;
            }
            points.add(g.vertex(edge.first));
            if (labelSetIds != null && labelSetLexicon != null) {
                fetcher.fetch(edge_id, labels)
                val set_id = labelSetLexicon.add(labels)
                labelSetIds.add(set_id)
            }
        }
    }

}

