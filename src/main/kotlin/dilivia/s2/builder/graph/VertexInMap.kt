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
package dilivia.s2.builder.graph

import dilivia.s2.builder.EdgeId

// A class that maps vertices to their incoming edge ids.  Example usage:
//   VertexInMap in(g)
//   for (Graph::EdgeId e : in.edge_ids(v)) { ... }
class VertexInMap(val g: Graph) {

    // sorted vector of all incoming edges (see GetInEdgeIds).
    val inEdgeIds: List<EdgeId> = g.getInEdgeIds()
    private val inEdgeBegins: ArrayList<EdgeId> = ArrayList(g.numVertices + 1)

    init {
        var e: EdgeId = 0
        for (v: VertexId in 0..g.numVertices) {
            while (e < g.numEdges && g.edge(inEdgeIds[e]).second < v) ++e
            inEdgeBegins.add(e)
        }
    }

    fun degree(v: VertexId): Int = edgeIds(v).size

    fun edgeIds(v: VertexId): VertexInEdgeIds = inEdgeIds.subList(inEdgeBegins[v], inEdgeBegins[v + 1])


    override fun toString(): String {
        return "VertexInMap: $inEdgeBegins"
    }
}
