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

import dilivia.collections.equalRange
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId

// A class that maps vertices to their outgoing edge ids.  Example usage:
//   VertexOutMap out(g)
//   for (Graph::EdgeId e : out.edge_ids(v)) { ... }
//   for (const Graph::Edge& edge : out.edges(v)) { ... }
class VertexOutMap(val g: Graph) {

    private val edgeBegins: ArrayList<EdgeId> = ArrayList(g.numVertices + 1)

    init {
        var e: EdgeId = 0
        for (v: VertexId in 0 .. g.numVertices) {
            while (e < g.numEdges && g.edge(e).first < v) ++e
            edgeBegins.add(e)
        }
    }

    fun edges(v: VertexId): VertexOutEdges = g.edges.subList(edgeBegins[v], edgeBegins[v + 1])

    fun edgeIds(v: VertexId): VertexOutEdgeIds = edgeBegins[v] until edgeBegins[v + 1]

    // Return the edges (or edge ids) between a specific pair of vertices.
    fun edges(v0: VertexId, v1: VertexId): VertexOutEdges {
        val range = g.edges.equalRange(edgeBegins[v0], edgeBegins[v0 + 1], Edge(v0, v1))
        return g.edges.subList(range.first, range.second)
    }

    fun edgeIds(v0: VertexId, v1: VertexId): VertexOutEdgeIds {
        val range = g.edges.equalRange(edgeBegins[v0], edgeBegins[v0 + 1], Edge(v0, v1))
        if (range.first == g.edges.size) return IntRange.EMPTY
        return range.first until range.second
    }

    fun degree(v: VertexId): Int {
        val edgeIds = edgeIds(v)
        return edgeIds.last - edgeIds.first + 1
    }

    override fun toString(): String {
        return "VertexOutMap: $edgeBegins"
    }

}
