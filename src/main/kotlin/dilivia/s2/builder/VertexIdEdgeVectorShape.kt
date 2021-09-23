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
package dilivia.s2.builder

import dilivia.s2.S2Point
import dilivia.s2.shape.S2Shape

/**
 * An S2Shape used to represent the entire collection of S2Builder input edges.
 *
 * @property vertices The vertices of the S2Builder.
 * @property edges The edges represented by pair of index in the vertices list.
 */
class VertexIdEdgeVectorShape(private val edges: List<Pair<Int, Int>>, private val vertices: List<S2Point>) : S2Shape() {

    fun vertex0(edgeId: Int): S2Point = vertex(edges[edgeId].first)
    fun vertex1(edgeId: Int): S2Point = vertex(edges[edgeId].second)

    // S2Shape interface:

    override val dimension: Int = 1
    override val numEdges: Int = edges.size
    override fun edge(edgeId: Int): dilivia.s2.shape.Edge = dilivia.s2.shape.Edge(vertices[edges[edgeId].first], vertices[edges[edgeId].second])

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int = edges.size
    override fun chain(chainId: Int): Chain = Chain(chainId, 1)
    override fun chainEdge(chainId: Int, offset: Int): dilivia.s2.shape.Edge = edge(chainId)

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    private fun vertex(i: Int): S2Point = vertices[i]

}
