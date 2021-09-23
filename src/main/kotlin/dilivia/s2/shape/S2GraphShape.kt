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
package dilivia.s2.shape

import dilivia.PreConditions.requireEQ
import dilivia.s2.builder.graph.Graph


// An S2Shape representing the edges in an S2Builder::Graph.
class S2GraphShape(val g: Graph) : S2Shape() {

    override val numEdges: Int
        get() = g.numEdges

    override fun edge(edgeId: Int): Edge {
        val graphEdge = g.edge(edgeId)
        return Edge(g.vertex(graphEdge.first), g.vertex(graphEdge.second))
    }

    override val dimension: Int = 1

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int
        get() = g.numEdges

    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireEQ(offset, 0)
        return edge(chainId)
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

}
