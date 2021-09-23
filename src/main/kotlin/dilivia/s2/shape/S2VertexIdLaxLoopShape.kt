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

import dilivia.PreConditions
import dilivia.s2.S2Point
import org.apache.commons.math3.util.FastMath.min

// S2VertexIdLaxLoopShape is just like S2LaxLoopShape, except that vertices are
// specified as indices into a vertex array.  This representation can be more
// compact when many loops are arranged in a mesh structure.
class S2VertexIdLaxLoopShape : S2Shape {

    private val vertices: Array<S2Point>
    private val vertexIds: IntArray

    // Constructs an empty loop.
    constructor() : super() {
        vertices = emptyArray()
        vertexIds = intArrayOf()
    }

    // Constructs the shape from the given vertex array and indices.
    // "vertex_ids" is a vector of indices into "vertex_array".
    //
    // ENSURES:  loop->vertex(i) == (*vertex_array)[vertex_ids[i]]
    // REQUIRES: "vertex_array" persists for the lifetime of this object.
    constructor(vertexIds: IntArray, vertices: Array<S2Point>) : super() {
        this.vertexIds = vertexIds
        this.vertices = vertices
    }

    // Returns the number of vertices in the loop.
    fun numVertices(): Int = this.vertexIds.size
    fun vertexId(i: Int): Int = this.vertexIds[i]
    fun vertex(i: Int): S2Point = vertices[vertexId(i)]

    // S2Shape interface:

    override val numEdges: Int
        get() = numVertices()

    override fun edge(edgeId: Int): Edge {
        PreConditions.requireLT(edgeId, numEdges)
        var e1 = edgeId + 1
        if (e1 == numVertices()) e1 = 0
        return Edge(vertex(edgeId), vertex(e1))
    }

    override val dimension: Int = 2

    override fun getReferencePoint(): ReferencePoint {
        // GetReferencePoint interprets a loop with no vertices as "full".
        if (numVertices() == 0) return ReferencePoint(contained = false)
        return getReferencePoint(this)
    }

    override val numChains: Int
        get() = min(1, numVertices())

    override fun chain(chainId: Int): Chain = Chain(0, numVertices())

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        PreConditions.requireEQ(chainId, 0)
        PreConditions.requireLT(offset, numEdges)
        val k = if(offset + 1 == numVertices()) 0 else offset + 1
        return Edge(vertex(offset), vertex(k))
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

}
