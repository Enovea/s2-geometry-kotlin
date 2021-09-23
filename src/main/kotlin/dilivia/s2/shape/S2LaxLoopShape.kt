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

import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireLT
import dilivia.s2.S2Point
import dilivia.s2.region.S2Loop
import org.apache.commons.math3.util.FastMath.min

//
// This file defines various S2Shape types representing loops:
//
// S2LaxLoopShape
//   - like S2Loop::Shape but allows duplicate vertices & edges, more compact
//     representation, and faster to initialize.
//
// S2LaxClosedPolylineShape
//   - like S2LaxLoopShape, but defines a loop that does not have an interior
//     (a closed polyline).
//
// S2VertexIdLaxLoopShape
//   - like S2LaxLoopShape, but vertices are specified as indices into an
//     existing vertex array.


// S2LaxLoopShape represents a closed loop of edges surrounding an interior
// region.  It is similar to S2Loop::Shape except that this class allows
// duplicate vertices and edges.  Loops may have any number of vertices,
// including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
// consisting of a single point.)
//
// Note that S2LaxLoopShape is faster to initialize and more compact than
// S2Loop::Shape, but does not support the same operations as S2Loop.
open class S2LaxLoopShape : S2Shape {

    // For clients that have many small loops, we save some memory by
    // representing the vertices as an array rather than using std::vector.
    private lateinit var vertices: Array<S2Point>

    // Constructs an empty loop.
    constructor() {
        vertices = emptyArray()
    }

    // Constructs an S2LaxLoopShape with the given vertices.
    constructor(vertices: List<S2Point>) {
        init(vertices)
    }

    // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
    constructor(loop: S2Loop) {
        init(loop)
    }

    // Initializes an S2LaxLoopShape with the given vertices.
    fun init(vertices: List<S2Point>) {
        this.vertices = vertices.toTypedArray()
    }

    // Initializes an S2LaxLoopShape from the given S2Loop, by copying its data.
    //
    // REQUIRES: !loop->is_full()
    //           [Use S2LaxPolygonShape if you need to represent a full loop.]
    fun init(loop: S2Loop) {
        requireArgument({ !loop.isFull() }, { "Full loops not supported; use S2LaxPolygonShape" })
        if (loop.isEmpty()) {
            vertices = emptyArray()
        } else {
            vertices = loop.verticesSpan().points.toTypedArray()
        }
    }

    fun numVertices(): Int = vertices.size
    fun vertex(i: Int): S2Point = vertices[i]

    override val numEdges: Int
        get() = numVertices()

    override fun edge(edgeId: Int): Edge {
        requireLT(edgeId, numEdges)
        var e = edgeId + 1
        if (e == numVertices()) e = 0
        return Edge(vertices[edgeId], vertices[e])
    }

    override val dimension: Int = 2

    override fun getReferencePoint(): ReferencePoint = S2Shape.getReferencePoint(this)

    override val numChains: Int
        get() = min(1, numVertices())

    override fun chain(chain_id: Int): Chain = Chain(0, numVertices())

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireEQ(chainId, 0)
        requireLT(offset, numEdges)
        val k = if(offset + 1 == numVertices()) 0 else offset + 1
        return Edge(vertices[offset], vertices[k])
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

    override fun toString(): String {
        return "${this.javaClass.simpleName}(id=$id, vertices=${vertices.contentToString()})"
    }


}


