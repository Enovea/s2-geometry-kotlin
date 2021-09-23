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
import dilivia.PreConditions.requireLT
import dilivia.s2.S2Point
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.TypeTags.kLaxPolylineTypeTag
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min

// S2LaxPolylineShape represents a polyline.  It is similar to
// S2Polyline::Shape except that duplicate vertices are allowed, and the
// representation is slightly more compact.
//
// Polylines may have any number of vertices, but note that polylines with
// fewer than 2 vertices do not define any edges.  (To create a polyline
// consisting of a single degenerate edge, either repeat the same vertex twice
// or use S2LaxClosedPolylineShape defined in s2_lax_loop_shape.h.)
// Constructs an S2LaxPolylineShape with the given vertices.
class S2LaxPolylineShape(val vertices: List<S2Point>) : S2Shape() {

    // Constructs an empty polyline.
    constructor() : this(emptyList()) {}

    // Constructs an S2LaxPolylineShape from the given S2Polyline, by copying
    // its data.
    constructor(polyline: S2Polyline) : this(polyline.vertices())

    // Initializes an S2LaxPolylineShape with the given vertices.
    init {
        if (vertices.size == 1) {
            logger.warn { "S2LaxPolylineShape with one vertex has no edges" }
        }
    }

    fun numVertices(): Int = vertices.size
    fun vertex(i: Int): S2Point = vertices[i]

    // S2Shape interface:
    override val numEdges: Int = max(0, numVertices() - 1)
    override fun edge(edgeId: Int): Edge {
        requireLT(edgeId, numEdges)
        return Edge(vertex(edgeId), vertex(edgeId + 1))
    }

    override val dimension: Int = 1
    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)
    override val numChains: Int = min(1, numEdges)
    override fun chain(chain_id: Int): Chain = Chain(0, numEdges)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireEQ(chainId, 0)
        requireLT(offset, numEdges)
        return Edge(vertex(offset), vertex(offset + 1))
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

    override val typeTag: TypeTag = kLaxPolylineTypeTag

    companion object {
        val logger = KotlinLogging.logger(S2LaxPolylineShape::class.java.name)
    }
}


