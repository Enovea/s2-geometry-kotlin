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
import dilivia.s2.S2Point

// S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is reprsented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.
class S2PointVectorShape(id: Int = -1, val points: List<S2Point> = emptyList()) : S2Shape(id) {

  fun numPoints(): Int = points.size

  fun point(i: Int): S2Point = points[i]

  // S2Shape interface:

    override val numEdges: Int
        get() = numPoints()

    override fun edge(edgeId: Int): Edge = Edge(points[edgeId], points[edgeId])

    override val dimension: Int = 0

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int
        get() = numPoints()

    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireEQ(offset, 0)
        return Edge(points[chainId], points[chainId])
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    override val typeTag: TypeTag = TypeTags.kPointVectorTypeTag

}


class S2PointShape(id: Int = -1, val point: S2Point): S2Shape(id) {

    override val numEdges: Int = 1

    override fun edge(edgeId: Int): Edge {
        requireEQ(edgeId, 0)
        return Edge(point, point)
    }

    override val dimension: Int = 0

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int = 1

    override fun chain(chainId: Int): Chain {
        requireEQ(chainId, 0)
        return Chain(chainId, 1)
    }

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireEQ(chainId, 0)
        requireEQ(offset, 0)
        return Edge(point, point)
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    override val typeTag: TypeTag = TypeTags.kPointTypeTag
}
