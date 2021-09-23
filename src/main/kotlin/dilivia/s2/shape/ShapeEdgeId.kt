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

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32, int32> except that it has named fields.
// It should be passed and returned by value.
data class ShapeEdgeId(val shapeId: Int = -1, val edgeId: Int = -1): Comparable<ShapeEdgeId> {

    override fun compareTo(other: ShapeEdgeId): Int {
        val shapeIdComparison = shapeId.compareTo(other.shapeId)
        return if (shapeIdComparison != 0) shapeIdComparison else edgeId.compareTo(other.edgeId)
    }

}
