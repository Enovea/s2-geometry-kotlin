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

// A class representing a ShapeEdgeId together with the two endpoints of that
// edge.  It should be passed by reference.
data class ShapeEdge(val id: ShapeEdgeId, val edge: Edge) {
    constructor(shapeId: Int, edgeId: Int, edge: Edge) : this(ShapeEdgeId(shapeId, edgeId), edge)
    constructor(shape: S2Shape, edge_id: Int) : this(shape.id, edge_id, shape.edge(edge_id))

    val v0 = edge.v0
    val v1 = edge.v1

}
