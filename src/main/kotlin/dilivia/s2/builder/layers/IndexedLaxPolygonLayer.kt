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
package dilivia.s2.builder.layers

import dilivia.s2.S2Error
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.shape.S2LaxPolygonShape

// Like LaxPolygonLayer, but adds the polygon to a MutableS2ShapeIndex (if the
// polygon is non-empty).
class IndexedLaxPolygonLayer(
    val index: MutableS2ShapeIndex,
    val options: LaxPolygonLayer.Options = LaxPolygonLayer.Options()
) : Layer() {

    val polygon = S2LaxPolygonShape()
    val layer = LaxPolygonLayer(polygon, options = options)

    override fun graphOptions(): GraphOptions = layer.graphOptions()

    override fun build(g: Graph, error: S2Error) {
        layer.build(g, error)
        if (error.isOk() && !polygon.isEmpty()) {
            index.add(polygon)
        }
    }
}
