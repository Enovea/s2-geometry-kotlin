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
import dilivia.s2.region.S2Polyline

// Like S2PolylineVectorLayer, but adds the points to a MutableS2ShapeIndex
class IndexedS2PolylineVectorLayer(
    private val index: MutableS2ShapeIndex,
    options: S2PolylineVectorLayer.Options = S2PolylineVectorLayer.Options()
) : Layer() {

    private val polylines: MutableList<S2Polyline> = mutableListOf()
    private val layer: S2PolylineVectorLayer = S2PolylineVectorLayer(polylines, options = options)

    override fun graphOptions(): GraphOptions = layer.graphOptions()

    override fun build(g: Graph, error: S2Error) {
        layer.build(g, error)
        if (error.isOk()) {
            polylines.forEach { polyline -> index.add(S2Polyline.Shape(polyline = polyline)) }
        }
    }

}
