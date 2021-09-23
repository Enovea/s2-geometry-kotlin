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
package dilivia.s2.index.shape

import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap
import dilivia.s2.shape.S2PointVectorShape

// Generates a cloud of points that approximately fills the given S2Cap.
class PointCloudShapeIndexFactory : ShapeIndexFactory {

    override fun addEdges(index_cap: S2Cap, num_edges: Int, index: MutableS2ShapeIndex) {
        val points = mutableListOf<S2Point>()
        repeat (num_edges) {
            points.add(S2Random.samplePoint(index_cap))
        }
        index.add(S2PointVectorShape(points = points))
    }

}
