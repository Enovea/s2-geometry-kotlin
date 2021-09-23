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
package dilivia.s2.index.point

import dilivia.math.matrix.Matrix3x3Double
import dilivia.s2.Fractal
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap

// Generates the vertices of a fractal whose convex hull approximately
// matches the given cap.
class FractalPointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val fractal = Fractal()
        fractal.setLevelForApproxMaxEdges(numPoints)
        fractal.setDimension(1.5)
        val loop = fractal.makeLoop(Matrix3x3Double.fromCols(S2Random.randomFrameAt(indexCap.center)), indexCap.radius())
        for (i in 0 until loop.numVertices) {
            index.add(loop.vertex(i), i)
        }
    }
}
