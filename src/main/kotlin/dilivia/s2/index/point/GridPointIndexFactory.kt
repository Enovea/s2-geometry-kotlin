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
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap
import org.apache.commons.math3.util.FastMath.ceil
import org.apache.commons.math3.util.FastMath.sqrt
import org.apache.commons.math3.util.FastMath.tan

// Generates vertices on a square grid that includes the entire given cap.
class GridPointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val sqrtNumPoints = ceil(sqrt(numPoints.toDouble())).toInt()
        val frame = S2Random.randomFrameAt(indexCap.center)
        val radius = indexCap.radius().radians
        val spacing = 2 * radius / sqrtNumPoints;
        for (i in 0 until sqrtNumPoints) {
            for (j in 0 until sqrtNumPoints) {
                val point = S2Point(tan((i + 0.5) * spacing - radius), tan((j + 0.5) * spacing - radius), 1.0)
                index.add(S2PointUtil.fromFrame(Matrix3x3Double.fromCols(frame), point.normalize()), i * sqrtNumPoints + j)
            }
        }
    }
}
