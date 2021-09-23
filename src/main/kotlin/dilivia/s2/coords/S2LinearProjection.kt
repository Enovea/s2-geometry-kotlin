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
package dilivia.s2.coords

import dilivia.math.M_PI
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * S2 Linear projection.
 *
 * This is the fastest transformation, but also produces the least uniform cell sizes. Cell areas vary by a factor of
 * about 5.2, with the largest cells at the center of each face and the smallest cells in the corners.
 *
 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2LinearProjection : S2Projection {

    override fun stToUv(s: Double): Double = 2 * s - 1
    override fun uvToSt(u: Double): Double = 0.5 * (u + 1)

    override val kMinAngleSpan: LengthMetric = LengthMetric(1.0)                     // 1.000
    override val kMaxAngleSpan: LengthMetric = LengthMetric(2.0)                     // 2.000
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(sqrt(2.0 / 3.0))                // 0.816
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.411459345844456965)        // 1.411

    override val kMinEdge: LengthMetric = LengthMetric(2 * sqrt(2.0) / 3.0)       // 0.943
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.440034192955603643)         // 1.440

    override val kMinDiag: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)     // 0.943
    override val kMaxDiag: LengthMetric = LengthMetric(2.0 * sqrt(2.0))           // 2.828
    override val kAvgDiag: LengthMetric = LengthMetric(2.031817866418812674)         // 2.032

    override val kMinArea: AreaMetric = AreaMetric(4.0 / (3.0 * sqrt(3.0)))       // 0.770
    override val kMaxArea: AreaMetric = AreaMetric(4.0)                              // 4.000
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6)                     // 2.094

    override val kMaxEdgeAspect: Double = sqrt(2.0)                                     // 1.414
    override val kMaxDiagAspect: Double = sqrt(3.0)                                     // 1.732

}
