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
 * S2 quadratic projection.
 *
 * This is an approximation of the tangent projection that is much faster and produces cells that are almost as uniform
 * in size. It is about 3 times faster than the tangent projection for converting cell ids to points or vice versa.
 * Cell areas vary by a maximum ratio of about 2.1.

 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2QuadraticProjection : S2Projection {

    override fun stToUv(s: Double): Double = if (s >= 0.5)
            (1.0 / 3.0) * (4 * s * s - 1)
    else
            (1.0 / 3.0) * (1 - 4 * (1 - s) * (1 - s))

    override fun uvToSt(u: Double): Double = if (u >= 0)
        0.5 * sqrt(1 + 3 * u)
    else
        1 - 0.5 * sqrt(1 - 3 * u)

    override val kMinAngleSpan: LengthMetric = LengthMetric(4.0 / 3.0)                  // 1.333
    override val kMaxAngleSpan: LengthMetric = LengthMetric(1.704897179199218452)       // 1.705
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                   // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)       // 0.943
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.434523672886099389)           // 1.435

    override val kMinEdge: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)        // 0.943
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.459213746386106062)            // 1.459

    override val kMinDiag: LengthMetric = LengthMetric(8.0 * sqrt(2.0) / 9.0)        // 1.257
    override val kMaxDiag: LengthMetric = LengthMetric(2.438654594434021032)            // 2.439
    override val kAvgDiag: LengthMetric = LengthMetric(2.060422738998471683)            // 2.060

    override val kMinArea: AreaMetric = AreaMetric(8.0 * sqrt(2.0) / 9.0)            // 1.257
    override val kMaxArea: AreaMetric = AreaMetric(2.635799256963161491)                // 2.636
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6)                        // 2.094

    override val kMaxEdgeAspect: Double = 1.442615274452682920                                // 1.443
    override val kMaxDiagAspect: Double = sqrt(3.0)                                        // 1.732

}
