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

import dilivia.math.M_1_PI
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import org.apache.commons.math3.util.FastMath.atan
import org.apache.commons.math3.util.FastMath.sqrt
import org.apache.commons.math3.util.FastMath.tan

/**
 * S2 tangente projection
 *
 * Transforming the coordinates via atan() makes the cell sizes more uniform. The areas vary by a maximum ratio of 1.4
 * as opposed to a maximum ratio of 5.2. However, each call to atan() is about as expensive as all of the other
 * calculations combined when converting from points to cell ids, i.e. it reduces performance by a factor of 3.
 *
 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2TanProjection : S2Projection {

    override fun stToUv(s: Double): Double {
        // Unfortunately, tan(M_PI_4) is slightly less than 1.0.  This isn't due to a flaw in the implementation of
        // tan(), it's because the derivative of tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
        // point numbers on either side of the infinite-precision value of pi/4 have tangents that are slightly below
        // and slightly above 1.0 when rounded to the nearest double-precision result.

        val t = tan(M_PI_2 * s - M_PI_4)
        return t + (1.0 / (1L shl 53)) * t
    }

    override fun uvToSt(u: Double): Double {
        val a = atan(u)
        return (2 * M_1_PI) * (a + M_PI_4)
    }

    override val kMinAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                         // 1.571
    override val kMaxAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                         // 1.571
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                         // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(M_PI / (2 * sqrt(2.0)))            // 1.111
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.437318638925160885)                 // 1.437

    override val kMinEdge: LengthMetric = LengthMetric(M_PI / (2 * sqrt(2.0)))             // 1.111
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.461667032546739266)                  // 1.462

    override val kMinDiag: LengthMetric = LengthMetric(M_PI * sqrt(2.0) / 3)               // 1.481
    override val kMaxDiag: LengthMetric = LengthMetric(M_PI * sqrt(2.0 / 3.0))             // 2.565
    override val kAvgDiag: LengthMetric = LengthMetric(2.063623197195635753)                  // 2.064

    override val kMinArea: AreaMetric = AreaMetric((M_PI * M_PI) / (4.0 * sqrt(2.0)))      // 1.745
    override val kMaxArea: AreaMetric = AreaMetric(M_PI * M_PI / 4.0)                         // 2.467
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6)                              // 2.094

    override val kMaxEdgeAspect: Double = sqrt(2.0)                                              // 1.414
    override val kMaxDiagAspect: Double = sqrt(3.0)                                              // 1.732

}
