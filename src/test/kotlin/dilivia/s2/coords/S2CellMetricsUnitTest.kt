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

import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow


// Note: obviously, I could have defined a bundle of metrics like this in the
// S2 class itself rather than just for testing.  However, it's not clear that
// this is useful other than for testing purposes, and I find
//  S2Coords.projection.kMinWidth.GetLevelForMinValue(width) to be slightly more readable than
// than  S2Coords.projection.kWidth.min().GetLevelForMinValue(width).  Also, there is no
// fundamental reason that we need to analyze the minimum, maximum, and average
// values of every metric; it would be perfectly reasonable to just define
// one of these.
data class MetricBundle(val dim: Int, val min: S2CellMetric, val max: S2CellMetric, val avg: S2CellMetric)

class S2CellMetricsTest {

    fun checkMinMaxAvg(bundle: MetricBundle) {
        assertThat(bundle.min.deriv <= bundle.avg.deriv).isTrue()
        assertThat(bundle.avg.deriv <= bundle.max.deriv).isTrue()
    }

    fun checkLessOrEqual(a: MetricBundle, b: MetricBundle) {
        assertThat(a.min.deriv <= b.min.deriv).isTrue()
        assertThat(a.max.deriv <= b.max.deriv).isTrue()
        assertThat(a.avg.deriv <= b.avg.deriv).isTrue()
    }

    @Test
    fun test() {
        val angleSpanBundle = MetricBundle(
            1,
            S2Coords.projection.kMinAngleSpan,
            S2Coords.projection.kMaxAngleSpan,
            S2Coords.projection.kAvgAngleSpan
        )
        val widthBundle =
            MetricBundle(1, S2Coords.projection.kMinWidth, S2Coords.projection.kMaxWidth, S2Coords.projection.kAvgWidth)
        val edgeBundle =
            MetricBundle(1, S2Coords.projection.kMinEdge, S2Coords.projection.kMaxEdge, S2Coords.projection.kAvgEdge)
        val diagBundle =
            MetricBundle(1, S2Coords.projection.kMinDiag, S2Coords.projection.kMaxDiag, S2Coords.projection.kAvgDiag)
        val areaBundle =
            MetricBundle(2, S2Coords.projection.kMinArea, S2Coords.projection.kMaxArea, S2Coords.projection.kAvgArea)

        // First, check that min <= avg <= max for each metric.
        checkMinMaxAvg(angleSpanBundle)
        checkMinMaxAvg(widthBundle)
        checkMinMaxAvg(edgeBundle)
        checkMinMaxAvg(diagBundle)
        checkMinMaxAvg(areaBundle)

        // Check that the maximum aspect ratio of an individual cell is consistent
        // with the global minimums and maximums.
        assertThat(S2Coords.projection.kMaxEdgeAspect >= 1).isTrue()
        assertThat(S2Coords.projection.kMaxEdgeAspect <= S2Coords.projection.kMaxEdge.deriv / S2Coords.projection.kMinEdge.deriv).isTrue()
        assertThat(S2Coords.projection.kMaxDiagAspect >= 1).isTrue()
        assertThat(S2Coords.projection.kMaxDiagAspect <= S2Coords.projection.kMaxDiag.deriv / S2Coords.projection.kMinDiag.deriv).isTrue()

        // Check various conditions that are provable mathematically.
        checkLessOrEqual(widthBundle, angleSpanBundle)
        checkLessOrEqual(widthBundle, edgeBundle)
        checkLessOrEqual(edgeBundle, diagBundle)

        assertThat(S2Coords.projection.kMinArea.deriv >= S2Coords.projection.kMinWidth.deriv * S2Coords.projection.kMinEdge.deriv - 1e-15).isTrue()
        assertThat(S2Coords.projection.kMaxArea.deriv <= S2Coords.projection.kMaxWidth.deriv * S2Coords.projection.kMaxEdge.deriv + 1e-15).isTrue()

        // GetLevelForMaxValue() and friends have built-in assertions, we just need
        // to call these functions to test them.
        //
        // We don't actually check that the metrics are correct here, e.g. that
        // GetMinWidth(10) is a lower bound on the width of cells at level 10.
        // It is easier to check these properties in s2cell_test, since
        // S2Cell has methods to compute the cell vertices, etc.

        for (level in -2..(S2Coords.kMaxCellLevel + 3)) {
            var width = S2Coords.projection.kMinWidth.deriv * 2.0.pow((-level).toDouble())
            if (level >= S2Coords.kMaxCellLevel + 3) width = 0.0
            // Check boundary cases (exactly equal to a threshold value).
            val expectedLevel = max(0, min(S2Coords.kMaxCellLevel, level))
            assertThat(S2Coords.projection.kMinWidth.getLevelForMaxValue(width)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinWidth.getLevelForMinValue(width)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinWidth.getClosestLevel(width)).isEqualTo(expectedLevel)

            // Also check non-boundary cases.
            assertThat(S2Coords.projection.kMinWidth.getLevelForMaxValue(1.2 * width)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinWidth.getLevelForMinValue(0.8 * width)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinWidth.getClosestLevel(1.2 * width)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinWidth.getClosestLevel(0.8 * width)).isEqualTo(expectedLevel)

            // Same thing for area.
            var area = S2Coords.projection.kMinArea.deriv * 4.0.pow((-level).toDouble())
            if (level <= -3) area = 0.0
            assertThat(S2Coords.projection.kMinArea.getLevelForMaxValue(area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getLevelForMinValue(area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getClosestLevel(area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getLevelForMaxValue(1.2 * area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getLevelForMinValue(0.8 * area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getClosestLevel(1.2 * area)).isEqualTo(expectedLevel)
            assertThat(S2Coords.projection.kMinArea.getClosestLevel(0.8 * area)).isEqualTo(expectedLevel)
        }
    }
}
