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

import dilivia.PreConditions
import dilivia.math.M_SQRT2
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min

/**
 * Defines a cell metric of the given dimension (1 == length, 2 == area).
 *
 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
open class S2CellMetric
/**
 * @param dim Metric dimension
 * @param deriv deriv value.
 */(
    /** The metric dimension (1 == length, 2 == area). */
    val dim: Int,
    /**
     * The "deriv" value of a metric is a derivative, and must be multiplied by a length or area in
     * (s,t)-space to get a useful value. */
    val deriv: Double
) {

    // Return the value of a metric for cells at the given level. The value is
    // either a length or an area on the unit sphere, depending on the
    // particular metric.
    fun getValue(level: Int): Double = StrictMath.scalb(deriv, -dim * level)

    // Return the level at which the metric has approximately the given
    // value.  For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the
    // level at which the average cell edge length is approximately 0.1.
    // The return value is always a valid level.
    fun getClosestLevel(value: Double): Int = getLevelForMaxValue((if (dim == 1) M_SQRT2 else 2.0) * value)

    // Return the minimum level such that the metric is at most the given
    // value, or S2CellId::kMaxLevel if there is no such level.  For example,
    // S2::kMaxDiag.GetLevelForMaxValue(0.1) returns the minimum level such
    // that all cell diagonal lengths are 0.1 or smaller.  The return value
    // is always a valid level.
    fun getLevelForMaxValue(value: Double): Int {
        if (value <= 0) return S2Coords.kMaxCellLevel

        // This code is equivalent to computing a floating-point "level"
        // value and rounding up.  ilogb() returns the exponent corresponding to a
        // fraction in the range [1,2).
        val exponent = StrictMath.getExponent(value / deriv)
        val level = max(0, min(S2Coords.kMaxCellLevel, -(exponent shr (dim - 1))))
        PreConditions.checkState(
            { level == S2Coords.kMaxCellLevel || getValue(level) <= value },
            { "getLevelForMaxValue($value): deriv = $deriv => level = $level, getValue(level) = ${getValue(level)}" })
        PreConditions.checkState { level == 0 || getValue(level - 1) > value }
        return level
    }

    // Return the maximum level such that the metric is at least the given
    // value, or zero if there is no such level.  For example,
    // S2::kMinWidth.GetLevelForMinValue(0.1) returns the maximum level such
    // that all cells have a minimum width of 0.1 or larger.  The return value
    // is always a valid level.
    fun getLevelForMinValue(value: Double): Int {
        if (value <= 0) return S2Coords.kMaxCellLevel

        // This code is equivalent to computing a floating-point "level"
        // value and rounding down.
        val exponent = StrictMath.getExponent(deriv / value)
        val level = max(0, min(S2Coords.kMaxCellLevel, exponent shr (dim - 1)))
        PreConditions.checkState { level == 0 || getValue(level) >= value }
        PreConditions.checkState { level == S2Coords.kMaxCellLevel || getValue(level + 1) < value }
        return level
    }

}

class LengthMetric(deriv: Double) : S2CellMetric(1, deriv)
class AreaMetric(deriv: Double) : S2CellMetric(2, deriv)
