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
package dilivia.math.vectors

import dilivia.PreConditions
import dilivia.math.DoubleType
import dilivia.math.NumberType
import org.apache.commons.math3.util.FastMath.abs

/**
 * Abstract base class for vector of doubles implementations.
 *
 * @param V The implementation class of the vector.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
abstract class AbstractRDoubleVector<V : AbstractRDoubleVector<V>>(val coords: DoubleArray) : RVector<Double, V> {

    override val numberType: NumberType<Double> = DoubleType

    override val size: Int
        get() = this.coords.size

    override operator fun get(index: Int): Double {
        PreConditions.requireBetween(index, 0, this.size) {
            check_message_component_idx_out_of_range.format(
                index,
                size - 1
            )
        }
        return this.coords[index]
    }

    override operator fun set(index: Int, value: Double) {
        PreConditions.requireBetween(index, 0, this.size) {
            check_message_component_idx_out_of_range.format(
                index,
                size - 1
            )
        }
        this.coords[index] = value
    }

    override fun normalize(): V {
        val norm = norm()
        PreConditions.checkState({ norm != 0.0 }, { "|$this| = 0" })
        for (i in 0 .. coords.lastIndex) {
            coords[i] /= norm
        }
        @Suppress("UNCHECKED_CAST")
        return this as V
    }

    override fun equals(other: Any?): Boolean {

        return this::class.isInstance(other) &&
                @Suppress("UNCHECKED_CAST")
                (this.size == (other as V).size) &&
                (0 until size).all { i ->
                    val ci = this[i]
                    val coi = other[i]
                    ci == coi || (abs(ci) == 0.0 && abs(coi) == 0.0) || (ci.isNaN() && coi.isNaN())
                }
    }

    /**
     * Calculates hashcode based on stored coordinates. Since we want +0.0 and
     * -0.0 to be treated the same, we ignore the sign of the coordinates.
     */
    override fun hashCode(): Int {
        var value: Long = 17
        coords.forEach { c ->
            value += 37 * value + java.lang.Double.doubleToLongBits(abs(c))
        }
        return (value xor (value ushr 32)).toInt()
    }

    override fun toString(): String {
        return coords.joinToString(prefix = "(", postfix = ")")
    }

}
