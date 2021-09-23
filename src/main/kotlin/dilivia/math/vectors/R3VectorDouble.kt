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

import dilivia.PreConditions.requireArgument
import dilivia.math.Decimal128Type
import dilivia.math.ExactFloatType
import org.apache.commons.math3.util.FastMath

/**
 * 3d-vector of Doubles.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class R3VectorDouble constructor(coords: DoubleArray) : AbstractRDoubleVector<R3VectorDouble>(coords), R3Vector<Double, R3VectorDouble> {

    init {
        requireArgument({ coords.size == 3 }, { "Must have exactly 3 coordinates" })
    }

    constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0) : this(doubleArrayOf(x, y, z))

    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    constructor(coord: List<Double>) : this(coord.toDoubleArray())

    constructor(vector: R3VectorExactFloat) : this(
        x = vector.x.toDouble(),
        y = vector.y.toDouble(),
        z = vector.z.toDouble()
    )

    override fun newInstance(coords: Array<Double>): R3VectorDouble = R3VectorDouble(coords.toDoubleArray())

    override fun newInstance(init: (Int) -> Double): R3VectorDouble = newInstance(Array(size, init))

    override fun largestAbsComponent(): Int {
        return if (kotlin.math.abs(coords[0]) > kotlin.math.abs(coords[1])) {
            if (kotlin.math.abs(coords[0]) > kotlin.math.abs(coords[2])) 0 else 2
        } else {
            if (kotlin.math.abs(coords[1]) > kotlin.math.abs(coords[2])) 1 else 2
        }
    }

    override fun plus(other: R3VectorDouble): R3VectorDouble = R3VectorDouble(doubleArrayOf(
        x + other.x,
        y + other.y,
        z + other.z
    ))

    override fun plus(other: Double): R3VectorDouble  = R3VectorDouble(doubleArrayOf(
        x + other,
        y + other,
        z + other
    ))

    override fun plusAssign(other: R3VectorDouble) {
        x += other.x
        y += other.y
        z += other.z
    }

    override fun plusAssign(other: Double) {
        x += other
        y += other
        z += other
    }

    override fun minus(other: R3VectorDouble): R3VectorDouble = R3VectorDouble(doubleArrayOf(
        x - other.x,
        y - other.y,
        z - other.z
    ))

    override fun minus(other: Double): R3VectorDouble  = R3VectorDouble(doubleArrayOf(
        x - other,
        y - other,
        z - other
    ))

    override fun minusAssign(other: R3VectorDouble) {
        x -= other.x
        y -= other.y
        z -= other.z
    }

    override fun minusAssign(other: Double) {
        x -= other
        y -= other
        z -= other
    }

    override fun crossProd(other: R3Vector<Double, R3VectorDouble>): R3VectorDouble {
        val thisX = coords[0]
        val thisY = coords[1]
        val thisZ = coords[2]
        val otherX = other.x
        val otherY = other.y
        val otherZ = other.z
        return R3VectorDouble(
            doubleArrayOf(thisY * otherZ - thisZ * otherY, thisZ * otherX - thisX * otherZ, thisX * otherY - thisY * otherX)
        )
    }

    override fun dotProd(other: RVector<Double, R3VectorDouble>): Double {
        return x * other[0] + y * other[1] + z * other[2]
    }

    fun toDecimal128(): R3VectorDecimal128 = R3VectorDecimal128(
        x = Decimal128Type.cast(x),
        y = Decimal128Type.cast(y),
        z = Decimal128Type.cast(z)
    )

    fun toExactFloat(): R3VectorExactFloat = R3VectorExactFloat(
        x = ExactFloatType.cast(x),
        y = ExactFloatType.cast(y),
        z = ExactFloatType.cast(z)
    )


    override fun equals(other: Any?): Boolean {
        if (other !is R3VectorDouble) return false
        val thisX = coords[0]
        val thisY = coords[1]
        val thisZ = coords[2]
        val otherX = other.coords[0]
        val otherY = other.coords[1]
        val otherZ = other.coords[2]
        return  (thisX == otherX || (FastMath.abs(thisX) == 0.0 && FastMath.abs(otherX) == 0.0) || (thisX.isNaN() && otherX.isNaN())) &&
                (thisY == otherY || (FastMath.abs(thisY) == 0.0 && FastMath.abs(otherY) == 0.0) || (thisY.isNaN() && otherY.isNaN())) &&
                (thisZ == otherZ || (FastMath.abs(thisZ) == 0.0 && FastMath.abs(otherZ) == 0.0) || (thisZ.isNaN() && otherZ.isNaN()))
    }

    override fun hashCode(): Int {
        var value: Long = 17
        coords.forEach { c ->
            value += 37 * value + java.lang.Double.doubleToLongBits(FastMath.abs(c))
        }
        return (value xor (value ushr 32)).toInt()
    }


}

operator fun Double.plus(other: R3VectorDouble): R3VectorDouble = other + this
operator fun Double.times(other: R3VectorDouble): R3VectorDouble = other * this
operator fun Int.times(other: R3VectorDouble): R3VectorDouble = other * this.toDouble()
