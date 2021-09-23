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

/**
 * 2d-vector of Doubles.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class R2VectorDouble constructor(coords: DoubleArray) : AbstractRDoubleVector<R2VectorDouble>(coords), R2Vector<Double, R2VectorDouble> {

    init {
        requireArgument({ coords.size == 2 }, { "Must have exactly 2 coordinates" })
    }

    constructor(x: Double = 0.0, y: Double = 0.0) : this(doubleArrayOf(x, y))

    constructor(x: Int, y: Int) : this(x.toDouble(), y.toDouble())

    constructor(coord: List<Double>) : this(coord.toDoubleArray())

    override fun newInstance(coords: Array<Double>): R2VectorDouble = R2VectorDouble(coords.toDoubleArray())

    override fun newInstance(init: (Int) -> Double): R2VectorDouble = newInstance(Array(size, init))

    override fun dotProd(other: RVector<Double, R2VectorDouble>): Double {
        return coords[0] * other[0] + coords[1] * other[1]
    }

    override fun minus(other: R2VectorDouble): R2VectorDouble {
        return R2VectorDouble(x = coords[0] - other.coords[0], y = coords[1] - other.coords[1])
    }

    override fun ortho(): R2VectorDouble {
        return R2VectorDouble(-y, x)
    }

}

typealias R2Point = R2VectorDouble

operator fun Double.plus(other: R2VectorDouble): R2VectorDouble = other + this
operator fun Double.times(other: R2VectorDouble): R2VectorDouble = other * this
