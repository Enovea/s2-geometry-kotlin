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
import dilivia.math.ExactFloatType
import dilivia.math.NumberType
import java.math.BigDecimal

/**
 * 3d-vector of exact floats.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class R3VectorExactFloat constructor(coords: Array<BigDecimal>) : AbstractRBigDecimalVector<R3VectorExactFloat>(coords), R3Vector<BigDecimal, R3VectorExactFloat> {

    init {
        requireArgument({ coords.size == 3 }, { "Must have exactly 3 coordinates" })
    }

    override val numberType: NumberType<BigDecimal> = ExactFloatType

    constructor(x: BigDecimal = BigDecimal.ZERO, y: BigDecimal = BigDecimal.ZERO, z: BigDecimal = BigDecimal.ZERO) : this(arrayOf(x, y, z))

    constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0) : this(ExactFloatType.cast(x), ExactFloatType.cast(y), ExactFloatType.cast(z))

    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    constructor(coord: List<BigDecimal>) : this(coord.toTypedArray())

    constructor(coord: DoubleArray) : this(coord.map(ExactFloatType::cast).toTypedArray())

    override fun newInstance(coords: Array<BigDecimal>): R3VectorExactFloat = R3VectorExactFloat(coords)

    override fun newInstance(init: (Int) -> BigDecimal): R3VectorExactFloat = newInstance(Array(size, init))

    /**
     * Gets the dot product of this vector with the other one.
     *
     * @param other A vector.
     * @return The dot product.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    fun dotProd(other: R3VectorExactFloat): BigDecimal =
        numberType.sum((0 until size).map { i: Int -> numberType.times(this[i], other[i]) })

    fun toDouble(): R3VectorDouble = R3VectorDouble(
        x = x.toDouble(),
        y = y.toDouble(),
        z = z.toDouble()
    )


    fun toDecimal128(): R3VectorDecimal128 = R3VectorDecimal128(
        x = x,
        y = y,
        z = z
    )
}
