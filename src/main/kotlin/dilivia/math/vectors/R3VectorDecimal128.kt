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
import dilivia.math.NumberType
import java.math.BigDecimal

/**
 * 3d-vector of decimal 128 numbers.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class R3VectorDecimal128 constructor(coords: Array<BigDecimal>) : AbstractRBigDecimalVector<R3VectorDecimal128>(coords), R3Vector<BigDecimal, R3VectorDecimal128> {

    init {
        requireArgument({ coords.size == 3 }, { "Must have exactly 3 coordinates" })
    }

    override val numberType: NumberType<BigDecimal> = Decimal128Type

    constructor(x: BigDecimal = BigDecimal.ZERO, y: BigDecimal = BigDecimal.ZERO, z: BigDecimal = BigDecimal.ZERO) : this(arrayOf(x, y, z))

    constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0) : this(ExactFloatType.cast(x), ExactFloatType.cast(y), ExactFloatType.cast(z))

    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    constructor(coord: List<BigDecimal>) : this(coord.toTypedArray())

    constructor(coord: DoubleArray) : this(coord.map(ExactFloatType::cast).toTypedArray())

    override fun newInstance(coords: Array<BigDecimal>): R3VectorDecimal128 = R3VectorDecimal128(coords)

    override fun newInstance(init: (Int) -> BigDecimal): R3VectorDecimal128 = newInstance(Array(size, init))

    fun toDouble(): R3VectorDouble = R3VectorDouble(x = x.toDouble(), y = y.toDouble(), z = z.toDouble())

    fun toExactFloat(): R3VectorExactFloat = R3VectorExactFloat(x = x, y = y, z = z)

}
