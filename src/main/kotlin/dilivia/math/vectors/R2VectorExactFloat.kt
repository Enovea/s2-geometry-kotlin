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
 * 2d-vector of exact floats.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class R2VectorExactFloat constructor(coords: Array<BigDecimal>) : AbstractRBigDecimalVector<R2VectorExactFloat>(coords), R2Vector<BigDecimal, R2VectorExactFloat> {

    init {
        requireArgument({ coords.size == 2 }, { "Must have exactly 2 coordinates" })
    }

    override val numberType: NumberType<BigDecimal> = ExactFloatType

    constructor(x: BigDecimal = BigDecimal.ZERO, y: BigDecimal = BigDecimal.ZERO) : this(arrayOf(x, y))

    constructor(x: Double = 0.0, y: Double = 0.0) : this(ExactFloatType.cast(x), ExactFloatType.cast(y))

    constructor(x: String, y: String) : this(BigDecimal(x, ExactFloatType.mathContext), BigDecimal(y, ExactFloatType.mathContext))

    constructor(x: Int, y: Int) : this(x.toDouble(), y.toDouble())

    constructor(coord: List<BigDecimal>) : this(coord.toTypedArray())

    constructor(coord: DoubleArray) : this(coord.map(ExactFloatType::cast).toTypedArray())

    override fun newInstance(coords: Array<BigDecimal>): R2VectorExactFloat = R2VectorExactFloat(coords)

    override fun newInstance(init: (Int) -> BigDecimal): R2VectorExactFloat = newInstance(Array(size, init))

}


operator fun BigDecimal.plus(other: R2VectorExactFloat): R2VectorExactFloat = other + this
operator fun BigDecimal.times(other: R2VectorExactFloat): R2VectorExactFloat = other * this
