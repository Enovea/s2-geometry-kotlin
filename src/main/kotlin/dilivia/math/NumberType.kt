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
package dilivia.math

import ch.obermuhlner.math.big.BigDecimalMath
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.FastMath.pow
import java.math.BigDecimal
import java.math.MathContext


fun ldexp(x: Double, exp: Int) = x * pow(2.0, exp)

interface NumberType<T : Number> : Comparator<T> {
    val epsilon: T
    val roundingEpsilon: T
    fun cast(v: Double): T
    fun plus(v1: T, v2: T): T
    fun plus(v1: T, v2: T, v3: T): T
    fun plus(v: Array<T>): T
    fun minus(v1: T, v2: T): T
    fun unaryMinus(v: T): T
    fun times(v1: T, v2: T): T
    fun div(v1: T, v2: T): T
    fun sum(values: Iterable<T>): T
    fun sqrt(v: T): T
    fun abs(v: T): T
    fun atan2(y: T, x: T): T
    fun one(): T
    fun zero(): T
    fun sign(x: T): Int

}

object DoubleType : NumberType<Double> {

    override val epsilon: Double = pow(2.0, -52)

    override val roundingEpsilon = pow(2.0, -53)

    override fun cast(v: Double): Double = v

    override fun plus(v1: Double, v2: Double): Double = v1 + v2

    override fun plus(v1: Double, v2: Double, v3: Double): Double = v1 + v2 + v3

    override fun plus(v: Array<Double>): Double = v.fold(0.0) { v1, v2 -> v1 + v2 }

    override fun minus(v1: Double, v2: Double): Double = v1 - v2

    override fun unaryMinus(v: Double): Double = if (v != 0.0) -v else v

    override fun times(v1: Double, v2: Double): Double = v1 * v2

    override fun div(v1: Double, v2: Double): Double = v1 / v2

    override fun sum(values: Iterable<Double>): Double = values.sum()

    override fun sqrt(v: Double): Double = FastMath.sqrt(v)

    override fun abs(v: Double): Double = FastMath.abs(v)

    override fun atan2(y: Double, x: Double): Double = FastMath.atan2(y, x)

    override fun one(): Double = 1.0

    override fun zero(): Double = 0.0

    override fun sign(x: Double): Int = FastMath.signum(x).toInt()

    fun epsilonForDigits(digits: Int): Double {
        return if (digits < 64) pow(2.0, -digits) else epsilonForDigits(digits - 63) / (1.toULong() shl 63).toDouble()
    }

    override fun compare(o1: Double, o2: Double): Int = if (o1 == o2) 0 else o1.compareTo(o2)

}

open class BigDecimalType(val mathContext: MathContext) : NumberType<BigDecimal> {

    override val epsilon: BigDecimal =
        if (mathContext.precision > 0) BigDecimal("1e-${mathContext.precision - 1}", mathContext) else BigDecimal.ZERO

    override val roundingEpsilon: BigDecimal =
        if (mathContext.precision > 0) BigDecimal("5e-${mathContext.precision}", mathContext) else BigDecimal.ZERO

    override fun cast(v: Double): BigDecimal = BigDecimal(v, mathContext)

    override fun plus(v1: BigDecimal, v2: BigDecimal): BigDecimal = v1.add(v2, mathContext)

    override fun plus(v1: BigDecimal, v2: BigDecimal, v3: BigDecimal): BigDecimal =
        v1.add(v2, mathContext).add(v3, mathContext)

    override fun plus(v: Array<BigDecimal>): BigDecimal = v.fold(BigDecimal.ZERO) { v1, v2 -> plus(v1, v2) }

    override fun minus(v1: BigDecimal, v2: BigDecimal): BigDecimal = v1.subtract(v2, mathContext)

    override fun unaryMinus(v: BigDecimal): BigDecimal = v.negate(mathContext)

    override fun times(v1: BigDecimal, v2: BigDecimal): BigDecimal = v1.multiply(v2, mathContext)

    override fun div(v1: BigDecimal, v2: BigDecimal): BigDecimal {
        return v1.divide(v2, mathContext)
    }

    override fun sum(values: Iterable<BigDecimal>): BigDecimal =
        values.fold(BigDecimal.ZERO) { acc: BigDecimal, v: BigDecimal -> acc.add(v, mathContext) }

    override fun sqrt(v: BigDecimal): BigDecimal = v.sqrt(mathContext)

    override fun abs(v: BigDecimal): BigDecimal = v.abs(mathContext)

    override fun atan2(y: BigDecimal, x: BigDecimal): BigDecimal = BigDecimalMath.atan2(y, x, mathContext)

    override fun one(): BigDecimal = BigDecimal.ONE

    override fun zero(): BigDecimal = BigDecimal.ZERO

    override fun sign(x: BigDecimal): Int = x.signum()

    override fun compare(o1: BigDecimal, o2: BigDecimal): Int = o1.compareTo(o2)

}

object ExactFloatType : BigDecimalType(MathContext.UNLIMITED)
object Decimal64Type : BigDecimalType(MathContext.DECIMAL64)
object Decimal128Type : BigDecimalType(MathContext.DECIMAL128)

fun Int.toDecimal128() = BigDecimal(this, Decimal128Type.mathContext)
fun Double.toDecimal128() = BigDecimal(this, Decimal128Type.mathContext)
fun Double.toExactFloat() = BigDecimal(this, ExactFloatType.mathContext)
