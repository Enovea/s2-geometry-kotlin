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

import dilivia.math.ExactFloatType
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.api.Assertions.assertThatThrownBy
import org.junit.jupiter.api.Disabled
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.EnumSource
import java.math.BigDecimal
import kotlin.math.sqrt

@DisplayName("R2VectorExactFloat Tests")
internal class R2VectorExactFloatUnitTest {

    @Disabled()
    @DisplayName("Operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R2VectorExactFloatOperationsTest::class)
    fun operations(test: R2VectorExactFloatOperationsTest) {
        assertThat(test.v1 + test.v2).isEqualTo(test.sum)
        assertThat(test.v1 - test.v2).isEqualTo(test.subtraction)
        assertThat(test.v1 * test.v2).isEqualTo(test.multiplication)
        if (test.division != null)
            assertThat(test.v1 / test.v2).isEqualTo(test.division)
        else
            assertThatThrownBy { test.v1 / test.v2 }.isInstanceOf(ArithmeticException::class.java).hasMessage("Division by zero")
        assertThat(test.v1.crossProd(test.v2)).isEqualTo(test.crossProduct)
        assertThat(test.v1.dotProd(test.v2)).isEqualTo(test.dotProduct)
        assertThat(test.v1 == test.v2).isEqualTo(test.equals)
        assertThat(test.v1.compareTo(test.v2)).isEqualTo(test.compareTo)
        assertThat(test.v1.approxEquals(test.v2, test.margin)).isEqualTo(test.approxEquals)
        assertThat(test.v1.angle(test.v2)).isEqualTo(test.angle)

        val sum = test.v1.clone()
        sum += test.v2
        assertThat(sum).isEqualTo(test.sum)

        val subtraction = test.v1.clone()
        subtraction -= test.v2
        assertThat(subtraction).isEqualTo(test.subtraction)

        val multiplication = test.v1.clone()
        multiplication *= test.v2
        assertThat(multiplication).isEqualTo(test.multiplication)

        if (test.division != null) {
            val division = test.v1.clone()
            division /= test.v2
            assertThat(division).isEqualTo(test.division)
        }
    }

    enum class R2VectorExactFloatOperationsTest(
            val v1: R2VectorExactFloat,
            val v2: R2VectorExactFloat,
            val equals: Boolean,
            val compareTo: Int,
            val margin: Double = 1e-13,
            val approxEquals: Boolean,
            val sum: R2VectorExactFloat,
            val subtraction: R2VectorExactFloat,
            val multiplication: R2VectorExactFloat,
            val division: R2VectorExactFloat?,
            val crossProduct: BigDecimal,
            val dotProduct: BigDecimal,
            val angle: BigDecimal
    ) {

        TEST1(
                v1 = R2VectorExactFloat(3.0, 2.0),
                v2 = R2VectorExactFloat(2.0, 3.0),
                equals = false,
                compareTo = 1,
                approxEquals = false,
                sum = R2VectorExactFloat(5.0, 5.0),
                subtraction = R2VectorExactFloat(1.0, -1.0),
                multiplication = R2VectorExactFloat(6.0, 6.0),
                division = R2VectorExactFloat("1.5", "0.66666666666666666666666666666666666666666666666666666666666666666666666666666667"),
                crossProduct = BigDecimal("5"),
                dotProduct = BigDecimal("12"),
                angle = BigDecimal("0.39479111969976151674009953038958058689517020757570420303537788048206793995648757")
        ),

        TEST2(
                v1 = R2VectorExactFloat(3.0, 0.0),
                v2 = R2VectorExactFloat(2.0, 3.0),
                equals = false,
                compareTo = 1,
                approxEquals = false,
                sum = R2VectorExactFloat(5.0, 3.0),
                subtraction = R2VectorExactFloat(1.0, -3.0),
                multiplication = R2VectorExactFloat(6.0, 0.0),
                division = R2VectorExactFloat(3.0 / 2.0, 0.0),
                crossProduct = ExactFloatType.cast(3.0 * 3.0 - 2.0 * 0.0),
                dotProduct = ExactFloatType.cast(3.0 * 2.0),
                angle = BigDecimal("0.98279372324732906798571061101466601449687745363162855676142508831798807154979604")
        ),

        TEST3(
                v1 = R2VectorExactFloat(1.0, 2.0),
                v2 = R2VectorExactFloat(2.0, 0.0),
                equals = false,
                compareTo = -1,
                approxEquals = false,
                sum = R2VectorExactFloat(3.0, 2.0),
                subtraction = R2VectorExactFloat(-1.0, 2.0),
                multiplication = R2VectorExactFloat(2.0, 0.0),
                division = null,
                crossProduct = ExactFloatType.cast(1.0 * 0.0 - 2.0 * 2.0),
                dotProduct = ExactFloatType.cast(1.0 * 2.0),
                angle = BigDecimal("-1.1071487177940905030170654601785370400700476454014326466765392074337103389773628")
        ),

        TEST4(
                v1 = R2VectorExactFloat(3.0, 2.0),
                v2 = R2VectorExactFloat(3.0, 2.0),
                equals = true,
                compareTo = 0,
                approxEquals = true,
                sum = R2VectorExactFloat(6.0, 4.0),
                subtraction = R2VectorExactFloat(0.0, 0.0),
                multiplication = R2VectorExactFloat(9.0, 4.0),
                division = R2VectorExactFloat(1.0, 1.0),
                crossProduct = BigDecimal.ZERO,
                dotProduct = ExactFloatType.cast(9.0 + 4.0),
                angle = BigDecimal.ZERO
        ),

        TEST5(
                v1 = R2VectorExactFloat(3.0, 2.0),
                v2 = R2VectorExactFloat(3.0 + 1e-14, 2.0 + 1e-14),
                equals = false,
                compareTo = -1,
                approxEquals = true,
                sum = R2VectorExactFloat(
                    "6.000000000000010214051826551440171897411346435546875",
                    "4.000000000000010214051826551440171897411346435546875"),
                subtraction = R2VectorExactFloat("-1.0214051826551440171897411346435546875E-14", "-1.0214051826551440171897411346435546875E-14"),
                multiplication = R2VectorExactFloat(
                    "9.000000000000030642155479654320515692234039306640625",
                    "4.000000000000020428103653102880343794822692871093750"
                ),
                division = R2VectorExactFloat(
                    "0.99999999999999659531605781619820124027571546104723434604727748507151251316954019",
                    "0.99999999999999489297408672430599576497319635184300015506495720418654504912587904"
                ),
                crossProduct = BigDecimal("1.0214051826551440171897411346435546875E-14"),
                dotProduct = BigDecimal("13.000000000000051070259132757200859487056732177734375"),
                angle = BigDecimal("7.8569629435010769586031224918122892528923865728453946006841147910947856327636826E-16")
        ),
    }



    @DisplayName("Scalar operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R2VectorExactFloatScalarOperationsTest::class)
    fun scalarOperations(test: R2VectorExactFloatScalarOperationsTest) {
        assertThat(test.v + test.s).isEqualTo(test.sum)
        assertThat(test.s + test.v).isEqualTo(test.sum)
        assertThat(test.v - test.s).isEqualTo(test.subtraction)
        assertThat(test.v * test.s).isEqualTo(test.multiplication)
        assertThat(test.s * test.v).isEqualTo(test.multiplication)
        if (test.division != null) {
            assertThat(test.v / test.s).isEqualTo(test.division)
        }
        else assertThatThrownBy { test.v / test.s }.isInstanceOf(ArithmeticException::class.java).hasMessage("Division by zero")

        val sum = test.v.clone()
        sum += test.s
        assertThat(sum).isEqualTo(test.sum)

        val subtraction = test.v.clone()
        subtraction -= test.s
        assertThat(subtraction).isEqualTo(test.subtraction)

        val multiplication = test.v.clone()
        multiplication *= test.s
        assertThat(multiplication).isEqualTo(test.multiplication)

        if (test.division != null) {
            val division = test.v.clone()
            division /= test.s
            assertThat(division).isEqualTo(test.division)
        }
    }


    enum class R2VectorExactFloatScalarOperationsTest(
            val v: R2VectorExactFloat,
            val s: BigDecimal,
            val sum: R2VectorExactFloat,
            val subtraction: R2VectorExactFloat,
            val multiplication: R2VectorExactFloat,
            val division: R2VectorExactFloat?,
    ) {

        TEST1(
                v = R2VectorExactFloat(3.0, 2.0),
                s = ExactFloatType.cast(2.5),
                sum = R2VectorExactFloat(5.5, 4.5),
                subtraction = R2VectorExactFloat(0.5, -0.5),
                multiplication = R2VectorExactFloat("7.5", "5.0"),
                division = R2VectorExactFloat("1.2", "0.8"),
        ),

        TEST2(
                v = R2VectorExactFloat(3.0, 2.0),
                s = ExactFloatType.cast(0.0),
                sum = R2VectorExactFloat(3.0, 2.0),
                subtraction = R2VectorExactFloat(3.0, 2.0),
                multiplication = R2VectorExactFloat(0.0, 0.0),
                division = null,
        )

    }


    @Disabled()
    @DisplayName("Unary operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R2VectorExactFloatUnaryOperationsTest::class)
    fun unaryOperations(test: R2VectorExactFloatUnaryOperationsTest) {
        assertThat(test.v.sqrt()).isEqualTo(test.sqrt)
        assertThat(test.v.abs()).isEqualTo(test.abs)
        assertThat(test.v.norm()).isEqualTo(test.norm)
        assertThat(test.v.norm2()).isEqualTo(test.norm2)
        if (test.normalized != null)
            assertThat(test.v.normalized()).isEqualTo(test.normalized)
        else assertThatThrownBy { test.v.normalized() }.isInstanceOf(IllegalStateException::class.java)
        assertThat(-test.v).isEqualTo(R2VectorExactFloat(-test.v.x, -test.v.y))
        assertThat(test.v.ortho()).isEqualTo(test.ortho)

        test.v.x = ExactFloatType.cast(10.0)
        test.v.y = ExactFloatType.cast(10.0)
        assertThat(test.v).isEqualTo(R2VectorExactFloat(10.0, 10.0))
    }

    enum class R2VectorExactFloatUnaryOperationsTest(
            val v: R2VectorExactFloat,
            val sqrt: R2VectorExactFloat,
            val abs: R2VectorExactFloat,
            val norm: BigDecimal,
            val norm2: BigDecimal,
            val normalized: R2VectorExactFloat?,
            val ortho: R2VectorExactFloat
    ) {

        TEST1(
                v = R2VectorExactFloat(3.0, 2.0),
                sqrt = R2VectorExactFloat(
                    "1.7320508075688772935274463415058723669428052538103806280558069794519330169088",
                    "1.414213562373095048801688724209698078569671875376948073176679737990732478462107"),
                abs = R2VectorExactFloat(3.0, 2.0),
                norm = BigDecimal("3.6055512754639892931192212674704959462512965738452462127104530562271669482930104", ExactFloatType.mathContext),
                norm2 = ExactFloatType.cast(3.0 * 3.0 + 2.0 * 2.0),
                normalized = R2VectorExactFloat(
                    "0.83205029433784368302751260018549906451952997857967220293318147451396160345223319",
                    "0.55470019622522912201834173345699937634635331905311480195545431634264106896815546"
                ),
                ortho = R2VectorExactFloat(-2.0, 3.0)
        ),

        TEST2(
                v = R2VectorExactFloat(0.0, 2.0),
                sqrt = R2VectorExactFloat("0", "1.414213562373095048801688724209698078569671875376948073176679737990732478462107"),
                abs = R2VectorExactFloat(0.0, 2.0),
                norm = BigDecimal("2"),
                norm2 = BigDecimal("4"),
                normalized = R2VectorExactFloat(0.0, 2.0 / sqrt(2.0 * 2.0)),
                ortho = R2VectorExactFloat(-2.0, 0.0)
        ),

        TEST3(
                v = R2VectorExactFloat(0.0, 0.0),
                sqrt = R2VectorExactFloat("0", "0"),
                abs = R2VectorExactFloat("0", "0"),
                norm = BigDecimal("0"),
                norm2 = BigDecimal("0"),
                normalized = null,
                ortho = R2VectorExactFloat(-0.0, 0.0)
        )

    }

}
