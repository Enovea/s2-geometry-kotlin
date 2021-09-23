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

import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.EnumSource
import kotlin.Double.Companion.NaN

@DisplayName("R3VectorDouble Tests")
internal class R3VectorDoubleUnitTest {

    @DisplayName("Operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R3VectorDoubleOperationsTest::class)
    fun operations(test: R3VectorDoubleOperationsTest) {
        assertThat(test.v1 + test.v2).isEqualTo(test.sum)
        assertThat(test.v1 - test.v2).isEqualTo(test.subtraction)
        assertThat(test.v1 * test.v2).isEqualTo(test.multiplication)
        assertThat(test.v1 / test.v2).isEqualTo(test.division)
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

        val division = test.v1.clone()
        division /= test.v2
        assertThat(division).isEqualTo(test.division)
    }

    enum class R3VectorDoubleOperationsTest(
            val v1: R3VectorDouble,
            val v2: R3VectorDouble,
            val equals: Boolean,
            val compareTo: Int,
            val margin: Double = 1e-13,
            val approxEquals: Boolean,
            val sum: R3VectorDouble,
            val subtraction: R3VectorDouble,
            val multiplication: R3VectorDouble,
            val division: R3VectorDouble,
            val crossProduct: R3VectorDouble,
            val dotProduct: Double,
            val angle: Double
    ) {

        TEST1(
                v1 = R3VectorDouble(1.0, -1.0, 0.0),
                v2 = R3VectorDouble(1.0, 1.0,0.0),
                equals = false,
                compareTo = -1,
                approxEquals = false,
                sum = R3VectorDouble(2.0, 0.0, 0.0),
                subtraction = R3VectorDouble(0.0, -2.0, 0.0),
                multiplication = R3VectorDouble(1.0, -1.0,  0.0),
                division = R3VectorDouble(1.0, -1.0, NaN),
                crossProduct = R3VectorDouble(-0.0,  0.0, 2.0),
                dotProduct = 0.0,
                angle = 1.5707963267948966
        ),

    }
/*
    @DisplayName("Scalar operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R2VectorDoubleScalarOperationsTest::class)
    fun scalarOperations(test: R2VectorDoubleScalarOperationsTest) {
        assertThat(test.v + test.s).isEqualTo(test.sum)
        assertThat(test.s + test.v).isEqualTo(test.sum)
        assertThat(test.v - test.s).isEqualTo(test.subtraction)
        assertThat(test.v * test.s).isEqualTo(test.multiplication)
        assertThat(test.s * test.v).isEqualTo(test.multiplication)
        assertThat(test.v / test.s).isEqualTo(test.division)

        val sum = test.v.clone()
        sum += test.s
        assertThat(sum).isEqualTo(test.sum)

        val subtraction = test.v.clone()
        subtraction -= test.s
        assertThat(subtraction).isEqualTo(test.subtraction)

        val multiplication = test.v.clone()
        multiplication *= test.s
        assertThat(multiplication).isEqualTo(test.multiplication)

        val division = test.v.clone()
        division /= test.s
        assertThat(division).isEqualTo(test.division)
    }


    enum class R2VectorDoubleScalarOperationsTest(
            val v: R2VectorDouble,
            val s: Double,
            val sum: R2VectorDouble,
            val subtraction: R2VectorDouble,
            val multiplication: R2VectorDouble,
            val division: R2VectorDouble,
    ) {

        TEST1(
                v = R2VectorDouble(3.0, 2.0),
                s = 2.5,
                sum = R2VectorDouble(5.5, 4.5),
                subtraction = R2VectorDouble(0.5, -0.5),
                multiplication = R2VectorDouble(3.0 * 2.5, 2.0 * 2.5),
                division = R2VectorDouble(3.0 / 2.5, 2.0 / 2.5),
        ),

        TEST2(
                v = R2VectorDouble(3.0, 2.0),
                s = 0.0,
                sum = R2VectorDouble(3.0, 2.0),
                subtraction = R2VectorDouble(3.0, 2.0),
                multiplication = R2VectorDouble(0.0, 0.0),
                division = R2VectorDouble(POSITIVE_INFINITY, POSITIVE_INFINITY),
        )

    }

    @DisplayName("Unary operations.")
    @ParameterizedTest(name = "{index}: {argumentsWithNames}")
    @EnumSource(R2VectorDoubleUnaryOperationsTest::class)
    fun unaryOperations(test: R2VectorDoubleUnaryOperationsTest) {
        assertThat(test.v.sqrt()).isEqualTo(test.sqrt)
        assertThat(test.v.abs()).isEqualTo(test.abs)
        assertThat(test.v.norm()).isEqualTo(test.norm)
        assertThat(test.v.norm2()).isEqualTo(test.norm2)
        if (test.normalized != null)
            assertThat(test.v.normalized()).isEqualTo(test.normalized)
        else assertThatThrownBy { test.v.normalized() }.isInstanceOf(IllegalStateException::class.java)
        assertThat(-test.v).isEqualTo(R2VectorDouble(-test.v.x, -test.v.y))
        assertThat(test.v.ortho()).isEqualTo(test.ortho)

        test.v.x = 10.0
        test.v.y = 10.0
        assertThat(test.v).isEqualTo(R2VectorDouble(10.0, 10.0))
    }

    enum class R2VectorDoubleUnaryOperationsTest(
            val v: R2VectorDouble,
            val sqrt: R2VectorDouble,
            val abs: R2VectorDouble,
            val norm: Double,
            val norm2: Double,
            val normalized: R2VectorDouble?,
            val ortho: R2VectorDouble
    ) {

        TEST1(
                v = R2VectorDouble(3.0, 2.0),
                sqrt = R2VectorDouble(sqrt(3.0), sqrt(2.0)),
                abs = R2VectorDouble(3.0, 2.0),
                norm = sqrt(3.0 * 3.0 + 2.0 * 2.0),
                norm2 = 3.0 * 3.0 + 2.0 * 2.0,
                normalized = R2VectorDouble(3.0 / sqrt(3.0 * 3.0 + 2.0 * 2.0), 2.0 / sqrt(3.0 * 3.0 + 2.0 * 2.0)),
                ortho = R2VectorDouble(-2.0, 3.0)
        ),

        TEST2(
                v = R2VectorDouble(0.0, 2.0),
                sqrt = R2VectorDouble(sqrt(0.0), sqrt(2.0)),
                abs = R2VectorDouble(0.0, 2.0),
                norm = sqrt(2.0 * 2.0),
                norm2 = 2.0 * 2.0,
                normalized = R2VectorDouble(0.0, 2.0 / sqrt(2.0 * 2.0)),
                ortho = R2VectorDouble(-2.0, 0.0)
        ),

        TEST3(
                v = R2VectorDouble(0.0, 0.0),
                sqrt = R2VectorDouble(sqrt(0.0), sqrt(0.0)),
                abs = R2VectorDouble(0.0, 0.0),
                norm = 0.0,
                norm2 = 0.0,
                normalized = null,
                ortho = R2VectorDouble(-0.0, 0.0)
        )

    }
*/
}


