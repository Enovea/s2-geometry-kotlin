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

import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.pow

internal class DoubleTypeUnitTest {

    @Test
    fun epsilonForDigitsRecursion() {
        assertThat(DoubleType.epsilonForDigits(0)).isEqualTo(1.0)
        assertThat(DoubleType.epsilonForDigits(24)).isEqualTo(2.0.pow(-24))
        assertThat(DoubleType.epsilonForDigits(53)).isEqualTo(2.0.pow(-53))
        assertThat(DoubleType.epsilonForDigits(64)).isEqualTo(2.0.pow(-64))
        assertThat(DoubleType.epsilonForDigits(106)).isEqualTo(2.0.pow(-106))
        assertThat(DoubleType.epsilonForDigits(113)).isEqualTo(2.0.pow(-113))
    }

    @Test
    fun testRoundingEpsilonVSNumericLimits() {
        // Check that rounding_epsilon<T>() returns the expected value for "float"
        // and "double".  We explicitly do not test "long double" since if this type
        // is implemented using double-double arithmetic then the numeric_limits
        // epsilon() value is completely unrelated to the maximum rounding error.
        assertThat(DoubleType.roundingEpsilon).isEqualTo(0.5 * DoubleType.epsilon)
    }
}

