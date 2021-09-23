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

import dilivia.math.vectors.R2VectorDouble
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.EnumSource


@DisplayName("R1Interval")
internal class R1IntervalUnitTest {

    @DisplayName("Accessors")
    @Test
    fun testAccessors() {
        assertThat(unit.lo).isEqualTo(0.0)
        assertThat(unit.hi).isEqualTo(1.0)
        assertThat(negunit.lo).isEqualTo(-1.0)
        assertThat(negunit.hi).isEqualTo(0.0)
    }

    @DisplayName("Mutations")
    @Test
    fun testMutations() {
        val ten = R1Interval(0, 0)
        ten.hi = 10.0
        assertThat(R1Interval(0, 10)).isEqualTo(ten)
        ten[0] = -10.0
        assertThat(R1Interval(-10, 10)).isEqualTo(ten)
        ten[1] = 0.0
        assertThat(R2VectorDouble(-10.0, 0.0)).isEqualTo(ten.bounds)
        ten.bounds = R2VectorDouble(0, 10)
        assertThat(R1Interval(0, 10)).isEqualTo(ten)
    }

    @DisplayName("Is Empty ? / Is not empty ?")
    @Test
    fun testIsEmpty() {
        assertThat(unit.isNotEmpty).isTrue()
        assertThat(unit.isEmpty).isFalse()
        assertThat(half.isNotEmpty).isTrue()
        assertThat(half.isEmpty).isFalse()
        assertThat(empty.isEmpty).isTrue()
        assertThat(empty.isNotEmpty).isFalse()
    }

    @DisplayName("Center")
    @Test
    fun testCenter() {
        assertThat(unit.center).isEqualTo(0.5)
        assertThat(half.center).isEqualTo(0.5)
    }

    @DisplayName("Length")
    @Test
    fun testLength() {
        assertThat(negunit.length).isEqualTo(1.0)
        assertThat(half.length).isEqualTo(0.0)
        assertThat(empty.length < 0).isTrue()
    }

    @DisplayName("Equality")
    @Test
    fun testEquality() {
        assertThat(empty == empty).isTrue()
        assertThat(unit == unit).isTrue()
        assertThat(unit != empty).isTrue()
        assertThat(R1Interval(1, 2) != R1Interval(1, 3)).isTrue()
    }

    @DisplayName("Default constructor is empty")
    @Test
    fun testDefaultConstructor() {
        // Check that the default R1Interval is identical to Empty().
        assertThat(defaultEmpty.isEmpty).isTrue()
        assertThat(empty.lo).isEqualTo(defaultEmpty.lo)
        assertThat(empty.hi).isEqualTo(defaultEmpty.hi)
    }

    @DisplayName("Contains")
    @Test
    fun testContains() {
        assertThat(unit.contains(0.5)).isTrue()
        assertThat(unit.interiorContains(0.5)).isTrue()
        assertThat(unit.contains(0.0)).isTrue()
        assertThat(!unit.interiorContains(0.0)).isTrue()
        assertThat(unit.contains(1.0)).isTrue()
        assertThat(!unit.interiorContains(1.0)).isTrue()
    }

    @DisplayName("Add point")
    @Test
    fun testAddPoint() {
        val r: R1Interval = empty.clone()
        r.addPoint(5.0)
        assertThat(r.lo == 5.0 && r.hi == 5.0).isTrue()
        r.addPoint(-1.0)
        assertThat(r.lo == -1.0 && r.hi == 5.0).isTrue()
        r.addPoint(0.0)
        assertThat(r.lo == -1.0 && r.hi == 5.0).isTrue()
    }

    @DisplayName("Project")
    @Test
    fun testProject() {
        assertThat(R1Interval(0.1, 0.4).project(0.3)).isEqualTo(0.3)
        assertThat(R1Interval(0.1, 0.4).project(-7.0)).isEqualTo(0.1)
        assertThat(R1Interval(0.1, 0.4).project(0.6)).isEqualTo(0.4)
    }

    @DisplayName("Interval operations")
    @ParameterizedTest(name = "{index} : {0}")
    @EnumSource(IntervalOpsTest::class)
    fun testIntervalOps(intervalOps: IntervalOpsTest) {
        val x = intervalOps.x
        val y = intervalOps.y
        val expectedRelation = intervalOps.expectedRelation
        assertThat(x.contains(y)).isEqualTo(expectedRelation[0] == 'T')
        assertThat(x.interiorContains(y))
            .withFailMessage("$x.interiorContains($y) != ${expectedRelation[1]}")
            .isEqualTo(expectedRelation[1] == 'T')
        assertThat(x.intersects(y)).isEqualTo(expectedRelation[2] == 'T')
        assertThat(x.interiorIntersects(y)).isEqualTo(expectedRelation[3] == 'T')
        assertThat(x.contains(y)).isEqualTo(x.union(y) == x)
        assertThat(x.intersects(y)).isEqualTo(!x.intersection(y).isEmpty)

        val z = x.clone()
        z.addInterval(y)
        assertThat(x.union(y)).isEqualTo(z)
    }

    enum class IntervalOpsTest(
        val displayName: String,
        val x: R1Interval,
        val y: R1Interval,
        val expectedRelation: String
    ) {
        TEST1(displayName = "Empty x Empty", x = empty, y = empty, expectedRelation = "TTFF"),
        TEST2(displayName = "Empty x Unit", x = empty, y = unit, expectedRelation = "FFFF"),
        TEST3(displayName = "Unit x Half", x = unit, y = half, expectedRelation = "TTTT"),
        TEST4(displayName = "Unit x Unit", x = unit, y = unit, expectedRelation = "TFTT"),
        TEST5(displayName = "Unit x Empty", x = unit, y = empty, expectedRelation = "TTFF"),
        TEST6(displayName = "Unit x NegUnit", x = unit, y = negunit, expectedRelation = "FFTF"),
        TEST7(displayName = "Unit x [0.0;0.5]", x = unit, y = R1Interval(0.0, 0.5), expectedRelation = "TFTT"),
        TEST8(displayName = "Half x [0.0;0.5]", x = half, y = R1Interval(0.0, 0.5), expectedRelation = "FFTF");

        override fun toString(): String = displayName
    }

    @DisplayName("From point pair")
    @Test
    fun testFromPointPair() {
        assertThat(R1Interval.fromPointPair(4.0, 4.0)).isEqualTo(R1Interval(4, 4))
        assertThat(R1Interval.fromPointPair(-1.0, -2.0)).isEqualTo(R1Interval(-2, -1))
        assertThat(R1Interval.fromPointPair(-5.0, 3.0)).isEqualTo(R1Interval(-5, 3))
    }

    @DisplayName("Expanded")
    @Test
    fun testExpanded() {
        assertThat(empty.expanded(0.45)).isEqualTo(empty)
        assertThat(unit.expanded(0.5)).isEqualTo(R1Interval(-0.5, 1.5))
        assertThat(R1Interval(0.5, 0.5)).isEqualTo(unit.expanded(-0.5))
        assertThat(unit.expanded(-0.51).isEmpty).isTrue()
        assertThat(unit.expanded(-0.51).expanded(0.51).isEmpty).isTrue()
    }

    @DisplayName("Expands")
    @Test
    fun testExpands() {
        assertThat(empty.clone().expands(0.45)).isEqualTo(empty)
        assertThat(unit.clone().expands(0.5)).isEqualTo(R1Interval(-0.5, 1.5))
        assertThat(unit.clone().expands(-0.5)).isEqualTo(R1Interval(0.5, 0.5))
        assertThat(unit.clone().expands(-0.51).isEmpty).isTrue()
        assertThat(unit.clone().expands(-0.51).expands(0.51).isEmpty).isTrue()
    }

    @DisplayName("Union")
    @Test
    fun testUnion() {
        assertThat(R1Interval(99, 100).union(empty)).isEqualTo(R1Interval(99, 100))
        assertThat(empty.union(R1Interval(99, 100))).isEqualTo(R1Interval(99, 100))
        assertThat(R1Interval(5, 3).union(R1Interval(0, -2)).isEmpty).isTrue()
        assertThat(R1Interval(0, -2).union(R1Interval(5, 3)).isEmpty).isTrue()
        assertThat(unit.union(unit)).isEqualTo(unit)
        assertThat(unit.union(negunit)).isEqualTo(R1Interval(-1, 1))
        assertThat(negunit.union(unit)).isEqualTo(R1Interval(-1, 1))
        assertThat(half.union(unit) == unit).isTrue()
    }

    @DisplayName("Intersection")
    @Test
    fun testIntersection() {
        assertThat(unit.intersection(half)).isEqualTo(half)
        assertThat(unit.intersection(negunit)).isEqualTo(R1Interval(0, 0))
        assertThat(negunit.intersection(half).isEmpty).isTrue()
        assertThat(unit.intersection(empty).isEmpty).isTrue()
        assertThat(empty.intersection(unit).isEmpty).isTrue()
    }

    @DisplayName("Approx equals")
    @Test
    fun testApproxEquals() {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        val kLo = 4 * DoubleType.epsilon  // < max_error default
        val kHi = 6 * DoubleType.epsilon  // > max_error default

        // Empty intervals.
        val empty = R1Interval.empty()
        assertThat(empty.approxEquals(empty)).isTrue()
        assertThat(R1Interval(0, 0).approxEquals(empty)).isTrue()
        assertThat(empty.approxEquals(R1Interval(0, 0))).isTrue()
        assertThat(R1Interval(1, 1).approxEquals(empty)).isTrue()
        assertThat(empty.approxEquals(R1Interval(1, 1))).isTrue()
        assertThat(empty.approxEquals(R1Interval(0, 1))).isFalse()
        assertThat(empty.approxEquals(R1Interval(1.0, 1 + 2 * kLo))).isTrue()
        assertThat(empty.approxEquals(R1Interval(1.0, 1 + 2 * kHi))).isFalse()

        // Singleton intervals.
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1, 1))).isTrue()
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1 - kLo, 1 - kLo))).isTrue()
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1 + kLo, 1 + kLo))).isTrue()
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1.0 - kHi, 1.0))).isFalse()
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1.0, 1 + kHi))).isFalse()
        assertThat(R1Interval(1, 1).approxEquals(R1Interval(1 - kLo, 1 + kLo))).isTrue()
        assertThat(R1Interval(0, 0).approxEquals(R1Interval(1, 1))).isFalse()

        // Other intervals.
        assertThat(R1Interval(1 - kLo, 2 + kLo).approxEquals(R1Interval(1, 2))).isTrue()
        assertThat(R1Interval(1 + kLo, 2 - kLo).approxEquals(R1Interval(1, 2))).isTrue()
        assertThat(R1Interval(1 - kHi, 2 + kLo).approxEquals(R1Interval(1, 2))).isFalse()
        assertThat(R1Interval(1 + kHi, 2 - kLo).approxEquals(R1Interval(1, 2))).isFalse()
        assertThat(R1Interval(1 - kLo, 2 + kHi).approxEquals(R1Interval(1, 2))).isFalse()
        assertThat(R1Interval(1 + kLo, 2 - kHi).approxEquals(R1Interval(1, 2))).isFalse()
    }

    companion object {
        val unit = R1Interval(0, 1)
        val negunit = R1Interval(-1.0, 0.0)
        val half = R1Interval(0.5, 0.5)
        val empty = R1Interval.empty()
        val defaultEmpty = R1Interval()
    }

}
