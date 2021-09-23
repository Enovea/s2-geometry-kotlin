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
package dilivia.s2

import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.vectors.R2VectorDouble
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.abs

class S1IntervalUnitTest {

    // "Quadrants" are numbered as follows:
    // quad1 == [0, Pi/2]
    // quad2 == [Pi/2, Pi]
    // quad3 == [-Pi, -Pi/2]
    // quad4 == [-Pi/2, 0]

    val empty = S1Interval.empty()
    val full = S1Interval.full()

    val zero = S1Interval(0.0, 0.0)
    val pi2 = S1Interval(M_PI_2, M_PI_2)
    val pi = S1Interval(M_PI, M_PI)
    val mipi = S1Interval(-M_PI, -M_PI)
    val mipi2 = S1Interval(-M_PI_2, -M_PI_2)

    val quad1 = S1Interval(0.0, M_PI_2)
    val quad2 = S1Interval(M_PI_2, -M_PI)
    val quad3 = S1Interval(M_PI, -M_PI_2)
    val quad4 = S1Interval(-M_PI_2, 0.0)

    val quad12 = S1Interval(0.0, -M_PI)
    val quad23 = S1Interval(M_PI_2, -M_PI_2) // inverted
    val quad34 = S1Interval(-M_PI, 0.0)
    val quad41 = S1Interval(-M_PI_2, M_PI_2)

    val quad123 = S1Interval(0.0, -M_PI_2)
    val quad234 = S1Interval(M_PI_2, 0.0)
    val quad341 = S1Interval(M_PI, M_PI_2)
    val quad412 = S1Interval(-M_PI_2, -M_PI)

    val mid12 = S1Interval(M_PI_2 - 0.01, M_PI_2 + 0.02)
    val mid23 = S1Interval(M_PI - 0.01, -M_PI + 0.02)
    val mid34 = S1Interval(-M_PI_2 - 0.01, -M_PI_2 + 0.02)
    val mid41 = S1Interval(-0.01, 0.02)

    val quad1lo = S1Interval(quad12.lo, mid41.hi)
    val quad2lo = S1Interval(quad23.lo, mid12.hi)
    val quad2hi = S1Interval(mid23.lo, quad12.hi)
    val quad3hi = S1Interval(mid34.lo, quad23.hi)
    val quad12eps = S1Interval(quad12.lo, mid23.hi)
    val quadeps12 = S1Interval(mid41.lo, quad12.hi)
    val quad123eps = S1Interval(quad12.lo, mid34.hi)
    val quadeps123 = S1Interval(mid41.lo, quad23.hi)
    val quad23eps = S1Interval(quad23.lo, mid34.hi)
    val quadeps23 = S1Interval(mid12.lo, quad23.hi)
    val quad412eps = S1Interval(mid34.lo, quad12.hi)

    @Test
    fun constructorsAndAccessors() {
        // Spot-check the constructors and accessors.
        assertThat(quad12.lo).isZero()
        assertThat(quad12.hi).isEqualTo(M_PI)
        assertThat(quad34[0]).isEqualTo(M_PI)
        assertThat(quad34[1]).isEqualTo(0.0)
        assertThat(quad34.bounds).isEqualTo(R2VectorDouble(M_PI, 0.0))
        assertThat(pi.lo).isEqualTo(M_PI)
        assertThat(pi.hi).isEqualTo(M_PI)

        // Check that [-Pi, -Pi] is normalized to [Pi, Pi].
        assertThat(mipi.lo).isEqualTo(M_PI)
        assertThat(mipi.hi).isEqualTo(M_PI)
        assertThat(quad23.lo).isEqualTo(M_PI_2)
        assertThat(quad23.hi).isEqualTo(-M_PI_2)

        // Check that the default S1Interval is identical to Empty().
        val defaultEmpty = S1Interval()
        assertThat(defaultEmpty.isValid).isTrue()
        assertThat(defaultEmpty.isEmpty).isTrue()
        assertThat(empty.lo).isEqualTo(defaultEmpty.lo)
        assertThat(empty.hi).isEqualTo(defaultEmpty.hi)
    }

    @Test
    fun simplePredicates() {
        // is_valid(), is_empty(), is_full(), is_inverted()
        assertThat(zero.isValid && !zero.isEmpty && !zero.isFull).isTrue()
        assertThat(empty.isValid && empty.isEmpty && !empty.isFull).isTrue()
        assertThat(empty.isInverted).isTrue()
        assertThat(full.isValid && !full.isEmpty && full.isFull).isTrue()
        assertThat(!quad12.isEmpty && !quad12.isFull && !quad12.isInverted).isTrue()
        assertThat(!quad23.isEmpty && !quad23.isFull && quad23.isInverted).isTrue()
        assertThat(pi.isValid && !pi.isEmpty && !pi.isInverted).isTrue()
        assertThat(mipi.isValid && !mipi.isEmpty && !mipi.isInverted).isTrue()
    }

    @Test
    fun almostEmptyOrFull() {
        // Test that rounding errors don't cause intervals that are almost empty or
        // full to be considered empty or full.  The following value is the greatest
        // representable value less than Pi.
        val kAlmostPi = M_PI - 2 * DoubleType.epsilon
        assertThat(S1Interval(-kAlmostPi, M_PI).isFull).isFalse()
        assertThat(S1Interval(-M_PI, kAlmostPi).isFull).isFalse()
        assertThat(S1Interval(M_PI, -kAlmostPi).isEmpty).isFalse()
        assertThat(S1Interval(kAlmostPi, -M_PI).isEmpty).isFalse()
    }

    @Test
    fun center() {
        assertThat(quad12.center).isEqualTo(M_PI_2)
        assertThat(S1Interval(3.1, 2.9).center).isEqualTo(3.0 - M_PI)
        assertThat(S1Interval(-2.9, -3.1).center).isEqualTo(M_PI - 3.0)
        assertThat(S1Interval(2.1, -2.1).center).isEqualTo(M_PI)
        assertThat(pi.center).isEqualTo(M_PI)
        assertThat(mipi.center).isEqualTo(M_PI)
        assertThat(abs(quad23.center)).isEqualTo(M_PI)
        assertThat(quad123.center).isEqualTo(0.75 * M_PI)
    }

    @Test
    fun length() {
        assertThat(quad12.length).isEqualTo(M_PI)
        assertThat(pi.length).isEqualTo(0.0)
        assertThat(mipi.length).isEqualTo(0.0)
        assertThat(quad123.length).isEqualTo(1.5 * M_PI)
        assertThat(abs(quad23.length)).isEqualTo(M_PI)
        assertThat(full.length).isEqualTo(2 * M_PI)
        assertThat(empty.length < 0).isTrue()
    }

    @Test
    fun complement() {
        assertThat(empty.complement.isFull).isTrue()
        assertThat(full.complement.isEmpty).isTrue()
        assertThat(pi.complement.isFull).isTrue()
        assertThat(mipi.complement.isFull).isTrue()
        assertThat(zero.complement.isFull).isTrue()
        assertThat(quad12.complement.approxEquals(quad34)).isTrue()
        assertThat(quad34.complement.approxEquals(quad12)).isTrue()
        assertThat(quad123.complement.approxEquals(quad4)).isTrue()
    }

    @Test
    fun contains() {
        // Contains(double), InteriorContains(double)
        assertThat(!empty.contains(0) && !empty.contains(M_PI) && !empty.contains(-M_PI)).isTrue()
        assertThat(!empty.interiorContains(M_PI) && !empty.interiorContains(-M_PI)).isTrue()
        assertThat(full.contains(0) && full.contains(M_PI) && full.contains(-M_PI)).isTrue()
        assertThat(full.interiorContains(M_PI) && full.interiorContains(-M_PI)).isTrue()
        assertThat(quad12.contains(0) && quad12.contains(M_PI) && quad12.contains(-M_PI)).isTrue()
        assertThat(quad12.interiorContains(M_PI_2) && !quad12.interiorContains(0)).isTrue()
        assertThat(!quad12.interiorContains(M_PI) && !quad12.interiorContains(-M_PI)).isTrue()
        assertThat(quad23.contains(M_PI_2) && quad23.contains(-M_PI_2)).isTrue()
        assertThat(quad23.contains(M_PI) && quad23.contains(-M_PI)).isTrue()
        assertThat(!quad23.contains(0)).isTrue()
        assertThat(!quad23.interiorContains(M_PI_2) && !quad23.interiorContains(-M_PI_2)).isTrue()
        assertThat(quad23.interiorContains(M_PI) && quad23.interiorContains(-M_PI)).isTrue()
        assertThat(!quad23.interiorContains(0)).isTrue()
        assertThat(pi.contains(M_PI) && pi.contains(-M_PI) && !pi.contains(0)).isTrue()
        assertThat(!pi.interiorContains(M_PI) && !pi.interiorContains(-M_PI)).isTrue()
        assertThat(mipi.contains(M_PI) && mipi.contains(-M_PI) && !mipi.contains(0)).isTrue()
        assertThat(!mipi.interiorContains(M_PI) && !mipi.interiorContains(-M_PI)).isTrue()
        assertThat(zero.contains(0) && !zero.interiorContains(0)).isTrue()
    }

    private fun testIntervalOps(
        x: S1Interval, y: S1Interval,
        expectedRelation: String, expectedUnion: S1Interval, expectedIntersection: S1Interval
    ) {
        // Test all of the interval operations on the given pair of intervals.
        // "expected_relation" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.
        assertThat(x.contains(y)).isEqualTo(expectedRelation[0] == 'T')
        assertThat(x.interiorContains(y)).isEqualTo(expectedRelation[1] == 'T')
        assertThat(x.intersects(y)).isEqualTo(expectedRelation[2] == 'T')
        assertThat(x.interiorIntersects(y)).isEqualTo(expectedRelation[3] == 'T')

        // bounds() returns a const reference to a member variable, so we need to
        // make a copy when invoking it on a temporary object.
        assertThat(expectedUnion).isEqualTo(x.union(y))
        assertThat(expectedIntersection).isEqualTo(x.intersection(y))
        assertThat(x.contains(y)).isEqualTo(x.union(y) === x)
        assertThat(x.intersects(y)).isEqualTo(!x.intersection(y).isEmpty)
        if (y.lo == y.hi) {
            val r = x.clone().addPoint(y.lo)
            assertThat(expectedUnion).isEqualTo(r)
        }
    }

    @Test
    fun intervalOps() {
        // Contains(S1Interval), InteriorContains(S1Interval),
        // Intersects(), InteriorIntersects(), Union(), Intersection()
        testIntervalOps(empty, empty, "TTFF", empty, empty)
        testIntervalOps(empty, full, "FFFF", full, empty)
        testIntervalOps(empty, zero, "FFFF", zero, empty)
        testIntervalOps(empty, pi, "FFFF", pi, empty)
        testIntervalOps(empty, mipi, "FFFF", mipi, empty)

        testIntervalOps(full, empty, "TTFF", full, empty)
        testIntervalOps(full, full, "TTTT", full, full)
        testIntervalOps(full, zero, "TTTT", full, zero)
        testIntervalOps(full, pi, "TTTT", full, pi)
        testIntervalOps(full, mipi, "TTTT", full, mipi)
        testIntervalOps(full, quad12, "TTTT", full, quad12)
        testIntervalOps(full, quad23, "TTTT", full, quad23)

        testIntervalOps(zero, empty, "TTFF", zero, empty)
        testIntervalOps(zero, full, "FFTF", full, zero)
        testIntervalOps(zero, zero, "TFTF", zero, zero)
        testIntervalOps(zero, pi, "FFFF", S1Interval(0.0, M_PI), empty)
        testIntervalOps(zero, pi2, "FFFF", quad1, empty)
        testIntervalOps(zero, mipi, "FFFF", quad12, empty)
        testIntervalOps(zero, mipi2, "FFFF", quad4, empty)
        testIntervalOps(zero, quad12, "FFTF", quad12, zero)
        testIntervalOps(zero, quad23, "FFFF", quad123, empty)

        testIntervalOps(pi2, empty, "TTFF", pi2, empty)
        testIntervalOps(pi2, full, "FFTF", full, pi2)
        testIntervalOps(pi2, zero, "FFFF", quad1, empty)
        testIntervalOps(pi2, pi, "FFFF", S1Interval(M_PI_2, M_PI), empty)
        testIntervalOps(pi2, pi2, "TFTF", pi2, pi2)
        testIntervalOps(pi2, mipi, "FFFF", quad2, empty)
        testIntervalOps(pi2, mipi2, "FFFF", quad23, empty)
        testIntervalOps(pi2, quad12, "FFTF", quad12, pi2)
        testIntervalOps(pi2, quad23, "FFTF", quad23, pi2)

        testIntervalOps(pi, empty, "TTFF", pi, empty)
        testIntervalOps(pi, full, "FFTF", full, pi)
        testIntervalOps(pi, zero, "FFFF", S1Interval(M_PI, 0.0), empty)
        testIntervalOps(pi, pi, "TFTF", pi, pi)
        testIntervalOps(pi, pi2, "FFFF", S1Interval(M_PI_2, M_PI), empty)
        testIntervalOps(pi, mipi, "TFTF", pi, pi)
        testIntervalOps(pi, mipi2, "FFFF", quad3, empty)
        testIntervalOps(pi, quad12, "FFTF", S1Interval(0.0, M_PI), pi)
        testIntervalOps(pi, quad23, "FFTF", quad23, pi)

        testIntervalOps(mipi, empty, "TTFF", mipi, empty)
        testIntervalOps(mipi, full, "FFTF", full, mipi)
        testIntervalOps(mipi, zero, "FFFF", quad34, empty)
        testIntervalOps(mipi, pi, "TFTF", mipi, mipi)
        testIntervalOps(mipi, pi2, "FFFF", quad2, empty)
        testIntervalOps(mipi, mipi, "TFTF", mipi, mipi)
        testIntervalOps(mipi, mipi2, "FFFF", S1Interval(-M_PI, -M_PI_2), empty)
        testIntervalOps(mipi, quad12, "FFTF", quad12, mipi)
        testIntervalOps(mipi, quad23, "FFTF", quad23, mipi)

        testIntervalOps(quad12, empty, "TTFF", quad12, empty)
        testIntervalOps(quad12, full, "FFTT", full, quad12)
        testIntervalOps(quad12, zero, "TFTF", quad12, zero)
        testIntervalOps(quad12, pi, "TFTF", quad12, pi)
        testIntervalOps(quad12, mipi, "TFTF", quad12, mipi)
        testIntervalOps(quad12, quad12, "TFTT", quad12, quad12)
        testIntervalOps(quad12, quad23, "FFTT", quad123, quad2)
        testIntervalOps(quad12, quad34, "FFTF", full, quad12)

        testIntervalOps(quad23, empty, "TTFF", quad23, empty)
        testIntervalOps(quad23, full, "FFTT", full, quad23)
        testIntervalOps(quad23, zero, "FFFF", quad234, empty)
        testIntervalOps(quad23, pi, "TTTT", quad23, pi)
        testIntervalOps(quad23, mipi, "TTTT", quad23, mipi)
        testIntervalOps(quad23, quad12, "FFTT", quad123, quad2)
        testIntervalOps(quad23, quad23, "TFTT", quad23, quad23)
        testIntervalOps(quad23, quad34, "FFTT", quad234, S1Interval(-M_PI, -M_PI_2))

        testIntervalOps(quad1, quad23, "FFTF", quad123, S1Interval(M_PI_2, M_PI_2))
        testIntervalOps(quad2, quad3, "FFTF", quad23, mipi)
        testIntervalOps(quad3, quad2, "FFTF", quad23, pi)
        testIntervalOps(quad2, pi, "TFTF", quad2, pi)
        testIntervalOps(quad2, mipi, "TFTF", quad2, mipi)
        testIntervalOps(quad3, pi, "TFTF", quad3, pi)
        testIntervalOps(quad3, mipi, "TFTF", quad3, mipi)

        testIntervalOps(quad12, mid12, "TTTT", quad12, mid12)
        testIntervalOps(mid12, quad12, "FFTT", quad12, mid12)

        testIntervalOps(quad12, mid23, "FFTT", quad12eps, quad2hi)
        testIntervalOps(mid23, quad12, "FFTT", quad12eps, quad2hi)

        // This test checks that the union of two disjoint intervals is the smallest
        // interval that contains both of them.  Note that the center of "mid34"
        // slightly CCW of -Pi/2 so that there is no ambiguity about the result.
        testIntervalOps(quad12, mid34, "FFFF", quad412eps, empty)
        testIntervalOps(mid34, quad12, "FFFF", quad412eps, empty)

        testIntervalOps(quad12, mid41, "FFTT", quadeps12, quad1lo)
        testIntervalOps(mid41, quad12, "FFTT", quadeps12, quad1lo)

        testIntervalOps(quad23, mid12, "FFTT", quadeps23, quad2lo)
        testIntervalOps(mid12, quad23, "FFTT", quadeps23, quad2lo)
        testIntervalOps(quad23, mid23, "TTTT", quad23, mid23)
        testIntervalOps(mid23, quad23, "FFTT", quad23, mid23)
        testIntervalOps(quad23, mid34, "FFTT", quad23eps, quad3hi)
        testIntervalOps(mid34, quad23, "FFTT", quad23eps, quad3hi)
        testIntervalOps(quad23, mid41, "FFFF", quadeps123, empty)
        testIntervalOps(mid41, quad23, "FFFF", quadeps123, empty)
    }

    @Test
    fun addPoint() {
        assertThat(empty.clone().addPoint(0)).isEqualTo(zero)
        assertThat(empty.clone().addPoint(M_PI)).isEqualTo(pi)
        assertThat(empty.clone().addPoint(-M_PI)).isEqualTo(mipi)
        assertThat(empty.clone().addPoint(M_PI).addPoint(-M_PI)).isEqualTo(pi)
        assertThat(empty.clone().addPoint(-M_PI).addPoint(M_PI)).isEqualTo(mipi)
        assertThat(empty.clone().addPoint(mid12.lo).addPoint(mid12.hi)).isEqualTo(mid12)
        assertThat(empty.clone().addPoint(mid23.lo).addPoint(mid23.hi)).isEqualTo(mid23)
        assertThat(quad1.clone().addPoint(-0.9 * M_PI).addPoint(-M_PI_2)).isEqualTo(quad123)
        assertThat(full.clone().addPoint(0).isFull).isTrue()
        assertThat(full.clone().addPoint(M_PI).isFull).isTrue()
        assertThat(full.clone().addPoint(-M_PI).isFull).isTrue()
    }

    @Test
    fun project() {
        var r = S1Interval(-M_PI, -M_PI)
        assertThat(r.project(-M_PI)).isEqualTo(M_PI)
        assertThat(r.project(0)).isEqualTo(M_PI)
        r = S1Interval(0.0, M_PI)
        assertThat(r.project(0.1)).isEqualTo(0.1)
        assertThat(r.project(-M_PI_2 + 1e-15)).isEqualTo(0.0)
        assertThat(r.project(-M_PI_2 - 1e-15)).isEqualTo(M_PI)
        r = S1Interval(M_PI - 0.1, -M_PI + 0.1)
        assertThat(r.project(M_PI)).isEqualTo(M_PI)
        assertThat(r.project(1e-15)).isEqualTo(M_PI - 0.1)
        assertThat(r.project(-1e-15)).isEqualTo(-M_PI + 0.1)
        assertThat(full.project(0)).isEqualTo(0.0)
        assertThat(full.project(M_PI)).isEqualTo(M_PI)
        assertThat(full.project(-M_PI)).isEqualTo(M_PI)
    }

    @Test
    fun fromPointPair() {
        assertThat(S1Interval.fromPointPair(-M_PI, M_PI)).isEqualTo(pi)
        assertThat(S1Interval.fromPointPair(M_PI, -M_PI)).isEqualTo(pi)
        assertThat(S1Interval.fromPointPair(mid34.hi, mid34.lo)).isEqualTo(mid34)
        assertThat(S1Interval.fromPointPair(mid23.lo, mid23.hi)).isEqualTo(mid23)
    }


    @Test
    fun expanded() {
        assertThat(empty.expanded(1)).isEqualTo(empty)
        assertThat(full.expanded(1)).isEqualTo(full)
        assertThat(zero.expanded(1)).isEqualTo(S1Interval(-1.0, 1.0))
        assertThat(mipi.expanded(0.01)).isEqualTo(S1Interval(M_PI - 0.01, -M_PI + 0.01))
        assertThat(pi.expanded(27)).isEqualTo(full)
        assertThat(pi.expanded(M_PI_2)).isEqualTo(quad23)
        assertThat(pi2.expanded(M_PI_2)).isEqualTo(quad12)
        assertThat(mipi2.expanded(M_PI_2)).isEqualTo(quad34)
        assertThat(empty.expanded(-1)).isEqualTo(empty)
        assertThat(full.expanded(-1)).isEqualTo(full)
        assertThat(quad123.expanded(-27)).isEqualTo(empty)
        assertThat(quad234.expanded(-27)).isEqualTo(empty)
        assertThat(quad123.expanded(-M_PI_2)).isEqualTo(quad2)
        assertThat(quad341.expanded(-M_PI_2)).isEqualTo(quad4)
        assertThat(quad412.expanded(-M_PI_2)).isEqualTo(quad1)
    }

    @Test
    fun approxEquals() {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        val kLo = 4 * DoubleType.epsilon;  // < max_error default
        val kHi = 6 * DoubleType.epsilon;  // > max_error default

        // Empty intervals.
        assertThat(empty.approxEquals(empty)).isTrue()
        assertThat(zero.approxEquals(empty) && empty.approxEquals(zero)).isTrue()
        assertThat(pi.approxEquals(empty) && empty.approxEquals(pi)).isTrue()
        assertThat(mipi.approxEquals(empty) && empty.approxEquals(mipi)).isTrue()
        assertThat(empty.approxEquals(full)).isFalse()
        assertThat(empty.approxEquals(S1Interval(1.0, 1 + 2 * kLo))).isTrue()
        assertThat(empty.approxEquals(S1Interval(1.0, 1 + 2 * kHi))).isFalse()
        assertThat(S1Interval(M_PI - kLo, -M_PI + kLo).approxEquals(empty)).isTrue()

        // Full intervals.
        assertThat(full.approxEquals(full)).isTrue()
        assertThat(full.approxEquals(empty)).isFalse()
        assertThat(full.approxEquals(zero)).isFalse()
        assertThat(full.approxEquals(pi)).isFalse()
        assertThat(full.approxEquals(S1Interval(kLo, -kLo))).isTrue()
        assertThat(full.approxEquals(S1Interval(2 * kHi, 0.0))).isFalse()
        assertThat(S1Interval(-M_PI + kLo, M_PI - kLo).approxEquals(full)).isTrue()
        assertThat(S1Interval(-M_PI, M_PI - 2 * kHi).approxEquals(full)).isFalse()

        // Singleton intervals.
        assertThat(pi.approxEquals(pi) && mipi.approxEquals(pi)).isTrue()
        assertThat(pi.approxEquals(S1Interval(M_PI - kLo, M_PI - kLo))).isTrue()
        assertThat(pi.approxEquals(S1Interval(M_PI - kHi, M_PI - kHi))).isFalse()
        assertThat(pi.approxEquals(S1Interval(M_PI - kLo, -M_PI + kLo))).isTrue()
        assertThat(pi.approxEquals(S1Interval(M_PI - kHi, -M_PI))).isFalse()
        assertThat(zero.approxEquals(pi)).isFalse()
        assertThat(pi.union(mid12).union(zero).approxEquals(quad12)).isTrue()
        assertThat(quad2.intersection(quad3).approxEquals(pi)).isTrue()
        assertThat(quad3.intersection(quad2).approxEquals(pi)).isTrue()

        // Intervals whose corresponding endpoints are nearly the same but where the
        // endpoints are in opposite order (i.e., inverted intervals).
        assertThat(S1Interval(0.0, kLo).approxEquals(S1Interval(kLo, 0.0))).isFalse()
        assertThat(
            S1Interval(M_PI - 0.5 * kLo, -M_PI + 0.5 * kLo).approxEquals(
                S1Interval(
                    -M_PI + 0.5 * kLo,
                    M_PI - 0.5 * kLo
                )
            )
        ).isFalse()

        // Other intervals.
        assertThat(S1Interval(1 - kLo, 2 + kLo).approxEquals(S1Interval(1, 2))).isTrue()
        assertThat(S1Interval(1 + kLo, 2 - kLo).approxEquals(S1Interval(1, 2))).isTrue()
        assertThat(S1Interval(2 - kLo, 1 + kLo).approxEquals(S1Interval(2, 1))).isTrue()
        assertThat(S1Interval(2 + kLo, 1 - kLo).approxEquals(S1Interval(2, 1))).isTrue()
        assertThat(S1Interval(1 - kHi, 2 + kLo).approxEquals(S1Interval(1, 2))).isFalse()
        assertThat(S1Interval(1 + kHi, 2 - kLo).approxEquals(S1Interval(1, 2))).isFalse()
        assertThat(S1Interval(2 - kHi, 1 + kLo).approxEquals(S1Interval(2, 1))).isFalse()
        assertThat(S1Interval(2 + kHi, 1 - kLo).approxEquals(S1Interval(2, 1))).isFalse()
        assertThat(S1Interval(1 - kLo, 2 + kHi).approxEquals(S1Interval(1, 2))).isFalse()
        assertThat(S1Interval(1 + kLo, 2 - kHi).approxEquals(S1Interval(1, 2))).isFalse()
        assertThat(S1Interval(2 - kLo, 1 + kHi).approxEquals(S1Interval(2, 1))).isFalse()
        assertThat(S1Interval(2 + kLo, 1 - kHi).approxEquals(S1Interval(2, 1))).isFalse()
    }

    @Test
    fun getDirectedHausdorffDistance() {
        assertThat(empty.getDirectedHausdorffDistance(empty)).isEqualTo(0.0)
        assertThat(empty.getDirectedHausdorffDistance(mid12)).isEqualTo(0.0)
        assertThat(mid12.getDirectedHausdorffDistance(empty)).isEqualTo(M_PI)

        assertThat(quad12.getDirectedHausdorffDistance(quad123)).isEqualTo(0.0)
        val interval = S1Interval(3.0, -3.0)  // an interval whose complement center is 0.
        assertThat(S1Interval(-0.1, 0.2).getDirectedHausdorffDistance(interval)).isEqualTo(3.0)
        assertThat(S1Interval(0.1, 0.2).getDirectedHausdorffDistance(interval)).isEqualTo(3.0 - 0.1)
        assertThat(S1Interval(-0.2, -0.1).getDirectedHausdorffDistance(interval)).isEqualTo(3.0 - 0.1)
    }

}
