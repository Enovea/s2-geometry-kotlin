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

import dilivia.math.vectors.R2Point
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.EnumSource

@DisplayName("R2Rect")
class R2RectUnitTest {

    @DisplayName("Empty Rectangles")
    @Test
    fun testEmptyRectangles() {
        // Test basic properties of empty rectangles.
        val empty = R2Rect.empty()
        assertThat(empty.isValid).isTrue()
        assertThat(empty.isEmpty).isTrue()
    }

    @DisplayName("Constructors and Accessors")
    @Test
    fun testConstructorsAndAccessors() {
        // Check various constructors and accessor methods.
        val r = R2Rect(R2Point(0.1, 0.0), R2Point(0.25, 1.0))
        assertThat(r.x.lo).isEqualTo(0.1)
        assertThat(r.x.hi).isEqualTo(0.25)
        assertThat(r.y.lo).isEqualTo(0.0)
        assertThat(r.y.hi).isEqualTo(1.0)

        assertThat(r[0][0]).isEqualTo(0.1)
        assertThat(r[0][1]).isEqualTo(0.25)
        assertThat(r[1][0]).isEqualTo(0.0)
        assertThat(r[1][1]).isEqualTo(1.0)

        assertThat(R1Interval(0.1, 0.25)).isEqualTo(r.x)
        assertThat(R1Interval(0, 1)).isEqualTo(r.y)

        assertThat(R1Interval(0.1, 0.25)).isEqualTo(r[0])
        assertThat(R1Interval(0, 1)).isEqualTo(r[1])

        val r2 = R2Rect()
        assertThat(r2.isEmpty).isTrue()
    }

    @DisplayName("FromCenterSize")
    @Test
    fun testFromCenterSize() {
        // FromCenterSize()
        assertThat(
            R2Rect.fromCenterSize(R2Point(0.3, 0.5), R2Point(0.2, 0.4))
                .approxEquals(R2Rect(R2Point(0.2, 0.3), R2Point(0.4, 0.7)))
        ).isTrue()
        assertThat(
            R2Rect.fromCenterSize(R2Point(1.0, 0.1), R2Point(0, 2))
                .approxEquals(R2Rect(R2Point(1.0, -0.9), R2Point(1.0, 1.1)))
        ).isTrue()
    }

    @DisplayName("FromPoint")
    @Test
    fun testFromPoint() {
        // FromPoint(), FromPointPair()
        val d1 = R2Rect(R2Point(0.1, 0.0), R2Point(0.25, 1.0))
        assertThat(R2Rect.fromPoint(d1.lo)).isEqualTo(R2Rect(d1.lo, d1.lo))
        assertThat(R2Rect.fromPointPair(R2Point(0.15, 0.9), R2Point(0.35, 0.3))).isEqualTo(
            R2Rect(
                R2Point(0.15, 0.3),
                R2Point(0.35, 0.9)
            )
        )
        assertThat(R2Rect.fromPointPair(R2Point(0.83, 0.0), R2Point(0.12, 0.5))).isEqualTo(
            R2Rect(
                R2Point(0.12, 0.0),
                R2Point(0.83, 0.5)
            )
        )
    }

    @DisplayName("SimplePredicates")
    @Test
    fun testSimplePredicates() {
        // GetCenter(), GetVertex(), Contains(R2Point), InteriorContains(R2Point).
        val sw1 = R2Point(0.0, 0.25)
        val ne1 = R2Point(0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)

        assertThat(r1.center).isEqualTo(R2Point(0.25, 0.5))
        assertThat(r1.getVertex(0)).isEqualTo(R2Point(0.0, 0.25))
        assertThat(r1.getVertex(1)).isEqualTo(R2Point(0.5, 0.25))
        assertThat(r1.getVertex(2)).isEqualTo(R2Point(0.5, 0.75))
        assertThat(r1.getVertex(3)).isEqualTo(R2Point(0.0, 0.75))
        assertThat(r1.contains(R2Point(0.2, 0.4))).isTrue()
        assertThat(r1.contains(R2Point(0.2, 0.8))).isFalse()
        assertThat(r1.contains(R2Point(-0.1, 0.4))).isFalse()
        assertThat(r1.contains(R2Point(0.6, 0.1))).isFalse()
        assertThat(r1.contains(sw1)).isTrue()
        assertThat(r1.contains(ne1)).isTrue()
        assertThat(r1.interiorContains(sw1)).isFalse()
        assertThat(r1.interiorContains(ne1)).isFalse()

        // Make sure that GetVertex() returns vertices in CCW order.
        for (k in 0 until 4) {
            val a = r1.getVertex(k - 1)
            val b = r1.getVertex(k)
            val c = r1.getVertex(k + 1)
            assertThat((b - a).ortho().dotProd(c - a) > 0).isTrue()
        }
    }

    companion object {
        val empty = R2Rect.empty()
        val sw1 = R2Point(0.0, 0.25)
        val ne1 = R2Point(0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)
        val r1_mid = R2Rect(R2Point(0.25, 0.5), R2Point(0.25, 0.5))
        val r_sw1 = R2Rect(sw1, sw1)
        val r_ne1 = R2Rect(ne1, ne1)
    }

    @DisplayName("Rectangle operations")
    @ParameterizedTest(name = "{index} : {0}")
    @EnumSource(RectangleOpsTest::class)
    fun testRectangleOps(test: RectangleOpsTest) {
        val x: R2Rect = test.x
        val y: R2Rect = test.y
        val expectedRexion: String = test.expectedRexion
        val expectedUnion: R2Rect = test.expectedUnion
        val expectedIntersection: R2Rect = test.expectedIntersection
        // Test all of the interval operations on the given pair of intervals.
        // "expected_rexion" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.

        assertThat(x.contains(y))
            .withFailMessage("$x.contains($y) != ${expectedRexion[0]}")
            .isEqualTo(expectedRexion[0] == 'T')
        assertThat(x.interiorContains(y)).isEqualTo(expectedRexion[1] == 'T')
        assertThat(x.intersects(y)).isEqualTo(expectedRexion[2] == 'T')
        assertThat(x.interiorIntersects(y)).isEqualTo(expectedRexion[3] == 'T')

        assertThat(x.contains(y)).isEqualTo(x.union(y) == x)
        assertThat(x.intersects(y)).isEqualTo(!x.intersection(y).isEmpty)

        assertThat(x.union(y)).isEqualTo(expectedUnion)
        assertThat(x.intersection(y)).isEqualTo(expectedIntersection)

        var r = x.clone()
        r = r.addRect(y)
        assertThat(r).isEqualTo(expectedUnion)
        if (y.size == R2Point(0, 0)) {
            r = x.clone()
            r = r.addPoint(y.lo)
            assertThat(r).isEqualTo(expectedUnion)
        }
    }

    enum class RectangleOpsTest(
        val displayName: String,
        val x: R2Rect, val y: R2Rect,
        val expectedRexion: String, val expectedUnion: R2Rect, val expectedIntersection: R2Rect
    ) {

        TEST1("r1 x r1_mid", r1, r1_mid, "TTTT", r1, r1_mid),
        TEST2("r1 x r_sw1", r1, r_sw1, "TFTF", r1, r_sw1),
        TEST3("r1 x r_ne1", r1, r_ne1, "TFTF", r1, r_ne1),

        TEST4(
            "r1 x R2Rect(R2Point(0.45, 0.1), R2Point(0.75, 0.3))",
            r1,
            R2Rect(R2Point(0.45, 0.1), R2Point(0.75, 0.3)),
        "FFTT",
            R2Rect(R2Point(0.0, 0.1), R2Point(0.75, 0.75)),
            R2Rect(R2Point(0.45, 0.25), R2Point(0.5, 0.3))
        ),
        TEST5(
            "r1 x R2Rect(R2Point(0.5, 0.1), R2Point(0.7, 0.3))",
            r1,
            R2Rect(R2Point(0.5, 0.1), R2Point(0.7, 0.3)),
            "FFTF",
            R2Rect(R2Point(0.0, 0.1), R2Point(0.7, 0.75)),
            R2Rect(R2Point(0.5, 0.25), R2Point(0.5, 0.3))
        ),
        TEST6(
            "r1 x R2Rect(R2Point(0.45, 0.1), R2Point(0.7, 0.25))",
            r1,
            R2Rect(R2Point(0.45, 0.1), R2Point(0.7, 0.25)),
            "FFTF",
            R2Rect(R2Point(0.0, 0.1), R2Point(0.7, 0.75)),
            R2Rect(R2Point(0.45, 0.25), R2Point(0.5, 0.25))
        ),

        TEST7(
            "R2Rect(R2Point(0.1, 0.2), R2Point(0.1, 0.3)) x R2Rect(R2Point(0.15, 0.7), R2Point(0.2, 0.8))",
            R2Rect(R2Point(0.1, 0.2), R2Point(0.1, 0.3)),
            R2Rect(R2Point(0.15, 0.7), R2Point(0.2, 0.8)),
            "FFFF",
            R2Rect(R2Point(0.1, 0.2), R2Point(0.2, 0.8)),
            empty
        ),

        // Check that the intersection of two rectangles that overlap in x but not y
        // is valid, and vice versa.
        TEST8(
            "R2Rect(R2Point(0.1, 0.2), R2Point(0.4, 0.5)) x R2Rect(R2Point(0, 0), R2Point(0.2, 0.1))",
            R2Rect(R2Point(0.1, 0.2), R2Point(0.4, 0.5)),
            R2Rect(R2Point(0, 0), R2Point(0.2, 0.1)),
            "FFFF",
            R2Rect(R2Point(0, 0), R2Point(0.4, 0.5)), empty
        ),
        TEST9(
            "R2Rect(R2Point(0, 0), R2Point(0.1, 0.3)) x R2Rect(R2Point(0.2, 0.1), R2Point(0.3, 0.4))",
            R2Rect(R2Point(0, 0), R2Point(0.1, 0.3)),
            R2Rect(R2Point(0.2, 0.1), R2Point(0.3, 0.4)),
            "FFFF",
            R2Rect(R2Point(0, 0), R2Point(0.3, 0.4)), empty
        );

        override fun toString(): String = displayName
    }

    @DisplayName("AddPoint")
    @Test
    fun testAddPoint() {
        // AddPoint()
        val sw1 = R2Point(0.0, 0.25)
        val ne1 = R2Point(0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)

        val r2 = R2Rect.empty()
        r2.addPoint(R2Point(0.0, 0.25))
        r2.addPoint(R2Point(0.5, 0.25))
        r2.addPoint(R2Point(0.0, 0.75))
        r2.addPoint(R2Point(0.1, 0.4))
        assertThat(r2).isEqualTo(r1)
    }

    @DisplayName("Project")
    @Test
    fun testProject() {
        val r1 = R2Rect(R1Interval(0.0, 0.5), R1Interval(0.25, 0.75))
        assertThat(r1.project(R2Point(-0.01, 0.24))).isEqualTo(R2Point(0.0, 0.25))
        assertThat(r1.project(R2Point(-5.0, 0.48))).isEqualTo(R2Point(0.0, 0.48))
        assertThat(r1.project(R2Point(-5.0, 2.48))).isEqualTo(R2Point(0.0, 0.75))
        assertThat(r1.project(R2Point(0.19, 2.48))).isEqualTo(R2Point(0.19, 0.75))
        assertThat(r1.project(R2Point(6.19, 2.48))).isEqualTo(R2Point(0.5, 0.75))
        assertThat(r1.project(R2Point(6.19, 0.53))).isEqualTo(R2Point(0.5, 0.53))
        assertThat(r1.project(R2Point(6.19, -2.53))).isEqualTo(R2Point(0.5, 0.25))
        assertThat(r1.project(R2Point(0.33, -2.53))).isEqualTo(R2Point(0.33, 0.25))
        assertThat(r1.project(R2Point(0.33, 0.37))).isEqualTo(R2Point(0.33, 0.37))
    }

    @DisplayName("Expanded")
    @Test
    fun testExpanded() {
        // Expanded()
        assertThat(R2Rect.empty().expanded(R2Point(0.1, 0.3)).isEmpty).isTrue()
        assertThat(R2Rect.empty().expanded(R2Point(-0.1, -0.3)).isEmpty).isTrue()
        assertThat(
            R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, 0.3))
                .approxEquals(R2Rect(R2Point(0.1, 0.1), R2Point(0.4, 1.0)))
        ).isTrue()
        assertThat(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(-0.1, 0.3)).isEmpty).isTrue()
        assertThat(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, -0.2)).isEmpty).isTrue()
        assertThat(
            R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, -0.1))
                .approxEquals(R2Rect(R2Point(0.1, 0.5), R2Point(0.4, 0.6)))
        ).isTrue()
        assertThat(
            R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(0.1)
                .approxEquals(R2Rect(R2Point(0.1, 0.3), R2Point(0.4, 0.8)))
        ).isTrue()
    }

}
