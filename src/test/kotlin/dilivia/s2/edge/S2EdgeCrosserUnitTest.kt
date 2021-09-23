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
package dilivia.s2.edge

import dilivia.PreConditions
import dilivia.math.vectors.times
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.nextTowards
import kotlin.math.pow

class S2EdgeCrosserUnitTest {

// In non-debug builds, check that default-constructed and/or NaN S2Point
// arguments don't cause crashes, especially on the very first method call
// (since S2CopyingEdgeCrosser checks whether the first vertex of each edge is
// the same as the last vertex of the previous edged when deciding whether or
// not to call Restart).

    fun testCrossingSignInvalid(point: S2Point, expected: Int) {
        val crosser = S2EdgeCrosser(point, point)
        assertThat(crosser.crossingSign(point, point)).isEqualTo(expected)
    }

    fun testEdgeOrVertexCrossingInvalid(point: S2Point, expected: Boolean) {
        val crosser = S2EdgeCrosser(point, point)
        assertThat(crosser.edgeOrVertexCrossing(point, point)).isEqualTo(expected)
    }

    @Test
    fun invalidDefaultPoints() {
        if (!PreConditions.enabled) {
            // Check that default-constructed S2Point arguments don't cause crashes.
            val point = S2Point(0, 0, 0)
            testCrossingSignInvalid(point, 0)
            testEdgeOrVertexCrossingInvalid(point, false)
        }
    }

    @Test
    fun invalidNanPoints() {
        if (!PreConditions.enabled) {
            // Check that NaN S2Point arguments don't cause crashes.
            val nan = Double.NaN
            val point = S2Point(nan, nan, nan)
            testCrossingSignInvalid(point, -1)
            testEdgeOrVertexCrossingInvalid(point, false)
        }
    }

    private fun testCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point, robust: Int, edge_or_vertex: Boolean) {
        var r = robust
        // Modify the expected result if two vertices from different edges match.
        if (a == c || a == d || b == c || b == d) r = 0
        assertThat(S2EdgeCrossings.crossingSign(a, b, c, d)).isEqualTo(r)
        val crosser = S2EdgeCrosser(a, b, c)
        assertThat(crosser.crossingSign(d))         .isEqualTo(r)
        assertThat(crosser.crossingSign(c))         .isEqualTo(r)
        assertThat(crosser.crossingSign(d, c))      .isEqualTo(r)
        assertThat(crosser.crossingSign(c, d))      .isEqualTo(r)

        assertThat(S2EdgeCrossings.edgeOrVertexCrossing(a, b, c, d))   .isEqualTo(edge_or_vertex)
        crosser.restartAt(c)
        assertThat(crosser.edgeOrVertexCrossing(d))            .isEqualTo(edge_or_vertex)
        assertThat(crosser.edgeOrVertexCrossing(c))            .isEqualTo(edge_or_vertex)
        assertThat(crosser.edgeOrVertexCrossing(d, c))         .isEqualTo(edge_or_vertex)
        assertThat(crosser.edgeOrVertexCrossing(c, d))         .isEqualTo(edge_or_vertex)

        // Check that the crosser can be re-used.
        crosser.init(c, d)
        crosser.restartAt(a)
        assertThat(crosser.crossingSign(b)) .isEqualTo(r)
        assertThat(crosser.crossingSign(a)) .isEqualTo(r)
    }

    private fun testCrossings(a: S2Point, b: S2Point, c: S2Point, d: S2Point, robust: Int, edge_or_vertex: Boolean) {
        a.normalize()
        b.normalize()
        c.normalize()
        d.normalize()
        testCrossing(a, b, c, d, robust, edge_or_vertex)
        testCrossing(b, a, c, d, robust, edge_or_vertex)
        testCrossing(a, b, d, c, robust, edge_or_vertex)
        testCrossing(b, a, d, c, robust, edge_or_vertex)
        testCrossing(a, a, c, d, -1, false)
        testCrossing(a, b, c, c, -1, false)
        testCrossing(a, a, c, c, -1, false)
        testCrossing(a, b, a, b, 0, true)
        testCrossing(c, d, a, b, robust, edge_or_vertex != (robust == 0))
    }

    @Test
    fun crossings() {
        // The real tests of edge crossings are in s2{loop,polygon}_test,
        // but we do a few simple tests here.

        // Two regular edges that cross.
        testCrossings(
            S2Point(1, 2, 1),
            S2Point(1.0, -3.0, 0.5),
            S2Point(1.0, -0.5, -3.0),
            S2Point(0.1, 0.5, 3.0),
            1,
            true
        )

        // Two regular edges that intersect antipodal points.
        testCrossings(
            S2Point(1, 2, 1),
            S2Point(1.0, -3.0, 0.5),
            S2Point(-1.0, 0.5, 3.0),
            S2Point(-0.1, -0.5, -3.0),
            -1,
            false
        )

        // Two edges on the same great circle that start at antipodal points.
        testCrossings(S2Point(0, 0, -1), S2Point(0, 1, 0), S2Point(0, 0, 1), S2Point(0, 1, 1), -1, false)

        // Two edges that cross where one vertex is S2::Origin().
        testCrossings(S2Point(1, 0, 0), S2PointUtil.origin(), S2Point(1.0, -0.1, 1.0), S2Point(1.0, 1.0, -0.1), 1, true)

        // Two edges that intersect antipodal points where one vertex is
        // S2::Origin().
        testCrossings(S2Point(1, 0, 0), S2PointUtil.origin(), S2Point(-1.0, 0.1, -1.0), S2Point(-1.0, -1.0, 0.1), -1, false)

        // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
        // and edge CD is further CCW around (2,3,4) than AB.
        testCrossings(S2Point(2, 3, 4), S2Point(-1, 2, 5), S2Point(7, -2, 3), S2Point(2, 3, 4), 0, false)

        // Two edges that barely cross each other near the middle of one edge.  The
        // edge AB is approximately in the x=y plane, while CD is approximately
        // perpendicular to it and ends exactly at the x=y plane.
        testCrossings(
            S2Point(1, 1, 1),
            S2Point(1.0, 1.0.nextTowards(0.0), -1.0),
            S2Point(11, -12, -1),
            S2Point(10, 10, 1),
            1,
            true
        )

        // In this version, the edges are separated by a distance of about 1e-15.
        testCrossings(
            S2Point(1, 1, 1),
            S2Point(1.0, 1.0.nextTowards(2.0), -1.0),
            S2Point(1, -1, 0),
            S2Point(1, 1, 0),
            -1,
            false
        )

        // Two edges that barely cross each other near the end of both edges.  This
        // example cannot be handled using regular double-precision arithmetic due
        // to floating-point underflow.
        testCrossings(
            S2Point(0, 0, 1),
            S2Point(2.0, -1e-323, 1.0),
            S2Point(1, -1, 1),
            S2Point(1e-323, 0.0, 1.0),
            1,
            true
        )

        // In this version, the edges are separated by a distance of about 1e-640.
        testCrossings(
            S2Point(0, 0, 1),
            S2Point(2.0, 1e-323, 1.0),
            S2Point(1, -1, 1),
            S2Point(1e-323, 0.0, 1.0),
            -1,
            false
        )

        // Two edges that barely cross each other near the middle of one edge.
        // Computing the exact determinant of some of the triangles in this test
        // requires more than 2000 bits of precision.
        testCrossings(
            S2Point(1.0, -1e-323, -1e-323),
            S2Point(1e-323, 1.0, 1e-323),
            S2Point(1.0, -1.0, 1e-323),
            S2Point(1, 1, 0),
            1,
            true
        )

        // In this version, the edges are separated by a distance of about 1e-640.
        testCrossings(
            S2Point(1.0, 1e-323, -1e-323),
            S2Point(-1e-323, 1.0, 1e-323),
            S2Point(1.0, -1.0, 1e-323),
            S2Point(1, 1, 0),
            -1,
            false
        )
    }

    @Test
    fun collinearEdgesThatDontTouch() {
        val kIters = 500
        repeat(kIters) {
            val a = S2Random.randomPoint()
            val d = S2Random.randomPoint()
            val b = S2EdgeDistances.interpolate(0.05, a, d)
            val c = S2EdgeDistances.interpolate(0.95, a, d)
            assertThat(S2EdgeCrossings.crossingSign(a, b, c, d)).isLessThanOrEqualTo(0)
            assertThat(S2EdgeCrossings.crossingSign(a, b, c, d)).isLessThanOrEqualTo(0)
            val crosser = S2EdgeCrosser(a, b, c)
            assertThat(crosser.crossingSign(d)).isLessThanOrEqualTo(0)
            assertThat(crosser.crossingSign(c)).isLessThanOrEqualTo(0)
        }
    }

    @Test
    fun coincidentZeroLengthEdgesThatDontTouch() {
        // It is important that the edge primitives can handle vertices that exactly
        // exactly proportional to each other, i.e. that are not identical but are
        // nevertheless exactly coincident when projected onto the unit sphere.
        // There are various ways that such points can arise.  For example,
        // Normalize() itself is not idempotent: there exist distinct points A,B
        // such that Normalize(A) == B  and Normalize(B) == A.  Another issue is
        // that sometimes calls to Normalize() are skipped when the result of a
        // calculation "should" be unit length mathematically (e.g., when computing
        // the cross product of two orthonormal vectors).
        //
        // This test checks pairs of edges AB and CD where A,B,C,D are exactly
        // coincident on the sphere and the norms of A,B,C,D are monotonically
        // increasing.  Such edge pairs should never intersect.  (This is not
        // obvious, since it depends on the particular symbolic perturbations used
        // by s2pred::Sign().  It would be better to replace this with a test that
        // says that the CCW results must be consistent with each other.)
        val kIters = 1000
        repeat(kIters) {
            // Construct a point P where every component is zero or a power of 2.
            val p = S2Point()
            for (i in 0..2) {
                val binary_exp = S2Random.skewed(11)
                p[i] = if (binary_exp > 1022) 0.0 else 2.0.pow(-binary_exp.toDouble())
            }
            // If all components were zero, try again.  Note that normalization may
            // convert a non-zero point into a zero one due to underflow (!)
            if (p.norm2() == 0.0) return@repeat
            p.normalize()
            if (p[0] == 0.0 && p[1] == 0.0 && p[2] == 0.0) return@repeat

            // Now every non-zero component should have exactly the same mantissa.
            // This implies that if we scale the point by an arbitrary factor, every
            // non-zero component will still have the same mantissa.  Scale the points
            // so that they are all distinct and are still very likely to satisfy
            // S2::IsUnitLength (which allows for a small amount of error in the norm).
            val a = (1 - 3e-16) * p
            val b = (1 - 1e-16) * p
            val c = p
            val d = (1 + 2e-16) * p
            if (!S2PointUtil.isUnitLength(a) || !S2PointUtil.isUnitLength(d)) return@repeat
            // Verify that the expected edges do not cross.
            assertThat(S2EdgeCrossings.crossingSign(a, b, c, d)).isLessThanOrEqualTo(0)
            val crosser = S2EdgeCrosser(a, b, c)
            assertThat(crosser.crossingSign(d)).isLessThanOrEqualTo(0)
            assertThat(crosser.crossingSign(c)).isLessThanOrEqualTo(0)
        }
    }


}

