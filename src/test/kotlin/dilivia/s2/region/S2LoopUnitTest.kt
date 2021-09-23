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
package dilivia.s2.region

import dilivia.PreConditions
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.R1Interval
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.region.S2LoopUnitTest.Companion.RelationFlags.CONTAINED
import dilivia.s2.region.S2LoopUnitTest.Companion.RelationFlags.CONTAINS
import dilivia.s2.region.S2LoopUnitTest.Companion.RelationFlags.COVERS
import dilivia.s2.region.S2LoopUnitTest.Companion.RelationFlags.DISJOINT
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import java.util.*
import kotlin.math.abs
import kotlin.math.acos
import kotlin.math.asin
import kotlin.math.cos
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.tan


class S2LoopUnitTest {

    companion object {

        private val logger = KotlinLogging.logger(S2LoopUnitTest::class.java.name)

        // The set of all loops declared below.
        val allLoops = mutableListOf<S2Loop>()

        // Some standard loops to use in the tests (see descriptions below).
        // The empty loop.
        val empty: S2Loop = addLoop(S2Loop(S2Loop.kEmpty))

        // The full loop.
        val full: S2Loop = addLoop(S2Loop(S2Loop.kFull))

        // The northern hemisphere, defined using two pairs of antipodal points.
        val north_hemi: S2Loop = addLoop("0:-180, 0:-90, 0:0, 0:90")

        // The northern hemisphere, defined using three points 120 degrees apart.
        val north_hemi3: S2Loop = addLoop("0:-180, 0:-60, 0:60")

        // The southern hemisphere, defined using two pairs of antipodal points.
        val south_hemi: S2Loop = addLoop("0:90, 0:0, 0:-90, 0:-180")

        // The western hemisphere, defined using two pairs of antipodal points.
        val west_hemi: S2Loop = addLoop("0:-180, -90:0, 0:0, 90:0")

        // The eastern hemisphere, defined using two pairs of antipodal points.
        val east_hemi: S2Loop = addLoop("90:0, 0:0, -90:0, 0:-180")

        // The "near" hemisphere, defined using two pairs of antipodal points.
        val near_hemi: S2Loop = addLoop("0:-90, -90:0, 0:90, 90:0")

        // The "far" hemisphere, defined using two pairs of antipodal points.
        val far_hemi: S2Loop = addLoop("90:0, 0:90, -90:0, 0:-90")

        // A spiral stripe that slightly over-wraps the equator.
        val candy_cane: S2Loop = addLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")

        // A small clockwise loop in the northern & eastern hemisperes.
        val small_ne_cw: S2Loop = addLoop("35:20, 45:20, 40:25")

        // Loop around the north pole at 80 degrees.
        val arctic_80: S2Loop = addLoop("80:-150, 80:-30, 80:90")

        // Loop around the south pole at 80 degrees.
        val antarctic_80: S2Loop = addLoop("-80:120, -80:0, -80:-120")

        // A completely degenerate triangle along the equator that Sign()
        // considers to be CCW.
        val line_triangle: S2Loop = addLoop("0:1, 0:2, 0:3")

        // A nearly-degenerate CCW chevron near the equator with very long sides
        // (about 80 degrees).  Its area is less than 1e-640, which is too small
        // to represent in double precision.
        val skinny_chevron: S2Loop = addLoop("0:0, -1e-320:80, 0:1e-320, 1e-320:80")

        // A diamond-shaped loop around the point 0:180.
        val loop_a: S2Loop = addLoop("0:178, -1:180, 0:-179, 1:-180")

        // Another diamond-shaped loop around the point 0:180.
        val loop_b: S2Loop = addLoop("0:179, -1:180, 0:-178, 1:-180")

        // The intersection of A and B.
        val a_intersect_b: S2Loop = addLoop("0:179, -1:180, 0:-179, 1:-180")

        // The union of A and B.
        val a_union_b: S2Loop = addLoop("0:178, -1:180, 0:-178, 1:-180")

        // A minus B (concave).
        val a_minus_b: S2Loop = addLoop("0:178, -1:180, 0:179, 1:-180")

        // B minus A (concave).
        val b_minus_a: S2Loop = addLoop("0:-179, -1:180, 0:-178, 1:-180")

        // A shape gotten from A by adding a triangle to one edge, and
        // subtracting a triangle from the opposite edge.
        val loop_c: S2Loop = addLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")

        // A shape gotten from A by adding a triangle to one edge, and
        // adding another triangle to the opposite edge.
        val loop_d: S2Loop = addLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")

        //   3------------2
        //   |            |               ^
        //   |  7-8  b-c  |               |
        //   |  | |  | |  |      Latitude |
        //   0--6-9--a-d--1               |
        //   |  | |       |               |
        //   |  f-e       |               +----------->
        //   |            |                 Longitude
        //   4------------5
        //
        // Important: It is not okay to skip over collinear vertices when
        // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
        // uses symbolic perturbations to ensure that no three vertices are
        // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
        // collinear).  In other words, it is unpredictable (modulo knowing the
        // details of the symbolic perturbations) whether 0123 contains 06123,
        // for example.
        //
        // Loop E:  0,6,9,a,d,1,2,3
        // Loop F:  0,4,5,1,d,a,9,6
        // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
        // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
        // Loop I:  7,6,f,e,9,8
        val loop_e: S2Loop = addLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30")
        val loop_f: S2Loop = addLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34")
        val loop_g: S2Loop = addLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
        val loop_h: S2Loop = addLoop("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
        val loop_i: S2Loop = addLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36")

        // Like loop_a, but the vertices are at leaf cell centers.
        val snapped_loop_a: S2Loop = addLoop(
            S2Loop(
                listOf<S2Point>(
                    S2CellId.fromPoint(S2TextParser.makePoint("0:178")).toPoint(),
                    S2CellId.fromPoint(S2TextParser.makePoint("-1:180")).toPoint(),
                    S2CellId.fromPoint(S2TextParser.makePoint("0:-179")).toPoint(),
                    S2CellId.fromPoint(S2TextParser.makePoint("1:-180")).toPoint()
                )
            )
        )

        val kRectError = S2LatLngRectBounder.maxErrorForTests()


        fun addLoop(str: String): S2Loop {
            logger.trace { "Add loop: $str" }
            return addLoop(S2TextParser.makeLoop(str))
        }

        fun addLoop(loop: S2Loop): S2Loop {
            logger.trace { "Add loop: $loop" }
            allLoops.add(loop)
            return loop
        }

        fun rotate(loop: S2Loop): S2Loop {
            val vertices = mutableListOf<S2Point>()
            for (i in 1..loop.numVertices) {
                vertices.add(loop.vertex(i))
            }
            return S2Loop(vertices)
        }

        // Check that the curvature is *identical* when the vertex order is
        // rotated, and that the sign is inverted when the vertices are reversed.
        fun checkCurvatureInvariants(loop: S2Loop) {
            val expected = loop.curvature
            var loop_copy = loop.clone()
            repeat(loop.numVertices) {
                loop_copy.invert()
                assertThat(loop_copy.curvature).isEqualTo(-expected)
                loop_copy.invert()
                loop_copy = rotate(loop_copy)
                assertThat(loop_copy.curvature).isEqualTo(expected)
            }
        }

        // Checks that if a loop is normalized, it doesn't contain a
        // point outside of it, and vice versa.
        fun checkNormalizeAndContains(loop: S2Loop) {
            val p: S2Point = S2TextParser.makePoint("40:40")

            val flip = loop.clone()
            flip.invert()
            assertThat(loop.isNormalized() xor loop.contains(p)).isTrue()
            assertThat(flip.isNormalized() xor flip.contains(p)).isTrue()

            assertThat(loop.isNormalized() xor flip.isNormalized()).isTrue()

            flip.normalize()
            assertThat(flip.contains(p)).isFalse()
        }

        // Given a pair of loops where A contains B, check various identities.
        fun yestOneNestedPair(a: S2Loop, b: S2Loop) {
            assertThat(a.contains(b)).isTrue()
            assertThat(b.contains(a)).isEqualTo(a.boundaryEquals(b))
            assertThat(a.intersects(b)).isEqualTo(!b.isEmpty())
            assertThat(b.intersects(a)).isEqualTo(!b.isEmpty())
        }

        // Given a pair of disjoint loops A and B, check various identities.
        fun testOneDisjointPair(a: S2Loop, b: S2Loop) {
            assertThat(a.intersects(b)).isFalse()
            assertThat(b.intersects(a)).isFalse()
            assertThat(a.contains(b)).isEqualTo(b.isEmpty())
            assertThat(b.contains(a)).isEqualTo(a.isEmpty())
        }

        // Given loops A and B whose union covers the sphere, check various identities.
        fun testOneCoveringPair(a: S2Loop, b: S2Loop) {
            assertThat(a.contains(b)).isEqualTo(a.isFull())
            assertThat(b.contains(a)).isEqualTo(b.isFull())
            val a1 = (a.clone())
            a1.invert()
            val complementary = a1.boundaryEquals(b)
            assertThat(a.intersects(b)).isEqualTo(!complementary)
            assertThat(b.intersects(a)).isEqualTo(!complementary)
        }

        // Given loops A and B such that both A and its complement intersect both B
        // and its complement, check various identities.
        fun testOneOverlappingPair(a: S2Loop, b: S2Loop) {
            assertThat(a.contains(b)).isFalse()
            assertThat(b.contains(a)).isFalse()
            assertThat(a.intersects(b)).isTrue()
            assertThat(b.intersects(a)).isTrue()
        }

        // Given a pair of loops where A contains B, test various identities
        // involving A, B, and their complements.
        fun testNestedPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            val b1 = (b.clone())
            a1.invert()
            b1.invert()
            yestOneNestedPair(a, b)
            yestOneNestedPair(b1, a1)
            testOneDisjointPair(a1, b)
            testOneCoveringPair(a, b1)
        }

        // Given a pair of disjoint loops A and B, test various identities
        // involving A, B, and their complements.
        fun testDisjointPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            a1.invert()
            testNestedPair(a1, b)
        }

        // Given loops A and B whose union covers the sphere, test various identities
        // involving A, B, and their complements.
        fun testCoveringPair(a: S2Loop, b: S2Loop) {
            val b1 = (b.clone())
            b1.invert()
            testNestedPair(a, b1)
        }

        // Given loops A and B such that both A and its complement intersect both B
        // and its complement, test various identities involving these four loops.
        fun testOverlappingPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            val b1 = (b.clone())
            a1.invert()
            b1.invert()
            testOneOverlappingPair(a, b)
            testOneOverlappingPair(a1, b1)
            testOneOverlappingPair(a1, b)
            testOneOverlappingPair(a, b1)
        }

        object RelationFlags {
            val CONTAINS = 0x01  // A contains B
            val CONTAINED = 0x02  // B contains A
            val DISJOINT = 0x04  // A and B are disjoint (intersection is empty)
            val COVERS = 0x08  // (A union B) covers the entire sphere
        }

        // Verify the relationship between two loops A and B.  "flags" is the set of
        // RelationFlags that apply.  "shared_edge" means that the loops share at
        // least one edge (possibly reversed).
        fun testRelationWithDesc(a: S2Loop, b: S2Loop, flags: Int, shared_edge: Boolean, test_description: String) {
            logger.info { test_description }
            if (flags and CONTAINS != 0) {
                testNestedPair(a, b)
            }
            if (flags and CONTAINED != 0) {
                testNestedPair(b, a)
            }
            if (flags and COVERS != 0) {
                testCoveringPair(a, b)
            }
            if (flags and DISJOINT != 0) {
                testDisjointPair(a, b)
            } else if ((flags and (CONTAINS or CONTAINED or COVERS)) == 0) {
                testOverlappingPair(a, b)
            }
            if (!shared_edge && (flags and (CONTAINS or CONTAINED or DISJOINT)) != 0) {
                assertThat(a.containsNested(b)).isEqualTo(a.contains(b))
            }
            // A contains the boundary of B if either A contains B, or the two loops
            // contain each other's boundaries and there are no shared edges (since at
            // least one such edge must be reversed, and therefore is not considered to
            // be contained according to the rules of CompareBoundary).
            var comparison = 0
            if ((flags and CONTAINS) != 0 || ((flags and COVERS) != 0 && !shared_edge)) {
                comparison = 1
            }
            // Similarly, A excludes the boundary of B if either A and B are disjoint,
            // or B contains A and there are no shared edges (since A is considered to
            // contain such edges according to the rules of CompareBoundary).
            if ((flags and DISJOINT) != 0 || ((flags and CONTAINED) != 0 && !shared_edge)) {
                comparison = -1
            }
            // CompareBoundary requires that neither loop is empty.
            if (!a.isEmpty() && !b.isEmpty()) {
                assertThat(a.compareBoundary(b)).isEqualTo(comparison)
            }
        }

        fun testRelation(a: S2Loop, b: S2Loop, flags: Int, shared_edge: Boolean) =
            testRelationWithDesc(a, b, flags, shared_edge, "Test $flags (sharedEdge = $shared_edge) with args $a , $b")

        fun makeCellLoop(begin: S2CellId, end: S2CellId): S2Loop {
            // Construct a CCW polygon whose boundary is the union of the cell ids
            // in the range [begin, end).  We add the edges one by one, removing
            // any edges that are already present in the opposite direction.
            val edges = TreeMap<S2Point, MutableSet<S2Point>>()
            var id = begin
            while (id != end) {
                val cell = S2Cell(id)
                for (k in 0..3) {
                    val a = cell.getVertex(k)
                    val b = cell.getVertex(k + 1)
                    val edgesB = edges[b]
                    val edgesA = edges[a]
                    if (edgesB == null || !edgesB.remove(a)) {
                        if (edgesA == null) edges[a] = mutableSetOf(b)
                        else edgesA.add(b)
                    } else if (edgesB.isEmpty()) {
                        edges.remove(b)
                    }
                }
                id = id.next()
            }

            // The remaining edges form a single loop.  We simply follow it starting
            // at an arbitrary vertex and build up a list of vertices.
            val vertices = mutableListOf<S2Point>()
            var p = edges.firstKey()
            while (!edges.isEmpty()) {
                assertThat(edges.getValue(p).size).isEqualTo(1)
                val next = edges.getValue(p).first()
                vertices.add(p)
                edges.remove(p)
                p = next
            }

            return S2Loop(vertices)
        }

        fun testNear(a_str: String, b_str: String, max_error: S1Angle, expected: Boolean) {
            val a = S2TextParser.makeLoop(a_str)
            val b = S2TextParser.makeLoop(b_str)
            assertThat(a.boundaryNear(b, max_error)).isEqualTo(expected)
            assertThat(b.boundaryNear(a, max_error)).isEqualTo(expected)
        }

        fun checkIdentical(loop: S2Loop, loop2: S2Loop) {
            assertThat(loop2.depth).isEqualTo(loop.depth)
            assertThat(loop2.numVertices).isEqualTo(loop.numVertices)
            for (i in 0 until loop.numVertices) {
                assertThat(loop2.vertex(i)).isEqualTo(loop.vertex(i))
            }
            assertThat(loop2.isEmpty()).isEqualTo(loop.isEmpty())
            assertThat(loop2.isFull()).isEqualTo(loop.isFull())
            assertThat(loop2.depth).isEqualTo(loop.depth)
            assertThat(loop2.isNormalized()).isEqualTo(loop.isNormalized())
            assertThat(loop2.contains(S2PointUtil.origin())).isEqualTo(loop.contains(S2PointUtil.origin()))
            assertThat(loop2.rectBound).isEqualTo(loop.rectBound)
        }
    }

    @Test
    fun getRectBound() {
        assertThat(empty.rectBound.isEmpty).isTrue()
        assertThat(full.rectBound.isFull).isTrue()
        assertThat(candy_cane.rectBound.lng.isFull).isTrue()
        assertThat(candy_cane.rectBound.latLo().degrees() <= -20).isTrue()
        assertThat(candy_cane.rectBound.latHi().degrees() >= 10).isTrue()
        assertThat(small_ne_cw.rectBound.isFull).isTrue()
        assertThat(
            arctic_80.rectBound.approxEquals(
                S2LatLngRect(
                    S2LatLng.fromDegrees(80, -180),
                    S2LatLng.fromDegrees(90, 180)
                ), kRectError
            )
        ).isTrue()
        assertThat(
            antarctic_80.rectBound.approxEquals(
                S2LatLngRect(
                    S2LatLng.fromDegrees(-90, -180),
                    S2LatLng.fromDegrees(-80, 180)
                ), kRectError
            )
        ).isTrue()

        // Create a loop that contains the complement of the "arctic_80" loop.
        val arctic_80inv = arctic_80.clone()
        arctic_80inv.invert()
        // The highest latitude of each edge is attained at its midpoint.
        val mid = (arctic_80inv.vertex(0) + arctic_80inv.vertex(1)) * 0.5
        assertThat(S2LatLng.fromPoint(mid).lat().radians).isCloseTo(
            arctic_80inv.rectBound.latHi().radians,
            Offset.offset(kRectError.lat().radians)
        )

        assertThat(south_hemi.rectBound.lng.isFull).isTrue()
        assertThat(south_hemi.rectBound.lat.approxEquals(R1Interval(-M_PI_2, 0.0), kRectError.lat().radians)).isTrue()
    }

    @Test
    fun areaConsistentWithCurvature() {
        // Check that the area computed using GetArea() is consistent with the
        // curvature of the loop computed using GetTurnangle().  According to
        // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
        // minus its curvature.
        for (loop in allLoops) {
            val area = loop.area
            val gauss_area = 2 * M_PI - loop.curvature
            // The error bound is sufficient for current tests but not guaranteed.
            assertThat(abs(area - gauss_area) <= 1e-14)
                .withFailMessage("Failed loop: $loop\nArea = $area, Gauss Area = $gauss_area")
                .isTrue()
        }
    }

    @Test
    fun getAreaConsistentWithSign() {
        // Test that GetArea() returns an area near 0 for degenerate loops that
        // contain almost no points, and an area near 4*Pi for degenerate loops that
        // contain almost all points.
        val kMaxVertices = 6
        repeat(50) { i ->
            S2Random.reset(i + 11)
            val num_vertices = 3 + S2Random.randomInt(kMaxVertices - 3 + 1)
            // Repeatedly choose N vertices that are exactly on the equator until we
            // find some that form a valid loop.
            var loop: S2Loop
            do {
                val vertices = mutableListOf<S2Point>()
                repeat(num_vertices) {
                    // We limit longitude to the range [0, 90] to ensure that the loop is
                    // degenerate (as opposed to following the entire equator).
                    vertices.add(S2LatLng.fromRadians(0.0, S2Random.randomDouble() * M_PI_2).toPoint())
                }
                loop = S2Loop(vertices, debugOverride = S2Debug.DISABLE)
            } while (!loop.isValid())
            val ccw = loop.isNormalized()
            val area = loop.area
            assertThat(area)
                .withFailMessage("Failed loop $i: ccw = $ccw, area = $area\n${S2TextParser.toString(loop)}")
                .isCloseTo(if (ccw) 0.0 else 4 * M_PI, Offset.offset(1e-15))
            assertThat(loop.contains(S2Point(0, 0, 1))).isEqualTo(!ccw)
        }
    }

    fun testGetAreaAccuracy() {
        // TODO(ericv): Test that GetArea() has an accuracy significantly better
        // than 1e-15 on loops whose area is small.
    }

    @Test
    fun getAreaAndCentroid() {
        assertThat(empty.area).isEqualTo(0.0)
        assertThat(full.area).isEqualTo(4 * M_PI)
        assertThat(empty.centroid).isEqualTo(S2Point(0, 0, 0))
        assertThat(full.centroid).isEqualTo(S2Point(0, 0, 0))

        assertThat(2 * M_PI).isEqualTo(north_hemi.area)
        assertThat(east_hemi.area).isCloseTo(2 * M_PI, Offset.offset(1e-15))

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices.  Then check that the area and centroid are
        // correct.

        repeat(50) { i ->
            S2Random.reset(i)
            // Choose a coordinate frame for the spherical cap.
            val (x, y, z) = S2Random.randomFrame()

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            val kMaxDist = 1e-6
            val height = 2 * S2Random.randomDouble()
            val phi = asin(1 - height)
            var max_dtheta = 2 * acos(tan(abs(phi)) / tan(abs(phi) + kMaxDist))
            max_dtheta = min(M_PI, max_dtheta);  // At least 3 vertices.

            val vertices = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                vertices.add(x * (cos(theta) * cos(phi)) + y * (sin(theta) * cos(phi)) + z * sin(phi))
                theta += S2Random.randomDouble() * max_dtheta
            }
            val loop = S2Loop(vertices)
            val area = loop.area
            val centroid = loop.centroid
            val expectedArea = 2 * M_PI * height
            assertThat(abs(area - expectedArea) <= 2 * M_PI * kMaxDist).isTrue()
            val expectedCentroid = z * (expectedArea * (1 - 0.5 * height))
            assertThat((centroid - expectedCentroid).norm() <= 2 * kMaxDist).isTrue()
        }
    }

    @Test
    fun getCurvature() {
        assertThat(empty.curvature).isEqualTo(2 * M_PI)
        assertThat(full.curvature).isEqualTo(-2 * M_PI)

        assertThat(north_hemi3.curvature).isCloseTo(0.0, Offset.offset(1e-15))
        checkCurvatureInvariants(north_hemi3)

        assertThat(west_hemi.curvature).isCloseTo(0.0, Offset.offset(1e-15))
        checkCurvatureInvariants(west_hemi)

        // We don't have an easy way to estimate the curvature of this loop, but
        // we can still check that the expected invariants hold.
        checkCurvatureInvariants(candy_cane)

        assertThat(line_triangle.curvature).isCloseTo(2 * M_PI, Offset.offset(1e-15))
        checkCurvatureInvariants(line_triangle)

        assertThat(skinny_chevron.curvature).isCloseTo(2 * M_PI, Offset.offset(1e-15))
        checkCurvatureInvariants(skinny_chevron)

        // Build a narrow spiral loop starting at the north pole.  This is designed
        // to test that the error in GetCurvature is linear in the number of
        // vertices even when the partial sum of the curvatures gets very large.
        // The spiral consists of two "arms" defining opposite sides of the loop.
        val kArmPoints = 10000;    // Number of vertices in each "arm"
        val kArmRadius = 0.01;  // Radius of spiral.
        val vertices = Array<S2Point>(2 * kArmPoints) { S2Point() }
        vertices[kArmPoints] = S2Point(0, 0, 1)
        for (i in 0 until kArmPoints) {
            val angle = (2 * M_PI / 3) * i
            val x = cos(angle)
            val y = sin(angle)
            val r1 = i * kArmRadius / kArmPoints
            val r2 = (i + 1.5) * kArmRadius / kArmPoints
            vertices[kArmPoints - i - 1] = S2Point(r1 * x, r1 * y, 1.0).normalize()
            vertices[kArmPoints + i] = S2Point(r2 * x, r2 * y, 1.0).normalize()
        }
        // This is a pathological loop that contains many long parallel edges, and
        // takes tens of seconds to validate in debug mode.
        val spiral = S2Loop(vertices.toMutableList(), debugOverride = S2Debug.DISABLE)

        // Check that GetCurvature() is consistent with GetArea() to within the
        // error bound of the former.  We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that.  The error in GetArea() can be ignored for the
        // purposes of this test since it is generally much smaller.
        assertThat(spiral.curvature).isCloseTo(
            2 * M_PI - spiral.area,
            Offset.offset(0.01 * spiral.getCurvatureMaxError())
        )
    }

    @Test
    fun normalizedCompatibleWithContains() {
        checkNormalizeAndContains(line_triangle)
        checkNormalizeAndContains(skinny_chevron)
    }

    @Test
    fun contains() {
        // Check the full and empty loops have the correct containment relationship
        // with the special "vertex" that defines them.
        assertThat(empty.contains(S2Loop.kEmpty[0])).isFalse()
        assertThat(full.contains(S2Loop.kFull[0])).isTrue()

        assertThat(candy_cane.contains(S2LatLng.fromDegrees(5, 71).toPoint())).isTrue()

        // Create copies of these loops so that we can change the vertex order.
        var north_copy = (north_hemi.clone())
        var south_copy = (south_hemi.clone())
        var west_copy = (west_hemi.clone())
        var east_copy = (east_hemi.clone())
        repeat(4) {
            assertThat(north_copy.contains(S2Point(0, 0, 1))).isTrue()
            assertThat(north_copy.contains(S2Point(0, 0, -1))).isFalse()
            assertThat(south_copy.contains(S2Point(0, 0, 1))).isFalse()
            assertThat(south_copy.contains(S2Point(0, 0, -1))).isTrue()
            assertThat(west_copy.contains(S2Point(0, 1, 0))).isFalse()
            assertThat(west_copy.contains(S2Point(0, -1, 0))).isTrue()
            assertThat(east_copy.contains(S2Point(0, 1, 0))).isTrue()
            assertThat(east_copy.contains(S2Point(0, -1, 0))).isFalse()
            north_copy = rotate(north_copy)
            south_copy = rotate(south_copy)
            east_copy = rotate(east_copy)
            west_copy = rotate(west_copy)
        }

        // This code checks each cell vertex is contained by exactly one of
        // the adjacent cells.
        for (level in 0..2) {
            val loops = mutableListOf<S2Loop>()
            val loop_vertices = mutableListOf<S2Point>()
            val points = mutableSetOf<S2Point>()
            var id = S2CellId.begin(level)
            while (id != S2CellId.end(level)) {
                val cell = S2Cell(id)
                points.add(cell.getCenter())
                for (k in 0..3) {
                    loop_vertices.add(cell.getVertex(k))
                    points.add(cell.getVertex(k))
                }
                loops.add(S2Loop(loop_vertices))
                loop_vertices.clear()
                id = id.next()
            }
            for (point in points) {
                var count = 0
                for (loop in loops) {
                    if (loop.contains(point)) ++count
                }
                assertThat(1).isEqualTo(count)
            }
        }
    }

    @Test
    fun containsMatchesCrossingSign() {
        // This test demonstrates a former incompatibility between CrossingSign()
        // and Contains(const S2Point&).  It constructs an S2Cell-based loop L and
        // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
        // previously, Contains() returned false for both endpoints of E.
        //
        // The reason for the bug was that the loop bound was sometimes too tight.
        // The Contains() code for a0 bailed out early because a0 was found not to
        // be inside the bound of L.

        // Start with a cell that ends up producing the problem.
        val cell_id = S2CellId.fromPoint(S2Point(1, 1, 1)).parent(21)

        val children = Array(4) { S2Cell() }
        S2Cell(cell_id).subdivide(children)

        val points = Array<S2Point>(4) { S2Point() }
        for (i in 0..3) {
            // Note extra normalization. GetCenter() is already normalized.
            // The test results will no longer be inconsistent if the extra
            // Normalize() is removed.
            points[i] = children[i].getCenter().normalize()
        }

        val loop = S2Loop(points.toMutableList())

        // Get a vertex from a grandchild cell.
        // +---------------+---------------+
        // |               |               |
        // |    points[3]  |   points[2]   |
        // |       v       |       v       |
        // |       +-------+------ +       |
        // |       |       |       |       |
        // |       |       |       |       |
        // |       |       |       |       |
        // +-------+-------+-------+-------+
        // |       |       |       |       |
        // |       |    <----------------------- grandchild_cell
        // |       |       |       |       |
        // |       +-------+------ +       |
        // |       ^       |       ^       | <-- cell
        // | points[0]/a0  |     points[1] |
        // |               |               |
        // +---------------+---------------+
        val grandchild_cell = S2Cell(cell_id.child(0).child(2))
        val a0 = grandchild_cell.getVertex(0)

//        // If this doesn't hold, the rest of the test is pointless.
//        assertThat(points[0] != a0)
//            .withFailMessage("This test depends on rounding errors that should make a0 slightly different from points[0]\npoints[0]: ${points[0]}\n       a0: $a0")
//            .isTrue()
//
//        // The edge from a0 to the origin crosses one boundary.
//        assertThat(S2EdgeCrossings.crossingSign(a0, S2PointUtil.origin(), loop.vertex(0), loop.vertex(1))).isEqualTo(-1)
//        assertThat(S2EdgeCrossings.crossingSign(a0, S2PointUtil.origin(), loop.vertex(1), loop.vertex(2))).isEqualTo(1)
//        assertThat(S2EdgeCrossings.crossingSign(a0, S2PointUtil.origin(), loop.vertex(2), loop.vertex(3))).isEqualTo(-1)
//        assertThat(S2EdgeCrossings.crossingSign(a0, S2PointUtil.origin(), loop.vertex(3), loop.vertex(4))).isEqualTo(-1)
//
//        // Contains should return false for the origin, and true for a0.
//        assertThat(loop.contains(S2PointUtil.origin())).isFalse()
//        assertThat(loop.contains(a0)).isTrue()
//
//        // Since a0 is inside the loop, it should be inside the bound.
//        val bound = loop.rectBound
//        assertThat(bound.contains(a0)).isTrue()
    }

    @Test
    fun loopRelations() {
        // Check full and empty relationships with normal loops and each other.
        testRelation(full, full, CONTAINS or CONTAINED or COVERS, true)
        testRelation(full, north_hemi, CONTAINS or COVERS, false)
        testRelation(full, empty, CONTAINS or DISJOINT or COVERS, false)
        testRelation(north_hemi, full, CONTAINED or COVERS, false)
        testRelation(north_hemi, empty, CONTAINS or DISJOINT, false)
        testRelation(empty, full, CONTAINED or DISJOINT or COVERS, false)
        testRelation(empty, north_hemi, CONTAINED or DISJOINT, false)
        testRelation(empty, empty, CONTAINS or CONTAINED or DISJOINT, false)

        testRelation(north_hemi, north_hemi, CONTAINS or CONTAINED, true)
        testRelation(north_hemi, south_hemi, DISJOINT or COVERS, true)
        testRelation(north_hemi, east_hemi, 0, false)
        testRelation(north_hemi, arctic_80, CONTAINS, false)
        testRelation(north_hemi, antarctic_80, DISJOINT, false)
        testRelation(north_hemi, candy_cane, 0, false)

        // We can't compare north_hemi3 vs. north_hemi or south_hemi because the
        // result depends on the "simulation of simplicity" implementation details.
        testRelation(north_hemi3, north_hemi3, CONTAINS or CONTAINED, true)
        testRelation(north_hemi3, east_hemi, 0, false)
        testRelation(north_hemi3, arctic_80, CONTAINS, false)
        testRelation(north_hemi3, antarctic_80, DISJOINT, false)
        testRelation(north_hemi3, candy_cane, 0, false)

        testRelation(south_hemi, north_hemi, DISJOINT or COVERS, true)
        testRelation(south_hemi, south_hemi, CONTAINS or CONTAINED, true)
        testRelation(south_hemi, far_hemi, 0, false)
        testRelation(south_hemi, arctic_80, DISJOINT, false)
        testRelation(south_hemi, antarctic_80, CONTAINS, false)
        testRelation(south_hemi, candy_cane, 0, false)

        testRelation(candy_cane, north_hemi, 0, false)
        testRelation(candy_cane, south_hemi, 0, false)
        testRelation(candy_cane, arctic_80, DISJOINT, false)
        testRelation(candy_cane, antarctic_80, DISJOINT, false)
        testRelation(candy_cane, candy_cane, CONTAINS or CONTAINED, true)

        testRelation(near_hemi, west_hemi, 0, false)

        testRelation(small_ne_cw, south_hemi, CONTAINS, false)
        testRelation(small_ne_cw, west_hemi, CONTAINS, false)

        testRelation(small_ne_cw, north_hemi, COVERS, false)
        testRelation(small_ne_cw, east_hemi, COVERS, false)

        testRelation(loop_a, loop_a, CONTAINS or CONTAINED, true)
        testRelation(loop_a, loop_b, 0, false)
        testRelation(loop_a, a_intersect_b, CONTAINS, true)
        testRelation(loop_a, a_union_b, CONTAINED, true)
        testRelation(loop_a, a_minus_b, CONTAINS, true)
        testRelation(loop_a, b_minus_a, DISJOINT, true)

        testRelation(loop_b, loop_a, 0, false)
        testRelation(loop_b, loop_b, CONTAINS or CONTAINED, true)
        testRelation(loop_b, a_intersect_b, CONTAINS, true)
        testRelation(loop_b, a_union_b, CONTAINED, true)
        testRelation(loop_b, a_minus_b, DISJOINT, true)
        testRelation(loop_b, b_minus_a, CONTAINS, true)

        testRelation(a_intersect_b, loop_a, CONTAINED, true)
        testRelation(a_intersect_b, loop_b, CONTAINED, true)
        testRelation(a_intersect_b, a_intersect_b, CONTAINS or CONTAINED, true)
        testRelation(a_intersect_b, a_union_b, CONTAINED, false)
        testRelation(a_intersect_b, a_minus_b, DISJOINT, true)
        testRelation(a_intersect_b, b_minus_a, DISJOINT, true)

        testRelation(a_union_b, loop_a, CONTAINS, true)
        testRelation(a_union_b, loop_b, CONTAINS, true)
        testRelation(a_union_b, a_intersect_b, CONTAINS, false)
        testRelation(a_union_b, a_union_b, CONTAINS or CONTAINED, true)
        testRelation(a_union_b, a_minus_b, CONTAINS, true)
        testRelation(a_union_b, b_minus_a, CONTAINS, true)

        testRelation(a_minus_b, loop_a, CONTAINED, true)
        testRelation(a_minus_b, loop_b, DISJOINT, true)
        testRelation(a_minus_b, a_intersect_b, DISJOINT, true)
        testRelation(a_minus_b, a_union_b, CONTAINED, true)
        testRelation(a_minus_b, a_minus_b, CONTAINS or CONTAINED, true)
        testRelation(a_minus_b, b_minus_a, DISJOINT, false)

        testRelation(b_minus_a, loop_a, DISJOINT, true)
        testRelation(b_minus_a, loop_b, CONTAINED, true)
        testRelation(b_minus_a, a_intersect_b, DISJOINT, true)
        testRelation(b_minus_a, a_union_b, CONTAINED, true)
        testRelation(b_minus_a, a_minus_b, DISJOINT, false)
        testRelation(b_minus_a, b_minus_a, CONTAINS or CONTAINED, true)
    }

    // Make sure the relations are correct if the loop crossing happens on
    // two ends of a shared boundary segment.
    @Test
    fun loopRelationsWhenSameExceptPiecesStickingOutAndIn() {
        testRelation(loop_a, loop_c, 0, true)
        testRelation(loop_c, loop_a, 0, true)
        testRelation(loop_a, loop_d, CONTAINED, true)
        testRelation(loop_d, loop_a, CONTAINS, true)
        testRelation(loop_e, loop_f, DISJOINT, true)
        testRelation(loop_e, loop_g, CONTAINS, true)
        testRelation(loop_e, loop_h, 0, true)
        testRelation(loop_e, loop_i, 0, false)
        testRelation(loop_f, loop_g, DISJOINT, true)
        testRelation(loop_f, loop_h, 0, true)
        testRelation(loop_f, loop_i, 0, false)
        testRelation(loop_g, loop_h, CONTAINED, true)
        testRelation(loop_h, loop_g, CONTAINS, true)
        testRelation(loop_g, loop_i, DISJOINT, true)
        testRelation(loop_h, loop_i, CONTAINS, true)
    }

    @Test
    fun loopRelations2() {
        // Construct polygons consisting of a sequence of adjacent cell ids
        // at some fixed level.  Comparing two polygons at the same level
        // ensures that there are no T-vertices.
        repeat(1000) {
            var begin = S2CellId((S2Random.randomLong() or 1L).toULong())
            if (!begin.isValid) return@repeat
            begin = begin.parent(S2Random.randomInt(S2CellId.kMaxLevel))
            val a_begin = begin.advance(S2Random.skewed(6).toLong())
            val a_end = a_begin.advance(S2Random.skewed(6).toLong() + 1L)
            val b_begin = begin.advance(S2Random.skewed(6).toLong())
            val b_end = b_begin.advance(S2Random.skewed(6).toLong() + 1L)
            if (!a_end.isValid || !b_end.isValid) return@repeat

            val a = (makeCellLoop(a_begin, a_end))
            val b = (makeCellLoop(b_begin, b_end))
            val contained = (a_begin <= b_begin && b_end <= a_end)
            val intersects = (a_begin < b_end && b_begin < a_end)
            logger.trace { "Checking ${a.numVertices} vs. ${b.numVertices}, contained = $contained, intersects = $intersects" }

            logger.trace { "a = ${a.toDebugString()}" }
            logger.trace { "b = ${b.toDebugString()}" }
            logger.trace { "a contains b: ${a.contains(b)}" }
            assertThat(contained).isEqualTo(a.contains(b))
            assertThat(intersects).isEqualTo(a.intersects(b))
        }
    }

    @Test
    fun boundsForLoopContainment() {
        // To reliably test whether one loop contains another, the bounds of the
        // outer loop are expanded slightly.  This test constructs examples where
        // this expansion is necessary and verifies that it is sufficient.
        repeat(1000) {
            // We construct a triangle ABC such that A,B,C are nearly colinear, B is
            // the point of maximum latitude, and the edge AC passes very slightly
            // below B (i.e., ABC is CCW).
            val b = (S2Random.randomPoint() + S2Point(0, 0, 1)).normalize()
            val v = b.crossProd(S2Point(0, 0, 1)).normalize()
            val a = S2EdgeDistances.interpolate(S2Random.randomDouble(), -v, b)
            val c = S2EdgeDistances.interpolate(S2Random.randomDouble(), b, v)
            if (S2Predicates.sign(a, b, c) < 0) return@repeat
            // Now construct another point D directly below B, and create two loops
            // ABCD and ACD.
            val d = S2Point(b.x, b.y, 0.0).normalize()
            val vertices = arrayOf(c, d, a, b)  // Reordered for convenience
            val outer = S2Loop(vertices.toMutableList())
            val inner = S2Loop(vertices.slice(0..2).toMutableList())
            // Now because the bounds calculation is less accurate when the maximum is
            // attained along an edge (rather than at a vertex), sometimes the inner
            // loop will have a *larger* bounding box than the outer loop.  We look
            // only for those cases.
            if (outer.rectBound.contains(inner.rectBound)) return@repeat
            assertThat(outer.contains(inner)).isTrue()
        }
    }

    fun debugDumpCrossings(loop: S2Loop) {
        // This function is useful for debugging.
        println(
            """
            |Ortho(v1): ${S2PointUtil.ortho(loop.vertex(1))}
            |Contains(kOrigin): ${loop.contains(S2PointUtil.origin())}
        """.trimMargin()
        )
        for (i in 1..loop.numVertices) {
            val a = S2PointUtil.ortho(loop.vertex(i))
            val b = loop.vertex(i - 1)
            val c = loop.vertex(i + 1)
            val o = loop.vertex(i)

            println(
                "Vertex %d: [%.17g, %.17g, %.17g], %d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %d".format(
                    i, loop.vertex(i).x, loop.vertex(i).y, loop.vertex(i).z,
                    i - 1, i, S2Predicates.sign(b, o, a),
                    i + 1, i, i - 1, S2Predicates.sign(c, o, b),
                    i, i + 1, S2Predicates.sign(a, o, c),
                    S2Predicates.orderedCCW(a, b, c, o)
                )
            )
        }
        for (i in 0 until loop.numVertices + 2) {
            var orig = S2PointUtil.origin()
            val dest: S2Point
            if (i < loop.numVertices) {
                dest = loop.vertex(i)
                println("Origin.%d crosses:".format(i))
            } else {
                dest = S2Point(0, 0, 1)
                if (i == loop.numVertices + 1) orig = loop.vertex(1)
                println("Case %d:".format(i))
            }
            for (j in 0 until loop.numVertices) {
                println(
                    " %d".format(
                        S2EdgeCrossings.edgeOrVertexCrossing(
                            orig,
                            dest,
                            loop.vertex(j),
                            loop.vertex(j + 1)
                        )
                    )
                )
            }
            println()
        }
        for (i in 0..2 step 2) {
            println("Origin.v1 crossing v%d.v1: ".format(i))
            val a = S2PointUtil.ortho(loop.vertex(1))
            val b = loop.vertex(i)
            val c = S2PointUtil.origin()
            val o = loop.vertex(1)
            println(
                "%d1R=%d, M1%d=%d, R1M=%d, crosses: %d".format(
                    i, S2Predicates.sign(b, o, a), i, S2Predicates.sign(c, o, b), S2Predicates.sign(a, o, c),
                    S2EdgeCrossings.edgeOrVertexCrossing(c, o, b, a)
                )
            )
        }
    }

    @Test
    fun boundaryNear() {
        val degree = S1Angle.degrees(1)

        testNear("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", degree * 0.5, true)
        testNear("0:0, 0:3, 0:7, 0:10, 3:7, 5:5", "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1", S1Angle.radians(1e-3), true)

        // All vertices close to some edge, but not equivalent.
        testNear("0:0, 0:2, 2:2, 2:0", "0:0, 1.9999:1, 0:2, 2:2, 2:0", degree * 0.5, false)

        // Two triangles that backtrack a bit on different edges.  A simple
        // greedy matching algorithm would fail on this example.
        val t1 = "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, 2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2"
        val t2 = "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, 0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4"
        testNear(t1, t2, degree * 1.5, true)
        testNear(t1, t2, degree * 0.5, false)
    }

    fun testEmptyFullSnapped(loop: S2Loop, level: Int) {
        PreConditions.requireArgument { loop.isEmptyOrFull() }
        val cellid = S2CellId.fromPoint(loop.vertex(0)).parent(level)
        val vertices = listOf(cellid.toPoint())
        val loop2 = S2Loop(vertices)
        assertThat(loop.boundaryEquals(loop2)).isTrue()
        assertThat(loop.boundaryApproxEquals(loop2)).isTrue()
        assertThat(loop.boundaryNear(loop2)).isTrue()
    }

    // Test converting the empty/full loops to S2LatLng representations.  (We
// don't bother testing E5/E6/E7 because that test is less demanding.)
    fun testEmptyFullLatLng(loop: S2Loop) {
        PreConditions.requireArgument { loop.isEmptyOrFull() }
        val vertices = listOf(S2LatLng.fromPoint(loop.vertex(0)).toPoint())
        val loop2 = S2Loop(vertices)
        assertThat(loop.boundaryEquals(loop2)).isTrue()
        assertThat(loop.boundaryApproxEquals(loop2)).isTrue()
        assertThat(loop.boundaryNear(loop2)).isTrue()
    }

    fun testEmptyFullConversions(loop: S2Loop) {
        testEmptyFullSnapped(loop, S2CellId.kMaxLevel)
        testEmptyFullSnapped(loop, 1);  // Worst case for approximation
        testEmptyFullSnapped(loop, 0)
        testEmptyFullLatLng(loop)
    }

    @Test
    fun emptyFullLossyConversions() {
        // Verify that the empty and full loops can be encoded lossily.
        val empty = S2Loop(S2Loop.kEmpty)
        testEmptyFullConversions(empty)

        val full = S2Loop(S2Loop.kFull)
        testEmptyFullConversions(full)
    }

// This test checks that S2Loops created directly from S2Cells behave
// identically to S2Loops created from the vertices of those cells; this
// previously was not the case, because S2Cells calculate their bounding
// rectangles slightly differently, and S2Loops created from them just copied
// the S2Cell bounds.
    @Test
    fun cellConstructorAndContains() {
        val cell = S2Cell(S2CellId.fromLatLng(S2LatLng.fromE6(40565459, -74645276)))
        val cell_as_loop = S2Loop(cell)

        val vertices = mutableListOf<S2Point>()
        for (i in 0 until cell_as_loop.numVertices) {
            vertices.add(cell_as_loop.vertex(i))
        }
        val loop_copy = S2Loop(vertices)
        assertThat(loop_copy.contains(cell_as_loop)).isTrue()
        assertThat(cell_as_loop.contains(loop_copy)).isTrue()

        // Demonstrates the reason for this test; the cell bounds are more
        // conservative than the resulting loop bounds.
        assertThat(loop_copy.rectBound.contains(cell.rectBound)).isFalse()
    }

    // Construct a loop using s2textformat::MakeLoop(str) and check that it
// produces a validation error that includes "snippet".
    fun checkLoopIsInvalid(str: String, snippet: String) {
        val loop = S2TextParser.makeLoop(str, S2Debug.DISABLE)
        val error = S2Error()
        assertThat(loop.findValidationError(error)).isTrue()
        assertThat(error.text).contains(snippet)
    }

    fun checkLoopIsInvalid(points: List<S2Point>, snippet: String) {
        val l = S2Loop(points, debugOverride = S2Debug.DISABLE)
        val error = S2Error()
        assertThat(l.findValidationError(error)).isTrue()
        assertThat(error.text).contains(snippet)
    }

    @Test
    fun isValidDetectsInvalidLoops() {
        // Not enough vertices.  Note that all single-vertex loops are valid; they
        // are interpreted as being either empty or full.
        checkLoopIsInvalid("", "at least 3 vertices")
        checkLoopIsInvalid("20:20, 21:21", "at least 3 vertices")

        // There is a degenerate edge
        checkLoopIsInvalid("20:20, 20:20, 20:21", "degenerate")
        checkLoopIsInvalid("20:20, 20:21, 20:20", "degenerate")

        // There is a duplicate vertex
        checkLoopIsInvalid("20:20, 21:21, 21:20, 20:20, 20:21", "duplicate vertex")

        // Some edges cross
        checkLoopIsInvalid("20:20, 21:21, 21:20.5, 21:20, 20:21", "crosses")

        // Points with non-unit length (triggers S2_DCHECK failure in debug)
        val savedPreConditionsState = PreConditions.enabled
        PreConditions.enabled = false
        checkLoopIsInvalid(
            listOf(S2Point(2, 0, 0), S2Point(0, 1, 0), S2Point(0, 0, 1)),
            "unit length"
        )
        PreConditions.enabled = savedPreConditionsState

        // Adjacent antipodal vertices
        checkLoopIsInvalid(
            listOf(S2Point(1, 0, 0), S2Point(-1, 0, 0), S2Point(0, 0, 1)),
            "antipodal"
        )
    }

    // Helper function for testing the distance methods.  "boundary_x" is the
// expected result of projecting "x" onto the loop boundary.  For convenience
// it can be set to S2Point() to indicate that (boundary_x == x).
    fun testDistanceMethods(loop: S2Loop, x: S2Point, boundary_x: S2Point) {
        var boundary_x = boundary_x
        // This error is not guaranteed by the implementation but is okay for tests.
        val kMaxError = S1Angle.radians(1e-15)

        if (boundary_x == S2Point()) boundary_x = x
        assertThat(S1Angle(boundary_x, loop.projectToBoundary(x))).isLessThanOrEqualTo(kMaxError)

        if (loop.isEmptyOrFull()) {
            assertThat(loop.distanceToBoundary(x)).isEqualTo(S1Angle.infinity())
        } else {
            // assertEquals only works with doubles.
            assertThat(loop.distanceToBoundary(x).degrees()).isCloseTo(
                S1Angle(x, boundary_x).degrees(),
                Offset.offset(kMaxError.degrees())
            )
        }
        if (loop.contains(x)) {
            assertThat(loop.distance(x)).isEqualTo(S1Angle.zero())
            assertThat(loop.project(x)).isEqualTo(x)
        } else {
            assertThat(loop.distance(x)).isEqualTo(loop.distanceToBoundary(x))
            assertThat(loop.project(x)).isEqualTo(loop.projectToBoundary(x))
        }
    }

    @Test
    fun distanceMethods() {
        // S2ClosestEdgeQuery is already tested, so just do a bit of sanity checking.

        // The empty and full loops don't have boundaries.
        testDistanceMethods(empty, S2Point(0, 1, 0), S2Point())
        testDistanceMethods(full, S2Point(0, 1, 0), S2Point())

        // A CCW square around the S2LatLng point (0,0).  Note that because lines of
        // latitude are curved on the sphere, it is not straightforward to project
        // points onto any edge except along the equator.  (The equator is the only
        // line of latitude that is also a geodesic.)
        val square = S2TextParser.makeLoop("-1:-1, -1:1, 1:1, 1:-1")
        assertThat(square.isNormalized()).isTrue()

        // A vertex.
        testDistanceMethods(square, S2LatLng.fromDegrees(1, -1).toPoint(), S2Point())
        // A point on one of the edges.
        testDistanceMethods(square, S2LatLng.fromDegrees(0.5, 1.0).toPoint(), S2Point())
        // A point inside the square.
        testDistanceMethods(square, S2LatLng.fromDegrees(0.0, 0.5).toPoint(), S2LatLng.fromDegrees(0, 1).toPoint())
        // A point outside the square that projects onto an edge.
        testDistanceMethods(square, S2LatLng.fromDegrees(0, -2).toPoint(), S2LatLng.fromDegrees(0, -1).toPoint())
        // A point outside the square that projects onto a vertex.
        testDistanceMethods(square, S2LatLng.fromDegrees(3, 4).toPoint(), S2LatLng.fromDegrees(1, 1).toPoint())
    }

    @Test
    fun makeRegularLoop() {
        val center = S2LatLng.fromDegrees(80, 135).toPoint()
        val radius = S1Angle.degrees(20)
        val loop = (S2Loop.makeRegularLoop(center, radius, 4))

        assertThat(loop.numVertices).isEqualTo(4)
        val p0 = loop.vertex(0)
        val p1 = loop.vertex(1)
        val p2 = loop.vertex(2)
        val p3 = loop.vertex(3)
        // Make sure that the radius is correct.
        assertThat(S2LatLng.fromPoint(center).getDistance(S2LatLng.fromPoint(p0)).degrees()).isCloseTo(20.0, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(center).getDistance(S2LatLng.fromPoint(p1)).degrees()).isCloseTo(20.0, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(center).getDistance(S2LatLng.fromPoint(p2)).degrees()).isCloseTo(20.0, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(center).getDistance(S2LatLng.fromPoint(p3)).degrees()).isCloseTo(20.0, Offset.offset(1e-14))
        // Make sure that all angles of the polygon are the same.
        assertThat((p1 - p0).angle(p3 - p0)).isEqualTo(M_PI_2)
        assertThat((p2 - p1).angle(p0 - p1)).isEqualTo(M_PI_2)
        assertThat((p3 - p2).angle(p1 - p2)).isCloseTo(M_PI_2, Offset.offset(1e-14))
        assertThat((p0 - p3).angle(p2 - p3)).isCloseTo(M_PI_2, Offset.offset(1e-14))
        // Make sure that all edges of the polygon have the same length.
        assertThat(S2LatLng.fromPoint(p0).getDistance(S2LatLng.fromPoint(p1)).degrees()).isCloseTo(27.990890717782829, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(p1).getDistance(S2LatLng.fromPoint(p2)).degrees()).isCloseTo(27.990890717782829, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(p2).getDistance(S2LatLng.fromPoint(p3)).degrees()).isCloseTo(27.990890717782829, Offset.offset(1e-14))
        assertThat(S2LatLng.fromPoint(p3).getDistance(S2LatLng.fromPoint(p0)).degrees()).isCloseTo(27.990890717782829, Offset.offset(1e-14))

        // Check actual coordinates. This may change if we switch the algorithm
        // intentionally.
        assertThat(S2LatLng.fromPoint(p0).lat().degrees()).isEqualTo(62.162880741097204)
        assertThat(S2LatLng.fromPoint(p0).lng().degrees()).isCloseTo(103.11051028343407 , Offset.offset(2e-14))
        assertThat(S2LatLng.fromPoint(p1).lat().degrees()).isEqualTo(61.955157772928345)
        assertThat(S2LatLng.fromPoint(p1).lng().degrees()).isEqualTo(165.25681963683536)
        assertThat(S2LatLng.fromPoint(p2).lat().degrees()).isEqualTo(75.139812547718478)
        assertThat(S2LatLng.fromPoint(p2).lng().degrees()).isEqualTo(-119.13042521187423)
        assertThat(S2LatLng.fromPoint(p3).lat().degrees()).isEqualTo(75.524190079054392)
        assertThat(S2LatLng.fromPoint(p3).lng().degrees()).isEqualTo(26.392175948257943)
    }

    @Test
    fun loopShapeBasic() {
        val loop = S2TextParser.makeLoop("0:0, 0:1, 1:0")
        val shape = S2Loop.Shape(loop = loop)
        assertThat(shape.loop).isEqualTo(loop)
        assertThat(shape.numEdges).isEqualTo(3)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.chain(0).start).isEqualTo(0)
        assertThat(shape.chain(0).length).isEqualTo(3)
        val edge2 = shape.edge(2)
        assertThat(S2TextParser.toString(edge2.v0)).isEqualTo("1:0")
        assertThat(S2TextParser.toString(edge2.v1)).isEqualTo("0:0")
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun loopShapeEmptyLoop() {
        val loop = S2Loop(S2Loop.kEmpty)
        val shape = S2Loop.Shape(loop = loop)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun loopShapeFullLoop() {
        val loop = S2Loop(S2Loop.kFull)
        val shape = S2Loop.Shape(loop = loop)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isTrue()
        assertThat(shape.getReferencePoint().contained).isTrue()
    }


    fun test() {
        val a =
            S2TextParser.makeLoop("50.303911598656235:-73.86564257495967, 50.28788331840551:-73.75332124455377, 50.20985629631176:-73.7995017621896, 50.193784232475195:-73.68742455910771, 50.17758772848922:-73.57530791436109, 50.16126667018284:-73.46315304269496, 50.23904872109849:-73.41612184648233, 50.316905061298975:-73.36888021009284, 50.300293378356415:-73.2561227535556, 50.37813992369764:-73.20838518387683, 50.36131819042698:-73.0953096522237, 50.4391530567744:-73.04707361235798, 50.45605937628929:-73.16043366279766, 50.472838688392926:-73.27376052499469, 50.48949109563151:-73.38705294925887, 50.50601670439513:-73.5003096868368, 50.42784686954798:-73.54740369844646, 50.34975119353782:-73.59428666639539, 50.365985863534945:-73.70693320281322, 50.38209508952241:-73.81954030436384")
        val b =
            S2TextParser.makeLoop("50.303911598656235:-73.86564257495967, 50.28788331840551:-73.75332124455377, 50.20985629631176:-73.7995017621896, 50.193784232475195:-73.68742455910771, 50.17758772848922:-73.57530791436109, 50.25545199470335:-73.5285596579851, 50.271730156068024:-73.6409598649805, 50.34975119353782:-73.59428666639539, 50.365985863534945:-73.70693320281322, 50.38209508952241:-73.81954030436384")

        println("a contains b: ${a.contains(b)}")
    }
}

