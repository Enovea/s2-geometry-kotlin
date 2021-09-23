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
package dilivia.s2.index.shape

import dilivia.collections.isSorted
import dilivia.math.DoubleType
import dilivia.math.vectors.R2Point
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2CellId
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2Factory.makeRegularPoints
import dilivia.s2.S2PaddedCell
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomInt
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2Coords.faceUvToXyz
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.CrossingType
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.S2EdgeVectorShape
import dilivia.s2.shape.ShapeEdge
import dilivia.s2.shape.ShapeEdgeId
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.lang.Math.pow
import kotlin.math.nextTowards


typealias TestEdge = Pair<S2Point, S2Point>

class S2CrossingEdgeQueryUnitTest {

    fun perturbAtDistance(distance: S1Angle, a0: S2Point, b0: S2Point): S2Point {
        val x = S2EdgeDistances.interpolateAtDistance(distance, a0, b0)
        if (S2Random.oneIn(2)) {
            for (i in 0..2) {
                x[i] = x[i].nextTowards(if (S2Random.oneIn(2)) -1.0 else 1.0)
            }
            x.normalize()
        }
        return x
    }

    // Generate sub-edges of some given edge (a0,b0).  The length of the sub-edges
// is distributed exponentially over a large range, and the endpoints may be
// slightly perturbed to one side of (a0,b0) or the other.
    fun getPerturbedSubEdges(a0: S2Point, b0: S2Point, count: Int, edges: MutableList<TestEdge>) {
        edges.clear()
        a0.normalize()
        b0.normalize()
        val length0 = S1Angle(a0, b0)
        repeat(count) {
            val length = length0 * pow(1e-15, randomDouble())
            val offset = (length0 - length) * randomDouble()
            edges.add(TestEdge(perturbAtDistance(offset, a0, b0), perturbAtDistance(offset + length, a0, b0)))
        }
    }

    // Generate edges whose center is randomly chosen from the given S2Cap, and
// whose length is randomly chosen up to "max_length".
    fun getCapEdges(center_cap: S2Cap, max_length: S1Angle, count: Int, edges: MutableList<TestEdge>) {
        edges.clear()
        repeat(count) {
            val center = S2Random.samplePoint(center_cap)
            val edge_cap = S2Cap.fromCenterAngle(center, 0.5 * max_length)
            val p1 = S2Random.samplePoint(edge_cap)
            // Compute p1 reflected through "center", and normalize for good measure.
            val p2 = ((2.0 * p1.dotProd(center)) * center - p1).normalize()
            edges.add(TestEdge(p1, p2))
        }
    }

    // Project ShapeEdges to ShapeEdgeIds.  Useful because
// ShapeEdge does not have operator==, but ShapeEdgeId does.
    fun getShapeEdgeIds(shape_edges: List<ShapeEdge>): List<ShapeEdgeId> {
        val shape_edge_ids = mutableListOf<ShapeEdgeId>()
        for (shape_edge in shape_edges) {
            shape_edge_ids.add(shape_edge.id)
        }
        return shape_edge_ids
    }

    fun testAllCrossings(edges: List<TestEdge>) {
        val shape = S2EdgeVectorShape()  // raw pointer since "shape" used below
        for (edge in edges) {
            shape.add(edge.first, edge.second)
        }
        // Force more subdivision than usual to make the test more challenging.
        val options = AbstractMutableS2ShapeIndex.Options(maxEdgesPerCell = 1)
        val index = MutableS2ShapeIndex(options)
        val shape_id = index.add(shape)
        assertThat(shape_id).isEqualTo(0)
        // To check that candidates are being filtered reasonably, we count the
        // total number of candidates that the total number of edge pairs that
        // either intersect or are very close to intersecting.
        var num_candidates = 0
        var num_nearby_pairs = 0
        var i = 0
        for (edge in edges) {
            val a = edge.first
            val b = edge.second
            val query = S2CrossingEdgeQuery(index)
            val candidates = query.getCandidates(a, b, shape)

            // Verify that the second version of GetCandidates returns the same result.
            val edge_candidates = query.getCandidates(a, b)
            assertThat(edge_candidates).isEqualTo(candidates)
            assertThat(!candidates.isEmpty()).isTrue()

            // Now check the actual candidates.
            assertThat(candidates.isSorted()).isTrue()
            assertThat(candidates.last().shapeId).isEqualTo(0)  // Implies all shape_ids are 0.
            assertThat(candidates.first().edgeId).isGreaterThanOrEqualTo(0)
            assertThat(candidates.last().edgeId).isLessThan(shape.numEdges)
            num_candidates += candidates.size
            var missing_candidates = ""
            val expected_crossings = mutableListOf<ShapeEdgeId>()
            var expected_interior_crossings = mutableListOf<ShapeEdgeId>()
            for (i in 0 until shape.numEdges) {
                val edge = shape.edge(i)
                val c = edge.v0
                val d = edge.v1
                val sign = S2EdgeCrossings.crossingSign(a, b, c, d)
                if (sign >= 0) {
                    expected_crossings.add(ShapeEdgeId(0, i))
                    if (sign > 0) {
                        expected_interior_crossings.add(ShapeEdgeId(0, i))
                    }
                    ++num_nearby_pairs
                    if (candidates.binarySearch(ShapeEdgeId(0, i)) < 0) {
                        missing_candidates += " $i"
                    }
                } else {
                    val kMaxDist = S2Coords.projection.kMaxDiag.getValue(S2Coords.kMaxCellLevel)
                    if (S2EdgeDistances.getDistance(a, c, d).radians < kMaxDist ||
                        S2EdgeDistances.getDistance(b, c, d).radians < kMaxDist ||
                        S2EdgeDistances.getDistance(c, a, b).radians < kMaxDist ||
                        S2EdgeDistances.getDistance(d, a, b).radians < kMaxDist
                    ) {
                        ++num_nearby_pairs
                    }
                }
            }
            assertThat(missing_candidates.isEmpty())
                .withFailMessage(missing_candidates)
                .isTrue()

            // Test that GetCrossings() returns only the actual crossing edges.
            val actual_crossings = query.getCrossingEdges(a, b, shape, CrossingType.ALL)
            assertThat(getShapeEdgeIds(actual_crossings)).isEqualTo(expected_crossings)

            // Verify that the second version of GetCrossings returns the same result.
            val actual_edge_crossings = query.getCrossingEdges(a, b, CrossingType.ALL)
            assertThat(getShapeEdgeIds(actual_edge_crossings)).isEqualTo(expected_crossings)

            // Verify that CrossingType::INTERIOR returns only the interior crossings.
            val actual_interior_crossings = query.getCrossingEdges(a, b, shape, CrossingType.INTERIOR)
            assertThat(getShapeEdgeIds(actual_interior_crossings)).isEqualTo(expected_interior_crossings)
        }
        // There is nothing magical about this particular ratio; this check exists
        // to catch changes that dramatically increase the number of candidates.
        assertThat(num_candidates).isLessThanOrEqualTo(3 * num_nearby_pairs)
    }

    // Test edges that lie in the plane of one of the S2 cube edges.  Such edges
// may lie on the boundary between two cube faces, or pass through a cube
// vertex, or follow a 45 diagonal across a cube face toward its center.
//
// This test is sufficient to demonstrate that padding the cell boundaries is
// necessary for correctness.  (It fails if MutableS2ShapeIndex::kCellPadding
// is set to zero.)
    @Test
    fun perturbedCubeEdges() {
        val edges = mutableListOf<TestEdge>()
        repeat(10) {
            val face = randomInt(6)
            val scale = pow(1e-15, randomDouble())
            val uv = R2Point(2 * randomInt(2) - 1, 2 * randomInt(2) - 1)  // vertex
            val a0 = faceUvToXyz(face, scale * uv)
            val b0 = a0 - 2 * S2Coords.norm(face)
            // TODO(ericv): This test is currently slow because *every* crossing test
            // needs to invoke s2pred::ExpensiveSign().
            getPerturbedSubEdges(a0, b0, 30, edges)
            testAllCrossings(edges)
        }
    }

    // Test edges that lie in the plane of one of the S2 cube face axes.  These
// edges are special because one coordinate is zero, and they lie on the
// boundaries between the immediate child cells of the cube face.
    @Test
    fun perturbedCubeFaceAxes() {
        val edges = mutableListOf<TestEdge>()
        repeat(5) {
            val face = randomInt(6)
            val scale = pow(1e-15, randomDouble())
            val axis = S2Coords.uvwAxis(face, randomInt(2))
            val a0 = scale * axis + S2Coords.norm(face)
            val b0 = scale * axis - S2Coords.norm(face)
            getPerturbedSubEdges(a0, b0, 30, edges)
            testAllCrossings(edges)
        }
    }

    @Test
    fun capEdgesNearCubeVertex() {
        // Test a random collection of edges near the S2 cube vertex where the
        // Hilbert curve starts and ends.
        val edges = mutableListOf<TestEdge>()
        getCapEdges(
            S2Cap.fromCenterAngle(S2Point(-1, -1, 1).normalize(), S1Angle.radians(1e-3)),
            S1Angle.radians(1e-4), 1000, edges
        )
        testAllCrossings(edges)
    }

    @Test
    fun degenerateEdgeOnCellVertexIsItsOwnCandidate() {
        repeat(100) {
            val edges = mutableListOf<TestEdge>()
            val cell = S2Cell(S2Random.randomCellId())
            edges.add(TestEdge(cell.getVertex(0), cell.getVertex(0)))
            testAllCrossings(edges)
        }
    }

    @Test
    fun collinearEdgesOnCellBoundaries() {
        val kNumEdgeIntervals = 8  // 9*8/2 = 36 edges
        for (level in 0..S2CellId.kMaxLevel) {
            val cell = S2Cell(S2Random.randomCellId(level))
            val t = randomInt(4)
            val p1 = cell.getVertexRaw(t)
            val p2 = cell.getVertexRaw(t + 1)
            val delta = (p2 - p1).div(kNumEdgeIntervals.toDouble())
            val edges = mutableListOf<TestEdge>()
            for (i in 0..kNumEdgeIntervals) {
                for (j in 0 until i) {
                    edges.add(TestEdge((p1 + i * delta).normalize(), (p1 + j * delta).normalize()))
                }
            }
            testAllCrossings(edges);
        }
    }

    // This is the example from the header file, with a few extras.
    fun testPolylineCrossings(index: S2ShapeIndex, a0: S2Point, a1: S2Point) {
        val query = S2CrossingEdgeQuery(index)
        val edges = query.getCrossingEdges(a0, a1, CrossingType.ALL)
        if (edges.isEmpty()) return
        for (edge in edges) {
            assertThat(S2EdgeCrossings.crossingSign(a0, a1, edge.v0, edge.v1)).isGreaterThanOrEqualTo(0)
        }
        // Also test that no edges are missing.
        for (i in 0 until index.nextNewShapeId()) {
            val shape = index.shape(i) as S2Polyline.Shape
            val polyline = shape.polyline
            for (e in 0 until (polyline.numVertices - 1)) {
                if (S2EdgeCrossings.crossingSign(a0, a1, polyline.vertex(e), polyline.vertex(e + 1)) >= 0) {
                    assertThat(edges.count { edge -> edge.id == ShapeEdgeId(i, e) }).isOne()
                }
            }
        }
    }

    @Test
    fun PolylineCrossings() {
        val index = MutableS2ShapeIndex()
        // Three zig-zag lines near the equator.
        index.add(S2Polyline.Shape(polyline = makePolyline("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")))
        index.add(S2Polyline.Shape(polyline = makePolyline("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")))
        index.add(S2Polyline.Shape(polyline = makePolyline("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")))
        testPolylineCrossings(index, makePoint("1:0"), makePoint("1:4"))
        testPolylineCrossings(index, makePoint("5:5"), makePoint("6:6"))
    }

    @Test
    fun shapeIdsAreCorrect() {
        // This tests that when some index cells contain only one shape, the
        // intersecting edges are returned with the correct shape id.
        val index = MutableS2ShapeIndex()
        index.add(S2Polyline.Shape(polyline = S2Polyline(makeRegularPoints(makePoint("0:0"), S1Angle.degrees(5), 100))))
        index.add(
            S2Polyline.Shape(
                polyline = S2Polyline(
                    makeRegularPoints(
                        makePoint("0:20"),
                        S1Angle.degrees(5),
                        100
                    )
                )
            )
        )
        testPolylineCrossings(index, makePoint("1:-10"), makePoint("1:30"));
    }

    // Verifies that when VisitCells() is called with a specified root cell and a
// query edge that barely intersects that cell, that at least one cell is
// visited.  (At one point this was not always true, because when the query edge
// is clipped to the index cell boundary without using any padding then the
// result is sometimes empty, i.e., the query edge appears not to intersect the
// specifed root cell.  The code now uses an appropriate amount of padding,
// i.e. S2::kFaceClipErrorUVCoord.)
    @Test
    fun visitCellsQueryEdgeOnFaceBoundary() {
        val kIters = 100
        repeat(kIters) {

            // Choose an edge AB such that B is nearly on the edge between two S2 cube
            // faces, and such that the result of clipping AB to the face that nominally
            // contains B (according to S2::GetFace) is empty when no padding is used.
            var a_face: Int
            var b_face: Int
            var a: S2Point
            var b: S2Point
            var a_uv: R2Point
            var b_uv: R2Point
            do {
                a_face = randomInt(6)
                a = faceUvToXyz(a_face, randomDouble(-1.0, 1.0), randomDouble(-1.0, 1.0)).normalize()
                b_face = S2Coords.uvwFace(a_face, 0, 1)  // Towards positive u-axis
                b = faceUvToXyz(
                    b_face,
                    1 - randomInt(2) * 0.5 * DoubleType.epsilon,
                    randomDouble(-1.0, 1.0)
                ).normalize()
                val clippedEdge = S2EdgeClipping.clipToFace(a, b, b_face)
            } while (S2Coords.face(b) != b_face || clippedEdge != null)

            // Verify that the clipping result is non-empty when a padding of
            // S2::kFaceClipErrorUVCoord is used instead.
            val clipped = S2EdgeClipping.clipToPaddedFace(a, b, b_face, S2EdgeClipping.kFaceClipErrorUVCoord)
            assertThat(clipped).isNotNull()

            // Create an S2ShapeIndex containing a single edge BC, where C is on the
            // same S2 cube face as B (which is different than the face containing A).
            val c = faceUvToXyz(b_face, randomDouble(-1.0, 1.0), randomDouble(-1.0, 1.0)).normalize()
            val index = MutableS2ShapeIndex()
            index.add(S2Polyline.Shape(polyline = S2Polyline(vertices = listOf(b, c))))

            // Check that the intersection between AB and BC is detected when the face
            // containing BC is specified as a root cell.  (Note that VisitCells()
            // returns false only if the CellVisitor returns false, and otherwise
            // returns true.)
            val query = S2CrossingEdgeQuery(index)
            val root = S2PaddedCell(S2CellId.fromFace(b_face), 0.0)
            assertThat(query.visitCells(a, b, root, object : CellVisitor {
                override fun visit(cell: S2ShapeIndexCell): Boolean {
                    return false
                }
            })).isFalse()
        }
    }

}  // namespace
