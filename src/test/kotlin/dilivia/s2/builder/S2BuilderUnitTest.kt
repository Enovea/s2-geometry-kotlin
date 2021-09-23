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
@file:Suppress("UsePropertyAccessSyntax")

package dilivia.s2.builder

import com.google.common.base.Stopwatch
import dilivia.math.matrix.Matrix3x3Double
import dilivia.s2.Fractal
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Debug
import dilivia.s2.S2Earth
import dilivia.s2.S2Error
import dilivia.s2.S2Factory.makePolyline
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomFrame
import dilivia.s2.S2TextParser
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.S2TextParser.makePolygon
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.layers.Layer
import dilivia.s2.builder.layers.S2PolygonLayer
import dilivia.s2.builder.layers.S2PolylineLayer
import dilivia.s2.builder.layers.S2PolylineVectorLayer
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.builder.snap.S2CellIdSnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.edge.S2EdgeCrossings.kIntersectionError
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Loop.Companion.makeRegularLoop
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.util.concurrent.TimeUnit
import kotlin.math.min
import kotlin.math.pow

// Iteration multiplier for randomized tests
const val iteration_multiplier = 1

// A set of (edge string, vector<InputEdgeId>) pairs representing the
// InputEdgeIds attached to the edges of a graph.  Edges are in
// S2TextParser.toString() format, such as "1:3, 4:5".
typealias EdgeInputEdgeIds = List<Pair<String, List<Int>>>

class S2BuilderUnitTest {

    val logger = KotlinLogging.logger {}

    private fun expectPolygonsEqual(expected: S2Polygon, actual: S2Polygon) {
        assertThat(actual)
            .withFailMessage(
                """
      |
      |Expected:
      |${S2TextParser.toString(expected)}
      |Actual:
      |${S2TextParser.toString(actual)}
    """.trimMargin()
            )
            .isEqualTo(expected)
    }

    private fun expectPolygonsApproxEqual(expected: S2Polygon, actual: S2Polygon, tolerance: S1Angle) {
        assertThat(expected.boundaryApproxEquals(actual, tolerance))
            .withFailMessage(
                """
      |
      |Expected:  ${S2TextParser.toString(expected)}
      |Actual:    ${S2TextParser.toString(actual)}
      |Tolerance: ${tolerance.degrees()}
      """.trimMargin()
            )
            .isTrue()
    }

    private fun expectPolylinesEqual(expected: S2Polyline, actual: S2Polyline) {
        assertThat(actual)
            .withFailMessage(
                """
      |
      |Expected:  ${S2TextParser.toString(expected)}
      |Actual:    ${S2TextParser.toString(actual)}
    """.trimMargin()
            )
            .isEqualTo(expected)
    }

    @Test
    fun addShape() {
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        val input = makePolygon("0:0, 0:5, 5:5, 5:0; 1:1, 1:4, 4:4, 4:1")
        builder.addShape(input.index.shape(0)!!)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        expectPolygonsEqual(input, output)
    }

    @Test
    fun simpleVertexMerging() {
        // When IdentitySnapFunction is used (i.e., no special requirements on
        // vertex locations), check that vertices closer together than the snap
        // radius are merged together.

        val snapRadius = S1Angle.degrees(0.5)
        val builder = S2Builder(S2Builder.Options(IdentitySnapFunction(snapRadius)))
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        val input = makePolygon("0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9")
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = makePolygon("0:0, 0:1, 1:0.9")
        expectPolygonsApproxEqual(expected, output, snapRadius)
    }

    @Test
    fun simpleS2CellIdSnapping() {
        // When S2CellIdSnapFunction is used, check that all output vertices are the
        // centers of S2CellIds at the specified level level.

        val level = S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.degrees(1))
        val snapFunction = S2CellIdSnapFunction(level)
        val builder = S2Builder(S2Builder.Options(snapFunction))
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        val input = makePolygon("2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3")
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output.numLoops()).isEqualTo(1)
        val loop = output.loop(0)
        for (i in 0 until loop.numVertices) {
            assertThat(loop.vertex(i)).isEqualTo(S2CellId.fromPoint(loop.vertex(i)).parent(level).toPoint())
        }
        expectPolygonsApproxEqual(input, output, snapFunction.snapRadius)
    }

    @Test
    fun simpleIntLatLngSnapping() {
        val builder = S2Builder(S2Builder.Options(IntLatLngSnapFunction(0)))  // E0 coords
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        val input =
            makePolygon("2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, 5.22:3.88, 5.55:2.49, 4.49:2.51")
        val expected = makePolygon("2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3")
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output.numLoops()).isEqualTo(1)
        expectPolygonsEqual(expected, output)
    }

    @Test
    fun verticesMoveLessThanSnapRadius() {
        // Check that chains of closely spaced vertices do not collapse into a
        // single vertex.

        val snapRadius = S1Angle.degrees(1)
        val builder = S2Builder(S2Builder.Options(IdentitySnapFunction(snapRadius)))
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        // The spacing between input vertices is about 2*pi*20/1000 = 0.125 degrees.
        // The output vertices are spaced between 1 and 2 degrees apart; the average
        // spacing is about 1.33 degrees.
        val input = S2Polygon(makeRegularLoop(S2Point(1, 0, 0), S1Angle.degrees(20), 1000))
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output.numLoops()).isEqualTo(1)
        assertThat(output.loop(0).numVertices).isGreaterThanOrEqualTo(90).isLessThanOrEqualTo(100)
        assertThat(output.boundaryNear(input, snapRadius)).isTrue()
    }

    @Test
    fun minEdgeVertexSeparation() {
        // Check that edges are separted from non-incident vertices by at least
        // min_edge_vertex_separation().  This requires adding new vertices (not
        // present in the input) in some cases.

        // The input is a skinny right triangle with two legs of length 10 and 1,
        // and whose diagonal is subdivided into 10 short edges.  Using a snap
        // radius of 0.5, about half of the long leg is snapped onto the diagonal
        // (which causes that part of the polygon to be removed).  But the real
        // problem is that the remaining part of the long leg gets too close to the
        // remaining vertices on the diagonal, i.e. it would violate the minimum
        // edge-vertex separation guarantee.  S2Builder handles this by creating at
        // least one vertex along the original long leg, to keep the snapped edge
        // far enough away from the diagonal.
        val input = makePolygon("0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0")
        val expected = makePolygon("0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 4.00021862252687:0")
        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(0.5)))
        val builder = S2Builder(options)
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        expectPolygonsApproxEqual(expected, output, S1Angle.radians(1e-15))
    }

    @Test
    fun idempotencySnapsInadequatelySeparatedVertices() {
        // This test checks that when vertices are closer together than
        // min_vertex_separation() then they are snapped together even when
        // options.idempotent() is true.
        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(1.0)))
        val builder = S2Builder(options)
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output))
        builder.addPolyline(makePolyline("0:0, 0:0.9, 0:2"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = "0:0, 0:2"
        assertThat(S2TextParser.toString(output)).isEqualTo(expected)
    }

    @Test
    fun idempotencySnapsIdenticalVerticesWithZeroSnapRadius() {
        // This test checks that even when the snap radius is zero, identical
        // vertices are snapped together.
        val builder = S2Builder(S2Builder.Options())
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        builder.addPolyline(makePolyline("0:1, 1:0"))
        builder.addPolyline(makePolyline("0:0, 0:1"))
        builder.addEdge(makePoint("0:1"), makePoint("0:1"))
        builder.addPolyline(makePolyline("1:0, 0:0"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = "0:0, 0:1, 1:0"
        assertThat(S2TextParser.toString(output)).isEqualTo(expected)
    }

    @Test
    fun idempotencySnapsIdenticalVerticesWithZeroSnapRadiusEdgeSplitting() {
        // This test checks that identical vertices are snapped together even when
        // the snap radius is zero and options.split_crossing_edges() is true.
        val options = S2Builder.Options()
        options.splitCrossingEdges = true
        val builder = S2Builder(options)
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        builder.addPolyline(makePolyline("0:1, 1:0"))
        builder.addPolyline(makePolyline("0:0, 0:1"))
        builder.addEdge(makePoint("0:1"), makePoint("0:1"))
        builder.addPolyline(makePolyline("1:0, 0:0"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = "0:0, 0:1, 1:0"
        assertThat(S2TextParser.toString(output)).isEqualTo(expected)
    }

    @Test
    fun idempotencySnapsUnsnappedVertices() {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that S2Builder detects vertices
        // that are not at a valid location returned by the given snap function.

        // In this example we snap two vertices to integer lat/lng coordinates.  The
        // two vertices are far enough apart (more than min_vertex_separation) so
        // that they might be the result of a previous snapping operation, but one
        // of the two vertices does not have integer lat/lng coordinates.  We use
        // internal knowledge of how snap sites are chosen (namely, that candidates
        // are considered in S2CellId order) to construct two different cases, one
        // where the snapped vertex is processed first and one where the unsnapped
        // vertex is processed first.  This exercises two different code paths.
        val snapFunction = IntLatLngSnapFunction(0)
        assertThat(snapFunction.snapRadius).isGreaterThanOrEqualTo(S1Angle.degrees(0.7))
        assertThat(snapFunction.minVertexSeparation()).isLessThanOrEqualTo(S1Angle.degrees(0.35))
        val builder = S2Builder(S2Builder.Options(snapFunction))

        // In this example, the snapped vertex (0, 0) is processed first and is
        // selected as a Voronoi site (i.e., output vertex).  The second vertex is
        // closer than min_(), therefore it is snapped to the first vertex
        // and the polyline becomes degenerate.
        val a = S2LatLng.fromDegrees(0, 0).toPoint()
        val b = S2LatLng.fromDegrees(0.01, 0.6).toPoint()
        assertThat(S2CellId.fromPoint(a)).isLessThan(S2CellId.fromPoint(b))
        val input1 = S2Polyline(listOf(a, b))
        val output1 = S2Polyline()
        builder.startLayer(S2PolylineLayer(output1))
        builder.addPolyline(input1)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output1)).isEqualTo("0:0, 0:1")

        // In this example the unsnapped vertex is processed first and is snapped to
        // (0, 0).  The second vertex is further than snap_radius() away, so it is
        // also snapped (which does nothing) and is left at (0, 1).
        val c = S2LatLng.fromDegrees(0.01, 0.4).toPoint()
        val d = S2LatLng.fromDegrees(0, 1).toPoint()
        assertThat(S2CellId.fromPoint(c)).isLessThan(S2CellId.fromPoint(d))
        val input2 = S2Polyline(listOf(c, d))
        val output2 = S2Polyline()
        builder.startLayer(S2PolylineLayer(output2))
        builder.addPolyline(input2)
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output2)).isEqualTo("0:0, 0:1")
    }

    @Test
    fun idempotencySnapsEdgesWithTinySnapRadius() {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that S2Builder detects edges that
        // are too close to vertices even when the snap radius is very small
        // (e.g., S2.kIntersectionError).
        //
        // Previously S2Builder used a conservative approximation to decide whether
        // edges were too close to vertices; unfortunately this meant that when the
        // snap radius was very small then no snapping would be done at all, because
        // even an edge/vertex distance of zero was considered far enough apart.
        //
        // This tests that the current code (which uses exact predicates) handles
        // this situation correctly (i.e., that an edge separated from a
        // non-incident vertex by a distance of zero cannot be the output of a
        // previous snapping operation).
        val options = S2Builder.Options()
        options.snapFunction = IdentitySnapFunction(kIntersectionError)
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.duplicateEdges = DuplicateEdges.MERGE
        val builder = S2Builder(options)
        val output = mutableListOf<S2Polyline>()
        builder.startLayer(S2PolylineVectorLayer(output, options = layerOptions))
        builder.addPolyline(makePolyline("0:0, 0:10"))
        builder.addPolyline(makePolyline("0:5, 0:7"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output.size).isEqualTo(1)
        assertThat(S2TextParser.toString(output[0])).isEqualTo("0:0, 0:5, 0:7, 0:10")
    }

    @Test
    fun idempotencyDoesNotSnapAdequatelySeparatedEdges() {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that when an edge is further away
        // than min_edge_vertex_separation() then no snapping is done.
        val options = S2Builder.Options(IntLatLngSnapFunction(0))
        options.idempotent = true  // Test fails if this is "false".
        val builder = S2Builder(options)
        val output1 = S2Polygon()
        val output2 = S2Polygon()
        builder.startLayer(S2PolygonLayer(output1))
        builder.addPolygon(makePolygon("1.49:0, 0:2, 0.49:3"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected = "1:0, 0:2, 0:3"
        assertThat(S2TextParser.toString(output1)).isEqualTo(expected)
        builder.startLayer(S2PolygonLayer(output2))
        builder.addPolygon(output1)
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output2)).isEqualTo(expected)
    }

    @Test
    fun kMaxSnapRadiusCanSnapAtLevel0() {
        // Verify that kMaxSnapRadius will allow snapping at S2CellId level 0.
        assertThat(S2CellIdSnapFunction.minSnapRadiusForLevel(0)).isLessThanOrEqualTo(SnapFunction.kMaxSnapRadius())
    }

    //@Disabled("infinite loop")
    @Test
    fun s2CellIdSnappingAtAllLevels() {
        val input = makePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0")
        for (level in 0..S2CellId.kMaxLevel) {
            val snapFunction = S2CellIdSnapFunction(level)
            val builder = S2Builder(S2Builder.Options(snapFunction))
            val output = S2Polygon()
            builder.startLayer(S2PolygonLayer(output))
            builder.addPolygon(input)
            val error = S2Error()
            assertThat(builder.build(error)).isTrue()
            assertThat(output.isValid()).isTrue()
            // The ApproxContains calls below are not guaranteed to succeed in general
            // because ApproxContains works by snapping both polygons together using
            // the given tolerance and then checking for containment.  Since
            // ApproxContains snaps to an arbitrary subset of the input vertices
            // rather than to S2CellId centers at the current level, this means that
            // corresponding vertices in "input" and "output" can snap to different
            // sites, which causes the containment test to fail.  Nevertheless, by
            // using a larger tolerance of 2 * snap_radius, all calls in this test
            // succeed (and would be likely to succeed in other similar tests).
            // (To guarantee correctness we would need to use S2CellIdSnapFunction
            // within the ApproxContains implementation.)
            val tolerance = minOf(2.times(snapFunction.snapRadius), SnapFunction.kMaxSnapRadius())
            assertThat(output.approxContains(input, tolerance)).isTrue()
            assertThat(input.approxContains(output, tolerance)).isTrue()
        }
    }

    @Test
    fun snappingDoesNotRotateVertices() {
        // This is already tested extensively elsewhere.
        val input = makePolygon(
            "49.9305505:-124.8345463, 49.9307448:-124.8299657, 49.9332101:-124.8301996, 49.9331224:-124.8341368; " +
                    "49.9311087:-124.8327042, 49.9318176:-124.8312621, 49.9318866:-124.8334451"
        )
        val options = S2Builder.Options((S2CellIdSnapFunction()))
        val builder = S2Builder(options)
        val output1 = S2Polygon()
        val output2 = S2Polygon()
        builder.startLayer(S2PolygonLayer(output1))
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        // This checks that the vertices are in the same cyclic order, and that
        // vertices have not moved by more than "snap_radius".
        expectPolygonsApproxEqual(input, output1, options.snapFunction.snapRadius())

        // Check that snapping twice doesn't rotate the vertices.  This also
        // verifies that S2Builder can be used again after Build() is called.
        builder.startLayer(S2PolygonLayer(output2))
        builder.addPolygon(output1)
        assertThat(builder.build(error)).isTrue()
        expectPolygonsEqual(output1, output2)
    }

    @Test
    fun selfIntersectingPolyline() {
        // Check that when two edges of a polyline cross, the intersection point is
        // added to both edges.

        val options = S2Builder.Options()
        val snapFunction = IntLatLngSnapFunction(1)  // Snap to E1 coordinates
        options.snapFunction = snapFunction
        options.splitCrossingEdges = true
        val builder = S2Builder(options)
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output))
        val input = makePolyline("3:1, 1:3, 1:1, 3:3")
        val expected = makePolyline("3:1, 2:2, 1:3, 1:1, 2:2, 3:3")
        builder.addPolyline(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        expectPolylinesEqual(expected, output)
    }

    @Test
    fun selfIntersectingPolygon() {
        // Check that when two edge of a polygon cross, the intersection point is
        // added to both edges, and that the resulting (undirected) edges can be
        // assembled into a valid polygon.

        val snapFunction = IntLatLngSnapFunction(1)  // Snap to E1 coordinates
        val options = S2Builder.Options()
        options.snapFunction = snapFunction
        options.splitCrossingEdges = true
        val builder = S2Builder(options)
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output, options = S2PolygonLayer.Options(EdgeType.UNDIRECTED)))
        val input = makePolyline("3:1, 1:3, 1:1, 3:3, 3:1")
        val expected = makePolygon("1:1, 1:3, 2:2; 3:3, 3:1, 2:2")
        builder.addPolyline(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        expectPolygonsEqual(expected, output)
    }

    @Test
    fun tieBreakingIsConsistent() {
        // Check that when an edge passes between two equally distant vertices, that
        // the choice of which one to snap to does not depend on the edge direction.

        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(2)))
        options.idempotent = false
        val builder = S2Builder(options)
        builder.forceVertex(S2LatLng.fromDegrees(1, 0).toPoint())
        builder.forceVertex(S2LatLng.fromDegrees(-1, 0).toPoint())
        val output1 = S2Polyline()
        val output2 = S2Polyline()
        builder.startLayer(S2PolylineLayer(output1))
        builder.addPolyline(makePolyline("0:-5, 0:5"))
        builder.startLayer(S2PolylineLayer(output2))
        builder.addPolyline(makePolyline("0:5, 0:-5"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output1.numVertices).isEqualTo(3)
        assertThat(output2.numVertices).isEqualTo(3)
        for (i in 0 until 3) {
            assertThat(output2.vertex(2 - i)).isEqualTo(output1.vertex(i))
        }
    }

    // Verifies that two graphs have the same vertices and edges.
    private fun expectGraphsEqual(expected: Graph, actual: Graph) {
        assertThat(actual.vertices).isEqualTo(expected.vertices)
        assertThat(actual.edges).isEqualTo(expected.edges)
        assertThat(actual.inputEdgeIdSetIds).isEqualTo(expected.inputEdgeIdSetIds)
    }

    // This layer makes both a shallow and a deep copy of the Graph object passed
// to its Build() method and appends them to two vectors.  Furthermore, it
// verifies that the shallow and deep copies of any graphs previously appended
// to those vectors are still identical.
    inner class GraphPersistenceLayer(private val graphOptions: GraphOptions) : Layer() {

        private val shallowCopies: MutableList<Graph> = mutableListOf()
        private val deepCopies: MutableList<Graph> = mutableListOf()

        override fun graphOptions(): GraphOptions = graphOptions

        override fun build(g: Graph, error: S2Error) {
            for (i in shallowCopies.indices) {
                expectGraphsEqual(shallowCopies[i], deepCopies[i])
            }
            shallowCopies.add(g)
            deepCopies.add(Graph(g))
        }
    }

    @Test
    fun graphPersistence() {
        // Ensure that the Graph objects passed to S2Builder.Layer.Build() methods
        // remain valid until all layers have been built.
        val builder = S2Builder(S2Builder.Options())
        repeat(20) {
            builder.startLayer(GraphPersistenceLayer(GraphOptions()))
            repeat(S2Random.randomInt(10)) {
                builder.addEdge(S2Random.randomPoint(), S2Random.randomPoint())
            }
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
    }

    private fun testPolylineLayers(
        input_strs: List<String>,
        expected_strs: List<String>,
        layer_options: S2PolylineLayer.Options,
        builder_options: S2Builder.Options = S2Builder.Options()
    ) {
        val builder = S2Builder(builder_options)
        val output = mutableListOf<S2Polyline>()
        for (input_str in input_strs) {
            output.add(S2Polyline())
            builder.startLayer(S2PolylineLayer(output.last(), options = layer_options))
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val outputStrs = mutableListOf<String>()
        for (polyline in output) {
            outputStrs.add(S2TextParser.toString(polyline))
        }
        assertThat(outputStrs.joinToString("; ")).isEqualTo(expected_strs.joinToString("; "))
    }

    private fun testPolylineVector(
        input_strs: List<String>,
        expected_strs: List<String>,
        layer_options: S2PolylineVectorLayer.Options,
        builder_options: S2Builder.Options = S2Builder.Options()
    ) {
        val builder = S2Builder(builder_options)
        val output = mutableListOf<S2Polyline>()
        builder.startLayer(S2PolylineVectorLayer(output, options = layer_options))
        for (input_str in input_strs) {
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val outputStrs = mutableListOf<String>()
        for (polyline in output) {
            outputStrs.add(S2TextParser.toString(polyline))
        }
        assertThat(outputStrs.joinToString("; ")).isEqualTo(expected_strs.joinToString("; "))
    }

    private fun testPolylineLayersBothEdgeTypes(
        input_strs: List<String>,
        expected_strs: List<String>,
        layer_options: S2PolylineLayer.Options,  // by value
        builder_options: S2Builder.Options = S2Builder.Options()
    ) {
        layer_options.edgeType = EdgeType.DIRECTED
        testPolylineLayers(input_strs, expected_strs, layer_options, builder_options)
        layer_options.edgeType = EdgeType.UNDIRECTED
        testPolylineLayers(input_strs, expected_strs, layer_options, builder_options)
    }

    @Test
    fun simplifyOneEdge() {
        // Simplify a perturbed edge chain into a single edge.

        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(1)))
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf("0:0, 1:0.5, 2:-0.5, 3:0.5, 4:-0.5, 5:0"),
            listOf("0:0, 5:0"),
            S2PolylineLayer.Options(), options
        )
    }

    @Test
    fun simplifyTwoLayers() {
        // Construct two layers, each containing a polyline that could be simplified
        // to a single edge on its own.  However the two polylines actually cross,
        // so make sure that the output still contains the intersection vertex.

        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(0.5)))
        options.splitCrossingEdges = true
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf("-2:-1, -1:0, 1:0, 2:1", "1:-2, 0:-1, 0:1, -1:2"),
            listOf("-2:-1, 0:0, 2:1", "1:-2, 0:0, -1:2"),
            S2PolylineLayer.Options(), options
        )
    }

    @Test
    fun simplifyOneLoop() {
        // Simplify a regular loop with 1000 vertices and a radius of 20 degrees.
        // Turning on edge chain simplification yields a dramatically smaller number
        // of vertices than snapping alone (10 vertices vs 95 vertices using a snap
        // radius of 1 degree).  This is because snapping alone yields vertices that
        // stay within 1 degree of the input *vertices*, while simplifying edge
        // chains yields edges that stay within 1 degree of the input *edges*.

        for (i in 0..1) {
            val edgeType = EdgeType.values()[i]
            val snapRadius = S1Angle.degrees(1)
            val options = S2Builder.Options((IdentitySnapFunction(snapRadius)))
            options.setSimplifyEdgeChains(true)
            val builder = S2Builder(options)
            val output = S2Polygon()
            builder.startLayer(S2PolygonLayer(output, options = S2PolygonLayer.Options(edgeType)))
            // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
            val input = S2Polygon(makeRegularLoop(S2Point(1, 0, 0), S1Angle.degrees(20), 1000))
            builder.addPolygon(input)
            val error = S2Error()
            assertThat(builder.build(error)).isTrue()
            assertThat(output.numLoops()).isEqualTo(1)

            logger.debug { "initial: ${S2TextParser.toString(input)}" }
            logger.debug { "output: ${S2TextParser.toString(output)}" }
            logger.debug { "num vertices: ${output.loop(0).numVertices}" }

            assertThat(output.loop(0).numVertices).isGreaterThanOrEqualTo(10).isLessThanOrEqualTo(12)
            assertThat(output.boundaryNear(input, snapRadius)).isTrue()
        }
    }

    @Test
    fun simplifyOppositeDirections() {
        // We build two layers with two polylines that follow the same circular arc
        // in opposite directions, and verify that they are snapped identically.
        // (The snap radius is adjusted so that the arc is simplified into a long
        // edge and a short edge, and therefore we would get a different result if
        // the two layers followed the edge chain in different directions.)

        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(0.5)))
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf(
                "-4:0.83, -3:0.46, -2:0.2, -1:0.05, 0:0, 1:0.5, 2:0.2, 3:0.46, 4:0.83",
                "4:.83, 3:.46, 2:.2, 1:.05, 0:0, -1:.5, -2:.2, -3:.46, -4:.83"
            ),
            listOf("-4:0.83, -2:0.2, 4:0.83", "4:0.83, -2:0.2, -4:0.83"),
            S2PolylineLayer.Options(), options
        )
    }

    @Test
    fun simplifyKeepsEdgeVertexSeparation() {
        // We build two layers each containing a polyline, such that the polyline in
        // the first layer could be simplified to a straight line except that then
        // it would create an intersection with the second polyline.

        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(1.0)))
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf("0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"),
            listOf("0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"),
            S2PolylineLayer.Options(), options
        )
    }

    @Test
    fun simplifyBacktrackingEdgeChain() {
        // Test simplifying an edge chain that backtracks on itself.
        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.degrees(0.5)))
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf("0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 4:0, 3:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0"),
            listOf("0:0, 2:0, 5:0, 2:0, 5:0, 7:0"),
            S2PolylineLayer.Options(), options
        )
    }

    @Test
    fun simplifyLimitsEdgeDeviation() {
        // Make sure that simplification does not create long edges such that the
        // midpoint of the edge might be further than max_edge_deviation() from an
        // input edge.  In the example below, vertices are snapped to integer
        // lat/lng coordinates, and the snap radius is approximately 0.707 degrees.
        // Snapping moves the input vertices perpendicular to the input edge by just
        // slightly less than the snap radius (0.693 degrees).  Now the midpoint of
        // the snapped edge is about 0.98 degrees from the input edge, which causes
        // an extra site to be added at the midpoint of the original edge.
        //
        // When simplify_edge_chains() is enabled, then usually an extra site like
        // this would be simplified away (because the simplified edge would still be
        // within snap_radius() of all the input vertices) except that there is an
        // explicit check in S2Builder that prevents this.  (If the check is removed
        // then this test fails.)

        val options = S2Builder.Options(IntLatLngSnapFunction(0))  // E0 coordinates
        options.setSimplifyEdgeChains(true)
        testPolylineLayersBothEdgeTypes(
            listOf("-30.49:-29.51, 29.51:30.49"), listOf("-29.99999999999999:-30, -1:1, 29.99999999999999:30"),
            S2PolylineLayer.Options(), options
        )
    }

    //@Disabled("Infinite loop")
    @Test
    fun simplifyPreservesTopology() {
        // Crate several nested concentric loops, and verify that the loops are
        // still nested after simplification.

        val kNumLoops = 5 // TODO: 20
        val kNumVerticesPerLoop = 1000
        val kBaseRadius = S1Angle.degrees(5)
        val kSnapRadius = S1Angle.degrees(0.1)
        val options = S2Builder.Options((IdentitySnapFunction(kSnapRadius)))
        options.setSimplifyEdgeChains(true)
        val builder = S2Builder(options)
        val input = mutableListOf<S2Polygon>()
        val output = mutableListOf<S2Polygon>()
        repeat(kNumLoops) { j ->
            // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
            val radius = kBaseRadius + 0.7 * j * j / kNumLoops * kSnapRadius
            input.add(S2Polygon(makeRegularLoop(S2Point(1, 0, 0), radius, kNumVerticesPerLoop)))
            output.add(S2Polygon())
            builder.startLayer(S2PolygonLayer(output.last()))
            builder.addPolygon(input.last())
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        repeat(kNumLoops) { j ->
            assertThat(output[j].boundaryNear(input[j], kSnapRadius)).isTrue()
            if (j > 0) assertThat(output[j].contains(output[j - 1])).isTrue()
        }
    }

    @Test
    fun simplifyRemovesSiblingPairs() {
        val options = S2Builder.Options(IntLatLngSnapFunction(0))  // E0 coords
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.setSiblingPairs(SiblingPairs.DISCARD)

        // Check that there is no sibling pair without simplification.
        testPolylineVector(
            listOf("0:0, 0:10", "0:10, 0.6:5, 0:0"),
            listOf("0:0, 0:10, 1:5, 0:0"), layerOptions, options
        )

        // Now check that (1) simplification produces a sibling pair,
        // and (2) the sibling pair is removed (since we requested it).
        options.setSimplifyEdgeChains(true)
        testPolylineVector(
            listOf("0:0, 0:10", "0:10, 0.6:5, 0:0"),
            emptyList(), layerOptions, options
        )
    }

    @Test
    fun simplifyMergesDuplicateEdges() {
        val options = S2Builder.Options(IntLatLngSnapFunction(0))  // E0 coords
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.duplicateEdges = DuplicateEdges.MERGE

        // Check that there are no duplicate edges without simplification.
        testPolylineVector(
            listOf("0:0, 0:10", "0:0, 0.6:5, 0:10"),
            listOf("0:0, 0:10", "0:0, 1:5, 0:10"), layerOptions, options
        )

        // Now check that (1) simplification produces a duplicate edge pair,
        // and (2) the duplicate pair is merged (since we requested it).
        options.setSimplifyEdgeChains(true)
        testPolylineVector(
            listOf("0:0, 0:10", "0:0, 0.6:5, 0:10"),
            listOf("0:0, 0:10"), layerOptions, options
        )
    }

    @Test
    fun simplifyKeepsForcedVertices() {
        val options = S2Builder.Options(IdentitySnapFunction(S1Angle.radians(1e-15)))
        options.setSimplifyEdgeChains(true)
        val builder = S2Builder(options)
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output))
        builder.addPolyline(makePolyline("0:0, 0:1, 0:2, 0:3"))
        builder.forceVertex(makePoint("0:1"))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(S2TextParser.toString(output)).isEqualTo("0:0, 0:1, 0:3")
    }


    class InputEdgeIdCheckingLayer(val expected: EdgeInputEdgeIds, private val graphOptions: GraphOptions) : Layer() {

        private val INPUT_EDGE_ID_MISMATCH = S2Error.USER_DEFINED_START

        override fun graphOptions(): GraphOptions = graphOptions

        override fun build(g: Graph, error: S2Error) {
            val actual = mutableListOf<Pair<String, List<Int>>>()
            val vertices = mutableListOf<S2Point>()
            for (e in 0 until g.numEdges) {
                vertices.clear()
                vertices.add(g.vertex(g.edge(e).first))
                vertices.add(g.vertex(g.edge(e).second))
                val edge = S2TextParser.toString(
                    listOf(
                        g.vertex(g.edge(e).first),
                        g.vertex(g.edge(e).second)
                    )
                )
                val ids = g.inputEdgeIds(e)
                actual.add(Pair(edge, ids.toList()))
            }
            // This comparison doesn't consider multiplicity, but that's fine.
            var missing = ""
            var extra = ""
            for (p in expected) {
                if (actual.filter { e -> e == p }.count() > 0) continue
                missing += toString(p)
            }
            for (p in actual) {
                if (expected.filter { a -> a == p }.count() > 0) continue
                extra += toString(p)
            }
            if (missing.isNotEmpty() || extra.isNotEmpty()) {
                error.init(INPUT_EDGE_ID_MISMATCH, String.format("Missing:\n%sExtra:\n%s\n", missing, extra))
            }
        }

        fun toString(p: Pair<String, List<Int>>): String {
            var r = "  (${p.first})={"
            if (p.second.isNotEmpty()) {
                for (id in p.second) {
                    r += "$id,"
                }
                r = r.removeSuffix(", ")
            }
            r += "}\n"
            return r
        }


    }

    private fun testInputEdgeIds(
        input_strs: List<String>,
        expected: EdgeInputEdgeIds,
        graph_options: GraphOptions,
        options: S2Builder.Options
    ) {
        val builder = S2Builder(options)
        builder.startLayer(InputEdgeIdCheckingLayer(expected, graph_options))
        for (input_str in input_strs) {
            builder.addPolyline(makePolyline(input_str))
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
    }

    @Test
    fun inputEdgeIdAssignment() {
        // Check that input edge ids are assigned in order.
        testInputEdgeIds(
            listOf("0:0, 0:1, 0:2"),
            listOf(Pair("0:0, 0:1", listOf(0)), Pair("0:1, 0:2", listOf(1))), GraphOptions(), S2Builder.Options()
        )
    }

    @Test
    fun undirectedSiblingsDontHaveInputEdgeIds() {
        // Check that the siblings of undirected edges do not have InputEdgeIds.
        val graphOptions = GraphOptions()
        graphOptions.edgeType = EdgeType.UNDIRECTED
        testInputEdgeIds(
            listOf("0:0, 0:1, 0:2"),
            listOf(
                Pair("0:0, 0:1", listOf(0)), Pair("0:1, 0:2", listOf(1)),
                Pair("0:1, 0:0", emptyList()), Pair("0:2, 0:1", emptyList())
            ),
            graphOptions, S2Builder.Options()
        )
    }

    @Test
    fun createdSiblingsDontHaveInputEdgeIds() {
        // Check that edges created by SiblingPairs.CREATE do not have
        // InputEdgeIds.
        val graphOptions = GraphOptions()
        graphOptions.siblingPairs = SiblingPairs.CREATE
        testInputEdgeIds(
            listOf("0:0, 0:1, 0:2"), listOf(Pair("0:0, 0:1", listOf(0)), Pair("0:1, 0:2", listOf(1))),
            GraphOptions(), S2Builder.Options()
        )
    }

    @Test
    fun edgeMergingDirected() {
        // Tests that input edge ids are merged when directed edges are merged.
        val graphOptions = GraphOptions()
        graphOptions.duplicateEdges = DuplicateEdges.MERGE
        testInputEdgeIds(
            listOf("0:0, 0:1", "0:0, 0:1"), listOf(Pair("0:0, 0:1", listOf(0, 1))),
            graphOptions, S2Builder.Options()
        )
    }

    @Test
    fun edgeMergingUndirected() {
        // Tests that input edge ids are merged when undirected edges are merged.
        val graphOptions = GraphOptions()
        graphOptions.duplicateEdges = DuplicateEdges.MERGE
        graphOptions.siblingPairs = SiblingPairs.KEEP
        testInputEdgeIds(
            listOf("0:0, 0:1, 0:2", "0:0, 0:1", "0:2, 0:1"),
            listOf(Pair("0:0, 0:1", listOf(0, 2)), Pair("0:1, 0:2", listOf(1)), Pair("0:2, 0:1", listOf(3))),
            graphOptions, S2Builder.Options()
        )
    }

    @Test
    fun simplifyDegenerateEdgeMergingEasy() {
        // Check that when an input edge is snapped to a chain that includes
        // degenerate edges, and the edge chain is simplified, that the InputEdgeIds
        // attached to those degenerate edges are transferred to the simplified
        // edge.  For example (using integers for vertices), an edge chain 1.2,
        // 2.2, 2.3 that is simplified to 1.3 should get the InputEdgeIds
        // associated with all three original edges.  (This ensures that the labels
        // attached to those edges are also transferred.)
        //
        // This also tests that degenerate edges at the start and end of the
        // simplified chain are *not* merged.  (It's up to the output layer to
        // decide what to do with these edges.  The only reason we merge degenerate
        // edges in the interior of the interior of the simplified edge is because
        // those edges are being removed from the graph.)
        val graphOptions = GraphOptions()
        graphOptions.degenerateEdges = DegenerateEdges.KEEP
        val options = S2Builder.Options(IntLatLngSnapFunction(0))
        options.setSimplifyEdgeChains(true)
        testInputEdgeIds(
            listOf("0:0, 0:0.1, 0:1.1, 0:1, 0:0.9, 0:2, 0:2.1"), listOf(
                Pair("0:0, 0:0", listOf(0)), Pair("0:0, 0:2", listOf(1, 2, 3, 4)), Pair("0:2, 0:2", listOf(5))
            ), graphOptions, options
        )
    }

    @Test
    fun simplifyDegenerateEdgeMergingHard() {
        // This is a harder version of the test above.  Now there are several edge
        // chains that overlap each other in both directions, and several degenerate
        // edges at that middle vertex.  This tests that if exactly one edge chain
        // contains a degenerate edge in input edge order (e.g., the input order was
        // AB, BB, BC), then the degenerate edge is assigned to that edge chain.
        // Otherwise the edge is assigned to an arbitrary chain.
        val graphOptions = GraphOptions()  // Default options keep everything.
        val options = S2Builder.Options(IntLatLngSnapFunction(0))
        options.setSimplifyEdgeChains(true)
        val input = listOf(
            "0:1, 0:1.1", "0:0, 0:1, 0:2",  // Degenerate edge defined before chain
            "0:0, 0:0.9, 0:1, 0:1.1, 0:2",  // Degenerate edge defined in chain
            "0:2, 0:1, 0:0.9, 0:0",         // Defined in chain, chain reversed
            "0:2, 0:1, 0:0", "0:1.1, 0:1", "0:1, 0:1.1",  // Defined after chain
        )
        var expected: EdgeInputEdgeIds = listOf(
            Pair("0:0, 0:2", listOf(0, 1, 2)), Pair("0:0, 0:2", listOf(3, 4, 5, 6)),
            Pair("0:2, 0:0", listOf(7, 8, 9)), Pair("0:2, 0:0", listOf(10, 11, 12, 13))
        )
        testInputEdgeIds(input, expected, graphOptions, options)

        // Now try the same test with undirected edges.  This results in four more
        // simplified edges that are not labelled with any input edge ids.
        expected = expected + listOf(
            Pair("0:0, 0:2", emptyList()),
            Pair("0:0, 0:2", emptyList()),
            Pair("0:2, 0:0", emptyList()),
            Pair("0:2, 0:0", emptyList())
        )

        graphOptions.edgeType = EdgeType.UNDIRECTED
        testInputEdgeIds(input, expected, graphOptions, options)
    }

    @Test
    fun simplifyDegenerateEdgeMergingMultipleLayers() {
        // Check that degenerate edges are assigned to an edge in the correct layer
        // when multiple edge chains in different layers are simplified in the same
        // way (i.e., yielding a set of identical or reversed edges in different
        // layers).
        val graphOptions = GraphOptions()  // Default options keep everything.
        val options = S2Builder.Options(IntLatLngSnapFunction(0))
        options.setSimplifyEdgeChains(true)

        // Note below that the edge chains in different layers have different vertex
        // locations, different number of interior vertices, different degenerate
        // edges, etc, and yet they can all be simplified together.
        val input: List<List<String>> = listOf(
            listOf(
                "0.1:5, 0:5.2", "0.1:0, 0:9.9",   // Defined before chain
                "0:10.1, 0:0.1", "0:3.1, 0:2.9",  // Defined after chain
            ),
            listOf(
                "0.1:3, 0:3.2", "-0.1:0, 0:4.1, 0:9.9",  // Defined before chain
                "0.1:9.9, 0:7, 0.1:6.9, 0.1:0.2",        // Defined inside chain
            ),
            listOf(
                "0.2:0.3, 0.1:6, 0:5.9, 0.1:10.2",       // Defined inside chain
                "0.1:0.1, 0:9.8", "0.1:2, 0:2.1",        // Defined after chain
            )
        )
        val expected: List<EdgeInputEdgeIds> = listOf(
            listOf(Pair("0:0, 0:10", listOf(0, 1)), Pair("0:10, 0:0", listOf(2, 3))),
            listOf(Pair("0:0, 0:10", listOf(4, 5, 6)), Pair("0:10, 0:0", listOf(7, 8, 9))),
            listOf(Pair("0:0, 0:10", listOf(10, 11, 12)), Pair("0:0, 0:10", listOf(13, 14)))
        )
        val builder = S2Builder(options)
        for (i in input.indices) {
            builder.startLayer(InputEdgeIdCheckingLayer(expected[i], graphOptions))
            for (input_str in input[i]) {
                builder.addPolyline(makePolyline(input_str))
            }
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
    }

    @Test
    fun highPrecisionPredicates() {
        // To produce correct output in this example, the algorithm needs fall back
        // to high precision predicates when the output of the normal predicates is
        // uncertain.
        val vertices = listOf(
            S2Point(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
            S2Point(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),
            S2Point(-0.10531192039116592, -0.80522217309701472, 0.58354661457028933),
        )
        val input = S2Polyline(vertices)
        val snapRadius = S2EdgeCrossings.kIntersectionMergeRadius
        val options = S2Builder.Options((IdentitySnapFunction(snapRadius)))
        options.idempotent = false
        val builder = S2Builder(options)
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output))
        builder.forceVertex(S2Point(-0.10531192039134191, -0.80522217309705857, 0.58354661457019719))
        builder.addPolyline(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
    }

    // Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
    private fun choosePoint(): S2Point {
        val x = S2Random.randomPoint()
        for (i in 0..2) {
            if (S2Random.oneIn(3)) {
                x[i] *= 1e-50.pow(randomDouble())
            }
        }
        return x.normalize()
    }

    @Test
    fun highPrecisionStressTest() {
        // This test constructs many small, random inputs such that the output is
        // likely to be inconsistent unless high-precision predicates are used.

        val snapRadius = S2EdgeCrossings.kIntersectionMergeRadius
        // Some S2Builder calculations use an upper bound that takes into account
        // S1ChordAngle errors.  We sometimes try perturbing points by very close to
        // that distance in an attempt to expose errors.
        val ca = S1ChordAngle(snapRadius)
        val snapRadiusWithError =
            ca.plusError(ca.getS1AngleConstructorMaxError() + S2EdgeDistances.getUpdateMinDistanceMaxError(ca))
                .toAngle()

        var nonDegenerate = 0
        val kIters = 8000 * iteration_multiplier
        repeat(kIters) { iter ->
            // TODO(ericv): This test fails with a random seed of 96.  Change this
            // back to "iter + 1" once all the exact predicates are implemented.
            S2Random.reset(iter + 1) // Easier to reproduce a specific case.

            // We construct a nearly degenerate triangle where one of the edges is
            // sometimes very short.  Then we add a forced vertex somewhere near the
            // shortest edge.  Then after snapping, we check that (1) the edges still
            // form a loop, and (2) if the loop is non-degenerate, then it has the
            // same orientation as the original triangle.
            //
            // v1 is located randomly.  (v0,v1) is the longest of the three edges.
            // v2 is located along (v0,v1) but is perturbed by up to 2 * snap_radius.
            val v1 = choosePoint()
            val v0Dir = choosePoint()
            val d0 = 1e-16.pow(randomDouble())
            var v0 = S2EdgeDistances.interpolateAtDistance(S1Angle.radians(d0), v1, v0Dir)
            val d2 = 0.5 * d0 * 1e-16.pow(randomDouble().pow(2.0))
            var v2 = S2EdgeDistances.interpolateAtDistance(S1Angle.radians(d2), v1, v0Dir)
            v2 = S2Random.samplePoint(S2Cap.fromCenterAngle(v2, 2 * snapRadius))
            // Vary the edge directions by randomly swapping v0 and v2.
            if (S2Random.oneIn(2)) {
                val tmp = v0
                v0 = v2
                v2 = tmp
            }

            // The forced vertex (v3) is either located near the (v1, v2) edge.
            // We perturb it either in a random direction from v1 or v2, or
            // perpendicular to (v1, v2) starting from an interior edge point.
            var d3 = if (S2Random.oneIn(2)) snapRadius else snapRadiusWithError
            if (S2Random.oneIn(3)) d3 = 1.5 * randomDouble() * d3
            var v3: S2Point
            if (S2Random.oneIn(5)) {
                v3 = if (S2Random.oneIn(2)) v1 else v2
                v3 = S2EdgeDistances.interpolateAtDistance(d3, v3, choosePoint())
            } else {
                v3 = S2EdgeDistances.interpolate(1e-16.pow(randomDouble()), v1, v2)
                v3 = S2EdgeDistances.interpolateAtDistance(d3, v3, v1.crossProd(v2).normalize())
            }
            val options = S2Builder.Options((IdentitySnapFunction(snapRadius)))
            options.idempotent = false
            val builder = S2Builder(options)
            val output = S2Polygon()
            output.debugOverride = S2Debug.DISABLE
            builder.startLayer(S2PolygonLayer(output))
            builder.forceVertex(v3)
            builder.addEdge(v0, v1)
            builder.addEdge(v1, v2)
            builder.addEdge(v2, v0)
            val error = S2Error()
            if (!builder.build(error)) {
                logger.error { "d0=$d0, d2=$d2, d3=$d3" }
            }
            if (error.isOk() && !output.isEmpty()) {
                assertThat(output.numLoops()).isEqualTo(1)
                if (output.numLoops() == 1) {
                    assertThat(output.isValid()).isTrue()
                    assertThat(output.loop(0).isNormalized())
                        .withFailMessage("d0=$d0, d2=$d2, d3=$d3")
                        .isEqualTo(S2Predicates.sign(v0, v1, v2) > 0)
                    ++nonDegenerate
                }
            }
        }
        logger.debug { "$nonDegenerate non-degenerate out of $kIters" }
        assertThat(nonDegenerate).isGreaterThanOrEqualTo(kIters / 10)
    }

    @Test
    fun selfIntersectionStressTest() {
        val kIters = 50 * iteration_multiplier
        repeat(kIters) { iter ->
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val timer = Stopwatch.createStarted()

            // The minimum radius is about 36cm on the Earth's surface.  The
            // performance is reduced for radii much smaller than this because
            // S2ShapeIndex only indexes regions down to about 1cm across.
            val cap = S2Random.randomCap(1e-14, 1e-2)

            val options = S2Builder.Options()
            options.splitCrossingEdges = true
            if (S2Random.oneIn(2)) {
                val radius = cap.radius()
                val minExp = IntLatLngSnapFunction.exponentForMaxSnapRadius(radius)
                val exponent = min(IntLatLngSnapFunction.kMaxExponent, minExp + S2Random.randomInt(5))
                options.snapFunction = IntLatLngSnapFunction(exponent)
            }
            val builder = S2Builder(options)

            // Note that the number of intersections (and the running time) is
            // quadratic in the number of vertices.  With 200 input vertices, the
            // output consists of about 2300 loops and 9000 vertices.
            val output = S2Polygon()
            builder.startLayer(S2PolygonLayer(output, options = S2PolygonLayer.Options(EdgeType.UNDIRECTED)))
            val verticesCount = 50 // google.DEBUG_MODE ? 50 : 200
            val vertices = mutableListOf<S2Point>()
            repeat(verticesCount) {
                vertices.add(S2Random.samplePoint(cap))
            }
            vertices.add(vertices.first().clone())
            val input = S2Polyline(vertices)
            builder.addPolyline(input)
            val error = S2Error()
            assertThat(builder.build(error)).isTrue()
            assertThat(output.findValidationError(error)).isFalse()
            if (iter == -1) {
                logger.debug { "S2Polyline: " + S2TextParser.toString(input) }
                logger.debug { "S2Polygon: " + S2TextParser.toString(output) }
            }
            if (iter < 50) {
                logger.debug(
                    String.format(
                        "iter=%4d: ms=%4d, radius=%8.3g, loops=%d, vertices=%d",
                        iter, timer.elapsed(TimeUnit.MILLISECONDS),
                        cap.radius().radians, output.numLoops(),
                        output.numVertices()
                    )
                )
            }
        }
    }

    @Test
    fun fractalStressTest() {
        val kIters = 100 * iteration_multiplier // (google.DEBUG_MODE ? 100 : 1000)
        repeat(kIters) { iter ->
            logger.debug { "fractalStressTest | iteration $iter" }
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val fractal = Fractal()
            fractal.setLevelForApproxMaxEdges(800) // google.DEBUG_MODE ? 800 : 12800
            fractal.setLevelForApproxMinEdges(12)
            fractal.setDimension(1.5 + 0.5 * randomDouble())
            val input = S2Polygon(fractal.makeLoop(Matrix3x3Double.fromCols(randomFrame()), S1Angle.degrees(20)))

            logger.debug { "Input = ${S2TextParser.toString(input)}" }

            val options = S2Builder.Options()
            when {
                S2Random.oneIn(3) -> {
                    val exponent = S2Random.randomInt(11)
                    options.snapFunction = IntLatLngSnapFunction(exponent)
                }
                S2Random.oneIn(2) -> {
                    val level = S2Random.randomInt(20)
                    options.snapFunction = S2CellIdSnapFunction(level)
                }
                else -> {
                    options.snapFunction = IdentitySnapFunction(S1Angle.degrees(10 * 1e-4.pow(randomDouble())))

                }
            }

            logger.debug { "Builder options = $options" }
            val builder = S2Builder(options)
            val output = S2Polygon()
            builder.startLayer(S2PolygonLayer(output))
            builder.addPolygon(input)
            val error = S2Error()
            assertThat(builder.build(error)).isTrue()
            assertThat(output.findValidationError(error)).isFalse()
            if (iter == -1) {
                logger.debug { "S2Polygon: " + S2TextParser.toString(input) }
                logger.debug { "S2Polygon: " + S2TextParser.toString(output) }
            }
            if (iter < 50) {
                logger.debug(
                    String.format(
                        "iter=%4d: in_vertices=%d, out_vertices=%d",
                        iter,
                        input.numVertices(),
                        output.numVertices()
                    )
                )
            }
        }
    }

    private fun testSnappingWithForcedVertices(
        input_str: String,
        snap_radius: S1Angle,
        vertices_str: String,
        expected_str: String
    ) {
        val builder = S2Builder(S2Builder.Options(IdentitySnapFunction(snap_radius)))
        val vertices = S2TextParser.parsePoints(vertices_str)
        for (vertex in vertices) {
            builder.forceVertex(vertex)
        }
        val output = S2Polyline()
        builder.startLayer(S2PolylineLayer(output))
        builder.addPolyline(makePolyline(input_str))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()

        logger.debug { "testSnappingWithForcedVertices $input_str -> ${S2TextParser.toString(output)}" }
        logger.debug { output.vertices() }
        logger.debug { output.vertices().map { S2LatLng.fromPoint(it) }.map { ll -> "${ll.latRadians}:${ll.lngRadians}" } }

        assertThat(S2TextParser.toString(output)).isEqualTo(expected_str)
    }

    @Test
    fun adjacentCoverageIntervalsSpanMoreThan90Degrees() {
        // The test for whether one Voronoi site excludes another along a given
        // input edge boils down to a test of whether two angle intervals "a" and
        // "b" overlap.  Let "ra" and "rb" be the semi-widths of the two intervals,
        // and let "d" be the angle between their centers.  Then "a" contains "b" if
        // (rb + d <= ra), and "b" contains "a" if (rb - d >= ra).  However the
        // actual code uses the sines of the angles, e.g.  sin(rb + d) <= sin(ra).
        // This works fine most of the time, but the first condition (rb + d <= ra)
        // also needs to check that rb + d < 90 degrees.  This test verifies that
        // case.

        // The following 3 tests have d < 90, d = 90, and d > 90 degrees, but in all
        // 3 cases rb + d > 90 degrees.
        testSnappingWithForcedVertices(
            "0:0, 0:80", S1Angle.degrees(60),
            "0:0, 0:70", "0:0, 0:70"
        )
        testSnappingWithForcedVertices(
            "0:0, 0:80", S1Angle.degrees(60),
            "0:0, 0:90", "0:0, 0:90"
        )
        testSnappingWithForcedVertices(
            "0:0, 0:80", S1Angle.degrees(60),
            "0:0, 0:110", "0:0, 0:110"
        )

        // This test has d = 180 degrees, i.e. the two sites project to points that
        // are 180 degrees apart along the input edge.  The snapped edge doesn't
        // stay within max_edge_deviation() of the input edge, so an extra site is
        // added and it is snapped again (yielding two edges).  The case we are
        // testing here is the first call to SnapEdge() before adding the site.
        testSnappingWithForcedVertices(
            "0:10, 0:170", S1Angle.degrees(50),
            "47:0, 49:180", "47:0, 0:90, 49.00000000000001:180"
        )

        // This test has d = 220 degrees, i.e. when the input edge is snapped it
        // goes the "wrong way" around the sphere.  Again, the snapped edge is too
        // far from the input edge so an extra site is added and it is resnapped.
        testSnappingWithForcedVertices(
            "0:10, 0:170", S1Angle.degrees(70),
            "0:-20, 0:-160", "0:-20, 0:90, 0:-160"
        )

        // Without using forced vertices, the maximum angle between the coverage
        // interval centers is d = 300 degrees.  This would use an edge 180 degrees
        // long, and then place two sites 60 degrees past either endpoint.  With
        // forced vertices we can increase the snap radius to 70 degrees and get an
        // angle of up to d = 320 degrees, but the sites are only 40 degrees apart
        // (which is why it requires forced vertices).  The test below is an
        // approximation of this situation with d = 319.6 degrees.
        testSnappingWithForcedVertices(
            "0:0.1, 0:179.9", S1Angle.degrees(70),
            "0:-69.8, 0:-110.2",
            "0:-69.8, 0:90, 0:-110.2"
        )
    }

    @Test
    fun oldS2PolygonBuilderBug() {
        // This is a polygon that caused the obsolete S2PolygonBuilder class to
        // generate an invalid output polygon (duplicate edges).
        val input = makePolygon(
            "32.2983095:72.3416582, 32.2986281:72.3423059, " +
                    "32.2985238:72.3423743, 32.2987176:72.3427807, " +
                    "32.2988174:72.3427056, 32.2991269:72.3433480, " +
                    "32.2991881:72.3433077, 32.2990668:72.3430462, " +
                    "32.2991745:72.3429778, 32.2995078:72.3436725, " +
                    "32.2996075:72.3436269, 32.2985465:72.3413832, " +
                    "32.2984558:72.3414530, 32.2988015:72.3421839, " +
                    "32.2991552:72.3429416, 32.2990498:72.3430073, " +
                    "32.2983764:72.3416059"
        )
        assertThat(input.isValid()).isTrue()

        val snapRadius = S1Angle.radians(S2Earth.metersToRadians(20 / 0.866))
        val builder = S2Builder(S2Builder.Options(IdentitySnapFunction(snapRadius)))
        val output = S2Polygon()
        builder.startLayer(S2PolygonLayer(output))
        builder.addPolygon(input)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(output.isValid()).isTrue()
        val expected = makePolygon(
            "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; " +
                    "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, " +
                    "32.2985238:72.3423743, 32.2987176:72.3427807"
        )
        expectPolygonsEqual(expected, output)
    }


}
