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

import dilivia.s2.S1Angle
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.layers.IndexedLaxPolygonLayer
import dilivia.s2.builder.layers.IndexedS2PointVectorLayer
import dilivia.s2.builder.layers.IndexedS2PolylineVectorLayer
import dilivia.s2.builder.layers.LaxPolygonLayer
import dilivia.s2.builder.layers.Layer
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.index.shape.S2BooleanOperation.OpType
import dilivia.s2.index.shape.S2BooleanOperation.PolygonModel
import dilivia.s2.index.shape.S2BooleanOperation.PolylineModel
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2LaxPolygonShape
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

typealias EdgeVector = List<Edge>

class IndexMatchingLayer(val index: S2ShapeIndex, val dimension: Int) : Layer() {

    override fun graphOptions(): GraphOptions = GraphOptions(
        edgeType = EdgeType.DIRECTED,
        degenerateEdges = DegenerateEdges.KEEP,
        duplicateEdges = DuplicateEdges.KEEP,
        siblingPairs = SiblingPairs.KEEP
    )

    override fun build(g: Graph, error: S2Error) {
        val actual = mutableListOf<Edge>()
        val expected = mutableListOf<Edge>()
        for (e in 0 until g.numEdges) {
            val edge = g.edge(e)
            actual.add(Edge(g.vertex(edge.first), g.vertex(edge.second)))
        }
        for (shape in index) {
            if (shape == null || shape.dimension != dimension) {
                continue
            }
            var e = shape.numEdges
            while (--e >= 0) {
                expected.add(shape.edge(e))
            }
        }
        actual.sort()
        expected.sort()

        // The edges are a multiset, so we can't use std::set_difference.
        val missing = mutableListOf<Edge>()
        val extra = mutableListOf<Edge>()
        var ai = 0
        var ei = 0
        while (ai != actual.size || ei != expected.size) {
            if (ei == expected.size || (ai != actual.size && actual[ai] < expected[ei])) {
                extra.add(actual[ai++])
            } else if (ai == actual.size || expected[ei] < actual[ai]) {
                missing.add(expected[ei++])
            } else {
                ++ai
                ++ei
            }
        }
        if (missing.isNotEmpty() || extra.isNotEmpty()) {
            // There may be errors in more than one dimension, so we append to the
            // existing error text.
            error.init(
                S2Error.USER_DEFINED_START,
                "%sDimension %d: Missing edges: %s Extra edges: %s\n".format(
                    error.text, dimension, toString(missing), toString(extra)
                )
            )
        }
    }

    companion object {

        fun toString(edges: EdgeVector): String {
            var msg = ""
            for (edge in edges) {
                val vertices = listOf(edge.v0, edge.v1)
                msg += S2TextParser.toString(vertices)
                msg += "; "
            }
            return msg
        }

    }

}

class S2BooleanOperationUnitTest {

    companion object {

        // The intersections in the "expected" data below were computed in lat-lng
        // space (i.e., the rectangular projection), while the actual intersections
        // are computed using geodesics.  We can compensate for this by rounding the
        // intersection points to a fixed precision in degrees (e.g., 2 decimals).
        fun roundToE(exp: Int): S2BooleanOperation.Options {
            val options = S2BooleanOperation.Options()
            options.snapFunction = IntLatLngSnapFunction(exp)
            return options
        }

        // The polygon used in the polyline/polygon vertex tests below.
        fun kVertexTestPolygonStr(): String = "0:0, 0:1, 0:2, 0:3, 0:4, 0:5, 5:5, 5:4, 5:3, 5:2, 5:1, 5:0"

        fun expectResult(
            op_type: OpType,
            options: S2BooleanOperation.Options,
            a_str: String,
            b_str: String,
            expected_str: String
        ) {
            val a = S2TextParser.makeIndex(a_str)
            val b = S2TextParser.makeIndex(b_str)
            val expected = S2TextParser.makeIndex(expected_str)
            val layers = mutableListOf<Layer>()
            for (dim in 0..2) {
                layers.add(IndexMatchingLayer(expected, dim))
            }
            val op = S2BooleanOperation(op_type, layers, options)
            val error = S2Error()
            assertThat(op.build(a, b, error))
                .withFailMessage("$op_type <failed:\nExpected result: $expected_str\n$error")
                .isTrue()

            // Now try the same thing with boolean output.
            assertThat(S2BooleanOperation.isEmpty(op_type, a, b, options)).isEqualTo(expected.nextNewShapeId() == 0)
        }


        // Performs the given operation and compares the result to "expected_str".  All
// arguments are in s2textformat::MakeLaxPolygon() format.
        fun expectPolygon(op_type: OpType, a_str: String, b_str: String, expected_str: String) {
            val a = S2TextParser.makeIndex("# # $a_str")
            val b = S2TextParser.makeIndex("# # $b_str")
            val polygonOptions = LaxPolygonLayer.Options()
            polygonOptions.degenerateBoundaries = LaxPolygonLayer.DegenerateBoundaries.DISCARD
            val output = S2LaxPolygonShape()
            val op = S2BooleanOperation(
                op_type,
                LaxPolygonLayer(output, options = polygonOptions),
                S2BooleanOperation.Options(IdentitySnapFunction(S1Angle.degrees(1.1)))
            )
            val error = S2Error()
            assertThat(op.build(a, b, error))
                .withFailMessage(error.toString())
                .isTrue()
            assertThat(S2TextParser.toString(output)).isEqualTo(expected_str)
        }
    }

// TODO(ericv): Clean up or remove these notes.
//
// Options to test:
//   polygon_model:                   OPEN, SEMI_OPEN, CLOSED
//   polyline_model:                  OPEN, SEMI_OPEN, CLOSED
//   polyline_loops_have_boundaries:  true, false
//   conservative:                    true, false
//
// Geometry combinations to test:
//
// Point/point:
//  - disjoint, coincident
// Point/polyline:
//  - Start vertex, end vertex, interior vertex, degenerate polyline
//  - With polyline_loops_have_boundaries: start/end vertex, degenerate polyline
// Point/polygon:
//  - Polygon interior, exterior, vertex
//  - Vertex of degenerate sibling pair shell, hole
//  - Vertex of degenerate single point shell, hole
// Polyline/polyline:
//  - Vertex intersection:
//    - Start, end, interior, degenerate, loop start/end, degenerate loop
//    - Test cases where vertex is not emitted because an incident edge is.
//  - Edge/edge: interior crossing, duplicate, reversed, degenerate
//  - Test that degenerate edges are ignored unless polyline has a single edge.
//    (For example, AA has one edge but AAA has no edges.)
// Polyline/polygon:
//  - Vertex intersection: polyline vertex cases already covered, but test
//    polygon normal vertex, sibling pair shell/hole, single vertex shell/hole
//    - Also test cases where vertex is not emitted because an edge is.
//  - Edge/edge: interior crossing, duplicate, reversed
//  - Edge/interior: polyline edge in polygon interior, exterior
// Polygon/polygon:
//  - Vertex intersection:
//    - normal vertex, sibling pair shell/hole, single vertex shell/hole
//    - Also test cases where vertex is not emitted because an edge is.
//    - Test that polygons take priority when there is a polygon vertex and
//      also isolated polyline vertices.  (There should not be any points.)
//  - Edge/edge: interior crossing, duplicate, reversed
//  - Interior/interior: polygons in interior/exterior of other polygons

    @Test
    fun degeneratePolylines() {
        // Verify that degenerate polylines are preserved under all boundary models.
        val options = S2BooleanOperation.Options()
        val a = "# 0:0, 0:0 #"
        val b = "# #"
        options.polylineModel = PolylineModel.OPEN
        expectResult(OpType.UNION, options, a, b, a)
        options.polylineModel = PolylineModel.SEMI_OPEN
        expectResult(OpType.UNION, options, a, b, a)
        options.polylineModel = PolylineModel.CLOSED
        expectResult(OpType.UNION, options, a, b, a)
    }

    @Test
    fun degeneratePolygons() {
        // Verify that degenerate polygon features (single-vertex and sibling pair
        // shells and holes) are preserved under all boundary models.
        val options = S2BooleanOperation.Options()
        val a = "# # 0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 3:3; 6:6; 7:7, 8:8"
        val b = "# #"
        options.polygonModel = PolygonModel.OPEN
        expectResult(OpType.UNION, options, a, b, a)
        options.polygonModel = PolygonModel.SEMI_OPEN
        expectResult(OpType.UNION, options, a, b, a)
        options.polygonModel = PolygonModel.CLOSED
        expectResult(OpType.UNION, options, a, b, a)
    }

    @Test
    fun pointPoint() {
        val options = S2BooleanOperation.Options()
        val a = "0:0 | 1:0 # #"
        val b = "0:0 | 2:0 # #"
        // Note that these results have duplicates, which is correct.  Clients can
        // eliminated the duplicates with the appropriate GraphOptions.
        expectResult(OpType.UNION, options, a, b, "0:0 | 0:0 | 1:0 | 2:0 # #")
        expectResult(OpType.INTERSECTION, options, a, b, "0:0 | 0:0 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "1:0 # #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "1:0 | 2:0 # #")
    }

    @Test
    fun pointOpenPolyline() {
        // Tests operations between an open polyline and its vertices.
        //
        // The polyline "3:0, 3:0" consists of a single degenerate edge and contains
        // no points (since polyline_model() is OPEN).  Since S2BooleanOperation
        // preserves degeneracies, this means that the union includes *both* the
        // point 3:0 and the degenerate polyline 3:0, since they do not intersect.
        //
        // This test uses Options::polyline_loops_have_boundaries() == true, which
        // means that the loop "4:0, 5:0, 4:0" does not contain the vertex "4:0".
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.OPEN
        val a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #"
        val b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
        expectResult(OpType.UNION, options, a, b, "0:0 | 2:0 | 3:0 | 4:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #")
        expectResult(OpType.INTERSECTION, options, a, b, "1:0 | 5:0 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "0:0 | 2:0 | 3:0 | 4:0 # #")
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE,
            options,
            a,
            b,
            "0:0 | 2:0 | 3:0 | 4:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
        )
    }

    @Test
    fun pointOpenPolylineLoopBoundariesFalse() {
        // With Options::polyline_loops_have_boundaries() == false, the loop
        // "4:0, 5:0, 4:0" has two vertices, both of which are contained.
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.OPEN
        options.polylineLoopsHaveBoundaries = false
        val a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #"
        val b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
        expectResult(OpType.UNION, options, a, b, "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #")
        expectResult(OpType.INTERSECTION, options, a, b, "1:0 | 4:0 | 5:0 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "0:0 | 2:0 | 3:0 # #")
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE,
            options,
            a,
            b,
            "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
        )
    }

    @Test
    fun pointSemiOpenPolyline() {
        // Degenerate polylines are defined not contain any points under the
        // SEMI_OPEN model either, so again the point 3:0 and the degenerate
        // polyline "3:0, 3:0" do not intersect.
        //
        // The result does not depend on Options::polyline_loops_have_boundaries().
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.SEMI_OPEN
        for (bool_value in listOf(false, true)) {
            options.polylineLoopsHaveBoundaries = bool_value
            val a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #"
            val b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
            expectResult(OpType.UNION, options, a, b, "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #")
            expectResult(OpType.INTERSECTION, options, a, b, "0:0 | 1:0 | 4:0 | 5:0 # #")
            expectResult(OpType.DIFFERENCE, options, a, b, "2:0 | 3:0 # #")
            expectResult(
                OpType.SYMMETRIC_DIFFERENCE,
                options,
                a,
                b,
                "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
            )
        }
    }

    @Test
    fun pointClosedPolyline() {
        // Under the CLOSED model, the degenerate polyline 3:0 does contain its
        // vertex.  Since polylines take precedence over points, the union of the
        // point 3:0 and the polyline 3:0 is the polyline only.  Similarly, since
        // subtracting a point from a polyline has no effect, the symmetric
        // difference includes only the polyline objects.
        //
        // The result does not depend on Options::polyline_loops_have_boundaries().
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.CLOSED
        for (bool_value in listOf(false, true)) {
            options.polylineLoopsHaveBoundaries = bool_value
            val a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #"
            val b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #"
            expectResult(OpType.UNION, options, a, b, "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #")
            expectResult(OpType.INTERSECTION, options, a, b, "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #")
            expectResult(OpType.DIFFERENCE, options, a, b, "# #")
            expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #")
        }
    }

    @Test
    fun pointPolygonInterior() {
        val options = S2BooleanOperation.Options()  // PolygonModel is irrelevant.
        // One interior point and one exterior point.
        val a = "1:1 | 4:4 # #"
        val b = "# # 0:0, 0:3, 3:0"
        expectResult(OpType.UNION, options, a, b, "4:4 # # 0:0, 0:3, 3:0")
        expectResult(OpType.INTERSECTION, options, a, b, "1:1 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "4:4 # #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "4:4 # # 0:0, 0:3, 3:0")
    }

    @Test
    fun pointOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN
        // See notes about the two vertices below.
        val a = "0:1 | 1:0 # #"
        val b = "# # 0:0, 0:1, 1:0"
        expectResult(OpType.UNION, options, a, b, "0:1 | 1:0 # # 0:0, 0:1, 1:0")
        expectResult(OpType.INTERSECTION, options, a, b, "# #")
        expectResult(OpType.DIFFERENCE, options, a, b, "0:1 | 1:0 # #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "0:1 | 1:0 # # 0:0, 0:1, 1:0")
    }

    @Test
    fun pointSemiOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.SEMI_OPEN
        // The two vertices are chosen such that the polygon contains one vertex but
        // not the other under PolygonModel::SEMI_OPEN.  (The same vertices are used
        // for all three PolygonModel options.)
        val polygon = S2TextParser.makePolygon("0:0, 0:1, 1:0")
        assertThat(polygon.contains(S2TextParser.makePoint("0:1"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("1:0"))).isFalse()
        val a = "0:1 | 1:0 # #"
        val b = "# # 0:0, 0:1, 1:0"
        expectResult(OpType.UNION, options, a, b, "1:0 # # 0:0, 0:1, 1:0")
        expectResult(OpType.INTERSECTION, options, a, b, "0:1 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "1:0 # #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "1:0 # # 0:0, 0:1, 1:0")
    }

    @Test
    fun pointClosedPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.CLOSED
// See notes about the two vertices above.
        val a = "0:1 | 1:0 # #"
        val b = "# # 0:0, 0:1, 1:0"
        expectResult(OpType.UNION, options, a, b, "# # 0:0, 0:1, 1:0")
        expectResult(OpType.INTERSECTION, options, a, b, "0:1 | 1:0 # #")
        expectResult(OpType.DIFFERENCE, options, a, b, "# #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "# # 0:0, 0:1, 1:0")
    }

    @Test
    fun PolylineVertexOpenPolylineVertex() {
        // Test first, last, and middle vertices of both polylines.  Also test
        // first/last and middle vertices of two polyline loops.
        //
        // Degenerate polylines are tested in PolylineEdgePolylineEdgeOverlap below.
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.OPEN
        val a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        val b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )

        // The output consists of the portion of each input polyline that intersects
        // the opposite region, so the intersection vertex is present twice.  This
        // allows reassembling the individual polylins that intersect, if desired.
        // (Otherwise duplicates can be removed using DuplicateEdges::MERGE.)
        expectResult(OpType.INTERSECTION, options, a, b, "# 0:1, 0:1 | 0:1, 0:1 #")

        // Note that all operations are defined such that subtracting a
        // lower-dimensional subset of an object has no effect.  In this case,
        // subtracting the middle vertex of a polyline has no effect.
        expectResult(OpType.DIFFERENCE, options, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #")
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )
    }


    @Test
    fun PolylineVertexOpenPolylineVertexLoopBoundariesFalse() {
        // With Options::polyline_loops_have_boundaries() == false, the 3 polyline
        // loops each have two vertices, both of which are contained.
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.OPEN
        options.polylineLoopsHaveBoundaries = false
        val a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        val b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )

        // Note that the polyline "0:3, 0:4, 0:3" only has two vertices, not three.
        // This means that 0:3 is emitted only once for that polyline, plus once for
        // the other polyline, for a total of twice.
        expectResult(
            OpType.INTERSECTION,
            options,
            a,
            b,
            "# 0:1, 0:1 | 0:1, 0:1 | 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #"
        )

        expectResult(OpType.DIFFERENCE, options, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #")
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )
    }


    @Test
    fun PolylineVertexSemiOpenPolylineVertex() {
        // The result does not depend on Options::polyline_loops_have_boundaries().
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.SEMI_OPEN
        for (bool_value in listOf(false, true)) {
            options.polylineLoopsHaveBoundaries = bool_value
            val a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
            val b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
            expectResult(
                OpType.UNION, options, a, b,
                "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
            )
            expectResult(
                OpType.INTERSECTION, options, a, b,
                "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 | 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #"
            )
            expectResult(
                OpType.DIFFERENCE, options, a, b,
                "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
            )
            expectResult(
                OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
            )
        }
    }


    @Test
    fun PolylineVertexClosedPolylineVertex() {
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.CLOSED
        val a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        val b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )

        // Since Options::polyline_loops_have_boundaries() == true, the polyline
        // "0:3, 0:4, 0:3" has three vertices.  Therefore 0:3 is emitted twice for
        // that polyline, plus once for the other polyline, for a total of thrice.
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 " +
                    "| 0:2, 0:2 | 0:2, 0:2 " +
                    "| 0:3, 0:3 | 0:3, 0:3 | 0:3, 0:3 " +
                    "| 0:4, 0:4 | 0:4, 0:4 | 0:4, 0:4 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                    "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )
    }


    @Test
    fun PolylineVertexClosedPolylineVertexLoopBoundariesFalse() {
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.CLOSED
        options.polylineLoopsHaveBoundaries = false
        val a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        val b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )

        // Since Options::polyline_loops_have_boundaries() == false, the polyline
        // "0:3, 0:4, 0:3" has two vertices.  Therefore 0:3 is emitted once for
        // that polyline, plus once for the other polyline, for a total of twice.
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 " +
                    "| 0:2, 0:2 | 0:2, 0:2 " +
                    "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #"
        )
    }


    @Test
    fun TestSemiOpenPolygonVerticesContained() {
        // Verify whether certain vertices of the test polygon are contained under
        // the semi-open boundary model (for use in the tests below).
        val polygon = S2TextParser.makePolygon(kVertexTestPolygonStr())
        assertThat(polygon.contains(S2TextParser.makePoint("0:1"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("0:2"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("0:3"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("0:4"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("5:1"))).isFalse()
        assertThat(polygon.contains(S2TextParser.makePoint("5:2"))).isFalse()
        assertThat(polygon.contains(S2TextParser.makePoint("5:3"))).isFalse()
        assertThat(polygon.contains(S2TextParser.makePoint("5:4"))).isFalse()
    }

    // Don't bother testing every PolylineModel with every PolygonModel for vertex
// intersection, since we have already tested the PolylineModels individually
// above.  It is sufficient to use PolylineModel::CLOSED with the various
// PolygonModel options.

    @Test
    fun PolylineVertexOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN

        // Define some constants to reduce code duplication.
        // Test all combinations of polylines that start or end on a polygon vertex,
        // where the polygon vertex is open or closed using semi-open boundaries,
        // and where the incident edge is inside or outside the polygon.
        val a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #")
        val b = "# # " + kVertexTestPolygonStr()

        val kDifferenceResult =
            "# 0:1, 0:1 | 0:2, 0:2 | -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #"
        expectResult(OpType.UNION, options, a, b, kDifferenceResult + kVertexTestPolygonStr())
        expectResult(OpType.INTERSECTION, options, a, b, "# 1:1, 0:1 | 0:2, 1:2 | 4:3, 5:3 | 5:4, 4:4 #")
        expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult)
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, kDifferenceResult + kVertexTestPolygonStr())
    }

    // Like the test above, except that every polygon vertex is also incident to a
// closed polyline vertex.  This tests that when an open vertex and a closed
// vertex coincide with each other, the result is considered closed.

    @Test
    fun PolylineVertexOpenPolygonClosedPolylineVertex() {
        val kTestGeometrySuffix =
            "-2:0, 0:1 | -2:1, 0:2 | -2:2, 0:3 | -2:3, 0:4 | 7:0, 5:1 | 7:1, 5:2 | 7:2, 5:3 | 7:3, 5:4 # " + kVertexTestPolygonStr()

        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN
        val a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #")
        val b = ("# $kTestGeometrySuffix")

        val kDifferencePrefix = "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2"
        expectResult(
            OpType.UNION, options, a, b,
            "$kDifferencePrefix | 0:1, 0:1 | 0:2, 0:2 | 5:3, 5:3 | 5:4, 5:4 | $kTestGeometrySuffix"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4" +
                    "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4" +
                    "| 0:1, 0:1 | 0:2, 0:2 | 0:3, 0:3 | 0:4, 0:4" +
                    "| 5:1, 5:1 | 5:2, 5:2 | 5:3, 5:3 | 5:4, 5:4 #"
        )
        expectResult(OpType.DIFFERENCE, options, a, b, "$kDifferencePrefix #")
        expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "$kDifferencePrefix | $kTestGeometrySuffix")
    }


    @Test
    fun PolylineVertexSemiOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.SEMI_OPEN
        // Test all combinations of polylines that start or end on a polygon vertex,
        // where the polygon vertex is open or closed using semi-open boundaries,
        // and where the incident edge is inside or outside the polygon.
        //
        // The vertices at latitude 0 used below are all closed while the vertices
        // at latitude 5 are all open (see TestSemiOpenPolygonVerticesContained).
        val a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #")
        val b = "# # " + kVertexTestPolygonStr()
        val kDifferenceResult =
            "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #"
        expectResult(
            OpType.UNION, options, a, b,
            kDifferenceResult + kVertexTestPolygonStr()
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4 | 4:3, 5:3 | 5:4, 4:4 #"
        )
        expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult)
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            kDifferenceResult + kVertexTestPolygonStr()
        )
    }


    @Test
    fun PolylineVertexClosedPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.CLOSED
        // Test all combinations of polylines that start or end on a polygon vertex,
        // where the polygon vertex is open or closed using semi-open boundaries,
        // and where the incident edge is inside or outside the polygon.
        val a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #")
        val b = "# # " + kVertexTestPolygonStr()
        val kDifferenceResult = "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 #"
        expectResult(
            OpType.UNION, options, a, b,
            kDifferenceResult + kVertexTestPolygonStr()
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4 | 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4 #"
        )
        expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult)
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            kDifferenceResult + kVertexTestPolygonStr()
        )
    }

    @Test
    fun PolylineEdgePolylineEdgeCrossing() {
        // Two polyline edges that cross at a point interior to both edges.
        val options = roundToE(1)
        val a = "# 0:0, 2:2 #"
        val b = "# 2:0, 0:2 #"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 1:1 | 1:1, 1:1 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 0:0, 2:2 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #"
        )
    }


    @Test
    fun PolylineEdgePolylineEdgeOverlap() {
        // The PolylineModel does not affect this calculation.  In particular the
        // intersection of a degenerate polyline edge with itself is non-empty, even
        // though the edge contains no points in the OPEN and SEMI_OPEN models.
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN
        // Test edges in the same and reverse directions, and degenerate edges.
        val a = "# 0:0, 1:0, 2:0, 2:5 | 3:0, 3:0 | 6:0, 5:0, 4:0 #"
        val b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0 #"
        // As usual, the expected output includes the relevant portions of *both*
        // input polylines.  Duplicates can be removed using GraphOptions.
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 1:0, 2:0, 2:5 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 | 6:0, 5:0, 4:0 | 4:0, 5:0 #"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 0:0, 1:0, 2:0 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 | 5:0, 4:0 | 4:0, 5:0 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 2:0, 2:5 | 6:0, 5:0 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 2:0, 2:5 | 6:0, 5:0 #"
        )
    }


    @Test
    fun PolylineEdgeOpenPolygonEdgeOverlap() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN
        // A polygon and two polyline edges that coincide with the polygon boundary,
        // one in the same direction and one in the reverse direction.
        val a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # "
        val b = "# # 1:1, 1:3, 3:3, 3:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 1:1, 1:3, 3:3 | 3:3, 1:3 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1"
        )
    }


    @Test
    fun PolylineEdgeSemiOpenPolygonEdgeOverlap() {
        val polygon = S2TextParser.makePolygon("1:1, 1:3, 3:3, 3:1")
        assertThat(polygon.contains(S2TextParser.makePoint("1:1"))).isFalse()
        assertThat(polygon.contains(S2TextParser.makePoint("1:3"))).isTrue()
        assertThat(polygon.contains(S2TextParser.makePoint("3:3"))).isFalse()
        assertThat(polygon.contains(S2TextParser.makePoint("3:1"))).isFalse()
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.SEMI_OPEN
        val a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # "
        val b = "# # 1:1, 1:3, 3:3, 3:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:3, 1:3 | 1:1, 1:3, 3:3 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1"
        )
    }


    @Test
    fun PolylineEdgeClosedPolygonEdgeOverlap() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.CLOSED
        val a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # "
        val b = "# # 1:1, 1:3, 3:3, 3:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 1:1, 1:3, 3:3, 3:1"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 1:3, 3:3 | 3:3, 1:3 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 1:1, 1:3, 3:3, 3:1"
        )
    }

    @Test
    fun PolygonVertexMatching() {
        // This test shows that CrossingProcessor::ProcessEdgeCrossings() must set
        // a0_matches_polygon and a1_matches_polygon correctly even when (a0, a1)
        // itself is a polygon edge (or its sibling).  (It requires degenerate
        // polygon geometry to demonstrate this.)
        val options = S2BooleanOperation.Options()
        options.polylineModel = PolylineModel.CLOSED
        options.polygonModel = PolygonModel.CLOSED
        val a = "# 0:0, 1:1 # "
        val b = "# # 0:0, 1:1"
        expectResult(OpType.UNION, options, a, b, "# # 0:0, 1:1")
    }


    @Test
    fun PolylineEdgePolygonInterior() {
        val options = S2BooleanOperation.Options()  // PolygonModel is irrelevant.
        // One normal and one degenerate polyline edge in the polygon interior, and
        // similarly for the polygon exterior.
        val a = "# 1:1, 2:2 | 3:3, 3:3 | 6:6, 7:7 | 8:8, 8:8 # "
        val b = "# # 0:0, 0:5, 5:5, 5:0"
        expectResult(
            OpType.UNION, options, a, b,
            "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 2:2 | 3:3, 3:3 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 6:6, 7:7 | 8:8, 8:8 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0"
        )
    }


    @Test
    fun PolygonVertexOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.OPEN
        val a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        val b = "# # 0:0, 5:3, 5:2"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
    }


    @Test
    fun PolygonVertexSemiOpenPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.SEMI_OPEN
        val a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        val b = "# # 0:0, 5:3, 5:2"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
    }


    @Test
    fun PolygonVertexClosedPolygonVertex() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.CLOSED
        val a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        val b = "# # 0:0, 5:3, 5:2"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 0:0"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5"
        )
        expectResult(
            OpType.DIFFERENCE, options, b, a,
            "# # 0:0, 5:3, 5:2"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2"
        )
    }


    @Test
    fun PolygonEdgePolygonEdgeCrossing() {
        // Two polygons whose edges cross at points interior to both edges.
        val options = roundToE(2)
        val a = "# # 0:0, 0:2, 2:2, 2:0"
        val b = "# # 1:1, 1:3, 3:3, 3:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:3, 3:3, 3:1, 2:1, 2:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 1:1, 1:2, 2:2, 2:1"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0; 1:2, 1:3, 3:3, 3:1, 2:1, 2:2"
        )
    }


    @Test
    fun PolygonEdgeOpenPolygonEdgeOverlap() {
        val options = S2BooleanOperation.Options()
        // One shape is a rectangle, the other consists of one triangle inside the
        // rectangle and one triangle outside the rectangle, where each triangle
        // shares one edge with the rectangle.  This implies that the edges are in
        // the same direction in one case and opposite directions in the other case.
        options.polygonModel = PolygonModel.OPEN
        val a = "# # 0:0, 0:4, 2:4, 2:0"
        val b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0; 0:4, 1:5, 2:4"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 0:0, 1:1, 2:0"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4"
        )
    }


    @Test
    fun PolygonEdgeSemiOpenPolygonEdgeOverlap() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.SEMI_OPEN
        val a = "# # 0:0, 0:4, 2:4, 2:0"
        val b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:4, 1:5, 2:4, 2:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 0:0, 1:1, 2:0"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1"
        )
        // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
        // normalized, i.e. the output could contain siblings pairs (which can be
        // discarded using S2Builder::GraphOptions).
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4"
        )
    }


    @Test
    fun PolygonEdgeClosedPolygonEdgeOverlap() {
        val options = S2BooleanOperation.Options()
        options.polygonModel = PolygonModel.CLOSED
        val a = "# # 0:0, 0:4, 2:4, 2:0"
        val b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:4, 1:5, 2:4, 2:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 0:0, 1:1, 2:0; 0:4, 2:4"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1"
        )
        // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
        // normalized, i.e. the output could contain siblings pairs.
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4"
        )
    }


    @Test
    fun PolygonPolygonInterior() {
        val options = S2BooleanOperation.Options()  // PolygonModel is irrelevant.
        // One loop in the interior of another polygon and one loop in the exterior.
        val a = "# # 0:0, 0:4, 4:4, 4:0"
        val b = "# # 1:1, 1:2, 2:2, 2:1; 5:5, 5:6, 6:6, 6:5"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:4, 4:4, 4:0; 5:5, 5:6, 6:6, 6:5"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 1:1, 1:2, 2:2, 2:1"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1; 5:5, 5:6, 6:6, 6:5"
        )
    }

    @Test
    fun PolygonEdgesDegenerateAfterSnapping() {
        val options = roundToE(0)
        val a = "# # 0:-1, 0:1, 0.1:1, 0.1:-1"
        val b = "# # -1:0.1, 1:0.1, 1:0, -1:0"
        // When snapping causes an output edge to become degenerate, it is still
        // emitted (since otherwise loops that contract to a single point would be
        // lost).  If the output layer doesn't want such edges, they can be removed
        // via DegenerateEdges::DISCARD or DISCARD_EXCESS.
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | -1:0, -1:0, 0:0, 1:0, 1:0, 0:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 0:0, 0:0, 0:0, 0:0"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | 0:0, 0:0"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | -1:0, -1:0, 0:0, 1:0, 1:0, 0:0 | 0:0, 0:0, 0:0, 0:0"
        )
    }

///////////////////////////////////////////////////////////////////////////
// The remaining tests are intended to cover combinations of features or
// interesting special cases.


    @Test
    fun ThreeOverlappingBars() {
        // Two vertical bars and a horizontal bar that overlaps both of the other
        // bars and connects them.

        // Round intersection points to E2 precision because the expected results
        // were computed in lat/lng space rather than using geodesics.
        val options = roundToE(2)
        val a = "# # 0:0, 0:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3"
        val b = "# # 1:1, 1:4, 2:4, 2:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 1:1, 1:2, 2:2, 2:1; 1:3, 1:4, 2:4, 2:3"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3; 1:2, 1:3, 2:3, 2:2"
        )
    }


    @Test
    fun FourOverlappingBars() {
        // Two vertical bars and two horizontal bars.

        // Round intersection points to E2 precision because the expected results
        // were computed in lat/lng space rather than using geodesics.
        val options = roundToE(2)
        val a = "# # 1:88, 1:93, 2:93, 2:88; -1:88, -1:93, 0:93, 0:88"
        val b = "# # -2:89, -2:90, 3:90, 3:89; -2:91, -2:92, 3:92, 3:91"
        expectResult(
            OpType.UNION, options, a, b,
            "# # -1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, " +
                    "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, " +
                    "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88; " +
                    "0:90, 1:90, 1:91, 0:91" /*CW*/
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # 1:89, 1:90, 2:90, 2:89; 1:91, 1:92, 2:92, 2:91; " +
                    "-1:89, -1:90, 0:90, 0:89; -1:91, -1:92, 0:92, 0:91"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # 1:88, 1:89, 2:89, 2:88; 1:90, 1:91, 2:91, 2:90; " +
                    "1:92, 1:93, 2:93, 2:92; -1:88, -1:89, 0:89, 0:88; " +
                    "-1:90, -1:91, 0:91, 0:90; -1:92, -1:93, 0:93, 0:92"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # 1:88, 1:89, 2:89, 2:88; -1:88, -1:89, 0:89, 0:88; " +
                    "1:90, 1:91, 2:91, 2:90; -1:90, -1:91, 0:91, 0:90; " +
                    "1:92, 1:93, 2:93, 2:92; -1:92, -1:93, 0:93, 0:92; " +
                    "-2:89, -2:90, -1:90, -1:89; -2:91, -2:92, -1:92, -1:91; " +
                    "0:89, 0:90, 1:90, 1:89; 0:91, 0:92, 1:92, 1:91; " +
                    "2:89, 2:90, 3:90, 3:89; 2:91, 2:92, 3:92, 3:91"
        )
    }


    @Test
    fun OverlappingDoughnuts() {
        // Two overlapping square doughnuts whose holes do not overlap.
        // This means that the union polygon has only two holes rather than three.

        // Round intersection points to E2 precision because the expected results
        // were computed in lat/lng space rather than using geodesics.
        val options = roundToE(1)
        val a = "# # -1:-93, -1:-89, 3:-89, 3:-93; " +
                "0:-92, 2:-92, 2:-90, 0:-90" /*CW*/
        val b = "# # -3:-91, -3:-87, 1:-87, 1:-91; " +
                "-2:-90, 0:-90, 0:-88, -2:-88" /*CW*/
        expectResult(
            OpType.UNION, options, a, b,
            "# # -1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93; " +
                    "0:-92, 2:-92, 2:-90, 1:-90, 1:-91, 0:-91; " /*CW */ +
                    "-2:-90, -1:-90, -1:-89, 0:-89, 0:-88, -2:-88" /* CW */
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# # -1:-91, -1:-90, 0:-90, 0:-91; " +
                    "0:-90, 0:-89, 1:-89, 1:-90"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, " +
                    "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; " +
                    "-1:-90, -1:-89, 0:-89, 0:-90"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, " +
                    "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; " +
                    "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88,-2:-88,-2:-90,-1:-90,-1:-91; " +
                    "-1:-90, -1:-89, 0:-89, 0:-90; " +
                    "1:-91, 0:-91, 0:-90, 1:-90"
        )
    }


    @Test
    fun PolylineEnteringRectangle() {
        // A polyline that enters a rectangle very close to one of its vertices.
        val options = roundToE(1)
        val a = "# 0:0, 2:2 #"
        val b = "# # 1:1, 1:3, 3:3, 3:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 1:1, 2:2 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 0:0, 1:1 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1"
        )
    }


    @Test
    fun PolylineCrossingRectangleTwice() {
        // A polyline that crosses a rectangle in one direction, then moves to a
        // different side and crosses the rectangle in the other direction.  Note
        // that an extra vertex is added where the two polyline edges cross.
        val options = roundToE(1)
        val a = "# 0:-5, 0:5, 5:0, -5:0 #"
        val b = "# # 1:1, 1:-1, -1:-1, -1:1"
        expectResult(
            OpType.UNION, options, a, b,
            "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 " +
                    "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1"
        )
        expectResult(
            OpType.INTERSECTION, options, a, b,
            "# 0:-1, 0:0, 0:1 | 1:0, 0:0, -1:0 #"
        )
        expectResult(
            OpType.DIFFERENCE, options, a, b,
            "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 #"
        )
        expectResult(
            OpType.SYMMETRIC_DIFFERENCE, options, a, b,
            "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 " +
                    "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1"
        )
    }

    // Subtracts a degenerate loop along the 180 degree meridian from the given
// input geometry, and compares the result to "expected_str".  The inputs should
// be in the format expected by s2textformat::MakeIndex().
    private fun testMeridianSplitting(input_str: String, expected_str: String) {
        val input = S2TextParser.makeIndex(input_str)
        val meridian = MutableS2ShapeIndex()
        val loops = listOf(
            listOf(
                S2Point(0, 0, -1), S2Point(-1, 0, 0), S2Point(0, 0, 1),
                S2Point(-1, 0, 0)
            )
        )
        meridian.add(S2LaxPolygonShape(loops))
        val output = MutableS2ShapeIndex()
        val layers = mutableListOf(
            IndexedS2PointVectorLayer(output),
            // TODO(ericv): Implement s2builderutil::IndexedS2LaxPolylineVectorLayer.
            IndexedS2PolylineVectorLayer(output),
            IndexedLaxPolygonLayer(output)
        )
        val op = S2BooleanOperation(OpType.DIFFERENCE, layers)
        val error = S2Error()
        assertThat(op.build(input, meridian, error))
            .withFailMessage(error.toString())
            .isTrue()
        assertThat(S2TextParser.toString(output)).isEqualTo(expected_str)
    }

    // This test demonstrated that S2 geometry can easily be transformed such that
// no edge crosses the 180 degree meridian, as required by formats such as
// GeoJSON, by simply subtracting a degenerate loop that follows the 180 degree
// meridian.  This not only splits polylines along the meridian, it also inserts
// the necessary extra vertices at the north/south poles.  (The only extra step
// is that the vertices along the 180 degree meridian or at the poles may need
// to be "doubled" into two vertices, one at longitude 180 and one at longitude
// -180, in order to match the longitudes of the adjacent vertices.)

    @Test
    fun meridianSplitting() {
        // A line along the equator crossing the 180 degree meridian.
        testMeridianSplitting("# 0:-160, 0:170 #", "# 0:-160, 0:180, 0:170 #")

        // The northern hemisphere.
        testMeridianSplitting(
            "# # 0:0, 0:120, 0:-120",
            "# # 90:0, 0:180, 0:-119.99999999999999, 0:0, 0:119.99999999999999, 0:180"
        )

        // A small square that crosses the 180th meridian.  Notice that one input
        // loop is split into two output loops.
        testMeridianSplitting(
            "# # 9:179, 9:-179, 10:-179, 10:179",
            "# # 9.00134850712993:180, 9:-179, 10:-179, 10.00149252698408:180; " +
            "10.00149252698408:180, 10:179, 9:179, 9.00134850712993:180");

        // An annulus that crosses the 180th meridian.  This turns into two shells.
        testMeridianSplitting(
            "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
            "# # 10.00149252698408:180, 10:-179, 9:-179, 9.00134850712993:180, " +
                    "8.00481316618607:180, 8:-178, 11:-178, 11.00654129428001:180; " +
                    "9.00134850712993:180, 9:179, 10:179, 10.00149252698408:180, " +
                    "11.00654129428001:180, 11:178, 8:178, 8.00481316618607:180"
        )

        // An annulus that crosses the 180th meridian.  This turns into two shells.
        testMeridianSplitting(
            "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
            "# # 10.00149252698408:180, 10:-179, 9:-179, 9.00134850712993:180, 8.00481316618607:180, " +
                    "8:-178, 11:-178, 11.00654129428001:180; 9.00134850712993:180, 9:179, 10:179, " +
                    "10.00149252698408:180, 11.00654129428001:180, 11:178, 8:178, 8.00481316618607:180"
        )

        // The whole world except for a small square that crosses the 180th meridian.
        // This is a single loop that visits both poles.  The result is correct
        // except that (1) +180 or -180 needs to be chosen consistently with the
        // adjacent points, and (2) each pole needs to be duplicated (once with
        // longitude -180 and once with longitude 180).
        testMeridianSplitting(
            "# # 9:-179, 9:179, 10:179, 10:-179",
            "# # 0:180, 9.00134850712993:180, 9:179, 10:179, 10.00149252698408:180, " +
                    "90:0, 10.00149252698408:180, 10:-179, 9:-179, 9.00134850712993:180, 0:180, -90:0"
        )
    }

    // This test exercises the "special case" documented in
// GraphEdgeClipper::GetCrossedVertexIndex().

    @Test
    @Strictfp
    fun getCrossedVertexIndexBug() {
        // The first two edges (a0, a1) and (b0, b1) of the following polygons cross
        // such that after snapping, the corresponding edge chains are:
        //
        //   a0 a1 -> a0 b0 b1 x a1
        //   b0 b1 -> b0 x b1
        //
        // where "x" is the computed intersection point of (a0, a1) and (b0, b1).
        // Previously there was a bug such that the two edge chains did not choose
        // the same vertex to represent the point where the two chains cross: the
        // (a0, a1) chain chose "x" as the crossing point while the (b0, b1) chain
        // chose "b0".  This has been fixed such that both chains now choose "x".
        // (Both "x" and "b1" happen to be valid choices in this example, but it is
        // essential that both subchains make the same choice.)

        // S2LatLng coordinates are not accurate enough to reproduce this example.
        val aLoops = listOf(
            listOf(
                // 51.5131559470858:-0.130381523356724
                S2Point(0.62233331065911901, -0.0014161759526823048, 0.78275107466533156),
                // 51.5131892038956:-0.130404244210776
                S2Point(0.6223328557578689, -0.0014164217071954736, 0.78275143589379825),
                S2TextParser.makePoint("51.51317:-0.1306")
            )
        )
        val bLoops = listOf(
            listOf(
                // 51.5131559705551:-0.13038153939079
                S2Point(0.62233331033809591, -0.001416176126110953, 0.78275107492024998),
                // 51.5131559705551:-0.130381539390786
                S2Point(0.62233331033809591, -0.0014161761261109063, 0.78275107492025009),
                S2TextParser.makePoint("51.52:-0.12"),
                S2TextParser.makePoint("51.52:-0.14")
            )
        )
        val a = MutableS2ShapeIndex()
        val b = MutableS2ShapeIndex()
        a.add(S2LaxPolygonShape(aLoops))
        b.add(S2LaxPolygonShape(bLoops))
        val actual = S2LaxPolygonShape()
        val options = LaxPolygonLayer.Options()
        options.degenerateBoundaries = LaxPolygonLayer.DegenerateBoundaries.DISCARD
        val op = S2BooleanOperation(OpType.UNION, LaxPolygonLayer(actual, options = options))
        val error = S2Error()
        assertThat(op.build(a, b, error)).withFailMessage(error.toString()).isTrue()
        assertThat(S2TextParser.toString(actual)).isEqualTo(
            "51.51318713547798:-0.13042532888806, " +
                    "51.51317:-0.1306, " +
                    "51.5131559470858:-0.13038152335672, " +
                    "51.51315597055508:-0.13038153939079, " +
                    "51.51315597055508:-0.13038153939079, " +
                    "51.52:-0.12, " +
                    "51.52:-0.14"
        )
    }


    @Test
    fun FullAndEmptyResults() {
        // The following constants are all in s2textformat::MakeLaxPolygon() format.
        val kEmpty = ""
        val kFull = "full"

        // Two complementary shell/hole pairs, together with alternative shells that
        // are slightly smaller or larger than the original.
        val kShell1 = "10:0, 10:10, 20:10"
        val kHole1 = "10:0, 20:10, 10:10"
        val kShell1Minus = "11:2, 11:9, 18:9"
        val kShell1Plus = "9:-2, 9:11, 22:11"
        val kShell2 = "10:20, 10:30, 20:30"
        val kHole2 = "10:20, 20:30, 10:30"

        // The northern and southern hemispheres.
        val kNorthHemi = "0:0, 0:120, 0:-120"
        val kSouthHemi = "0:0, 0:-120, 0:120"
        // These edges deviate from kSouthHemi by slightly more than 1 degree.
        val kSouthHemiPlus = "0.5:0, 0.5:-120, 0.5:120"

        // A shell and hole that cover complementary hemispheres, such that each
        // hemisphere intersects all six S2 cube faces.  There are also alternative
        // shells that are slightly smaller or larger than the original.
        val k6FaceShell1 = "0:-45, 45:0, 45:90, 0:135, -45:180, -45:-90"
        val k6FaceHole1 = "0:-45, -45:-90, -45:180, 0:135, 45:90, 45:0"
        val k6FaceShell1Minus = "-1:-45, 44:0, 44:90, -1:135, -46:180, -46:-90"
        val k6FaceShell1Plus = "1:-45, 46:0, 46:90, 1:135, -44:180, -44:-90"

        // Two complementary shell/hole pairs that are small enough so that they will
        // disappear when the snap radius chosen above is used.
        val kAlmostEmpty1 = "2:0, 2:10, 3:0"
        val kAlmostFull1 = "2:0, 3:0, 2:10"
        val kAlmostEmpty2 = "4:0, 4:10, 5:0"
        val kAlmostFull2 = "4:0, 5:0, 4:10"

        // A polygon that intersects all 6 faces such but snaps to an empty polygon.
        val k6FaceAlmostEmpty1 = "$k6FaceShell1Minus; $k6FaceHole1"

        // Test empty UNION results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.UNION, kEmpty, kEmpty, kEmpty)
        //  - Empty due to snapping, union does not intersect all 6 cube faces.
        expectPolygon(OpType.UNION, kAlmostEmpty1, kAlmostEmpty2, kEmpty)
        //  - Empty due to snapping, union intersects all 6 cube faces.
        expectPolygon(OpType.UNION, k6FaceAlmostEmpty1, k6FaceAlmostEmpty1, kEmpty)

        // Test full UNION results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.UNION, kEmpty, kFull, kFull)
        expectPolygon(OpType.UNION, kEmpty, kFull, kFull)
        expectPolygon(OpType.UNION, kFull, kFull, kFull)
        //  - Exact result, some input edges.
        expectPolygon(OpType.UNION, kFull, kShell1, kFull)
        expectPolygon(OpType.UNION, kHole1, kHole2, kFull)
        expectPolygon(OpType.UNION, kHole1, kShell1, kFull)
        //  - Full due to snapping, almost complementary polygons.
        expectPolygon(OpType.UNION, kHole1, kShell1Minus, kFull)
        expectPolygon(OpType.UNION, k6FaceHole1, k6FaceShell1Minus, kFull)

        // Test empty INTERSECTION results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.INTERSECTION, kEmpty, kEmpty, kEmpty)
        expectPolygon(OpType.INTERSECTION, kEmpty, kFull, kEmpty)
        expectPolygon(OpType.INTERSECTION, kFull, kEmpty, kEmpty)
        //  - Exact result, inputs do not both intersect all 6 cube faces.
        expectPolygon(OpType.INTERSECTION, kEmpty, kHole1, kEmpty)
        expectPolygon(OpType.INTERSECTION, kShell1, kShell2, kEmpty)
        expectPolygon(OpType.INTERSECTION, kShell1, kHole1, kEmpty)
        //  - Exact result, inputs both intersect all 6 cube faces.
        expectPolygon(OpType.INTERSECTION, k6FaceShell1, k6FaceHole1, kEmpty)
        //  - Empty due to snapping, inputs do not both intersect all 6 cube faces.
        expectPolygon(OpType.INTERSECTION, kShell1Plus, kHole1, kEmpty)
        //  - Empty due to snapping, inputs both intersect all 6 cube faces.
        expectPolygon(OpType.INTERSECTION, k6FaceShell1Plus, k6FaceHole1, kEmpty)

        // Test full INTERSECTION results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.INTERSECTION, kFull, kFull, kFull)
        //  - Full due to snapping, almost full input polygons.
        expectPolygon(OpType.INTERSECTION, kAlmostFull1, kAlmostFull2, kFull)

        // Test empty DIFFERENCE results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.DIFFERENCE, kEmpty, kEmpty, kEmpty)
        expectPolygon(OpType.DIFFERENCE, kEmpty, kFull, kEmpty)
        expectPolygon(OpType.DIFFERENCE, kFull, kFull, kEmpty)
        //  - Exact result, first input does not intersect all 6 cube faces.
        expectPolygon(OpType.DIFFERENCE, kEmpty, kShell1, kEmpty)
        expectPolygon(OpType.DIFFERENCE, kShell1, kFull, kEmpty)
        expectPolygon(OpType.DIFFERENCE, kShell1, kShell1, kEmpty)
        expectPolygon(OpType.DIFFERENCE, kShell1, kHole2, kEmpty)
        //  - Exact result, first input intersects all 6 cube faces.
        expectPolygon(OpType.DIFFERENCE, k6FaceShell1, k6FaceShell1Plus, kEmpty)
        //  - Empty due to snapping, first input does not intersect all 6 cube faces.
        expectPolygon(OpType.DIFFERENCE, kShell1Plus, kShell1, kEmpty)
        //  - Empty due to snapping, first input intersect all 6 cube faces.
        expectPolygon(OpType.DIFFERENCE, k6FaceShell1Plus, k6FaceShell1, kEmpty)

        // Test full DIFFERENCE results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.DIFFERENCE, kFull, kEmpty, kFull)
        //  - Full due to snapping, almost full/empty input polygons.
        expectPolygon(OpType.DIFFERENCE, kAlmostFull1, kAlmostEmpty2, kFull)

        // Test empty SYMMETRIC_DIFFERENCE results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kEmpty, kEmpty)
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kFull, kEmpty)
        //  - Exact result, union does not intersect all 6 cube faces.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kShell1, kEmpty)
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kNorthHemi, kEmpty)
        //  - Exact result, union intersects all 6 cube faces.  This case is only
        //    handled correctly due to the kBiasTowardsEmpty heuristic.
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceShell1,
            kEmpty
        )
        //  - Empty due to snapping, union does not intersect all 6 cube faces.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kShell1, kEmpty)
        //  - Empty due to snapping, union intersects all 6 cube faces.  This case is
        //    only handled correctly due to the kBiasTowardsEmpty heuristic.
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceShell1,
            kEmpty
        )
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceShell1,
            kEmpty
        )

        // Test full SYMMETRIC_DIFFERENCE results.
        //  - Exact result, no input edges.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kEmpty, kFull)
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kFull, kFull)
        //  - Exact result, complementary input polygons.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kHole1, kFull)
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, kAlmostEmpty1, kAlmostFull1,
            kFull
        )
        //  - Full due to snapping, almost complementary input polygons.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kHole1, kFull)
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, kAlmostFull1, kAlmostEmpty2,
            kFull
        )
        //  - Exact result, complementary hemispheres, at least one input does not
        //    intersect all 6 cube faces.
        expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemi, kFull)
        //  - Exact result, almost complementary hemispheres, at least one input does
        //    not intersect all 6 cube faces.
        expectPolygon(
            OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemiPlus,
            kFull
        )

        // TODO(ericv): The following case is not currently implemented.
        //  - Full result, complementary (to within the snap radius) input polygons
        //    each with an area of approximately 2*Pi, and both polygons intersect all
        //    6 cube faces.
/*
  expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceHole1, kFull)
  expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceHole1,
                kFull)
  expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceHole1,
                kFull)
*/
    }

    // Tests whether the two S2ShapeIndexes are equal according to
    // S2BooleanOperation::Equals().
    private fun testEqual(a_str: String, b_str: String): Boolean {
        val a = S2TextParser.makeIndex(a_str)
        val b = S2TextParser.makeIndex(b_str)
        return S2BooleanOperation.equals(a, b)
    }

    // Tests S2BooleanOperation::Equals, which computes the symmetric difference
    // between two geometries and tests whether the result is empty.
    //
    // This also indirectly tests IsEmpty(), which is used to implement Contains()
    // and Intersects().
    @Test
    fun Equals() {
        assertThat(testEqual("# #", "# #")).isTrue()
        assertThat(testEqual("# # full", "# # full")).isTrue()

        assertThat(testEqual("# #", "# # full")).isFalse()
        assertThat(testEqual("0:0 # #", "# #")).isFalse()
        assertThat(testEqual("0:0 # #", "# # full")).isFalse()
        assertThat(testEqual("# 0:0, 1:1 #", "# #")).isFalse()
        assertThat(testEqual("# 0:0, 1:1 #", "# # full")).isFalse()
        assertThat(testEqual("# # 0:0, 0:1, 1:0 ", "# #")).isFalse()
        assertThat(testEqual("# # 0:0, 0:1, 1:0 ", "# # full")).isFalse()
    }

    // Tests Contains() on empty and full geometries.
    @Test
    fun ContainsEmptyAndFull() {
        val empty = S2TextParser.makeIndex("# #")
        val full = S2TextParser.makeIndex("# # full")
        assertThat(S2BooleanOperation.contains(empty, empty)).isTrue()
        assertThat(S2BooleanOperation.contains(empty, full)).isFalse()
        assertThat(S2BooleanOperation.contains(full, empty)).isTrue()
        assertThat(S2BooleanOperation.contains(full, full)).isTrue()
    }

    // Tests Intersects() on empty and full geometries.
    @Test
    fun IntersectsEmptyAndFull() {
        val empty = S2TextParser.makeIndex("# #")
        val full = S2TextParser.makeIndex("# # full")
        assertThat(S2BooleanOperation.intersects(empty, empty)).isFalse()
        assertThat(S2BooleanOperation.intersects(empty, full)).isFalse()
        assertThat(S2BooleanOperation.intersects(full, empty)).isFalse()
        assertThat(S2BooleanOperation.intersects(full, full)).isTrue()
    }

}
