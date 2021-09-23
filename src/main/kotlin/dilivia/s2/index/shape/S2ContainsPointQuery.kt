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

import dilivia.s2.S2Point
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.index.shape.S2ShapeIndex.CellIterator
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdge
import mu.KotlinLogging

// Defines whether shapes are considered to contain their vertices.  Note that
// these definitions differ from the ones used by S2BooleanOperation.
//
//  - In the OPEN model, no shapes contain their vertices (not even points).
//    Therefore Contains(S2Point) returns true if and only if the point is
//    in the interior of some polygon.
//
//  - In the SEMI_OPEN model, polygon point containment is defined such that
//    if several polygons tile the region around a vertex, then exactly one of
//    those polygons contains that vertex.  Points and polylines still do not
//    contain any vertices.
//
//  - In the CLOSED model, all shapes contain their vertices (including points
//    and polylines).
//
// Note that points other than vertices are never contained by polylines.
// If you want need this behavior, use S2ClosestEdgeQuery::IsDistanceLess()
// with a suitable distance threshold instead.
enum class S2VertexModel { OPEN, SEMI_OPEN, CLOSED };

// This class defines the options supported by S2ContainsPointQuery.
data class S2ContainsPointQueryOptions(

        // Controls whether shapes are considered to contain their vertices (see
        // definitions above).  By default the SEMI_OPEN model is used.
        //
        // DEFAULT: S2VertexModel::SEMI_OPEN
        var vertex_model: S2VertexModel = S2VertexModel.SEMI_OPEN) {

}

// S2ContainsPointQuery determines whether one or more shapes in an
// S2ShapeIndex contain a given S2Point.  The S2ShapeIndex may contain any
// number of points, polylines, and/or polygons (possibly overlapping).
// Shape boundaries may be modeled as OPEN, SEMI_OPEN, or CLOSED (this affects
// whether or not shapes are considered to contain their vertices).
//
// Example usage:
//   auto query = MakeS2ContainsPointQuery(&index, S2VertexModel::CLOSED);
//   return query.Contains(point);
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
//
// However, note that if you need to do a large number of point containment
// tests, it is more efficient to re-use the S2ContainsPointQuery object
// rather than constructing a new one each time.
class S2ContainsPointQuery<T : S2ShapeIndex> {

    private lateinit var index: T
    private lateinit var options: S2ContainsPointQueryOptions
    private lateinit var iter: CellIterator

    private val logger = KotlinLogging.logger {  }

    // Default constructor; requires Init() to be called.
    constructor() {
    }

    // Rather than calling this constructor, which requires specifying the
    // IndexType template argument explicitly, the preferred idiom is to call
    // MakeS2ContainsPointQuery() instead.  For example:
    //
    //   return MakeS2ContainsPointQuery(&index).Contains(p);
    constructor(index: T, options: S2ContainsPointQueryOptions = S2ContainsPointQueryOptions()) {
        init(index, options)
    }

    // Convenience constructor that accepts the S2VertexModel directly.
    constructor(index: T, vertex_model: S2VertexModel) {
        init(index, S2ContainsPointQueryOptions(vertex_model))
    }

    fun index() = index
    fun options() = options

    // Equivalent to the two-argument constructor above.
    fun init(index: T, options: S2ContainsPointQueryOptions = S2ContainsPointQueryOptions()) {
        this.index = index
        this.options = options
        this.iter = index.cellIterator()
    }

    // Returns true if any shape in the given index() contains the point "p"
    // under the vertex model specified (OPEN, SEMI_OPEN, or CLOSED).
    fun contains(p: S2Point): Boolean {
        if (!iter.locate(p)) return false

        val cell = iter.cell()
        val numClipped = cell.numClipped
        for (s in 0 until numClipped) {
            if (shapeContains(iter, cell.clipped(s), p)) return true
        }
        return false
    }

    // Returns true if the given shape contains the point "p" under the vertex
    // model specified (OPEN, SEMI_OPEN, or CLOSED).
    //
    // REQUIRES: "shape" belongs to index().
    fun shapeContains(shape: S2Shape, p: S2Point): Boolean {
        if (!iter.locate(p)) return false
        val clipped = iter.cell().findClipped(shape.id) ?: return false
        return shapeContains(iter, clipped, p)
    }

    // Visits all shapes in the given index() that contain the given point "p",
    // terminating early if the given ShapeVisitor function returns false (in
    // which case VisitContainingShapes returns false as well).  Each shape is
    // visited at most once.
    //
    // Note that the API allows non-const access to the visited shapes.
    fun interface ShapeVisitor {
        fun visit(shape: S2Shape): Boolean
    }

    fun visitContainingShapes(p: S2Point, visitor: ShapeVisitor): Boolean {
        // This function returns "false" only if the algorithm terminates early
        // because the "visitor" function returned false.
        if (!iter.locate(p)) return true

        val cell = iter.cell()
        val numClipped = cell.numClipped
        for (s in 0 until numClipped) {
            val clipped = cell.clipped(s)
            val shape = index.shape(clipped.shapeId)
            val contains = shape?.let { shapeContains(iter, clipped, p) } ?: false
            logger.trace { "visitContainingShapes | test shape: $shape: contains = $contains" }
            if (contains && !visitor.visit(shape!!)) {
                return false
            }
        }
        return true
    }

    // Convenience function that returns all the shapes that contain the given
    // point "p".
    fun getContainingShapes(p: S2Point): List<S2Shape> {
        val results = mutableListOf<S2Shape>()
        visitContainingShapes(p, object : ShapeVisitor {
            override fun visit(shape: S2Shape): Boolean {
                results.add(shape)
                return true
            }
        })
        return results
    }

    // Visits all edges in the given index() that are incident to the point "p"
    // (i.e., "p" is one of the edge endpoints), terminating early if the given
    // EdgeVisitor function returns false (in which case VisitIncidentEdges
    // returns false as well).  Each edge is visited at most once.
    fun interface EdgeVisitor {
        fun visit(edge: ShapeEdge): Boolean
    }

    fun visitIncidentEdges(p: S2Point, visitor: EdgeVisitor): Boolean {
        // This function returns "false" only if the algorithm terminates early
        // because the "visitor" function returned false.
        if (!iter.locate(p)) return true

        val cell = iter.cell()
        val numClipped = cell.numClipped
        for (s in 0 until numClipped) {
            val clipped = cell.clipped(s)
            val numEdges = clipped.numEdges
            if (numEdges == 0) continue
            val shape = index.shape(clipped.shapeId) ?: continue
            for (i in 0 until numEdges) {
                val edgeId = clipped.edge(i)
                val edge = shape.edge(edgeId)
                if ((edge.v0 == p || edge.v1 == p) && !visitor.visit(ShapeEdge(shape.id, edgeId, edge))) {
                    return false
                }
            }
        }
        return true
    }

    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use.

    // Returns a pointer to the iterator used internally by this class, in order
    // to avoid the need for clients to create their own iterator.  Clients are
    // allowed to reposition this iterator arbitrarily between method calls.
    fun mutableIter() = iter

    // Low-level helper method that returns true if the given S2ClippedShape
    // referred to by an S2ShapeIndex::Iterator contains the point "p".
    fun shapeContains(iter: CellIterator, clipped: S2ClippedShape, p: S2Point): Boolean {
        var inside = clipped . containsCenter
        val num_edges = clipped.numEdges
        if (num_edges > 0) {
            val shape = index.shape(clipped.shapeId)!!
            if (shape.dimension < 2) {
                // Points and polylines can be ignored unless the vertex model is CLOSED.
                if (options.vertex_model != S2VertexModel.CLOSED) return false

                // Otherwise, the point is contained if and only if it matches a vertex.
                for (i in 0 until num_edges) {
                    val edge = shape . edge (clipped.edge(i))
                    if (edge.v0 == p || edge.v1 == p) return true
                }
                return false
            }
            // Test containment by drawing a line segment from the cell center to the
            // given point and counting edge crossings.
            val crosser = S2EdgeCrosser(iter.center(), p)
            for (i in 0 until num_edges) {
                val edge = shape . edge (clipped.edge(i))
                var sign = crosser . crossingSign (edge.v0, edge.v1)
                if (sign < 0) continue
                if (sign == 0) {
                    // For the OPEN and CLOSED models, check whether "p" is a vertex.
                    if (options.vertex_model != S2VertexModel.SEMI_OPEN && (edge.v0 == p || edge.v1 == p)) {
                        return (options.vertex_model == S2VertexModel.CLOSED)
                    }
                    sign = if (S2EdgeCrossings.vertexCrossing(crosser.a(), crosser.b(), edge.v0, edge.v1)) 1 else 0
                }
                inside = inside xor (sign > 0)
            }
        }
        return inside
    }

    companion object {
        // Returns an S2ContainsPointQuery for the given S2ShapeIndex.  Note that
        // it is efficient to return S2ContainsPointQuery objects by value.
        fun <T : S2ShapeIndex> makeS2ContainsPointQuery(index: T, options: S2ContainsPointQueryOptions = S2ContainsPointQueryOptions()): S2ContainsPointQuery<T> {
            return S2ContainsPointQuery<T>(index, options)
        }
    }

}
