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

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireEQ
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.math.R2Rect
import dilivia.math.vectors.R2Point
import dilivia.s2.S2CellId
import dilivia.s2.S2PaddedCell
import dilivia.s2.S2Point
import dilivia.s2.edge.FaceSegmentList
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.edge.S2EdgeClipping.kFaceClipErrorUVCoord
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.index.CrossingType
import dilivia.s2.index.shape.S2ShapeIndex.CellIterator
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdge
import dilivia.s2.shape.ShapeEdgeId


// A function that is called with each candidate intersecting edge.  The
// function may return false in order to request that the algorithm should
// be terminated, i.e. no further crossings are needed.
@FunctionalInterface
interface ShapeEdgeIdVisitor {
    fun visit(id: ShapeEdgeId): Boolean
}

// A function that is called with each S2ShapeIndexCell that might contain
// edges intersecting the given query edge.  The function may return false
// in order to request that the algorithm should be terminated, i.e. no
// further crossings are needed.
@FunctionalInterface
interface CellVisitor {
    fun visit(cell: S2ShapeIndexCell): Boolean
}

// S2CrossingEdgeQuery is used to find edges or shapes that are crossed by
// an edge.  Here is an example showing how to index a set of polylines,
// and then find the polylines that are crossed by a given edge AB:
//
// void Test(const vector<S2Polyline*>& polylines,`
//           const S2Point& a0, const S2Point &a1) {
//   MutableS2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(absl::make_unique<S2Polyline::Shape>(polyline));
//   }
//   S2CrossingEdgeQuery query(&index);
//   for (const auto& edge : query.GetCrossingEdges(a, b, CrossingType::ALL)) {
//     S2_CHECK_GE(S2::CrossingSign(a0, a1, edge.v0(), edge.v1()), 0);
//   }
// }
//
// Note that if you need to query many edges, it is more efficient to declare
// a single S2CrossingEdgeQuery object and reuse it so that temporary storage
// does not need to be reallocated each time.
//
// If you want to find *all* pairs of crossing edges, use
// s2shapeutil::VisitCrossingEdgePairs() instead.
// REQUIRES: "index" is not modified after this method is called.
class S2CrossingEdgeQuery(val index: S2ShapeIndex) {

    private val iter: CellIterator = index.cellIterator()

    //////////// Temporary storage used while processing a query ///////////
    lateinit var a0_: R2Point
    lateinit var a1_: R2Point

    // Avoids repeated allocation when methods are called many times.
    private val tmp_candidates: MutableList<ShapeEdgeId> = mutableListOf()

    // Returns all edges that intersect the given query edge (a0,a1) and that
    // have the given CrossingType (ALL or INTERIOR).  Edges are sorted and
    // unique.
    fun getCrossingEdges(a0: S2Point, a1: S2Point, type: CrossingType): List<ShapeEdge> {
        val edges = mutableListOf<ShapeEdge>();
        getCrossingEdges(a0, a1, type, edges)
        return edges;
    }

    // A specialized version of GetCrossingEdges() that only returns the edges
    // that belong to a particular S2Shape. 
    fun getCrossingEdges(a0: S2Point, a1: S2Point, shape: S2Shape, type: CrossingType): List<ShapeEdge> {
        val edges = mutableListOf<ShapeEdge>()
        getCrossingEdges(a0, a1, shape, type, edges)
        return edges
    }

    // These versions can be more efficient when they are called many times,
    // since they do not require allocating a new vector on each call.
    fun getCrossingEdges(a0: S2Point, a1: S2Point, type: CrossingType, edges: MutableList<ShapeEdge>): Unit {
        edges.clear()
        getCandidates(a0, a1, tmp_candidates)
        val min_sign = if (type == CrossingType.ALL) 0 else 1
        val crosser = S2EdgeCrosser(a0, a1)
        var shape_id = -1
        lateinit var shape: S2Shape
        for (candidate in tmp_candidates) {
            if (candidate.shapeId != shape_id) {
                shape_id = candidate.shapeId
                shape = index.shape(shape_id)!!
            }
            val edge_id = candidate.edgeId
            val b = shape.edge(edge_id)
            if (crosser.crossingSign(b.v0, b.v1) >= min_sign) {
                edges.add(ShapeEdge(shape_id, edge_id, b))
            }
        }
    }

    fun getCrossingEdges(
        a0: S2Point,
        a1: S2Point,
        shape: S2Shape,
        type: CrossingType,
        edges: MutableList<ShapeEdge>
    ): Unit {
        edges.clear()
        getCandidates(a0, a1, shape, tmp_candidates)
        val min_sign = if (type == CrossingType.ALL) 0 else 1
        val crosser = S2EdgeCrosser(a0, a1)
        for (candidate in tmp_candidates) {
            val edge_id = candidate.edgeId
            val b = shape.edge(edge_id)
            if (crosser.crossingSign(b.v0, b.v1) >= min_sign) {
                edges.add(ShapeEdge(shape.id, edge_id, b))
            }
        }
    }

    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.

    // Returns a superset of the edges that intersect a query edge (a0, a1).
    // This method is useful for clients that want to test intersections in some
    // other way, e.g. using S2::EdgeOrVertexCrossing().
    fun getCandidates(a0: S2Point, a1: S2Point): List<ShapeEdgeId> {
        val edges = mutableListOf<ShapeEdgeId>()
        getCandidates(a0, a1, edges)
        return edges
    }

    // A specialized version of GetCandidates() that only returns the edges that
    // belong to a particular S2Shape.
    fun getCandidates(a0: S2Point, a1: S2Point, shape: S2Shape): List<ShapeEdgeId> {
        val edges = mutableListOf<ShapeEdgeId>()
        getCandidates(a0, a1, shape, edges)
        return edges
    }

    // These versions can be more efficient when they are called many times,
    // since they do not require allocating a new vector on each call.
    fun getCandidates(a0: S2Point, a1: S2Point, edges: MutableList<ShapeEdgeId>) {
        edges.clear()
        val num_edges = S2CountEdges.countEdgesUpTo(index, kMaxBruteForceEdges + 1);
        visitRawCandidates(a0, a1, object : ShapeEdgeIdVisitor {
            override fun visit(id: ShapeEdgeId): Boolean {
                edges.add(id)
                return true
            }
        })
        edges.sortAndRemoveDuplicates()
    }

    fun getCandidates(a0: S2Point, a1: S2Point, shape: S2Shape, edges: MutableList<ShapeEdgeId>) {
        edges.clear()
        val num_edges = shape.numEdges
        visitRawCandidates(a0, a1, shape, object : ShapeEdgeIdVisitor {
            override fun visit(id: ShapeEdgeId): Boolean {
                edges.add(id)
                return true
            }
        })
        edges.sortAndRemoveDuplicates()
    }

    // Visits a superset of the edges that intersect the query edge (a0, a1),
    // terminating early if the given ShapeEdgeIdVisitor returns false (in which
    // case this function returns false as well).
    //
    // CAVEAT: Edges may be visited more than once.
    fun visitRawCandidates(a0: S2Point, a1: S2Point, visitor: ShapeEdgeIdVisitor): Boolean {
        val num_edges = S2CountEdges.countEdgesUpTo(index, kMaxBruteForceEdges + 1);
        if (num_edges <= kMaxBruteForceEdges) {
            val num_shape_ids = index.nextNewShapeId()
            for (s in 0 until num_shape_ids) {
                val shape = index.shape(s) ?: continue
                val num_shape_edges = shape.numEdges
                for (e in 0 until num_shape_edges) {
                    if (!visitor.visit(ShapeEdgeId(s, e))) return false
                }
            }
            return true;
        }
        return visitCells(a0, a1, object : CellVisitor {
            override fun visit(cell: S2ShapeIndexCell): Boolean {
                for (s in 0 until cell.numClipped) {
                    val clipped = cell.clipped(s)
                    for (j in 0 until clipped.numEdges) {
                        if (!visitor.visit(ShapeEdgeId(clipped.shapeId, clipped.edge(j)))) {
                            return false
                        }
                    }
                }
                return true
            }
        })
    }

    fun visitRawCandidates(a0: S2Point, a1: S2Point, shape: S2Shape, visitor: ShapeEdgeIdVisitor): Boolean {
        val num_edges = shape.numEdges
        if (num_edges <= kMaxBruteForceEdges) {
            for (e in 0 until num_edges) {
                if (!visitor.visit(ShapeEdgeId(shape.id, e))) return false
            }
            return true;
        }
        return visitCells(a0, a1, object : CellVisitor {
            override fun visit(cell: S2ShapeIndexCell): Boolean {
                val clipped = cell.findClipped(shape.id) ?: return true
                for (j in 0 until clipped.numEdges) {
                    if (!visitor.visit(ShapeEdgeId(shape.id, clipped.edge(j)))) return false
                }
                return true
            }
        })
    }


    // Visits all S2ShapeIndexCells that might contain edges intersecting the
    // given query edge (a0, a1), terminating early if the given CellVisitor
    // returns false (in which case this function returns false as well).
    //
    // NOTE: Each candidate cell is visited exactly once.
    fun visitCells(a0: S2Point, a1: S2Point, visitor: CellVisitor): Boolean {
        val segments: FaceSegmentList = mutableListOf()
        S2EdgeClipping.getFaceSegments(a0, a1, segments)
        for (segment in segments) {
            a0_ = segment.a;
            a1_ = segment.b;

            // Optimization: rather than always starting the recursive subdivision at
            // the top level face cell, instead we start at the smallest S2CellId that
            // contains the edge (the "edge root cell").  This typically lets us skip
            // quite a few levels of recursion since most edges are short.
            val edge_bound = R2Rect.fromPointPair(a0_, a1_)
            var pcell = S2PaddedCell(S2CellId.fromFace(segment.face), 0.0)
            val edge_root = pcell.shrinkToFit(edge_bound)

            // Now we need to determine how the edge root cell is related to the cells
            // in the spatial index (cell_map_).  There are three cases:
            //
            //  1. edge_root is an index cell or is contained within an index cell.
            //     In this case we only need to look at the contents of that cell.
            //  2. edge_root is subdivided into one or more index cells.  In this case
            //     we recursively subdivide to find the cells intersected by a0a1.
            //  3. edge_root does not intersect any index cells.  In this case there
            //     is nothing to do.
            val relation = iter.locate(edge_root)
            if (relation == CellRelation.INDEXED) {
                // edge_root is an index cell or is contained by an index cell (case 1).
                checkState { iter.id().contains(edge_root) }
                if (!visitor.visit(iter.cell())) return false
            } else if (relation == CellRelation.SUBDIVIDED) {
                // edge_root is subdivided into one or more index cells (case 2).  We
                // find the cells intersected by a0a1 using recursive subdivision.
                if (!edge_root.isFace) pcell = S2PaddedCell(edge_root, 0.0);
                if (!visitCells(pcell, edge_bound, visitor)) return false
            }
        }
        return true
    }

    // Visits all S2ShapeIndexCells within "root" that might contain edges
    // intersecting the given query edge (a0, a1), terminating early if the
    // given CellVisitor returns false (in which case this function returns
    // false as well).
    //
    // NOTE: Each candidate cell is visited exactly once.
    //
    // REQUIRES: root.padding() == 0
    //   [This low-level method does not support padding; the argument is supplied
    //    as an S2PaddedCell in order to avoid constructing it repeatedly when
    //    this method is called using different query edges with the same root.]
    fun visitCells(a0: S2Point, a1: S2Point, root: S2PaddedCell, visitor: CellVisitor): Boolean {
        requireEQ(root.padding, 0.0)
        // We use padding when clipping to ensure that the result is non-empty
        // whenever the edge (a0, a1) intersects the given root cell.
        val clippedEdge = S2EdgeClipping.clipToPaddedFace(a0, a1, root.id.face(), kFaceClipErrorUVCoord)
        if (clippedEdge != null) {
            a0_ = clippedEdge.first
            a1_ = clippedEdge.second
            val edge_bound = R2Rect.fromPointPair(a0_, a1_)
            if (root.bound.intersects(edge_bound)) {
                return visitCells(root, edge_bound, visitor)
            }
        }
        return true
    }

    // Given a query edge AB and a cell "root", returns all S2ShapeIndex cells
    // within "root" that might contain edges intersecting AB.
    //
    // REQUIRES: root.padding() == 0 (see above)
    fun getCells(a0: S2Point, a1: S2Point, root: S2PaddedCell, cells: MutableList<S2ShapeIndexCell>) {
        cells.clear()
        visitCells(a0, a1, root, object : CellVisitor {
            override fun visit(cell: S2ShapeIndexCell): Boolean {
                cells.add(cell)
                return true
            }
        })
    }

    // Computes the index cells intersected by the current edge that are
    // descendants of "pcell" and calls visitor_ for each one.
    //
    // WARNING: This function is recursive with a maximum depth of 30.  The frame
    // size is about 2K in versions of GCC prior to 4.7 due to poor overlapping
    // of storage for temporaries.  This is fixed in GCC 4.7, reducing the frame
    // size to about 350 bytes (i.e., worst-case total stack usage of about 10K).
    private fun visitCells(pcell: S2PaddedCell, edgeBound: R2Rect, visitor: CellVisitor): Boolean {
        // This code uses S2PaddedCell because it has the methods we need for
        // efficient splitting, however the actual padding is required to be zero.
        requireEQ(pcell.padding, 0.0)

        iter.seek(pcell.id.rangeMin())
        if (iter.done() || iter.id() > pcell.id.rangeMax()) {
            // The index does not contain "pcell" or any of its descendants.
            return true
        }
        if (iter.id() == pcell.id) {
            return visitor.visit(iter.cell())
        }

        // Otherwise, split the edge among the four children of "pcell".
        val center = pcell.middle().lo
        if (edgeBound[0].hi < center[0]) {
            // Edge is entirely contained in the two left children.
            return clipVAxis(edgeBound, center[1], 0, pcell, visitor);
        } else if (edgeBound[0].lo >= center[0]) {
            // Edge is entirely contained in the two right children.
            return clipVAxis(edgeBound, center[1], 1, pcell, visitor);
        } else {
            val childBounds = Array(2) { R2Rect.empty() }
            splitUBound(edgeBound, center[0], childBounds);
            if (edgeBound[1].hi < center[1]) {
                // Edge is entirely contained in the two lower children.
                return (visitCells(S2PaddedCell(pcell, 0, 0), childBounds[0], visitor) &&
                        visitCells(S2PaddedCell(pcell, 1, 0), childBounds[1], visitor))
            } else if (edgeBound[1].lo >= center[1]) {
                // Edge is entirely contained in the two upper children.
                return (visitCells(S2PaddedCell(pcell, 0, 1), childBounds[0], visitor) &&
                        visitCells(S2PaddedCell(pcell, 1, 1), childBounds[1], visitor))
            } else {
                // The edge bound spans all four children.  The edge itself intersects
                // at most three children (since no padding is being used).
                return (clipVAxis(childBounds[0], center[1], 0, pcell, visitor) &&
                        clipVAxis(childBounds[1], center[1], 1, pcell, visitor));
            }
        }
    }

    // Given either the left (i=0) or right (i=1) side of a padded cell "pcell",
    // determine whether the current edge intersects the lower child, upper child,
    // or both children, and call VisitCells() recursively on those children.
    // "center" is the v-coordinate at the center of "pcell".
    private fun clipVAxis(
        edgeBound: R2Rect,
        center: Double,
        i: Int,
        pcell: S2PaddedCell,
        visitor: CellVisitor
    ): Boolean {
        if (edgeBound[1].hi < center) {
            // Edge is entirely contained in the lower child.
            return visitCells(S2PaddedCell(pcell, i, 0), edgeBound, visitor);
        } else if (edgeBound[1].lo >= center) {
            // Edge is entirely contained in the upper child.
            return visitCells(S2PaddedCell(pcell, i, 1), edgeBound, visitor);
        } else {
            // The edge intersects both children.
            val childBounds = Array<R2Rect>(2) { R2Rect.empty() }
            splitVBound(edgeBound, center, childBounds);
            return (visitCells(S2PaddedCell(pcell, i, 0), childBounds[0], visitor) &&
                    visitCells(S2PaddedCell(pcell, i, 1), childBounds[1], visitor));
        }
    }

    // Split the current edge into two child edges at the given u-value "u" and
    // return the bound for each child.
    fun splitUBound(edgeBound: R2Rect, u: Double, childBounds: Array<R2Rect>) {
        // See comments in MutableS2ShapeIndex::ClipUBound.
        val v = edgeBound[1].project(S2EdgeClipping.interpolateDouble(u, a0_[0], a1_[0], a0_[1], a1_[1]))

        // "diag_" indicates which diagonal of the bounding box is spanned by a0a1:
        // it is 0 if a0a1 has positive slope, and 1 if a0a1 has negative slope.
        val diag = if ((a0_[0] > a1_[0]) != (a0_[1] > a1_[1])) 1 else 0
        splitBound(edgeBound, 0, u, diag, v, childBounds);
    }

    // Split the current edge into two child edges at the given v-value "v" and
    // return the bound for each child.
    fun splitVBound(edge_bound: R2Rect, v: Double, child_bounds: Array<R2Rect>) {
        val u = edge_bound[0].project(S2EdgeClipping.interpolateDouble(v, a0_[1], a1_[1], a0_[0], a1_[0]))
        val diag = if ((a0_[0] > a1_[0]) != (a0_[1] > a1_[1])) 1 else 0
        splitBound(edge_bound, diag, u, 0, v, child_bounds)
    }

    companion object {

        // For small loops it is faster to use brute force.  The threshold below was
        // determined using the benchmarks in the unit test.
        val kMaxBruteForceEdges = 27

        // Split the current edge into two child edges at the given point (u,v) and
        // return the bound for each child.  "u_end" and "v_end" indicate which bound
        // endpoints of child 1 will be updated.
        fun splitBound(edge_bound: R2Rect, u_end: Int, u: Double, v_end: Int, v: Double, child_bounds: Array<R2Rect>) {
            val childBound0 = edge_bound.clone()
            childBound0[0][1 - u_end] = u;
            childBound0[1][1 - v_end] = v;
            child_bounds[0] = childBound0
            checkState { !child_bounds[0].isEmpty }
            checkState { edge_bound.contains(child_bounds[0]) }

            val childBound1 = edge_bound.clone()
            childBound1[0][u_end] = u;
            childBound1[1][v_end] = v;
            child_bounds[1] = childBound1;
            checkState { !child_bounds[1].isEmpty }
            checkState { edge_bound.contains(child_bounds[1]) }
        }
    }

}
