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

import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.index.shape.CellRelation
import dilivia.s2.index.shape.S2ClippedShape
import dilivia.s2.index.shape.S2ContainsPointQuery
import dilivia.s2.index.shape.S2ShapeIndex

// This class wraps an S2ShapeIndex object with the additional methods needed
// to implement the S2Region API, in order to allow S2RegionCoverer to compute
// S2CellId coverings of arbitrary collections of geometry.
//
// These methods could conceivably be made part of S2ShapeIndex itself, but
// there are several advantages to having a separate class:
//
//  - The class can be templated in order to avoid virtual calls and memory
//    allocation (for dilivia.iterators) when the concrete S2ShapeIndex type is known.
//
//  - Implementing these methods efficiently requires an S2ShapeIndex iterator,
//    and this design allows a single iterator to be allocated and reused.
//
//  - S2Region::Clone() is not a good fit for the S2ShapeIndex API because
//    it can't be implemented for some subtypes (e.g., EncodedS2ShapeIndex).
//
// Example usage:
//
// S2CellUnion GetCovering(const S2ShapeIndex& index) {
//   S2RegionCoverer coverer;
//   coverer.mutable_options()->set_max_cells(20);
//   S2CellUnion covering;
//   coverer.GetCovering(MakeS2ShapeIndexRegion(&index), &covering);
//   return covering;
// }
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
class S2ShapeIndexRegion<T : S2ShapeIndex>(val index: T) : S2Region {

    // This class is not thread-safe!
    private val containsQuery: S2ContainsPointQuery<T> = S2ContainsPointQuery(index)

    // Optimization: rather than declaring our own iterator, instead we reuse
    // the iterator declared by S2ContainsPointQuery.  (This improves benchmark
    // times significantly for classes that create a new S2ShapeIndexRegion
    // object on every call to Contains/MayIntersect(S2Cell).
    private val iter = containsQuery.mutableIter()

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    override fun clone() = S2ShapeIndexRegion(index)

    override val capBound: S2Cap
        get() {
            val covering = mutableListOf<S2CellId>()
            getCellUnionBound(covering)
            return S2CellUnion(covering).capBound
        }

    override val rectBound: S2LatLngRect
        get() {
            val covering = mutableListOf<S2CellId>()
            getCellUnionBound(covering)
            return S2CellUnion(covering).rectBound
        }


    // This method currently returns at most 4 cells, unless the index spans
    // multiple faces in which case it may return up to 6 cells.
    override fun getCellUnionBound(cellIds: MutableList<S2CellId>) {
        // We find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  There are
        // two cases:
        //
        //  - If the index intersects two or more faces, then for each intersected
        //    face we add one cell to the covering.  Rather than adding the entire
        //    face, instead we add the smallest S2Cell that covers the S2ShapeIndex
        //    cells within that face.
        //
        //  - If the index intersects only one face, then we first find the smallest
        //    cell S that contains the index cells (just like the case above).
        //    However rather than using the cell S itself, instead we repeat this
        //    process for each of its child cells.  In other words, for each
        //    child cell C we add the smallest S2Cell C' that covers the index cells
        //    within C.  This extra step is relatively cheap and produces much
        //    tighter coverings when the S2ShapeIndex consists of a small region
        //    near the center of a large S2Cell.
        //
        // The following code uses only a single Iterator object because creating an
        // Iterator may be relatively expensive for some S2ShapeIndex types (e.g.,
        // it may involve memory allocation).
        cellIds.clear()
        if (cellIds is ArrayList) {
            cellIds.ensureCapacity(6)
        }

        // Find the last S2CellId in the index.
        iter.finish()
        if (!iter.prev()) return  // Empty index.
        val lastIndexId = iter.id()
        iter.begin()
        if (iter.id() != lastIndexId) {
            // The index has at least two cells.  Choose an S2CellId level such that
            // the entire index can be spanned with at most 6 cells (if the index
            // spans multiple faces) or 4 cells (it the index spans a single face).
            val level = iter.id().getCommonAncestorLevel(lastIndexId) + 1

            // For each cell C at the chosen level, we compute the smallest S2Cell
            // that covers the S2ShapeIndex cells within C.
            val lastId = lastIndexId.parent(level)
            var id = iter.id().parent(level)
            while (id != lastId) {
                // If the cell C does not contain any index cells, then skip it.
                if (id.rangeMax() < iter.id()) {
                    id = id.next()
                    continue
                }

                // Find the range of index cells contained by C and then shrink C so
                // that it just covers those cells.
                val first = iter.id()
                iter.seek(id.rangeMax().next())
                iter.prev()
                coverRange(first, iter.id(), cellIds)
                iter.next()
                id = id.next()
            }
        }
        coverRange(iter.id(), lastIndexId, cellIds)
    }

    // Returns true if "target" is contained by any single shape.  If the cell
    // is covered by a union of different shapes then it may return false.
    //
    // The implementation is conservative but not exact; if a shape just barely
    // contains the given cell then it may return false.  The maximum error is
    // less than 10 * DBL_EPSILON radians (or about 15 nanometers).
    override fun contains(cell: S2Cell): Boolean {
        val relation = iter.locate(cell.id())

        // If the relation is DISJOINT, then "target" is not contained.  Similarly if
        // the relation is SUBDIVIDED then "target" is not contained, since index
        // cells are subdivided only if they (nearly) intersect too many edges.
        if (relation != CellRelation.INDEXED) return false

        // Otherwise, the iterator points to an index cell containing "target".
        // If any shape contains the target cell, we return true.
        check(iter.id().contains(cell.id()))
        val currentShapeIndexCell = iter.cell()
        for (s in 0 until currentShapeIndexCell.numClipped) {
            val clipped = currentShapeIndexCell.clipped(s)
            // The shape contains the target cell iff the shape contains the cell
            // center and none of its edges intersects the (padded) cell interior.
            if (iter.id() == cell.id()) {
                if (clipped.numEdges == 0 && clipped.containsCenter) return true
            } else {
                // It is faster to call AnyEdgeIntersects() before Contains().
                val shape = index.shape(clipped.shapeId) ?: continue
                if (shape.dimension == 2 && !anyEdgeIntersects(clipped, cell) && containsQuery.shapeContains(iter, clipped, cell.getCenter())) {
                    return true
                }
            }
        }
        return false
    }

    // Returns true if any shape intersects "target".
    //
    // The implementation is conservative but not exact; if a shape is just
    // barely disjoint from the given cell then it may return true.  The maximum
    // error is less than 10 * DBL_EPSILON radians (or about 15 nanometers).
    override fun mayIntersect(cell: S2Cell): Boolean {
        val relation = iter.locate(cell.id())

        // If "target" does not overlap any index cell, there is no intersection.
        if (relation == CellRelation.DISJOINT) return false

        // If "target" is subdivided into one or more index cells, then there is an
        // intersection to within the S2ShapeIndex error bound.
        if (relation == CellRelation.SUBDIVIDED) return true

        // Otherwise, the iterator points to an index cell containing "target".
        //
        // If "target" is an index cell itself, there is an intersection because index
        // cells are created only if they have at least one edge or they are
        // entirely contained by the loop.
        check(iter.id().contains(cell.id()))
        if (iter.id() == cell.id()) return true

        // Test whether any shape intersects the target cell or contains its center.
        val c = iter.cell()
        for (s in 0 until c.numClipped) {
            val clipped = c.clipped(s)
            if (anyEdgeIntersects(clipped, cell)) return true
            if (containsQuery.shapeContains(iter, clipped, cell.getCenter())) {
                return true
            }
        }
        return false
    }

    // Returns true if the given point is contained by any two-dimensional shape
    // (i.e., polygon).  Boundaries are treated as being semi-open (i.e., the
    // same rules as S2Polygon).  Zero and one-dimensional shapes are ignored by
    // this method (if you need more flexibility, see S2BooleanOperation).

    override fun contains(p: S2Point): Boolean {
        if (iter.locate(p)) {
            val cell = iter.cell()
            for (s in 0 until cell.numClipped) {
                if (containsQuery.shapeContains(iter, cell.clipped(s), p)) {
                    return true
                }
            }
        }
        return false
    }

    // Returns true if the indexed shape "clipped" in the indexed cell "id"
    // contains the point "p".
    //
    // REQUIRES: id.contains(S2CellId(p))
    private fun contains(id: S2CellId, clipped: S2ClippedShape, p: S2Point): Boolean = TODO()

    // Returns true if any edge of the indexed shape "clipped" intersects the
    // cell "target".  It may also return true if an edge is very close to
    // "target"; the maximum error is less than 10 * DBL_EPSILON radians (about
    // 15 nanometers).
    fun anyEdgeIntersects(clipped: S2ClippedShape, target: S2Cell): Boolean {
        val kMaxError = (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist)
        val bound = target.boundUV().expanded(kMaxError)
        val face = target.face()
        val shape = index.shape(clipped.shapeId) ?: return false
        val numEdges = clipped.numEdges
        for (i in 0 until numEdges) {
            val edge = shape.edge(clipped.edge(i))
            val clippedEdge = S2EdgeClipping.clipToPaddedFace(edge.v0, edge.v1, face, kMaxError)
            if (clippedEdge != null && S2EdgeClipping.intersectsRect(clippedEdge.first, clippedEdge.second, bound)) {
                return true
            }
        }
        return false
    }


    companion object {


        // Computes the smallest S2Cell that covers the S2Cell range (first, last) and
        // adds this cell to "cell_ids".
        //
        // REQUIRES: "first" and "last" have a common ancestor.
        private fun coverRange(first: S2CellId, last: S2CellId, cellIds: MutableList<S2CellId>) {
            if (first == last) {
                // The range consists of a single index cell.
                cellIds.add(first)
            } else {
                // Add the lowest common ancestor of the given range.
                val level = first.getCommonAncestorLevel(last)
                check(level >= 0)
                cellIds.add(first.parent(level))
            }
        }

        // Returns an S2ShapeIndexRegion that wraps the given S2ShapeIndex.  Note that
        // it is efficient to return S2ShapeIndexRegion objects by value.
        fun <T : S2ShapeIndex> makeS2ShapeIndexRegion(index: T): S2ShapeIndexRegion<T> = S2ShapeIndexRegion(index)
    }

}
