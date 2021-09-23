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

import dilivia.PreConditions
import dilivia.s2.S2CellId
import dilivia.s2.S2PaddedCell
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.index.CrossingType
import dilivia.s2.shape.ShapeEdge
import mu.KotlinLogging

// IndexCrosser is a helper class for finding the edge crossings between a
// pair of S2ShapeIndexes.  It is instantiated twice, once for the index pair
// (A,B) and once for the index pair (B,A), in order to be able to test edge
// crossings in the most efficient order.
// @constructor
// If "swapped" is true, the loops A and B have been swapped.  This affects
// how arguments are passed to the given loop relation, since for example
// A.Contains(B) is not the same as B.Contains(A).
class IndexCrosser(val aIndex: S2ShapeIndex, val bIndex: S2ShapeIndex, type: CrossingType, val visitor: EdgePairVisitor, val swapped: Boolean) {

    private val minCrossingSign: Int = if (type == CrossingType.INTERIOR) 1 else 0

    // Temporary data declared here to avoid repeated memory allocations.
    private val bQuery: S2CrossingEdgeQuery = S2CrossingEdgeQuery(bIndex)
    private val bCells = mutableListOf<S2ShapeIndexCell>()
    private val aShapeEdges: ShapeEdgeVector = mutableListOf()
    private val bShapeEdges: ShapeEdgeVector = mutableListOf()

    // Given two dilivia.iterators positioned such that ai->id().Contains(bi->id()),
    // visits all crossings between edges of A and B that intersect a->id().
    // Terminates early and returns false if visitor_ returns false.
    // Advances both dilivia.iterators past ai->id().
    fun visitCrossings(ai: S2ShapeIndex.RangeIterator, bi: S2ShapeIndex.RangeIterator): Boolean {
        logger.trace { "visitCrossings | ai = ${ai.id()}, bi = ${bi.id()}" }
        PreConditions.requireArgument { ai.id().contains(bi.id()) }
        if (ai.cell().numEdges() == 0) {
            // Skip over the cells of B using binary search.
            bi.seekBeyond(ai)
        } else {
            // If ai->id() intersects many edges of B, then it is faster to use
            // S2CrossingEdgeQuery to narrow down the candidates.  But if it
            // intersects only a few edges, it is faster to check all the crossings
            // directly.  We handle this by advancing "bi" and keeping track of how
            // many edges we would need to test.
            val kEdgeQueryMinEdges = 23
            var bEdges = 0
            bCells.clear()
            do {
                val cellEdges = bi.cell().numEdges()
                if (cellEdges > 0) {
                    bEdges += cellEdges
                    if (bEdges >= kEdgeQueryMinEdges) {
                        // There are too many edges, so use an S2CrossingEdgeQuery.
                        if (!visitSubcellCrossings(ai.cell(), ai.id())) return false
                        bi.seekBeyond(ai)
                        return true
                    }
                    bCells.add(bi.cell())
                }
                bi.next()
            } while (bi.id() <= ai.rangeMax())
            if (bCells.isNotEmpty()) {
                // Test all the edge crossings directly.
                S2CrossingEdgePairsScanner.getShapeEdges(aIndex, ai.cell(), aShapeEdges)
                S2CrossingEdgePairsScanner.getShapeEdges(bIndex, bCells, bShapeEdges)
                if (!visitEdgesEdgesCrossings(aShapeEdges, bShapeEdges)) {
                    return false
                }
            }
        }
        ai.next()
        return true
    }

    // Given two index cells, visits all crossings between edges of those cells.
    // Terminates early and returns false if visitor_ returns false.
    fun visitCellCellCrossings(aCell: S2ShapeIndexCell, bCell: S2ShapeIndexCell): Boolean {
        // Test all edges of "a_cell" against all edges of "b_cell".
        S2CrossingEdgePairsScanner.getShapeEdges(aIndex, aCell, aShapeEdges)
        S2CrossingEdgePairsScanner.getShapeEdges(bIndex, bCell, bShapeEdges)
        return visitEdgesEdgesCrossings(aShapeEdges, bShapeEdges)
    }

    private fun visitEdgePair(a: ShapeEdge, b: ShapeEdge, isInterior: Boolean): Boolean {
        if (swapped) {
            return visitor.visit(b, a, isInterior)
        } else {
            return visitor.visit(a, b, isInterior)
        }
    }

    // Visits all crossings of the current edge with all edges of the given index
    // cell of B.  Terminates early and returns false if visitor_ returns false.
    private fun visitEdgeCellCrossings(a: ShapeEdge, bCell: S2ShapeIndexCell): Boolean {
        // Test the current edge of A against all edges of "b_cell".

        // Note that we need to use a new S2EdgeCrosser (or call Init) whenever we
        // replace the contents of b_shape_edges_, since S2EdgeCrosser requires that
        // its S2Point arguments point to values that persist between Init() calls.
        S2CrossingEdgePairsScanner.getShapeEdges(bIndex, bCell, bShapeEdges)
        val crosser = S2EdgeCrosser(a.v0, a.v1)
        for (b in bShapeEdges) {
            if (crosser.c() == null || crosser.c() != b.v0) {
                crosser.restartAt(b.v0)
            }
            val sign = crosser.crossingSign(b.v1)
            if (sign >= minCrossingSign && !visitEdgePair(a, b, sign == 1)) return false
        }
        return true
    }

    // Visits all crossings of any edge in "a_cell" with any index cell of B that
    // is a descendant of "b_id".  Terminates early and returns false if
    // visitor_ returns false.
    private fun visitSubcellCrossings(aCell: S2ShapeIndexCell, bId: S2CellId): Boolean {
        // Test all edges of "a_cell" against the edges contained in B index cells
        // that are descendants of "b_id".
        S2CrossingEdgePairsScanner.getShapeEdges(aIndex, aCell, aShapeEdges)
        val bRoot = S2PaddedCell(bId, 0.0)
        for (a in aShapeEdges) {
            // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
            // of B that might contain crossing edges.
            if (!bQuery.visitCells(a.v0, a.v1, bRoot, object : CellVisitor {
                        override fun visit(cell: S2ShapeIndexCell): Boolean {
                            return visitEdgeCellCrossings(a, cell)
                        }

                    })) {
                return false
            }
        }
        return true
    }

    // Visits all crossings of any edge in "a_edges" with any edge in "b_edges".
    private fun visitEdgesEdgesCrossings(aEdges: ShapeEdgeVector, bEdges: ShapeEdgeVector): Boolean {
        // Test all edges of "a_edges" against all edges of "b_edges".
        for (a in aEdges) {
            val crosser = S2EdgeCrosser(a.v0, a.v1)
            for (b in bEdges) {
                if (crosser.c() == null || crosser.c() != b.v0) {
                    crosser.restartAt(b.v0)
                }
                val sign = crosser.crossingSign(b.v1)
                if (sign >= minCrossingSign && !visitEdgePair(a, b, sign == 1)) return false
            }
        }
        return true
    }

    companion object {
        private val logger = KotlinLogging.logger(IndexCrosser::class.java.name)
    }

}
