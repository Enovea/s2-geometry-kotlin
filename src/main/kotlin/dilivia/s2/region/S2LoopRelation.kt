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

import dilivia.PreConditions.requireArgument
import dilivia.s2.S2CellId
import dilivia.s2.S2PaddedCell
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2WedgeRelations
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.index.shape.S2ClippedShape
import dilivia.s2.index.shape.S2CrossingEdgeQuery
import dilivia.s2.index.shape.S2ShapeIndex.RangeIterator
import dilivia.s2.index.shape.S2ShapeIndexCell

// LoopRelation is an abstract class that defines a relationship between two
// loops (Contains, Intersects, or CompareBoundary).
abstract class LoopRelation {

    // Optionally, a_target() and b_target() can specify an early-exit condition
    // for the loop relation.  If any point P is found such that
    //
    //   A.Contains(P) == a_crossing_target() &&
    //   B.Contains(P) == b_crossing_target()
    //
    // then the loop relation is assumed to be the same as if a pair of crossing
    // edges were found.  For example, the Contains() relation has
    //
    //   a_crossing_target() == 0
    //   b_crossing_target() == 1
    //
    // because if A.Contains(P) == 0 (false) and B.Contains(P) == 1 (true) for
    // any point P, then it is equivalent to finding an edge crossing (i.e.,
    // since Contains() returns false in both cases).
    //
    // Loop relations that do not have an early-exit condition of this form
    // should return -1 for both crossing targets.
    abstract fun aCrossingTarget(): Int
    abstract fun bCrossingTarget(): Int

    // Given a vertex "ab1" that is shared between the two loops, return true if
    // the two associated wedges (a0, ab1, b2) and (b0, ab1, b2) are equivalent
    // to an edge crossing.  The loop relation is also allowed to maintain its
    // own internal state, and can return true if it observes any sequence of
    // wedges that are equivalent to an edge crossing.
    abstract fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean

}

// Loop relation for Contains().
class ContainsRelation : LoopRelation() {
    private var foundSharedVertex: Boolean = false

    // If A.Contains(P) == false && B.Contains(P) == true, it is equivalent to
    // having an edge crossing (i.e., Contains returns false).
    override fun aCrossingTarget(): Int = 0
    override fun bCrossingTarget(): Int = 1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        foundSharedVertex = true
        return !S2WedgeRelations.wedgeContains(a0, ab1, a2, b0, b2)
    }

    fun foundSharedVertex(): Boolean = foundSharedVertex

}

// Loop relation for Intersects().
class IntersectsRelation : LoopRelation() {
    private var foundSharedVertex: Boolean = false

    // If A.Contains(P) == true && B.Contains(P) == true, it is equivalent to
    // having an edge crossing (i.e., Intersects returns true).
    override fun aCrossingTarget(): Int = 1
    override fun bCrossingTarget(): Int = 1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        foundSharedVertex = true;
        return S2WedgeRelations.wedgeIntersects(a0, ab1, a2, b0, b2)
    }

    fun foundSharedVertex(): Boolean = foundSharedVertex

}


// Loop relation for CompareBoundary().
class CompareBoundaryRelation(reversedB: Boolean) : LoopRelation() {

    val reverseB: Boolean = reversedB       // True if loop B should be reversed.
    private var foundSharedVertex: Boolean = false  // True if any wedge was processed.
    private var containsEdge: Boolean = false        // True if any edge of B is contained by A.
    private var excludesEdge: Boolean = false        // True if any edge of B is excluded by A.

    // The CompareBoundary relation does not have a useful early-exit condition,
    // so we return -1 for both crossing targets.
    //
    // Aside: A possible early exit condition could be based on the following.
    //   If A contains a point of both B and ~B, then A intersects Boundary(B).
    //   If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
    //   So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
    //   the return value is 0, i.e., Boundary(A) intersects Boundary(B).
    // Unfortunately it isn't worth detecting this situation because by the
    // time we have seen a point in all four intersection regions, we are also
    // guaranteed to have seen at least one pair of crossing edges.
    override fun aCrossingTarget(): Int = -1
    override fun bCrossingTarget(): Int = -1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        // Because we don't care about the interior of B, only its boundary, it is
        // sufficient to check whether A contains the semiwedge (ab1, b2).
        foundSharedVertex = true
        if (wedgeContainsSemiwedge(a0, ab1, a2, b2, reverseB)) {
            containsEdge = true
        } else {
            excludesEdge = true
        }
        return containsEdge and excludesEdge
    }

    fun foundSharedVertex(): Boolean = foundSharedVertex

    fun containsEdge(): Boolean = containsEdge

    fun excludesEdge(): Boolean = excludesEdge

}


// Returns true if the wedge (a0, ab1, a2) contains the "semiwedge" defined as
// any non-empty open set of rays immediately CCW from the edge (ab1, b2).  If
// "reverse_b" is true, then substitute "clockwise" for "CCW"; this simulates
// what would happen if the direction of loop B was reversed.
fun wedgeContainsSemiwedge(a0: S2Point, ab1: S2Point, a2: S2Point, b2: S2Point, reverseB: Boolean): Boolean {
    if (b2 == a0 || b2 == a2) {
        // We have a shared or reversed edge.
        return (b2 == a0) == reverseB
    } else {
        return S2Predicates.orderedCCW(a0, a2, b2, ab1)
    }
}


// LoopCrosser is a helper class for determining whether two loops cross.
// It is instantiated twice for each pair of loops to be tested, once for the
// pair (A,B) and once for the pair (B,A), in order to be able to process
// edges in either loop nesting order.
class LoopCrosser(val a: S2Loop, val b: S2Loop, val relation: LoopRelation, val swapped: Boolean) {

    val aCrossingTarget: Int = if (swapped) relation.bCrossingTarget() else relation.aCrossingTarget()
    val bCrossingTarget: Int = if (swapped) relation.aCrossingTarget() else relation.bCrossingTarget()

    // State maintained by StartEdge() and EdgeCrossesCell().
    private val crosser: S2EdgeCrosser = S2EdgeCrosser()
    private var aj: Int = -1
    private var bjPrev: Int = -1

    // Temporary data declared here to avoid repeated memory allocations.
    private val bQuery: S2CrossingEdgeQuery = S2CrossingEdgeQuery(b.index)
    private val bCells = mutableListOf<S2ShapeIndexCell>()

    // Given two dilivia.iterators positioned such that ai->id().Contains(bi->id()),
    // return true if there is a crossing relationship anywhere within ai->id().
    // Specifically, this method returns true if there is an edge crossing, a
    // wedge crossing, or a point P that matches both "crossing targets".
    // Advances both dilivia.iterators past ai->id().
    fun hasCrossingRelation(ai: RangeIterator, bi: RangeIterator): Boolean {
        requireArgument { ai.id().contains(bi.id()) }
        if (ai.numEdges(0) == 0) {
            if (ai.containsCenter(0) && aCrossingTarget == 1) {
                // All points within ai->id() satisfy the crossing target for A, so it's
                // worth iterating through the cells of B to see whether any cell
                // centers also satisfy the crossing target for B.
                do {
                    if (bi.containsCenter(0) && bCrossingTarget == 1) return true
                    bi.next()
                } while (bi.id() <= ai.rangeMax())
            } else {
                // The crossing target for A is not satisfied, so we skip over the cells
                // of B using binary search.
                bi.seekBeyond(ai)
            }
        } else {
            // The current cell of A has at least one edge, so check for crossings.
            if (hasCrossing(ai, bi)) return true
        }
        ai.next()
        return false
    }

    // Given two index cells, return true if there are any edge crossings or
    // wedge crossings within those cells.
    fun cellCrossesCell(aClipped: S2ClippedShape, bClipped: S2ClippedShape): Boolean {
        // Test all edges of "a_clipped" against all edges of "b_clipped".
        val aNumEdges = aClipped.numEdges
        for (i in 0 until aNumEdges) {
            startEdge(aClipped.edge(i))
            if (edgeCrossesCell(bClipped)) return true
        }
        return false
    }

    // Given two dilivia.iterators positioned such that ai->id().Contains(bi->id()),
    // return true if there is an edge crossing or wedge crosssing anywhere
    // within ai->id().  Advances "bi" (only) past ai->id().
    private fun hasCrossing(ai: RangeIterator, bi: RangeIterator): Boolean {
        requireArgument { ai.id().contains(bi.id()) }
        // If ai->id() intersects many edges of B, then it is faster to use
        // S2CrossingEdgeQuery to narrow down the candidates.  But if it intersects
        // only a few edges, it is faster to check all the crossings directly.
        // We handle this by advancing "bi" and keeping track of how many edges we
        // would need to test.
        val kEdgeQueryMinEdges = 20;  // Tuned using benchmarks.
        var totalEdges = 0
        bCells.clear()
        do {
            if (bi.numEdges(0) > 0) {
                totalEdges += bi.cell().numEdges()
                if (totalEdges >= kEdgeQueryMinEdges) {
                    // There are too many edges to test them directly, so use
                    // S2CrossingEdgeQuery.
                    if (cellCrossesAnySubcell(ai.clipped(0), ai.id())) return true
                    bi.seekBeyond(ai)
                    return false;
                }
                bCells.add(bi.cell())
            }
            bi.next()
        } while (bi.id() <= ai.rangeMax())

        // Test all the edge crossings directly.
        for (bCell in bCells) {
            if (cellCrossesCell(ai.clipped(0), bCell.clipped(0))) {
                return true
            }
        }
        return false
    }

    // Given an index cell of A, return true if there are any edge or wedge
    // crossings with any index cell of B contained within "b_id".
    private fun cellCrossesAnySubcell(aClipped: S2ClippedShape, bId: S2CellId): Boolean {
        // Test all edges of "a_clipped" against all edges of B.  The relevant B
        // edges are guaranteed to be children of "b_id", which lets us find the
        // correct index cells more efficiently.
        val bRoot = S2PaddedCell(bId, 0.0)
        val aNumEdges = aClipped.numEdges
        for (i in 0 until aNumEdges) {
            val aj = aClipped.edge(i)
            // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
            // of B that might contain crossing edges.
            bQuery.getCells(a.vertex(aj), a.vertex(aj + 1), bRoot, bCells)
            if (bCells.isEmpty()) continue
            startEdge(aj)
            for (b_cell in bCells) {
                if (edgeCrossesCell(b_cell.clipped(0))) return true
            }
        }
        return false
    }

    // Prepare to check the given edge of loop A for crossings.
    private fun startEdge(aj: Int) {
        // Start testing the given edge of A for crossings.
        crosser.init(a.vertex(aj), a.vertex(aj + 1));
        this.aj = aj;
        bjPrev = -2;
    }

    // Check the current edge of loop A for crossings with all edges of the
    // given index cell of loop B.
    private fun edgeCrossesCell(bClipped: S2ClippedShape): Boolean {
        // Test the current edge of A against all edges of "b_clipped".
        val bNumEdges = bClipped.numEdges
        for (j in 0 until bNumEdges) {
            val bj = bClipped.edge(j)
            if (bj != bjPrev + 1) crosser.restartAt(b.vertex(bj))
            bjPrev = bj
            val crossing = crosser.crossingSign(b.vertex(bj + 1))
            if (crossing < 0) continue
            if (crossing > 0) return true
            // We only need to check each shared vertex once, so we only
            // consider the case where a_vertex(aj_+1) == b_.vertex(bj+1).

            val vertexAjPlus1 = a.vertex(aj + 1)
            val vertexBjPlus1 = b.vertex(bj + 1)
            if (vertexAjPlus1 == vertexBjPlus1) {
                val vertexAj = a.vertex(aj)
                val vertexAjPlus2 = a.vertex(aj + 2)
                val vertexBj = b.vertex(bj)
                val vertexBjPlus2 = b.vertex(bj + 2)
                if (swapped) {
                    if (relation.wedgesCross(vertexBj, vertexBjPlus1, vertexBjPlus2, vertexAj, vertexAjPlus2)) {
                        return true
                    }
                } else if (relation.wedgesCross(vertexAj, vertexAjPlus1, vertexAjPlus2, vertexBj, vertexBjPlus2)) {
                    return true
                }
            }
        }
        return false
    }

}
