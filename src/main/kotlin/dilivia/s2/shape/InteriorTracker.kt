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
package dilivia.s2.shape

import dilivia.PreConditions.checkState
import dilivia.collections.lowerBound
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeCrosser


// Given a set of shapes, InteriorTracker keeps track of which shapes contain
// a particular point (the "focus").  It provides an efficient way to move the
// focus from one point to another and incrementally update the set of shapes
// which contain it.  We use this to compute which shapes contain the center
// of every S2CellId in the index, by advancing the focus from one cell center
// to the next.
//
// Initially the focus is at the start of the S2CellId space-filling curve.
// We then visit all the cells that are being added to the MutableS2ShapeIndex
// in increasing order of S2CellId.  For each cell, we draw two edges: one
// from the entry vertex to the center, and another from the center to the
// exit vertex (where "entry" and "exit" refer to the points where the
// space-filling curve enters and exits the cell).  By counting edge crossings
// we can incrementally compute which shapes contain the cell center.  Note
// that the same set of shapes will always contain the exit point of one cell
// and the entry point of the next cell in the index, because either (a) these
// two points are actually the same, or (b) the intervening cells in S2CellId
// order are all empty, and therefore there are no edge crossings if we follow
// this path from one cell to the other.
class InteriorTracker {

    // Constructs the InteriorTracker.  You must call AddShape() for each shape
    // that will be tracked before calling MoveTo() or DrawTo().
    constructor() {

    }

    private var isActive: Boolean = false
    private var a: S2Point = S2Point()
    private var b: S2Point = origin()
    private var nextCellid: S2CellId = S2CellId.begin(S2CellId.kMaxLevel)
    private val crosser = S2EdgeCrosser()
    private val shapeIds = mutableListOf<Int>()

    // Shape ids saved by SaveAndClearStateBefore().  The state is never saved
    // recursively so we don't need to worry about maintaining a stack.
    private val savedIds = mutableSetOf<Int>()

    // Returns the current focus point (see above).
    fun focus(): S2Point { return b }

    // Returns true if any shapes are being tracked.
    fun isActive(): Boolean = isActive

    // Adds a shape whose interior should be tracked.  "is_inside" indicates
    // whether the current focus point is inside the shape.  Alternatively, if
    // the focus point is in the process of being moved (via MoveTo/DrawTo), you
    // can also specify "is_inside" at the old focus point and call TestEdge()
    // for every edge of the shape that might cross the current DrawTo() line.
    // This updates the state to correspond to the new focus point.
    //
    // REQUIRES: shape->dimension() == 2
    fun addShape(shapeId: Int, isInside: Boolean) {
        isActive = true;
        if (isInside) {
            toggleShape(shapeId);
        }
    }

    // Moves the focus to the given point.  This method should only be used when
    // it is known that there are no edge crossings between the old and new
    // focus locations; otherwise use DrawTo().
    fun moveTo(b: S2Point) { this.b = b }

    // Moves the focus to the given point.  After this method is called,
    // TestEdge() should be called with all edges that may cross the line
    // segment between the old and new focus locations.
    fun drawTo(point: S2Point) {
        a = this.b
        this.b = point
        crosser.init(this.a, this.b)
    }

    // Indicates that the given edge of the given shape may cross the line
    // segment between the old and new focus locations (see DrawTo).
    // REQUIRES: shape->dimension() == 2
    fun testEdge(shapeId: Int, edge: Edge) {
        if (crosser.edgeOrVertexCrossing(edge.v0, edge.v1)) {
            toggleShape(shapeId)
        }
    }

    // The set of shape ids that contain the current focus.
    fun shapeIds(): List<Int> { return shapeIds }

    // Indicates that the last argument to MoveTo() or DrawTo() was the entry
    // vertex of the given S2CellId, i.e. the tracker is positioned at the start
    // of this cell.  By using this method together with at_cellid(), the caller
    // can avoid calling MoveTo() in cases where the exit vertex of the previous
    // cell is the same as the entry vertex of the current cell.
    fun setNextCellId(nextCellid: S2CellId) {
        this.nextCellid = nextCellid.rangeMin()
    }

    // Returns true if the focus is already at the entry vertex of the given
    // S2CellId (provided that the caller calls set_next_cellid() as each cell
    // is processed).
    fun atCellId(cellid: S2CellId): Boolean {
        return cellid.rangeMin() == nextCellid
    }

    // Makes an internal copy of the state for shape ids below the given limit,
    // and then clear the state for those shapes.  This is used during
    // incremental updates to track the state of added and removed shapes
    // separately.
    fun saveAndClearStateBefore(limitShapeId: Int): Unit {
        checkState { savedIds.isEmpty() }
        val limit = shapeIds.lowerBound(0, shapeIds.size, limitShapeId)
        savedIds.addAll(shapeIds.subList(0, limit))
        repeat(limit) { shapeIds.removeAt(0) }
    }

    // Restores the state previously saved by SaveAndClearStateBefore().  This
    // only affects the state for shape_ids below "limit_shape_id".
    fun restoreStateBefore(limitShapeId: Int) {
        val limit = shapeIds.lowerBound(0, shapeIds.size, limitShapeId)
        repeat(limit) { shapeIds.removeAt(0)}
        shapeIds.addAll(0, savedIds);
        savedIds.clear()
    }

    // Removes "shape_id" from shape_ids_ if it exists, otherwise insert it.
    private fun toggleShape(shapeId: Int) {
        // Since shape_ids_.size() is typically *very* small (0, 1, or 2), it turns
        // out to be significantly faster to maintain a sorted array rather than
        // using an STL set or btree_set.
        if (shapeIds.isEmpty()) {
            shapeIds.add(shapeId);
        } else if (shapeIds.first() == shapeId) {
            shapeIds.removeAt(0)
        } else {
            var pos = 0
            while (pos < shapeIds.size && shapeIds[pos] < shapeId) {
                if (++pos == shapeIds.size) {
                    shapeIds.add(shapeId)
                    return
                }
            }
            if (shapeIds[pos] == shapeId) {
                shapeIds.removeAt(pos)
            } else {
                shapeIds.add(pos, shapeId)
            }
        }
    }

    // Returns a pointer to the first entry "x" where x >= shape_id.
    private fun lowerBound(shapeId: Int): Iterator<Int> {
        val pos = shapeIds.listIterator()
        while (pos.hasNext() && pos.next() < shapeId) { pos.next() }
        if (pos.hasPrevious()) pos.previous()
        return pos.iterator();
    }

    override fun toString(): String {
        return "InteriorTracker(isActive=$isActive, a=$a, b=$b, nextCellid=$nextCellid, shapeIds=$shapeIds, savedIds=$savedIds)"
    }


    companion object {

        // Returns the initial focus point when the InteriorTracker is constructed
        // (corresponding to the start of the S2CellId space-filling curve).
        fun origin(): S2Point {
            // The start of the S2CellId space-filling curve.
            return S2Coords.faceUvToXyz(0, -1.0, -1.0).normalize()
        }

    }




}
