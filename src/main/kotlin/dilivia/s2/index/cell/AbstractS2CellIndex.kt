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
package dilivia.s2.index.cell

import dilivia.ComparisonChain
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.s2.S2CellId
import dilivia.s2.index.lexicon.Label
import dilivia.s2.region.S2CellUnion

abstract class AbstractS2CellIndex {

    abstract fun cellIterator(): CellIterator

    /**
     * Represents a node in the (cell_id, value) tree. Cells are organized in a tree such that the ancestors of a given
     * node contain that node.
     *
     * @property cellId Cell id of this node.
     * @property label Associated value.
     * @property parent Index of the ancestor node.
     */
    data class CellNode(
        val cellId: S2CellId = S2CellId.none,
        val label: Label = S2CellIndex.kDoneContents,
        val parent: Int = -1
    )

    /**
     * A RangeNode represents a range of leaf S2CellIds. The range starts at "startId" (a leaf cell) and ends at the
     * "startId" field of the next RangeNode. "contents" points to the node of cellTree representing the cells that
     * overlap this range.
     *
     * @property startId First leaf cell contained by this range.
     * @property contents Contents of this node (an index within cellTree).
     */
    data class RangeNode(
        val startId: S2CellId,
        val contents: Int
    ) : Comparable<Any?> {

        override fun compareTo(other: Any?): Int {
            return when (other) {
                is S2CellId -> compareToCell(other)
                is RangeNode -> compareToRange(other)
                else -> throw IllegalArgumentException("")
            }

        }

        fun compareToCell(other: S2CellId): Int = startId.compareTo(other)

        fun compareToRange(other: RangeNode): Int = ComparisonChain.start()
            .compare(startId, other.startId)
            .compare(contents, other.contents)
            .result()

    }

    // To build the cell tree and leaf cell ranges, we maintain a stack of (cell_id, value) pairs that contain the
    // current leaf cell.  This class represents an instruction to push or pop a (cell_id, value) pair.
    //
    // If label >= 0, the (cell_id, value) pair is pushed on the stack.
    // If cell_id == S2CellId::Sentinel(), a pair is popped from the stack.
    // Otherwise the stack is unchanged but a RangeNode is still emitted.
    data class Delta(val startId: S2CellId, val cellId: S2CellId, val value: Label, val idx: Int = 0) : Comparable<Delta> {

        // Deltas are sorted first by start_id, then in reverse order by cell_id,
        // and then by label.  This is necessary to ensure that
        //   (1) larger cells are pushed on the stack before smaller cells, and
        //   (2) cells are popped off the stack before any new cells are added.

        override fun compareTo(other: Delta): Int {
            return ComparisonChain.start()
                .compare(startId, other.startId)
                .compare(other.cellId, cellId)
                .compare(value, other.value)
                .result()
        }

    }

    /**
     * Clears the index so that it can be re-used.
     */
    abstract fun clear()

    /**
     * Visits all (cellId, value) pairs in the given index that intersect the given S2CellUnion "target", terminating
     * early if the given CellVisitor function returns false (in which case visitIntersectingCells returns false
     * as well). Each (cellId, value) pair in the index is visited at most once. (If the index contains duplicates,
     * then each copy is visited.)
     *
     * @param target Target cell union.
     * @param visitor Cell visitor.
     * @return true if all the intersecting cells have been visited and false if the visitor exit the processing.
     */
    abstract fun visitIntersectingCells(target: S2CellUnion, visitor: CellVisitor): Boolean

    /**
     * Convenience function that returns the labels of all indexed cells that intersect the given S2CellUnion "target".
     * The output contains each label at most once, but is not sorted.
     *
     * @param target A cell union.
     * @return The values that intersects the specified target.
     */
    fun getIntersectingValues(target: S2CellUnion): List<Label> = getIntersectingValues(target, mutableListOf())

    /**
     * Convenience function that returns the labels of all indexed cells that intersect the given S2CellUnion "target".
     *
     * This version can be more efficient when it is called many times, since it does not require allocating a new list
     * on each call.
     *
     * @param target A cell union.
     * @param values The list to populate
     * @return The values that intersects the specified target.
     */
    fun getIntersectingValues(target: S2CellUnion, values: MutableList<Label>): MutableList<Label> {
        values.clear()
        visitIntersectingCells(target, object : CellVisitor {

            override fun visit(cellId: S2CellId, label: Label): Boolean {
                values.add(label)
                return true
            }
        })
        values.sortAndRemoveDuplicates()
        return values
    }

    interface CellIterator {

        // The S2CellId of the current (cell_id, value) pair.
        // REQUIRES: !done()
        fun cellId(): S2CellId

        // The Label of the current (cell_id, value) pair.
        // REQUIRES: !done()
        fun value(): Label

        // Returns the current (cell_id, value) pair.
        fun valuedCell(): ValuedCell<Label>

        // Returns true if all (cell_id, value) pairs have been visited.
        fun done(): Boolean

        // Advances the iterator to the next (cell_id, value) pair.
        // REQUIRES: !done()
        fun next()

    }
}
