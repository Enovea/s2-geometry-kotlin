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

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.collections.upperBound
import dilivia.s2.S2CellId
import dilivia.s2.index.lexicon.Label
import dilivia.s2.region.S2CellUnion
import mu.KotlinLogging

/**
 * S2CellIndex stores a collection of (cellId, value) pairs. The S2CellIds may be overlapping or contain duplicate
 * values.  For example, an S2CellIndex could store a collection of S2CellUnions, where each S2CellUnion has its own
 * value.
 *
 * Value can be 32-bit non-negative integers, and be typically used to map the results of queries back to client data
 * structures. See ValueLexicon, which maintains a set of distinct labels (Int) and maps them to sequentially numbered
 * integers.  For example, the following code uses strings as labels:
 *
 *   val myLabelLexicon = ValueLexicon<string>
 *   var labelStr = ...
 *   var label = cellIndex.add(cellId, myLabelLexicon.add(labelStr))
 *   ...
 *   var label = ...
 *   var labelStr = myLabelLexicon.value(label)
 *
 * To build an S2CellIndex, call add() for each (cellId, value) pair, and then call the build() method.
 * For example:
 *
 *  val contents = listOf<S2CellId>(...)
 *  for (int i = 0; i < contents.size(); ++i) {
 *      index.add(contents[i], i /*label*/)
 *  }
 *  index.build()
 *
 * There is also a convenience method that adds an S2CellUnion:
 *
 *  index.add(cellUnion, value)
 *
 * Note that the index is not dynamic; the contents of the index cannot be changed once it has been built.
 *
 * There are several options for retrieving data from the index. The simplest is to use a built-in method such as
 * getIntersectingLabels (which returns the labels of all cells that intersect a given target S2CellUnion):
 *
 *  val values: List<T> = index.getIntersectingValues(targetUnion)
 *
 * Alternatively, you can use an external class such as S2ClosestCellQuery, which computes the cell(s) that are closest
 * to a given target geometry.
 * For example, here is how to find all cells that are closer than "distanceLimit" to a given target point:
 *
 *   val query = S2ClosestCellQuery(index)
 *   query.options.setMaxDistance(distanceLimit)
 *   val target = S2ClosestCellQuery.PointTarget(target_point)
 *   for (result in query.findClosestCells(target)) {
 *     // result.distance is the distance to the target.
 *     // result.cellId is the indexed S2CellId.
 *     // result.label is the integer label associated with the S2CellId.
 *     doSomething(targetPoint, result)
 *   }
 *
 * Finally, you can access the index contents directly. Internally, the index consists of a set of non-overlapping
 * leaf cell ranges that subdivide the sphere and such that each range intersects a particular set of (cell_id, label)
 * pairs.
 * Data is accessed using the following iterator types:
 *
 *  RangeIterator:
 *    - used to seek and iterate over the non-overlapping leaf cell ranges.
 *  NonEmptyRangeIterator:
 *    - like RangeIterator, but skips ranges whose contents are empty.
 *  ContentsIterator:
 *    - iterates over the (cell_id, value) pairs that intersect a given range.
 *  CellIterator:
 *    - iterates over the entire set of (cell_id, value) pairs.
 *
 * Note that these are low-level, efficient types intended mainly for implementing new query classes. Most clients
 * should use either the built-in methods such as visitIntersectingCells and getIntersectingLabels,
 * or a helper such as S2ClosestCellQuery or S2Closest*Query*.CellUnionTarget.
 *
 * This class is a port of the S2CellIndex class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T The type of indexed data.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2CellIndex : AbstractS2CellIndex() {

    /**
     * A tree of (cellId, value) pairs such that if X is an ancestor of Y, then X.cell_id contains Y.cell_id.
     * The contents of a given range of leaf cells can be represented by pointing to a node of this tree.
     */
    private val cellTree = mutableListOf<CellNode>()

    /**
     * The last element of rangeNodes is a sentinel value, which is necessary in order to represent the range covered
     * by the previous element.
     */
    private val rangeNodes = ArrayList<RangeNode>()

    /** The number of (cell_id, value) pairs in the index. */
    val numCells: Int
        get() = cellTree.size

    override fun cellIterator(): CellIterator = CellIterator(this)

    /**
     * Adds the given (cell_id, value) pair to the index. Note that the index is not valid until build() is called.
     *
     * The S2CellIds in the index may overlap (including duplicate values). Duplicate (cell_id, value) pairs are also
     * allowed, although be aware that S2ClosestCellQuery will eliminate such duplicates anyway.
     *
     * REQUIRES: cellId.isValid
     *
     * @param cellId A valid cell id.
     * @param value The data to index.
     */
    fun add(cellId: S2CellId, value: Label) {
        requireArgument { cellId.isValid }
        cellTree.add(CellNode(cellId, value, -1))
    }

    /**
     * Convenience function that adds a collection of cells with the same label.
     *
     * @param cellIds A valid cell union.
     * @param value The data to index.
     */
    fun add(cellIds: S2CellUnion, value: Label) {
        cellIds.forEach { cellId -> add(cellId, value) }
    }

    /**
     * Constructs the index.
     * This method may only be called once. No dilivia.iterators may be used until the index is built.
     */
    fun build() {
        val startTime = System.currentTimeMillis()
        logger.debug { "Start building index" }

        val deltas = ArrayList<Delta>(2 * cellTree.size + 2)
        // Create two deltas for each (cell_id, label) pair: one to add the pair to
        // the stack (at the start of its leaf cell range), and one to remove it from
        // the stack (at the end of its leaf cell range).
        for (node in cellTree) {
            deltas.add(Delta(node.cellId.rangeMin(), node.cellId, node.label))
            deltas.add(Delta(node.cellId.rangeMax().next(), S2CellId.sentinel, -1))
        }
        // We also create two special deltas to ensure that a RangeNode is emitted at
        // the beginning and end of the S2CellId range.
        deltas.add(Delta(S2CellId.begin(S2CellId.kMaxLevel), S2CellId.none, -1))
        deltas.add(Delta(S2CellId.end(S2CellId.kMaxLevel), S2CellId.none, -1))
        deltas.sort()
        logger.trace { """
            |Deltas: ${deltas.size} elements
            |-----------------------------
            |${deltas.joinToString("\n")}
            |""".trimMargin() }

        // Now walk through the deltas to build the leaf cell ranges and cell tree
        // (which is essentially a permanent form of the "stack" described above).
        cellTree.clear()
        rangeNodes.ensureCapacity(deltas.size)
        var contents = kDoneContents
        var i = 0
        while (i < deltas.size) {
            val startId = deltas[i].startId
            // Process all the deltas associated with the current start_id.
            while (i < deltas.size && deltas[i].startId == startId) {
                if (deltas[i].value >= 0) {
                    cellTree.add(CellNode(cellId = deltas[i].cellId, label = deltas[i].value, parent = contents))
                    contents = cellTree.size - 1
                } else if (deltas[i].cellId == S2CellId.sentinel) {
                    contents = cellTree[contents].parent
                }
                ++i
            }
            rangeNodes.add(RangeNode(startId = startId, contents = contents))
        }

        logger.trace { """Cell tree: ${cellTree.size} elements
            |-----------------------------
            |${cellTree.joinToString("\n")}
            |Range nodes: ${rangeNodes.size} elements
            |-----------------------------
            |${rangeNodes.joinToString("\n")}
            |""".trimMargin()
        }
        logger.debug { "Index build in ${System.currentTimeMillis() - startTime} ms" }
    }

    /**
     * Clears the index so that it can be re-used.
     */
    override fun clear() {
        cellTree.clear()
        rangeNodes.clear()
    }

    fun toDebugString(): String {
        val sep = "\n\t- "
        val cells = cellTree.joinToString(sep, prefix = sep)
        val ranges = rangeNodes.joinToString(sep, prefix = sep)

        return """LevelDbS2CellIndex:
            |Cells:$cells
            |Ranges:$ranges
        """.trimMargin()
    }

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
    override fun visitIntersectingCells(target: S2CellUnion, visitor: CellVisitor): Boolean {
        logger.trace { "--> visitIntersectingCells(target = $target)" }
        if (target.isEmpty()) return true
        val targetIterator = target.listIterator()
        var it = targetIterator.next()
        val contents = ContentsIterator(this)
        val range = RangeIterator(this)
        range.begin()
        do {
            logger.trace {
                """
                |
                |Current target cell: $it [ range min = ${it.rangeMin()} ; range max = ${it.rangeMax()} ]
                |=====================================
                |Current range: 
                |  - start id:  ${range.startId()}
                |  - limit id : ${range.limitId()}
                |---------------------------------------
            """.trimMargin()
            }

            if (range.limitId() <= it.rangeMin()) {
                logger.trace { "Range limit id <= current cell range min id => seek to cell range min = ${it.rangeMin()}" }
                range.seek(it.rangeMin())  // Only seek when necessary.
            }
            while (range.startId() <= it.rangeMax()) {
                logger.trace { "Range start id = ${range.startId()} <= target cell range max = ${it.rangeMax()} => visit range contents" }
                contents.startUnion(range)
                if (!visitRangeContents(contents, visitor)) return false
                range.next()
            }
            // Check whether the next target cell is also contained by the leaf cell
            // range that we just processed.  If so, we can skip over all such cells
            // using binary search.  This speeds up benchmarks by between 2x and 10x
            // when the average number of intersecting cells is small (< 1).
            if (!range.done() && targetIterator.hasNext()) {
                it = targetIterator.next()
                if (it.rangeMax() < range.startId()) {
                    logger.trace { "Next cell ($it) range max = ${it.rangeMax()} < range start id = ${range.startId()}" }
                    // Skip to the first target cell that extends past the previous range.

                    //it = std::lower_bound(it + 1, target.end(), range.start_id());
                    //if ((it - 1)->range_max() >= range.start_id())--it;
                    while (targetIterator.hasNext() && target[targetIterator.nextIndex()] < range.startId()) {
                        it = targetIterator.next()
                    }
                    logger.trace { "Current cell = $it (idx = ${targetIterator.previousIndex() + 1}) after Skip to the first target cell that extends past the previous range" }
                    if (target[targetIterator.previousIndex()].rangeMax() >= range.startId()) {
                        logger.trace { "Previous cell = ${target[targetIterator.previousIndex()]} range max >= range start id = ${range.startId()} => previous" }
                        it = targetIterator.previous()
                    }
                }
            } else break
        } while (true)
        return true
    }

    private fun visitRangeContents(contents: ContentsIterator, visitor: CellVisitor): Boolean {
        while (!contents.done()) {
            logger.trace { "Visit cell ${ValuedCell(cellId = contents.cellId(), label = contents.value())}" }
            if (!visitor.visit(contents.cellId(), contents.value())) {
                return false
            }
            contents.next()
        }
        return true
    }


    override fun toString(): String {
        val sep = "\n\t- "
        return """S2CellIndex:
            |Cells:${
            cellTree.mapIndexed { i, c -> i to c }.joinToString(sep, prefix = sep) { "${it.first}: ${it.second}" }
        }
            |Ranges:${rangeNodes.joinToString(sep, prefix = sep)}
        """.trimMargin()
    }

    /**
     * An iterator that visits the entire set of indexed (cell_id, value) pairs in an unspecified order.
     * Initializes a CellIterator for the given S2CellIndex, positioned at the first cell (if any).
     */
    class CellIterator(val index: S2CellIndex): AbstractS2CellIndex.CellIterator {

        // NOTE(ericv): There is a potential optimization that would require this
        // class to iterate over both cell_tree_ *and* range_nodes_.
        private var cellIdx: Int = -1

        init {
            cellIdx = 0
            checkState({ index.rangeNodes.isNotEmpty() }, { "Call build() first." })
        }

        // The S2CellId of the current (cell_id, value) pair.
        // REQUIRES: !done()
        override fun cellId(): S2CellId {
            checkState({ !done() }, { "End of cell iterator is reached. " })
            return index.cellTree[cellIdx].cellId
        }

        // The Label of the current (cell_id, value) pair.
        // REQUIRES: !done()
        override fun value(): Label {
            checkState({ !done() }, { "End of cell iterator is reached. " })
            return index.cellTree[cellIdx].label
        }

        // Returns the current (cell_id, value) pair.
        override fun valuedCell(): ValuedCell<Label> {
            check(!done()) { "End of cell iterator is reached. " }
            val cellNode = index.cellTree[cellIdx]
            return ValuedCell(cellId = cellNode.cellId, label = cellNode.label)
        }

        // Returns true if all (cell_id, value) pairs have been visited.
        override fun done(): Boolean = cellIdx == index.cellTree.size

        // Advances the iterator to the next (cell_id, value) pair.
        // REQUIRES: !done()
        override fun next(): Unit {
            check(!done()) { "End of cell iterator is reached. " }
            ++cellIdx
        }

    }

    /**
     * An iterator that seeks and iterates over a set of non-overlapping leaf cell ranges that cover the entire sphere.
     * The indexed (s2cell_id, value) pairs that intersect the current leaf cell range can be visited using
     * ContentsIterator (see below). Initializes a RangeIterator for the given S2CellIndex.  The iterator is
     * initially *unpositioned*; you must call a positioning method such as begin() or seek() before accessing its
     * contents.
     */
    open class RangeIterator(val index: S2CellIndex) {

        private var rangeNodeIdx: Int

        init {
            check(index.rangeNodes.isNotEmpty()) { "Call build() first." }
            rangeNodeIdx = kUninitialized
        }

        constructor(iterator: RangeIterator) : this(iterator.index) {
            this.rangeNodeIdx = iterator.rangeNodeIdx
        }

        // The start of the current range of leaf S2CellIds.
        //
        // If done() is true, returns S2CellId.end(S2CellId.kMaxLevel).  This
        // property means that most loops do not need to test done() explicitly.
        fun startId(): S2CellId = index.rangeNodes[rangeNodeIdx].startId

        // The (non-inclusive) end of the current range of leaf S2CellIds.
        // REQUIRES: !done()
        fun limitId(): S2CellId {
            check(!done())
            return index.rangeNodes[rangeNodeIdx + 1].startId
        }

        // Returns true if the iterator is positioned beyond the last valid range.
        fun done(): Boolean {
            check(rangeNodeIdx != kUninitialized) { "Call begin() or seek() first." }

            // Note that the last element of range_nodes_ is a sentinel value.
            return rangeNodeIdx >= index.rangeNodes.lastIndex
        }

        // Positions the iterator at the first range of leaf cells (if any).
        open fun begin() {
            rangeNodeIdx = 0
        }

        // Positions the iterator so that done() is true.
        fun finish(): Unit {
            rangeNodeIdx = index.rangeNodes.lastIndex
        }

        // Advances the iterator to the next range of leaf cells.
        // REQUIRES: !done()
        open fun next() {
            check(!done())
            ++rangeNodeIdx
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        open fun prev(): Boolean {
            if (rangeNodeIdx == 0) return false
            --rangeNodeIdx
            return true
        }

        // Positions the iterator at the first range with start_id() >= target.
        // (Such an entry always exists as long as "target" is a valid leaf cell.
        // Note that it is valid to access start_id() even when done() is true.)
        //
        // REQUIRES: target.is_leaf()
        open fun seek(target: S2CellId) {
            requireArgument { target.isLeaf }
            rangeNodeIdx = index.rangeNodes.upperBound(value = target) - 1
            // std::upper_bound(range_nodes_->begin(), range_nodes_->end(), target)-1;
        }

        // Returns true if no (s2cell_id, label) pairs intersect this range.
        // Also returns true if done() is true.
        fun isEmpty(): Boolean = index.rangeNodes[rangeNodeIdx].contents == kDoneContents

        // If advancing the iterator "n" times would leave it positioned on a
        // valid range, does so and returns true.  Otherwise leaves the iterator
        // unmodified and returns false.
        fun advance(n: Int): Boolean {
            // Note that the last element of range_nodes_ is a sentinel value.
            if (n >= index.rangeNodes.lastIndex - 1 - rangeNodeIdx) return false
            rangeNodeIdx += n
            return true
        }

        fun contents(): Int = index.rangeNodes[rangeNodeIdx].contents

        companion object {
            // A special value used to indicate that the RangeIterator has not yet
            // been initialized by calling Begin() or Seek().
            // Note that since the last element of range_nodes_ is a sentinel value,
            // it_ will never legitimately be positioned at range_nodes_->end().
            const val kUninitialized = Int.MAX_VALUE
        }

    }

    // Like RangeIterator, but only visits leaf cell ranges that overlap at
    // least one (cell_id, label) pair.
    // Initializes a NonEmptyRangeIterator for the given S2CellIndex.
    // The iterator is initially *unpositioned*; you must call a positioning
    // method such as Begin() or Seek() before accessing its contents.
    class NonEmptyRangeIterator : RangeIterator {

        constructor(index: S2CellIndex) : super(index)

        constructor(iterator: NonEmptyRangeIterator) : super(iterator)

        // Positions the iterator at the first non-empty range of leaf cells.
        override fun begin() {
            super.begin()
            while (isEmpty() && !done()) next()
        }

        // Advances the iterator to the next non-empty range of leaf cells.
        // REQUIRES: !done()
        override fun next() {
            do {
                super.next()
            } while (isEmpty() && !done())
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        override fun prev(): Boolean {
            while (super.prev()) {
                if (!isEmpty()) return true
            }
            // Return the iterator to its original position.
            if (isEmpty() && !done()) next()
            return false
        }

        // Positions the iterator at the first non-empty range with
        // start_id() >= target.
        //
        // REQUIRES: target.is_leaf()
        override fun seek(target: S2CellId) {
            super.seek(target)
            while (isEmpty() && !done()) super.next()
        }

    }

    /**
     * An iterator that visits the (cell_id, value) pairs that cover a set of leaf cell ranges (see RangeIterator).
     * Note that when multiple leaf cell ranges are visited, this class only guarantees that each result will be
     * reported at least once, i.e. duplicate values may be suppressed.  If you want duplicate values to be reported
     * again, be sure to call clear() first.
     *
     * [In particular, the implementation guarantees that when multiple leaf cell ranges are visited in monotonically
     * increasing order, then each (cell_id, value) pair is reported exactly once.]
     */
    class ContentsIterator {

        private var index: S2CellIndex? = null

        // The value of it.start_id() from the previous call to StartUnion().
        // This is used to check whether these values are monotonically
        // increasing.
        private var prevStartId: S2CellId = S2CellId.none

        // The maximum index within the cell_tree_ vector visited during the
        // previous call to StartUnion().  This is used to eliminate duplicate
        // values when StartUnion() is called multiple times.
        private var nodeCutoff: Int = -1

        // The maximum index within the cell_tree_ vector visited during the
        // current call to StartUnion().  This is used to update node_cutoff_.
        private var nextNodeCutoff: Int = -1

        // A copy of the current node in the cell tree.
        private var node: AbstractS2CellIndex.CellNode = AbstractS2CellIndex.CellNode()

        // Default constructor; must be followed by a call to Init().
        constructor() {

        }

        // Convenience constructor that calls Init().
        constructor(index: S2CellIndex) {
            init(index)
        }

        // Initializes the iterator.  Should be followed by a call to UnionWith()
        // to visit the contents of each desired leaf cell range.
        fun init(index: S2CellIndex) {
            this.index = index
            clear()
        }

        // Clears all state with respect to which range(s) have been visited.
        fun clear() {
            prevStartId = S2CellId.none
            nodeCutoff = -1
            nextNodeCutoff = -1
            setDone()
        }

        // Positions the ContentsIterator at the first (cell_id, value) pair that
        // covers the given leaf cell range.  Note that when multiple leaf cell
        // ranges are visited using the same ContentsIterator, duplicate values
        // may be suppressed.  If you don't want this behavior, call Clear() first.
        fun startUnion(range: RangeIterator) {
            check(index != null)
            if (range.startId() < prevStartId) {
                nodeCutoff = -1;  // Can't automatically eliminate duplicates.
            }
            prevStartId = range.startId()

            // TODO(ericv): Since RangeNode only uses 12 of its 16 bytes, we could add a
            // "label" field without using any extra space.  Then we could store a leaf
            // node of cell_tree_ directly in each RangeNode, where the cell_id is
            // implicitly defined as the one that covers the current leaf cell range.
            // This would save quite a bit of space; e.g. if the given cells are
            // non-overlapping, then cell_tree_ would be empty (since every node is a
            // leaf node and could therefore be stored directly in a RangeNode).  It
            // would also be faster because cell_tree_ would rarely be accessed.
            val contents = range.contents()
            if (contents <= nodeCutoff) {
                setDone()
            } else {
                node = index!!.cellTree[contents]
            }

            // When visiting ancestors, we can stop as soon as the node index is smaller
            // than any previously visited node index.  Because indexes are assigned
            // using a preorder traversal, such nodes are guaranteed to have already
            // been reported.
            nextNodeCutoff = contents

        }

        // The S2CellId of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun cellId(): S2CellId {
            check(!done())
            return node.cellId
        }

        // The value of the current (cell_id, value) pair.
        // REQUIRES: !done()
        fun value(): Label {
            check(!done())
            return node.label
        }

        // Returns the current (cell_id, value) pair.
        // REQUIRES: !done()
        fun valuedCell(): ValuedCell<Label> {
            check(!done())
            return ValuedCell(node.cellId, node.label)
        }

        // Returns true if all (cell_id, value) pairs have been visited.
        fun done(): Boolean = node.label == kDoneContents

        // Advances the iterator to the next (cell_id, value) pair covered by the
        // current leaf cell range.
        // REQUIRES: !done()
        fun next() {
            check(index != null)
            check(!done())
            if (node.parent <= nodeCutoff) {
                // We have already processed this node and its ancestors.
                nodeCutoff = nextNodeCutoff
                setDone()
            } else {
                node = index!!.cellTree[node.parent].copy()
            }

        }

        // node_.label == kDoneContents indicates that done() is true.
        private fun setDone() {
            node = AbstractS2CellIndex.CellNode()
        }

    }

    companion object {

        private val logger = KotlinLogging.logger(S2CellIndex::class.java.name)

        const val kDoneContents: Label = -1;

    }

}
