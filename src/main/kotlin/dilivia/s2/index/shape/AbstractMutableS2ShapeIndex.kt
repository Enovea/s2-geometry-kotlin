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

import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkLE
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireEQ
import dilivia.math.R1Interval
import dilivia.math.R2Rect
import dilivia.math.vectors.R2Point
import dilivia.s2.S2CellId
import dilivia.s2.S2PaddedCell
import dilivia.s2.S2Point
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.index.shape.InitialPosition.BEGIN
import dilivia.s2.index.shape.InitialPosition.UNPOSITIONED
import dilivia.s2.region.S2CellUnion
import dilivia.s2.shape.Edge
import dilivia.s2.shape.InteriorTracker
import dilivia.s2.shape.S2Shape
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.pow
import java.time.Duration
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicReference
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import kotlin.concurrent.withLock

interface S2ShapeIndexState {
    fun pendingAdditionsBegin(): Int
    fun pendingAdditionsBegin(value: Int)
    fun pendingRemovals(): List<RemovedShape>
    fun addRemovedShape(shape: RemovedShape)
    fun clearRemovedShapes()
    fun status(status: S2ShapeIndexStatus)
    fun status(): S2ShapeIndexStatus
    fun clear() {
        pendingAdditionsBegin(0)
        clearRemovedShapes()
        status(S2ShapeIndexStatus.FRESH)
    }
}

class S2ShapeIndexMapState(
    private val stateMap: MutableMap<String, DefaultS2ShapeIndexState>,
    indexName: String,
) : S2ShapeIndexState {
    companion object {
        val stateLocks = ConcurrentHashMap<String, Lock>()
    }

    private val indexName: String = indexName.intern()

    init {
        stateLocks.computeIfAbsent(indexName) { ReentrantLock(true) }
    }

    private fun getIndexState(): DefaultS2ShapeIndexState {
        val state: DefaultS2ShapeIndexState
        if (!stateMap.containsKey(indexName)) {
            state = DefaultS2ShapeIndexState()
            stateMap[indexName] = state
        } else {
            state = stateMap.getValue(indexName)
        }
        return state
    }

    override fun pendingAdditionsBegin(): Int {
        val state: S2ShapeIndexState = getIndexState()
        return state.pendingAdditionsBegin()
    }

    override fun pendingAdditionsBegin(value: Int): Unit {
        val state = getIndexState()
        state.pendingAdditionsBegin(value)
        stateMap[indexName] = state
    }

    override fun pendingRemovals(): List<RemovedShape> {
        val state = getIndexState()
        return state.pendingRemovals()
    }

    override fun addRemovedShape(shape: RemovedShape): Unit {
        val state = getIndexState()
        state.addRemovedShape(shape)
        stateMap[indexName] = state
    }

    override fun clearRemovedShapes(): Unit {
        val state = getIndexState()
        state.clearRemovedShapes()
        stateMap[indexName] = state
    }

    override fun status(): S2ShapeIndexStatus = stateLocks.getValue(indexName).withLock {
        val state = getIndexState()
        state.status()
    }

    override fun status(status: S2ShapeIndexStatus) = stateLocks.getValue(indexName).withLock {
        val state = getIndexState()
        state.status(status)
        stateMap[indexName] = state
    }

    override fun clear(): Unit = synchronized(indexName) {
        val state = getIndexState()
        state.clear()
        stateMap[indexName] = state
    }

}

data class DefaultS2ShapeIndexState(
    private var pendingAdditionsBegin: Int = 0,
    private val pendingRemovals: MutableList<RemovedShape> = mutableListOf(),
    private var indexStatus: AtomicReference<S2ShapeIndexStatus> = AtomicReference(S2ShapeIndexStatus.FRESH)
) : S2ShapeIndexState {
    override fun pendingAdditionsBegin(): Int = pendingAdditionsBegin

    override fun pendingAdditionsBegin(value: Int) {
        pendingAdditionsBegin = value
    }

    override fun pendingRemovals(): List<RemovedShape> = pendingRemovals

    override fun addRemovedShape(shape: RemovedShape) {
        pendingRemovals.add(shape)
    }

    override fun clearRemovedShapes() {
        pendingRemovals.clear()
    }

    override fun status(status: S2ShapeIndexStatus) {
        indexStatus.set(status)
    }

    override fun status(): S2ShapeIndexStatus = indexStatus.get()

}


/**
 * @property options The index options.
 * @constructor Create a MutableS2ShapeIndex with the given options. Option values may be changed by calling Init().
 * @param options The options supplied for this index.
 *
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
abstract class AbstractMutableS2ShapeIndex(protected var options: Options = Options()) : S2ShapeIndex() {

    private var updateState: UpdateState? = null

    //    /** The id of the first shape that has been queued for addition but not processed yet. */
//    protected abstract var pendingAdditionsBegin: Int
//
//    /** The set of shapes that have been queued for removal but not processed yet. */
//    private val pendingRemovals = mutableListOf<RemovedShape>()
//
//    /** The current index status. Reads and writes to this field are guarded by "lock". */
//    protected abstract var indexStatus: S2ShapeIndexStatus
    protected abstract val state: S2ShapeIndexState
    val pendingAdditionsBegin: Int
        get() = state.pendingAdditionsBegin()
    val pendingRemovals: List<RemovedShape>
        get() = state.pendingRemovals()
    val status: S2ShapeIndexStatus
        get() = state.status()

    /**
     * Additions and removals are queued and processed on the first subsequent query. There are several reasons to do
     * this:
     *
     * - It is significantly more efficient to process updates in batches.
     * - Often the index will never be queried, in which case we can save both the time and memory required to build it.
     *   Examples:
     *   + S2Loops that are created simply to pass to an S2Polygon. (We don't need the S2Loop index, because S2Polygon
     *     builds its own index.)
     *   + Applications that load a database of geometry and then query only a small fraction of it.
     *   + Applications that only read and write geometry (Decode/Encode).
     *
     * The main drawback is that we need to go to some extra work to ensure that methods are still thread-safe. Note
     * that the goal is *not* to make this class thread-safe in general, but simply to hide the fact that we defer some
     * of the indexing work until query time.
     */
    private val lock: ReentrantLock = ReentrantLock()

    /**
     * Options that affect construction of the MutableS2ShapeIndex.
     *
     * @property maxEdgesPerCell The maximum number of edges per cell. If a cell has more than this many edges that are
     * not considered "long" relative to the cell size, then it is subdivided.
     *
     * Values between 10 and 50 represent a reasonable balance between memory usage, construction time, and query time.
     * Small values make queries faster, while large values make construction faster and use less memory.
     * Values higher than 50 do not save significant additional memory, and query times can increase substantially,
     * especially for algorithms that visit all pairs of potentially intersecting edges (such as polygon validation),
     * since this is quadratic in the number of edges per cell.
     *
     * Note that the *average* number of edges per cell is generally slightly less than half of the maximum value
     * defined here.
     *
     * Defaults = 10
     */
    data class Options(val maxEdgesPerCell: Int = kIndexDefaultMaxEdgesPerCell)

    fun options(): Options = options

    /**
     * Index cell iterator. Used by S2ShapeIndex class to provided a S2ShapeIndexCellIterator.
     *
     * @constructor  Default constructor; must be followed by a call to init().
     */
    class Iterator() : IteratorBase() {

        /** The iterate index. */
        private var index: AbstractMutableS2ShapeIndex? = null

        /** Current iterator position. S2CellId.sentinel() means that the the iterator is at the end of the index. */
        private var iter: S2CellId = S2CellId.sentinel

        /**
         * Constructs an iterator positioned as specified. By default dilivia.iterators are unpositioned, since this avoids an
         * extra seek in this situation where one of the seek methods (such as Locate) is immediately called.
         *
         * If you want to position the iterator at the beginning, e.g. in order to loop through the entire index,
         * do this:
         *
         * <pre>
         *     val iter = MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.BEGIN)
         *     while (!iter.done()) {
         *          ...
         *         iter.next()
         *     }
         * </pre>
         *
         * @param index The shape index to iterate.
         * @param pos The inital position of the iterator. (default = UNPOSITIONED)
         */
        constructor(index: AbstractMutableS2ShapeIndex, pos: InitialPosition = UNPOSITIONED) : this() {
            init(index, pos)
        }

        /**
         * Initializes an iterator for the given MutableS2ShapeIndex. This method may also be called in order to
         * restore an iterator to a valid state after the underlying index has been updated (although it is usually
         * easier just to declare a new iterator whenever required, since iterator construction is cheap).
         *
         * @param index The shape index to iterate.
         * @param pos The inital position of the iterator. (default = UNPOSITIONED)
         */
        fun init(index: AbstractMutableS2ShapeIndex, pos: InitialPosition) {
            index.maybeApplyUpdates()
            initStale(index, pos)
        }

        /**
         * Initialize an iterator for the given MutableS2ShapeIndex without applying any pending updates.
         * This can be used to observe the actual current state of the index without modifying it in any way.
         *
         * @param index The shape index to iterate.
         * @param pos The initial position of the iterator. (default = UNPOSITIONED)
         */
        fun initStale(index: AbstractMutableS2ShapeIndex, pos: InitialPosition = UNPOSITIONED) {
            this.index = index
            this.iter = when {
                index.isEmpty() -> S2CellId.sentinel
                pos == BEGIN -> index.firstCellId()!!
                else -> S2CellId.sentinel
            }
            refresh()
        }

        override fun cell(): S2ShapeIndexCell {
            // Since MutableS2ShapeIndex always sets the "cell" field, we can skip the logic in the base class that
            // conditionally calls getCell()
            val cell = rawCell()
            check(cell != null)
            return cell
        }

        override fun begin() {
            // Make sure that the iterator is initialized
            val index = this.index ?: throw IllegalStateException("The iterator must be initialized before calling refresh.")
            // Make sure that the index has not been modified since Init() was called.
            check(index.isFresh())
            // Position the iterator on the first cell id and refresh.
            iter = if (index.isEmpty()) S2CellId.sentinel else index.firstCellId()!!
            refresh()
        }

        override fun finish() {
            iter = S2CellId.sentinel
            refresh()
        }

        override fun next() {
            checkState { !done() }
            // Make sure that the iterator is initialized
            val index = this.index ?: throw IllegalStateException("The iterator must be initialized before calling refresh.")
            // Position the iterator on the next cell. If we are on the last cell, position the iterator on the sentinel
            // to indicate that we are at the end.
            iter = index.higherCellId(iter) ?: S2CellId.sentinel
            refresh()
        }

        override fun prev(): Boolean {
            // Make sure that the iterator is initialized
            val index = this.index ?: throw IllegalStateException("The iterator must be initialized before calling refresh.")
            // If the index is empty or if we are on the begin of the index, do nothing and return false
            if (index.isEmpty() || iter == index.firstCellId()) return false
            // Position the iterator on the previous cell.
            iter = index.lowerCellId(iter)
                ?: throw IllegalStateException("The iterator is already on the beginning of the index. (iter = $iter, first=${index.firstCellId()})")
            refresh()
            return true
        }

        override fun seek(target: S2CellId) {
            // Make sure that the iterator is initialized
            val index = this.index ?: throw IllegalStateException("The iterator must be initialized before calling refresh.")

            iter = index.ceilingCellId(target) ?: S2CellId.sentinel
            refresh()
        }

        override fun clone(): Iterator {
            val clone = super.clone() as Iterator
            clone.index = index
            clone.iter = this.iter
            return clone
        }

        /**
         * Updates the iterator base fields.
         */
        private fun refresh() {
            val index = this.index
                ?: throw IllegalStateException("The iterator must be initialized before calling refresh.")
            if (iter == S2CellId.sentinel) {
                setFinished()
            } else {
                setState(iter, index.getShapeIndexCell(iter))
            }
        }

        override fun getCell(): S2ShapeIndexCell? {
            throw IllegalStateException("Should never be called.")
        }

        companion object {
            private val logger = KotlinLogging.logger(Iterator::class.java.name)
        }

    }

    abstract fun isEmpty(): Boolean

    abstract fun firstCellId(): S2CellId?

    abstract fun higherCellId(id: S2CellId): S2CellId?

    abstract fun lowerCellId(id: S2CellId): S2CellId?

    abstract fun ceilingCellId(id: S2CellId): S2CellId?

    /**
     * Adds a shape to the index. Also assigns a unique id to the shape (shape.id) and returns that id. Shape ids are
     * assigned sequentially starting from 0 in the order shapes are added.
     * Invalidates all dilivia.iterators and their associated data by setting the index status to STALE.
     *
     * @param shape The shape to add.
     * @return The id affected to the shape.
     */
    @Synchronized
    fun add(shape: S2Shape): Int {
        // Additions are processed lazily by ApplyUpdates().
        addShapeToStore(shape)
        state.status(S2ShapeIndexStatus.STALE)
        return shape.id
    }

    protected abstract fun addShapeToStore(shape: S2Shape)

    /**
     * Removes the given shape from the index and return the shape.
     * Invalidates all dilivia.iterators and their associated data by setting the index status to STALE.
     *
     * @param shapeId The identifier of the shape to remove.
     * @return The removed shape.
     */
    @Synchronized
    fun remove(shapeId: Int): S2Shape {
        // This class updates itself lazily, because it is much more efficient to process additions and removals in
        // batches.
        val shape = removeShapeFromStore(shapeId)
        check(shape != null)
        if (shapeId >= state.pendingAdditionsBegin()) {
            // We are removing a shape that has not yet been added to the index,
            // so there is nothing else to do.
        } else {
            // We build the new RemovedShape
            val numEdges = shape.numEdges
            val removedEdges = mutableListOf<Edge>()
            for (e in 0 until numEdges) {
                removedEdges.add(shape.edge(e))
            }
            state.addRemovedShape(
                RemovedShape(
                    shapeId = shape.id,
                    hasInterior = shape.dimension == 2,
                    containsTrackerOrigin = S2Shape.containsBruteForce(shape, kInteriorTrackerOrigin),
                    edges = removedEdges
                )
            )
        }
        state.status(S2ShapeIndexStatus.STALE)
        return shape
    }

    protected abstract fun removeShapeFromStore(shapeId: Int): S2Shape?

    /**
     * Resets the index to its original state and returns all shapes to the caller. This method is much more efficient
     * than removing all shapes one at a time.
     *
     * @return The removed shapes.
     */
    @Synchronized
    fun removeAll() {
        check(updateState == null)
        state.clear()
        removeAllShapesFromStore()
    }

    protected abstract fun removeAllShapesFromStore()

    protected abstract fun getShapeIndexCell(id: S2CellId): S2ShapeIndexCell

    protected abstract fun setShapeIndexCell(id: S2CellId, cell: S2ShapeIndexCell)

    protected abstract fun removeEmptyShapeIndexCells()

    protected abstract fun removeShapeIndexCell(id: S2CellId)

    /**
     * Forces the index build. Calls to add() and remove() are normally queued and processed on the first subsequent
     * query (in a thread-safe way). This has many advantages, the most important of which is that sometimes there *is*
     * no subsequent query, which lets us avoid building the index completely.
     *
     * This method forces any pending updates to be applied immediately. Calling this method is rarely a good idea.
     * (One valid reason is to exclude the cost of building the index from benchmark results.)
     */
    @Synchronized
    fun forceBuild() {
        logger.debug { "forceBuild | status = ${state.status()}" }
        // No locks required because this is not a const method. It is the client's responsibility to ensure correct
        // thread synchronization.
        if (state.status() != S2ShapeIndexStatus.FRESH) {
            applyUpdatesInternal()
            state.status(S2ShapeIndexStatus.FRESH)
            logger.trace {
                """
                |
                |==============================================
                | Shape index updated
                |----------------------------------------------
                |${toDebugString()}
                |----------------------------------------------
                | """.trimMargin()
            }
        }
    }

    /**
     * Indicates if there are no pending updates that need to be applied. This can be useful to avoid building the
     * index unnecessarily, or for choosing between two different algorithms depending on whether the index is
     * available.
     *
     * The returned index status may be slightly out of date if the index was built in a different thread. This is fine
     * for the intended use (as an efficiency hint), but it should not be used by internal methods
     * (see MaybeApplyUpdates).
     *
     * @return true if the index is fresh and false otherwise.
     */
    fun isFresh(): Boolean = state.status() == S2ShapeIndexStatus.FRESH

    /////////////////////////////////////////////////////////////////////
    //                       Internal Methods                          //
    /////////////////////////////////////////////////////////////////////

    /**
     * Constructs an iterator positioned as specified. By default dilivia.iterators are unpositioned, since this avoids an
     * extra seek in this situation where one of the seek methods (such as Locate) is immediately called.
     *
     * If you want to position the iterator at the beginning, e.g. in order to loop through the entire index, do this:
     *
     * <pre>
     *     val iter = MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.BEGIN)
     *     while (!iter.done()) {
     *          ...
     *         iter.next()
     *     }
     * </pre>
     *
     * @param pos The initial position
     * @return A new iterator instance.
     */
    override fun newIterator(pos: InitialPosition): IteratorBase = Iterator(this, pos)

    /**
     * Indicates if this is the first update to the index.
     *
     * @return true if it is the first update and false otherwise.
     */
    private fun isFirstUpdate(): Boolean {
        // Note that it is not sufficient to check whether cellMap is empty, since entries are added during the update
        // process.
        return state.pendingAdditionsBegin() == 0
    }

    /**
     * Given that the given shape is being updated, indicates if it is being removed (as opposed to being added).
     * NOTE: This method does not check if the shape is being updated.
     *
     * @param shapeId The shape identifier to check.
     * @return true if the shape is being removed.
     */
    private fun isShapeBeingRemoved(shapeId: Int): Boolean {
        // All shape ids being removed are less than all shape ids being added.
        return shapeId < state.pendingAdditionsBegin()
    }

    /**
     * Ensure that any pending updates have been applied. This method must be called before accessing the cellMap field,
     * even if the indexStatus appears to be FRESH, because a memory barrier is required in order to ensure that all
     * the index updates are visible if the updates were done in another thread.
     */
    protected fun maybeApplyUpdates() {
        // To avoid acquiring and releasing the spinlock on every query, we use
        // atomic operations when testing whether the status is FRESH and when
        // updating the status to be FRESH.  This guarantees that any thread that
        // sees a status of FRESH will also see the corresponding index updates.
        if (state.status() != S2ShapeIndexStatus.FRESH) {
            applyUpdatesThreadSafe()
        }
    }

    /** Apply any pending updates in a thread-safe way. */
    private fun applyUpdatesThreadSafe() {
        lock.lock()
        when (state.status()) {
            S2ShapeIndexStatus.FRESH -> lock.unlock()
            S2ShapeIndexStatus.UPDATING -> {
                val updateState = updateState
                check(updateState != null)
                // Wait until the updating thread is finished. We do this by attempting to lock a mutex that is held by
                // the updating thread.  When this mutex is unlocked the indexStatus is guaranteed to be FRESH.
                ++updateState.numWaiting
                lock.unlock()
                updateState.waitMutex.lock()
                lock.lock()
                --updateState.numWaiting
                unlockAndSignal()  // Notify other waiting threads.
            }
            else -> {
                checkEQ(S2ShapeIndexStatus.STALE, state.status())
                state.status(S2ShapeIndexStatus.UPDATING)
                // Allocate the extra state needed for thread synchronization. We keep the spinlock held while doing
                // this, because
                // (1) memory allocation is fast, so the chance of a context switch while holding the lock is low
                // (2) by far the most common situation is that there is no contention, and this saves an extra lock and
                //     unlock step;
                // (3) even in the rare case where there is contention, the main side effect is that some other thread
                //     will burn a few CPU cycles rather than sleeping.
                val updateState = UpdateState()
                this.updateState = updateState
                // lock.Lock wait_mutex *before* calling Unlock() to ensure that all other threads will block on it.
                updateState.waitMutex.lock()
                // Release the spinlock before doing any real work.
                lock.unlock()
                applyUpdatesInternal()
                lock.lock()
                // indexStatus can be updated to FRESH only while locked *and* using an atomic store operation, so that
                // maybeApplyUpdates() can check whether the index is FRESH without acquiring the spinlock.
                state.status(S2ShapeIndexStatus.FRESH)

                logger.trace {
                    """
                    |
                    |==============================================
                    | Shape index updated
                    |----------------------------------------------
                    |${toDebugString()}
                    |----------------------------------------------
                    | """.trimMargin()
                }

                unlockAndSignal()  // Notify any waiting threads.
            }
        }
    }

    /**
     * Releases lock and wakes up any waiting threads by releasing waitMutex. If this was the last waiting thread,
     * also deletes updateState.
     * REQUIRES: lock is held.
     * REQUIRES: updateState.waitMutex is held.
     */
    private fun unlockAndSignal() {
        checkEQ(S2ShapeIndexStatus.FRESH, state.status())
        val updateState = updateState
        check(updateState != null)
        val numWaiting = updateState.numWaiting
        lock.unlock()
        // Allow another waiting thread to proceed. Note that no new threads can start waiting because the indexStatus
        // is now FRESH, and the caller is required to prevent any new mutations from occurring while these
        // methods are running.
        //
        // We need to unlock waitMutex before destroying it even if there are no waiting threads.
        updateState.waitMutex.unlock()
        if (numWaiting == 0) {
            this.updateState = null
        }
    }

    /**
     * This method updates the index by applying all pending additions and removals. It does *not* update indexStatus
     * (see ApplyUpdatesThreadSafe).
     */
    private fun applyUpdatesInternal() {
        val startTime = System.nanoTime()
        // Check whether we have so many edges to process that we should process them in multiple batches to save
        // memory. Building the index can use up to 20x as much memory (per edge) as the final index size.
        val removals = state.pendingRemovals()
        val removalsCount = removals.size
        val batches = getUpdateBatches(state.pendingAdditionsBegin(), removals)
        val shapeAdditions = if (batches.isEmpty()) 0 else (batches.last().additionsEnd - state.pendingAdditionsBegin())

        logger.trace {
            """Update batches:
            |--------------------------------------
            |${batches.joinToString("\n")}
            |--------------------------------------
        """.trimMargin()
        }

        var numEdges: Long = 0L
        batches.forEachIndexed { index, batch ->
            if (batches.size > 1) {
                print("\rProcess batch ${index + 1} / ${batches.size}")
            }
            val allEdges: Array<MutableList<FaceEdge>> = Array(6) { mutableListOf() }
            val tracker = InteriorTracker()
            if (removals.isNotEmpty()) {
                logger.trace { "Start removing shapes: ${removals.joinToString { s -> s.shapeId.toString() }}" }
                // The first batch implicitly includes all shapes being removed.
                for (pendingRemoval in removals) {
                    removeShape(pendingRemoval, allEdges, tracker)
                }
                state.clearRemovedShapes()
            }
            for (id in state.pendingAdditionsBegin() until batch.additionsEnd) {
                addShape(id, allEdges, tracker)
            }
            numEdges += allEdges.map { edges -> edges.size }.sum()
            for (face in 0..5) {
                updateFaceEdges(face, allEdges[face], tracker)
                // Save memory by clearing vectors after we are done with them.
                allEdges[face].clear()
            }
            state.pendingAdditionsBegin(batch.additionsEnd)
        }

        removeEmptyShapeIndexCells()

        // It is the caller's responsibility to update index_status_.
        val duration = Duration.ofNanos(System.nanoTime() - startTime)
        if (duration.toMillis() >= 1000L) {
            logger.info { "Index update took $duration for ${batches.size} batched, $shapeAdditions shapes added, $removalsCount shapes removed, $numEdges face edges" }
        }
    }

    // Count the number of edges being updated, and break them into several
    // batches if necessary to reduce the amount of memory needed.  (See the
    // documentation for FLAGS_s2shape_index_tmp_memory_budget_mb.)
    private fun getUpdateBatches(pendingAdditionsBegin: Int, pendingRemovals: List<RemovedShape>): List<BatchDescriptor> {
        val batches = mutableListOf<BatchDescriptor>()
        // Count the edges being removed and added.
        var numEdgesRemoved = 0
        if (pendingRemovals.isNotEmpty()) {
            for (pending_removal in pendingRemovals) {
                numEdgesRemoved += pending_removal.edges.size
            }
        }
        var numEdgesAdded = 0
        for (id in pendingAdditionsBegin until nextNewShapeId()) {
            val shape = shape(id) ?: continue
            numEdgesAdded += shape.numEdges
        }
        var numEdges = numEdgesRemoved + numEdgesAdded

//        if (numEdges * kTmpBytesPerEdge <= kTmpMemoryBudgetBytes) {
        if (nextNewShapeId() > 1) {
            logger.debug { "Num edges = $numEdges" }
        }
        if (numEdges <= 400_000) {
            // We can update all edges at once without exceeding kTmpMemoryBudgetBytes.
            batches.add(BatchDescriptor(nextNewShapeId(), numEdges))
            return batches
        }
        // Otherwise, break the updates into up to several batches, where the size
        // of each batch is chosen so that all batches use approximately the same
        // high-water memory.  GetBatchSizes() returns the recommended number of
        // edges in each batch.
        val batchSizes = mutableListOf<Int>()
        var remainingsEdgeCount = numEdges
        while (remainingsEdgeCount > 0) {
            val size = min(remainingsEdgeCount, 400_000)
            batchSizes.add(size)
            remainingsEdgeCount -= size
        }
//        val batchSizes = getBatchSizes(
//                numEdges,
//                kMaxUpdateBatches,
//                kFinalBytesPerEdge.toDouble(),
//                kTmpBytesPerEdge.toDouble(),
//                kTmpMemoryBudgetBytes.toDouble()
//        )

        // We always process removed edges in a single batch, since (1) they already
        // take up a lot of memory because we have copied all their edges, and (2)
        // AbsorbIndexCell() uses (shapes_[id] == nullptr) to detect when a shape is
        // being removed, so in order to split the removals into batches we would
        // need a different approach (e.g., temporarily add fake entries to shapes_
        // and restore them back to nullptr as shapes are actually removed).
        numEdges = 0
        if (pendingRemovals.isNotEmpty()) {
            numEdges += numEdgesRemoved
            if (numEdges >= batchSizes[0]) {
                batches.add(BatchDescriptor(pendingAdditionsBegin, numEdges))
                numEdges = 0
            }
        }
        // Keep adding shapes to each batch until the recommended number of edges
        // for that batch is reached, then move on to the next batch.
        for (id in pendingAdditionsBegin until nextNewShapeId()) {
            val shape = shape(id) ?: continue
            numEdges += shape.numEdges
            if (numEdges >= batchSizes[batches.size]) {
                batches.add(BatchDescriptor(id + 1, numEdges))
                numEdges = 0
            }
        }
        // Some shapes have no edges.  If a shape with no edges is the last shape to
        // be added or removed, then the final batch may not include it, so we fix
        // that problem here.
        batches.last().additionsEnd = nextNewShapeId()
        checkLE(batches.size, kMaxUpdateBatches)
        return batches
    }

    // Given "num_items" items, each of which uses "tmp_bytes_per_item" while it
    // is being updated but only "final_bytes_per_item" in the end, divide the
    // items into batches that have approximately the same *total* memory usage
    // consisting of the temporary memory needed for the items in the current
    // batch plus the final size of all the items that have already been
    // processed.  Use the fewest number of batches (but never more than
    // "max_batches") such that the total memory usage does not exceed the
    // combined final size of all the items plus "tmp_memory_budget_bytes".
    /* static */
    private fun getBatchSizes(
        numItems: Int,
        maxBatches: Int,
        finalBytesPerItem: Double,
        tmpBytesPerItem: Double,
        tmpMemoryBudgetBytes: Double
    ): List<Int> {
        val batchSizes = mutableListOf<Int>()
        var remainingItems = numItems
        // This code tries to fit all the data into the same memory space
        // ("total_budget_bytes") at every iteration.  The data consists of some
        // number of processed items (at "final_bytes_per_item" each), plus some
        // number being updated (at "tmp_bytes_per_item" each).  The space occupied
        // by the items being updated is the "free space".  At each iteration, the
        // free space is multiplied by (1 - final_bytes_per_item/tmp_bytes_per_item)
        // as the items are converted into their final form.
        val finalBytes = remainingItems * finalBytesPerItem
        val finalBytesRatio = finalBytesPerItem / tmpBytesPerItem
        val freeSpaceMultiplier = 1 - finalBytesRatio

        // The total memory budget is the greater of the final size plus the allowed
        // temporary memory, or the minimum amount of memory required to limit the
        // number of batches to "max_batches".
        val totalBudgetBytes = max(finalBytes + tmpMemoryBudgetBytes, finalBytes / (1 - pow(freeSpaceMultiplier, maxBatches.toDouble())))

        // "max_batch_items" is the number of items in the current batch.
        var maxBatchItems = totalBudgetBytes / tmpBytesPerItem
        var i = 0
        while (i + 1 < maxBatches && remainingItems > 0) {
            val batchItems = min(remainingItems, (maxBatchItems + 1).toInt())
            batchSizes.add(batchItems)
            remainingItems -= batchItems
            maxBatchItems *= freeSpaceMultiplier
            ++i
        }
        checkLE(batchSizes.size, maxBatches)
        return batchSizes
    }

    // Clip all edges of the given shape to the six cube faces, add the clipped
    // edges to "all_edges", and start tracking its interior if necessary.
    private fun addShape(id: Int, all_edges: Array<MutableList<FaceEdge>>, tracker: InteriorTracker) {
        val shape = shape(id) ?: return  // This shape has already been removed.
        // Construct a template for the edges to be added.
        val faceEdge = FaceEdge(
            shapeId = id,
            hasInterior = shape.dimension == 2
        )
        if (faceEdge.hasInterior) {
            val containsFocus = S2Shape.containsBruteForce(shape, tracker.focus())
            tracker.addShape(id, containsFocus)
        }
        val numEdges = shape.numEdges
        for (e in 0 until numEdges) {
            val edge = shape.edge(e)
            addFaceEdge(faceEdge.copy(edgeId = e, edge = edge, maxLevel = getEdgeMaxLevel(edge)), all_edges)
        }

//        logger.trace {
//            """Add shape $id
//            |---------------------------------------
//            |Interior tracker: $tracker
//            |${
//                all_edges
//                        .mapIndexed { index, faceEdgeList ->
//                            "$index: \n${
//                                faceEdgeList.joinToString(",\n") { faceEdge ->
//                                    "   - (sid=${faceEdge.shapeId}, eid=${faceEdge.edgeId}, edge=${faceEdge.edge})"
//                                }
//                            }"
//                        }.joinToString("")
//            }
//            |---------------------------------------
//        """.trimMargin()
//        }
    }

    private fun removeShape(removed: RemovedShape, allEdges: Array<MutableList<FaceEdge>>, tracker: InteriorTracker) {
        val faceEdge = FaceEdge(
            shapeId = removed.shapeId,
            edgeId = -1,  // Not used or needed for removed edges.
            hasInterior = removed.hasInterior,
        )
        if (faceEdge.hasInterior) {
            tracker.addShape(faceEdge.shapeId, removed.containsTrackerOrigin)
        }
        for (removedEdge in removed.edges) {
            addFaceEdge(faceEdge.copy(edge = removedEdge, maxLevel = getEdgeMaxLevel(removedEdge)), allEdges)
        }
    }

    private fun addFaceEdge(faceEdge: FaceEdge, allEdges: Array<MutableList<FaceEdge>>) {
        var edge = faceEdge.copy()
        // Fast path: both endpoints are on the same face, and are far enough from
        // the edge of the face that don't intersect any (padded) adjacent face.
        val aFace = S2Coords.face(edge.edge.v0)
        if (aFace == S2Coords.face(edge.edge.v1)) {
            edge.a = S2Coords.validFaceXyzToUv(aFace, edge.edge.v0)
            edge.b = S2Coords.validFaceXyzToUv(aFace, edge.edge.v1)
            val kMaxUV = 1 - kCellPadding
            if (abs(edge.a[0]) <= kMaxUV && abs(edge.a[1]) <= kMaxUV && abs(edge.b[0]) <= kMaxUV && abs(edge.b[1]) <= kMaxUV) {
                allEdges[aFace].add(edge)
                return
            }
        }
        // Otherwise we simply clip the edge to all six faces.
        for (face in 0..5) {
            val abClippedUV = S2EdgeClipping.clipToPaddedFace(edge.edge.v0, edge.edge.v1, face, kCellPadding)
            if (abClippedUV != null) {
                edge = faceEdge.copy(a = abClippedUV.first, b = abClippedUV.second)
                allEdges[face].add(edge)
            }
        }

    }

    // Return the first level at which the edge will *not* contribute towards
    // the decision to subdivide.
    fun getEdgeMaxLevel(edge: Edge): Int {
        // Compute the maximum cell size for which this edge is considered "long".
        // The calculation does not need to be perfectly accurate, so we use Norm()
        // rather than Angle() for speed.
        val cellSize = ((edge.v0 - edge.v1).norm() * kIndexCellSizeToLongEdgeRatio)
        // Now return the first level encountered during subdivision where the
        // average cell size is at most "cell_size".
        return S2Coords.projection.kAvgEdge.getLevelForMaxValue(cellSize)
    }

    // Given a face and a vector of edges that intersect that face, add or remove
    // all the edges from the index.  (An edge is added if shapes_[id] is not
    // nullptr, and removed otherwise.)
    private fun updateFaceEdges(face: Int, faceEdges: List<FaceEdge>, tracker: InteriorTracker) {
        logger.trace { "Update face $face with ${faceEdges.size} edges." }
        val numEdges = faceEdges.size
        if (numEdges == 0 && tracker.shapeIds().isEmpty()) return

        // Create the initial ClippedEdge for each FaceEdge.  Additional clipped
        // edges are created when edges are split between child cells.  We create
        // two arrays, one containing the edge data and another containing pointers
        // to those edges, so that during the recursion we only need to copy
        // pointers in order to propagate an edge to the correct child.
        val clippedEdges = mutableListOf<ClippedEdge>()
        val bound = R2Rect()
        for (e in 0 until numEdges) {
            val clipped = ClippedEdge(
                faceEdge = faceEdges[e],
                bound = R2Rect.fromPointPair(faceEdges[e].a, faceEdges[e].b)
            )
            //clipped.face_edge = &face_edges[e]
            //clipped.bound = R2Rect::FromPointPair(face_edges[e].a, face_edges[e].b)
            clippedEdges.add(clipped)
            bound.addRect(clipped.bound)
        }
        logger.trace { "Updated bound: $bound" }


        // Construct the initial face cell containing all the edges, and then update
        // all the edges in the index recursively.
        //EdgeAllocator alloc
        val faceId = S2CellId.fromFace(face)
        var pcell = S2PaddedCell(faceId, kCellPadding)

        // "disjoint_from_index" means that the current cell being processed (and
        // all its descendants) are not already present in the index.
        val disjointFromIndex = isFirstUpdate()

        if (numEdges > 0) {
            val shrunkId = shrinkToFit(pcell, bound)
            logger.trace { "Shrink id = $shrunkId (pcell = $pcell)" }
            if (shrunkId != pcell.id) {
                // All the edges are contained by some descendant of the face cell.  We
                // can save a lot of work by starting directly with that cell, but if we
                // are in the interior of at least one shape then we need to create
                // index entries for the cells we are skipping over.
                skipCellRange(faceId.rangeMin(), shrunkId.rangeMin(), tracker, disjointFromIndex)
                pcell = S2PaddedCell(shrunkId, kCellPadding)
                updateEdges(pcell, clippedEdges, tracker, disjointFromIndex)
                skipCellRange(shrunkId.rangeMax().next(), faceId.rangeMax().next(), tracker, disjointFromIndex)
                return
            }
        }
        // Otherwise (no edges, or no shrinking is possible), subdivide normally.
        updateEdges(pcell, clippedEdges, tracker, disjointFromIndex)
    }


    private fun shrinkToFit(pcell: S2PaddedCell, bound: R2Rect): S2CellId {
        var shrunkId = pcell.shrinkToFit(bound)
        if (!isFirstUpdate() && shrunkId != pcell.id) {
            // Don't shrink any smaller than the existing index cells, since we need
            // to combine the new edges with those cells.
            // Use InitStale() to avoid applying updated recursively.
            val iter = Iterator()
            iter.initStale(this)
            val r = iter.locate(shrunkId)
            if (r == CellRelation.INDEXED) {
                shrunkId = iter.id()
            }
        }
        return shrunkId
    }

    // Skip over the cells in the given range, creating index cells if we are
    // currently in the interior of at least one shape.
    private fun skipCellRange(
        begin: S2CellId,
        end: S2CellId,
        tracker: InteriorTracker/*, EdgeAllocator* alloc*/,
        disjoint_from_index: Boolean
    ) {
        // If we aren't in the interior of a shape, then skipping over cells is easy.
        if (tracker.shapeIds().isEmpty()) return

        // Otherwise generate the list of cell ids that we need to visit, and create
        // an index entry for each one.
        for (skipped_id in S2CellUnion.fromBeginEnd(begin, end)) {
            updateEdges(S2PaddedCell(skipped_id, kCellPadding), mutableListOf(), tracker, disjoint_from_index)
        }
    }

    // Given a cell and a set of ClippedEdges whose bounding boxes intersect that
    // cell, add or remove all the edges from the index.  Temporary space for
    // edges that need to be subdivided is allocated from the given EdgeAllocator.
    // "disjoint_from_index" is an optimization hint indicating that cell_map_
    // does not contain any entries that overlap the given cell.
    private fun updateEdges(pcell: S2PaddedCell, edges: MutableList<ClippedEdge>, tracker: InteriorTracker, disjointFromIndex: Boolean) {
        logger.trace { "Update edges: cell = ${pcell.id}, ${edges.size} edges, disjoint from index = $disjointFromIndex" }
        var disjointFromIndex = disjointFromIndex
        // Cases where an index cell is not needed should be detected before this.
        checkState { (edges.isNotEmpty() || tracker.shapeIds().isNotEmpty()) }

        // This function is recursive with a maximum recursion depth of 30
        // (S2CellId::kMaxLevel).  Note that using an explicit stack does not seem
        // to be any faster based on profiling.

        // Incremental updates are handled as follows.  All edges being added or
        // removed are combined together in "edges", and all shapes with interiors
        // are tracked using "tracker".  We subdivide recursively as usual until we
        // encounter an existing index cell.  At this point we "absorb" the index
        // cell as follows:
        //
        //   - Edges and shapes that are being removed are deleted from "edges" and
        //     "tracker".
        //   - All remaining edges and shapes from the index cell are added to
        //     "edges" and "tracker".
        //   - Continue subdividing recursively, creating new index cells as needed.
        //   - When the recursion gets back to the cell that was absorbed, we
        //     restore "edges" and "tracker" to their previous state.
        //
        // Note that the only reason that we include removed shapes in the recursive
        // subdivision process is so that we can find all of the index cells that
        // contain those shapes efficiently, without maintaining an explicit list of
        // index cells for each shape (which would be expensive in terms of memory).
        var indexCellAbsorbed = false
        if (!disjointFromIndex) {
            // There may be existing index cells contained inside "pcell".  If we
            // encounter such a cell, we need to combine the edges being updated with
            // the existing cell contents by "absorbing" the cell.
            // Use InitStale() to avoid applying updated recursively.
            val iter = Iterator()
            iter.initStale(this)
            when (val r = iter.locate(pcell.id)) {
                CellRelation.DISJOINT -> disjointFromIndex = true
                CellRelation.INDEXED -> {
                    // Absorb the index cell by transferring its contents to "edges" and
                    // deleting it.  We also start tracking the interior of any new shapes.
                    absorbIndexCell(pcell, iter, edges, tracker)
                    indexCellAbsorbed = true
                    disjointFromIndex = true
                }
                else -> checkEQ(CellRelation.SUBDIVIDED, r)
            }
        }

        // If there are existing index cells below us, then we need to keep
        // subdividing so that we can merge with those cells.  Otherwise,
        // MakeIndexCell checks if the number of edges is small enough, and creates
        // an index cell if possible (returning true when it does so).
        if (!disjointFromIndex || !makeIndexCell(pcell, edges, tracker)) {
            // Reserve space for the edges that will be passed to each child.  This is
            // important since otherwise the running time is dominated by the time
            // required to grow the vectors.  The amount of memory involved is
            // relatively small, so we simply reserve the maximum space for every child.
            val childEdges = Array(2) { Array(2) { mutableListOf<ClippedEdge>() } }  // [i][j]
            val numEdges = edges.size
            //for (i in 0..1) {
            //    for (j in 0..1) {
            //    child_edges[i][j].reserve(num_edges)
            //}
            //}

            // Remember the current size of the EdgeAllocator so that we can free any
            // edges that are allocated during edge splitting.
            //size_t alloc_size = alloc->size()

            // Compute the middle of the padded cell, defined as the rectangle in
            // (u,v)-space that belongs to all four (padded) children.  By comparing
            // against the four boundaries of "middle" we can determine which children
            // each edge needs to be propagated to.
            val middle = pcell.middle()

            // Build up a vector edges to be passed to each child cell.  The (i,j)
            // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
            // Note that the vast majority of edges are propagated to a single child.
            // This case is very fast, consisting of between 2 and 4 floating-point
            // comparisons and copying one pointer.  (ClipVAxis is inline.)
            for (e in 0 until numEdges) {
                val edge = edges[e]
                when {
                    edge.bound[0].hi <= middle[0].lo -> {
                        // Edge is entirely contained in the two left children.
                        clipVAxis(edge, middle[1], childEdges[0])
                    }
                    edge.bound[0].lo >= middle[0].hi -> {
                        // Edge is entirely contained in the two right children.
                        clipVAxis(edge, middle[1], childEdges[1])
                    }
                    edge.bound[1].hi <= middle[1].lo -> {
                        // Edge is entirely contained in the two lower children.
                        childEdges[0][0].add(clipUBound(edge, 1, middle[0].hi))
                        childEdges[1][0].add(clipUBound(edge, 0, middle[0].lo))
                    }
                    edge.bound[1].lo >= middle[1].hi -> {
                        // Edge is entirely contained in the two upper children.
                        childEdges[0][1].add(clipUBound(edge, 1, middle[0].hi))
                        childEdges[1][1].add(clipUBound(edge, 0, middle[0].lo))
                    }
                    else -> {
                        // The edge bound spans all four children.  The edge itself intersects
                        // either three or four (padded) children.
                        val left = clipUBound(edge, 1, middle[0].hi)
                        clipVAxis(left, middle[1], childEdges[0])
                        val right = clipUBound(edge, 0, middle[0].lo)
                        clipVAxis(right, middle[1], childEdges[1])
                    }
                }
            }
            // Free any memory reserved for children that turned out to be empty.  This
            // step is cheap and reduces peak memory usage by about 10% when building
            // large indexes (> 10M edges).
            //for (int i = 0; i < 2; ++i) {
            //    for (int j = 0; j < 2; ++j) {
            //    if (child_edges[i][j].empty()) {
            //        vector<const ClippedEdge*>().swap(child_edges[i][j])
            //    }
            //}
            //}

            // Now recursively update the edges in each child.  We call the children in
            // increasing order of S2CellId so that when the index is first constructed,
            // all insertions into cell_map_ are at the end (which is much faster).
            for (pos in 0..3) {
                val (i, j) = pcell.getChildIJ(pos)
                if (childEdges[i][j].isNotEmpty() || tracker.shapeIds().isNotEmpty()) {
                    updateEdges(S2PaddedCell(pcell, i, j), childEdges[i][j], tracker, disjointFromIndex)
                }
            }
            // Free any temporary edges that were allocated during clipping.
            //alloc->Reset(alloc_size)
        }
        if (indexCellAbsorbed) {
            // Restore the state for any edges being removed that we are tracking.
            tracker.restoreStateBefore(state.pendingAdditionsBegin())
        }
    }

    // Given an edge and an interval "middle" along the v-axis, clip the edge
    // against the boundaries of "middle" and add the edge to the corresponding
    // children.
    /* static */
    private fun clipVAxis(edge: ClippedEdge, middle: R1Interval, child_edges: Array<MutableList<ClippedEdge>>) {
        when {
            edge.bound[1].hi <= middle.lo -> child_edges[0].add(edge) // Edge is entirely contained in the lower child.
            edge.bound[1].lo >= middle.hi -> child_edges[1].add(edge)  // Edge is entirely contained in the upper child.
            else -> {
                // The edge bound spans both children.
                child_edges[0].add(clipVBound(edge, 1, middle.hi))
                child_edges[1].add(clipVBound(edge, 0, middle.lo))
            }
        }
    }

    // Given an edge, clip the given endpoint (lo=0, hi=1) of the u-axis so that
    // it does not extend past the given value.
    /* static */
    private fun clipUBound(edge: ClippedEdge, u_end: Int, u: Double): ClippedEdge {
        // First check whether the edge actually requires any clipping.  (Sometimes
        // this method is called when clipping is not necessary, e.g. when one edge
        // endpoint is in the overlap area between two padded child cells.)
        if (u_end == 0) {
            if (edge.bound[0].lo >= u) return edge
        } else {
            if (edge.bound[0].hi <= u) return edge
        }
        // We interpolate the new v-value from the endpoints of the original edge.
        // This has two advantages: (1) we don't need to store the clipped endpoints
        // at all, just their bounding box; and (2) it avoids the accumulation of
        // roundoff errors due to repeated interpolations.  The result needs to be
        // clamped to ensure that it is in the appropriate range.
        val e = edge.faceEdge
        val v = edge.bound[1].project(S2EdgeClipping.interpolateDouble(u, e.a[0], e.b[0], e.a[1], e.b[1]))

        // Determine which endpoint of the v-axis bound to update.  If the edge
        // slope is positive we update the same endpoint, otherwise we update the
        // opposite endpoint.
        val vEnd = u_end xor (if ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1])) 1 else 0)
        return updateBound(edge, u_end, u, vEnd, v)
    }

    // Given an edge, clip the given endpoint (lo=0, hi=1) of the v-axis so that
    // it does not extend past the given value.
    /* static */
    private fun clipVBound(edge: ClippedEdge, v_end: Int, v: Double): ClippedEdge {
        // See comments in ClipUBound.
        if (v_end == 0) {
            if (edge.bound[1].lo >= v) return edge
        } else {
            if (edge.bound[1].hi <= v) return edge
        }
        val e = edge.faceEdge
        val u = edge.bound[0].project(S2EdgeClipping.interpolateDouble(v, e.a[1], e.b[1], e.a[0], e.b[0]))
        val uEnd = v_end xor (if ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1])) 1 else 0)
        return updateBound(edge, uEnd, u, v_end, v)
    }

    // Given an edge and two bound endpoints that need to be updated, allocate and
    // return a new edge with the updated bound.
    /* static */
    private fun updateBound(edge: ClippedEdge, u_end: Int, u: Double, v_end: Int, v: Double): ClippedEdge {
        val bound = edge.bound.clone()
        bound[0][u_end] = u
        bound[1][v_end] = v
        bound[0][1 - u_end] = bound[0][1 - u_end]
        bound[1][1 - v_end] = bound[1][1 - v_end]
        val clipped = ClippedEdge(faceEdge = edge.faceEdge, bound = bound)
        checkState { !clipped.bound.isEmpty }
        checkState { edge.bound.contains(clipped.bound) }
        return clipped
    }

    // Absorb an index cell by transferring its contents to "edges" and/or
    // "tracker", and then delete this cell from the index.  If "edges" includes
    // any edges that are being removed, this method also updates their
    // InteriorTracker state to correspond to the exit vertex of this cell, and
    // saves the InteriorTracker state by calling SaveAndClearStateBefore().  It
    // is the caller's responsibility to restore this state by calling
    // RestoreStateBefore() when processing of this cell is finished.
    private fun absorbIndexCell(pcell: S2PaddedCell, iter: Iterator, edges: MutableList<ClippedEdge>, tracker: InteriorTracker) {
        requireEQ(pcell.id, iter.id())

        // When we absorb a cell, we erase all the edges that are being removed.
        // However when we are finished with this cell, we want to restore the state
        // of those edges (since that is how we find all the index cells that need
        // to be updated).  The edges themselves are restored automatically when
        // UpdateEdges returns from its recursive call, but the InteriorTracker
        // state needs to be restored explicitly.
        //
        // Here we first update the InteriorTracker state for removed edges to
        // correspond to the exit vertex of this cell, and then save the
        // InteriorTracker state.  This state will be restored by UpdateEdges when
        // it is finished processing the contents of this cell.
        if (tracker.isActive() && edges.isNotEmpty() && isShapeBeingRemoved(edges[0].faceEdge.shapeId)) {
            // We probably need to update the InteriorTracker.  ("Probably" because
            // it's possible that all shapes being removed do not have interiors.)
            if (!tracker.atCellId(pcell.id)) {
                tracker.moveTo(pcell.getEntryVertex())
            }
            tracker.drawTo(pcell.getExitVertex())
            tracker.setNextCellId(pcell.id.next())
            for (edge in edges) {
                val faceEdge = edge.faceEdge
                if (!isShapeBeingRemoved(faceEdge.shapeId)) {
                    break  // All shapes being removed come first.
                }
                if (faceEdge.hasInterior) {
                    tracker.testEdge(faceEdge.shapeId, faceEdge.edge)
                }
            }
        }
        // Save the state of the edges being removed, so that it can be restored
        // when we are finished processing this cell and its children.  We don't
        // need to save the state of the edges being added because they aren't being
        // removed from "edges" and will therefore be updated normally as we visit
        // this cell and its children.
        tracker.saveAndClearStateBefore(state.pendingAdditionsBegin())

        // Create a FaceEdge for each edge in this cell that isn't being removed.
        val faceEdges = mutableListOf<FaceEdge>()
        var trackerMoved = false
        val cell = iter.cell()
        for (s in 0 until cell.numClipped) {
            val clipped = cell.clipped(s)
            val shapeId = clipped.shapeId
            val shape = shape(shapeId) ?: continue // This shape is being removed.
            val numEdges = clipped.numEdges

            // If this shape has an interior, start tracking whether we are inside the
            // shape.  UpdateEdges() wants to know whether the entry vertex of this
            // cell is inside the shape, but we only know whether the center of the
            // cell is inside the shape, so we need to test all the edges against the
            // line segment from the cell center to the entry vertex.
            val edge = FaceEdge(
                shapeId = shape.id,
                hasInterior = (shape.dimension == 2)
            )
            if (edge.hasInterior) {
                tracker.addShape(shapeId, clipped.containsCenter)
                // There might not be any edges in this entire cell (i.e., it might be
                // in the interior of all shapes), so we delay updating the tracker
                // until we see the first edge.
                if (!trackerMoved && numEdges > 0) {
                    tracker.moveTo(pcell.getCenter())
                    tracker.drawTo(pcell.getEntryVertex())
                    tracker.setNextCellId(pcell.id)
                    trackerMoved = true
                }
            }
            for (i in 0 until numEdges) {
                val e = clipped.edge(i)
                edge.edgeId = e
                edge.edge = shape.edge(e)
                edge.maxLevel = getEdgeMaxLevel(edge.edge)
                if (edge.hasInterior) tracker.testEdge(shapeId, edge.edge)
                val clippedToFace = S2EdgeClipping.clipToPaddedFace(edge.edge.v0, edge.edge.v1, pcell.id.face(), kCellPadding)
                if (clippedToFace == null) {
                    logger.error { "Invariant failure in MutableS2ShapeIndex" }
                }
                faceEdges.add(edge)
            }
        }
        // Now create a ClippedEdge for each FaceEdge, and put them in "new_edges".
        val newEdges = mutableListOf<ClippedEdge>()
        for (face_edge in faceEdges) {
            val clipped = ClippedEdge(
                faceEdge = face_edge,
                bound = S2EdgeClipping.getClippedEdgeBound(face_edge.a, face_edge.b, pcell.bound)
            )
            newEdges.add(clipped)
        }
        // Discard any edges from "edges" that are being removed, and append the
        // remainder to "new_edges".  (This keeps the edges sorted by shape id.)
        for (i in edges.indices) {
            val clipped = edges[i]
            if (!isShapeBeingRemoved(clipped.faceEdge.shapeId)) {
                newEdges.addAll(edges.subList(i, edges.size))
                break
            }
        }
        // Update the edge list and delete this cell from the index.
        edges.clear()
        edges.addAll(newEdges)
        removeShapeIndexCell(pcell.id)
    }

    /**
     * Attempt to build an index cell containing the given edges, and return true if successful.
     * (Otherwise the edges should be subdivided further.)
     *
     * @param pcell
     * @param edges The edges to index.
     * @param tracker The shape interior tracker.
     *
     * @return true if a cell has been created and false otherwise.
     */
    private fun makeIndexCell(pcell: S2PaddedCell, edges: List<ClippedEdge>, tracker: InteriorTracker): Boolean {
        logger.trace { "makeIndexCell(pcell = $pcell, edges = [${edges.joinToString { e -> "(${e.faceEdge.shapeId}, ${e.faceEdge.edgeId})" }}], tracker = ${tracker.shapeIds()})" }
        if (edges.isEmpty() && tracker.shapeIds().isEmpty()) {
            // No index cell is needed. (In most cases this situation is detected before we get to this point, but this
            // can happen when all shapes in a cell are removed.)
            logger.trace { "Cell $pcell: edges is empty and tracker shape ids is empty. makeIndexCell = true" }
            return true
        }

        // Count the number of edges that have not reached their maximum level yet.
        // Return false if there are too many such edges.
        var count = 0
        for (edge in edges) {
            count += if (pcell.level < edge.faceEdge.maxLevel) 1 else 0
            if (count > options.maxEdgesPerCell) {
                logger.trace { "Cell $pcell contains more than ${options.maxEdgesPerCell} edges. makeIndexCell = false" }
                return false
            }
        }

        // Possible optimization: Continue subdividing as long as exactly one child of "pcell" intersects the given
        // edges.  This can be done by finding the bounding box of all the edges and calling ShrinkToFit():
        //
        // val cellid = pcell.shrinkToFit(getRectBound(edges))
        //
        // Currently this is not beneficial; it slows down construction by 4-25% (mainly computing the union of the
        // bounding rectangles) and also slows down queries (since more recursive clipping is required to get down to
        // the level of a spatial index cell).  But it may be worth trying again once "contains_center" is computed and
        // all algorithms are modified to take advantage of it.

        // We update the InteriorTracker as follows.  For every S2Cell in the index we construct two edges: one edge
        // from entry vertex of the cell to its center, and one from the cell center to its exit vertex.  Here "entry"
        // and "exit" refer the S2CellId ordering, i.e. the order in which points are encountered along the S2
        // space-filling curve.  The exit vertex then becomes the entry vertex for the next cell in the index, unless
        // there are one or more empty intervening cells, in which case the InteriorTracker state is unchanged because
        // the intervening cells have no edges.

        // Shift the InteriorTracker focus point to the center of the current cell.
        if (tracker.isActive() && edges.isNotEmpty()) {
            if (!tracker.atCellId(pcell.id)) {
                tracker.moveTo(pcell.getEntryVertex())
            }
            tracker.drawTo(pcell.getCenter())
            testAllEdges(edges, tracker)
        }
        logger.trace { "tracker = ${tracker.shapeIds()})" }


        // Allocate and fill a new index cell.  To get the total number of shapes we need to merge the shapes
        // associated with the intersecting edges together with the shapes that happen to contain the cell center.
        val containingShapeIds: ShapeIdSet = tracker.shapeIds()
        val cellNumShapes = countShapes(edges, containingShapeIds)
        val cell = S2ShapeIndexCell()

        // To fill the index cell we merge the two sources of shapes: "edge shapes" (those that have at least one edge
        // that intersects this cell), and "containing shapes" (those that contain the cell center).  We keep track
        // of the index of the next intersecting edge and the next containing shape as we go along.  Both sets of shape
        // ids are already sorted.
        var currentEdgeIndex = 0
        var currentContainingShapeIndex = 0
        repeat(cellNumShapes) {
            val currentEdgeShapeId = if (currentEdgeIndex < edges.size) edges[currentEdgeIndex].faceEdge.shapeId else Int.MAX_VALUE
            val currentContainingShapeId =
                if (currentContainingShapeIndex < containingShapeIds.size) containingShapeIds[currentContainingShapeIndex] else Int.MAX_VALUE

            if (currentContainingShapeId < currentEdgeShapeId) {
                // The entire cell is in the shape interior.
                cell.addClipped(S2ClippedShape(currentContainingShapeId, emptyList(), true))
                ++currentContainingShapeIndex
            } else {
                val clippedShapeEdges = mutableListOf<Int>()
                while (currentEdgeIndex < edges.size && edges[currentEdgeIndex].faceEdge.shapeId == currentEdgeShapeId) {
                    clippedShapeEdges.add(edges[currentEdgeIndex].faceEdge.edgeId)
                    ++currentEdgeIndex
                }
                val containsCenter = if (currentContainingShapeId == currentEdgeShapeId) {
                    ++currentContainingShapeIndex
                    true
                } else false

                cell.addClipped(S2ClippedShape(currentEdgeShapeId, clippedShapeEdges, containsCenter))

            }
        }

        // UpdateEdges() visits cells in increasing order of S2CellId, so during
        // initial construction of the index all insertions happen at the end.  It
        // is much faster to give an insertion hint in this case.  Otherwise the
        // hint doesn't do much harm.  With more effort we could provide a hint even
        // during incremental updates, but this is probably not worth the effort.
        setShapeIndexCell(pcell.id, cell)

        // Shift the InteriorTracker focus point to the exit vertex of this cell.
        if (tracker.isActive() && edges.isNotEmpty()) {
            tracker.drawTo(pcell.getExitVertex())
            testAllEdges(edges, tracker)
            tracker.setNextCellId(pcell.id.next())
        }

        logger.trace { "MakeCellIndex: id = ${pcell.id} => $cell" }
        return true
    }

    companion object {

        private val logger = KotlinLogging.logger(AbstractMutableS2ShapeIndex::class.java.name)

        /** The cell size relative to the length of an edge at which it is first considered to be "long". Long edges do
         * not contribute toward the decision to subdivide a cell further.  For example, a value of 2.0 means that the
         * cell must be at least twice the size of the edge in order for that edge to be counted. There are two reasons
         * for not counting long edges:
         * (1) such edges typically need to be propagated to several children, which increases time and memory costs
         *     without much benefit, and
         * (2) in pathological cases, many long edges close together could force subdivision to continue all the way to
         * the leaf cell level.
         */
        var kIndexCellSizeToLongEdgeRatio = 1.0

        /**
         * The default maximum number of edges per cell (not counting "long" edges). If a cell has more than this many
         * edges, and it is not a leaf cell, then it is subdivided. This flag can be overridden via
         * MutableS2ShapeIndex.Options.
         * Reasonable values range from 10 to about 50 or so.
         */
        var kIndexDefaultMaxEdgesPerCell = 10

        /** The amount by which cells are "padded" to compensate for numerical errors when clipping line segments to
         * cell boundaries.
         * The total error when clipping an edge comes from two sources:
         * (1) Clipping the original spherical edge to a cube face (the "face edge").
         *     The maximum error in this step is S2EdgeClipping.kFaceClipErrorUVCoord.
         * (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
         *     The maximum error in this step is S2EdgeClipping.kEdgeClipErrorUVCoord.
         * Finally, since we encounter the same errors when clipping query edges, we double the total error so that we
         * only need to pad edges during indexing and not at query time. */
        val kCellPadding: Double = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kEdgeClipErrorUVCoord)

        /**
         * Defines the initial focus point of MutableS2ShapeIndex.InteriorTracker (the start of the S2CellId
         * space-filling curve).
         *
         * TODO(ericv): Move InteriorTracker here to avoid the need for this method.
         */
        val kInteriorTrackerOrigin: S2Point = S2Coords.faceUvToXyz(0, -1.0, -1.0).normalize()

        // Call tracker->TestEdge() on all edges from shapes that have interiors.
        private fun testAllEdges(edges: List<ClippedEdge>, tracker: InteriorTracker) {
            for (edge in edges) {
                val faceEdge = edge.faceEdge
                if (faceEdge.hasInterior) {
                    tracker.testEdge(faceEdge.shapeId, faceEdge.edge)
                }
            }
        }

        // Return the number of distinct shapes that are either associated with the
        // given edges, or that are currently stored in the InteriorTracker.
        private fun countShapes(edges: List<ClippedEdge>, cshape_ids: List<Int>): Int {
            val edgeIterator = edges.iterator()
            var cShapeIndex = 0
            var lastShapeId = -1
            var count = 0
            while (edgeIterator.hasNext()) {
                val clippedEdge = edgeIterator.next()
                if (clippedEdge.faceEdge.shapeId != lastShapeId) {
                    ++count
                    lastShapeId = clippedEdge.faceEdge.shapeId

                    while (cShapeIndex < cshape_ids.size && cshape_ids[cShapeIndex] <= lastShapeId) {
                        if (cshape_ids[cShapeIndex] < lastShapeId) ++count
                        ++cShapeIndex
                    }
                }

            }
            while (cShapeIndex < cshape_ids.size) {
                ++count
                ++cShapeIndex
            }

            return count
        }

        // TODO(fmeurisse) check memory estimation with kotlin objects.
        // The following memory estimates are based on heap profiling.
        //
        // The final size of a MutableS2ShapeIndex depends mainly on how finely the
        // index is subdivided, as controlled by Options::max_edges_per_cell() and
        // --s2shape_index_default_max_edges_per_cell. For realistic values of
        // max_edges_per_cell() and shapes with moderate numbers of edges, it is
        // difficult to get much below 8 bytes per edge.  [The minimum possible size
        // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
        // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
        // shapes_ vector.]
        //
        // The temporary memory consists mainly of the FaceEdge and ClippedEdge
        // structures plus a ClippedEdge pointer for every level of recursive
        // subdivision.  For very large indexes this can be 200 bytes per edge.
        const val kFinalBytesPerEdge = 8
        const val kTmpBytesPerEdge = 200
        const val kTmpMemoryBudgetBytes = 100 shl 20// static_cast<size_t>(FLAGS_s2shape_index_tmp_memory_budget_mb) << 20

        // We arbitrarily limit the number of batches just as a safety measure.
        // With the current default memory budget of 100 MB, this limit is not
        // reached even when building an index of 350 million edges.
        const val kMaxUpdateBatches = 100

    }

    /**
     * UpdateState holds temporary data related to thread synchronization. It is only allocated while updates are being
     * applied.
     */
    class UpdateState(
        // This mutex is used as a condition variable.  It is locked by the
        // updating thread for the entire duration of the update; other threads
        // lock it in order to wait until the update is finished.
        val waitMutex: ReentrantLock = ReentrantLock(),

        // The number of threads currently waiting on "wait_mutex_".  The
        // UpdateState can only be freed when this number reaches zero.
        //
        // Reads and writes to this field are guarded by "lock".
        var numWaiting: Int = 0
    ) {

        fun destroy() {
            checkEQ(0, numWaiting)
        }
    }

}

/**
 * The representation of a shape that has been queued for removal.
 *
 * @property shapeId The identifier of the removed shape.
 * @property hasInterior Indicates if the shape has an interior (dimension = 2)
 * @property containsTrackerOrigin Indicates if the shape contains the interior tracker origin.
 * @property edges The shape edges.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
data class RemovedShape(
    val shapeId: Int,
    val hasInterior: Boolean,
    val containsTrackerOrigin: Boolean,
    val edges: List<Edge>
)

/**
 * ClippedEdge represents the portion of that edge that has been clipped to a given S2Cell.
 *
 * @property faceEdge The original unclipped edge
 * @property bound Bounding box for the clipped portion
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
data class ClippedEdge(
    val faceEdge: FaceEdge = FaceEdge(),
    val bound: R2Rect = R2Rect()
)

/**
 * FaceEdge represents an edge that has been projected onto a given face.
 *
 * @property shapeId The shape that this edge belongs to
 * @property edgeId Edge id within that shape
 * @property maxLevel Not desirable to subdivide this edge beyond this level
 * @property hasInterior  Belongs to a shape of dimension 2.
 * @property a The start edge endpoint, clipped to a given face
 * @property b The end edge endpoint, clipped to a given face
 * @property edge The edge endpoints.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
data class FaceEdge(
    val shapeId: Int = -1,
    var edgeId: Int = 0,
    var maxLevel: Int = 0,
    val hasInterior: Boolean = false,
    var a: R2Point = R2Point(),
    var b: R2Point = R2Point(),
    var edge: Edge = Edge()
)

/**
 * A BatchDescriptor represents a set of pending updates that will be applied at the same time. The batch consists of
 * all updates with shape ids between the current value of "pendingAdditionsBegin" (inclusive) and  "additionsEnd"
 * (exclusive). The first batch to be processed also implicitly includes all shapes being removed.
 *
 * @property additionsEnd Exclusive end of the batch.
 * @property numEdges total number of edges that will be added or removed in this batch.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
data class BatchDescriptor(
    var additionsEnd: Int,
    val numEdges: Int
)

