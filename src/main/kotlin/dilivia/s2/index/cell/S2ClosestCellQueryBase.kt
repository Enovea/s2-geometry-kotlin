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

import dilivia.PreConditions.checkGE
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireGE
import dilivia.s2.S2CellId
import dilivia.s2.index.Delta
import dilivia.s2.index.Distance
import dilivia.s2.index.DistanceFactory
import dilivia.s2.index.S2DistanceTarget
import dilivia.s2.index.cell.S2CellIndex.NonEmptyRangeIterator
import dilivia.s2.index.cell.S2CellIndex.RangeIterator
import dilivia.s2.index.lexicon.Label
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Region
import dilivia.s2.region.S2RegionCoverer
import mu.KotlinLogging
import java.util.*

/**
 * S2ClosestCellQueryBase is a templatized class for finding the closest (cell_id, value) pairs in an S2CellIndex to a
 * given target.
 *
 * It is not intended to be used directly, but rather to serve as the implementation of various specialized classes
 * with more convenient APIs (such as S2ClosestCellQuery).  It is flexible enough so that it can be adapted to compute
 * maximum distances and even potentially Hausdorff distances.
 *
 * By using the appropriate options, this class can answer questions such as:
 *
 *  - Find the minimum distance between a cell collection A and a target B.
 *  - Find all cells in collection A that are within a distance D of target B.
 *  - Find the k cells of collection A that are closest to a given point P.
 *
 * The target is any class that implements the S2DistanceTarget interface. There are predefined targets for points,
 * edges, S2Cells, S2CellUnions, and S2ShapeIndexes (arbitrary collections of points, polylines, and polygons).
 *
 * The Distance template argument is used to represent distances.  Usually it is a thin wrapper around S1ChordAngle,
 * but another distance type may be used as long as it implements the Distance concept described in
 * S2DistanceTarget.  For example this can be used to measure maximum distances, to get more accuracy,
 * or to measure non-spheroidal distances.
 *
 * This class is a port of the S2ClosestCellQueryBase class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T Distance used to search cells.
 * @param V The type of data stored by the index.
 * @property distanceFactory Factory to build distance of type T.
 * @constructor requires init() to be called.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2ClosestCellQueryBase<T : Distance<T>>(private val distanceFactory: DistanceFactory<T>) {

    private val uid = UUID.randomUUID().toString()

    /////////////// Query parameters /////////////////

    /** Index on which the query is done. */
    private lateinit var index: S2CellIndex

    /**
     * Get the underlying S2CellIndex.
     *
     * @return the index.
     */
    fun index(): S2CellIndex = index

    /** Query options. */
    private lateinit var options: Options<T>

    /** Current query target. */
    private lateinit var target: S2DistanceTarget<T>

    /**
     * Indicates if maxError() must be subtracted from priority queue cell distances in order to ensure that such
     * distances are measured conservatively.  This is true only if the target takes advantage of maxError() in order to
     * return faster results, and 0 < maxError() < distanceLimit.
     */
    private var useConservativeCellDistance: Boolean = false

    /**
     * The distance beyond which we can safely ignore further candidate cells. (Candidates that are exactly at the limit
     * are ignored; this is more efficient for updateMinDistance() and should not affect clients since distance
     * measurements have a small amount of error anyway.)
     *
     * Initially this is the same as the maximum distance specified by the user, but it can also be updated by the
     * algorithm (see maybeAddResult).
     */
    private lateinit var distanceLimit: T

    /**
     * For the optimized algorithm we precompute the top-level S2CellIds that will be added to the priority queue.
     * There can be at most 6 of these cells. Essentially this is just a covering of the indexed cells.
     */
    private var indexCovering: ArrayList<S2CellId> = ArrayList(6)

    /////////////// Results /////////////////

    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.

    // The current result set is stored in one of three ways:

    /** If maxResults == 1, contains the best query result. */
    private lateinit var resultSingleton: S2ClosestCellQueryResult<T>

    /** If maxResults == kMaxMaxResults, results are appended to resultList and sorted/uniqued at the end.*/
    private var resultList: MutableList<S2ClosestCellQueryResult<T>> = mutableListOf()

    /**
     * If 1 < maxResults < kMaxMaxResults, results are kept in a Set so that we can progressively reduce the distance
     * limit once maxResults results have been found. (A priority queue is not sufficient because we need to be able to
     * check whether a candidate cell is already in the result set.
     */
    private var resultSet: TreeSet<S2ClosestCellQueryResult<T>> = TreeSet()

    /////////////// Internal /////////////////

    // Used to iterate over the contents of an S2CellIndex range.  It is defined
    // here to take advantage of the fact that when multiple ranges are visited
    // in increasing order, duplicates can automatically be eliminated.
    private lateinit var contentsIt: S2CellIndex.ContentsIterator

    /**
     * When the results are stored in a btree_set (see above), usually duplicates can be removed simply by inserting
     * candidate cells in the current result set.
     * However this is not true if Options.maxError > 0 and the Target subtype takes advantage of this by returning
     * suboptimal distances. This is because when updateMinDistance() is called with different "min_dist" parameters
     * (i.e., the distance to beat), the implementation may return a different distance for the same cell.
     * Since the set is keyed by (distance, cell_id, value) this can create duplicate results.
     *
     * The flag below is true when duplicates must be avoided explicitly. This is achieved by maintaining a separate
     * set keyed by (cell_id, value) only, and checking whether each edge is in that set before computing the distance
     * to it.
     *
     * TODO(ericv): Check whether it is faster to avoid duplicates by default (even when Options::max_results() == 1), rather than just when we need to.
     */
    private var avoidDuplicates: Boolean = true

    /** Maintains the tested cells in order to avoid duplicates when avoidDuplicates == true. */
    private val testedCells: MutableSet<ValuedCell<Label>> = HashSet()

    /**
     * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted in increasing order of distance from
     * the target.
     */
    private val queue: Queue<QueueEntry<T>> = PriorityQueue(16)

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private val maxDistanceCovering: MutableList<S2CellId> = mutableListOf()
    private val intersectionWithMaxDistance: MutableList<S2CellId> = mutableListOf()

    /////////////// Init /////////////////

    init {
        testedCells.add(ValuedCell(S2CellId.none))
    }

    /**
     * Convenience constructor that calls init().
     *
     * @param distanceFactory Factory to build distance of type T.
     * @param index The index on which the query is done.
     */
    constructor(distanceFactory: DistanceFactory<T>, index: S2CellIndex) : this(distanceFactory) {
        init(index)
    }

    /**
     * Initializes the query.
     * REQUIRES: reInit() must be called if "index" is modified.
     *
     * @param index The underlying index.
     */
    fun init(index: S2CellIndex) {
        this.index = index
        this.contentsIt = S2CellIndex.ContentsIterator(index)
        reInit()
    }

    /**
     *  Reinitializes the query.  This method must be called whenever the underlying index is modified.
     */
    fun reInit() {
        indexCovering.clear()
        distanceLimit = distanceFactory.infinity()
    }

    /////////////// Find cells /////////////////

    /**
     * Returns the closest (cell_id, value) pairs to the given target that satisfy the given options. This method may be
     * called multiple times.
     *
     * @param target The target distance
     * @param options Query options.
     */
    fun findClosestCells(target: S2DistanceTarget<T>, options: Options<T>): List<S2ClosestCellQueryResult<T>> {
        val results = mutableListOf<S2ClosestCellQueryResult<T>>()
        findClosestCells(target, options, results)
        return results
    }

    /**
     *  Returns the closest (cell_id, value) pairs to the given target that satisfy the given options.
     *
     * This version can be more efficient when this method is called many times, since it does not require allocating a
     * new list on each call.
     *
     * @param target The target distance
     * @param options Query options.
     * @param results The result list to populate.
     *
     * @return The populated results parameter.
     */
    fun findClosestCells(target: S2DistanceTarget<T>, options: Options<T>, results: MutableList<S2ClosestCellQueryResult<T>>): MutableList<S2ClosestCellQueryResult<T>> {
        logger.trace { """
            |
            |===========================================================================================================
            |findClosestCells: $uid
            |-----------------------------------------------------------------------------------------------------------
            |target: $target,
            |options: $options
            |===========================================================================================================
        """.trimMargin() }
        findClosestCellsInternal(target, options)
        results.clear()
        if (options.maxResults == 1) {
            if (!resultSingleton.isEmpty()) {
                results.add(resultSingleton)
            }
        } else if (options.maxResults == Options.kMaxMaxResults) {
            results.addAll(resultList)
            results.sort()
            resultList.clear()
        } else {
            results.addAll(resultSet)
            resultSet.clear()
        }

        logger.trace { """
            |
            |===========================================================================================================
            |findClosestCells: $uid - Results
            |-----------------------------------------------------------------------------------------------------------
            |${results.joinToString("\n")}
            |===========================================================================================================
        """.trimMargin() }
        return results
    }

    /**
     * Convenience method that returns exactly one (cell_id, label) pair. If no cells satisfy the given search criteria,
     * then a Result with distance() == Infinity() and isEmpty() == true is returned.
     *
     * REQUIRES: options.getMaxResults() == 1
     *
     * @param target The target distance
     * @param options Query options.
     *
     * @return The closest result.
    */
    fun findClosestCell(target: S2DistanceTarget<T>, options: Options<T>): S2ClosestCellQueryResult<T> {
        logger.trace { """
            |
            |===========================================================================================================
            |findClosestCells: $uid
            |-----------------------------------------------------------------------------------------------------------
            |target: $target,
            |options: $options
            |===========================================================================================================
        """.trimMargin() }
        requireEQ(options.maxResults, 1)
        findClosestCellsInternal(target, options)

        logger.trace { """
            |
            |===========================================================================================================
            |findClosestCells: $uid - Results
            |-----------------------------------------------------------------------------------------------------------
            |$resultSingleton
            |===========================================================================================================
        """.trimMargin() }
        return resultSingleton
    }

    /////////////// Options /////////////////

    /**
     * Options that control the set of cells returned. Note that by default *all* cells are returned, so you will always
     * want to set either the maxResults option or the maxDistance option (or both).
     *
     * @property maxResults Specifies that at most "maxResults" cells should be returned.
     *
     * @param T Distance used to search cells.
     */
    open class Options<T : Distance<T>>(
        distanceFactory: DistanceFactory<T>,
        maxResults: Int = kMaxMaxResults,
        maxDistance: T = distanceFactory.infinity(),
        maxError: Delta = Delta.zero(),
        region: S2Region? = null,
        useBruteForce: Boolean = false
    ) : Cloneable {

        /** Factory to build distance of type T. */
        @Suppress("CanBePrimaryConstructorProperty")
        val distanceFactory = distanceFactory

        /**
         * Specifies that at most "maxResults" cells should be returned.
         *
         * REQUIRES: max_results >= 1
         * DEFAULT: kMaxMaxResults
         */
        var maxResults = maxResults
            set(value) {
                requireGE(value, 1)
                field = value
            }

        /**
         * Specifies that only cells whose distance to the target is less than "maxDistance" should be returned.
         *
         * Note that cells whose distance is exactly equal to "max_distance" are not returned. In most cases this
         * doesn't matter (since distances are not computed exactly in the first place), but if such cells are needed
         * then you can retrieve them by specifying "max_distance" as the next largest representable Distance.
         * For example, if Distance is an S1ChordAngle then you can specify maxDistance.successor().
         *
         * DEFAULT: Distance.infinity
         */
        @Suppress("CanBePrimaryConstructorProperty")
        var maxDistance = maxDistance

        /**
         * Specifies that cells up to maxError further away than the true closest cells may be substituted in the
         * result set, as long as such cells satisfy all the remaining search criteria (such as maxDistance).
         *
         * This option only has an effect if maxResults is also specified otherwise all cells closer than maxDistance
         * will always be returned.
         *
         * Note that this does not affect how the distance between cells is computed; it simply gives the algorithm
         * permission to stop the search early as soon as the best possible improvement drops below maxError.
         *
         * This can be used to implement distance predicates efficiently. For example, to determine whether the minimum
         * distance is less than D, the isDistanceLess() method sets maxResults == 1 and maxDistance == maxError == D.
         *
         * This causes the algorithm to terminate as soon as it finds any cell whose distance is less than D, rather
         * than continuing to search for a cell that is even closer.
         *
         * DEFAULT: Delta.zero)
         */
        @Suppress("CanBePrimaryConstructorProperty")
        var maxError = maxError

        /**
         * Specifies that cells must intersect the given S2Region.  "region" is owned by the caller and must persist
         * during the lifetime of this object. The value may be changed between calls to findClosestPoints(), or reset
         * by calling region = null.
         *
         * Note that if you want to set the region to a disc around a target point, it is faster to use a PointTarget
         * with maxDistance) instead.  You can also call both methods, e.g. to set a maximum distance and also require
         * that cells lie within a given rectangle.
         */
        @Suppress("CanBePrimaryConstructorProperty")
        var region = region

        /**
         * Specifies that distances should be computed by examining every cell rather than using the S2ShapeIndex.
         * This is useful for testing, benchmarking, and debugging.
         *
         * DEFAULT: false
         */
        @Suppress("CanBePrimaryConstructorProperty")
        var useBruteForce = useBruteForce

        override fun clone(): Options<T> {
            return Options(
                distanceFactory, maxResults, maxDistance, maxError, region, useBruteForce
            )
        }

        override fun toString(): String {
            return "Options(maxResults=$maxResults, maxDistance=$maxDistance, maxError=$maxError, region=$region, useBruteForce=$useBruteForce)"
        }

        companion object {

            const val kMaxMaxResults = Int.MAX_VALUE

        }
    }

    companion object {

        private val logger = KotlinLogging.logger(S2ClosestCellQueryBase::class.java.name)

        // The minimum number of ranges that a cell must contain to enqueue it
        // rather than processing its contents immediately.
        const val kMinRangesToEnqueue = 6
    }

    /////////////// Internal methods /////////////////

    data class QueueEntry<T : Distance<T>>(
        // A lower bound on the distance from the target to "id".  This is the key of the priority queue.
        val distance: T,
        // The cell being queued.
        val id: S2CellId
    ) : Comparable<QueueEntry<T>> {

        // The priority queue returns the smallest elements first
        override fun compareTo(other: QueueEntry<T>): Int = distance.compareTo(other.distance)

    }

    private fun findClosestCellsInternal(target: S2DistanceTarget<T>, options: Options<T>) {
        this.target = target
        this.options = options

        testedCells.clear()
        contentsIt.clear()
        distanceLimit = options.maxDistance.clone()
        resultSingleton = S2ClosestCellQueryResult(distanceFactory.infinity())
        checkState { resultList.isEmpty() }
        checkState { resultSet.isEmpty() }
        requireGE(target.maxBruteForceIndexSize(), 0)

        if (distanceLimit == distanceFactory.zero()) {
            return
        }

        if (options.maxResults == Options.kMaxMaxResults &&
            options.maxDistance == distanceFactory.infinity() &&
            options.region == null
        ) {
            logger.warn { "Returning all cells (max_results/max_distance/region not set)" }
        }

        // If max_error() > 0 and the target takes advantage of this, then we may
        // need to adjust the distance estimates to the priority queue cells to
        // ensure that they are always a lower bound on the true distance.  For
        // example, suppose max_distance == 100, max_error == 30, and we compute the
        // distance to the target from some cell C0 as d(C0) == 80.  Then because
        // the target takes advantage of max_error(), the true distance could be as
        // low as 50.  In order not to miss edges contained by such cells, we need
        // to subtract max_error() from the distance estimates.  This behavior is
        // controlled by the use_conservative_cell_distance_ flag.
        //
        // However there is one important case where this adjustment is not
        // necessary, namely when max_distance() < max_error().  This is because
        // max_error() only affects the algorithm once at least max_edges() edges
        // have been found that satisfy the given distance limit.  At that point,
        // max_error() is subtracted from distance_limit_ in order to ensure that
        // any further matches are closer by at least that amount.  But when
        // max_distance() < max_error(), this reduces the distance limit to 0,
        // i.e. all remaining candidate cells and edges can safely be discarded.
        // (Note that this is how IsDistanceLess() and friends are implemented.)
        //
        // Note that Distance::Delta only supports operator==.
        val targetUsesMaxError = (options.maxError != Delta.zero() && target.setMaxError(options.maxError))

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a Delta and one is a Distance.  Instead we subtract them.
        useConservativeCellDistance = targetUsesMaxError &&
                (distanceLimit == distanceFactory.infinity() || distanceFactory.zero() < distanceLimit - options.maxError)

        // Use the brute force algorithm if the index is small enough.
        if (options.useBruteForce || index.numCells <= target.maxBruteForceIndexSize()) {
            avoidDuplicates = false
            findClosestCellsBruteForce()
        } else {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoidDuplicates = (targetUsesMaxError && options.maxResults > 1)
            findClosestCellsOptimized()
        }
    }

    private fun findClosestCellsBruteForce() {
        logger.trace { """
            |
            |-------------------------------------------------------
            |findClosestCellsBruteForce: $uid
            |-------------------------------------------------------
        """.trimMargin() }
        val iter = S2CellIndex.CellIterator(index)
        while (!iter.done()) {
            maybeAddResult(iter.cellId(), iter.value())
            iter.next()
        }

        logger.trace { """
            |
            |-------------------------------------------------------
            |findClosestCells: $uid - End
            |-------------------------------------------------------
        """.trimMargin() }
    }

    private fun findClosestCellsOptimized() {
        logger.trace { """
            |
            |-------------------------------------------------------
            |findClosestCellsOptimized: $uid
            |-------------------------------------------------------
        """.trimMargin() }
        initQueue()
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to remove
            // it before adding any new entries to the queue.
            val entry = queue.poll()
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            val distance = entry.distance
            if (distance >= distanceLimit) {
                queue.clear()  // Clear any remaining entries.
                break
            }
            var child = entry.id.childBegin()
            // We already know that it has too many cells, so process its children.
            // Each child may either be processed directly or enqueued again.  The
            // loop is optimized so that we don't seek unnecessarily.
            var seek = true
            val range = NonEmptyRangeIterator(index)
            for (i in 0..3) {
                seek = processOrEnqueue(child, range, seek)
                child = child.next()
            }
        }

        logger.trace { """
            |
            |-------------------------------------------------------
            |findClosestCellsOptimized: $uid - End
            |-------------------------------------------------------
        """.trimMargin() }
    }

    private fun initQueue() {
        logger.trace { """
            |
            |------------------------------
            | Init Queue: $uid
            |------------------------------
        """.trimMargin() }
        checkState { queue.isEmpty() }

        // Optimization: rather than starting with the entire index, see if we can
        // limit the search region to a small disc.  Then we can find a covering for
        // that disc and intersect it with the covering for the index.  This can
        // save a lot of work when the search region is small.
        val cap = target.getCapBound()
        logger.trace { "$uid: Cap = $cap" }

        if (cap.isEmpty) {
            logger.trace { "$uid: cap is empty, return" }
            return // Empty target.
        }
        if (options.maxResults == 1) {
            // If the user is searching for just the closest cell, we can compute an
            // upper bound on search radius by seeking to the center of the target's
            // bounding cap and looking at the contents of that leaf cell range.  If
            // the range intersects any cells, then the distance is zero.  Otherwise
            // we can still look at the two neighboring ranges, and use the minimum
            // distance to any cell in those ranges as an upper bound on the search
            // radius.  These cells may wind up being processed twice, but in general
            // this is still faster.
            //
            // First check the range containing or immediately following "center".
            val range = NonEmptyRangeIterator(index)
            val target = S2CellId.fromPoint(cap.center)
            range.seek(target)
            addRange(range)
            if (distanceLimit == distanceFactory.zero()) return

            // If the range immediately follows "center" (rather than containing it),
            // then check the previous non-empty range as well.
            if (range.startId() > target && range.prev()) {
                addRange(range)
                if (distanceLimit == distanceFactory.zero()) return
            }
        }

        // We start with a covering of the set of indexed cells, then intersect it
        // with the maximum search radius disc (if any).
        //
        // Note that unlike S2ClosestPointQuery, we can't also intersect with the
        // given region (if any).  This is because the index cells in the result are
        // only required to intersect the region.  This means that an index cell that
        // intersects the region's covering may be much closer to the target than the
        // covering itself, which means that we cannot use the region's covering to
        // restrict the search.
        //
        // TODO(ericv): If this feature becomes important, this could be fixed by
        // (1) computing a covering of the region, (2) looking up any index cells
        // that contain each covering cell by seeking to covering_cell.range_min(),
        // (3) replacing each covering cell by the largest such cell (if any), and
        // (4) normalizing the result.
        if (indexCovering.isEmpty()) initCovering()
        logger.trace { "$uid: index covering = $indexCovering" }

        var initialCells = indexCovering.toList()
        if (distanceLimit < distanceFactory.infinity()) {
            val coverer = S2RegionCoverer(maxCells = 4)
            val radius = cap.radius + distanceLimit.getChordAngleBound()
            val searchCap = S2Cap(cap.center, radius)
            coverer.getFastCovering(searchCap, maxDistanceCovering)
            S2CellUnion.getIntersection(initialCells, maxDistanceCovering, intersectionWithMaxDistance)
            initialCells = intersectionWithMaxDistance.toList()
        }

        logger.trace { "$uid: initial cells: $initialCells" }
        val range = NonEmptyRangeIterator(index)
        for (i in initialCells.indices) {
            val id = initialCells[i]
            val seek = (i == 0) || id.rangeMin() >= range.limitId()
            processOrEnqueue(id, range, seek)
            if (range.done()) break
        }

        logger.trace { """
            |
            |------------------------------
            | Init Queue: $uid - Initial queue
            |------------------------------
            |${queue.joinToString("\n")}
            |------------------------------
        """.trimMargin() }

    }

    private fun initCovering() {
        // Compute the "index covering", which is a small number of S2CellIds that
        // cover the indexed cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one covering cell
        // per spanned face, just big enough to cover the indexed cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the indexed cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a covering cell that is big enough to just fit
        // those indexed cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the covering cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.
        indexCovering.ensureCapacity(6)
        val it = NonEmptyRangeIterator(index)
        var last = NonEmptyRangeIterator(index)
        it.begin()
        last.finish()
        if (!last.prev()) return  // Empty index.
        val indexLastId = last.limitId().previous()
        if (it.startId() != last.startId()) {
            // The index contains at least two distinct S2CellIds (because otherwise
            // there would only be one non-empty range).  Choose a level such that the
            // entire index can be spanned with at most 6 cells (if the index spans
            // multiple faces) or 4 cells (it the index spans a single face).
            val level = it.startId().getCommonAncestorLevel(indexLastId) + 1

            // Visit each potential covering cell except the last (handled below).
            val startId = it.startId().parent(level)
            val lastId = indexLastId.parent(level)
            var id = startId
            while (id != lastId) {
                // Skip any covering cells that don't contain an indexed range.
                if (id.rangeMax() < it.startId()) {
                    id = id.next()
                    continue
                }

                // Find the indexed range contained by this covering cell and then
                // shrink the cell if necessary so that it just covers this range.
                val cellFirstId = it.startId()
                it.seek(id.rangeMax().next())
                // Find the last leaf cell covered by the previous non-empty range.
                last = NonEmptyRangeIterator(it)
                last.prev()
                addInitialRange(cellFirstId, last.limitId().previous())
                id = id.next()
            }
        }
        addInitialRange(it.startId(), indexLastId)
    }

    // Adds a cell to index_covering_ that covers the given inclusive range.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    fun addInitialRange(firstId: S2CellId, lastId: S2CellId) {
        // Add the lowest common ancestor of the given range.
        val level = firstId.getCommonAncestorLevel(lastId)
        checkGE(level, 0)
        indexCovering.add(firstId.parent(level))
    }

    // TODO(ericv): Consider having this method return false when distance_limit_
    // is reduced to zero, and terminating any calling loops early.
    fun maybeAddResult(cellId: S2CellId, value: Label) {

        if (avoidDuplicates && !testedCells.add(ValuedCell(cellId, value))) {
            logger.trace { "$uid:  maybeAddResult(cellId: $cellId, value: $value): duplicated result, do nothing" }
            return
        }

        // TODO(ericv): It may be relatively common to add the same S2CellId
        // multiple times with different labels.  This could be optimized by
        // remembering the last "cell_id" argument and its distance.  However this
        // may not be beneficial when Options::max_results() == 1, for example.
        val cell = S2Cell(cellId)
        val distance = distanceLimit.clone()
        if (!target.updateMinDistance(cell, distance)) {
            logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): distance from cell to target is > $distanceLimit, do nothing" }
            return
        }

        val region = options.region
        if (region != null && !region.mayIntersect(cell)) {
            logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): cell does not intersect region $region, do nothing" }
            return
        }

        val result = S2ClosestCellQueryResult(distance, cellId, value)
        if (options.maxResults == 1) {
            // Optimization for the common case where only the closest cell is wanted.
            resultSingleton = result
            distanceLimit = result.distance - options.maxError
            logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): (max result = 1) cell is nearer - new distanceLimit = $distanceLimit" }
        }
        else if (options.maxResults == Options.kMaxMaxResults) {
            resultList.add(result)  // Sort/unique at end.
            logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): (max result = infinity) add cell - distanceLimit = $distanceLimit" }
        } else {
            logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): (max result = ${options.maxResults}) add cell an remove result if necessary" }
            // Add this cell to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            resultSet.add(result)
            val size = resultSet.size
            if (size >= options.maxResults) {
                if (size > options.maxResults) {
                    val last = resultSet.pollLast()
                    logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): (max result = ${options.maxResults}) reach max result, remove last $last" }
                }
                distanceLimit = resultSet.last().distance - options.maxError
                logger.trace { "$uid: maybeAddResult(cellId: $cellId, value: $value): (max result = ${options.maxResults}) reach max result, update distance limit = $distanceLimit" }
            }
        }
    }

    // Either process the contents of the given cell immediately, or add it to the
    // queue to be subdivided.  If "seek" is false, then "iter" must be positioned
    // at the first non-empty range (if any) with start_id() >= id.range_min().
    //
    // Returns "true" if the cell was added to the queue, and "false" if it was
    // processed immediately, in which case "iter" is positioned at the first
    // non-empty range (if any) with start_id() > id.range_max().
    private fun processOrEnqueue(id: S2CellId, iter: NonEmptyRangeIterator, seek: Boolean): Boolean {
        logger.trace { """
            |
            |------------------------------
            | processOrEnqueue: $uid
            |------------------------------
            | id: $id
            | seek: $seek
            |------------------------------
        """.trimMargin() }
        if (seek) {
            iter.seek(id.rangeMin())
            logger.trace { "$uid: Seek iter to ${id.rangeMin()} => ${if (!iter.done()) iter.startId() else "done"}" }
        }
        val last: S2CellId = id.rangeMax()
        if (iter.startId() > last) {
            logger.trace { "$uid: iter.startId = ${iter.startId()} > last = $last. Return false." }
            return false  // No need to seek to next child.
        }

        // If this cell intersects at least "kMinRangesToEnqueue" leaf cell ranges
        // (including ranges whose contents are empty), then enqueue it.  We test
        // this by advancing (n - 1) ranges and checking whether that range also
        // intersects this cell.
        val maxIt: RangeIterator = RangeIterator(iter)
        if (maxIt.advance(kMinRangesToEnqueue - 1) && maxIt.startId() <= last) {
            // This cell intersects at least kMinRangesToEnqueue ranges, so enqueue it.
            val cell = S2Cell(id)
            var distance = distanceLimit.clone()
            // We check "region" second because it may be relatively expensive.
            val region = options.region
            if (target.updateMinDistance(cell, distance) && (region == null || region.mayIntersect(cell))) {
                if (useConservativeCellDistance) {
                    // Ensure that "distance" is a lower bound on distance to the cell.
                    distance -= options.maxError
                }
                logger.trace { "$uid: enqueue cell $id" }
                queue.offer(QueueEntry(distance, id))
            }
            return true  // Seek to next child.
        }

        // There were few enough ranges that we might as well process them now.
        while (iter.startId() <= last) {
            addRange(iter)
            iter.next()
        }
        return false  // No need to seek to next child.
    }

    fun addRange(range: RangeIterator) {
        logger.trace { """
            |
            |------------------------------
            | addRange: $uid
            |------------------------------
            | startId: ${range.startId()}
            | limitId: ${range.limitId()}
            |------------------------------
        """.trimMargin() }

        contentsIt.startUnion(range)
        while (!contentsIt.done()) {
            maybeAddResult(contentsIt.cellId(), contentsIt.value())
            contentsIt.next()
        }
    }

}

