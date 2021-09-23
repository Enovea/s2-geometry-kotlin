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
package dilivia.s2.index.point

import com.google.common.collect.ComparisonChain
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireGE
import dilivia.collections.isSorted
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.Delta
import dilivia.s2.index.Distance
import dilivia.s2.index.DistanceFactory
import dilivia.s2.index.S2DistanceTarget
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Region
import dilivia.s2.region.S2RegionCoverer
import mu.KotlinLogging
import java.util.*

/**
 * S2ClosestPointQueryBase is a templatized class for finding the closest point(s) to a given target.  It is not
 * intended to be used directly, but rather to serve as the implementation of various specialized classes with
 * more convenient APIs (such as S2ClosestPointQuery).  It is flexible enough so that it can be adapted to compute
 * maximum distances and even potentially Hausdorff distances.
 *
 * By using the appropriate options, this class can answer questions such as:
 *
 *  - Find the minimum distance between a point collection A and a target B.
 *  - Find all points in collection A that are within a distance D of target B.
 *  - Find the k points of collection A that are closest to a given point P.
 *
 * The target is any class that implements the S2DistanceTarget interface. There are predefined targets for points,
 * edges, S2Cells, and S2ShapeIndexes (arbitrary collections of points, polylines, and polygons).
 *
 * The Distance template argument is used to represent distances. Usually it is a thin wrapper around S1ChordAngle,
 * but another distance type may be used as long as it implements the Distance concept described in S2DistanceTargets.
 * For example this can be used to measure maximum distances, to get more accuracy, or to measure non-spheroidal
 * distances.
 *
 * This class is a port of the S2CellIndex class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T The distance type.
 * @param D The type of indexed data.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2ClosestPointQueryBase<T: Distance<T>, D : Comparable<D>> {

    ////////////////////////// Query parameters /////////////////////////////

    private val uid: String = UUID.randomUUID().toString()

    /** The distance factory. */
    private val distanceFactory: DistanceFactory<T>
    /** The point index. */
    private lateinit var index: S2PointIndex<D>
    /**
     * Get a reference to the underlying S2PointIndex.
     * @return The index of the query.
     */
    fun index(): S2PointIndex<D> = index
    /** The query options. */
    private lateinit var options: Options<T>
    /** The query target. */
    private lateinit var target: S2DistanceTarget<T>
    /**
     * True if maxEror() must be subtracted from priority queue cell distances in order to ensure that such distances
     * are measured conservatively. This is true only if the target takes advantage of maxError() in order to return
     * faster results, and 0 < maxError() < distanceLimit.
     */
    private var useConservativeCellDistance: Boolean = false

    /**
     * For the optimized algorihm we precompute the top-level S2CellIds that will be added to the priority queue. There
     * can be at most 6 of these cells. Essentially this is just a covering of the indexed points.
     */
    private val indexCovering = ArrayList<S2CellId>()

    // The distance beyond which we can safely ignore further candidate points.
    // (Candidates that are exactly at the limit are ignored; this is more
    // efficient for UpdateMinDistance() and should not affect clients since
    // distance measurements have a small amount of error anyway.)
    //
    // Initially this is the same as the maximum distance specified by the user,
    // but it can also be updated by the algorithm (see MaybeAddResult).
    private lateinit var distanceLimit: T

    ////////////////////////// Results /////////////////////////////

    // The current result set is stored in one of three ways:
    //
    //  - Otherwise results are kept in a priority queue so that we can
    //    progressively reduce the distance limit once max_results() results
    //    have been found.
    /** If maxResults == 1, keep the best result. */
    private var resultSingleton: Result<T, D>
    /** If maxResults == "infinity", results are appended to resultVector and sorted/uniqued at the end. */
    private val resultVector = mutableListOf<Result<T, D>>()
    /**
     * If 1 < maxResults < "infinity", results are kept in a priority queue so that we can progressively reduce the
     * distance limit once maxResults results have been found.
     */
    private val resultSet: PriorityQueue<Result<T, D>> = PriorityQueue(16) { d1, d2 -> -d1.compareTo(d2)}

    ////////////////////////// Internal /////////////////////////////

    private var queue: Queue<QueueEntry<T>> = PriorityQueue(16)

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private lateinit var iter: S2PointIndex.Iterator<D>
    private val regionCovering = mutableListOf<S2CellId>()
    private val maxDistanceCovering = mutableListOf<S2CellId>()
    private val intersectionWithRegion = ArrayList<S2CellId>()
    private val intersectionWithMaxDistance = ArrayList<S2CellId>()
    private val tmpPointData: Array<PointData<D>?> = Array(kMinPointsToEnqueue - 1) { null }

    ////////////////////////// Init /////////////////////////////

    /**
     * Default constructor;
     * requires init() to be called.
     *
     * @param distanceFactory The distance factory for the distance type of the query.
     */
    constructor(distanceFactory: DistanceFactory<T>) {
        this.distanceFactory = distanceFactory
        this.resultSingleton = Result(distanceFactory.infinity())
    }

    /**
     * Convenience constructor that calls init().
     *
     * @param distanceFactory The distance factory for the distance type of the query.
     * @param index The index on which the query is executed.
     */
    constructor(distanceFactory: DistanceFactory<T>, index: S2PointIndex<D>): this(distanceFactory)  {
        init(index)
    }

    /**
     * Initializes the query.
     *
     * REQUIRES: ReInit() must be called if "index" is modified.
     * @param index The index on which the query is executed.
     */
    fun init(index: S2PointIndex<D>) {
        this.index = index
        this.iter = S2PointIndex.Iterator(index)
        reInit()
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    fun reInit() {
        iter.init(index)
        indexCovering.clear()
    }

    ////////////////////////// Queries /////////////////////////////
    
    /**
     * Get the closest points to the given target that satisfy the given options.
     * This method may be called multiple times.
     *
     * @param target The target of the query.
     * @param options The query options.
     * @return The query results.
     */
    fun findClosestPoints(target: S2DistanceTarget<T>, options: Options<T>): List<Result<T, D>> {
        val results = mutableListOf<Result<T, D>>()
        findClosestPoints(target, options, results)
        return results
    }

    /**
     * Get the closest points to the given target that satisfy the given options. This version can be more efficient
     * when this method is called many times, since it does not require allocating a new list on each call.
     *
     * @param target The target of the query.
     * @param options The query options.
     * @param results The list to populate with the query results.
     */
    fun findClosestPoints(target: S2DistanceTarget<T>, options: Options<T>, results: MutableList<Result<T, D>>) {
        findClosestPointsInternal(target, options);
        results.clear()
        if (options.getMaxResult() == 1) {
            val result = resultSingleton
            if (!result.isEmpty()) {
                results.add(result)
            }
        } else if (options.getMaxResult() == Options.kMaxMaxResults) {
            results.addAll(resultVector)
            results.sort()
            resultVector.clear()
        } else {
            if (results is ArrayList) results.ensureCapacity(resultSet.size)
            while(!resultSet.isEmpty()) {
                results.add(resultSet.poll())
            }
            // The priority queue returns the largest elements first.
            results.reverse()
            checkState { results.isSorted() }
        }

        logger.trace { """
            |
            |===========================================
            | findClosestPointsInternal: $uid
            |-------------------------------------------
            | ${results.joinToString("\n")}
            |===========================================
        """.trimMargin() }
    }

    /**
     * Convenience method that returns exactly one point. If no points satisfy the given search criteria, then a
     * Result with distance() == infinity() and isEmpty() == true is returned.
     *
     * REQUIRES: options.maxResults == 1
     *
     * @param target The target of the query.
     * @param options The query options.
     * @return The closest result.
     */
    fun findClosestPoint(target: S2DistanceTarget<T>, options: Options<T>): Result<T, D> {
        requireEQ(options.getMaxResult(), 1)
        findClosestPointsInternal(target, options)

        logger.trace { """
            |
            |===========================================
            | findClosestPointsInternal: $uid
            |-------------------------------------------
            | result singleton: $resultSingleton
            |===========================================
        """.trimMargin() }
        return resultSingleton
    }

    ////////////////////////// Internal /////////////////////////////

    private fun findClosestPointsInternal(target: S2DistanceTarget<T>, options: Options<T>) {

        logger.trace { """
            |
            |===========================================
            | findClosestPointsInternal: $uid
            |-------------------------------------------
            | target: $target
            | options: $options
            |===========================================
        """.trimMargin() }

        this.target = target
        this.options = options

        distanceLimit = options.maxDistance.clone()
        resultSingleton = Result(distanceFactory.infinity())
        requireArgument { resultVector.isEmpty() }
        requireArgument { resultSet.isEmpty() }
        requireGE(target.maxBruteForceIndexSize(), 0)
        if (distanceLimit == distanceFactory.zero()) return

        if (options.getMaxResult() == Options.kMaxMaxResults &&
                options.maxDistance == distanceFactory.infinity() &&
                options.region == null) {
            logger.warn { "Returning all points (max_results/max_distance/region not set)" }
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
        // max_error() only affects the algorithm once at least max_results() edges
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

        // Note that given point is processed only once (unlike S2ClosestEdgeQuery),
        // and therefore we don't need to worry about the possibility of having
        // duplicate points in the results.
        if (options.useBruteForce || index.numPoints() <= target.maxBruteForceIndexSize()) {
            logger.trace { "$uid: use brute force" }
            findClosestPointsBruteForce()
        } else {
            logger.trace { "$uid: use optimized" }
            findClosestPointsOptimized()
        }
    }

    private fun findClosestPointsBruteForce() {
        iter.begin()
        while (!iter.done()) {
            maybeAddResult(iter.pointData())
            iter.next()
        }
    }

    private fun findClosestPointsOptimized() {
        initQueue()
        logger.trace { "Intialized queue: $queue" }
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to remove
            // it before adding any new entries to the queue.
            val entry = queue.poll()
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            val distance = entry.distance
            if (distance >= distanceLimit) {
                queue = PriorityQueue(16)  // Clear any remaining entries.
                break
            }
            var child = entry.id.childBegin()
            // We already know that it has too many points, so process its children.
            // Each child may either be processed directly or enqueued again.  The
            // loop is optimized so that we don't seek unnecessarily.
            var seek = true
            for (i in 0..3) {
                seek = processOrEnqueue(child, iter, seek)
                child = child.next()
            }
        }

    }

    private fun initQueue() {
        requireArgument { queue.isEmpty() }

        // Optimization: rather than starting with the entire index, see if we can
        // limit the search region to a small disc.  Then we can find a covering for
        // that disc and intersect it with the covering for the index.  This can
        // save a lot of work when the search region is small.
        val cap = target.getCapBound()
        if (cap.isEmpty) return  // Empty target.

        if (options.getMaxResult() == 1) {
            // If the user is searching for just the closest point, we can compute an
            // upper bound on search radius by seeking to the center of the target's
            // bounding cap and looking at the adjacent index points (in S2CellId
            // order).  The minimum distance to either of these points is an upper
            // bound on the search radius.
            //
            // TODO(ericv): The same strategy would also work for small values of
            // max_results() > 1, e.g. max_results() == 20, except that we would need to
            // examine more neighbors (at least 20, and preferably 20 in each
            // direction).  It's not clear whether this is a common case, though, and
            // also this would require extending MaybeAddResult() so that it can
            // remove duplicate entries.  (The points added here may be re-added by
            // ProcessOrEnqueue(), but this is okay when max_results() == 1.)
            iter.seek(S2CellId.fromPoint(cap.center))
            if (!iter.done()) {
                maybeAddResult(iter.pointData())
            }
            if (iter.prev()) {
                maybeAddResult(iter.pointData())
            }
            // Skip the rest of the algorithm if we found a matching point.
            if (distanceLimit == distanceFactory.zero()) return
        }
        // We start with a covering of the set of indexed points, then intersect it
        // with the given region (if any) and maximum search radius disc (if any).
        if (indexCovering.isEmpty()) initCovering()
        var initialCells = indexCovering
        val region = options.region
        if (region != null) {
            val coverer = S2RegionCoverer(maxCells = 4)
            coverer.getCovering(region, regionCovering)
            regionCovering.sort()
            S2CellUnion.getIntersection(indexCovering, regionCovering, intersectionWithRegion)
            initialCells = intersectionWithRegion
        }
        if (distanceLimit < distanceFactory.infinity()) {
            val coverer = S2RegionCoverer(maxCells = 4)
            val radius = cap.radius + distanceLimit.getChordAngleBound()
            val searchCap = S2Cap(cap.center, radius)
            coverer.getFastCovering(searchCap, maxDistanceCovering)
            S2CellUnion.getIntersection(initialCells, maxDistanceCovering, intersectionWithMaxDistance)
            initialCells = intersectionWithMaxDistance
        }
        iter.begin()
        var i = 0
        while (i < initialCells.size && !iter.done()) {
            val id = initialCells[i]
            processOrEnqueue(id, iter, id.rangeMin() > iter.id() /*seek*/)
            ++i
        }
    }
    
    private fun initCovering() {
        // Compute the "index covering", which is a small number of S2CellIds that
        // cover the indexed points.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one covering cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a covering cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the covering cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.
        indexCovering.ensureCapacity(6)
        iter.finish()
        if (!iter.prev()) return  // Empty index.
        val indexLastId = iter.id()
        iter.begin()
        if (iter.id() != indexLastId) {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            val level = iter.id().getCommonAncestorLevel(indexLastId) + 1

            // Visit each potential covering cell except the last (handled below).
            val lastId = indexLastId.parent(level)
            var id = iter.id().parent(level)
            while (id != lastId) {
                // Skip any covering cells that don't contain any index cells.
                if (id.rangeMax() < iter.id()) {
                    id = id.next()
                    continue
                }

                // Find the range of index cells contained by this covering cell and
                // then shrink the cell if necessary so that it just covers them.
                val cellFirstId = iter.id()
                iter.seek(id.rangeMax().next())
                iter.prev()
                val cellLastId = iter.id()
                iter.next()
                addInitialRange(cellFirstId, cellLastId)
                id = id.next()
            }
        }
        addInitialRange(iter.id(), indexLastId)
    }

    // Adds a cell to index_covering_ that covers the given inclusive range.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private fun addInitialRange(firstId: S2CellId, lastId: S2CellId) {
        // Add the lowest common ancestor of the given range.
        val level = firstId.getCommonAncestorLevel(lastId)
        requireGE(level, 0)
        indexCovering.add(firstId.parent(level))
    }
    
    private fun maybeAddResult(pointData: PointData<D>) {
        //logger.trace { "Maybe add result ? $pointData" }
        val distance = distanceLimit.clone()
        if (!target.updateMinDistance(pointData.point, distance)) {
            //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> distance ${target.distance(pointData.point)} >= $distanceLimit: skip point" }
            return
        }

        val region = options.region
        if (region != null && !region.contains(pointData.point)) {
            //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> Point is not in the region of the query: skip point" }
            return
        }

        val result = Result(distance, pointData)
        if (options.getMaxResult() == 1) {
            // Optimization for the common case where only the closest point is wanted.
            resultSingleton = result
            distanceLimit = result.distance - options.maxError
            //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> set point as best result. New distance limit = $distanceLimit" }
        } else if (options.getMaxResult() == Options.kMaxMaxResults) {
            resultVector.add(result)  // Sort/unique at end.
            //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> add point to result list." }
        } else {
            // Add this point to result_set_.  Note that with the current algorithm
            // each candidate point is considered at most once (except for one special
            // case where max_results() == 1, see InitQueue for details), so we don't
            // need to worry about possibly adding a duplicate entry here.
            if (resultSet.size >= options.getMaxResult()) {
                val removedResult = resultSet.poll()  // Replace the furthest result point.
                //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> max result reached, remove $removedResult." }
            }
            resultSet.add(result)
            //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> add point to result." }
            if (resultSet.size >= options.getMaxResult()) {
                distanceLimit = resultSet.peek().distance - options.maxError
                //logger.trace { "$uid: Point maybeAddResult(pointData: $pointData) -> New distance limit = $distanceLimit, results=$resultSet" }
            }
        }
    }

    // Either process the contents of the given cell immediately, or add it to the
    // queue to be subdivided.  If "seek" is false, then "iter" must already be
    // positioned at the first indexed point within or after this cell.
    //
    // Returns "true" if the cell was added to the queue, and "false" if it was
    // processed immediately, in which case "iter" is left positioned at the next
    // cell in S2CellId order.
    private fun processOrEnqueue(id: S2CellId, iter: S2PointIndex.Iterator<D>, seek: Boolean): Boolean {
        if (seek) iter.seek(id.rangeMin())
        if (id.isLeaf) {
            // Leaf cells can't be subdivided.
            while (!iter.done() && iter.id() == id) {
                maybeAddResult(iter.pointData())
                iter.next()
            }
            return false  // No need to seek to next child.
        }
        val last = id.rangeMax()
        var numPoints = 0
        while (!iter.done() && iter.id() <= last) {
            if (numPoints == kMinPointsToEnqueue - 1) {
                // This cell has too many points (including this one), so enqueue it.
                val cell = S2Cell(id)
                var distance = distanceLimit.clone()
                // We check "region_" second because it may be relatively expensive.
                val region = options.region
                if (target.updateMinDistance(cell, distance) && (region == null || region.mayIntersect(cell))) {
                    if (useConservativeCellDistance) {
                        // Ensure that "distance" is a lower bound on distance to the cell.
                        distance -= options.maxError
                    }
                    queue.offer(QueueEntry(distance, id))
                }
                return true;  // Seek to next child.
            }
            tmpPointData[numPoints++] = iter.pointData()
            iter.next()
        }
        // There were few enough points that we might as well process them now.
        for (i in 0 until numPoints) {
            maybeAddResult(tmpPointData[i]!!)
        }
        return false  // No need to seek to next child.
    }

    // Options that control the set of points returned.  Note that by default
    // *all* points are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    //
    // This class is also available as S2ClosestPointQueryBase<Data>::Options.
    // (It is defined here to avoid depending on the "Data" template argument.)
    //
    // The Distance template argument is described below.
    open class Options<T : Distance<T>>(

            val distanceFactory: DistanceFactory<T>,

            // Specifies that at most "max_results" points should be returned.
            //
            // REQUIRES: max_results >= 1
            // DEFAULT: numeric_limits<int>::max()
            private var maxResult: Int = kMaxMaxResults,

            // Specifies that only points whose distance to the target is less than
            // "max_distance" should be returned.
            //
            // Note that points whose distance is exactly equal to "max_distance" are
            // not returned.  In most cases this doesn't matter (since distances are
            // not computed exactly in the first place), but if such points are needed
            // then you can retrieve them by specifying "max_distance" as the next
            // largest representable Distance.  For example, if Distance is an
            // S1ChordAngle then you can specify max_distance.Successor().
            //
            // DEFAULT: Distance::Infinity()
            var maxDistance: T = distanceFactory.infinity(),

            // Specifies that points up to max_error() further away than the true
            // closest points may be substituted in the result set, as long as such
            // points satisfy all the remaining search criteria (such as max_distance).
            // This option only has an effect if max_results() is also specified;
            // otherwise all points closer than max_distance() will always be returned.
            //
            // Note that this does not affect how the distance between points is
            // computed; it simply gives the algorithm permission to stop the search
            // early as soon as the best possible improvement drops below max_error().
            //
            // This can be used to implement distance predicates efficiently.  For
            // example, to determine whether the minimum distance is less than D, the
            // IsDistanceLess() method sets max_results() == 1 and max_distance() ==
            // max_error() == D.  This causes the algorithm to terminate as soon as it
            // finds any point whose distance is less than D, rather than continuing to
            // search for a point that is even closer.
            //
            // DEFAULT: Distance::Delta::Zero()
            var maxError: Delta = Delta.zero(),

            // Specifies that points must be contained by the given S2Region.  "region"
            // is owned by the caller and must persist during the lifetime of this
            // object.  The value may be changed between calls to FindClosestPoints(),
            // or reset by calling set_region(nullptr).
            //
            // Note that if you want to set the region to a disc around a target point,
            // it is faster to use a PointTarget with set_max_distance() instead.  You
            // can also call both methods, e.g. to set a maximum distance and also
            // require that points lie within a given rectangle.
            var region: S2Region? = null,

            // Specifies that distances should be computed by examining every point
            // rather than using the S2ShapeIndex.  This is useful for testing,
            // benchmarking, and debugging.
            //
            // DEFAULT: false
            var useBruteForce: Boolean = false
    ) {

        fun setMaxResult(value: Int) {
            requireGE(value, 1)
            this.maxResult = value
        }

        fun getMaxResult(): Int = maxResult
        override fun toString(): String {
            return "Options(maxResult=$maxResult, maxDistance=$maxDistance, maxError=$maxError, region=$region, useBruteForce=$useBruteForce)"
        }


        companion object {

            const val kMaxMaxResults = Int.MAX_VALUE
        }

    }

    // Each "Result" object represents a closest point.
    data class Result<T: Distance<T>, D : Comparable<D>>(val distance: T, val pointData: PointData<D>? = null): Comparable<Result<T, D>> {

        // Returns true if this Result object does not refer to any data point.
        // (The only case where an empty Result is returned is when the
        // FindClosestPoint() method does not find any points that meet the
        // specified criteria.)
        fun isEmpty(): Boolean = pointData == null

        // The point itself.
        fun point(): S2Point = pointData!!.point

        // The client-specified data associated with this point.
        fun data(): D = pointData!!.data

        // Compares two Result objects first by distance, then by point_data().
        override fun compareTo(other: Result<T, D>): Int = ComparisonChain.start()
                .compare(distance, other.distance)
                .compare(pointData, other.pointData) { a, b  ->
                    when {
                        a != null && b != null -> a.compareTo(b)
                        a == null && b != null -> -1
                        a != null && b == null -> 1
                        else -> 0
                    }
                }
                .result()

        override fun toString(): String {
            return "($distance, $pointData)"
        }


    }

    // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    data class QueueEntry<T : Distance<T>>(
            // A lower bound on the distance from the target to "id".  This is the key
            // of the priority queue.
            val distance: T,

            // The cell being queued.
            val id: S2CellId
    ) : Comparable<QueueEntry<T>> {

        // The priority queue returns the largest elements first, so we want the
        // "largest" entry to have the smallest distance.
        override fun compareTo(other: QueueEntry<T>): Int = distance.compareTo(other.distance)

    }
    
    companion object {

        private val logger = KotlinLogging.logger(S2ClosestPointQueryBase::class.java.name)

        // The minimum number of points that a cell must contain to enqueue it
        // rather than processing its contents immediately.
        const val kMinPointsToEnqueue = 13
    }

}
