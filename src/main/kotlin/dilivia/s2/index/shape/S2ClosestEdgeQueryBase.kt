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

import dilivia.PreConditions.checkIsEmpty
import dilivia.PreConditions.requireGE
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.Delta
import dilivia.s2.index.Distance
import dilivia.s2.index.DistanceFactory
import dilivia.s2.index.S2DistanceTarget
import dilivia.s2.index.shape.S2ShapeIndex.CellIterator
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2RegionCoverer
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import java.util.*

/**
 * S2ClosestEdgeQueryBase is a templatized class for finding the closest edge(s) between two geometries.  It is not
 * intended to be used directly, but rather to serve as the implementation of various specialized classes with more
 * convenient APIs (such as S2ClosestEdgeQuery).  It is flexible enough so that it can be adapted to compute maximum
 * distances and even potentially Hausdorff distances.
 *
 * By using the appropriate options, this class can answer questions such as:
 *
 *  - Find the minimum distance between two geometries A and B.
 *  - Find all edges of geometry A that are within a distance D of geometry B.
 *  - Find the k edges of geometry A that are closest to a given point P.
 *
 * You can also specify whether polygons should include their interiors (i.e., if a point is contained by a polygon,
 * should the distance be zero or should it be measured to the polygon boundary?)
 *
 * The input geometries may consist of any number of points, polylines, and polygons (collectively referred to as
 * "shapes").  Shapes do not need to be disjoint; they may overlap or intersect arbitrarily.  The implementation is
 * designed to be fast for both simple and complex geometries.
 *
 * The Distance template argument is used to represent distances. Usually it is a thin wrapper around S1ChordAngle,
 * but another distance type may be used as long as it implements the Distance concept described in S2DistanceTarget.
 * For example this can be used to measure maximum distances, to get more accuracy, or to measure
 * non-spheroidal distances.
 *
 * This class is a port of the S2ClosestEdgeQueryBase class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T The distance type.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2ClosestEdgeQueryBase<T : Distance<T>> {

    /**
     * Options that control the set of edges returned. Note that by default *all* edges are returned, so you will
     * always want to set either the maxResults option or the maxDistance option (or both).
     *
     * @param T The distance type.
     */
    open class Options<T : Distance<T>>(

        /** The distance factory used to create Distance instance. */
        val distanceFactory: DistanceFactory<T>,

        maxResults: Int = kMaxMaxResults,

        /**
         * Specifies that only edges whose distance to the target is less than "maxDistance" should be returned.
         *
         * Note that edges whose distance is exactly equal to "maxDistance" are not returned. In most cases this
         * doesn't matter (since distances are not computed exactly in the first place), but if such edges are needed
         * then you can retrieve them by specifying "maxDistance" as the next largest representable Distance.
         * For example, if Distance is an S1ChordAngle then you can specify max_distance.successor().
         *
         * DEFAULT: Distance.infinity()
         */
        var maxDistance: T = distanceFactory.infinity(),

        /**
         * Specifies that edges up to max_error() further away than the true closest edges may be substituted in the
         * result set, as long as such edges satisfy all the remaining search criteria (such as max_distance).
         * This option only has an effect if max_results() is also specified; otherwise all edges closer than
         * maxDistance will always be returned.
         *
         * Note that this does not affect how the distance between edges is computed; it simply gives the algorithm
         * permission to stop the search early as soon as the best possible improvement drops below maxError.
         *
         * This can be used to implement distance predicates efficiently.  For example, to determine whether the minimum
         * distance is less than D, set maxResults == 1 and maxDistance == maxError == D.  This causes the algorithm to
         * terminate as soon as it finds any edge whose distance is less than D, rather than continuing to search for
         * an edge that is even closer.
         *
         * DEFAULT: Distance.Delta.zero()
         */
        var maxError: Delta = Delta.zero(),

        /**
         * Specifies that polygon interiors should be included when measuring distances. In other words, polygons that
         * contain the target should have a distance of zero.  (For targets consisting of multiple connected
         * components, the distance is zero if any component is contained.)  This is indicated in the results by
         * returning a (shapeId, edgeId) pair with edgeId == -1, i.e. this value denotes the polygons's interior.
         *
         * Note that for efficiency, any polygon that intersects the target may or may not have an (edgeId == -1)
         * result.  Such results are optional because in that case the distance to the polygon is already zero.
         *
         * DEFAULT: true
         */
        var includeInteriors: Boolean = true,

        /**
         * Specifies that distances should be computed by examining every edge rather than using the S2ShapeIndex.
         * This is useful for testing, benchmarking, and debugging.
         *
         * DEFAULT: false
         */
        var useBruteForce: Boolean = false

    ) : Cloneable {

        /**
         * Specifies that at most "maxResults" edges should be returned.
         *
         * REQUIRES: max_results >= 1
         * DEFAULT: kMaxMaxResults
         */
        var maxResults: Int = maxResults
            set(value) {
                requireGE(value, 1)
                field = value
            }

        override fun clone(): Options<T> {
            return Options(
                distanceFactory, maxResults, maxDistance, maxError, includeInteriors, useBruteForce
            )
        }

        companion object {
            const val kMaxMaxResults = Int.MAX_VALUE
        }

    }

    private data class QueueEntry<T : Distance<T>>(
        // A lower bound on the distance from the target to "id".  This is the key
        // of the priority queue.
        val distance: T,

        // The cell being queued.
        val id: S2CellId,

        // If "id" belongs to the index, this field stores the corresponding
        // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
        // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
        // field is to avoid an extra Seek() when the queue entry is processed.
        val indexCell: S2ShapeIndexCell? = null

    ) : Comparable<QueueEntry<T>> {

        // The priority queue returns the smallest elements first
        override fun compareTo(other: QueueEntry<T>): Int = distance.compareTo(other.distance)

    }

    private val uid = UUID.randomUUID().toString()

    /////////////// Query parameters /////////////////

    /** The distance factory. */
    val distanceFactory: DistanceFactory<T>

    /** The index of the query. */
    private lateinit var index: S2ShapeIndex

    /** The query options. */
    private lateinit var options: Options<T>

    /** The query target. */
    private lateinit var target: S2DistanceTarget<T>

    /**
     * True if maxError must be subtracted from priority queue cell distances in order to ensure that such distances
     * are measured conservatively.  This is true only if the target takes advantage of max_error() in order to
     * return faster results, and 0 < maxError < distanceLimit.
     */
    private var useConservativeCellDistance: Boolean = false

    /**
     * For the optimized algorihm we precompute the top-level S2CellIds that will be added to the priority queue.
     * There can be at most 6 of these cells.  Essentially this is just a covering of the indexed edges, except
     * that we also store pointers to the corresponding S2ShapeIndexCells to reduce the number of index seeks required.
     */
    private var indexCovering: ArrayList<S2CellId> = ArrayList(6)
    private var indexCells: ArrayList<S2ShapeIndexCell?> = ArrayList(6)

    // The decision about whether to use the brute force algorithm is based on counting the total number of edges in
    // the index.  However if the index contains a large number of shapes, this in itself might take too long.
    // So instead we only count edges up to (maxBruteForceIndexSize() + 1) for the current target type
    // (stored as indexNumEdgesLimit).
    private var indexNumEdges: Int = 0
    private var indexNumEdgesLimit: Int = 0

    /**
     * The distance beyond which we can safely ignore further candidate edges. (Candidates that are exactly at the limit
     * are ignored; this is more efficient for updateMinDistance() and should not affect clients since distance
     * measurements have a small amount of error anyway.)
     *
     * Initially this is the same as the maximum distance specified by the user, but it can also be updated by the
     * algorithm (see maybeAddResult).
     */
    private lateinit var distanceLimit: T

    /////////////// Results /////////////////

    // The current result set is stored in one of three ways:

    //  - Otherwise results are kept in a btree_set so that we can progressively
    //    reduce the distance limit once max_results() results have been found.
    //    (A priority queue is not sufficient because we need to be able to
    //    check whether a candidate edge is already in the result set.)
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.

    /** If maxResults == 1, keep the best result */
    private lateinit var resultSingleton: S2ClosestEdgeQueryResult<T>

    /** If maxResults() == "infinity", results are appended to this list and sorted/uniqued at the end. */
    private var resultVector: MutableList<S2ClosestEdgeQueryResult<T>> = mutableListOf()

    /**
     * Otherwise results are kept in a TreeSet so that we can progressively reduce the distance limit once maxResults
     * results have been found. (A priority queue is not sufficient because we need to be able to check whether a
     * candidate edge is already in the result set.)
     */
    private var resultSet: TreeSet<S2ClosestEdgeQueryResult<T>> = TreeSet()

    /**
     * When the result edges are stored in a btree_set (see above), usually duplicates can be removed simply by
     * inserting candidate edges in the current set. However this is not true if Options.maxError() > 0 and the Target
     * subtype takes advantage of this by returning suboptimal distances. This is because when updateMinDistance()
     * is called with different "minDist" parameters (i.e., the distance to beat), the implementation may return a
     * different distance for the same edge.  Since the TreeSet is keyed by (distance, shape_id, edge_id) this can
     * create duplicate edges in the results.
     *
     * The flag below is true when duplicates must be avoided explicitly.  This is achieved by maintaining a separate
     * set keyed by (shape_id, edge_id) only, and checking whether each edge is in that set before computing the
     * distance to it.
     *
     * TODO(ericv): Check whether it is faster to avoid duplicates by default (even when Options::max_results() == 1),
     * rather than just when we need to.
     */
    private var avoidDuplicates: Boolean = true

    /** Keep already processed edges. */
    private val testedEdges: MutableSet<ShapeEdgeId> = HashSet()

    /**
     * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted in increasing order of distance from
     * the target.
     */
    private val queue: Queue<QueueEntry<T>> = PriorityQueue(16)

    // Temporaries, defined here to avoid multiple allocations / initializations.
    /** Index cell iterator. */
    private lateinit var iter: CellIterator
    private val maxDistanceCovering: MutableList<S2CellId> = mutableListOf()
    private val initialCells: MutableList<S2CellId> = mutableListOf()

    //////////////////////////////////// Init ////////////////////////////////////////

    /**
     * Default constructor; requires init() to be called.
     *
     * @param distanceFactory The distance factory.
     */
    constructor(distanceFactory: DistanceFactory<T>) {
        this.distanceFactory = distanceFactory
    }

    /**
     * Convenience constructor that calls init().
     *
     * @param distanceFactory The distance factory.
     * @param index The shape index on which the request is done.
     */
    constructor(distanceFactory: DistanceFactory<T>, index: S2ShapeIndex) {
        this.distanceFactory = distanceFactory
        init(index)
    }

    /**
     * Initializes the query.
     *
     * REQUIRES: reInit() must be called if "index" is modified.
     * @param index The index on which the query is done.
     */
    fun init(index: S2ShapeIndex) {
        this.index = index
        reInit()
    }

    /**
     * Reinitializes the query.  This method must be called whenever the underlying index is modified.
     */
    fun reInit() {
        indexNumEdges = 0
        indexNumEdgesLimit = 0
        indexCovering.clear()
        indexCells.clear()
        // We don't initialize iter_ here to make queries on small indexes a bit
        // faster (i.e., where brute force is used).
    }

    /** Gets the underlying S2ShapeIndex. */
    fun index(): S2ShapeIndex = index

    /////////////// Find edges /////////////////

    /**
     * Finds the closest edges to the given target that satisfy the given options. This method may be called multiple
     * times.
     *
     * Note that if options.includeInteriors is true, the result vector may include some entries with edge_id == -1.
     * This indicates that the target intersects the indexed polygon with the given shapeId.
     *
     * @param target The query target.
     * @param options The query options.
     * @return The query result ordered by increase distance.
     */
    fun findClosestEdges(target: S2DistanceTarget<T>, options: Options<T>): List<S2ClosestEdgeQueryResult<T>> {
        val results = mutableListOf<S2ClosestEdgeQueryResult<T>>()
        findClosestEdges(target, options, results)
        return results;
    }

    /**
     * Finds the closest edges to the given target that satisfy the given options. This version can be more efficient
     * when this method is called many times, since it does not require allocating a new list on each call.
     *
     * @param target The query target.
     * @param options The query options.
     * @param results The list to populate with the query result ordered by increase distance.
     */
    fun findClosestEdges(
        target: S2DistanceTarget<T>,
        options: Options<T>,
        results: MutableList<S2ClosestEdgeQueryResult<T>>
    ) {
        findClosestEdgesInternal(target, options)
        results.clear()
        if (options.maxResults == 1) {
            if (resultSingleton.shapeId >= 0) {
                results.add(resultSingleton)
            }
        } else if (options.maxResults == Options.kMaxMaxResults) {
            results.addAll(resultVector)
            results.sortAndRemoveDuplicates()
            resultVector.clear()
        } else {
            results.addAll(resultSet)
            resultSet.clear()
        }

        logger.trace {
            """
            |
            |=====================================================================================
            |findClosestEdge: $uid
            |-------------------------------------------------------------------------------------
            | result: $results
            |=====================================================================================
        """.trimMargin()
        }
    }

    /**
     *  Convenience method that returns exactly one edge. If no edges satisfy the given search criteria, then a
     *  Result with distance == infinity() and shapeId == edgeId == -1 is returned.
     *
     * Note that if options.includeInteriors is true, edge_id == -1 is also used to indicate that the target intersects
     * an indexed polygon (but in that case distance == Zero() and shapeId >= 0).
     *
     * REQUIRES: options.maxResults == 1
     *
     * @param target The query target.
     * @param options The query options.
     * @return The closest edge result.
     */
    fun findClosestEdge(target: S2DistanceTarget<T>, options: Options<T>): S2ClosestEdgeQueryResult<T> {
        check(options.maxResults == 1)
        findClosestEdgesInternal(target, options)
        logger.trace {
            """
            |
            |=====================================================================================
            |findClosestEdge: $uid
            |-------------------------------------------------------------------------------------
            | result: $resultSingleton
            |=====================================================================================
        """.trimMargin()
        }
        return resultSingleton
    }

    ///////////////////////////////// Internal ////////////////////////////

    private fun findClosestEdgesInternal(target: S2DistanceTarget<T>, options: Options<T>) {
        this.target = target
        this.options = options

        testedEdges.clear()
        distanceLimit = options.maxDistance.clone()
        resultSingleton = S2ClosestEdgeQueryResult(distance = distanceFactory.infinity())
        checkIsEmpty(resultVector)
        checkIsEmpty(resultSet)
        //checkGE(target.maxBruteForceIndexSize(), 0) { "Target $target : maxBruteForceIndexSize < 0" }

        logger.trace {
            """
            |
            |=====================================================================================
            |findClosestEdgesInternal: $uid
            |-------------------------------------------------------------------------------------
            | target: $target
            | options: $options
            | distance limit: $distanceLimit
            |=====================================================================================
        """.trimMargin()
        }

        if (distanceLimit == distanceFactory.zero()) return

        if (options.maxResults == Options.kMaxMaxResults && options.maxDistance == distanceFactory.infinity()) {
            logger.warn { "Returning all edges (max_results/max_distance not set)" }
        }

        // If includeInteriors == true. Process all the shapes contained by the target.
        if (options.includeInteriors) {
            val shapeIds = TreeSet<Int>()
            target.visitContainingShapes(index, object : S2DistanceTarget.ShapeVisitor {

                override fun visit(containing_shape: S2Shape, target_point: S2Point): Boolean {
                    shapeIds.add(containing_shape.id)
                    return shapeIds.size < options.maxResults
                }

            })

            logger.trace { "$uid: Contained shape = $shapeIds" }

            for (shape_id in shapeIds) {
                addResult(S2ClosestEdgeQueryResult(distanceFactory.zero(), shape_id, -1));
            }
            if (distanceLimit == distanceFactory.zero()) return
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

        // Use the brute force algorithm if the index is small enough.  To avoid
        // spending too much time counting edges when there are many shapes, we stop
        // counting once there are too many edges.  We may need to recount the edges
        // if we later see a target with a larger brute force edge threshold.
        val minOptimizedEdges = target.maxBruteForceIndexSize() + 1
        if (minOptimizedEdges > indexNumEdgesLimit && indexNumEdges >= indexNumEdgesLimit) {
            indexNumEdges = S2CountEdges.countEdgesUpTo(index, minOptimizedEdges)
            indexNumEdgesLimit = minOptimizedEdges
        }

        if (options.useBruteForce || indexNumEdges < minOptimizedEdges) {
            // The brute force algorithm considers each edge exactly once.
            avoidDuplicates = false
            logger.trace { "$uid: use brute force" }
            findClosestEdgesBruteForce()
        } else {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoidDuplicates = (targetUsesMaxError && options.maxResults > 1)
            logger.trace { "$uid: use optimized" }
            findClosestEdgesOptimized();
        }
    }

    private fun findClosestEdgesBruteForce() {
        for (shape in index) {
            if (shape == null) continue
            val numEdges = shape.numEdges
            for (e in 0 until numEdges) {
                maybeAddResult(shape, e)
            }
        }
    }

    private fun findClosestEdgesOptimized() {
        initQueue();
        // Repeatedly find the closest S2Cell to "target" and either split it into
        // its four children or process all of its edges.
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to
            // remove it before adding any new entries to the queue.
            val entry = queue.poll()
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            val distance = entry.distance
            if (!(distance < distanceLimit)) {
                queue.clear()  // Clear any remaining entries.
                break
            }
            // If this is already known to be an index cell, just process it.
            if (entry.indexCell != null) {
                processEdges(entry)
                continue
            }
            // Otherwise split the cell into its four children.  Before adding a
            // child back to the queue, we first check whether it is empty.  We do
            // this in two seek operations rather than four by seeking to the key
            // between children 0 and 1 and to the key between children 2 and 3.
            val id = entry.id;
            iter.seek(id.child(1).rangeMin())
            if (!iter.done() && iter.id() <= id.child(1).rangeMax()) {
                processOrEnqueue(id.child(1))
            }
            if (iter.prev() && iter.id() >= id.rangeMin()) {
                processOrEnqueue(id.child(0))
            }
            iter.seek(id.child(3).rangeMin())
            if (!iter.done() && iter.id() <= id.rangeMax()) {
                processOrEnqueue(id.child(3));
            }
            if (iter.prev() && iter.id() >= id.child(2).rangeMin()) {
                processOrEnqueue(id.child(2));
            }
        }
    }

    private fun initQueue() {
        check(queue.isEmpty())
        if (indexCovering.isEmpty()) {
            // We delay iterator initialization until now to make queries on very
            // small indexes a bit faster (i.e., where brute force is used).
            iter = index.cellIterator(InitialPosition.UNPOSITIONED)
        }

        // Optimization: if the user is searching for just the closest edge, and the
        // center of the target's bounding cap happens to intersect an index cell,
        // then we try to limit the search region to a small disc by first
        // processing the edges in that cell.  This sets distance_limit_ based on
        // the closest edge in that cell, which we can then use to limit the search
        // area.  This means that the cell containing "target" will be processed
        // twice, but in general this is still faster.
        //
        // TODO(ericv): Even if the cap center is not contained, we could still
        // process one or both of the adjacent index cells in S2CellId order,
        // provided that those cells are closer than distance_limit_.
        val cap = target.getCapBound()
        if (cap.isEmpty) return  // Empty target.
        if (options.maxResults == 1 && iter.locate(cap.center)) {
            processEdges(QueueEntry(distanceFactory.zero(), iter.id(), iter.cell()))
            // Skip the rest of the algorithm if we found an intersecting edge.
            if (distanceLimit == distanceFactory.zero()) return
        }
        if (indexCovering.isEmpty()) initCovering()
        if (distanceLimit == distanceFactory.infinity()) {
            // Start with the precomputed index covering.
            for (i in 0 until indexCovering.size) {
                processOrEnqueue(indexCovering[i], indexCells[i])
            }
        } else {
            // Compute a covering of the search disc and intersect it with the
            // precomputed index covering.
            val coverer = S2RegionCoverer()
            coverer.setMaxCells(4)
            val radius = cap.radius + distanceLimit.getChordAngleBound()
            val searchCap = S2Cap(cap.center, radius)
            coverer.getFastCovering(searchCap, maxDistanceCovering)
            S2CellUnion.getIntersection(indexCovering, maxDistanceCovering, initialCells)

            // Now we need to clean up the initial cells to ensure that they all
            // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
            // the index at all, while other may be descendants of an index cell.)
            var i = 0
            var j = 0
            while (i < initialCells.size) {
                val id_i = initialCells[i]
                // Find the top-level cell that contains this initial cell.
                while (indexCovering[j].rangeMax() < id_i) ++j
                val id_j = indexCovering[j]
                if (id_i == id_j) {
                    // This initial cell is one of the top-level cells.  Use the
                    // precomputed S2ShapeIndexCell pointer to avoid an index seek.
                    processOrEnqueue(id_j, indexCells[j])
                    ++i; ++j
                } else {
                    // This initial cell is a proper descendant of a top-level cell.
                    // Check how it is related to the cells of the S2ShapeIndex.
                    val r = iter.locate(id_i)
                    if (r == CellRelation.INDEXED) {
                        // This cell is a descendant of an index cell.  Enqueue it and skip
                        // any other initial cells that are also descendants of this cell.
                        processOrEnqueue(iter.id(), iter.cell())
                        val last_id = iter.id().rangeMax()
                        while (++i < initialCells.size && initialCells[i] <= last_id)
                            continue;
                    } else {
                        // Enqueue the cell only if it contains at least one index cell.
                        if (r == CellRelation.SUBDIVIDED) processOrEnqueue(id_i, null)
                        ++i;
                    }
                }
            }
        }
    }

    private fun initCovering() {
        // Find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  These are
        // the "top-level" cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one top-level cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a top-level cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the top-level cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.

        // Don't need to reserve index_cells_ since it is an InlinedVector.
        indexCovering.ensureCapacity(6)

        // TODO(ericv): Use a single iterator (iter_) below and save position
        // information using pair<S2CellId, const S2ShapeIndexCell*> type.
        val next = index.newIterator(InitialPosition.BEGIN)
        val last = index.newIterator(InitialPosition.END)
        last.prev()
        if (next.id() != last.id()) {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            val level = next.id().getCommonAncestorLevel(last.id()) + 1

            // Visit each potential top-level cell except the last (handled below).
            val lastId = last.id().parent(level)
            var id = next.id().parent(level)
            while (id != lastId) {
                // Skip any top-level cells that don't contain any index cells.
                if (id.rangeMax() < next.id()) {
                    id = id.next()
                    continue
                }

                // Find the range of index cells contained by this top-level cell and
                // then shrink the cell if necessary so that it just covers them.
                val cell_first = next.id() to next.cell()
                next.seek(id.rangeMax().next())
                next.prev()
                val cell_last = next.id() to next.cell()
                next.next()
                addInitialRange(cell_first, cell_last);
                id = id.next()
            }
        }
        addInitialRange(next.id() to next.cell(), last.id() to last.cell());
/*

        // Find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  These are
        // the "top-level" cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one top-level cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a top-level cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the top-level cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.

        // Don't need to reserve index_cells_ since it is an InlinedVector.
        indexCovering.ensureCapacity(6)

        // TODO(ericv): Use a single iterator (iter_) below and save position
        // information using pair<S2CellId, const S2ShapeIndexCell*> type.
        val next = index.cellIterator(InitialPosition.BEGIN)
        val last = index.cellIterator(InitialPosition.END)
        last.prev();
        if (next.id() != last.id()) {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            val level = next.id().getCommonAncestorLevel(last.id()) + 1

            // Visit each potential top-level cell except the last (handled below).
            val last_id = last.id().parent(level)
            var id = next.id().parent(level)
            while (id != last_id) {
                println("Init covering test cell: $id")
                // Skip any top-level cells that don't contain any index cells.
                if (id.rangeMax() < next.id()) continue

                // Find the range of index cells contained by this top-level cell and
                // then shrink the cell if necessary so that it just covers them.
                val cell_first = next.clone()
                next.seek(id.rangeMax().next())
                val cell_last = next.clone()
                cell_last.prev()
                addInitialRange(cell_first, cell_last);
                id = id.next()
            }
        }
        addInitialRange(next, last);

 */
    }

    // Add an entry to index_covering_ and index_cells_ that covers the given
    // inclusive range of cells.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private fun addInitialRange(first: Pair<S2CellId, S2ShapeIndexCell>, last: Pair<S2CellId, S2ShapeIndexCell>) {
        if (first.first == last.first) {
            // The range consists of a single index cell.
            indexCovering.add(first.first)
            indexCells.add(first.second)
        } else {
            // Add the lowest common ancestor of the given range.
            val level = first.first.getCommonAncestorLevel(last.first)
            check(level >= 0)
            indexCovering.add(first.first.parent(level))
            indexCells.add(null)
        }
    }

    private fun maybeAddResult(shape: S2Shape, edge_id: Int) {
        if (avoidDuplicates && !testedEdges.add(ShapeEdgeId(shape.id, edge_id))) {
            logger.trace { "$uid: maybeAddResult(shape: ${shape.id}, edge_id: $edge_id) = false, edge already tested." }
            return
        }
        val edge = shape.edge(edge_id)
        val distance = distanceLimit.clone()
        if (target.updateMinDistance(edge.v0, edge.v1, distance)) {
            addResult(S2ClosestEdgeQueryResult(distance, shape.id, edge_id))
        } else {
            logger.trace { "$uid: maybeAddResult(shape: ${shape.id}, edge_id: $edge_id) = false, edge is further that $distanceLimit." }
        }
    }

    private fun addResult(result: S2ClosestEdgeQueryResult<T>) {
        if (options.maxResults == 1) {
            // Optimization for the common case where only the closest edge is wanted.
            resultSingleton = result;
            distanceLimit = result.distance - options.maxError
            logger.trace { "$uid: set closest result = $result, distance limit = $distanceLimit" }
        } else if (options.maxResults == Options.kMaxMaxResults) {
            logger.trace { "$uid: add result = $result" }
            resultVector.add(result);  // Sort/unique at end.
        } else {
            // Add this edge to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            logger.trace { "$uid: add result = $result" }
            resultSet.add(result)
            val size = resultSet.size
            if (size >= options.maxResults) {
                if (size > options.maxResults) {
                    val removedResult = resultSet.pollLast()
                    logger.trace { "$uid: max result reached, remove further result = $removedResult" }
                }
                distanceLimit = resultSet.last().distance - options.maxError
                logger.trace { "$uid: max result reached, distance limit = $distanceLimit" }
            }
        }
    }

    private fun processEdges(entry: QueueEntry<T>) {
        val indexCell = entry.indexCell!!
        for (s in 0 until indexCell.numClipped) {
            val clipped = indexCell.clipped(s)
            val shape = index.shape(clipped.shapeId)
            if (shape != null) {
                for (j in 0 until clipped.numEdges) {
                    maybeAddResult(shape, clipped.edge(j))
                }
            }
        }
    }

    // Enqueue the given cell id.
    // REQUIRES: iter_ is positioned at a cell contained by "id".
    private fun processOrEnqueue(id: S2CellId) {
        check(id.contains(iter.id()))
        if (iter.id() == id) {
            processOrEnqueue(id, iter.cell())
        } else {
            processOrEnqueue(id, null)
        }
    }

    // Add the given cell id to the queue.  "index_cell" is the corresponding
    // S2ShapeIndexCell, or nullptr if "id" is not an index cell.
    //
    // This version is called directly only by InitQueue().
    private fun processOrEnqueue(id: S2CellId, index_cell: S2ShapeIndexCell?) {
        if (index_cell != null) {
            // If this index cell has only a few edges, then it is faster to check
            // them directly rather than computing the minimum distance to the S2Cell
            // and inserting it into the queue.
            val kMinEdgesToEnqueue = 10
            val num_edges = countEdges(index_cell)
            if (num_edges == 0) return
            if (num_edges < kMinEdgesToEnqueue) {
                // Set "distance" to zero to avoid the expense of computing it.
                processEdges(QueueEntry(distanceFactory.zero(), id, index_cell))
                return
            }
        }
        // Otherwise compute the minimum distance to any point in the cell and add
        // it to the priority queue.
        val cell = S2Cell(id)
        var distance = distanceLimit.clone()
        if (!target.updateMinDistance(cell, distance)) return
        if (useConservativeCellDistance) {
            // Ensure that "distance" is a lower bound on the true distance to the cell.
            distance -= options.maxError;  // operator-=() not defined.
        }
        queue.offer(QueueEntry(distance, id, index_cell))
    }

    companion object {
        private val logger = KotlinLogging.logger(S2ClosestEdgeQueryBase::class.java.name)

        // Return the number of edges in the given index cell.
        private fun countEdges(cell: S2ShapeIndexCell): Int {
            var count = 0;
            for (s in 0 until cell.numClipped) {
                count += cell.clipped(s).numEdges
            }
            return count
        }
    }


}

