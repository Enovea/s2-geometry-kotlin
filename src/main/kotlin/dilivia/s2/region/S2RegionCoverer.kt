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

import dilivia.PreConditions.checkLE
import dilivia.PreConditions.checkNE
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLE
import dilivia.collections.isSorted
import dilivia.collections.lowerBound
import dilivia.collections.upperBound
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import java.util.*

// An S2RegionCoverer is a class that allows arbitrary regions to be
// approximated as unions of cells (S2CellUnion).  This is useful for
// implementing various sorts of search and precomputation operations.
//
// Typical usage:
//
// S2RegionCoverer::Options options
// options.set_max_cells(5)
// S2RegionCoverer coverer(options)
// S2Cap cap(center, radius)
// S2CellUnion covering = coverer.GetCovering(cap)
//
// This yields a vector of at most 5 cells that is guaranteed to cover the
// given cap (a disc-shaped region on the sphere).
//
// The approximation algorithm is not optimal but does a pretty good job in
// practice.  The output does not always use the maximum number of cells
// allowed, both because this would not always yield a better approximation,
// and because max_cells() is a limit on how much work is done exploring the
// possible covering as well as a limit on the final output size.
//
// Because it is an approximation algorithm, one should not rely on the
// stability of the output.  In particular, the output of the covering algorithm
// may change across different versions of the library.
//
// One can also generate interior coverings, which are sets of cells which
// are entirely contained within a region.  Interior coverings can be
// empty, even for non-empty regions, if there are no cells that satisfy
// the provided constraints and are contained by the region.  Note that for
// performance reasons, it is wise to specify a max_level when computing
// interior coverings - otherwise for regions with small or zero area, the
// algorithm may spend a lot of time subdividing cells all the way to leaf
// level to try to find contained cells.
class S2RegionCoverer(
        maxCells: Int = kDefaultMaxCells,
        minLevel: Int = 0,
        maxLevel: Int = S2CellId.kMaxLevel,
        levelMod: Int = 1
) {

    private val logger = KotlinLogging.logger { }

    // Sets the desired maximum number of cells in the approximation.  Note
    // the following:
    //
    //  - For any setting of max_cells(), up to 6 cells may be returned if
    //    that is the minimum number required (e.g. if the region intersects
    //    all six cube faces).  Even for very tiny regions, up to 3 cells may
    //    be returned if they happen to be located at the intersection of
    //    three cube faces.
    //
    //  - min_level() takes priority over max_cells(), i.e. cells below the
    //    given level will never be used even if this causes a large number of
    //    cells to be returned.
    //
    //  - If max_cells() is less than 4, the area of the covering may be
    //    arbitrarily large compared to the area of the original region even
    //    if the region is convex (e.g. an S2Cap or S2LatLngRect).
    //
    // Accuracy is measured by dividing the area of the covering by the area
    // of the original region.  The following table shows the median and worst
    // case values for this area ratio on a test case consisting of 100,000
    // spherical caps of random size (generated using s2region_coverer_test):
    //
    //   max_cells:        3      4     5     6     8    12    20   100   1000
    //   median ratio:  5.33   3.32  2.73  2.34  1.98  1.66  1.42  1.11   1.01
    //   worst case:  215518  14.41  9.72  5.26  3.91  2.75  1.92  1.20   1.02
    //
    // The default value of 8 gives a reasonable tradeoff between the number
    // of cells used and the accuracy of the approximation.
    //
    // DEFAULT: kDefaultMaxCells
    private var maxCells: Int = kDefaultMaxCells
    fun setMaxCells(value: Int) {
        this.maxCells = value
    }

    fun getMaxCells(): Int = maxCells

    // Sets the minimum and maximum cell levels to be used.  The default is to
    // use all cell levels.
    //
    // To find the cell level corresponding to a given physical distance, use
    // the S2Cell metrics defined in s2metrics.h.  For example, to find the
    // cell level that corresponds to an average edge length of 10km, use:
    //
    //   int level =
    //       S2::kAvgEdge.GetClosestLevel(S2Earth::KmToRadians(length_km))
    //
    // Note that min_level() takes priority over max_cells(), i.e. cells below
    // the given level will never be used even if this causes a large number
    // of cells to be returned.  (This doesn't apply to interior coverings,
    // since interior coverings make no completeness guarantees -- the result
    // is simply a set of cells that covers as much of the interior as
    // possible while satisfying the given restrictions.)
    //
    // REQUIRES: min_level() <= max_level()
    // DEFAULT: 0
    private var minLevel: Int = 0
    fun setMinLevel(value: Int) {
        requireGE(value, 0)
        requireLE(value, S2CellId.kMaxLevel)
        minLevel = max(0, min(S2CellId.kMaxLevel, value))
    }

    fun getMinLevel(): Int = minLevel

    // REQUIRES: min_level() <= max_level()
    // DEFAULT: S2CellId::kMaxLevel
    private var maxLevel: Int = S2CellId.kMaxLevel
    fun setMaxLevel(value: Int) {
        requireGE(value, 0)
        requireLE(value, S2CellId.kMaxLevel)
        maxLevel = max(0, min(S2CellId.kMaxLevel, value))
    }

    fun getMaxLevel(): Int = maxLevel

    // Convenience function that sets both the maximum and minimum cell levels.
    fun setFixedLevel(level: Int) {
        setMinLevel(level)
        setMaxCells(level)
    }

    // If specified, then only cells where (level - min_level) is a multiple
    // of "level_mod" will be used (default 1).  This effectively allows the
    // branching factor of the S2CellId hierarchy to be increased.  Currently
    // the only parameter values allowed are 1, 2, or 3, corresponding to
    // branching factors of 4, 16, and 64 respectively.
    //
    // DEFAULT: 1
    private var levelMod: Int = 1
    fun setLevelMod(value: Int) {
        requireGE(value, 1)
        requireLE(value, 3)
        levelMod = max(1, min(3, value))
    }

    fun getLevelMod(): Int = levelMod

    init {
        checkLE(minLevel, maxLevel)
        setMaxCells(maxCells)
        setMinLevel(minLevel)
        setMaxLevel(maxLevel)
        setLevelMod(levelMod)
    }

    // We save a temporary copy of the pointer passed to GetCovering() in order
    // to avoid passing this parameter around internally.  It is only used (and
    // only valid) for the duration of a single GetCovering() call.
    private var region: S2Region? = null

    // The set of S2CellIds that have been added to the covering so far.
    private val result: MutableList<S2CellId> = mutableListOf()

    // We keep the candidates in a priority queue.  We specify a vector to hold
    // the queue entries since for some reason priority_queue<> uses a deque by
    // default.  We define our own own comparison function on QueueEntries in
    // order to make the results deterministic.  (Using the default
    // less<QueueEntry>, entries of equal priority would be sorted according to
    // the memory address of the candidate.)
    private val candidateQueue: PriorityQueue<S2RegionCovererQueueEntry> = PriorityQueue(10)

    // True if we're computing an interior covering.
    private var interiorCovering: Boolean = false

    // Counter of number of candidates created, for performance evaluation.
    private var candidatesCreatedCounter: Int = 0

    /**
     *
     * @property cell
     * @property isTerminal Cell should not be expanded further.
     * @property numChildren Number of children that intersect the region.
     * @property children Actual size may be 0, 4, 16, or 64 elements.
     */
    data class Candidate(val cell: S2Cell, var isTerminal: Boolean, var numChildren: Int, val children: MutableList<Candidate>) {

        constructor(cell: S2Cell, maxChildren: Int) : this(cell, maxChildren == 0, 0, mutableListOf())

    }

    // Convenience function that returns the maximum level such that
    //
    //   (level <= max_level()) && (level - min_level()) % level_mod() == 0.
    //
    // This is the maximum level that will actually be used in coverings.
    fun trueMaxLevel(): Int {
        if (levelMod == 1) return maxLevel
        return maxLevel - (maxLevel - minLevel) % levelMod
    }

    // Returns an S2CellUnion that covers (GetCovering) or is contained within
    // (GetInteriorCovering) the given region and satisfies the current options.
    //
    // Note that if options().min_level() > 0 or options().level_mod() > 1, the
    // by definition the S2CellUnion may not be normalized, i.e. there may be
    // groups of four child cells that can be replaced by their parent cell.
    fun getCovering(region: S2Region): S2CellUnion {
        interiorCovering = false
        getCoveringInternal(region)
        val union = S2CellUnion.fromVerbatim(result)
        result.clear()
        return union
    }

    fun getInteriorCovering(region: S2Region): S2CellUnion {
        interiorCovering = true
        getCoveringInternal(region)
        val union = S2CellUnion.fromVerbatim(result)
        result.clear()
        return union
    }

    // Like the methods above, but works directly with a vector of S2CellIds.
    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    fun getCovering(region: S2Region, covering: MutableList<S2CellId>) {
        interiorCovering = false
        getCoveringInternal(region)
        covering.addAll(result)
        result.clear()
    }

    fun getInteriorCovering(region: S2Region, covering: MutableList<S2CellId>) {
        interiorCovering = true
        getCoveringInternal(region)
        covering.addAll(result)
        result.clear()
    }

    // Like GetCovering(), except that this method is much faster and the
    // coverings are not as tight.  All of the usual parameters are respected
    // (max_cells, min_level, max_level, and level_mod), except that the
    // implementation makes no attempt to take advantage of large values of
    // max_cells().  (A small number of cells will always be returned.)
    //
    // This function is useful as a starting point for algorithms that
    // recursively subdivide cells.
    fun getFastCovering(region: S2Region, covering: MutableList<S2CellId>) {
        region.getCellUnionBound(covering)
        logger.trace { "getCellUnionBound(region = $region) = $covering" }
        canonicalizeCovering(covering)
        logger.trace { "getFastCovering(region = $region) = $covering" }
    }

    // Returns true if the given S2CellId vector represents a valid covering
    // that conforms to the current covering parameters.  In particular:
    //
    //  - All S2CellIds must be valid.
    //
    //  - S2CellIds must be sorted and non-overlapping.
    //
    //  - S2CellId levels must satisfy min_level(), max_level(), and level_mod().
    //
    //  - If covering.size() > max_cells(), there must be no two cells with
    //    a common ancestor at min_level() or higher.
    //
    //  - There must be no sequence of cells that could be replaced by an
    //    ancestor (i.e. with level_mod() == 1, the 4 child cells of a parent).
    fun isCanonical(covering: S2CellUnion): Boolean = isCanonical(covering.cellIds())
    fun isCanonical(covering: List<S2CellId>): Boolean {
        // We check this on each call because of mutable_options().
        requireLE(minLevel, maxLevel)

        val tooManyCells = covering.size > maxCells
        var sameParentCount = 1
        var prevId = S2CellId.none
        for (id in covering) {
            if (!id.isValid) {
                logger.trace { "Cell $id is not valid. Covering is not canonical." }
                return false
            }

            // Check that the S2CellId level is acceptable.
            val level = id.level()
            if (level < minLevel || level > maxLevel) {
                logger.debug { "Cell $id: level $level is outside the range $minLevel .. $maxLevel. Covering is not canonical." }
                return false
            }
            if (levelMod > 1 && (level - minLevel) % levelMod != 0) {
                logger.debug { "Cell $id: levelMod (=$levelMod) > 1 and (level - minLevel) % levelMod = ${(level - minLevel) % levelMod } != 0. Covering is not canonical." }
                return false
            }

            if (prevId != S2CellId.none) {
                // Check that cells are sorted and non-overlapping.
                if (prevId.rangeMax() >= id.rangeMin()) {
                    logger.debug { "Cell $id; previous = $prevId: prevId.rangeMax() (=${prevId.rangeMax()}) >= id.rangeMin(=${id.rangeMin()}). Covering is not canonical." }
                    return false
                }

                // If there are too many cells, check that no pair of adjacent cells
                // could be replaced by an ancestor.
                if (tooManyCells) {
                    val commonAncestorLevel = id.getCommonAncestorLevel(prevId)
                    if (commonAncestorLevel >= minLevel) {
                        logger.debug { "Has too many cells and id(=$id).getCommonAncestorLevel(prevId(=$prevId)) = $commonAncestorLevel >= minLevel (=$minLevel)" }
                        return false
                    }
                }

                // Check that there are no sequences of (4 ** level_mod) cells that all
                // have the same parent (considering only multiples of "level_mod").
                val plevel = level - levelMod
                if (plevel < minLevel || level != prevId.level() || id.parent(plevel) != prevId.parent(plevel)) {
                    sameParentCount = 1
                } else if (++sameParentCount == (1 shl (2 * levelMod))) {
                    return false
                }
            }
            prevId = id
        }
        return true
    }

    // Modify "covering" if necessary so that it conforms to the current
    // covering parameters (max_cells, min_level, max_level, and level_mod).
    // There are no restrictions on the input S2CellIds (they may be unsorted,
    // overlapping, etc).
    fun canonicalizeCovering(covering: S2CellUnion): S2CellUnion {
        val ids = covering.cellIds().toMutableList()
        canonicalizeCovering(ids)
        return S2CellUnion(ids)
    }

    fun canonicalizeCovering(covering: MutableList<S2CellId>) {
        logger.trace { """
            |
            |============================================================
            | canonicalizeCovering($covering)
            |------------------------------------------------------------
        """.trimMargin() }
        // We check this on each call because of mutable_options().
        checkLE(minLevel, maxLevel)

        // Note that when the covering parameters have their default values, almost
        // all of the code in this function is skipped.

        // If any cells are too small, or don't satisfy level_mod(), then replace
        // them with ancestors.
        if (maxLevel < S2CellId.kMaxLevel || levelMod > 1) {
            covering.forEachIndexed { index, id ->
                val level = id.level()
                val newLevel = adjustLevel(min(level, maxLevel))
                if (newLevel != level) {
                    logger.trace { "Cell $id is too small, replace by parent at level $newLevel: ${id.parent(newLevel)}" }
                    covering[index] = id.parent(newLevel)
                }
            }
        }

        // Sort the cells and simplify them.
        S2CellUnion.normalize(covering)
        logger.trace { "Normalized order: $covering" }

        // Make sure that the covering satisfies min_level() and level_mod(),
        // possibly at the expense of satisfying max_cells().
        if (minLevel > 0 || levelMod > 1) {
            S2CellUnion.denormalize(covering, minLevel, levelMod, result)
            covering.clear()
            covering.addAll(result)
            result.clear()
        }

        // If there are too many cells and the covering is very large, use the
        // S2RegionCoverer to compute a new covering.  (This avoids possible O(n^2)
        // behavior of the simpler algorithm below.)
        val excess = covering.size - maxCells
        if (excess <= 0 || isCanonical(covering)) {
            return
        }
        if (excess * covering.size > 10000) {
            getCovering(S2CellUnion(covering), covering)
        } else {
            // Repeatedly replace two adjacent cells in S2CellId order by their lowest
            // common ancestor until the number of cells is acceptable.
            while (covering.size > maxCells) {
                logger.trace { "Covering size = ${covering.size} > $maxCells (covering = $covering)" }
                var bestIndex = -1
                var bestLevel = -1
                var i = 0
                while (i + 1 < covering.size) {
                    var level = covering[i].getCommonAncestorLevel(covering[i + 1])
                    logger.trace { "i = $i : Common ancestor level of covering[i] = ${covering[i]} and covering[i+1] = ${covering[i+1]} => $level" }
                    level = adjustLevel(level)
                    logger.trace { "i = $i : Adjusted level = $level" }
                    if (level > bestLevel) {
                        bestLevel = level
                        bestIndex = i
                        logger.trace { "i = $i : level $level is better => bestLevel = $level, bestIndex = $i" }
                    }
                    ++i
                }

                logger.trace { "BestLevel = $bestLevel, $bestIndex = $bestIndex" }

                if (bestLevel < minLevel) break

                // Replace all cells contained by the new ancestor cell.
                var id = covering[bestIndex].parent(bestLevel)
                replaceCellsWithAncestor(covering, id)

                // Now repeatedly check whether all children of the parent cell are
                // present, in which case we can replace those cells with their parent.
                while (bestLevel > minLevel) {
                    bestLevel -= levelMod
                    id = id.parent(bestLevel)
                    if (!containsAllChildren(covering, id)) break
                    replaceCellsWithAncestor(covering, id)
                }
            }
        }

        logger.trace { """
            |
            |------------------------------------------------------------
            | canonicalizeCovering => $covering
            |============================================================
        """.trimMargin() }
        checkState { isCanonical(covering) }
    }

    // If the cell intersects the given region, return a new candidate with no
    // children, otherwise return nullptr.  Also marks the candidate as "terminal"
    // if it should not be expanded further.
    private fun newCandidate(cell: S2Cell): Candidate? {
        val region = region
        check(region != null)
        if (!region.mayIntersect(cell)) return null

        var isTerminal = false
        if (cell.level() >= minLevel) {
            if (interiorCovering) {
                if (region.contains(cell)) {
                    isTerminal = true
                } else if (cell.level() + levelMod > maxLevel) {
                    return null
                }
            } else {
                if (cell.level() + levelMod > maxLevel || region.contains(cell)) {
                    isTerminal = true
                }
            }
        }
        ++candidatesCreatedCounter
        val maxChildren = if (isTerminal) 0 else 1 shl maxChildrenShift()
        return Candidate(cell, maxChildren)
    }

    // Returns the log base 2 of the maximum number of children of a candidate.
    private fun maxChildrenShift(): Int = 2 * levelMod

    // Processes a candidate by either adding it to the result_ vector or
    // expanding its children and inserting it into the priority queue.
    // Passing an argument of nullptr does nothing.
    private fun addCandidate(candidate: Candidate?) {
        if (candidate == null) return

        if (candidate.isTerminal) {
            result.add(candidate.cell.id())
            return
        }
        requireEQ(0, candidate.numChildren)

        // Expand one level at a time until we hit min_level() to ensure that we
        // don't skip over it.
        val numLevels = (if (candidate.cell.level() < minLevel) 1 else levelMod)
        val numTerminals = expandChildren(candidate, candidate.cell, numLevels)

        if (candidate.numChildren == 0) {
            // do nothing
        } else if (!interiorCovering &&
                numTerminals == 1 shl maxChildrenShift() &&
                candidate.cell.level() >= minLevel) {
            // Optimization: add the parent cell rather than all of its children.
            // We can't do this for interior coverings, since the children just
            // intersect the region, but may not be contained by it - we need to
            // subdivide them further.
            candidate.isTerminal = true
            addCandidate(candidate)
        } else {
            // We negate the priority so that smaller absolute priorities are returned
            // first.  The heuristic is designed to refine the largest cells first,
            // since those are where we have the largest potential gain.  Among cells
            // of the same size, we prefer the cells with the fewest children.
            // Finally, among cells with equal numbers of children we prefer those
            // with the smallest number of children that cannot be refined further.
            val priority = -((((candidate.cell.level() shl maxChildrenShift()) + candidate.numChildren) shl maxChildrenShift()) + numTerminals)
            candidateQueue.add(S2RegionCovererQueueEntry(priority, candidate))
            logger.trace { "Push: ${candidate.cell.id()} ($priority)" }
        }
    }

    // Populates the children of "candidate" by expanding the given number of
    // levels from the given cell.  Returns the number of children that were
    // marked "terminal".
    private fun expandChildren(candidate: Candidate, cell: S2Cell, num_levels: Int): Int {
        val region = region
        check(region != null)
        val levels = num_levels - 1
        val childCells = Array(4) { S2Cell() }
        cell.subdivide(childCells)
        var numTerminals = 0
        for (i in 0..3) {
            if (levels > 0) {
                if (region.mayIntersect(childCells[i])) {
                    numTerminals += expandChildren(candidate, childCells[i], levels)
                }
                continue
            }
            val child = newCandidate(childCells[i])
            if (child != null) {
                candidate.children.add(child)
                candidate.numChildren++
                if (child.isTerminal) ++numTerminals
            }
        }
        return numTerminals
    }

    // Computes a set of initial candidates that cover the given region.
    private fun getInitialCandidates() {
        val region = region
        check(region != null)
        // Optimization: start with a small (usually 4 cell) covering of the
        // region's bounding cap.
        val tmpCoverer = S2RegionCoverer(
                maxCells = min(4, maxCells),
                maxLevel = maxLevel
        )
        val cells = mutableListOf<S2CellId>()
        tmpCoverer.getFastCovering(region, cells)
        adjustCellLevels(cells)
        for (cellId in cells) {
            addCandidate(newCandidate(S2Cell(cellId)))
        }
        logger.debug { "Initial candidate = ${candidateQueue.map { it.candidate.cell.id() }}" }
    }

    // Generates a covering and stores it in result_.
    private fun getCoveringInternal(region: S2Region) {
        // We check this on each call because of mutable_options().
        checkLE(minLevel, maxLevel)

        // Strategy: Start with the 6 faces of the cube.  Discard any
        // that do not intersect the shape.  Then repeatedly choose the
        // largest cell that intersects the shape and subdivide it.
        //
        // result_ contains the cells that will be part of the output, while pq_
        // contains cells that we may still subdivide further.  Cells that are
        // entirely contained within the region are immediately added to the output,
        // while cells that do not intersect the region are immediately discarded.
        // Therefore pq_ only contains cells that partially intersect the region.
        // Candidates are prioritized first according to cell size (larger cells
        // first), then by the number of intersecting children they have (fewest
        // children first), and then by the number of fully contained children
        // (fewest children first).

        checkState { candidateQueue.isEmpty() }
        checkState { result.isEmpty() }
        this.region = region
        candidatesCreatedCounter = 0

        getInitialCandidates()
        while (!candidateQueue.isEmpty() && (!interiorCovering || result.size < maxCells)) {
            val candidate = candidateQueue.poll().candidate
            // For interior coverings we keep subdividing no matter how many children
            // the candidate has.  If we reach max_cells() before expanding all
            // children, we will just use some of them.  For exterior coverings we
            // cannot do this, because the result has to cover the whole region, so
            // all children have to be used.  The (candidate->num_children == 1) case
            // takes care of the situation when we already have more than max_cells()
            // in results (min_level is too high).  Subdividing the candidate with one
            // child does no harm in this case.
            if (interiorCovering ||
                    candidate.cell.level() < minLevel ||
                    candidate.numChildren == 1 ||
                    (result.size + candidateQueue.size + candidate.numChildren <= maxCells)) {
                // Expand this candidate into its children.
                for (i in 0 until candidate.numChildren) {
                    if (interiorCovering && result.size >= maxCells) {
                        // Do nothing
                    } else {
                        addCandidate(candidate.children[i])
                    }
                }
            } else {
                candidate.isTerminal = true
                addCandidate(candidate)
            }
        }
        logger.trace { "Created ${result.size} cells, $candidatesCreatedCounter candidates created, ${candidateQueue.size} left" }
        candidateQueue.clear()
        this.region = null

        // Rather than just returning the raw list of cell ids, we construct a cell
        // union and then denormalize it.  This has the effect of replacing four
        // child cells with their parent whenever this does not violate the covering
        // parameters specified (min_level, level_mod, etc).  This significantly
        // reduces the number of cells returned in many cases, and it is cheap
        // compared to computing the covering in the first place.
        S2CellUnion.normalize(result)
        if (minLevel > 0 || levelMod > 1) {
            val resultCopy = mutableListOf<S2CellId>()
            resultCopy.addAll(result)
            S2CellUnion.denormalize(resultCopy, minLevel, levelMod, result)
        }
        checkState({ isCanonical(result) }, { "getCoveringInternal(region = $region) is not canonical: $result" })
    }

    // If level > min_level(), then reduces "level" if necessary so that it also
    // satisfies level_mod().  Levels smaller than min_level() are not affected
    // (since cells at these levels are eventually expanded).
    private fun adjustLevel(level: Int): Int {
        return if (levelMod > 1 && level > minLevel) {
            level - (level - minLevel) % levelMod
        } else level
    }

    // Ensures that all cells with level > min_level() also satisfy level_mod(),
    // by replacing them with an ancestor if necessary.  Cell levels smaller
    // than min_level() are not modified (see AdjustLevel).  The output is
    // then normalized to ensure that no redundant cells are present.
    private fun adjustCellLevels(cells: MutableList<S2CellId>) {
        requireArgument { cells.isSorted() }
        if (levelMod == 1) return

        var out = 0
        for (id in cells) {
            var currentId = id
            val level = currentId.level()
            val newLevel = adjustLevel(level)
            if (newLevel != level) currentId = currentId.parent(newLevel)
            if (out > 0 && cells[out - 1].contains(currentId)) continue
            while (out > 0 && currentId.contains(cells[out - 1])) --out
            cells[out++] = currentId
        }
        while (cells.size > out) cells.removeLast()
    }

    // Returns true if "covering" contains all children of "id" at level
    // (id.level() + options_.level_mod()).
    private fun containsAllChildren(covering: List<S2CellId>, id: S2CellId): Boolean {
        var it = covering.indexOfFirst { cellId -> cellId >= id.rangeMin() }
        if (it == -1) it = covering.size
        val level = id.level() + levelMod
        var child = id.childBegin(level)
        while (child != id.childEnd(level)) {
            logger.trace { "it = $it, level = $level, child = $child " }
            if (it > covering.lastIndex || covering[it] != child) {
                logger.trace { "containsAllChildren(id = $id, covering = $covering) = false" }
                return false
            }
            ++it
            child = child.next()
        }
        logger.trace { "containsAllChildren(id = $id, covering = $covering) = true" }
        return true
    }

    // Replaces all descendants of "id" in "covering" with "id".
    // REQUIRES: "covering" contains at least one descendant of "id".
    private fun replaceCellsWithAncestor(covering: MutableList<S2CellId>, id: S2CellId) {
        val begin = covering.lowerBound(0, covering.size, id.rangeMin())
        val end = covering.upperBound(0, covering.size, id.rangeMax())
        checkNE(begin, end)
        logger.trace { """
            | ReplaceCellsWithAncestor(id = $id, covering = $covering)
            | --------------------------------------------------------------------
            | begin = $begin (${if (begin in covering.indices) covering[begin] else null})
            | end = $end (${if (end in covering.indices) covering[end] else null})
        """.trimMargin() }
        repeat(end - begin - 1) { covering.removeAt(begin + 1) }
        covering[begin] = id
        logger.trace { "ReplaceCellsWithAncestor => $covering" }
    }

    companion object {
        const val kDefaultMaxCells = 8

        // Given a connected region and a starting point on the boundary or inside the
        // region, returns a set of cells at the given level that cover the region.
        // The output cells are returned in arbitrary order.
        //
        // Note that this method is *not* faster than the regular GetCovering()
        // method for most region types, such as S2Cap or S2Polygon, and in fact it
        // can be much slower when the output consists of a large number of cells.
        // Currently it can be faster at generating coverings of long narrow regions
        // such as polylines, but this may change in the future, in which case this
        // method will most likely be removed.
        fun getSimpleCovering(region: S2Region, start: S2Point, level: Int, output: MutableList<S2CellId>): Unit = floodFill(region, S2CellId.fromPoint(start).parent(level), output)

        // Like GetSimpleCovering(), but accepts a starting S2CellId rather than a
        // starting point and cell level.  Returns all edge-connected cells at the
        // same level as "start" that intersect "region", in arbitrary order.
        fun floodFill(region: S2Region, start: S2CellId, output: MutableList<S2CellId>) {
            val all = mutableSetOf<S2CellId>()
            val frontier = mutableListOf<S2CellId>()
            output.clear()
            all.add(start)
            frontier.add(start)
            while (frontier.isNotEmpty()) {
                val id = frontier.removeLast()
                if (!region.mayIntersect(S2Cell(id))) continue
                output.add(id)

                val neighbors = id.getEdgeNeighbors()
                for (edge in 0..3) {
                    val nbr = neighbors[edge]
                    if (all.add(nbr)) {
                        frontier.add(nbr)
                    }
                }
            }
        }
    }
}

data class S2RegionCovererQueueEntry(val id: Int, val candidate: S2RegionCoverer.Candidate) : Comparable<S2RegionCovererQueueEntry> {
    /**
     * Compares this object with the specified object for order. Returns zero if this object is equal
     * to the specified [other] object, a negative number if it's less than [other], or a positive number
     * if it's greater than [other].
     */
    override fun compareTo(other: S2RegionCovererQueueEntry): Int = -id.compareTo(other.id)

}
