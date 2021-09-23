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

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLE
import dilivia.collections.Container
import dilivia.collections.isSorted
import dilivia.collections.lowerBound
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.coords.S2Coords
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import java.math.BigInteger


/**
 * A S2CellUnion is a region consisting of cells of various sizes.  Typically a cell union is used to approximate some
 * other shape.  There is a tradeoff between the accuracy of the approximation and how many cells are used. Unlike
 * polygons, cells have a fixed hierarchical structure.  This makes them more suitable for optimizations based on
 * preprocessing.
 *
 * A S2CellUnion is represented as a vector of sorted, non-overlapping S2CellIds. By default the vector is also
 * "normalized", meaning that groups of 4 child cells have been replaced by their parent cell whenever possible.
 * S2CellUnions are not required to be normalized, but certain operations will return different results if they are not
 * (e.g., Contains(S2CellUnion).)
 *
 * This class is a port of the S2CellUnion class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2CellUnion private constructor(private val cellIds: MutableList<S2CellId>, verbatim: Boolean) : S2Region,
    Iterable<S2CellId> by cellIds,
    Container<S2CellId> {

    init {
        if (!verbatim) {
            normalize()
        }
    }

    // Constructs a cell union with the given S2CellIds, then calls Normalize()
    // to sort them, remove duplicates, and merge cells when possible.  (See
    // FromNormalized if your vector is already normalized.)
    //
    // The argument is passed by value, so if you are passing a named variable
    // and have no further use for it, consider using std::move().
    //
    // A cell union containing a single S2CellId may be constructed like this:
    //
    //     S2CellUnion example({cell_id})
    constructor(cellIds: List<S2CellId>) : this(cellIds.toMutableList(), false)
    constructor(vararg cellId: S2CellId) : this(cellId.toMutableList(), false)

    // Creates an empty cell union.
    constructor() : this(emptyList())

    // Clears the contents of the cell union and minimizes memory usage.
    fun clear() {
        cellIds.clear()
    }

    override fun size(): Int = numCells()
    fun numCells(): Int = cellIds.size

    fun cellId(i: Int): S2CellId = cellIds[i]

    // Returns true if the cell union is valid, meaning that the S2CellIds are
    // valid, non-overlapping, and sorted in increasing order.
    fun isValid(): Boolean {
        if (numCells() > 0 && !cellId(0).isValid) return false
        for (i in 1 until numCells()) {
            if (!cellId(i).isValid) return false
            if (cellId(i - 1).rangeMax() >= cellId(i).rangeMin()) return false
        }
        return true
    }

    // Returns true if the cell union is normalized, meaning that it is
    // satisfies IsValid() and that no four cells have a common parent.
    // Certain operations such as Contains(S2CellUnion) will return a different
    // result if the cell union is not normalized.
    fun isNormalized(): Boolean {
        if (numCells() > 0 && !cellId(0).isValid) return false
        for (i in 1 until numCells()) {
            if (!cellId(i).isValid) return false
            if (cellId(i - 1).rangeMax() >= cellId(i).rangeMin()) return false
            if (i >= 3 && areSiblings(cellId(i - 3), cellId(i - 2), cellId(i - 1), cellId(i))) {
                return false
            }
        }
        return true
    }

    // Normalizes the cell union by discarding cells that are contained by other
    // cells, replacing groups of 4 child cells by their parent cell whenever
    // possible, and sorting all the cell ids in increasing order.
    //
    // Returns true if the number of cells was reduced.
    // TODO(ericv): Change this method to return void.
    fun normalize(): Boolean = normalize(cellIds)

    // Replaces "output" with an expanded version of the cell union where any
    // cells whose level is less than "min_level" or where (level - min_level)
    // is not a multiple of "level_mod" are replaced by their children, until
    // either both of these conditions are satisfied or the maximum level is
    // reached.
    //
    // This method allows a covering generated by S2RegionCoverer using
    // min_level() or level_mod() constraints to be stored as a normalized cell
    // union (which allows various geometric computations to be done) and then
    // converted back to the original list of cell ids that satisfies the
    // desired constraints.
    fun denormalize(min_level: Int, level_mod: Int): List<S2CellId> = denormalize(cellIds, min_level, level_mod)

    // Returns true if the cell union contains the given cell id.  Containment
    // is defined with respect to regions, e.g. a cell contains its 4 children.
    // This is a fast operation (logarithmic in the size of the cell union).
    //
    // CAVEAT: If you have constructed a non-normalized S2CellUnion using
    // FromVerbatim, note that groups of 4 child cells are *not* considered to
    // contain their parent cell.  To get this behavior you must use one of the
    // other constructors or call Normalize() explicitly.
    fun contains(id: S2CellId): Boolean {
        // This is an exact test.  Each cell occupies a linear span of the S2
        // space-filling curve, and the cell id is simply the position at the center
        // of this span.  The cell union ids are sorted in increasing order along
        // the space-filling curve.  So we simply find the pair of cell ids that
        // surround the given cell id (using binary search).  There is containment
        // if and only if one of these two cell ids contains this cell.

        var i = cellIds.indexOfFirst { cellId -> cellId >= id }
        if (i == -1) i = cellIds.size
        if (i != cellIds.size && cellIds[i].rangeMin() <= id) return true
        return i != 0 && cellIds[i - 1].rangeMax() >= id
    }

    // Returns true if the cell union intersects the given cell id.
    // This is a fast operation (logarithmic in the size of the cell union).
    fun intersects(id: S2CellId): Boolean {
        // This is an exact test; see the comments for Contains() above.
        val i = cellIds.lowerBound(0, cellIds.size, id)
        if (i != cellIds.size && cellIds[i].rangeMin() <= id.rangeMax()) return true
        return i != 0 && cellIds[i - 1].rangeMax() >= id.rangeMin()
    }

    // Returns true if this cell union contains the given other cell union.
    //
    // CAVEAT: If you have constructed a non-normalized S2CellUnion using
    // FromVerbatim, note that groups of 4 child cells are *not* considered to
    // contain their parent cell.  To get this behavior you must use one of the
    // other constructors or call Normalize() explicitly.
    fun contains(y: S2CellUnion): Boolean {
        // TODO(ericv): A divide-and-conquer or alternating-skip-search
        // approach may be sigificantly faster in both the average and worst case.
        for (y_id in y) {
            if (!contains(y_id)) return false
        }
        return true
    }

    // Returns true if this cell union intersects the given other cell union.
    fun intersects(y: S2CellUnion): Boolean {
        // TODO(ericv): A divide-and-conquer or alternating-skip-search
        // approach may be sigificantly faster in both the average and worst case.

        for (y_id in y) {
            if (intersects(y_id)) return true
        }
        return false
    }

    // Returns the union of the two given cell unions.
    fun union(y: S2CellUnion): S2CellUnion {
        return S2CellUnion(cellIds + y.cellIds)
    }

    // Returns the intersection of the two given cell unions.
    fun intersection(y: S2CellUnion): S2CellUnion {
        checkState { cellIds.isSorted() }
        requireArgument { y.cellIds.isSorted() }

        // This is a fairly efficient calculation that uses binary search to skip
        // over sections of both input vectors.  It takes logarithmic time if all the
        // cells of "x" come before or after all the cells of "y" in S2CellId order.
        val out = mutableListOf<S2CellId>()
        var i = 0
        var j = 0
        while (i != cellIds.size && j != y.cellIds.size) {
            val cellI = cellIds[i]
            val cellJ = y.cellIds[j]
            val imin = cellI.rangeMin()
            val jmin = cellJ.rangeMin()
            if (imin > jmin) {
                // Either j->contains(*i) or the two cells are disjoint.
                if (cellI <= cellJ.rangeMax()) {
                    out.add(cellI)
                    i++
                } else {
                    // Advance "j" to the first cell possibly contained by *i.
                    j = y.cellIds.lowerBound(j + 1, y.cellIds.size, imin)
                    // y.cellIds.subList(j + 1, y.cellIds.size).indexOfFirst { cellId -> cellId >= imin }
                    // std::lower_bound(j + 1, y.end(), imin)
                    if (j == -1) j = y.cellIds.size
                    // The previous cell *(j-1) may now contain *i.
                    if (cellI <= y.cellIds[j - 1].rangeMax()) --j
                }
            } else if (jmin > imin) {
                // Identical to the code above with "i" and "j" reversed.
                if (cellJ <= cellIds[i].rangeMax()) {
                    out.add(cellJ)
                    j++
                } else {
                    i = cellIds.lowerBound(i + 1, cellIds.size, jmin)
                    // cellIds.subList(i + 1, cellIds.size).indexOfFirst { cellId -> cellId >= jmin }
                    // std::lower_bound(i + 1, x.end(), jmin)
                    //if (i == -1) i = cellIds.size
                    if (cellJ <= cellIds[i - 1].rangeMax()) --i
                }
            } else {
                // "i" and "j" have the same range_min(), so one contains the other.
                if (cellI < cellJ) {
                    out.add(cellI)
                    i++
                } else {
                    out.add(cellJ)
                    j++
                }
            }
        }

        val intersection = S2CellUnion(out, true)
        // The output is generated in sorted order.
        checkState { out.isSorted() }
        // The output is normalized as long as at least one input is normalized.
        checkState { intersection.isNormalized() || (!isNormalized() && !intersection.isNormalized()) }

        return intersection
    }

    // Specialized version of GetIntersection() that returns the intersection of
    // a cell union with an S2CellId.  This can be useful for splitting a cell
    // union into pieces.
    fun intersection(id: S2CellId): S2CellUnion {
        val result = mutableListOf<S2CellId>()
        if (contains(id)) {
            result.add(id)
        } else {
            var i = cellIds.lowerBound(0, cellIds.size, id.rangeMin())
            val idMax = id.rangeMax()
            while (i != cellIds.size && cellIds[i] <= idMax) result.add(cellIds[i++])
        }
        val intersection = S2CellUnion(result, true)
        checkState { intersection.isNormalized() || !isNormalized() }
        return intersection
    }

    // Returns the difference of the two given cell unions.
    fun difference(y: S2CellUnion): S2CellUnion {
        // TODO(ericv): this is approximately O(N*log(N)), but could probably
        // use similar techniques as GetIntersection() to be more efficient.

        val result = mutableListOf<S2CellId>()
        for (id in this) {
            result.addAll(getDifferenceInternal(id, y))
        }
        // The output is normalized as long as the first argument is normalized.
        val difference = S2CellUnion(result, true)
        checkState { difference.isNormalized() || !isNormalized() }
        return difference
    }

    // Expands the cell union by adding a buffer of cells at "expand_level"
    // around the union boundary.
    //
    // For each cell "c" in the union, we add all neighboring cells at level
    // "expand_level" that are adjacent to "c".  Note that there can be many
    // such cells if "c" is large compared to "expand_level".  If "c" is smaller
    // than "expand_level", we first add the parent of "c" at "expand_level" and
    // then add all the neighbors of that cell.
    //
    // Note that the size of the output is exponential in "expand_level".  For
    // example, if expand_level == 20 and the input has a cell at level 10,
    // there will be on the order of 4000 adjacent cells in the output.  For
    // most applications the Expand(min_radius, max_level_diff) method below is
    // easier to use.
    fun expand(expand_level: Int) {
        val output = mutableListOf<S2CellId>()
        val levelLsb = S2CellId.lsbForLevel(expand_level)
        var i = numCells()
        while (--i >= 0) {
            var id = cellId(i)
            if (id.lsb() < levelLsb) {
                id = id.parent(expand_level)
                // Optimization: skip over any cells contained by this one.  This is
                // especially important when very small regions are being expanded.
                while (i > 0 && id.contains(cellId(i - 1))) --i
            }
            output.add(id)
            id.appendAllNeighbors(expand_level, output)
        }
        cellIds.clear()
        cellIds.addAll(output)
        normalize()
    }

    // Expands the cell union such that it contains all points whose distance to
    // the cell union is at most "min_radius", but do not use cells that are
    // more than "max_level_diff" levels higher than the largest cell in the
    // input.  The second parameter controls the tradeoff between accuracy and
    // output size when a large region is being expanded by a small amount
    // (e.g. expanding Canada by 1km).  For example, if max_level_diff == 4 the
    // region will always be expanded by approximately 1/16 the width of its
    // largest cell.  Note that in the worst case, the number of cells in the
    // output can be up to 4 * (1 + 2 ** max_level_diff) times larger than the
    // number of cells in the input.
    fun expand(min_radius: S1Angle, max_level_diff: Int) {
        var min_level = S2CellId.kMaxLevel
        for (id in this) {
            min_level = min(min_level, id.level())
        }
        // Find the maximum level such that all cells are at least "min_radius" wide.
        val radius_level = S2Coords.projection.kMinWidth.getLevelForMinValue(min_radius.radians)
        if (radius_level == 0 && min_radius.radians > S2Coords.projection.kMinWidth.getValue(0)) {
            // The requested expansion is greater than the width of a face cell.
            // The easiest way to handle this is to expand twice.
            expand(0)
        }
        expand(min(min_level + max_level_diff, radius_level))
    }

    // The number of leaf cells covered by the union.
    // This will be no more than 6*2^60 for the whole sphere.
    fun leafCellsCovered(): ULong {
        var numLeaves = BigInteger.ZERO
        for (id in this) {
            val invertedLevel = S2CellId.kMaxLevel - id.level()
            numLeaves += (BigInteger.ONE.shiftLeft(invertedLevel shl 1))
        }
        return numLeaves.toLong().toULong()

        /*
        var num_leaves = 0UL
        for (id in this) {
            val invertedLevel = S2CellId.kMaxLevel - id.level()
            num_leaves += (1UL shl (invertedLevel shl 1))
        }
        return num_leaves*/
    }

    // Approximates this cell union's area in steradians by summing the average
    // area of each contained cell's average area, using the AverageArea method
    // from the S2Cell class.  This is equivalent to the number of leaves covered,
    // multiplied by the average area of a leaf.  Note that AverageArea does not
    // take into account distortion of cell, and thus may be off by up to a
    // factor of up to 1.7.
    //
    // NOTE: Since this is proportional to LeafCellsCovered(), it is
    // always better to use that function if all you care about is
    // the relative average area between objects.
    fun averageBasedArea(): Double = S2Cell.averageArea(S2CellId.kMaxLevel) * leafCellsCovered().toDouble()

    // Calculates this cell union's area in steradians by summing the approximate
    // area for each contained cell, using the ApproxArea method from the S2Cell
    // class.
    fun approxArea(): Double {
        var area = 0.0
        for (id in this) {
            area += S2Cell(id).approxArea()
        }
        return area
    }

    // Calculates this cell union's area in steradians by summing the exact area
    // for each contained cell, using the Exact method from the S2Cell class.
    fun exactArea(): Double {
        var area = 0.0
        for (id in this) {
            area += S2Cell(id).exactArea()
        }
        return area
    }


    fun cellIds(): List<S2CellId> = cellIds.toList()

    fun listIterator(): ListIterator<S2CellId> = cellIds.listIterator()

    fun listIterator(index: Int): ListIterator<S2CellId> = cellIds.listIterator(index)

    fun begin(): S2CellId = cellIds.first()

    fun end(): S2CellId = cellIds.last()

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    override fun clone(): S2CellUnion = S2CellUnion(MutableList(cellIds.size) { i -> cellIds[i] }, true)
    override val capBound: S2Cap
        get() {
            // Compute the approximate centroid of the region.  This won't produce the
            // bounding cap of minimal area, but it should be close enough.
            if (cellIds.isEmpty()) return S2Cap.empty
            var centroid = S2Point(0, 0, 0)
            for (id in this) {
                val area = S2Cell.averageArea(id.level())
                centroid = centroid + area * id.toPoint()
            }
            if (centroid == S2Point(0, 0, 0)) {
                centroid = S2Point(1, 0, 0)
            } else {
                centroid = centroid.normalize()
            }

            // Use the centroid as the cap axis, and expand the cap angle so that it
            // contains the bounding caps of all the individual cells.  Note that it is
            // *not* sufficient to just bound all the cell vertices because the bounding
            // cap may be concave (i.e. cover more than one hemisphere).
            var cap = S2Cap.fromPoint(centroid)
            for (id in this) {
                val cellCapBound = S2Cell(id).capBound
                logger.trace { "Add cell cap bound $id: $cellCapBound" }
                cap = cap.addCap(cellCapBound)
                logger.trace { "Cap = $cap" }
            }
            logger.debug { "getCapBound = $cap" }
            return cap
        }

    override val rectBound: S2LatLngRect
        get() {
            var bound = S2LatLngRect.empty()
            for (id in this) {
                bound = bound.union(S2Cell(id).rectBound)
            }
            return bound
        }

    override fun contains(cell: S2Cell): Boolean = contains(cell.id())

    override fun mayIntersect(cell: S2Cell): Boolean = intersects(cell.id())

    override fun contains(p: S2Point): Boolean = contains(S2CellId.fromPoint(p))

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S2CellUnion) return false

        if (cellIds != other.cellIds) return false

        return true
    }

    override fun hashCode(): Int {
        return cellIds.hashCode()
    }

    fun isEmpty(): Boolean = cellIds.isEmpty()

    override operator fun get(i: Int): S2CellId = cellId(i)
    override fun toString(): String {
        return "S2CellUnion(cellIds=${cellIds.joinToString(" ; ")})"
    }


    companion object {

        private val logger = KotlinLogging.logger(S2CellUnion::class.java.name)

        fun getIntersection(x: List<S2CellId>, y: List<S2CellId>, out: MutableList<S2CellId>) {
            logger.trace {
                """ Get intersection of:
                |- x = $x
                |- y = $y 
            """.trimMargin()
            }
            requireArgument({ x.isSorted() }, { "x is not sorted: $x" })
            requireArgument { y.isSorted() }

            // This is a fairly efficient calculation that uses binary search to skip
            // over sections of both input vectors.  It takes logarithmic time if all the
            // cells of "x" come before or after all the cells of "y" in S2CellId order.

            out.clear()
            var iIdx = 0
            var jIdx = 0
            while (iIdx != x.size && jIdx != y.size) {
                var i = x[iIdx]
                var j = y[jIdx]
                val imin = i.rangeMin()
                val jmin = j.rangeMin()
                if (imin > jmin) {
                    // Either j->contains(*i) or the two cells are disjoint.
                    if (i <= j.rangeMax()) {
                        out.add(i)
                        ++iIdx
                    } else {
                        // Advance "j" to the first cell possibly contained by *i.
                        jIdx = y.lowerBound(jIdx + 1, y.size, imin)
                        // The previous cell *(j-1) may now contain *i.
                        if (i <= y[(jIdx - 1)].rangeMax()) {
                            --jIdx
                        }
                    }
                } else if (jmin > imin) {
                    // Identical to the code above with "i" and "j" reversed.
                    if (j <= i.rangeMax()) {
                        out.add(j)
                        ++jIdx
                    } else {
                        iIdx = x.lowerBound(iIdx + 1, x.size, jmin)
                        if (j <= x[iIdx - 1].rangeMax()) {
                            --iIdx
                        }
                    }
                } else {
                    // "i" and "j" have the same range_min(), so one contains the other.
                    if (i < j) {
                        out.add(i)
                        ++iIdx
                    } else {
                        out.add(j)
                        ++jIdx
                    }
                }
            }
            // The output is generated in sorted order.
            checkState { out.isSorted() }

        }

        // Converts a vector of uint64 to a vector of S2CellIds.
        private fun toS2CellIds(ids: List<ULong>): List<S2CellId> = ids.map { S2CellId(it) }

        // Returns true if the given four cells have a common parent.
        // REQUIRES: The four cells are distinct.
        private fun areSiblings(a: S2CellId, b: S2CellId, c: S2CellId, d: S2CellId): Boolean {
            // A necessary (but not sufficient) condition is that the XOR of the
            // four cells must be zero.  This is also very fast to test.
            if ((a.id xor b.id xor c.id) != d.id) return false

            // Now we do a slightly more expensive but exact test.  First, compute a
            // mask that blocks out the two bits that encode the child position of
            // "id" with respect to its parent, then check that the other three
            // children all agree with "mask".
            var mask = d.lsb() shl 1
            mask = (mask + (mask shl 1)).inv()
            val idMasked = (d.id and mask)
            return ((a.id and mask) == idMasked &&
                    (b.id and mask) == idMasked &&
                    (c.id and mask) == idMasked &&
                    !d.isFace)
        }

        // Convenience constructor that accepts a vector of uint64.  Note that
        // unlike the constructor above, this one makes a copy of "cell_ids".
        fun fromIds(ids: List<ULong>): S2CellUnion = S2CellUnion(toS2CellIds(ids))

        // Constructs a cell union for the whole sphere.
        fun wholeSphere(): S2CellUnion = S2CellUnion((0..5).map { S2CellId.fromFace(it) })

        // Constructs a cell union from S2CellIds that have already been normalized
        // (typically because they were extracted from another S2CellUnion).
        //
        // The argument is passed by value, so if you are passing a named variable
        // and have no further use for it, consider using std::move().
        //
        // REQUIRES: "cell_ids" satisfies the requirements of IsNormalized().
        fun fromNormalized(cellIds: List<S2CellId>): S2CellUnion {
            val result = S2CellUnion(cellIds)
            assert(result.isNormalized())
            return result
        }

        // Constructs a cell union from a vector of sorted, non-overlapping
        // S2CellIds.  Unlike the other constructors, FromVerbatim does not require
        // that groups of 4 child cells have been replaced by their parent cell.  In
        // other words, "cell_ids" must satisfy the requirements of IsValid() but
        // not necessarily IsNormalized().
        //
        // Note that if the cell union is not normalized, certain operations may
        // return different results (e.g., Contains(S2CellUnion)).
        //
        // REQUIRES: "cell_ids" satisfies the requirements of IsValid().
        fun fromVerbatim(cellIds: List<S2CellId>): S2CellUnion {
            val result = S2CellUnion(cellIds.toMutableList(), true)
            checkState { result.isValid() }
            return result
        }

        // Constructs a cell union that corresponds to a continuous range of cell
        // ids.  The output is a normalized collection of cell ids that covers the
        // leaf cells between "min_id" and "max_id" inclusive.
        //
        // REQUIRES: min_id.is_leaf(), max_id.is_leaf(), min_id <= max_id.
        fun fromMinMax(minId: S2CellId, maxId: S2CellId): S2CellUnion {
            requireArgument { maxId.isValid }
            return fromBeginEnd(minId, maxId.next())
        }

        // Like FromMinMax() except that the union covers the range of leaf cells
        // from "begin" (inclusive) to "end" (exclusive), as with Python ranges or
        // STL iterator ranges.  If (begin == end) the result is empty.
        //
        // REQUIRES: begin.is_leaf(), end.is_leaf(), begin <= end.
        fun fromBeginEnd(begin: S2CellId, end: S2CellId): S2CellUnion {
            requireArgument { begin.isLeaf }
            requireArgument { end.isLeaf }
            requireLE(begin, end)

            // We repeatedly add the largest cell we can.
            var id = begin.maximumTile(end)
            val cellIds = mutableListOf<S2CellId>()
            while (id != end) {
                cellIds.add(id)
                id = id.next().maximumTile(end)
            }
            // The output is already normalized.
            val output = S2CellUnion(cellIds, true)
            checkState { output.isNormalized() }
            return output
        }

        fun normalize(cellIds: MutableList<S2CellId>): Boolean {
            // Optimize the representation by discarding cells contained by other cells,
            // and looking for cases where all subcells of a parent cell are present.
            cellIds.sort()
            var out = 0
            var currentId: S2CellId
            for (id in cellIds) {
                currentId = id
                // Check whether this cell is contained by the previous cell.
                if (out > 0 && cellIds[out - 1].contains(currentId)) continue

                // Discard any previous cells contained by this cell.
                while (out > 0 && currentId.contains(cellIds[out - 1])) --out

                // Check whether the last 3 elements plus "id" can be collapsed into a
                // single parent cell.
                while (out >= 3 && areSiblings(cellIds[out - 3], cellIds[out - 2], cellIds[out - 1], currentId)) {
                    // Replace four children by their parent cell.
                    currentId = currentId.parent()
                    out -= 3
                }
                cellIds[out++] = currentId
            }
            if (cellIds.size == out) return false
            while (cellIds.size > out) cellIds.removeLast()
            return true
        }

        fun denormalize(cellIds: List<S2CellId>, min_level: Int, level_mod: Int): List<S2CellId> {
            val out = mutableListOf<S2CellId>()
            denormalize(cellIds, min_level, level_mod, out)
            return out
        }

        fun denormalize(cellIds: List<S2CellId>, min_level: Int, level_mod: Int, out: MutableList<S2CellId>) {
            requireGE(min_level, 0)
            requireLE(min_level, S2CellId.kMaxLevel)
            requireGE(level_mod, 1)
            requireLE(level_mod, 3)

            out.clear()
            for (id in cellIds) {
                val level = id.level()
                var new_level = max(min_level, level)
                if (level_mod > 1) {
                    // Round up so that (new_level - min_level) is a multiple of level_mod.
                    // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
                    new_level += (S2CellId.kMaxLevel - (new_level - min_level)) % level_mod
                    new_level = min(S2CellId.kMaxLevel, new_level)
                }
                if (new_level == level) {
                    out.add(id)
                } else {
                    val end = id.childEnd(new_level)
                    var cellId = id.childBegin(new_level)
                    while (cellId != end) {
                        out.add(cellId)
                        cellId = cellId.next()
                    }
                }
            }
        }

        private fun getDifferenceInternal(cell: S2CellId, y: S2CellUnion): List<S2CellId> {
            // Add the difference between cell and y to cell_ids.
            // If they intersect but the difference is non-empty, divide and conquer.
            val cell_ids = mutableListOf<S2CellId>()
            if (!y.intersects(cell)) {
                cell_ids.add(cell)
            } else if (!y.contains(cell)) {
                var child = cell.childBegin()
                for (i in 0..3) {
                    cell_ids.addAll(getDifferenceInternal(child, y))
                    if (i == 3) break  // Avoid unnecessary next() computation.
                    child = child.next()
                }
            }
            return cell_ids
        }
    }

    object S2CellUnionTestPeer {
        fun fromVerbatimNoChecks(cellIds: List<S2CellId>): S2CellUnion = S2CellUnion(cellIds.toMutableList(), true)
    }

}
