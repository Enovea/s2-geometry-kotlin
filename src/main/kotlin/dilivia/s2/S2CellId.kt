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
package dilivia.s2

import dilivia.Bits
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireBetween
import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLE
import dilivia.PreConditions.requireLT
import dilivia.math.DoubleType
import dilivia.math.R1Interval
import dilivia.math.R2Rect
import dilivia.math.vectors.R2Point
import dilivia.s2.S1Angle.Companion.sin
import dilivia.s2.coords.FaceIJ
import dilivia.s2.coords.FaceSiTi
import dilivia.s2.coords.LookupCellTables
import dilivia.s2.coords.LookupCellTables.kLookupBits
import dilivia.s2.coords.LookupCellTables.lookupIj
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2Coords.kInvertMask
import dilivia.s2.coords.S2Coords.kPosToOrientation
import dilivia.s2.coords.S2Coords.kSwapMask
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * An S2CellId is a 64-bit unsigned integer that uniquely identifies a cell in the S2 cell decomposition. It has the
 * following format:
 *
 *   id = [face][facePos]
 *
 *   face:     a 3-bit number (range 0..5) encoding the cube face.
 *
 *   face_pos: a 61-bit number encoding the position of the center of this cell along the Hilbert curve over this face
 *            (see the Wiki pages for details).
 *
 * Sequentially increasing cell ids follow a continuous space-filling curve over the entire sphere. They have the
 * following properties:
 *
 *  - The id of a cell at level k consists of a 3-bit face number followed by k bit pairs that recursively select one
 *    of the four children of each cell. The next bit is always 1, and all other bits are 0. Therefore, the level of a
 *    cell is determined by the position of its lowest-numbered bit that is turned on (for a cell at level k, this
 *    position is 2 * (kMaxLevel - k).)
 *
 *  - The id of a parent cell is at the midpoint of the range of ids spanned by its children (or by its descendants at
 *    any level).
 *
 * Leaf cells are often used to represent points on the unit sphere, and this class provides methods for converting
 * directly between these two representations. For cells that represent 2D regions rather than discrete point, it is
 * better to use the S2Cell class.
 *
 * This class is a port of the S2CellId class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @see S2Coords
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2CellId : Comparable<S2CellId> {

    /** The 64-bit unique identifier for this cell. */
    val id: ULong

    /**
     * Creates a cell id with the specified ULong value.
     *
     * @param id The id of the cell as ULong.
     */
    constructor(id: ULong) {
        this.id = id
    }

    /** The default constructor returns an invalid cell id (0UL). */
    constructor() : this(0UL)

    /**
     * Get the direction vector corresponding to the center of the given cell. The vector returned by ToPointRaw is
     * not necessarily unit length. This method returns the same result as S2Cell.getCenter().
     *
     * The maximum directional error in ToPoint() (compared to the exact mathematical result) is 1.5 * DBL_EPSILON
     * radians, and the maximum length error is 2 * DBL_EPSILON (the same as normalize).
     *
     * @return The direction vector of the cell center.
     */
    fun toPoint(): S2Point = toPointRaw().normalize()

    /**
     * Get the direction vector corresponding to the center of the cell. The vector returned by toPointRaw is not
     * necessarily unit length.
     *
     * @return The direction vector of the cell center (not necessarily normalized).
     */
    fun toPointRaw(): S2Point = S2Coords.faceSiTiToXyz(centerSiTi)

    /**
     * Get the S2LatLng corresponding to the center of the given cell.
     *
     * @return The center as S2LatLng.
     */
    fun toLatLng(): S2LatLng = S2LatLng.fromPoint(toPointRaw())

    /**
     * Center of the cell in (s,t) coordinates.
     * @see S2Coords
     */
    val centerST: R2Point
        get() {
            val (_, si, ti) = centerSiTi
            return R2Point(S2Coords.siTiToSt(si), S2Coords.siTiToSt(ti))
        }

    /**
     * The edge length of this cell in (s,t)-space. In the range [0, 1]
     */
    val sizeST: Double
        get() = getSizeST(level())

    /**
     * The bounds of this cell in (s,t)-space.
     */
    val boundST: R2Rect
        get() {
            val size = sizeST
            return R2Rect.fromCenterSize(centerST, R2Point(size, size))
        }

    /**
     * The center of the cell in (u,v) coordinates.
     * Note that the center of the cell is defined as the point at which it is recursively subdivided into four
     * children. in general, it is not at the midpoint of the (u,v) rectangle covered by the cell.
     */
    val centerUV: R2Point
        get() {
            val (_, si, ti) = centerSiTi
            return R2Point(
                S2Coords.stToUv(S2Coords.siTiToSt(si)),
                S2Coords.stToUv(S2Coords.siTiToSt(ti))
            )
        }

    /**
     * The bound of this cell in (u,v)-space.
     */
    val boundUV: R2Rect
        get() {
            val (_, i, j) = toFaceIJOrientation(false)
            return ijLevelToBoundUV(i, j, level())
        }

    /**
     * Get the (face, si, ti) coordinates of the center of the cell.
     * Note that although (si,ti) coordinates span the range [0,2**31] in general, the cell center coordinates are
     * always in the range [1,2**31-1] and therefore can be represented using a signed 32-bit integer.
     *
     * @return The (face, si, ti) coordinates of the cell.
     */
    val centerSiTi: FaceSiTi
        get() {
            // First we compute the discrete (i,j) coordinates of a leaf cell contained within the given cell.
            // Given that cells are represented by the Hilbert curve position corresponding at their center, it turns out
            // that the cell returned by ToFaceIJOrientation is always one of two leaf cells closest to the center of the
            // cell (unless the given cell is a leaf cell itself, in which case there is only one possibility).
            //
            // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin, jmin) be the coordinates of its
            // lower left-hand corner, the leaf cell returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
            // (imin + s/2 - 1, jmin + s/2 - 1).  The first case is the one we want.
            // We can distinguish these two cases by looking at the low bit of "i" or "j". In the second case the low bit
            // is one, unless s == 2 (i.e. the level just above leaf cells) in which case the low bit is zero.
            //
            // In the code below, the expression ((i ^ (int(id_) >> 2)) & 1) is true if we are in the second case described
            // above.
            val faceIJ = toFaceIJOrientation(false)
            val delta = if (isLeaf) 1U else if (((faceIJ.i xor (id.toInt() shr 2)) and 1) != 0) 2U else 0U

            // Note that (2 * {i,j} + delta) will never overflow a 32-bit integer.
            return FaceSiTi(
                face = faceIJ.face,
                si = 2U * faceIJ.i.toUInt() + delta,
                ti = 2U * faceIJ.j.toUInt() + delta
            )
        }

    /**
     * Indicates if the id represents a valid cell.
     *
     * All methods require isValid to be true unless otherwise specified (although not all methods enforce this).
     */
    val isValid: Boolean
        get() = (face() < kNumFaces && (lsb() and 0x1555555555555555UL) != 0UL)

    /**
     * Which cube face this cell belongs to, in the range 0..5.
     *
     * @return The face index.
     */
    fun face(): Int = (id shr kPosBits).toInt()

    /**
     * The position of the cell center along the Hilbert curve over this face, in
     * the range 0..(2**kPosBits-1).
     *
     * @return The cell position.
     */
    fun pos(): ULong = id and (0UL.inv() shr kFaceBits)

    /**
     * Get the subdivision level of the cell (range 0..kMaxLevel).
     *
     * @return The cell level.
     * @see kMaxLevel
     */
    fun level(): Int {
        // We can't just assert(isValid()) because we want level() to be
        // defined for end-dilivia.iterators, i.e. S2CellId::End(kLevel).  However there is
        // no good way to define S2CellId::None().level(), so we do prohibit that.
        assert(id != 0UL)

        // A special case for leaf cells is not worthwhile.
        return kMaxLevel - (Bits.findLSBSetNonZero64(id) shr 1)
    }

    /**
     * The edge length of this cell in (i,j)-space.
     */
    val sizeIJ: Int
        get() = getSizeIJ(level())

    /**
     * Indicate if this is a leaf cell (more efficient than checking whether level() == kMaxLevel).
     */
    val isLeaf: Boolean
        get() = id.toInt() and 1 != 0

    /**
     * Indicate if this is a top-level face cell (more efficient than checking whether level() == 0).
     */
    val isFace: Boolean
        get() = (id and (lsbForLevel(0) - 1UL)) == 0UL

    /**
     * Get the child position (0..3) of this cell's ancestor at the given level, relative to its parent. The argument
     * should be in the range 1..kMaxLevel. For example, childPosition(1) returns the position of this cell's
     * level-1 ancestor within its top-level face cell.
     *
     * @param level Level of the child.
     * @return The position at the given level.
     */
    fun childPosition(level: Int): Int {
        checkState { isValid }
        requireBetween(level, 1, this.level() + 1)
        return (id shr (2 * (kMaxLevel - level) + 1)).toInt() and 3
    }

    /**
     * Get the child position (0..3) of this cell within its parent.
     * REQUIRES: level() >= 1.
     * @return The position of this cell within its parent.
     */
    fun childPosition(): Int {
        val level = level()
        checkState { level >= 1 }
        return childPosition(level)
    }

    /**
     * Get the first id of the cells contained within this cell.
     *
     * These methods return the range of cell ids that are contained within this cell (including itself). The range
     * is *inclusive* (i.e. test using >= and <=) and the return values of both methods are valid leaf cell ids.
     * In other words, a.contains(b) if and only if
     *
     *     (b >= a.range_min() && b <= a.range_max())
     *
     * If you want to iterate through all the descendants of this cell at a particular level, use childBegin(level)
     * and childEnd(level) instead.
     * Also see maximumTile(), which can be used to iterate through a range of cells using S2CellIds at different
     * levels that are as large as possible.
     *
     * If you need to convert the range to a semi-open interval [min, limit) (e.g., in order to use a key-value store
     * that only supports semi-open range queries), do not attempt to define "limit" as range_max.next().
     * The problem is that leaf S2CellIds are 2 units apart, so the semi-open interval [min, limit) includes an
     * additional value (range_max.id() + 1) which is happens to be a valid S2CellId about one-third of the time and
     * is *never* contained by this cell.  (It always correpsonds to a cell that is larger than this one.)  You can
     * define "limit" as (range_max.id() + 1) if necessary (which is not always a valid S2CellId but can still be used
     * with fromToken/toToken), or you can convert rangeMax() to the key space of your key-value store and define
     * "limit" as successor(key).
     *
     * Note that Sentinel().range_min() == Sentinel.range_max() == Sentinel().
     *
     * @return This cell range start.
     * @see rangeMax
     */
    fun rangeMin(): S2CellId = S2CellId(id - (lsb() - 1UL))

    /**
     * Get the last id of the cells contained within this cell.
     *
     * @return This cell range end.
     * @see rangeMax
     */
    fun rangeMax(): S2CellId = S2CellId(id + (lsb() - 1UL))

    /**
     * Get an iterator over the cell range contained by this cell.
     *
     * @return An iterator over all the cells contained by this one.
     */
    fun rangeIterator(): Iterator<S2CellId> = S2CellIdIterator(parent = this)

    /**
     * Indicate if the given cell is contained within this one.
     * @param other A cell to test.
     * @return true if this cell contains the given one.
     */
    operator fun contains(other: S2CellId): Boolean {
        checkState { isValid }
        requireArgument { other.isValid }
        return other in rangeMin()..rangeMax()
    }

    /**
     * Indicate if the given cell intersects this one.
     *
     * @param other A cell id.
     * @return true if this cell intersects the given one.
     */
    fun intersects(other: S2CellId): Boolean {
        checkState { isValid }
        requireArgument { other.isValid }
        return other.rangeMin() <= rangeMax() && other.rangeMax() >= rangeMin()
    }

    /**
     * Get the cell at the previous level.
     *
     * @return The parent cell.
     */
    fun parent(): S2CellId {
        checkState { isValid }
        checkState { !isFace }
        val newLsb = lsb() shl 2
        return S2CellId((id and (newLsb.inv() + 1UL)) or newLsb)
    }

    /**
     * Get the cell at the given level (which must be less than or equal to the current level).
     * @param level The level.
     * @return The parent of this cell at the given level.
     */
    fun parent(level: Int): S2CellId {
        checkState { isValid }
        requireBetween(level, 0, level() + 1)
        val newLsb = lsbForLevel(level)
        return S2CellId((id and (newLsb.inv() + 1UL)) or newLsb)
    }

    /**
     * Get the immediate child of this cell at the given traversal order position (in the range 0 to 3).
     * This cell must not be a leaf cell.
     *
     * @param position The position (0..3)
     * @return The child of this cell at the given position.
     */
    fun child(position: Int): S2CellId {
        checkState { isValid }
        checkState { !isLeaf }
        // To change the level, we need to move the least-significant bit two
        // positions downward.  We do this by subtracting (4 * new_lsb) and adding
        // new_lsb.  Then to advance to the given child cell, we add
        // (2 * position * new_lsb).
        val newLsb = lsb() shr 2
        return S2CellId(id + (2UL * position.toULong() + 1UL - 4UL) * newLsb)
    }

    /**
     * Get the first child of this cell.
     *
     * Iterator-style methods for traversing the immediate children of a cell or
     * all of the children at a given level (greater than or equal to the current
     * level).  Note that the end value is exclusive, just like standard STL
     * dilivia.iterators, and may not even be a valid cell id.  You should iterate using
     * code like this:
     *
     *   var c = id.childBegin()
     *   while(c != id.childEnd()) {
     *     ...
     *     c = c.next()
     *   }
     *
     * The convention for advancing the iterator is "c = c.next()"
     *
     * @return The first immediate child.
     */
    fun childBegin(): S2CellId {
        checkState { isValid }
        checkState { !isLeaf }
        val oldLsb = lsb()
        return S2CellId(id - oldLsb + (oldLsb shr 2))
    }

    /**
     * Get the first child of this cell at the given level.
     *
     * @param level A level in range this.level()..kMaxLevel
     * @return The first immediate child.
     * @see childEnd
     */
    fun childBegin(level: Int): S2CellId {
        checkState { isValid }
        requireBetween(level, this.level(), kMaxLevel + 1)
        return S2CellId(id - lsb() + lsbForLevel(level))
    }

    /**
     * Get the last child of this cell.
     *
     * @return The last immediate child.
     * @see childBegin
     */
    fun childEnd(): S2CellId {
        checkState { isValid }
        checkState { !isLeaf }
        val oldLsb = lsb()
        return S2CellId(id + oldLsb + (oldLsb shr 2))
    }

    /**
     * Get the last child of this cell at the given level.
     *
     * @param level A level in range this.level()..kMaxLevel
     * @return The last immediate child.
     * @see childBegin
     */
    fun childEnd(level: Int): S2CellId {
        checkState { isValid }
        requireBetween(level, this.level(), kMaxLevel + 1)
        return S2CellId(id + lsb() + lsbForLevel(level))
    }

    /**
     * Return the next cell at the same level along the Hilbert curve. Works correctly when advancing from one face
     * to the next, but does *not* wrap around from the last face to the first or vice versa.
     *
     * @return The next cell.
     */
    operator fun next(): S2CellId = S2CellId(id + (lsb() shl 1))

    /**
     * Return the previous cell at the same level along the Hilbert curve. Works correctly when advancing from one face
     * to the next, but does *not* wrap around from the last face to the first or vice versa.
     *
     * @return The previous cell.
     */
    fun previous(): S2CellId = S2CellId(id - (lsb() shl 1))

    /**
     * This method advances or retreats the indicated number of steps along the Hilbert curve at the current level,
     * and returns the new position. The position is never advanced past end() or before begin().
     *
     * @param steps The number of steps (can be negative)
     * @return The cell at the advanced position.
     */
    fun advance(steps: Long): S2CellId {
        if (steps == 0L) return this

        // We clamp the number of steps if necessary to ensure that we do not
        // advance past the End() or before the Begin() of this level.  Note that
        // min_steps and max_steps always fit in a signed 64-bit integer.
        var s = steps
        val stepShift = 2 * (kMaxLevel - level()) + 1
        if (s < 0L) {
            val minSteps = -(id shr stepShift).toLong()
            if (s < minSteps) s = minSteps
        } else {
            val maxSteps = ((kWrapOffset + lsb() - id) shr stepShift).toLong()
            if (s > maxSteps) s = maxSteps
        }
        // If steps is negative, then shifting it left has undefined behavior.
        // Cast to uint64 for a 2's complement answer.
        return S2CellId(id + (s.toULong() shl stepShift))
    }

    /**
     * Gets the number of steps that this cell is from Begin(level()). The return value is always non-negative.
     * @return The number of steps from the begin range.
     */
    fun distanceFromBegin(): ULong {
        val stepShift = 2 * (kMaxLevel - level()) + 1
        return id shr stepShift
    }

    /**
     * Like next(), but wraps around from the last face to the first and vice versa. Should *not* be used for iteration
     * in conjunction with childBegin(), childEnd(), begin(), or end().
     *
     * @return The next cell.
     */
    fun nextWrap(): S2CellId {
        checkState { isValid }
        val n = next()
        return if (n.id < kWrapOffset) n
        else S2CellId(n.id - kWrapOffset)
    }

    /**
     * Like prev(), but wraps around from the last face to the first and vice versa. Should *not* be used for iteration
     * in conjunction with childBegin(), childEnd(), begin(), or end().
     *
     * @return The prev cell.
     */
    fun prevWrap(): S2CellId {
        checkState { isValid }
        val p = previous()
        return if (p.id < kWrapOffset) p
        else S2CellId(p.id + kWrapOffset)
    }

    /**
     * This method advances or retreats the indicated number of steps along the Hilbert curve at the current level,
     * and returns the new position.  The position wraps between the first and last faces as necessary.  The input must
     * be a valid cell id.
     *
     * @param steps The number of steps.
     * @return The cell at the new position.
     */
    fun advanceWrap(steps: Long): S2CellId {
        checkState { isValid }
        if (steps == 0L) return this

        var wrappedSteps = steps
        val stepShift: Int = 2 * (kMaxLevel - level()) + 1
        if (wrappedSteps < 0) {
            val minSteps = -(id shr stepShift).toLong()
            if (wrappedSteps < minSteps) {
                val stepWrap = (kWrapOffset shr stepShift).toLong()
                wrappedSteps %= stepWrap
                if (wrappedSteps < minSteps) wrappedSteps += stepWrap
            }
        } else {
            // Unlike advance(), we don't want to return End(level).
            val maxSteps = ((kWrapOffset - id) shr stepShift).toLong()
            if (wrappedSteps > maxSteps) {
                val stepWrap = (kWrapOffset shr stepShift).toLong()
                wrappedSteps %= stepWrap
                if (wrappedSteps > maxSteps) wrappedSteps -= stepWrap
            }
        }
        return S2CellId(id + (wrappedSteps.toULong() shl stepShift))
    }

    /**
     * Get the largest cell with the same rangeMin() and such that rangeMax() < limit.rangeMin(). Returns "limit" if no
     * such cell exists.
     * This method can be used to generate a small set of S2CellIds that covers a given range (a "tiling").
     * This example shows how to generate a tiling for a semi-open range of leaf cells [start, limit):
     *
     *   var id = start.maximumTile(limit)
     *   while (id != limit) {
     *      ...
     *      id = id.next().maximumTile(limit)
     *   }
     *
     * Note that in general the cells in the tiling will be of different sizes; they gradually get larger (near the
     * middle of the range) and then gradually get smaller (as "limit" is approached).
     *
     * @param limit The cell limit.
     * @return The largest cell
     */
    fun maximumTile(limit: S2CellId): S2CellId {
        var cellId = this
        val start = cellId.rangeMin()
        if (start >= limit.rangeMin()) return limit

        if (cellId.rangeMax() >= limit) {
            // The cell is too large.  Shrink it.  Note that when generating coverings
            // of S2CellId ranges, this loop usually executes only once.  Also because
            // id.range_min() < limit.range_min(), we will always exit the loop by the
            // time we reach a leaf cell.
            do {
                cellId = cellId.child(0); } while (cellId.rangeMax() >= limit)
            return cellId
        }
        // The cell may be too small.  Grow it if necessary.  Note that generally
        // this loop only iterates once.
        while (!cellId.isFace) {
            val parent = cellId.parent()
            if (parent.rangeMin() != start || parent.rangeMax() >= limit) break
            cellId = parent
        }
        return cellId
    }

    internal fun getCommonAncestorLevel(other: S2CellId): Int {
        // Basically we find the first bit position at which the two S2CellIds
        // differ and convert that to a level.  The max() below is necessary for the
        // case where one S2CellId is a descendant of the other.
        val bits: ULong = kotlin.math.max(id xor other.id, kotlin.math.max(lsb(), other.lsb()))
        checkState { bits != 0UL }  // Because lsb() is non-zero.

        // Compute the position of the most significant bit, and then map the bit
        // position as follows:
        // {0} -> 30, {1,2} -> 29, {3,4} -> 28, ... , {59,60} -> 0, {61,62,63} -> -1.
        return FastMath.max(60 - Bits.findMSBSetNonZero64(bits), -1) shr 1
    }

    /**
     * Encode cell ids to compact text strings suitable for display or indexing. Cells at lower levels
     * (i.e. larger cells) are encoded into fewer characters. The maximum token length is 16.
     *
     * Tokens preserve ordering, i.e. toToken(x) < toToken(y) if x < y.
     *
     * toToken() returns a string by value for convenience; the compiler does this without intermediate copying in most
     * cases.
     *
     * These methods guarantee that fromToken(toToken(x)) == x even when "x" is an invalid cell id. All tokens are
     * alphanumeric strings. FromToken() returns S2CellId::None() for malformed inputs.
     *
     * @return A string representation of the cell id.
     */
    fun toToken(): String {
        // Simple implementation: print the id in hex without trailing zeros.
        // Using hex has the advantage that the tokens are case-insensitive, all
        // characters are alphanumeric, no characters require any special escaping
        // in queries for most indexing systems, and it's easy to compare cell
        // tokens against the feature ids of the corresponding features.
        //
        // Using base 64 would produce slightly shorter tokens, but for typical cell
        // sizes used during indexing (up to level 15 or so) the average savings
        // would be less than 2 bytes per cell which doesn't seem worth it.

        // "0" with trailing 0s stripped is the empty string, which is not a
        // reasonable token.  Encode as "X".
        if (id == 0UL) return "X"
        val numZeroDigits = Bits.findLSBSetNonZero64(id) / 4
        return hexFormatString(id shr (4 * numZeroDigits), 16 - numZeroDigits)
    }

    /**
     * Creates a human readable debug string.  Used for << and available for direct usage as well.
     * The format is "f/dd..d" where "f" is a digit in the range [0-5] representing the S2CellId face, and "dd..d" is
     * a string of digits in the range [0-3] representing each child's position with respect to its parent.
     * (Note that the latter string may be empty.)
     *
     * For example "4/" represents S2CellId.fromFace(4), and "3/02" represents S2CellId.fromFace(3).child(0).child(2).
     *
     * @return A debug string representation of the cell id.
     */
    override fun toString(): String {
        if (!isValid) {
            return "Invalid: ${java.lang.Long.toHexString(id.toLong()).padStart(16, '0')}"
        }
        var out = "${face()}/"
        for (current_level in 1..level()) {
            out += "0123"[childPosition(current_level)]
        }
        return out
    }

    /**
     * Get the four cells that are adjacent across the cell's four edges. Neighbors are returned in the order defined
     * by S2Cell.getEdge(). All neighbors are guaranteed to be distinct.
     *
     * @return The four adjacent cells
     */
    fun getEdgeNeighbors(): Array<S2CellId> {
        val level = this.level()
        val size = getSizeIJ(level)
        val (face, i, j) = toFaceIJOrientation(false)


        // Edges 0, 1, 2, 3 are in the down, right, up, left directions.
        val downCell = fromFaceIJSame(face, i, j - size, j - size >= 0)
        val rightCell = fromFaceIJSame(face, i + size, j, i + size < kMaxSize)
        val upCell = fromFaceIJSame(face, i, j + size, j + size < kMaxSize)
        val leftCell = fromFaceIJSame(face, i - size, j, i - size >= 0)
        return arrayOf(
            downCell.parent(level), rightCell.parent(level),
            upCell.parent(level), leftCell.parent(level)
        )
    }

    /**
     * Get the neighbors of closest vertex to this cell at the given level, by appending them to "output".
     * Normally there are four neighbors, but the closest vertex may only have three neighbors if it is one of the 8
     * cube vertices.
     *
     * Requires: level < this.level(), so that we can determine which vertex is closest (in particular,
     * level == kMaxLevel is not allowed).
     *
     * @param level A level.
     * @param output The list to append to.
     */
    fun appendVertexNeighbors(level: Int, output: MutableList<S2CellId>) {
        // "level" must be strictly less than this cell's level so that we can
        // determine which vertex this cell is closest to.
        requireLT(level, level())
        val (face, i, j) = toFaceIJOrientation(false)

        // Determine the i- and j-offsets to the closest neighboring cell in each
        // direction.  This involves looking at the next bit of "i" and "j" to
        // determine which quadrant of this->parent(level) this cell lies in.
        val halfsize = getSizeIJ(level + 1)
        val size = halfsize shl 1
        val isame: Boolean
        val jsame: Boolean
        val ioffset: Int
        val joffset: Int
        if (i and halfsize != 0) {
            ioffset = size
            isame = (i + size) < kMaxSize
        } else {
            ioffset = -size
            isame = (i - size) >= 0
        }
        if (j and halfsize != 0) {
            joffset = size
            jsame = (j + size) < kMaxSize
        } else {
            joffset = -size
            jsame = (j - size) >= 0
        }

        output.add(parent(level))
        output.add(fromFaceIJSame(face, i + ioffset, j, isame).parent(level))
        output.add(fromFaceIJSame(face, i, j + joffset, jsame).parent(level))
        // If i- and j- edge neighbors are *both* on a different face, then this
        // vertex only has three neighbors (it is one of the 8 cube vertices).
        if (isame || jsame) {
            output.add(fromFaceIJSame(face, i + ioffset, j + joffset, isame && jsame).parent(level))
        }
    }

    /**
     * Append all neighbors of this cell at the given level to "output". Two cells X and Y are neighbors if their
     * boundaries intersect but their interiors do not. In particular, two cells that intersect at a single point are
     * neighbors.
     *
     * Requires: nbrLevel >= this.level(). Note that for cells adjacent to a face vertex, the same neighbor may be
     * appended more than once.
     *
     * @param nbrLevel A level
     * @param output The list to append to.
     */
    fun appendAllNeighbors(nbrLevel: Int, output: MutableList<S2CellId>) {
        requireGE(nbrLevel, level())
        var (face, i, j, _) = toFaceIJOrientation(false)

        // Find the coordinates of the lower left-hand leaf cell.  We need to
        // normalize (i,j) to a known position within the cell because nbr_level
        // may be larger than this cell's level.
        val size = sizeIJ
        i = i and -size
        j = j and -size

        val nbrSize = getSizeIJ(nbrLevel)
        requireLE(nbrSize, size)

        // We compute the top-bottom, left-right, and diagonal neighbors in one
        // pass.  The loop test is at the end of the loop to avoid 32-bit overflow.
        var k = -nbrSize
        while (true) {
            val sameFace: Boolean = when {
                k < 0 -> j + k >= 0
                k >= size -> (j + k) < kMaxSize
                else -> {
                    // Top and bottom neighbors.
                    output.add(fromFaceIJSame(face, i + k, j - nbrSize, j - size >= 0).parent(nbrLevel))
                    output.add(fromFaceIJSame(face, i + k, j + size, j + size < kMaxSize).parent(nbrLevel))
                    true
                }
            }
            // Left, right, and diagonal neighbors.
            output.add(fromFaceIJSame(face, i - nbrSize, j + k, sameFace && i - size >= 0).parent(nbrLevel))
            output.add(fromFaceIJSame(face, i + size, j + k, sameFace && i + size < kMaxSize).parent(nbrLevel))
            if (k >= size) break
            k += nbrSize
        }
    }

    /**
     * Return the (face, i, j) coordinates for the leaf cell corresponding to this cell id. Since cells are represented
     * by the Hilbert curve position at the center of the cell, the returned (i,j) for non-leaf cells will be a leaf
     * cell adjacent to the cell center. If "orientation" is non-NULL, also return the Hilbert curve orientation for
     * the current cell.
     *
     * @param computeOrientation Indicates if the orientation must be computed.
     * @return The (face, i, j, orientation) coordinates of the cell.
     */
    @JvmOverloads
    fun toFaceIJOrientation(computeOrientation: Boolean = false): FaceIJ {
        var i = 0
        var j = 0
        val face = face()
        var bits: UInt = ((face and kSwapMask).toUInt())

        // Each iteration maps 8 bits of the Hilbert curve position into
        // 4 bits of "i" and "j".  The lookup table transforms a key of the
        // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
        // letters [ijpo] represents bits of "i", "j", the Hilbert curve
        // position, and the Hilbert curve orientation respectively.
        //
        // On the first iteration we need to be careful to clear out the bits
        // representing the cube face.
        for (k in 7 downTo 0) {
            val nbits: Int = if (k == 7) (kMaxLevel - 7 * kLookupBits) else kLookupBits
            bits += (((id shr (k * 2 * kLookupBits + 1)).toInt() and ((1 shl (2 * nbits)) - 1)) shl 2).toUInt()
            bits = lookupIj[bits.toInt()]
            i += ((bits shr (kLookupBits + 2)) shl (k * kLookupBits)).toInt()
            j += (((bits shr 2) and (((1 shl kLookupBits) - 1).toUInt())) shl (k * kLookupBits)).toInt()
            bits = bits and ((kSwapMask or kInvertMask).toUInt())
        }

        var orientation: Int? = null
        if (computeOrientation) {
            // The position of a non-leaf cell at level "n" consists of a prefix of
            // 2*n bits that identifies the cell, followed by a suffix of
            // 2*(kMaxLevel-n)+1 bits of the form 10*.  If n==kMaxLevel, the suffix is
            // just "1" and has no effect.  Otherwise, it consists of "10", followed
            // by (kMaxLevel-n-1) repetitions of "00", followed by "0".  The "10" has
            // no effect, while each occurrence of "00" has the effect of reversing
            // the kSwapMask bit.
            checkState { 0 == kPosToOrientation[2] }
            checkState { kSwapMask == kPosToOrientation[0] }
            if ((lsb() and 0x1111111111111110UL) != 0UL) {
                bits = bits xor kSwapMask.toUInt()
            }
            orientation = bits.toInt()
        }
        return FaceIJ(face = face, i = i, j = j, orientation = orientation)

    }

    // Return the lowest-numbered bit that is on for this cell id, which is
    // equal to (uint64{1} << (2 * (kMaxLevel - level))).  So for example,
    // a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the
    // first test is more efficient.
    internal fun lsb(): ULong {
        return id and (id.inv() + 1UL)
    }

    override fun equals(other: Any?): Boolean {
        if (other !is S2CellId) {
            return false
        }
        return id == other.id
    }

    override fun hashCode(): Int {
        return id.hashCode()
    }

    override fun compareTo(other: S2CellId): Int = id.compareTo(other.id)

    companion object {

        // The extra position bit (61 rather than 60) let us encode each cell as its
        // Hilbert curve position at the cell center (which is halfway along the
        // portion of the Hilbert curve that fills that cell).
        /** Number of bits used to encode the cell face. */
        const val kFaceBits = 3
        /** Number of faces. */
        const val kNumFaces = 6
        /** Max level. Valid levels: 0..kMaxLevel */
        const val kMaxLevel = 30

        const val kPosBits = 2 * kMaxLevel + 1

        const val kMaxSize = 1 shl kMaxLevel

        // Constant related to unsigned long's
        private const val kMaxUnsigned = ULong.MAX_VALUE // Equivalent to 0xffffffffffffffffL

        // This is the offset required to wrap around from the beginning of the
        // Hilbert curve to the end or vice versa; see next_Wrap() and prevWrap().
        private val kWrapOffset = kNumFaces.toULong() shl kPosBits

        /**
         * An invalid cellId instance.
         * The default constructor returns an invalid cell id.
         */
        @JvmStatic
        val none: S2CellId = S2CellId()

        /**
         * An invalid cell id guaranteed to be larger than any valid cell id.
         * Useful for creating indexes.
         */
        @JvmStatic
        val sentinel: S2CellId = S2CellId(kMaxUnsigned) // -1

        /**
         * Get the cell corresponding to a given S2 cube face.
         *
         * @param face A face id in 0..5
         * @return A cell id instance that represents the corresponding S2 cube face.
         */
        @JvmStatic
        fun fromFace(face: Int): S2CellId = S2CellId((face.toULong() shl kPosBits) + lsbForLevel(0))

        /**
         * Get a cell given its face (range 0..5), 61-bit Hilbert curve position within that face, and level
         * (range 0..kMaxLevel).
         * The given position will be modified to correspond to the Hilbert curve position at the center of the
         * returned cell.
         *
         * @param face The face of the cell.
         * @param pos The position of the cell on tge Hilbert curve.
         * @param level Level of the cell.
         * @return The new cell id instance.
         */
        @JvmStatic
        fun fromFacePosLevel(face: Int, pos: ULong, level: Int): S2CellId {
            val cell = S2CellId((face.toULong() shl kPosBits) + (pos or 1UL))
            return cell.parent(level)
        }

        /**
         * Get the edge length in (s,t)-space of cells at the given level.
         *
         * @param level A level in range 0..kMaxLevel
         * @return The edge length.
         */
        @JvmStatic
        fun getSizeST(level: Int): Double = S2Coords.ijToStMin(getSizeIJ(level))

        /**
         * Like getSizeST, but return the size of cells at the given level.
         *
         * @param level A level in range 0..kMaxLevel
         * @return The size of the cells.
         */
        fun getSizeIJ(level: Int): Int = 1 shl (kMaxLevel - level)

        /**
         * Expand a rectangle in (u,v)-space so that it contains all points within the given distance of the boundary,
         * and return the smallest such rectangle.  If the distance is negative, then instead shrink this rectangle so
         * that it excludes all points within the given absolute distance of the boundary.
         *
         * Distances are measured *on the sphere*, not in (u,v)-space.  For example, you can use this method to expand
         * the (u,v)-bound of an S2CellId so that it contains all points within 5km of the original cell.  You can then
         * test whether a point lies within the expanded bounds like this:
         *
         * val uv = R2Point()
         * if (S2Coords.faceXyzToUv(face, point, uv) && bound.contains(uv)) { ... }
         *
         * Limitations:
         *
         *  - Because the rectangle is drawn on one of the six cube-face planes (i.e., {x,y,z} = +/-1), it can cover at
         *    most one hemisphere.  This limits the maximum amount that a rectangle can be expanded.  For example,
         *    S2CellId bounds can be expanded safely by at most 45 degrees (about 5000 km on the Earth's surface).
         *
         *  - The implementation is not exact for negative distances.  The resulting rectangle will exclude all points
         *    within the given distance of the boundary but may be slightly smaller than necessary.
         *
         * @param uv The rectangle to expand.
         * @param distance The expansion distance onto the sphere surface.
         * @return The expanded rectangle.
        */
        @JvmStatic
        fun expandedByDistanceUV(uv: R2Rect, distance: S1Angle): R2Rect {
            // Expand each of the four sides of the rectangle just enough to include all
            // points within the given distance of that side.  (The rectangle may be
            // expanded by a different amount in (u,v)-space on each side.)
            val u0 = uv[0][0]
            val u1 = uv[0][1]
            val v0 = uv[1][0]
            val v1 = uv[1][1]
            val maxU = max(abs(u0), abs(u1))
            val maxV = max(abs(v0), abs(v1))
            val sinDist = sin(distance)
            return R2Rect(
                R1Interval(expandEndpoint(u0, maxV, -sinDist), expandEndpoint(u1, maxV, sinDist)),
                R1Interval(expandEndpoint(v0, maxU, -sinDist), expandEndpoint(v1, maxU, sinDist))
            )
        }

        // This is a helper function for ExpandedByDistanceUV().
        //
        // Given an edge of the form (u,v0)-(u,v1), let max_v = max(abs(v0), abs(v1)).
        // This method returns a new u-coordinate u' such that the distance from the
        // line u=u' to the given edge (u,v0)-(u,v1) is exactly the given distance
        // (which is specified as the sine of the angle corresponding to the distance).
        @Strictfp
        private fun expandEndpoint(u: Double, maxV: Double, sinDist: Double): Double {
            // This is based on solving a spherical right triangle, similar to the
            // calculation in S2Cap::GetRectBound.
            val sinUShift = sinDist * sqrt((1 + u * u + maxV * maxV) / (1 + u * u))
            val cosUShift = sqrt(1 - sinUShift * sinUShift)
            // The following is an expansion of tan(atan(u) + asin(sin_u_shift)).
            return (cosUShift * u + sinUShift) / (cosUShift - sinUShift * u)
        }


        // Iterator-style methods for traversing all the cells along the Hilbert
        // curve at a given level (across all 6 faces of the cube).  Note that the
        // end value is exclusive (just like standard STL dilivia.iterators), and is not a
        // valid cell id.
        @JvmStatic
        fun begin(level: Int): S2CellId = fromFace(0).childBegin(level)

        @JvmStatic
        fun end(level: Int): S2CellId = fromFace(5).childEnd(level)

        /**
         * Decodes the cell id from a compact text string suitable for display or
         * indexing. Cells at lower levels (i.e. larger cells) are encoded into
         * fewer characters. The maximum token length is 16.
         *
         * @param token the token to decode
         * @return the S2CellId for that token
         * @throws NumberFormatException if the token is not formatted correctly
         */
        @JvmStatic
        fun fromToken(token: String): S2CellId {
            if (token.length > 16) return none
            var id: ULong = 0UL
            var pos = 60
            for (i in token.indices) {
                val d: ULong = when {
                    token[i] in '0'..'9' -> (token[i] - '0').toULong()
                    token[i] in 'a'..'f' -> (token[i] - 'a' + 10).toULong()
                    token[i] in 'A'..'F' -> (token[i] - 'A' + 10).toULong()
                    else -> return none
                }
                id = id or (d shl pos)
                pos -= 4
            }
            return S2CellId(id)
        }

        // Converts a string in the format returned by ToString() to an S2CellId.
        // Returns S2CellId::None() if the string could not be parsed.
        //
        // The method name includes "Debug" in order to avoid possible confusion
        // with FromToken() above.
        @JvmStatic
        fun fromDebugString(str: String): S2CellId {
            // This function is reasonably efficient, but is only intended for use in
            // tests.
            val level = str.length - 2
            if (level < 0 || level > kMaxLevel) return none
            val face = str[0] - '0'
            if (face < 0 || face > 5 || str[1] != '/') return none
            var id = fromFace(face)
            for (i in 2 until str.length) {
                val childPos = str[i] - '0'
                if (childPos < 0 || childPos > 3) return none
                id = id.child(childPos)
            }
            return id
        }

        // ///////////////////////////////////////////////////////////////////
        // Low-level methods.
        /**
         * Return a leaf cell given its cube face (range 0..5) and i- and
         * j-coordinates (see s2.h).
         */
        @JvmStatic
        fun fromFaceIJ(faceIJ: FaceIJ): S2CellId = fromFaceIJ(faceIJ.face, faceIJ.i, faceIJ.j)

        @JvmStatic
        fun fromFaceIJ(face: Int, i: Int, j: Int): S2CellId {
            // Optimization notes:
            //  - Non-overlapping bit fields can be combined with either "+" or "|".
            //    Generally "+" seems to produce better code, but not always.

            // Note that this value gets shifted one bit to the left at the end
            // of the function.
            var n: ULong = face.toULong() shl (kPosBits - 1)

            // Alternating faces have opposite Hilbert curve orientations; this
            // is necessary in order for all faces to have a right-handed
            // coordinate system.
            var bits: ULong = (face and kSwapMask).toULong()

            // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
            // curve position.  The lookup table transforms a 10-bit key of the form
            // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
            // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
            // Hilbert curve orientation respectively.
            for (k in 7 downTo 0) {
                val mask: Int = (1 shl kLookupBits) - 1
                bits += (((i shr (k * kLookupBits)) and mask) shl (kLookupBits + 2)).toUInt()
                bits += (((j shr (k * kLookupBits)) and mask) shl 2).toUInt()
                bits = LookupCellTables.lookupPos[bits.toInt()].toULong()
                n = n or ((bits shr 2) shl (k * 2 * kLookupBits))
                bits = bits and ((kSwapMask or kInvertMask).toULong())
            }

            return S2CellId(n * 2UL + 1UL)
        }

        /**
         * Return the lowest-numbered bit that is on for this cell id, which is equal
         * to (uint64(1) << (2 * (MAX_LEVEL - level))). So for example, a.lsb() <=
         * b.lsb() if and only if a.level() >= b.level(), but the first test is more
         * efficient.
         */
        @JvmStatic
        fun lsbForLevel(level: Int): ULong {
            return 1UL shl (2 * (kMaxLevel - level))
        }

        // Return the bound in (u,v)-space for the cell at the given level containing
        // the leaf cell with the given (i,j)-coordinates.
        @JvmStatic
        fun ijLevelToBoundUV(i: Int, j: Int, level: Int): R2Rect {
            val cellSize = getSizeIJ(level)
            val iLo = i and -cellSize
            val jLo = j and -cellSize
            return R2Rect(
                x = R1Interval(
                    lo = S2Coords.stToUv(S2Coords.ijToStMin(iLo)),
                    hi = S2Coords.stToUv(S2Coords.ijToStMin(iLo + cellSize))
                ),
                y = R1Interval(
                    lo = S2Coords.stToUv(S2Coords.ijToStMin(jLo)),
                    hi = S2Coords.stToUv(S2Coords.ijToStMin(jLo + cellSize))
                )
            )
        }

        // Given a face and a point (i,j) where either i or j is outside the valid
        // range [0..kMaxSize-1], this function first determines which neighboring
        // face "contains" (i,j), and then returns the leaf cell on that face which
        // is adjacent to the given face and whose distance from (i,j) is minimal.
        const val kScale = 1.0 / kMaxSize
        val kLimit = 1.0 + DoubleType.epsilon
        private fun fromFaceIJWrap(face: Int, i: Int, j: Int): S2CellId {
            // Convert i and j to the coordinates of a leaf cell just beyond the
            // boundary of this face.  This prevents 32-bit overflow in the case
            // of finding the neighbors of a face cell.
            val normalizedI = max(-1, min(kMaxSize, i))
            val normalizedJ = max(-1, min(kMaxSize, j))

            // We want to wrap these coordinates onto the appropriate adjacent face.
            // The easiest way to do this is to convert the (i,j) coordinates to (x,y,z)
            // (which yields a point outside the normal face boundary), and then call
            // S2::XYZtoFaceUV() to project back onto the correct face.
            //
            // The code below converts (i,j) to (si,ti), and then (si,ti) to (u,v) using
            // the linear projection (u=2*s-1 and v=2*t-1).  (The code further below
            // converts back using the inverse projection, s=0.5*(u+1) and t=0.5*(v+1).
            // Any projection would work here, so we use the simplest.)  We also clamp
            // the (u,v) coordinates so that the point is barely outside the
            // [-1,1]x[-1,1] face rectangle, since otherwise the reprojection step
            // (which divides by the new z coordinate) might change the other
            // coordinates enough so that we end up in the wrong leaf cell.

            // The arithmetic below is designed to avoid 32-bit integer overflows.
            checkState { 0 == kMaxSize % 2 }
            val u = max(-kLimit, min(kLimit, kScale * (2 * (normalizedI - kMaxSize / 2) + 1)))
            val v = max(-kLimit, min(kLimit, kScale * (2 * (normalizedJ - kMaxSize / 2) + 1)))

            // Find the leaf cell coordinates on the adjacent face, and convert
            // them to a cell id at the appropriate level.
            val p = S2Coords.faceUvToXyz(face, u, v)
            val faceUV = S2Coords.xyzToFaceUv(p)
            val s = 0.5 * (faceUV.u + 1)
            val t = 0.5 * (faceUV.v + 1)
            return fromFaceIJ(
                face = faceUV.face,
                i = S2Coords.stToIj(s),
                j = S2Coords.stToIj(t)
            )
        }

        /**
         * Public helper function that calls FromFaceIJ if sameFace is true, or
         * FromFaceIJWrap if sameFace is false.
         */
        fun fromFaceIJSame(face: Int, i: Int, j: Int, sameFace: Boolean): S2CellId {
            return if (sameFace) {
                fromFaceIJ(face, i, j)
            } else {
                fromFaceIJWrap(face, i, j)
            }
        }


        /**
         * Construct a leaf cell containing the given point "p".  Usually there is
         * is exactly one such cell, but for points along the edge of a cell, any
         * adjacent cell may be (deterministically) chosen.  This is because
         * S2CellIds are considered to be closed sets.  The returned cell will
         * always contain the given point, i.e.
         *
         *   S2Cell(S2CellId(p)).Contains(p)
         *
         * is always true.  The point "p" does not need to be normalized.
         *
         * If instead you want every point to be contained by exactly one S2Cell,
         * you will need to convert the S2CellIds to S2Loops (which implement point
         * containment this way).
         */
        @JvmStatic
        fun fromPoint(p: S2Point): S2CellId {
            val (face, u, v) = S2Coords.xyzToFaceUv(p)
            val i = S2Coords.stToIj(S2Coords.uvToSt(u))
            val j = S2Coords.stToIj(S2Coords.uvToSt(v))
            return fromFaceIJ(face, i, j)
        }

        /** Return the leaf cell containing the given S2LatLng.  */
        @JvmStatic
        fun fromLatLng(ll: S2LatLng): S2CellId {
            return fromPoint(ll.toPoint())
        }

        // Print the num_digits low order hex digits.
        @JvmStatic
        fun hexFormatString(value: ULong, numDigits: Int): String {
            var result = ""
            var i = numDigits
            var v = value
            while (i >= 1) {
                result = "0123456789abcdef"[(v and 0xF.toULong()).toInt()] + result
                v = v shr 4
                i--
            }
            return result
        }

    }

}


