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
package dilivia.s2.coords

import dilivia.s2.coords.S2Coords.kPosToIJ

/**
 * The following lookup tables are used to convert efficiently between an
 * (i,j) cell index and the corresponding position along the Hilbert curve.
 * "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
 * orientation of the current cell into 8 bits representing the order in which
 * that subcell is visited by the Hilbert curve, plus 2 bits indicating the
 * new orientation of the Hilbert curve within that subcell. (Cell
 * orientations are represented as combination of kSwapMask and kInvertMask.)
 *
 * "lookup_ij" is an inverted table used for mapping in the opposite
 * direction.
 *
 * We also experimented with looking up 16 bits at a time (14 bits of position
 * plus 2 of orientation) but found that smaller lookup tables gave better
 * performance. (2KB fits easily in the primary cache.)
 * Values for these constants are *declared* in the *.h file. Even though
 * the declaration specifies a value for the constant, that declaration
 * is not a *definition* of storage for the value. Because the values are
 * supplied in the declaration, we don't need the values here. Failing to
 * define storage causes link errors for any code that tries to take the
 * address of one of these values.
 */
@ExperimentalUnsignedTypes
object LookupCellTables {

    const val kLookupBits = 4
    val lookupPos = UIntArray(1 shl 2 * kLookupBits + 2)
    val lookupIj = UIntArray(1 shl 2 * kLookupBits + 2)

    init {
        initLookupCell(0, 0, 0, 0, 0, 0)
        initLookupCell(0, 0, 0, S2Coords.kSwapMask, 0, S2Coords.kSwapMask)
        initLookupCell(0, 0, 0, S2Coords.kInvertMask, 0, S2Coords.kInvertMask)
        initLookupCell(0, 0, 0, S2Coords.kSwapMask or S2Coords.kInvertMask, 0, (S2Coords.kSwapMask or S2Coords.kInvertMask))
    }

    private fun initLookupCell(level: Int, i: Int, j: Int, origOrientation: Int, pos: Int, orientation: Int) {
        var currentLevel = level
        var currentI = i
        var currentJ = j
        var currentPos = pos
        if (currentLevel == kLookupBits) {
            val ij = (currentI shl kLookupBits) + currentJ
            lookupPos[(ij shl 2) + origOrientation] = ((currentPos shl 2) + orientation).toUInt()
            lookupIj[(currentPos shl 2) + origOrientation] = ((ij shl 2) + orientation).toUInt()
        } else {
            currentLevel++
            currentI = currentI shl 1
            currentJ = currentJ shl 1
            currentPos = currentPos shl 2
            // Initialize each sub-cell recursively.
            for (subPos in 0..3) {
                val ij = posToIJ(orientation, subPos)
                val orientationMask = posToOrientation(subPos)
                initLookupCell(currentLevel, currentI + (ij ushr 1), currentJ + (ij and 1), origOrientation,
                        currentPos + subPos, orientation xor orientationMask)
            }
        }
    }

    /**
     * Returns an XOR bit mask indicating how the orientation of a child subcell
     * is related to the orientation of its parent cell. The returned value can
     * be XOR'd with the parent cell's orientation to give the orientation of
     * the child cell.
     *
     * @param position the position of the subcell in the Hilbert traversal, in
     * the range [0,3].
     * @return a bit mask containing some combination of [.SWAP_MASK] and
     * [.INVERT_MASK].
     * @throws IllegalArgumentException if position is out of bounds.
     */
    private fun posToOrientation(position: Int): Int {
        assert(position in 0..3)
        return S2Coords.kPosToOrientation[position]
    }

    /**
     * Return the IJ-index of the subcell at the given position in the Hilbert
     * curve traversal with the given orientation. This is the inverse of
     * [.ijToPos].
     *
     * @param orientation the subcell orientation, in the range [0,3].
     * @param position the position of the subcell in the Hilbert traversal, in
     * the range [0,3].
     * @return the IJ-index where `0->(0,0), 1->(0,1), 2->(1,0), 3->(1,1)`.
     * @throws IllegalArgumentException if either parameter is out of bounds.
     */
    fun posToIJ(orientation: Int, position: Int): Int {
        assert(orientation in 0..3)
        assert(position in 0..3)
        return kPosToIJ[orientation].get(position)
    }

}
