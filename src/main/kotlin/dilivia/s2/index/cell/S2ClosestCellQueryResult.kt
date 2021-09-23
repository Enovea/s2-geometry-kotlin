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
import dilivia.s2.S2CellId
import dilivia.s2.index.Distance
import dilivia.s2.index.lexicon.Label

/**
 * Result a result element of a S2ClosestCellQueryBase.
 * Each "Result" object represents a closest (cell_id, value) pair.
 *
 * @property distance The distance from the result value to the query target.
 * @property cellId The pair cell id.
 * @property value The pair value
 *
 * @param T Distance used to search cells.
 * @param V The type of data stored by the index.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
data class S2ClosestCellQueryResult<T : Distance<T>>(
    val distance: T,
    val cellId: S2CellId = S2CellId.none,
    val value: Label = -1
) : Comparable<S2ClosestCellQueryResult<T>> {

    /**
     * Indicates if the result is empty.
     *
     * @return true if this Result object does not refer to any cell. (The only case where an empty Result is returned
     * is when the findClosestCell() method does not find any cells that meet the specified criteria.)
     */
    fun isEmpty(): Boolean = cellId == S2CellId.none

    override fun compareTo(other: S2ClosestCellQueryResult<T>): Int {
        return ComparisonChain.start()
            .compare(distance, other.distance)
            .compare(cellId, other.cellId)
            .compare(value, other.value)
            .result()
    }

}
