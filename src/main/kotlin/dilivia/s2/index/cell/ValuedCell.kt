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

// Convenience class that represents a (cell_id, label) pair.
data class ValuedCell<T : Comparable<T>>(
    val cellId: S2CellId = S2CellId.none,
    val label: T? = null // -1
) : Comparable<ValuedCell<T>> {

    override fun compareTo(other: ValuedCell<T>): Int = ComparisonChain.start()
        .compare(cellId, other.cellId)
        .compare(label, other.label) { o1, o2 ->
            when {
                o1 == o2 -> 0
                o1 == null && o2 != null -> -1
                o1 != null && o2 == null -> 1
                else -> o1!!.compareTo(o2!!)
            }
        }
        .result()

}
