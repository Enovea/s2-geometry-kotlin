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

class S2CellIdIterator(val parent: S2CellId) : Iterator<S2CellId> {

    private var currentId: S2CellId = parent.rangeMin()
    private val max = parent.rangeMax()

    override fun hasNext(): Boolean = currentId < max

    override fun next(): S2CellId {
        val next = currentId
        currentId = currentId.next()
        return next
    }

}
