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

import dilivia.s2.S2CellId
import dilivia.s2.index.lexicon.Label

// A function that is called with each (cell_id, label) pair to be visited.
// The function may return false in order to indicate that no further
// (cell_id, label) pairs are needed.
interface CellVisitor {
    fun visit(cellId: S2CellId, label: Label): Boolean
}
