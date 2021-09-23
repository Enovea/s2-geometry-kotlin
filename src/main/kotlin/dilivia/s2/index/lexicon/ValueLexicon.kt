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
package dilivia.s2.index.lexicon

/**
 * Labels are 32-bit non-negative integers.  To support other label types, you can use ValueLexicon to map label values
 * to integers:
 *
 * val myLabelLexicon = ValueLexicon<MyLabel>()
 * index.add(cellId, myLabelLexicon.add(label))
 */
typealias Label = Int

interface ValueLexicon<T> {

    fun clear()

    fun add(value: T): Int

    fun size(): Int

    fun value(id: Int): T

}

