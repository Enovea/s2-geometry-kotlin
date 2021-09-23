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

import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantReadWriteLock
import kotlin.concurrent.withLock

/**
 * ValueLexicon is a class that maps distinct values to sequentially numbered integer identifiers. It automatically
 * eliminates duplicates.  See also SequenceLexicon.
 *
 * This class is threat-safe.
 *
 * Example usage:
 *
 *  val lexicon = ValueLexicon<String>()
 *  val catId = lexicon.add("cat");
 *  assertThat(lexicon.add("cat")).isEqualTo(cat_id)
 *  assertThat(lexicon.value(cat_id)).isEqualTo("cat")
 *
 * This class is a port of the ValueLexicon class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @see SequenceLexicon
 * @since 1.0
 */
class InMemoryValueLexicon<T> : ValueLexicon<T> {

    /** Lock. */
    private val lock = ReentrantReadWriteLock()
    private val writeLock: Lock
        get() = lock.writeLock()
    private val readLock: Lock
        get() = lock.readLock()

    /** Values in the lexicon. */
    private val values: MutableList<T> = mutableListOf()

    /** ValueKey map to value index. */
    private val idSet: MutableMap<ValueKey, Int> = mutableMapOf()

    // Clears all data from the lexicon.
    override fun clear() {
        lock.writeLock().lock()
        values.clear()
        idSet.clear()
        lock.writeLock().unlock()
    }

    // Add the given value to the lexicon if it is not already present, and
    // return its integer id.  Ids are assigned sequentially starting from zero.
    override fun add(value: T): Int {
        readLock.withLock {
            if (values.isNotEmpty() && value == values.last()) {
                return values.lastIndex
            }
        }

        writeLock.withLock {
            values.add(value)
            val valueKey = ValueKey(values.lastIndex)
            val id = idSet[valueKey]
            return if (id == null) {
                idSet[valueKey] = values.lastIndex
                values.lastIndex
            } else {
                values.removeLast()
                id
            }
        }

    }

    // Return the number of values in the lexicon.
    override fun size(): Int = readLock.withLock { values.size }

    // Return the value with the given id.
    override fun value(id: Int) = readLock.withLock { values[id] }

    private inner class ValueKey(val id: Int) {

        fun value(): T = readLock.withLock { values[id] }

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as InMemoryValueLexicon<*>.ValueKey

            val value = value()
            val otherValue = other.value() as T

            if (value != otherValue) return false

            return true
        }

        override fun hashCode(): Int = readLock.withLock { value().hashCode() }

    }
}
