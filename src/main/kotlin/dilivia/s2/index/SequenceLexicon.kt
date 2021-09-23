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
package dilivia.s2.index

import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantReadWriteLock
import kotlin.concurrent.withLock

/**
 * SequenceLexicon is a class for compactly representing sequences of values (e.g., tuples). It automatically eliminates
 * duplicates, and maps the remaining sequences to sequentially increasing integer ids. See also ValueLexicon and
 * IdSetLexicon.
 *
 * This class is thread-safety.
 *
 * Example usage:
 *
 *   val lexicon = SequenceLexicon<String>()
 *   val pets = listOf<String>("cat", "dog", "parrot")
 *   val petsId = lexicon.add(pets)
 *   assertThat(lexicon.add(pets)).isEqualTo(petsId)
 *   var values = ""
 *   for (pet in lexicon.sequence(petsId)) {
 *     values += pet
 *   }
 *   assertThat(values).isEqualTo("catdogparrot")
 *
 * This class is a port of the SequenceLexicon class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class SequenceLexicon<T>(
    private val values: MutableList<T> = mutableListOf(),
    private val begins: MutableList<Int> = mutableListOf(0),
    private val idSet: MutableMap<SequenceLexicon<T>.SequenceKey, Int> = mutableMapOf()
) : Cloneable {

    /** Lock. */
    private val lock = ReentrantReadWriteLock()
    private val writeLock: Lock
        get() = lock.writeLock()
    private val readLock: Lock
        get() = lock.readLock()


    override fun toString(): String {
        return """
            |values = $values
            |begins = $begins
            |idSet = $idSet
            |""".trimMargin()
    }

    // Clears all data from the lexicon.
    fun clear() {
        writeLock.withLock {
            values.clear()
            begins.clear()
            idSet.clear()
            begins.add(0)
        }
    }

    public override fun clone(): SequenceLexicon<T> {
        readLock.withLock {
            val clone = SequenceLexicon<T>(values.toMutableList(), begins.toMutableList())
            for (entry in idSet) {
                clone.idSet[clone.SequenceKey(entry.key.id)] = entry.value
            }
                return clone
        }
    }

    // Add the given sequence of values to the lexicon if it is not already
    // present, and return its integer id.  Ids are assigned sequentially
    // starting from zero.  This is a convenience method equivalent to
    // Add(std::begin(container), std::end(container)).
    fun add(container: Iterable<T>): Int {
        writeLock.withLock {
            values.addAll(container)
            begins.add(values.size)
            val sequenceKey = SequenceKey(begins.size - 2)
            var id = idSet[sequenceKey]
            if (id == null) {
                id = begins.size - 2
                idSet[sequenceKey] = id
                return id
            } else {
                begins.removeLast()
                while (values.size > begins.last()) values.removeLast()
                return id
            }
        }
    }

    // Return the number of value sequences in the lexicon.
    fun size(): Int = readLock.withLock { begins.size - 1 }

    // Return the value sequence with the given id.  This method can be used
    // with range-based for loops as follows:
    //   for (const auto& value : lexicon.sequence(id)) { ... }
    fun sequence(id: Int): Sequence = readLock.withLock { Sequence(begins[id], begins[id + 1]) }

    // A class representing a sequence of values.
    inner class Sequence(val beginIdx: Int, val endIdx: Int) {

        fun begin(): T = readLock.withLock { values[beginIdx] }
        fun end(): T = readLock.withLock { values[endIdx] }
        fun size() = endIdx - beginIdx

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as SequenceLexicon<*>.Sequence

            return values() == other.values()
        }

        override fun hashCode(): Int = readLock.withLock { values().hashCode() }

        fun values(): List<T> = readLock.withLock { values.slice(beginIdx until endIdx) }

    }

    inner class SequenceKey(val id: Int) {

        fun sequence(): Sequence = sequence(id)

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as SequenceLexicon<*>.SequenceKey

            return sequence() == other.sequence()
        }

        override fun hashCode(): Int {
            return sequence().hashCode()
        }

        override fun toString(): String {
            return "SequenceKey($id)"
        }

    }
}

