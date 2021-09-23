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
package dilivia.s2.builder

import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLE
import dilivia.PreConditions.requireNE
import dilivia.s2.index.SequenceLexicon

// IdSetLexicon is a class for compactly representing sets of non-negative
// integers such as array indices ("id sets").  It is especially suitable when
// either (1) there are many duplicate sets, or (2) there are many singleton
// or empty sets.  See also ValueLexicon and SequenceLexicon.
//
// Each distinct id set is mapped to a 32-bit integer.  Empty and singleton
// sets take up no additional space whatsoever; the set itself is represented
// by the unique id assigned to the set. Sets of size 2 or more occupy about
// 11 bytes per set plus 4 bytes per element (as compared to 24 bytes per set
// plus 4 bytes per element for std::vector).  Duplicate sets are
// automatically eliminated.  Note also that id sets are referred to using
// 32-bit integers rather than 64-bit pointers.
//
// This class is especially useful in conjunction with ValueLexicon<T>.  For
// example, suppose that you want to label objects with a set of strings.  You
// could use a ValueLexicon<string> to map the strings to "label ids" (32-bit
// integers), and then use IdSetLexicon to map each set of labels to a "label
// set id".  Each reference to that label set then takes up only 4 bytes.
//
// Example usage:
//
//   ValueLexicon<string> labels_;
//   IdSetLexicon label_sets_;
//
//   int32 GetLabelSet(const vector<string>& label_strings) {
//     vector<int32> label_ids;
//     for (const auto& str : label_strings) {
//       label_ids.push_back(labels_.Add(str));
//     }
//     return label_sets_.Add(label_ids);
//   }
//
//   int label_set_id = GetLabelSet(...);
//   for (auto id : label_sets_.id_set(label_set_id)) {
//     S2_LOG(INFO) << id;
//   }
//
// This class is similar to SequenceLexicon, except:
//
// 1. Empty and singleton sets are represented implicitly; they use no space.
// 2. Sets are represented rather than sequences; the ordering of values is
//    not important and duplicates are removed.
// 3. The values must be 32-bit non-negative integers (only).
class IdSetLexicon(
    private val idSets: SequenceLexicon<Int> = SequenceLexicon()
) : Cloneable {

    // This class represents a set of integers stored in the IdSetLexicon.
    class IdSet(val values: List<Int> = emptyList()): Iterable<Int> {

        constructor(v: Int): this(listOf(v))

        fun size(): Int = values.size

        override fun iterator(): Iterator<Int> = values.iterator()

    }

    private val tmp: MutableList<Int> = mutableListOf()  // temporary storage used during Add()


    // Clears all data from the lexicon.
    fun clear() {
        idSets.clear()
    }

    public override fun clone(): IdSetLexicon {
        return IdSetLexicon(idSets.clone())
    }

    // Add the given set of integers to the lexicon if it is not already
    // present, and return the unique id for this set.  "begin" and "end" are
    // forward dilivia.iterators over a sequence of values that can be converted to
    // non-negative 32-bit integers.  The values are automatically sorted and
    // duplicates are removed.  Returns a signed integer representing this set.
    //
    // REQUIRES: All values in [begin, end) are non-negative 32-bit integers.
    // Add the given set of integers to the lexicon if it is not already
    // present, and return the unique id for this set.  This is a convenience
    // method equivalent to Add(std::begin(container), std::end(container)).
    fun add(container: Iterable<Int>): Int {
        tmp.clear()
        container.forEach { value ->
            requireGE(value, 0)
            requireLE(value, Int.MAX_VALUE)
            tmp.add(value)
        }
        return addInternal(tmp)
    }

    // Convenience method that returns the unique id for a singleton set.
    // Note that because singleton sets take up no space, this method is
    // const.  Equivalent to calling Add(&id, &id + 1).
    fun addSingleton(id: Int): Int {
        requireGE(id, 0)
        requireLE(id, Int.MAX_VALUE)
        // Singleton sets are represented by their element.
        return id;
    }

    // Return the set of integers corresponding to an id returned by Add().
    fun idSet(id: Int): IdSet {
        if (id >= 0) {
            return IdSet(id)
        } else if (id == kEmptySetId) {
            return IdSet()
        } else {
            val sequence = idSets.sequence(id.inv())
            requireNE(0, sequence.size())
            return IdSet(sequence.values())
        }
    }

    private fun addInternal(ids: List<Int>): Int {
        if (ids.isEmpty()) {
            // Empty sets have a special id chosen not to conflict with other ids.
            return kEmptySetId;
        } else if (ids.size == 1) {
            // Singleton sets are represented by their element.
            return ids[0]
        } else {
            // Canonicalize the set by sorting and removing duplicates.
            // Non-singleton sets are represented by the bitwise complement of the id
            // returned by SequenceLexicon.
            return idSets.add(tmp.distinct().sorted()).inv();
        }
    }

    override fun toString(): String {
        return "$idSets"
    }

    companion object {

        // Choose kEmptySetId to be the last id that will ever be generated.
        // (Non-negative ids are reserved for singleton sets.)
        val kEmptySetId = Int.MIN_VALUE

        // Convenience method that returns the unique id for the empty set.  Note
        // that because the empty set takes up no space and has a fixed id, this
        // method is static.  Equivalent to calling Add() with an empty container.
        fun emptySetId(): Int = kEmptySetId

    }
}
