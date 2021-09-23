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
package dilivia.collections

import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLT
import org.apache.commons.math3.util.FastMath.min
import java.util.function.Predicate

fun <T> MutableList<T>.removeFirstThat(predicate: Predicate<T>): T? {
    val iterator = this.iterator()
    var removedValue: T? = null
    while (removedValue == null && iterator.hasNext()) {
        val value = iterator.next()
        if (predicate.test(value)) {
            iterator.remove()
            removedValue = value
        }
    }
    return removedValue
}

fun MutableList<Int>.iota(startIdx: Int, endIdx: Int, value: Int) {
    var current = value
    var index = startIdx
    while (index < endIdx) {
        this[index++] = current
        ++current
    }
}

fun <T> MutableList<T>.remove(fromIndex: Int, toIndex: Int) {
    require(fromIndex >= 0 && fromIndex <= this.lastIndex)
    require(toIndex >= fromIndex && toIndex <= this.size)
    val nbElementToRemove = toIndex - fromIndex
    repeat(nbElementToRemove) { this.removeAt(fromIndex) }
}

fun <T> MutableList<T>.assignWith(size: Int, allocator: () -> T) {
    if (this is ArrayList) {
        this.ensureCapacity(size)
    }
    (0 until min(this.size, size)).forEach { i -> this[i] = allocator() }
    while (this.size < size) this.add(allocator())
    while (this.size > size) this.removeLast()
}

fun <T> MutableList<T>.assign(size: Int, value: T) {
    if (this is ArrayList) {
        this.ensureCapacity(size)
    }
    (0 until min(this.size, size)).forEach { i -> this[i] = value }
    while (this.size < size) this.add(value)
    while (this.size > size) this.removeLast()
}


fun <T : Comparable<T>> MutableList<T>.sortAndRemoveDuplicates() {
    this.sort()
    var idx = this.lastIndex
    while (idx >= 1) {
        if (this[idx] == this[idx - 1]) this.removeAt(idx)
        --idx
    }
}

fun <T> MutableList<T>.sortAndRemoveDuplicatesWith(comparator: Comparator<T>) {
    this.sortWith(comparator)
    var idx = this.lastIndex
    while (idx >= 1) {
        if (this[idx] == this[idx - 1]) this.removeAt(idx)
        --idx
    }
}

fun <T> MutableList<T>.removeAll(value: T) {
    var idx = this.lastIndex
    while (idx >= 0) {
        if (this[idx] == value) this.removeAt(idx)
        --idx
    }
}

fun <T> MutableList<T>.reverse(startIdx: Int, endIdx: Int) {
    this.subList(startIdx, endIdx).reverse()
}

fun <T> MutableList<T>.rotate(startIdx: Int, newStartIdx: Int, endIdx: Int): Int {
    if (startIdx == newStartIdx) return endIdx
    if (newStartIdx == endIdx) return startIdx

    var read = newStartIdx
    var write = startIdx
    var nextRead = startIdx // read position for when "read" hits "last"

    while (read != endIdx) {
        if (write == nextRead) nextRead = read // track where "first" went
        val tmp = this[read]
        this[read] = this[write]
        this[write] = tmp
        ++read;++write
    }

    // rotate the remaining sequence into place
    this.rotate(write, nextRead, endIdx)
    return write
}


fun <T> MutableList<T>.resize(size: Int, value: T) {
    while (this.size > size) this.removeLast()
    while(this.size < size) this.add(value)
}


fun <T> MutableList<T>.resizeWith(size: Int, allocator: () -> T) {
    while (this.size > size) this.removeLast()
    while(this.size < size) this.add(allocator.invoke())
}


@JvmName("resizeInt")
fun MutableList<Int>.resize(size: Int) {
    this.resize(size, 0)
}

@JvmName("resizeBoolean")
fun MutableList<Boolean>.resize(size: Int) {
    this.resize(size, false)
}

fun <T> MutableList<T>.erase(startIdx: Int, endIdx: Int) {
    repeat(endIdx - startIdx) { this.removeAt(startIdx) }
}

fun <T> MutableList<T>.swap(i1: Int, i2: Int) {
    requireGE(i1, 0)
    requireLT(i1, this.size)
    requireGE(i2, 0)
    requireLT(i2, this.size)
    val tmp = this[i1]
    this[i1] = this[i2]
    this[i2] = tmp
}
