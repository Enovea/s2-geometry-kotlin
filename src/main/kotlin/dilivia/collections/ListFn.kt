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

import org.apache.commons.math3.util.FastMath.min


fun <T> List<T>.isSorted(): Boolean where T:Comparable<T>{
    if (this.size <= 1) return true
    val iter = this.iterator()
    var current: T
    var previous = iter.next()
    while (iter.hasNext()) {
        current = iter.next();
        if (previous > current) return false
        previous = current;
    }
    return true
}

fun <T> List<T>.isSorted(comparator: Comparator<T>): Boolean {
    if (this.size <= 1) return true
    val iter = this.iterator()
    var current: T
    var previous = iter.next()
    while (iter.hasNext()) {
        current = iter.next();
        if (comparator.compare(previous, current) > 0) return false
        previous = current;
    }
    return true
}


/**
 *
 */
fun <V, T: Comparable<V>> List<T>.lowerBound(first: Int = 0, last: Int = this.size, value: V): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        step = count / 2
        idx = f + step
        if (this[idx] < value) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}

/**
 *
 */
fun <T> List<T>.lowerBound(first: Int = 0, last: Int = this.size, value: T, comparator: Comparator<T>): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        step = count / 2
        idx = f + step
        if (comparator.compare(this[idx], value) < 0) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}

/**
 *
 */
fun <V, T: Comparable<V>> List<T>.upperBound(first: Int = 0, last: Int = this.size, value: V): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        idx = f
        step = count / 2
        idx += step
        if (this[idx] <= value) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}

/**
 *
 */
fun <T> List<T>.upperBound(first: Int = 0, last: Int = this.size, value: T, comp: Comparator<T>): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        idx = f
        step = count / 2
        idx += step
        if (comp.compare(this[idx], value) <= 0) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}

/**
 *
 */
fun <V, T: Comparable<V>> List<T>.equalRange(first: Int = 0, last: Int = this.size, value: V): Pair<Int, Int> = Pair(
        this.lowerBound(first, last, value),
        this.upperBound(first, last, value)
)


/**
 *
 */
fun <T> List<T>.equalRange(first: Int = 0, last: Int = this.size, value: T, comp: Comparator<T>): Pair<Int, Int> = Pair(
        this.lowerBound(first, last, value, comp),
        this.upperBound(first, last, value, comp)
)
