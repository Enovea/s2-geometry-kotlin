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

interface Container<T> {
    fun size(): Int
    operator fun get(index: Int):T
}

fun <V, T: Comparable<V>> Container<T>.lowerBound(first: Int, last: Int, value: V): Int {
    var f = min(first, this.size())
    val l = min(last, this.size())
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
fun <V, T: Comparable<V>> Container<T>.upperBound(first: Int = 0, last: Int, value: V): Int {
    var f = min(first, this.size())
    val l = min(last, this.size())
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


