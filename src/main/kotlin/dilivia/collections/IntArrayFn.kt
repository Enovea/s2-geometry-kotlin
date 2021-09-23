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


/**
 * Gets the index of the first element in index range [first, last) that is not less (i.e. greater or equal to) value,
 * org last if no such element is found.
 *
 * @param first The first index of the research range.
 * @param last The last index (exclusive) of the research range.
 * @param value The value to compare the elements to.
 * @return The index of the first not less element.
 */
fun IntArray.lowerBound(first: Int = 0, last: Int = this.size, value: Int): Int {
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
 * Gets the index of the first element in range [first, last) that is greater than value, or last if no such element
 * is found.
 *
 * @param first The first index of the research range.
 * @param last The last index (exclusive) of the research range.
 * @param value The value to compare the elements to.
 * @return The index of the first greater element.
 */
fun IntArray.upperBound(first: Int, last: Int, value: Int): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        idx = f
        step = count / 2
        idx += step
        if (value >= this[idx]) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}
