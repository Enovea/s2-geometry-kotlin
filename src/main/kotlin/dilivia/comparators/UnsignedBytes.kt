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
package dilivia.comparators

import kotlin.math.min

object UnsignedBytes: Comparator<Byte> {

    fun toInt(value: Byte): Int {
        return value.toInt() and 255
    }
    override fun compare(a: Byte, b: Byte): Int = toInt(a) - toInt(b)

    object LexicographicalComparator : Comparator<ByteArray> {

        override fun compare(left: ByteArray,  right: ByteArray): Int {
            val minLength: Int = min(left.size, right.size)

            for (i in 0 until minLength) {
                val result: Int = compare(left[i], right[i])
                if (result != 0) {
                    return result
                }
            }

            return left.size - right.size
        }

    }
}