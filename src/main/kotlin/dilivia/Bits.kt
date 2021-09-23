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
package dilivia

object Bits {

    fun findLSBSetNonZero(n: UInt): Int {
        var v = n
        var rc = 31
        var shift = 1 shl 4
        var i = 4
        while(i >= 0) {
            val x = v shl shift
            if (x != 0.toUInt()) {
                v = x
                rc -= shift
            }
            shift = shift shr 1
            --i
        }
        return rc
    }

    fun findLSBSetNonZero64(n: ULong): Int {
        val bottombits = n.toUInt()
        if (bottombits == 0U) {
            // Bottom bits are zero, so scan in top bits
            return 32 + findLSBSetNonZero((n shr 32).toUInt())
        } else {
            return findLSBSetNonZero(bottombits)
        }
    }

    fun findMSBSetNonZero(n: UInt): Int { return log2FloorNonZero(n) }

    fun findMSBSetNonZero64(n: ULong): Int { return log2FloorNonZero64(n) }

    // Log2FloorNonZero64() is defined in terms of Log2FloorNonZero32()
    fun log2FloorNonZero64(n: ULong): Int {
        val topbits: UInt = (n shr 32).toUInt()
        if (topbits == 0U) {
            // Top bits are zero, so scan in bottom bits
            return log2FloorNonZero(n.toUInt())
        } else {
            return 32 + log2FloorNonZero(topbits)
        }
    }

    fun log2FloorNonZero(n: UInt): Int {
        // Just use the common routine
        return log2Floor(n)
    }

    fun log2Floor(n: UInt): Int {
        if (n == 0U)
            return -1
        var log = 0
        var value: UInt = n
        var i = 4
        while (i >= 0) {
            val shift = (1 shl i)
            val x: UInt = value shr shift
            if (x != 0U) {
                value = x
                log += shift
            }
            --i
        }
        assert(value == 1U)
        return log
    }

}