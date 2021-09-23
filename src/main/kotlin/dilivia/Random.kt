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

import kotlin.math.abs
import kotlin.random.Random
import kotlin.random.nextUInt
import kotlin.random.nextULong

/**
 *
 */
object Random {

    private var random = Random(1)

    // Reset the generator state using the given seed.
    @JvmStatic
    fun reset(seed: Int) {
        random = Random(seed)
    }

    // Return a uniformly distributed 64-bit unsigned integer.
    @JvmStatic
    fun randomULong(): ULong = random.nextULong()

    // Return a uniformly distributed 32-bit unsigned integer.
    @JvmStatic
    fun randomUInt(): UInt = random.nextUInt()

    @JvmStatic
    fun randomUInt(min: UInt, limit: UInt): UInt = random.nextUInt(min, limit)

    @JvmStatic
    fun randomDouble(): Double = random.nextDouble()

    // Return a uniformly distributed integer in the range [0,n).
    @JvmStatic
    fun randomInt(n: Int): Int = random.nextInt(n)

    @JvmStatic
    fun randomInt(min: Int, limit: Int): Int = random.nextInt(min, limit)

    @JvmStatic
    fun randomInt(): Int = random.nextInt()

    @JvmStatic
    fun randomLong(): Long = random.nextLong()

    @JvmStatic
    fun randomLong(min: Long, limit: Long): Long = random.nextLong(min, limit)

    @JvmStatic
    fun randomLong(n: Long): Long = random.nextLong(n)

    // Return a uniformly distributed "double" in the range [min, limit).
    @JvmStatic
    fun randomDouble(min: Double, limit: Double): Double = random.nextDouble(min, limit)

    // Return true with probability 1 in n.
    @JvmStatic
    fun oneIn(n: Int): Boolean = randomInt(n) == 0

    // Skewed: pick "base" uniformly from range [0,max_log] and then
    // return "base" random bits.  The effect is to pick a number in the
    // range [0,2^max_log-1] with bias towards smaller numbers.

    // Pick "base" uniformly from range [0,maxLog] and then return
    // "base" random bits. The effect is to pick a number in the range
    // [0,2^maxLog-1] with bias towards smaller numbers.
    @JvmStatic
    fun skewed(maxLog: Int): Int {
        val base = abs(random.nextInt()) % (maxLog + 1)
        // if (!base) return 0; // if 0==base, we & with 0 below.
        //
        // this distribution differs slightly from ACMRandom's Skewed,
        // since 0 occurs approximately 3 times more than 1 here, and
        // ACMRandom's Skewed never outputs 0.
        return random.nextInt() and (1 shl base) - 1
    }


}
