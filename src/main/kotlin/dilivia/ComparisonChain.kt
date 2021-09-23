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

import java.lang.Boolean.compare as booleanCompare
import java.lang.Integer.compare as integerCompare
import java.lang.Long.compare as longCompare

/**
 */
abstract class ComparisonChain private constructor() {
    private class InactiveComparisonChain internal constructor(val result: Int) : ComparisonChain() {
        override fun <T> compare(left: Comparable<T>, right: T): ComparisonChain {
            return this
        }

        override fun <T> compare(left: T, right: T, comparator: Comparator<T>?): ComparisonChain {
            return this
        }

        override fun compare(left: Int, right: Int): ComparisonChain {
            return this
        }

        override fun compare(left: Long, right: Long): ComparisonChain {
            return this
        }

        override fun compare(left: Float, right: Float): ComparisonChain {
            return this
        }

        override fun compare(left: Double, right: Double): ComparisonChain {
            return this
        }

        override fun compareTrueFirst(left: Boolean, right: Boolean): ComparisonChain {
            return this
        }

        override fun compareFalseFirst(left: Boolean, right: Boolean): ComparisonChain {
            return this
        }

        override fun result(): Int {
            return result
        }
    }

    /**
     * Compares two comparable objects as specified by [Comparable.compareTo], *if* the
     * result of this comparison chain has not already been determined.
     */
    abstract fun <T> compare(left: Comparable<T>, right: T): ComparisonChain

    /**
     * Compares two objects using a comparator, *if* the result of this comparison chain has not
     * already been determined.
     */
    abstract fun <T> compare(left: T, right: T, comparator: Comparator<T>?): ComparisonChain

    /**
     * Compares two `int` values as specified by [Ints.compare], *if* the result of
     * this comparison chain has not already been determined.
     */
    abstract fun compare(left: Int, right: Int): ComparisonChain

    /**
     * Compares two `long` values as specified by [Longs.compare], *if* the result of
     * this comparison chain has not already been determined.
     */
    abstract fun compare(left: Long, right: Long): ComparisonChain

    /**
     * Compares two `float` values as specified by [Float.compare], *if* the result
     * of this comparison chain has not already been determined.
     */
    abstract fun compare(left: Float, right: Float): ComparisonChain

    /**
     * Compares two `double` values as specified by [Double.compare], *if* the result
     * of this comparison chain has not already been determined.
     */
    abstract fun compare(left: Double, right: Double): ComparisonChain

    /**
     * Compares two `boolean` values, considering `true` to be less than `false`,
     * *if* the result of this comparison chain has not already been determined.
     *
     */
    abstract fun compareTrueFirst(left: Boolean, right: Boolean): ComparisonChain

    /**
     * Compares two `boolean` values, considering `false` to be less than `true`,
     * *if* the result of this comparison chain has not already been determined.
     *
     */
    abstract fun compareFalseFirst(left: Boolean, right: Boolean): ComparisonChain

    /**
     * Ends this comparison chain and returns its result: a value having the same sign as the first
     * nonzero comparison result in the chain, or zero if every result was zero.
     */
    abstract fun result(): Int

    companion object {
        /** Begins a new chained comparison statement. See example in the class documentation.  */
        fun start(): ComparisonChain {
            return ACTIVE
        }

        private val ACTIVE: ComparisonChain = object : ComparisonChain() {

            override fun <T> compare(left: Comparable<T>, right: T): ComparisonChain {
                return classify(left.compareTo(right))
            }


            override fun <T> compare(left: T, right: T, comparator: Comparator<T>?): ComparisonChain {
                return classify(comparator!!.compare(left, right))
            }

            override fun compare(left: Int, right: Int): ComparisonChain {
                return classify(integerCompare(left, right))
            }

            override fun compare(left: Long, right: Long): ComparisonChain {
                return classify(longCompare(left, right))
            }

            override fun compare(left: Float, right: Float): ComparisonChain {
                return classify(java.lang.Float.compare(left, right))
            }

            override fun compare(left: Double, right: Double): ComparisonChain {
                return classify(java.lang.Double.compare(left, right))
            }

            override fun compareTrueFirst(left: Boolean, right: Boolean): ComparisonChain {
                return classify(booleanCompare(right, left)) // reversed
            }

            override fun compareFalseFirst(left: Boolean, right: Boolean): ComparisonChain {
                return classify(booleanCompare(left, right))
            }

            fun classify(result: Int): ComparisonChain {
                return if (result < 0) LESS else if (result > 0) GREATER else this
            }

            override fun result(): Int {
                return 0
            }
        }
        private val LESS: ComparisonChain = InactiveComparisonChain(-1)
        private val GREATER: ComparisonChain = InactiveComparisonChain(1)
    }
}
