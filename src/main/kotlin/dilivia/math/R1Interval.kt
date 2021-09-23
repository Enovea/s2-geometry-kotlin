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
package dilivia.math

import dilivia.PreConditions.checkState
import dilivia.math.vectors.R2VectorDouble
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min

/**
 * An R1Interval represents a closed, bounded interval on the real line. It is capable of representing the empty
 * interval (containing no points) and zero-length intervals (containing a single point).
 *
 * This class is a port of the R1Interval class of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @since 1.0
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 */
class R1Interval : Cloneable {

    /** The low bound of the interval. */
    var lo: Double

    /** The high bound of the interval. */
    var hi: Double

    /**
     * Interval constructor. If lo > hi, the interval is empty.
     *
     * @param lo The low bound of the interval.
     * @param hi The high bound of the interval. If hi < lo, the interval is empty.
     */
    constructor(lo: Double, hi: Double) {
        this.lo = lo
        this.hi = hi
    }

    /**
     * Builds a R1Interval instance from integer bounds.
     *
     * @param lo The low bound of the interval
     * @param hi The high bound of the interval
     */
    constructor(lo: Int, hi: Int) : this(lo.toDouble(), hi.toDouble())

    /**
     * The default constructor creates an empty interval. (Any interval where lo > hi is considered to be empty.)
     */
    constructor() : this(1, 0)

    public override fun clone(): R1Interval = R1Interval(lo, hi)

    /**
     * Gets the interval bound lo or hi as an array element.
     *
     * @param i Index of the bound (0 for lo and 1 for hi)
     * @return The interval bound at the given index.
     * @throws ArrayIndexOutOfBoundsException If the parameter i is not in range 0..1
     */
    @Throws(ArrayIndexOutOfBoundsException::class)
    operator fun get(i: Int) = when (i) {
        0 -> lo
        1 -> hi
        else -> throw ArrayIndexOutOfBoundsException("Index $i is out of bounds 0..1")
    }

    operator fun set(i: Int, v: Double) {
        when(i) {
            0 -> lo = v
            1 -> hi = v
            else -> throw ArrayIndexOutOfBoundsException("Index $i is out of bounds 0..1")
        }
    }

    fun set(other: R1Interval) {
        this.lo = other.lo
        this.hi = other.hi
    }

    /**
     * The interval bounds as an array of doubles.
     */
    var bounds: R2VectorDouble
        get() = R2VectorDouble(lo, hi)
        set(value) {
            lo = value[0]
            hi = value[1]
        }

    /**
     * Indicates if the interval is empty, i.e. it contains no points.
     */
    val isEmpty: Boolean
        get() = lo > hi

    /**
     * Indicates if the interval contains at least a point, i.e. it is not empty.
     */
    val isNotEmpty: Boolean
        get() = !isEmpty

    /**
     * The center of the interval. For empty intervals, the result is arbitrary.
     */
    val center: Double
        get() = 0.5 * (lo + hi)

    /**
     * The length of the interval. The length of an empty interval is negative.
     */
    val length: Double
        get() = hi - lo

    /**
     * Checks if the interval contains a given point.
     *
     * @param p A real point.
     * @return true if the point p is in the range lo..hi
     */
    operator fun contains(p: Double): Boolean = p in lo..hi

    /**
     * Checks if the interior (bounds excluded) of the interval containt a given point.
     *
     * @param p A point
     * @return true if the point p is > lo and < hi.
     */
    fun interiorContains(p: Double): Boolean = p > lo && p < hi

    /**
     * Check if this interval contains a given interval, i.e. all the point of the interval y are contained by this
     * interval.
     *
     * @param y An interval.
     * @return true if this interval contains the interval 'y'.
     */
    operator fun contains(y: R1Interval): Boolean = y.isEmpty || (y.lo >= lo && y.hi <= hi)

    /**
     * Checks if the interior of this interval contains a given interval.
     *
     * @param y an interval.
     * @return true if the interior of this interval contains the entire interval 'y' (including its boundary).
     */
    fun interiorContains(y: R1Interval): Boolean = y.isEmpty || (y.lo > lo && y.hi < hi)

    /**
     * Checks the a given interval intersects this one.
     *
     * @param y an interval.
     * @return true if this interval intersects the given interval, i.e. if they have any points in common.
     */
    fun intersects(y: R1Interval): Boolean {
        return if (lo <= y.lo) {
            y.lo <= hi && y.lo <= y.hi
        } else {
            lo <= y.hi && lo <= hi
        }
    }

    /**
     * Checks if the interior of this interval intersects an other one.
     *
     * @param y an interval.
     * @return true if the interior of this interval intersects any point of the given interval (including its boundary).
     */
    fun interiorIntersects(y: R1Interval): Boolean {
        return y.lo < hi && lo < y.hi && lo < hi && y.lo <= y.hi
    }

    /**
     * Computes the Hausdorff distance from this interval to another one.
     * For two R1Intervals x and y, this distance is defined as
     *     h(x, y) = max_{p in x} min_{q in y} d(p, q).
     * @param y an interval.
     * @return the Hausdorff distance to the given interval 'y'.
     */
    fun directedHausdorffDistance(y: R1Interval): Double {
        if (isEmpty) return 0.0
        if (y.isEmpty) return Double.MAX_VALUE
        return max(0.0, max(hi - y.hi, y.lo - lo))
    }

    /**
     * Expand the interval so that it contains the given point "p". Modify this interval such that it becomes the
     * minimum interval that contains this interval and the given point.
     *
     * @param p The point to add.
     * @return this.
     */
    fun addPoint(p: Double): R1Interval {
        when {
            isEmpty -> {
                lo = p
                hi = p
            }
            p < lo -> lo = p
            p > hi -> hi = p
        }
        return this
    }

    /**
     * Expand the interval so that it contains the given interval "y" that is to say the minimum interval that contains
     * this interval and the interval 'y'.
     *
     * @param y the interval to add.
     * @return this.
     */
    fun addInterval(y: R1Interval): R1Interval {
        if (y.isNotEmpty) {
            if (isEmpty) {
                lo = y.lo
                hi = y.hi
            }
            else {
                lo = min(lo, y.lo)
                hi = max(hi, y.hi)
            }
        }
        return this
    }

    /**
     * Project a point to this interval.
     *
     * @param p A point
     * @return the closest point in the interval to the given point "p". The interval must be non-empty.
     * @throws IllegalStateException if this interval is empty.
     */
    @Throws(IllegalStateException::class)
    fun project(p: Double): Double {
        checkState({isNotEmpty}, {"Trying to project point $p on an empty interval."})
        check(isNotEmpty)
        return max(lo, min(hi, p))
    }

    /**
     * Get a new interval instance that represents this interval expanded with the given margin.
     * Note that the expansion of an empty interval is always empty.
     *
     * @param margin The margin.
     * @return an interval that contains all points with a distance "margin" of a point in this interval.
     */
    fun expanded(margin: Double): R1Interval {
        return if (isEmpty) {
            this
        } else R1Interval(lo - margin, hi + margin)
    }

    /**
     * Expands this interval on each side by the given distance "margin".  If "margin" is negative, then shrink the
     * interval on each side by "margin" instead.  The resulting interval may be empty. Any expansion of an empty
     * interval remains empty.
     *
     * @param margin The margin.
     * @return this.
     */
    fun expands(margin: Double): R1Interval {
        if (isNotEmpty) {
            lo -= margin
            hi += margin
        }
        return this
    }

    /**
     * Compute the union of this interval and another one.
     *
     * @param y an interval.
     * @return the smallest interval that contains this interval and the given interval "y".
     */
    fun union(y: R1Interval): R1Interval = when {
        isEmpty -> y
        y.isEmpty -> this
        else -> R1Interval(min(lo, y.lo), max(hi, y.hi))
    }

    /**
     * Updates this interval with the union of this interval and "y".
     *
     * @param y an interval.
     * @return this.
     */
    fun setUnion(y: R1Interval): R1Interval {
        if (isEmpty) {
            lo = y.lo
            hi = y.hi
        }
        else if (y.isNotEmpty) {
            lo = min(lo, y.lo)
            hi = max(hi, y.hi)
        }
        return this
    }

    /**
     * Compute the intersection of this interval and another one.
     *
     * @param y An interval.
     * @return the intersection of this interval with the given interval.
     */
    //  Empty intervals do not need to be special-cased.
    fun intersection(y: R1Interval): R1Interval = R1Interval(max(lo, y.lo), min(hi, y.hi))

    /**
     * Updates this interval with the intersection of this interval and "y".
     *
     * @param y an interval.
     * @return this.
     */
    fun setIntersection(y: R1Interval): R1Interval {
        lo = max(lo, y.lo)
        hi = min(hi, y.hi)
        return this
    }

    /**
     * Sets this interval as an empty interval.
     *
     * @return this.
     */
    fun setEmpty(): R1Interval {
        lo = 1.0
        hi = 0.0
        return this
    }

    override fun equals(other: Any?): Boolean {
        if (other is R1Interval) {
            // Return true if two intervals contain the same set of points.
            return lo == other.lo && hi == other.hi || isEmpty && other.isEmpty
        }
        return false
    }

    override fun hashCode(): Int {
        if (isEmpty) {
            return 17
        }
        var value: Long = 17
        value = 37 * value + java.lang.Double.doubleToLongBits(lo)
        value = 37 * value + java.lang.Double.doubleToLongBits(hi)
        return (value xor (value ushr 32)).toInt()
    }

    /**
     * Check is this interval and another one are equals with an approximation of maximum maxError.
     * The empty interval is considered to be positioned arbitrarily on the real line, thus any interval with
     * (length <= 2*max_error) matches the empty interval.
     *
     * @param y an interval.
     * @param maxError The maximum error authorized to consider two interval equals.
     * @return true if this interval can be transformed into the given interval by moving each endpoint by at most
     * "max_error".
     */
    @JvmOverloads
    fun approxEquals(y: R1Interval, maxError: Double = 1e-15): Boolean {
        return when {
            isEmpty -> y.length <= 2 * maxError
            y.isEmpty -> length <= 2 * maxError
            else -> abs(y.lo - lo) <= maxError && abs(y.hi - hi) <= maxError
        }
    }

    override fun toString(): String {
        return "[$lo, $hi]"
    }

    companion object {

        /**
         * Returns an empty interval. (Any interval where lo > hi is considered empty.)
         *
         * @return An empty interval instance with lo = 1 and hi = 0
         */
        @JvmStatic
        fun empty(): R1Interval {
            return R1Interval(1, 0)
        }

        /**
         * Convenience method to construct an interval containing a single point.
         *
         * @param p The unique point in the interval
         * @return An interval instance that contains only the given point.
         */
        @JvmStatic
        fun fromPoint(p: Double): R1Interval {
            return R1Interval(p, p)
        }

        /**
         * Convenience method to construct the minimal interval containing the two given points. This is equivalent to
         * starting with an empty interval and calling addPoint() twice, but it is more efficient.
         *
         * @param p1 A point to be included.
         * @param p2 A point to be included.
         * @return The minimum interval that contains the two points.
         */
        @JvmStatic
        fun fromPointPair(p1: Double, p2: Double): R1Interval = when {
            p1 <= p2 -> R1Interval(p1, p2)
            else -> R1Interval(p2, p1)
        }

    }

}
