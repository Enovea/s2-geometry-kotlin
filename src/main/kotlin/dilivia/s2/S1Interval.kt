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
package dilivia.s2

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.math.vectors.R2VectorDouble
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max

/**
 * An S1Interval represents a closed interval on a unit circle (also known as a 1-dimensional sphere).  It is capable
 * of representing the empty interval (containing no points), the full interval (containing all points), and zero-length
 * intervals (containing a single point).
 *
 * Points are represented by the angle they make with the positive x-axis in the range [-Pi, Pi].  An interval is
 * represented by its lower and upper bounds (both inclusive, since the interval is closed).  The lower bound may
 * be greater than the upper bound, in which case the interval is "inverted" (i.e. it passes through the point (-1, 0)).
 *
 * Note that the point (-1, 0) has two valid representations, Pi and -Pi. The normalized representation of this point
 * internally is Pi, so that endpoints of normal intervals are in the range (-Pi, Pi].  However, we take advantage of
 * the point -Pi to construct two special intervals:
 *   - the full() interval is [-Pi, Pi],
 *   - and the empty() interval is [Pi, -Pi].
 *
 * This class is a port of the S1Interval class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S1Interval : Cloneable {

    /** The lower bound of the interval. */
    var lo: Double

    /** The upper bound of the interval. */
    var hi: Double

    /**
     * Internal constructor that assumes that both arguments are in the correct range, i.e. normalization from
     * -Pi to Pi is already done.
     *
     * @param lo The lower bound of the interval. Must be in the range [-Pi, Pi]
     * @param hi The upper bound of the interval. Must be in the range [-Pi, Pi]
     */
    private constructor(lo: Double, hi: Double, checked: Boolean) {
        var newLo = lo
        var newHi = hi
        if (!checked) {
            if (lo == -M_PI && hi != M_PI) {
                newLo = M_PI
            }
            if (hi == -M_PI && lo != M_PI) {
                newHi = M_PI
            }
        }
        this.lo = newLo
        this.hi = newHi
        checkState({ isValid }, { "Interval $this is not valid." })
    }

    /**
     * Builds a new S1Interval instance. Both endpoints must be in the range -Pi to Pi inclusive. The value -Pi is
     * converted internally to Pi except for the full() and empty() intervals.
     *
     * @param lo The lo bound of the interval. Must be  in the range [-Pi, Pi]
     * @param hi The hi bound of the interval. Must be  in the range [-Pi, Pi]
     */
    constructor(lo: Double, hi: Double) : this(lo, hi, false)

    constructor(lo: Int, hi: Int) : this(lo.toDouble(), hi.toDouble(), false)

    /**
     * Copy constructor. Assumes that the given interval is valid.
     *
     * @param interval The interval to copy.
     */
    constructor(interval: S1Interval) : this(interval.lo, interval.hi, true)

    /** The default constructor creates an empty interval. */
    constructor() : this(M_PI, -M_PI, true)

    public override fun clone(): S1Interval = S1Interval(this)

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
        when (i) {
            0 -> lo = v
            1 -> hi = v
            else -> throw ArrayIndexOutOfBoundsException("Index $i is out of bounds 0..1")
        }
    }

    fun set(other: S1Interval) {
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
     * Indicate if the interval is valid. An interval is valid if neither bound exceeds Pi in absolute value, and the
     * value -Pi appears only in the empty() and full() intervals.
     */
    val isValid: Boolean
        get() = abs(lo) <= M_PI
                && abs(hi) <= M_PI
                && !(lo == -M_PI && hi != M_PI)
                && !(hi == -M_PI && lo != M_PI)

    /** Indicate if the interval contains all points on the unit circle.  */
    val isFull: Boolean
        get() = lo == -M_PI && hi == M_PI

    /** Indicate if the interval is empty, i.e. it contains no points.  */
    val isEmpty: Boolean
        get() = lo == M_PI && hi == -M_PI

    /** Indicate if the interval is not empty.  */
    val isNotEmpty: Boolean
        get() = !isEmpty

    /* Indicate if the interval is inverted. = true if lo > hi. (This is true for empty intervals.) */
    val isInverted: Boolean
        get() = lo > hi

    /**
     * The midpoint of the interval. For full and empty intervals, the result is arbitrary.
     */
    val center: Double
        get() {
            val center = 0.5 * (lo + hi)
            return if (!isInverted) {
                center
            } else {
                // Return the center in the range (-Pi, Pi].
                if (center <= 0) center + M_PI else center - M_PI
                // Empty intervals have a negative length.
            }
        }

    /**
     * Length of the interval. The length of an empty interval is negative.
     */
    val length: Double
        get() {
            var length = hi - lo
            if (length >= 0) {
                return length
            }
            length += 2 * M_PI
            // Empty intervals have a negative length.
            return if (length > 0) length else -1.0
        }

    /**
     * the complement of the interior of the interval. An interval and its complement have the same boundary
     * but do not share any interior values. The complement operator is not a bijection, since the complement of a
     * singleton interval (containing a single value) is the same as the complement of an empty interval.
     */
    val complement: S1Interval
        get() = if (lo == hi) {
            full()
        } else S1Interval(hi, lo, true)

    /**
     * the midpoint of the complement of the interval. For full and empty intervals, the result is arbitrary. For a
     * singleton interval (containing a single point), the result is its antipodal point on S1.
     */
    val complementCenter: Double
        get() = when {
            lo != hi -> complement.center
            hi <= 0.0 -> hi + M_PI
            else -> hi - M_PI
        }

    /**
     * Check if the interval contains a given point.
     *
     * @param p A point in range [-PI, PI].
     * @return true if the interval (which is closed) contains the point 'p'.
     */
    operator fun contains(p: Double): Boolean {
        // Works for empty, full, and singleton intervals.
        requireArgument { abs(p) <= M_PI }
        return fastContains(if (p == -M_PI) M_PI else p)
    }

    /**
     * @see S1Interval.contains
     */
    operator fun contains(p: Int) = contains(p.toDouble())

    /**
     * Check if the interior of the interval contains a point.
     *
     * @param p A point in range [-PI, PI].
     * @return true if the interior of the interval contains the point 'p'.
     */
    fun interiorContains(p: Double): Boolean {
        // Works for empty, full, and singleton intervals.
        requireArgument { abs(p) <= M_PI }
        var cp = p
        if (cp == -M_PI) {
            cp = M_PI
        }
        return if (isInverted) {
            cp > lo || cp < hi
        } else {
            (cp > lo && cp < hi) || isFull
        }
    }

    /**
     * @see S1Interval.interiorContains
     */
    fun interiorContains(p: Int) = interiorContains(p.toDouble())

    /**
     * Check if this interval contains all points of another one. Works for empty, full, and singleton intervals.
     *
     * @param y An interval.
     * @return true if the interval contains the given interval 'y'.
     */
    operator fun contains(y: S1Interval): Boolean {
        // It might be helpful to compare the structure of these tests to
        // the simpler Contains(double) method above.
        return if (isInverted) {
            if (y.isInverted) {
                y.lo >= lo && y.hi <= hi
            } else (y.lo >= lo || y.hi <= hi) && !isEmpty
        } else {
            if (y.isInverted) {
                isFull || y.isEmpty
            } else y.lo >= lo && y.hi <= hi
        }
    }

    /**
     * Check if the interior of this interval contains another one.
     * Note that x.interiorContains(x) is true only when x is the empty or full interval, and
     * x.interiorContains(S1Interval(p,p)) is equivalent to x.interiorContains(p).
     *
     * @param y An interval.
     * @return true if the interior of this interval contains the entire interval 'y'.
     */
    fun interiorContains(y: S1Interval): Boolean {
        return if (isInverted) {
            if (!y.isInverted) {
                y.lo > lo || y.hi < hi
            } else y.lo > lo && y.hi < hi || y.isEmpty
        } else {
            if (y.isInverted) {
                isFull || y.isEmpty
            } else y.lo > lo && y.hi < hi || isFull
        }
    }

    /**
     * Check if this interval intersects another one.
     * Note that the point +/-Pi has two representations, so the intervals [-Pi,-3] and [2,Pi] intersect, for example.
     *
     * @param y An interval.
     * @return true if the two intervals contain any points in common.
     */
    fun intersects(y: S1Interval): Boolean {
        if (isEmpty || y.isEmpty) {
            return false
        }
        return if (isInverted) {
            // Every non-empty inverted interval contains Pi.
            y.isInverted || y.lo <= hi || y.hi >= lo
        } else {
            if (y.isInverted) {
                y.lo <= hi || y.hi >= lo
            } else y.lo <= hi && y.hi >= lo
        }
    }

    /**
     * Check if the interior of this interval intersects another one.
     * Works for empty, full, and singleton intervals.
     *
     * @param y An interval.
     * @return true if the interior of this interval contains any point of the interval 'y' (including its boundary).
     */
    fun interiorIntersects(y: S1Interval): Boolean {
        if (isEmpty || y.isEmpty || lo == hi) {
            return false
        }
        return if (isInverted) {
            y.isInverted || y.lo < hi || y.hi > lo
        } else {
            if (y.isInverted) {
                y.lo < hi || y.hi > lo
            } else y.lo < hi && y.hi > lo || isFull
        }
    }

    /**
     * Expand the interval by the minimum amount necessary so that it contains the given point "p"
     * (an angle in the range [-Pi, Pi]).
     *
     * @param p The point to add.
     * @return this after modification.
     */
    fun addPoint(p: Double): S1Interval {
        requireArgument { abs(p) <= M_PI }
        var cp = p
        if (cp == -M_PI) {
            cp = M_PI
        }
        if (!fastContains(cp)) {
            if (isEmpty) {
                lo = cp
                hi = cp
            } else {
                // Compute distance from p to each endpoint.
                val dlo = positiveDistance(cp, lo)
                val dhi = positiveDistance(hi, cp)
                if (dlo < dhi) {
                    lo = cp
                } else {
                    hi = cp
                }
                // Adding a point can never turn a non-full interval into a full one.
            }
        }
        return this
    }

    fun addPoint(p: Int): S1Interval = addPoint(p.toDouble())

    /**
     * Project a point to this interval. The interval must be non-empty.
     *
     * @param p The point to project.
     * @return the closest point in the interval to the given point "p".
     */
    fun project(p: Double): Double {
        checkState({isNotEmpty}, { "Interval is empty." })
        requireArgument({abs(p) <= M_PI}, { "" })
        val cp = if (p == -M_PI) M_PI else p
        if (fastContains(cp)) return cp
        // Compute distance from p to each endpoint.
        // Compute distance from p to each endpoint.
        val dlo: Double = positiveDistance(cp, lo)
        val dhi: Double = positiveDistance(hi, cp)
        return if (dlo < dhi) lo else hi
    }

    fun project(p: Int): Double = project(p.toDouble())

    /**
     * Expand this interval in both direction.
     * Note that the expansion of an empty interval is always empty.
     *
     * @param margin The expansion amount.
     * @return an interval that contains all points within a distance "radius" of a point in this interval.
     */
    fun expanded(margin: Double): S1Interval {
        if (margin >= 0) {
            if (isEmpty) return this
            // Check whether this interval will be full after expansion, allowing
            // for a 1-bit rounding error when computing each endpoint.
            if (length + 2 * margin + 2 * DoubleType.epsilon >= 2 * M_PI) return full()
        } else {
            if (isFull) return this
            // Check whether this interval will be empty after expansion, allowing
            // for a 1-bit rounding error when computing each endpoint.
            if (length + 2 * margin - 2 * DoubleType.epsilon <= 0) return empty()
        }
        var result = S1Interval(IEEEremainder(lo - margin, 2 * M_PI), IEEEremainder(hi + margin, 2 * M_PI))
        if (result.lo <= -M_PI) result = S1Interval(M_PI, result.hi)
        return result
    }

    fun expanded(radius: Int) = expanded(radius.toDouble())

    /**
     * Compute the union of this interval and the interval "y".
     *
     * @param y An interval.
     * @return the smallest interval that contains this interval and the given interval "y".
     */
    fun union(y: S1Interval): S1Interval {
        // The y.is_full() case is handled correctly in all cases by the code
        // below, but can follow three separate code paths depending on whether
        // this interval is inverted, is non-inverted but contains Pi, or neither.

        if (y.isEmpty) return this
        if (fastContains(y.lo)) {
            if (fastContains(y.hi)) {
                // Either this interval contains y, or the union of the two
                // intervals is the Full() interval.
                if (contains(y)) return this  // is_full() code path
                return full()
            }
            return S1Interval(lo, y.hi, true)
        }
        if (fastContains(y.hi)) return S1Interval(y.lo, hi, true)

        // This interval contains neither endpoint of y.  This means that either y
        // contains all of this interval, or the two intervals are disjoint.
        if (isEmpty || y.fastContains(lo)) return y

        // Check which pair of endpoints are closer together.
        val dlo = positiveDistance(y.hi, lo)
        val dhi = positiveDistance(hi, y.lo)
        if (dlo < dhi) {
            return S1Interval(y.lo, hi, true)
        } else {
            return S1Interval(lo, y.hi, true)
        }
    }

    /**
     * Compute the intersection of this interval and another one.
     * Note that the region of intersection may consist of two disjoint intervals.
     *
     * @param y An interval.
     * @return the smallest interval that contains the intersection of this interval with "y".
     */
    fun intersection(y: S1Interval): S1Interval {
        // The y.isFull case is handled correctly in all cases by the code below, but can follow three separate code
        // paths depending on whether this interval is inverted, is non-inverted but contains Pi, or neither.
        if (y.isEmpty) {
            return empty()
        }
        if (fastContains(y.lo)) {
            return if (fastContains(y.hi)) {
                // Either this interval contains y, or the region of intersection
                // consists of two disjoint subintervals. In either case, we want
                // to return the shorter of the two original intervals.
                if (y.length < length) {
                    y // is_full() code path
                } else this
            } else S1Interval(y.lo, hi, true)
        }
        if (fastContains(y.hi)) {
            return S1Interval(lo, y.hi, true)
        }

        // This interval contains neither endpoint of y. This means that either y
        // contains all of this interval, or the two intervals are disjoint.
        return if (y.fastContains(lo)) {
            this // is_empty() okay here
        } else {
            checkState({ !intersects(y) }, { "Wrong state." })
            empty()
        }
    }

    fun setEmpty(): S1Interval {
        lo = M_PI
        hi = -M_PI
        return this
    }

    /**
     * @return true if two intervals contains the same set of points.
     */
    override fun equals(other: Any?): Boolean {
        if (other is S1Interval) {
            return lo == other.lo && hi == other.hi
        }
        return false
    }

    override fun hashCode(): Int {
        var value: Long = 17
        value = 37 * value + java.lang.Double.doubleToLongBits(lo)
        value = 37 * value + java.lang.Double.doubleToLongBits(hi)
        return (value ushr 32 xor value).toInt()
    }

    /**
     * Check if this interval is equal to another one with a max error tolerance.
     * Empty and full intervals are considered to start at an arbitrary point on the unit circle, thus any interval
     * with (length <= 2*max_error) matches the empty interval, and any interval with (length >= 2*Pi - 2*max_error)
     * matches the full interval.
     *
     * @return true if this interval can be transformed into the given interval by moving each endpoint by at most
     * "maxError" (and without the endpoints crossing, which would invert the interval).
     */
    @JvmOverloads
    fun approxEquals(y: S1Interval, maxError: Double = 1e-15): Boolean {
        return when {
            isEmpty -> y.length <= 2 * maxError
            y.isEmpty -> length <= 2 * maxError
            isFull -> y.length >= 2 * (M_PI - maxError)
            y.isFull -> length >= 2 * (M_PI - maxError)
            else -> abs(IEEEremainder(y.lo - lo, 2 * M_PI)) <= maxError
                    && abs(IEEEremainder(y.hi - hi, 2 * M_PI)) <= maxError
                    && abs(length - y.length) <= 2 * maxError
        }
    }

    /**
     * Return true if the interval (which is closed) contains the point 'p'. Skips
     * the normalization of 'p' from -Pi to Pi.
     *
     */
    protected fun fastContains(p: Double): Boolean {
        return if (isInverted) {
            (p >= lo || p <= hi) && !isEmpty
        } else {
            p in lo..hi
        }
    }

    override fun toString(): String {
        return "[$lo, $hi]"
    }

    /**
     * Compute the Hausdorff distance with the interval 'y'.
     * For two S1Intervals x and y, this distance is defined by
     *     h(x, y) = max_{p in x} min_{q in y} d(p, q),
     * where d(.,.) is measured along S1.
     *
     * @param y An interval.
     * @return the Hausdorff distance to the given interval 'y'.
     */
    fun getDirectedHausdorffDistance(y: S1Interval): Double {
        if (y.contains(this)) return 0.0  // this includes the case this is empty
        if (y.isEmpty) return M_PI        // maximum possible distance on S1

        val yComplementCenter = y.complementCenter
        if (contains(yComplementCenter)) {
            return positiveDistance(y.hi, yComplementCenter)
        } else {
            // The Hausdorff distance is realized by either two hi) endpoints or two lo) endpoints, whichever is farther
            // apart.
            val hiHi = if (S1Interval(y.hi, yComplementCenter).contains(hi)) positiveDistance(y.hi, hi) else 0.0
            val loLo = if (S1Interval(yComplementCenter, y.lo).contains(lo)) positiveDistance(lo, y.lo) else 0.0
            checkState({hiHi > 0 || loLo > 0}, { "" })
            return max(hiHi, loLo)
        }
    }

    companion object {

        /** The empty interval.*/
        fun empty(): S1Interval = S1Interval(M_PI, -M_PI, true)

        /** the full interval. */
        fun full(): S1Interval = S1Interval(-M_PI, M_PI, true)

        /**
         * Convenience method to construct an interval containing a single point.
         *
         * @param p A point in [-PI, PI]
         * @return An interval that contains only the given point.
         */
        fun fromPoint(p: Double): S1Interval {
            var checkedPoint = p
            if (checkedPoint == -M_PI) {
                checkedPoint = M_PI
            }
            return S1Interval(checkedPoint, checkedPoint, true)
        }

        /**
         * Convenience method to construct the minimal interval containing the two given points. This is equivalent to
         * starting with an empty interval and calling addPoint() twice, but it is more efficient.
         *
         * @param p1 A point in [-PI, PI].
         * @param p2 A point in [-PI, PI].
         * @return The minimum interval that contains the 2 points.
         */
        @JvmStatic
        fun fromPointPair(p1: Double, p2: Double): S1Interval {
            requireArgument { abs(p1) <= M_PI }
            requireArgument { abs(p2) <= M_PI }
            val cp1 = if (p1 == -M_PI) M_PI else p1
            val cp2 = if (p2 == -M_PI) M_PI else p2
            return if (positiveDistance(cp1, cp2) <= M_PI) {
                S1Interval(cp1, cp2, true)
            } else {
                S1Interval(cp2, cp1, true)
            }
        }

        /**
         * Compute the distance from "a" to "b" in the range [0, 2*Pi). This is
         * equivalent to (remainder(b - a - S2.M_PI, 2 * S2.M_PI) + S2.M_PI), except that
         * it is more numerically stable (it does not lose precision for very small
         * positive distances).
         */
        fun positiveDistance(a: Double, b: Double): Double {
            // Compute the distance from "a" to "b" in the range [0, 2*Pi).
            // This is equivalent to (remainder(b - a - M_PI, 2 * M_PI) + M_PI),
            // except that it is more numerically stable (it does not lose
            // precision for very small positive distances).
            // Compute the distance from "a" to "b" in the range [0, 2*Pi).
            // This is equivalent to (remainder(b - a - M_PI, 2 * M_PI) + M_PI),
            // except that it is more numerically stable (it does not lose
            // precision for very small positive distances).
            val d = b - a
            return if (d >= 0) d
            // We want to ensure that if b == Pi and a == (-Pi + eps),
            // the return result is approximately 2*Pi and not zero.
            // We want to ensure that if b == Pi and a == (-Pi + eps),
            // the return result is approximately 2*Pi and not zero.
            else b + M_PI - (a - M_PI)
        }
    }
}
