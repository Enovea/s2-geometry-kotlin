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

import dilivia.math.M_PI
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.round
import java.util.*

/**
 * This class represents a one-dimensional angle (as opposed to a two-dimensional solid angle). It has methods for
 * converting angles to or from radians, degrees, and the e5/e6/e7 representations (i.e. degrees multiplied by
 * 1e5/1e6/1e7 and rounded to the nearest integer).
 *
 * The internal representation is a double-precision value in radians, so conversion to and from radians is exact.
 * Conversions between e5, e6, e7, and Degrees are not always exact; for example, degrees(3.1) is different from
 * e6(3100000) or e7(310000000).
 * However, the following properties are guaranteed for any integer "n", provided that "n" is in the input range of
 * both functions:
 *
 *     degrees(n) == e6(1000000 * n)
 *     degrees(n) == e7(10000000 * n)
 *          e6(n) == e7(10 * n)
 *
 * The corresponding properties are *not* true for e5, so if you use e5 then don't test for exact equality when
 * comparing to other formats such as Degrees or e7.
 *
 * The following conversions between degrees and radians are exact:
 *
 *          degrees(180) == radians(M_PI)
 *       degrees(45 * k) == radians(k * M_PI / 4)  for k == 0..8
 *
 * These identities also hold when the arguments are scaled up or down by any power of 2. Some similar identities are
 * also true, for example, degrees(60) == radians(M_PI / 3), but be aware that this type of identity does not hold in
 * general. For example, degrees(3) != radians(M_PI / 60).
 *
 * Similarly, the conversion to radians means that Angle::degrees(x).degrees() does not always equal "x".  For example,
 *
 *         S1Angle.degrees(45 * k).degrees() == 45 * k      for k == 0..8
 *   but       S1Angle.degrees(60).degrees() != 60.
 *
 * This means that when testing for equality, you should allow for numerical errors (EXPECT_DOUBLE_EQ) or convert to
 * discrete e5/e6/e7 values first.
 *
 * CAVEAT: All of the above properties depend on "double" being the usual 64-bit IEEE 754 type (which is true on almost
 * all modern platforms).
 *
 * This class is a port of the S1Angle class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S1Angle : Comparable<S1Angle>, Cloneable {

    /** Value of the angle in radian. */
    var radians: Double

    /**
     * This is internal to avoid ambiguity about which units are expected.
     */
    private constructor(radians: Double) {
        this.radians = radians
    }

    /**
     * The default constructor yields a zero angle.
     */
    constructor() : this(0.0)

    /**
     * Compute the angle between two points, which is also equal to the distance between these points on the unit
     * sphere. The points do not need to be normalized. This function has a maximum error of 3.25 * DBL_EPSILON (or
     * 2.5 * DBL_EPSILON for angles up to 1 radian). If either point is zero-length (e.g. an uninitialized S2Point), or
     * almost zero-length, the resulting angle will be zero.
     *
     * @param x point 1
     * @param y point 2
     */
    constructor(x: S2Point, y: S2Point) : this(x.angle(y))

    /**
     * Like the constructor with S2Point, but return the angle (i.e., distance) between two S2LatLng points. This
     * function has about 15 digits of accuracy for small distances but only about 8 digits of accuracy as the distance
     * approaches 180 degrees (i.e., nearly-antipodal points).
     *
     * @param x point 1
     * @param y point 2
     */
    constructor(x: S2LatLng, y: S2LatLng): this(x.getDistance(y).radians)

    public override fun clone(): S1Angle = radians(radians)

    /**
     * Gets the degrees value of the angle.
     *
     * @return degrees value of the angle.
     */
    fun degrees(): Double = radians * (180 / M_PI)

    /**
     * Gets the e5 degrees value of the angle.
     *
     * @return e5 degrees value of the angle.
     */
    fun e5(): Int = round(degrees() * 1e5).toInt()

    /**
     * Gets the e6 degrees value of the angle.
     *
     * @return e6 degrees value of the angle.
     */
    fun e6(): Int = round(degrees() * 1e6).toInt()

    /**
     * Gets the e7 degrees value of the angle.
     *
     * @return e7 degrees value of the angle.
     */
    fun e7(): Int = round(degrees() * 1e7).toInt()

    /**
     * Compute the absolute of this angle.
     *
     * @return the absolute value of an angle.
     */
    fun abs(): S1Angle = S1Angle(FastMath.abs(radians))

    override fun equals(other: Any?): Boolean = other is S1Angle && this.radians == other.radians

    override fun hashCode(): Int {
        val value = java.lang.Double.doubleToLongBits(radians)
        return (value xor (value ushr 32)).toInt()
    }

    /**
     * Defines the comparison of 2 angles.
     *
     * @param other The angle to compare.
     * @return -1 if this angle is less than other, 0 in case of equality and 1 otherwise.
     */
    override fun compareTo(other: S1Angle): Int = radians.compareTo(other.radians)

    /**
     * @return -this
     */
    operator fun unaryMinus(): S1Angle = S1Angle(-radians)

    /**
     * Defines the sum of 2 angles.
     *
     * @param other an angle
     * @return this + other
     */
    operator fun plus(other: S1Angle): S1Angle = S1Angle(radians + other.radians)

    operator fun plusAssign(other: S1Angle): Unit {
        radians += other.radians
    }

    /**
     * @param other an angle
     * @return this - other
     */
    operator fun minus(other: S1Angle): S1Angle = S1Angle(radians - other.radians)

    operator fun minusAssign(other: S1Angle): Unit {
        radians -= other.radians
    }

    /**
     * @param other an angle
     * @return this * other
     */
    operator fun times(other: S1Angle): S1Angle = S1Angle(radians * other.radians)

    operator fun timesAssign(other: S1Angle) {
        radians *= other.radians
    }

    /**
     * @param other an angle value in radians
     * @return this * other
     */
    operator fun times(other: Double): S1Angle = S1Angle(radians * other)

    operator fun timesAssign(other: Double) {
        radians *= other
    }

    /**
     * @param other an angle
     * @return this / other
     */
    operator fun div(other: S1Angle): S1Angle = S1Angle(radians / other.radians)

    operator fun divAssign(other: S1Angle) {
        radians /= other.radians
    }

    /**
     * @param other an angle value in radians
     * @return this / other
     */
    operator fun div(other: Double): S1Angle = S1Angle(radians / other)

    operator fun divAssign(other: Double) {
        radians /= other
    }

    /**
     * Compute the sine of this angle.
     *
     * @return The sine of this angle.
     */
    fun sin(): Double = FastMath.sin(radians)

    /**
     * Compute the cosine of this angle.
     *
     * @return The cosine of this angle.
     */
    fun cos(): Double = FastMath.cos(radians)

    /**
     * Compute the tangent of this angle.
     *
     * @return The tangent of this angle.
     */
    fun tan(): Double = FastMath.tan(radians)

    /**
     * Get a normalized instance of this angle.
     *
     * @return the angle normalized to the range (-180, 180] degrees.
     */
    fun normalized(): S1Angle {
        var rad = IEEEremainder(radians, 2.0 * M_PI)
        if (rad <= -M_PI) rad = M_PI
        return S1Angle(rad)
    }

    /**
     * Normalize this angle
     * @return this after normalization.
     */
    fun normalize(): S1Angle {
        var rad = IEEEremainder(radians, 2.0 * M_PI)
        if (rad <= -M_PI) rad = M_PI
        this.radians = rad
        return this
    }

    /**
     * Writes the angle in degrees with a "d" suffix, e.g. "17.3745d". By default 6 digits are printed; this can be
     * changed using setprecision(). Up to 17 digits are required to distinguish one angle from another.
     */
    override fun toString(): String {
        return "%.7fd".format(Locale.ENGLISH, degrees())
    }

    companion object {

        /**
         * Creates a S1Angle instance with this given measure in radians.
         *
         * @param radians The angle value in radians.
         * @return The S1Angle instance.
         */
        @JvmStatic
        fun radians(radians: Double): S1Angle = S1Angle(radians)
        fun radians(radians: Int): S1Angle = S1Angle(radians.toDouble())

        /**
         * Creates a S1Angle instance with this given measure in degrees.
         *
         * @param degrees The angle value in degrees.
         * @return The S1Angle instance.
         */
        @JvmStatic
        fun degrees(degrees: Double): S1Angle = S1Angle(degrees * (M_PI / 180))
        fun degrees(degrees: Int): S1Angle = degrees(degrees.toDouble())

        /**
         * Creates a S1Angle instance with this given measure in e5 degrees.
         *
         * @param e5 The angle value in degrees * 1e5.
         * @return The S1Angle instance.
         */
        @JvmStatic
        fun e5(e5: Int): S1Angle = degrees(e5 * 1e-5)

        /**
         * Creates a S1Angle instance with this given measure in e6 degrees.
         *
         * @param e6 The angle value in degrees * 1e6.
         * @return The S1Angle instance.
         */
        // Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6, but it's about 10 times faster and more than
        // accurate enough.
        @JvmStatic
        fun e6(e6: Int): S1Angle = degrees(e6 * 1e-6)

        /**
         * Creates a S1Angle instance with this given measure in e7 degrees.
         *
         * @param e7 The angle value in degrees * 1e7.
         * @return The S1Angle instance.
         */
        @JvmStatic
        fun e7(e7: Int): S1Angle = degrees(e7 * 1e-7)

        /**
         * Creates an infinity S1Angle instance.
         *
         * @return an angle larger than any finite angle.
         */
        fun infinity(): S1Angle = S1Angle(Double.MAX_VALUE)

        /**
         * An explicit shorthand for the default constructor.
         *
         * @return A 0.0 radians angle instance.
         */
        fun zero(): S1Angle = S1Angle(0.0)

        @JvmStatic
        fun max(left: S1Angle, right: S1Angle): S1Angle {
            return maxOf(left, right)
        }

        @JvmStatic
        fun min(left: S1Angle, right: S1Angle): S1Angle {
            return minOf(left, right)
        }

        @JvmStatic
        fun cos(angle: S1Angle) = angle.cos()

        @JvmStatic
        fun sin(angle: S1Angle)= angle.sin()

        @JvmStatic
        fun tan(angle: S1Angle) = angle.tan()

        operator fun Int.times(angle: S1Angle): S1Angle = angle * this.toDouble()

        operator fun UInt.times(angle: S1Angle): S1Angle = angle * this.toDouble()

        operator fun Double.times(angle: S1Angle): S1Angle = angle * this

        operator fun Double.div(angle: S1Angle): S1Angle = radians(this / angle.radians)

    }
}


