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
import dilivia.math.DoubleType.epsilon
import dilivia.math.M_PI
import dilivia.s2.S2PointUtil.isUnitLength
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.FastMath.asin
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * S1ChordAngle represents the angle subtended by a chord (i.e., the straight line segment connecting two points on the
 * sphere).  Its representation makes it very efficient for computing and comparing distances, but unlike S1Angle it is
 * only capable of representing angles between 0 and Pi radians. Generally, S1ChordAngle should only be used in loops
 * where many angles need to be calculated and compared. Otherwise it is simpler to use S1Angle.
 *
 * S1ChordAngle also loses some accuracy as the angle approaches Pi radians. Specifically, the representation of
 * (Pi - x) radians has an error of about (1e-15 / x), with a maximum error of about 2e-8 radians (about 13cm on the
 * Earth's surface).  For comparison, for angles up to 90 degrees (10000km)the worst-case representation error is about
 * 2e-16 radians (1 nanometer), which is about the same as S1Angle.
 *
 * S1ChordAngles are represented by the squared chord length, which can range from 0 to 4. infinity uses an infinite
 * squared length.
 *
 * This class is a port of the S1ChordAngle class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S1ChordAngle : Comparable<S1ChordAngle>, Cloneable {

    /** The squared length of the chord. */
    var length2: Double

    /**
     * Builds a S1ChordAngle instance.
     * @param length2 squared length of the chord.
     */
    constructor(length2: Double) {
        this.length2 = length2
        checkState({ isValid }, { "S1ChordAngle(length2 = $length2) is not valid." })
    }

    /** The default constructor yields a zero angle. */
    constructor() : this(0.0)

    /**
     * Construct the S1ChordAngle corresponding to the distance between the two
     * given points. The points must be unit length.
     *
     * @param x A unit length point.
     * @param y 1 unit length point.
     */
    constructor(x: S2Point, y: S2Point) : this(min(4.0, (x - y).norm2())) {
        checkState({ isUnitLength(x) }, { "Point x = $x is not unit length" })
        checkState({ isUnitLength(y) }, { "Point y = $y is not unit length" })
        // The squared distance may slightly exceed 4.0 due to roundoff errors.
        // The maximum error in the result is 2 * DBL_EPSILON * length2.
    }

    /**
     * Conversion from an S1Angle.
     *
     * Angles outside the range [0, Pi] are handled as follows:
     *   - infinity() is mapped to infinity(),
     *   - negative angles are mapped to negative(),
     *   - and finite angles larger than Pi are mapped to straight().
     *
     * Note that this operation is relatively expensive and should be avoided. To use S1ChordAngle effectively, you
     * should structure your code so that input arguments are converted to S1ChordAngles at the beginning of your
     * algorithm, and results are converted back to S1Angles only at the end.
     *
     * @param angle an angle.
     */
    constructor(angle: S1Angle) : this(
        when {
            angle.radians < 0.0 -> -1.0
            angle == S1Angle.infinity() -> Double.MAX_VALUE
            else -> {
                // The chord length is 2 * sin(angle / 2).
                val length = if (angle.radians.isNaN()) 2.0 else 2.0 * FastMath.sin(0.5 * min(M_PI, angle.radians))
                length * length
            }
        }
    )

    constructor(distance: S1ChordAngle) : this(distance.length2)

    public override fun clone(): S1ChordAngle = S1ChordAngle(this.length2)

    /**
     * Converts to an S1Angle.  Can be used just like an S1Angle constructor:
     *
     * infinity() is converted to S1Angle.infinity(), and negative() is converted to an unspecified negative S1Angle.
     *
     * Note that the conversion uses trigonometric functions and therefore should be avoided in inner loops.
     *
     * @return The angle represented by this chord.
     */
    fun toAngle(): S1Angle = when {
        isNegative() -> S1Angle.radians(-1)
        isInfinity() -> S1Angle.infinity()
        else -> S1Angle.radians(2.0 * asin(0.5 * sqrt(length2)))
    }

    // Convenience methods implemented by calling ToAngle() first.  Note that
    // because of the S1Angle conversion these methods are relatively expensive
    // (despite their lowercase names), so the results should be cached if they
    // are needed inside loops.
    fun radians(): Double = toAngle().radians
    fun degrees(): Double = toAngle().degrees()
    fun e5(): Int = toAngle().e5()
    fun e6(): Int = toAngle().e6()
    fun e7(): Int = toAngle().e7()

    /**
     * Defines the comparison of 2 chord.
     *
     * @param other a chord.
     * @return length2.compareTo(other.length2)
     */
    override fun compareTo(other: S1ChordAngle): Int = length2.compareTo(other.length2)

    /**
     * Indicates if this chord is zero length.
     *
     * @return true if length2 == 0.0
     */
    fun isZero(): Boolean = length2 == 0.0

    /**
     * Indicates if this chord is negative.
     *
     * @return true if length2 < 0.0
     */
    fun isNegative(): Boolean = length2 < 0.0

    /**
     * Indicates if this chord is infinite.
     *
     * @return true if length2 == Double.MAX_VALUE
     */
    fun isInfinity(): Boolean = length2 == Double.MAX_VALUE

    /**
     * Indicates if this chord is a special case (negative or infinity).
     *
     * @return true if isNegative() || isInfinity()
     */
    fun isSpecial(): Boolean = isNegative() || isInfinity()

    // Only addition and subtraction of S1ChordAngles is supported.  These
    // methods add or subtract the corresponding S1Angles, and clamp the result
    // to the range [0, Pi].  Both arguments must be non-negative and
    // non-infinite.
    @Throws(IllegalStateException::class, IllegalArgumentException::class)
    operator fun plus(b: S1ChordAngle): S1ChordAngle {
        // Note that this method is much more efficient than converting the chord
        // angles to S1Angles and adding those.  It requires only one square root
        // plus a few additions and multiplications.
        checkState({ !isSpecial() }, { "This angle $this is special" })
        requireArgument({ !b.isSpecial() }, { "The angle $b to add is special." })

        // Optimization for the common case where "b" is an error tolerance
        // parameter that happens to be set to zero.
        val a2 = length2
        val b2 = b.length2
        return when {
            b2 == 0.0 -> this

            // Clamp the angle sum to at most 180 degrees.
            a2 + b2 >= kMaxLength2 -> straight()

            // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
            // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
            // Then the formula below can be derived from c = 2 * sin(A+B) and the
            // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
            //                 cos(X) = sqrt(1 - sin^2(X)) .
            else -> {
                val x = a2 * (1 - 0.25 * b2)  // isValid) => non-negative
                val y = b2 * (1 - 0.25 * a2)  // isValid) => non-negative
                S1ChordAngle(min(kMaxLength2, x + y + 2 * sqrt(x * y)))
            }
        }
    }

    operator fun minus(b: S1ChordAngle): S1ChordAngle {
        // See comments in operator+().
        checkState({!isSpecial()}, { "This angle $this is special." })
        requireArgument({!b.isSpecial()}, { "Angle $b to subtract is special." })
        val a2 = length2
        val b2 = b.length2
        return when {
            b2 == 0.0 -> this
            a2 <= b2 -> zero()
            else -> {
                val x = a2 * (1 - 0.25 * b2)
                val y = b2 * (1 - 0.25 * a2)
                S1ChordAngle(max(0.0, x + y - 2 * sqrt(x * y)))
            }
        }
    }


    // Trigonmetric functions.

    /**
     * Compute the sine of the angle represented by this chord. (= sin(toAngle()))
     * It is more accurate and efficient to call these rather than first converting to an S1Angle.
     * @return the sine of this angle.
     */
    fun sin(): Double = sqrt(sin2())

    /**
     * Compute the cosine of the angle represented by this chord. (= cos(toAngle()))
     * It is more accurate and efficient to call these rather than first converting to an S1Angle.
     * @return the cosine of this angle.
     */
    fun cos(): Double {
        // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
        checkState({!isSpecial()}, { "This angle $this is special." })
        return 1 - 0.5 * length2
    }

    /**
     * Compute the tangente of the angle represented by this chord. (= tan(toAngle()))
     * It is more accurate and efficient to call these rather than first converting to an S1Angle.
     * @return the tangente of this angle.
     */
    fun tan(): Double {
        val tan = sin() / cos()
        return if (tan == -0.0) 0.0 else tan
    }

    /**
     * Compute the squared sine of the angle.
     *
     * @return sin(a)**2, but computed more efficiently.
     */
    @Throws(IllegalStateException::class)
    fun sin2(): Double {
        checkState({!isSpecial()}, { "This angle $this is special." })
        // Let "a" be the (non-squared) chord length, and let A be the corresponding
        // half-angle (a = 2*sin(A)).  The formula below can be derived from:
        //   sin(2*A) = 2 * sin(A) * cos(A)
        //   cos^2(A) = 1 - sin^2(A)
        // This is much faster than converting to an angle and computing its sine.
        return length2 * (1 - 0.25 * length2)
    }

    /**
     * Get the smallest representable S1ChordAngle larger than this object. This can be used to convert a "<" comparison
     * to a "<=" comparison. For example:
     *
     * val query = S2ClosestEdgeQuery(...)
     * val limit: S1ChordAngle = ...
     * if (query.isDistanceLess(target, limit.Successor())) {
     *   // Distance to "target" is less than or equal to "limit".
     * }
     *
     * Note the following special cases:
     *   - negative().successor() == zero()
     *   - straight().successor() == infinity()
     *   - infinity().successor() == infinity()
     *
     * @return the successor of this angle.
     */
    fun successor(): S1ChordAngle = when {
        length2 >= kMaxLength2 -> infinity()
        length2 < 0.0 -> zero()
        else -> S1ChordAngle(FastMath.nextUp(length2))
    }

    /**
     * Get the largest representable S1ChordAngle less than this object. This can be used to convert a ">" comparison
     * to a ">=" comparison.
     *
     * Note the following special cases:
     *   - negative().predecessor() == Negative()
     *   - straight().predecessor() == Negative()
     *   - infinity().predecessor() == Straight()
     *
     * @return the successor of this angle.
     */
    fun predecessor(): S1ChordAngle = when {
        length2 <= 0.0 -> negative()
        length2 > kMaxLength2 -> straight()
        else -> S1ChordAngle(FastMath.nextDown(length2))
    }

    // Returns a new S1ChordAngle that has been adjusted by the given error
    // bound (which can be positive or negative).  "error" should be the value
    // returned by one of the error bound methods below.  For example:
    //    S1ChordAngle a(x, y);
    //    S1ChordAngle a1 = a.plusError(a.getS2PointConstructorMaxError());
    fun plusError(error: Double): S1ChordAngle {
        // If angle is Negative() or Infinity(), don't change it. Otherwise clamp it to the valid range.
        return if (isSpecial()) this
        else S1ChordAngle(max(0.0, min(kMaxLength2, length2 + error)))
    }

    fun plusError(error: Int): S1ChordAngle = plusError(error.toDouble())

    // Return the maximum error in length2() for the S1ChordAngle(x, y)
    // constructor, assuming that "x" and "y" are normalized to within the
    // bounds guaranteed by S2Point::Normalize().  (The error is defined with
    // respect to the true distance after the points are projected to lie
    // exactly on the sphere.)
    fun getS2PointConstructorMaxError(): Double {
        // There is a relative error of 2.5 * DBL_EPSILON when computing the squared
        // distance, plus a relative error of 2 * DBL_EPSILON and an absolute error
        // of (16 * DBL_EPSILON**2) because the lengths of the input points may
        // differ from 1 by up to (2 * DBL_EPSILON) each.  (This is the maximum
        // length error in S2Point::Normalize.)
        return 4.5 * epsilon * length2 + 16 * epsilon * epsilon
    }

    // Return the maximum error in length2() for the S1Angle constructor.
    fun getS1AngleConstructorMaxError(): Double {
        // Assuming that an accurate math library is being used, the sin() call and
        // the multiply each have a relative error of 0.5 * DBL_EPSILON.
        return epsilon * length2
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S1ChordAngle) return false

        if (length2 != other.length2) return false

        return true
    }

    override fun hashCode(): Int {
        return length2.hashCode()
    }

    override fun toString(): String {
        return "S1ChordAngle(length2=$length2, angle=${toAngle()})"
    }

    /**
     * Indicates if this angle is valid. true if the internal representation is valid. negative() and infinity() are
     * both considered valid.
     */
    val isValid: Boolean
        get() = (length2 in 0.0..kMaxLength2) || isSpecial()

    companion object {

        const val kMaxLength2 = 4.0

        /** the zero chord angle. */
        fun zero() = S1ChordAngle()

        /** the chord angle of 90 degrees (a "right angle"). */
        fun right() = S1ChordAngle(2.0)

        /** The chord angle of 180 degrees (a "straight angle"). This is the maximum finite chord angle. */
        fun straight() = S1ChordAngle(4.0)

        // Return a chord angle larger than any finite chord angle.  The only valid
        // operations on Infinity() are comparisons, S1Angle conversions, and
        // Successor() / Predecessor().
        fun infinity() = S1ChordAngle(Double.MAX_VALUE)

        // Return a chord angle smaller than Zero().  The only valid operations on
        // Negative() are comparisons, S1Angle conversions, and Successor() /
        // Predecessor().
        fun negative() = S1ChordAngle(-1.0)


        // Convenience methods implemented by converting from an S1Angle.
        fun radians(radians: Double): S1ChordAngle = S1ChordAngle(S1Angle.radians(radians))
        fun degrees(degrees: Double): S1ChordAngle = S1ChordAngle(S1Angle.degrees(degrees))
        fun degrees(degrees: Int): S1ChordAngle = degrees(degrees.toDouble())
        fun e5(e5: Int): S1ChordAngle = S1ChordAngle(S1Angle.e5(e5))
        fun e6(e6: Int): S1ChordAngle = S1ChordAngle(S1Angle.e6(e6))
        fun e7(e7: Int): S1ChordAngle = S1ChordAngle(S1Angle.e7(e7))

        // Construct an S1ChordAngle that is an upper bound on the given S1Angle,
        // i.e. such that FastUpperBoundFrom(x).ToAngle() >= x.  Unlike the S1Angle
        // constructor above, this method is very fast, and the bound is accurate to
        // within 1% for distances up to about 3100km on the Earth's surface.
        // This method uses the distance along the surface of the sphere as an upper
        // bound on the distance through the sphere's interior.
        fun fastUpperBoundFrom(angle: S1Angle): S1ChordAngle = fromLength2(angle.radians * angle.radians)

        // Construct an S1ChordAngle from the squared chord length.  Note that the
        // argument is automatically clamped to a maximum of 4.0 to handle possible
        // roundoff errors.  The argument must be non-negative.
        fun fromLength2(length2: Double): S1ChordAngle = S1ChordAngle(min(4.0, length2))
        fun fromLength2(length2: Int) = fromLength2(length2.toDouble())

        fun between(p1: S2Point, p2: S2Point): S1ChordAngle = S1ChordAngle(p1, p2)

        fun sin(angle: S1ChordAngle): Double = angle.sin()
        fun cos(angle: S1ChordAngle): Double = angle.cos()
        fun tan(angle: S1ChordAngle): Double = angle.tan()
    }
}

