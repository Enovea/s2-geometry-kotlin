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

import dilivia.PreConditions.requireArgument
import dilivia.math.DoubleType
import dilivia.math.matrix.Matrix3x3Double
import dilivia.math.vectors.R3VectorDouble
import dilivia.math.vectors.times
import dilivia.s2.S1Angle.Companion.cos
import dilivia.s2.S1Angle.Companion.sin
import org.apache.commons.math3.util.FastMath.abs

typealias S2Point = R3VectorDouble

object S2PointUtil {

    var testDegeneracies: Boolean = false

    /**
     * Get a unique "origin" on the sphere for operations that need a fixed reference point.
     * In particular, this is the "point at infinity" used for point-in-polygon testing (by counting the number of
     * edge crossings).
     *
     * @return The S2Point origin.
     */
    fun origin(): S2Point {
        if (testDegeneracies) {
            // This value makes polygon operations much slower, because it greatly
            // increases the number of degenerate cases that need to be handled using
            // s2pred::ExpensiveSign().
            return dilivia.s2.S2Point(0, 0, 1)
        } else {
            // The origin should not be a point that is commonly used in edge tests in
            // order to avoid triggering code to handle degenerate cases.  (This rules
            // out the north and south poles.)  It should also not be on the boundary of
            // any low-level S2Cell for the same reason.
            //
            // The point chosen here is about 66km from the north pole towards the East
            // Siberian Sea.  See the unittest for more details.  It is written out
            // explicitly using floating-point literals because the optimizer doesn't
            // seem willing to evaluate Normalize() at compile time.
            return dilivia.s2.S2Point(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195)
        }
    }

    /**
     * Check if a S2Point is unit length.
     *
     * @return true if the given point is approximately unit length (this is mainly useful for assertions).
     */
    fun isUnitLength(point: S2Point): Boolean = abs(point.norm2() - 1) <= 5 * DoubleType.epsilon

    /**
     * Check if two points are approximately equals.
     * It is an error if either point is a zero-length vector (default S2Point), but this is only checked in debug mode.
     * In non-debug mode it will always return true.
     *
     * @return true if two points are within the given distance of each other (this is mainly useful for testing).
     */
    fun approxEquals(a: S2Point, b: S2Point, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        requireArgument { a != dilivia.s2.S2Point() }
        requireArgument { b != dilivia.s2.S2Point() }
        return S1Angle(a, b) <= maxError
    }

    fun ortho(a: S2Point): S2Point = if (testDegeneracies) {
        // Vector3.ortho() always returns a point on the X-Y, Y-Z, or X-Z planes.
        // This leads to many more degenerate cases in polygon operations.
        a.ortho()
    } else {
        var k = a.largestAbsComponent() - 1
        if (k < 0) k = 2
        val temp = dilivia.s2.S2Point(0.012, 0.0053, 0.00457)
        temp[k] = 1.0
        a.crossProd(temp).normalize()
    }

    /**
     * Get a vector "c" that is orthogonal to the given unit-length vectors "a" and "b".
     *
     * This function is similar to a.crossProd(b) except that it does a better job of ensuring orthogonality when "a"
     * is nearly parallel to "b", and it returns a non-zero result even when a == b or a == -b.
     *
     * It satisfies the following properties (RCP == robustCrossProd):
     *
     * (1) RCP(a,b) != 0 for all a, b
     * (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b
     * (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
     * (4) RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
     *
     * The result is not guaranteed to be unit length.
     *
     * @param a A S2 point.
     * @param b A S2 point.
     * @return An orthogonal vector of a and b.
     */
    fun robustCrossProd(a: S2Point, b: S2Point): S2Point {
        // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
        // approaches zero.  This leads to situations where a.CrossProd(b) is not
        // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
        // but we also want b.RobustCrossProd(a) == -a.RobustCrossProd(b).
        //
        // The easiest fix is to just compute the cross product of (b+a) and (b-a).
        // Mathematically, this cross product is exactly twice the cross product of
        // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
        // are always perpendicular (since "a" and "b" are unit length).  This
        // yields a result that is nearly orthogonal to both "a" and "b" even if
        // these two values differ only in the lowest bit of one component.

        requireArgument({ isUnitLength(a) }, { "Point a $a is not unit: ${a.norm()}" })
        requireArgument({ isUnitLength(b) }, { "Point a $b is not unit: ${b.norm()}" })
        val x = (b + a).crossProd(b - a)
        if (x != dilivia.s2.S2Point(0, 0, 0)) return x

        // The only result that makes sense mathematically is to return zero, but
        // we find it more convenient to return an arbitrary orthogonal vector.
        return ortho(a)
    }

    /**
     * Rotate the given point about the given axis by the given angle. "p" and "axis" must be unit length; "angle" has
     * no restrictions (e.g., it can be positive, negative, greater than 360 degrees, etc).
     *
     * @param p The point to rotate
     * @param axis The rotation axis.
     * @param angle The rotation angle.
     * @return The rotation result.
     */
    fun rotate(p: S2Point, axis: S2Point, angle: S1Angle): S2Point {
        requireArgument { isUnitLength(p) }
        requireArgument { isUnitLength(axis) }
        // Let M be the plane through P that is perpendicular to "axis", and let
        // "center" be the point where M intersects "axis".  We construct a
        // right-handed orthogonal frame (dx, dy, center) such that "dx" is the
        // vector from "center" to P, and "dy" has the same length as "dx".  The
        // result can then be expressed as (cos(angle)*dx + sin(angle)*dy + center).
        val center = p.dotProd(axis) * axis
        val dx = p - center
        val dy = axis.crossProd(p)
        // Mathematically the result is unit length, but normalization is necessary
        // to ensure that numerical errors don't accumulate.
        return (cos(angle) * dx + sin(angle) * dy + center).normalize()
    }

    /**
     * Extend the given point "z" on the unit sphere into a right-handed coordinate frame of unit-length column
     * vectors m = (x,y,z).  Note that the vectors (x,y) are an orthonormal frame for the tangent space at "z", while
     * "z" itself is an orthonormal frame for the normal space at "z".
     *
     * @param z a point.
     * @return The corresponding frame.
     */
    fun getFrame(z: S2Point): Matrix3x3Double {
        val m = Matrix3x3Double()
        getFrame(z, m)
        return m
    }

    /**
     * Extend the given point "z" on the unit sphere into a right-handed coordinate frame of unit-length column
     * vectors m = (x,y,z).  Note that the vectors (x,y) are an orthonormal frame for the tangent space at "z", while
     * "z" itself is an orthonormal frame for the normal space at "z".
     *
     * @param z a point.
     * @param m The corresponding frame.
     */
    fun getFrame(z: S2Point, m: Matrix3x3Double) {
        requireArgument { isUnitLength(z) }
        m.setCol(2, z)
        m.setCol(1, ortho(z))
        m.setCol(0, m.col(1).crossProd(z));  // Already unit-length.
    }

    // Given an orthonormal basis "m" of column vectors and a point "p", return
// the coordinates of "p" with respect to the basis "m".  The resulting
// point "q" satisfies the identity (m * q == p).
    fun toFrame(m: Matrix3x3Double, p: S2Point): S2Point {
        // The inverse of an orthonormal matrix is its transpose.
        return m.transpose() * p
    }

    // Given an orthonormal basis "m" of column vectors and a point "q" with
// respect to that basis, return the equivalent point "p" with respect to
// the standard axis-aligned basis.  The result satisfies (p == m * q).
    fun fromFrame(m: Matrix3x3Double, q: S2Point): S2Point {
        return m * q
    }

    /**
     * Return true if the points A, B, C are strictly counterclockwise. Return false if the points are clockwise or
     * collinear (i.e. if they are all contained on some great circle).
     *
     * Due to numerical errors, situations may arise that are mathematically impossible, e.g. ABC may be considered
     * strictly CCW while BCA is not. However, the implementation guarantees the following:
     *
     *   If simpleCCW(a,b,c), then !simpleCCW(c,b,a) for all a,b,c.
     *
     * @param a a point
     * @param b a point
     * @param c a point
     *
     * @return true if the points are strictly counterclockwise.
     */
    fun simpleCCW(a: S2Point, b: S2Point, c: S2Point): Boolean {
        // We compute the signed volume of the parallelepiped ABC.  The usual
        // formula for this is (AxB).C, but we compute it here using (CxA).B
        // in order to ensure that ABC and CBA are not both CCW.  This follows
        // from the following identities (which are true numerically, not just
        // mathematically):
        //
        //     (1) x.crossProd(y) == -(y.crossProd(x))
        //     (2) (-x).dotProd(y) == -(x.dotProd(y))

        return c.crossProd(a).dotProd(b) > 0
    }

    fun toDegreesString(p: S2Point): String {
        val s2LatLng = S2LatLng.fromPoint(p)
        return "(" + s2LatLng.latDegrees() + ", " + s2LatLng.lngDegrees() + ")"
    }

}
