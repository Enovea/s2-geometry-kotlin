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
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2PointUtil.robustCrossProd
import org.apache.commons.math3.util.FastMath.atan
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.sqrt
import org.apache.commons.math3.util.FastMath.tan

/**
 *
 */
object S2Measures {

    // Return the interior angle at the vertex B in the triangle ABC.  The
    // return value is always in the range [0, Pi].  All points should be
    // normalized.  Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
    //
    // The angle is undefined if A or C is diametrically opposite from B, and
    // becomes numerically unstable as the length of edge AB or BC approaches
    // 180 degrees.
    @JvmStatic
    fun angle(a: S2Point, b: S2Point, c: S2Point): Double {
        // RobustCrossProd() is necessary to get good accuracy when two of the input
        // points are very close together.
        return robustCrossProd(a, b).angle(robustCrossProd(c, b))
    }

    // Return the exterior angle at vertex B in the triangle ABC.  The return
    // value is positive if ABC is counterclockwise and negative otherwise.  If
    // you imagine an ant walking from A to B to C, this is the angle that the
    // ant turns at vertex B (positive = left = CCW, negative = right = CW).
    // This quantity is also known as the "geodesic curvature" at B.
    //
    // Ensures that TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all distinct
    // a,b,c. The result is undefined if (a == b || b == c), but is either
    // -Pi or Pi if (a == c).  All points should be normalized.
    @JvmStatic
    fun turnAngle(a: S2Point, b: S2Point, c: S2Point): Double {
        // We use RobustCrossProd() to get good accuracy when two points are very
        // close together, and Sign() to ensure that the sign is correct for
        // turns that are close to 180 degrees.
        //
        // Unfortunately we can't save RobustCrossProd(a, b) and pass it as the
        // optional 4th argument to Sign(), because Sign() requires a.CrossProd(b)
        // exactly (the robust version differs in magnitude).
        val angle = robustCrossProd(a, b).angle(robustCrossProd(b, c))

        // Don't return Sign() * angle because it is legal to have (a == c).
        return if (S2Predicates.sign(a, b, c) > 0) angle else -angle
    }

    // Return the area of triangle ABC.  This method combines two different
    // algorithms to get accurate results for both large and small triangles.
    // The maximum error is about 5e-15 (about 0.25 square meters on the Earth's
    // surface), the same as GirardArea() below, but unlike that method it is
    // also accurate for small triangles.  Example: when the true area is 100
    // square meters, Area() yields an error about 1 trillion times smaller than
    // GirardArea().
    //
    // All points should be unit length, and no two points should be antipodal.
    // The area is always positive.
    @JvmStatic
    fun area(a: S2Point, b: S2Point, c: S2Point): Double {
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        requireArgument { isUnitLength(c) }
        // This method is based on l'Huilier's theorem,
        //
        //   tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
        //
        // where E is the spherical excess of the triangle (i.e. its area),
        //       a, b, c, are the side lengths, and
        //       s is the semiperimeter (a + b + c) / 2 .
        //
        // The only significant source of error using l'Huilier's method is the
        // cancellation error of the terms (s-a), (s-b), (s-c).  This leads to a
        // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c).  This compares
        // to a relative error of about 1e-15 / E using Girard's formula, where E is
        // the true area of the triangle.  Girard's formula can be even worse than
        // this for very small triangles, e.g. a triangle with a true area of 1e-30
        // might evaluate to 1e-5.
        //
        // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
        // dmin = min(s-a, s-b, s-c).  This basically includes all triangles
        // except for extremely long and skinny ones.
        //
        // Since we don't know E, we would like a conservative upper bound on
        // the triangle area in terms of s and dmin.  It's possible to show that
        // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
        // Using this, it's easy to show that we should always use l'Huilier's
        // method if dmin >= k2 * s^5, where k2 is about 1e-2.  Furthermore,
        // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
        // k3 is about 0.1.  Since the best case error using Girard's formula
        // is about 1e-15, this means that we shouldn't even consider it unless
        // s >= 3e-4 or so.
        //
        // TODO(ericv): Implement rigorous error bounds (analysis already done).
        val sa = b.angle(c)
        val sb = c.angle(a)
        val sc = a.angle(b)
        val s = 0.5 * (sa + sb + sc)
        if (s >= 3e-4) {
            // Consider whether Girard's formula might be more accurate.
            val s2 = s * s
            val dmin = s - max(sa, max(sb, sc))
            if (dmin < 1e-2 * s * s2 * s2) {
                // This triangle is skinny enough to consider using Girard's formula.
                // We increase the area by the approximate maximum error in the Girard
                // calculation in order to ensure that this test is conservative.
                val area = girardArea(a, b, c)
                if (dmin < s * (0.1 * (area + 5e-15))) return area
            }
        }
        // Use l'Huilier's formula.
        return 4 * atan(sqrt(max(0.0, tan(0.5 * s) * tan(0.5 * (s - sa)) * tan(0.5 * (s - sb)) * tan(0.5 * (s - sc)))))
    }

    // Return the area of the triangle computed using Girard's formula.  All
    // points should be unit length, and no two points should be antipodal.
    //
    // This method is about twice as fast as Area() but has poor relative
    // accuracy for small triangles.  The maximum error is about 5e-15 (about
    // 0.25 square meters on the Earth's surface) and the average error is about
    // 1e-15.  These bounds apply to triangles of any size, even as the maximum
    // edge length of the triangle approaches 180 degrees.  But note that for
    // such triangles, tiny perturbations of the input points can change the
    // true mathematical area dramatically.
    @JvmStatic
    fun girardArea(a: S2Point, b: S2Point, c: S2Point): Double {
        // This is equivalent to the usual Girard's formula but is slightly more
        // accurate, faster to compute, and handles a == b == c without a special
        // case.  RobustCrossProd() is necessary to get good accuracy when two of
        // the input points are very close together.

        val ab = robustCrossProd(a, b)
        val bc = robustCrossProd(b, c)
        val ac = robustCrossProd(a, c)
        return max(0.0, ab.angle(ac) - ab.angle(bc) + bc.angle(ac));
    }

    // Like Area(), but returns a positive value for counterclockwise triangles
    // and a negative value otherwise.
    @JvmStatic
    fun signedArea(a: S2Point, b: S2Point, c: S2Point): Double = S2Predicates.sign(a, b, c) * area(a, b, c)

}
