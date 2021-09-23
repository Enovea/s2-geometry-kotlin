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
package dilivia.s2.builder.snap

import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLE
import dilivia.math.DoubleType
import dilivia.math.M_SQRT1_2
import dilivia.math.M_SQRT2
import dilivia.s2.S1Angle
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import org.apache.commons.math3.util.FastMath.ceil
import org.apache.commons.math3.util.FastMath.log10
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.round

// A SnapFunction that snaps vertices to S2LatLng E5, E6, or E7 coordinates.
// These coordinates are expressed in degrees multiplied by a power of 10 and
// then rounded to the nearest integer.  For example, in E6 coordinates the
// point (23.12345651, -45.65432149) would become (23123457, -45654321).
//
// The main argument of the SnapFunction is the exponent for the power of 10
// that coordinates should be multipled by before rounding.  For example,
// IntLatLngSnapFunction(7) is a function that snaps to E7 coordinates.  The
// exponent can range from 0 to 10.
//
// Each exponent has a corresponding minimum snap radius, which is simply the
// maximum distance that a vertex can move when snapped.  It is approximately
// equal to 1/sqrt(2) times the nominal point spacing; for example, for
// snapping to E7 the minimum snap radius is (1e-7 / sqrt(2)) degrees.
// You can also set the snap radius to any value larger than this; this can
// result in significant extra simplification (similar to using a larger
// exponent) but does not move vertices unnecessarily.
class IntLatLngSnapFunction(exponent: Int = 7) : SnapFunction() {

    // Snaps vertices to points whose (lat, lng) coordinates are integers after
    // converting to degrees and multiplying by 10 raised to the given exponent.
    // For example, (exponent == 7) yields E7 coordinates.  As a side effect,
    // this method also resets "snap_radius" to the minimum value allowed for
    // this exponent:
    //
    //   set_snap_radius(MinSnapRadiusForExponent(exponent))
    //
    // This means that if you want to use a larger snap radius than the minimum,
    // you must call set_snap_radius() *after* calling set_exponent().
    //
    // REQUIRES: kMinExponent <= exponent <= kMaxExponent
    var exponent: Int = exponent
        set(value) {
            requireGE(value, kMinExponent)
            requireLE(value, kMaxExponent)
            field = value
            snapRadius = minSnapRadiusForExponent(value)

            // Precompute the scale factors needed for snapping.  Note that these
            // calculations need to exactly match the ones in s1angle.h to ensure
            // that the same S2Points are generated.
            var power = 1.0
            repeat(value) { power *= 10.0 }
            fromDegrees = power
            toDegrees = 1.0 / power
        }

    // Defines the snap radius to be used (see s2builder.h).  The snap radius
    // must be at least the minimum value for the current exponent(), but larger
    // values can also be used (e.g., to simplify the geometry).
    //
    // REQUIRES: snap_radius >= MinSnapRadiusForExponent(exponent())
    // REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
    var snapRadius: S1Angle = minSnapRadiusForExponent(exponent)
        set(value) {
            requireGE(value, minSnapRadiusForExponent(exponent))
            requireLE(value, kMaxSnapRadius())
            field = value
        }

    override fun snapRadius(): S1Angle = snapRadius.clone()

    var fromDegrees: Double = 0.0
    var toDegrees: Double = 0.0

    init {
        this.exponent = exponent
    }

    // For IntLatLng snapping, the minimum separation between vertices depends on
    // exponent() and snap_radius().  It can vary between snap_radius()
    // and snap_radius().
    override fun minVertexSeparation(): S1Angle {
        // We have two bounds for the minimum vertex separation: one is proportional
        // to snap_radius, and one is equal to snap_radius minus a constant.  These
        // bounds give the best results for small and large snap radii respectively.
        // We return the maximum of the two bounds.
        //
        // 1. Proportional bound: It can be shown that in the plane, the worst-case
        //    configuration has a vertex separation of (sqrt(2) / 3) * snap_radius.
        //    This is verified in the unit test, except that on the sphere the ratio
        //    is slightly smaller (0.471337 vs. 0.471404).  We reduce that value a
        //    bit more below to be conservative.
        //
        // 2. Best asymptotic bound: This bound bound is derived by observing we
        //    only select a new site when it is at least snap_radius() away from all
        //    existing sites, and snapping a vertex can move it by up to
        //    ((1 / sqrt(2)) * to_degrees_) degrees.
        return maxOf(snapRadius * 0.471,        // sqrt(2) / 3 in the plane
                snapRadius - S1Angle.degrees(M_SQRT1_2 * toDegrees))
    }

    // For IntLatLng snapping, the minimum separation between edges and
    // non-incident vertices depends on level() and snap_radius().  It can
    // be as low as 0.222 * snap_radius(), but is typically 0.39 * snap_radius()
    // or more.
    override fun minEdgeVertexSeparation(): S1Angle {
        // Similar to min_vertex_separation(), in this case we have three bounds:
        // one is a constant bound, one is proportional to snap_radius, and one is
        // equal to snap_radius minus a constant.
        //
        // 1. Constant bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of ((1 / sqrt(13)) * to_degrees_) degrees.
        //    The unit test verifies this, except that on the sphere the ratio is
        //    slightly lower when small exponents such as E1 are used
        //    (0.2772589 vs 0.2773501).
        //
        // 2. Proportional bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of (2 / 9) * snap_radius (0.222222222222).  The
        //    unit test verifies this, except that on the sphere the bound can be
        //    slightly worse with large exponents (e.g., E9) due to small numerical
        //    errors (0.222222126756717).
        //
        // 3. Best asymptotic bound: If snap_radius() is large compared to the
        //    minimum snap radius, then the best bound is achieved by 3 sites on a
        //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
        //    apart (see S2CellIdSnapFunction::min_edge_vertex_separation).  This
        //    bound approaches 0.5 * snap_radius as the snap radius becomes large
        //    relative to the grid spacing.

        val vertexSep = minVertexSeparation()
        return maxOf(S1Angle.degrees(toDegrees) * 0.277,  // 1/sqrt(13) in the plane
                maxOf(snapRadius * 0.222,               // 2/9 in the plane
                        (vertexSep / snapRadius) * vertexSep * 0.5))
    }

    override fun snapPoint(point: S2Point): S2Point {
        requireGE(exponent, 0)  // Make sure the snap function was initialized.
        val input = S2LatLng.fromPoint(point)
        val lat = round(input.lat().degrees() * fromDegrees)
        val lng = round(input.lng().degrees() * fromDegrees)
        return S2LatLng.fromDegrees(lat * toDegrees, lng * toDegrees).toPoint()
    }

    override fun clone(): SnapFunction {
        val function = IntLatLngSnapFunction(exponent)
        function.snapRadius = snapRadius
        return function
    }

    companion object {

        // The minum exponent supported for snapping.
        const val kMinExponent = 0

        // The maximum exponent supported for snapping.
        const val kMaxExponent = 10

        // Returns the minimum allowable snap radius for the given exponent
        // (approximately equal to (pow(10, -exponent) / sqrt(2)) degrees).
        fun minSnapRadiusForExponent(exponent: Int): S1Angle {
            // snap_radius() needs to be an upper bound on the true distance that a
            // point can move when snapped, taking into account numerical errors.
            //
            // The maximum errors in latitude and longitude can be bounded as
            // follows (as absolute errors in terms of DBL_EPSILON):
            //
            //                                      Latitude      Longitude
            // Convert to S2LatLng:                    1.000          1.000
            // Convert to degrees:                     1.032          2.063
            // Scale by 10**exp:                       0.786          1.571
            // Round to integer: 0.5 * S1Angle::Degrees(to_degrees_)
            // Scale by 10**(-exp):                    1.375          2.749
            // Convert to radians:                     1.252          1.503
            // ------------------------------------------------------------
            // Total (except for rounding)             5.445          8.886
            //
            // The maximum error when converting the S2LatLng back to an S2Point is
            //
            //   sqrt(2) * (maximum error in latitude or longitude) + 1.5 * DBL_EPSILON
            //
            // which works out to (9 * sqrt(2) + 1.5) * DBL_EPSILON radians.  Finally
            // we need to consider the effect of rounding to integer coordinates
            // (much larger than the errors above), which can change the position by
            // up to (sqrt(2) * 0.5 * to_degrees_) radians.
            var power = 1.0
            repeat(exponent) { power *= 10.0 }
            return (S1Angle.degrees(M_SQRT1_2 / power) +
                    S1Angle.radians((9 * M_SQRT2 + 1.5) * DoubleType.epsilon))
        }

        // Returns the minimum exponent such that vertices will not move by more
        // than "snap_radius".  This can be useful when choosing an appropriate
        // exponent for snapping.  The return value is always a valid exponent
        // (out of range values are silently clamped).
        //
        // If you want to choose the exponent based on a distance, and then use
        // the minimum possible snap radius for that exponent, do this:
        //
        //   IntLatLngSnapFunction f(
        //       IntLatLngSnapFunction::ExponentForMaxSnapRadius(distance));
        fun exponentForMaxSnapRadius(snapRadius: S1Angle): Int {
            var radius = snapRadius
            // When choosing an exponent, we need to acount for the error bound of
            // (9 * sqrt(2) + 1.5) * DBL_EPSILON added by MinSnapRadiusForExponent().
            radius = radius - S1Angle.radians((9 * M_SQRT2 + 1.5) * DoubleType.epsilon)
            radius = maxOf(radius, S1Angle.radians(1e-30))
            val exponent = log10(M_SQRT1_2 / radius.degrees())

            // There can be small errors in the calculation above, so to ensure that
            // this function is the inverse of MinSnapRadiusForExponent() we subtract a
            // small error tolerance.
            return max(kMinExponent, min(kMaxExponent, ceil(exponent - 2 * DoubleType.epsilon).toInt()))
        }

    }

}
