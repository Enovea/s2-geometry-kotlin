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
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.max
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.coords.S2Coords


// A SnapFunction that snaps vertices to S2CellId centers.  This can be useful
// if you want to encode your geometry compactly using S2Polygon::Encode(),
// for example.  You can snap to the centers of cells at any level.
//
// Every snap level has a corresponding minimum snap radius, which is simply
// the maximum distance that a vertex can move when snapped.  It is
// approximately equal to half of the maximum diagonal length for cells at the
// chosen level.  You can also set the snap radius to a larger value; for
// example, you could snap to the centers of leaf cells (1cm resolution) but
// set the snap_radius() to 10m.  This would result in significant extra
// simplification, without moving vertices unnecessarily (i.e., vertices that
// are at least 10m away from all other vertices will move by less than 1cm).
class S2CellIdSnapFunction(level: Int = S2CellId.kMaxLevel) : SnapFunction() {

    // Snaps vertices to S2Cell centers at the given level.  As a side effect,
    // this method also resets "snap_radius" to the minimum value allowed at
    // this level:
    //
    //   set_snap_radius(MinSnapRadiusForLevel(level))
    //
    // This means that if you want to use a larger snap radius than the minimum,
    // you must call set_snap_radius() *after* calling set_level().
    var level: Int = level
        set(value) {
            requireGE(value, 0)
            requireLE(value, S2CellId.kMaxLevel)
            field = value
            snapRadius = minSnapRadiusForLevel(level)
        }

    // Defines the snap radius to be used (see s2builder.h).  The snap radius
    // must be at least the minimum value for the current level(), but larger
    // values can also be used (e.g., to simplify the geometry).
    //
    // REQUIRES: snap_radius >= MinSnapRadiusForLevel(level())
    // REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
    var snapRadius: S1Angle = minSnapRadiusForLevel(level)
        set(value) {
            requireGE(value, minSnapRadiusForLevel(level))
            requireLE(value, kMaxSnapRadius())
            field = value
        }

    override fun snapRadius(): S1Angle = snapRadius.clone()

    // For S2CellId snapping, the minimum separation between vertices depends on
    // level() and snap_radius().  It can vary between 0.5 * snap_radius()
    // and snap_radius().
    override fun minVertexSeparation(): S1Angle {
        // We have three different bounds for the minimum vertex separation: one is
        // a constant bound, one is proportional to snap_radius, and one is equal to
        // snap_radius minus a constant.  These bounds give the best results for
        // small, medium, and large snap radii respectively.  We return the maximum
        // of the three bounds.
        //
        // 1. Constant bound: Vertices are always separated by at least
        //    kMinEdge(level), the minimum edge length for the chosen snap level.
        //
        // 2. Proportional bound: It can be shown that in the plane, the worst-case
        //    configuration has a vertex separation of 2 / sqrt(13) * snap_radius.
        //    This is verified in the unit test, except that on the sphere the ratio
        //    is slightly smaller at cell level 2 (0.54849 vs. 0.55470).  We reduce
        //    that value a bit more below to be conservative.
        //
        // 3. Best asymptotic bound: This bound bound is derived by observing we
        //    only select a new site when it is at least snap_radius() away from all
        //    existing sites, and the site can move by at most 0.5 * kMaxDiag(level)
        //    when snapped.
        val minEdge = S1Angle.radians(S2Coords.projection.kMinEdge.getValue(level))
        val maxDiag = S1Angle.radians(S2Coords.projection.kMaxDiag.getValue(level))
        return max(minEdge, max(snapRadius * 0.548 /* 2 / sqrt(13) in the plane */, snapRadius - maxDiag * 0.5));
    }

    // For S2CellId snapping, the minimum separation between edges and
    // non-incident vertices depends on level() and snap_radius().  It can
    // be as low as 0.219 * snap_radius(), but is typically 0.5 * snap_radius()
    // or more.
    override fun minEdgeVertexSeparation(): S1Angle {
        // Similar to min_vertex_separation(), in this case we have four bounds: a
        // constant bound that holds only at the minimum snap radius, a constant
        // bound that holds for any snap radius, a bound that is proportional to
        // snap_radius, and a bound that is equal to snap_radius minus a constant.
        //
        // 1. Constant bounds:
        //
        //    (a) At the minimum snap radius for a given level, it can be shown that
        //    vertices are separated from edges by at least 0.5 * kMinDiag(level) in
        //    the plane.  The unit test verifies this, except that on the sphere the
        //    worst case is slightly better: 0.5652980068 * kMinDiag(level).
        //
        //    (b) Otherwise, for arbitrary snap radii the worst-case configuration
        //    in the plane has an edge-vertex separation of sqrt(3/19) *
        //    kMinDiag(level), where sqrt(3/19) is about 0.3973597071.  The unit
        //    test verifies that the bound is slighty better on the sphere:
        //    0.3973595687 * kMinDiag(level).
        //
        // 2. Proportional bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of 2 * sqrt(3/247) * snap_radius, which is
        //    about 0.2204155075.  The unit test verifies this, except that on the
        //    sphere the bound is slightly worse for certain large S2Cells: the
        //    minimum ratio occurs at cell level 6, and is about 0.2196666953.
        //
        // 3. Best asymptotic bound: If snap_radius() is large compared to the
        //    minimum snap radius, then the best bound is achieved by 3 sites on a
        //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
        //    apart.  An input edge passing just to one side of the center of the
        //    circle intersects the Voronoi regions of the two end sites but not the
        //    Voronoi region of the center site, and gives an edge separation of
        //    (min_vertex_separation ** 2) / (2 * snap_radius).  This bound
        //    approaches 0.5 * snap_radius for large snap radii, i.e.  the minimum
        //    edge-vertex separation approaches half of the minimum vertex
        //    separation as the snap radius becomes large compared to the cell size.
        val minDiag = S1Angle.radians(S2Coords.projection.kMinDiag.getValue(level))
        if (snapRadius == minSnapRadiusForLevel(level)) {
            // This bound only holds when the minimum snap radius is being used.
            return minDiag * 0.565            // 0.500 in the plane
        }
        // Otherwise, these bounds hold for any snap_radius().
        val vertexSep = minVertexSeparation()
        return max(minDiag * 0.397,          // sqrt(3 / 19) in the plane
                max(snapRadius * 0.219,      // 2 * sqrt(3 / 247) in the plane
                        (vertexSep / snapRadius) * vertexSep * 0.5))
    }

    override fun snapPoint(point: S2Point): S2Point = S2CellId.fromPoint(point).parent(level).toPoint()

    override fun clone(): SnapFunction {
        val function = S2CellIdSnapFunction(level)
        function.snapRadius = snapRadius
        return function
    }

    override fun toString(): String {
        return "S2CellIdSnapFunction(level=$level, snapRadius=$snapRadius)"
    }


    companion object {

        // Returns the minimum allowable snap radius for the given S2Cell level
        // (approximately equal to half of the maximum cell diagonal length).
        fun minSnapRadiusForLevel(level: Int): S1Angle {
            // snap_radius() needs to be an upper bound on the true distance that a
            // point can move when snapped, taking into account numerical errors.
            //
            // The maximum error when converting from an S2Point to an S2CellId is
            // S2::kMaxDiag.deriv() * DBL_EPSILON.  The maximum error when converting an
            // S2CellId center back to an S2Point is 1.5 * DBL_EPSILON.  These add up to
            // just slightly less than 4 * DBL_EPSILON.
            return S1Angle.radians(0.5 * S2Coords.projection.kMaxDiag.getValue(level) + 4 * DoubleType.epsilon)
        }

        // Returns the minimum S2Cell level (i.e., largest S2Cells) such that
        // vertices will not move by more than "snap_radius".  This can be useful
        // when choosing an appropriate level to snap to.  The return value is
        // always a valid level (out of range values are silently clamped).
        //
        // If you want to choose the snap level based on a distance, and then use
        // the minimum possible snap radius for the chosen level, do this:
        //
        //   S2CellIdSnapFunction f(
        //       S2CellIdSnapFunction::LevelForMaxSnapRadius(distance));
        fun levelForMaxSnapRadius(snapRadius: S1Angle): Int {
            // When choosing a level, we need to acount for the error bound of
            // 4 * DBL_EPSILON that is added by MinSnapRadiusForLevel().
            return S2Coords.projection.kMaxDiag.getLevelForMaxValue(2.0 * (snapRadius.radians - 4 * DoubleType.epsilon))
        }

    }

}
