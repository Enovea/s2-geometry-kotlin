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

import dilivia.PreConditions.checkLE
import dilivia.s2.S1Angle
import dilivia.s2.S2Point


// A SnapFunction restricts the locations of the output vertices.  For
// example, there are predefined snap functions that require vertices to be
// located at S2CellId centers or at E5/E6/E7 coordinates.  The SnapFunction
// can also specify a minimum spacing between vertices (the "snap radius").
//
// A SnapFunction defines the following methods:
//
// 1. The SnapPoint() method, which snaps a point P to a nearby point (the
//    "candidate snap site").  Any point may be returned, including P
//    itself (this is the "identity snap function").
//
// 2. "snap_radius", the maximum distance that vertices can move when
//    snapped.  The snap_radius must be at least as large as the maximum
//    distance between P and SnapPoint(P) for any point P.
//
// 3. "max_edge_deviation", the maximum distance that edges can move when
//    snapped.  It is slightly larger than "snap_radius" because when a
//    geodesic edge is snapped, the center of the edge moves further than
//    its endpoints.  This value is computed automatically by S2Builder.
//
// 4. "min_vertex_separation", the guaranteed minimum distance between
//    vertices in the output.  This is generally a fraction of
//    "snap_radius" where the fraction depends on the snap function.
//
// 5. A "min_edge_vertex_separation", the guaranteed minimum distance
//    between edges and non-incident vertices in the output.  This is
//    generally a fraction of "snap_radius" where the fraction depends on
//    the snap function.
//
// It is important to note that SnapPoint() does not define the actual
// mapping from input vertices to output vertices, since the points it
// returns (the candidate snap sites) are further filtered to ensure that
// they are separated by at least the snap radius.  For example, if you
// specify E7 coordinates (2cm resolution) and a snap radius of 10m, then a
// subset of points returned by SnapPoint will be chosen (the "snap sites"),
// and each input vertex will be mapped to the closest site.  Therefore you
// cannot assume that P is necessarily snapped to SnapPoint(P).
//
// S2Builder makes the following guarantees:
//
// 1. Every vertex is at a location returned by SnapPoint().
//
// 2. Vertices are within "snap_radius" of the corresponding input vertex.
//
// 3. Edges are within "max_edge_deviation" of the corresponding input edge
//    (a distance slightly larger than "snap_radius").
//
// 4. Vertices are separated by at least "min_vertex_separation"
//    (a fraction of "snap_radius" that depends on the snap function).
//
// 5. Edges and non-incident vertices are separated by at least
//    "min_edge_vertex_separation" (a fraction of "snap_radius").
//
// 6. Vertex and edge locations do not change unless one of the conditions
//    above is not already met (idempotency / stability).
//
// 7. The topology of the input geometry is preserved (up to the creation
//    of degeneracies).  This means that there exists a continuous
//    deformation from the input to the output such that no vertex
//    crosses an edge.
abstract class SnapFunction : Cloneable {

    // The maximum distance that vertices can move when snapped.
    //
    // If the snap radius is zero, then vertices are snapped together only if
    // they are identical.  Edges will not be snapped to any vertices other
    // than their endpoints, even if there are vertices whose distance to the
    // edge is zero, unless split_crossing_edges() is true.
    //
    // REQUIRES: snap_radius() <= kMaxSnapRadius
    abstract fun snapRadius(): S1Angle

    // The maximum distance that the center of an edge can move when snapped.
    // This is slightly larger than "snap_radius" because when a geodesic edge
    // is snapped, the center of the edge moves further than its endpoints.
    fun maxEdgeDeviation(): S1Angle {
        // We want max_edge_deviation() to be large enough compared to snap_radius()
        // such that edge splitting is rare.
        //
        // Using spherical trigonometry, if the endpoints of an edge of length L
        // move by at most a distance R, the center of the edge moves by at most
        // asin(sin(R) / cos(L / 2)).  Thus the (max_edge_deviation / snap_radius)
        // ratio increases with both the snap radius R and the edge length L.
        //
        // We arbitrarily limit the edge deviation to be at most 10% more than the
        // snap radius.  With the maximum allowed snap radius of 70 degrees, this
        // means that edges up to 30.6 degrees long are never split.  For smaller
        // snap radii, edges up to 49 degrees long are never split.  (Edges of any
        // length are not split unless their endpoints move far enough so that the
        // actual edge deviation exceeds the limit; in practice, splitting is rare
        // even with long edges.)  Note that it is always possible to split edges
        // when max_edge_deviation() is exceeded; see MaybeAddExtraSites().
        checkLE(snapRadius(), kMaxSnapRadius())
        return snapRadius() * kMaxEdgeDeviationRatio
    }

    // The guaranteed minimum distance between vertices in the output.
    // This is generally some fraction of "snap_radius".
    abstract fun minVertexSeparation(): S1Angle

    // The guaranteed minimum spacing between edges and non-incident vertices
    // in the output.  This is generally some fraction of "snap_radius".
    abstract fun minEdgeVertexSeparation(): S1Angle

    // Returns a candidate snap site for the given point.  The final vertex
    // locations are a subset of the snap sites returned by this function
    // (spaced at least "min_vertex_separation" apart).
    //
    // The only requirement is that SnapPoint(x) must return a point whose
    // distance from "x" is no greater than "snap_radius".
    abstract fun snapPoint(point: S2Point): S2Point

    // Returns a deep copy of this SnapFunction.
    abstract override fun clone(): SnapFunction

    companion object {
        private val kMaxEdgeDeviationRatio = 1.1

        // The maximum supported snap radius (equivalent to about 7800km).
        // The maximum snap radius is just large enough to support snapping to
        // S2CellId level 0.  It is equivalent to 7800km on the Earth's surface.
        // This value can't be larger than 85.7 degrees without changing the code
        // related to min_edge_length_to_split_ca_, and increasing it to 90 degrees
        // or more would most likely require significant changes to the algorithm.
        fun kMaxSnapRadius(): S1Angle = S1Angle.degrees(70)

    }
};
