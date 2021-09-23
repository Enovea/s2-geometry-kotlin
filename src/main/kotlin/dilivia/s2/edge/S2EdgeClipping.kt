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
package dilivia.s2.edge

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireNE
import dilivia.math.DoubleType
import dilivia.math.M_SQRT1_2
import dilivia.math.M_SQRT2
import dilivia.math.R1Interval
import dilivia.math.R2Rect
import dilivia.math.ldexp
import dilivia.math.vectors.R2Point
import dilivia.math.vectors.times
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.coords.S2Coords
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.signum

/**
 * S2PointUVW is used to document that a given S2Point is expressed in the (u,v,w) coordinates of some cube face.
 * */
typealias S2PointUVW = S2Point

/**
 * FaceSegment represents an edge AB clipped to an S2 cube face. It is represented by a face index and a pair of (u,v)
 * coordinates.
 */
data class FaceSegment(
    val face: Int,
    val a: R2Point,
    val b: R2Point
)

/** Mutable FaceSegment List. */
typealias FaceSegmentList = MutableList<FaceSegment>

/**
 * Defines a collection of functions for:
 *
 *  (1) Robustly clipping geodesic edges to the faces of the S2 biunit cube (see S2Coords), and
 *  (2) Robustly clipping 2D edges against 2D rectangles.
 *
 * These functions can be used to efficiently find the set of S2CellIds thatvare intersected by a geodesic edge
 * (e.g., see S2CrossingEdgeQuery).
 *
 * This class is a port of the s2edge_clipping methods of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @constructor Default constructor. must be followed by a call to init().
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2EdgeClipping {

    /** Logger. */
    private val logger = KotlinLogging.logger {  }

    // The maximum error in the vertices returned by GetFaceSegments and ClipToFace (compared to an exact calculation):

    /**
     * kFaceClipErrorRadians is the maximum angle between a returned vertex and the nearest point on the exact edge AB.
     * It is equal to the maximum directional error in S2::RobustCrossProd, plus the error when projecting points onto
     * a cube face.
     */
    val kFaceClipErrorRadians = 3 * DoubleType.epsilon

    /**
     * kFaceClipErrorDist is the same angle expressed as a maximum distance in (u,v)-space. In other words, a returned
     * vertex is at most this far from the exact edge AB projected into (u,v)-space.
     */
    val kFaceClipErrorUVDist = 9 * DoubleType.epsilon

    /**
     * kFaceClipErrorUVCoord is the same angle expressed as the maximum error in an individual u- or v-coordinate.
     * In other words, for each returned vertex there is a point on the exact edge AB whose u- and v-coordinates differ
     * from the vertex by at most this amount.
     */
    val kFaceClipErrorUVCoord = 9 * M_SQRT1_2 * DoubleType.epsilon

    /**
     * The maximum error in intersectRect. If some point of AB is inside the rectangle by at least this distance, the
     * result is guaranteed to be true; if all points of AB are outside the rectangle by at least this distance,
     * the result is guaranteed to be false. This bound assumes that "rect" is a subset of the rectangle [-1,1]x[-1,1]
     * or extends slightly outside it (e.g., by 1e-10 or less).
     */
    val kIntersectsRectErrorUVDist = 3 * M_SQRT2 * DoubleType.epsilon

    // The maximum error in the vertices generated by ClipEdge and the bounds generated by ClipEdgeBound (compared to
    // an exact calculation):

    /**
     * kEdgeClipErrorUVCoord is the maximum error in a u- or v-coordinate compared to the exact result, assuming that
     * the points A and B are in the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
     */
    val kEdgeClipErrorUVCoord = 2.25 * DoubleType.epsilon

    /**
     * kEdgeClipErrorUVDist is the maximum distance from a clipped point to the corresponding exact result. It is equal
     * to the error in a single coordinate because at most one coordinate is subject to error.
     */
    val kEdgeClipErrorUVDist = 2.25 * DoubleType.epsilon

    /**
     * Subdivides the given edge AB at every point where it crosses the boundary between two S2 cube faces and returns
     * the corresponding FaceSegments.
     * The segments are returned in order from A toward B. The input points must be unit length.
     *
     * This method guarantees that the returned segments form a continuous path from A to B, and that all vertices are
     * within kFaceClipErrorUVDist of the line AB. All vertices lie within the [-1,1]x[-1,1] cube face rectangles.
     * The results are consistent with S2Predicates.sign(), i.e. the edge is well-defined even its endpoints are
     * antipodal.
     *
     * TODO(ericv): Extend the implementation of S2::RobustCrossProd so that this statement is true.
     *
     * @param a The first point of the edge AB (require to be unit length).
     * @param b The second point of the edge AB (require to be unit length).
     * @param segments The container to fill with the face segments (it is cleared before computing face segments).
    */
    fun getFaceSegments(a: S2Point, b: S2Point, segments: FaceSegmentList) {
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        segments.clear()

        // Fast path: both endpoints are on the same face.
        var (aFace, aU, aV) = S2Coords.xyzToFaceUv(a)
        var (bFace, bU, bV) = S2Coords.xyzToFaceUv(b)
        val aUV = R2Point(aU, aV)
        val bUV = R2Point(bU, bV)
        var segment = FaceSegment(face = aFace, a = aUV.clone(), b = bUV.clone())
        if (aFace == bFace) {
            segments.add(segment)
            return
        }
        // Starting at A, we follow AB from face to face until we reach the face
        // containing B.  The following code is designed to ensure that we always
        // reach B, even in the presence of numerical errors.
        //
        // First we compute the normal to the plane containing A and B.  This normal
        // becomes the ultimate definition of the line AB; it is used to resolve all
        // questions regarding where exactly the line goes.  Unfortunately due to
        // numerical errors, the line may not quite intersect the faces containing
        // the original endpoints.  We handle this by moving A and/or B slightly if
        // necessary so that they are on faces intersected by the line AB.
        val ab = S2PointUtil.robustCrossProd(a, b)
        aFace = moveOriginToValidFace(aFace, a, ab, aUV)
        bFace= moveOriginToValidFace(bFace, b, -ab, bUV)

        // Now we simply follow AB from face to face until we reach B.
        segment = FaceSegment(face = aFace, a = aUV.clone(), b = bUV.clone())
        val bSaved = segment.b
        var face = aFace

        while (face != bFace) {
            // Complete the current segment by finding the point where AB exits the
            // current face.
            val n = S2Coords.faceXyzToUvw(face, ab)
            val exitAxis = getExitAxis(n)
            segment = segment.copy(b = getExitPoint(n, exitAxis))
            segments.add(segment)

            // Compute the next face intersected by AB, and translate the exit point
            // of the current segment into the (u,v) coordinates of the next face.
            // This becomes the first point of the next segment.
            val exitXyz = S2Coords.faceUvToXyz(face, segment.b)
            face = getNextFace(face, segment.b, exitAxis, n, bFace)
            val exitUvw = S2Coords.faceXyzToUvw(face, exitXyz)
            segment = segment.copy(face = face, a = R2Point(exitUvw[0], exitUvw[1]))
        }
        // Finish the last segment.
        segment = segment.copy(b = bSaved)
        segments.add(segment)
    }

    /**
     * Given an edge AB and a face, returns the (u,v) coordinates for the portion of AB that intersects that face.
     * This method guarantees that the clipped vertices lie within the [-1,1]x[-1,1] cube face rectangle and are within
     * kFaceClipErrorUVDist of the line AB, but the results may differ from those produced by GetFaceSegments.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param face The target face.
     * @return null if AB does not intersect the given face and the clipped edge otherwise.
    */
    fun clipToFace(a: S2Point, b: S2Point, face: Int): Pair<R2Point, R2Point>? = clipToPaddedFace(a, b, face, 0.0)

    /**
     * Like clipToFace, but rather than clipping to the square [-1,1]x[-1,1] in (u,v) space, this method clips to
     * [-R,R]x[-R,R] where R=(1+padding).
     *
     * @param aXyz The first endpoint of the edge AB to clip.
     * @param bXyz The second endpoint of the edge AB to clip.
     * @param face The target face.
     * @param padding The face padding.
     *
     * @return null if AB does not intersect the given face and the clipped edge otherwise.
     */
    fun clipToPaddedFace(aXyz: S2Point, bXyz: S2Point, face: Int, padding: Double): Pair<R2Point, R2Point>? {
        requireGE(padding, 0.0)

        val aUv: R2Point
        val bUv: R2Point
        // Fast path: both endpoints are on the given face.
        if (S2Coords.face(aXyz) == face && S2Coords.face(bXyz) == face) {
            aUv = S2Coords.validFaceXyzToUv(face, aXyz)
            bUv = S2Coords.validFaceXyzToUv(face, bXyz)
            return Pair(aUv, bUv)
        }
        // Convert everything into the (u,v,w) coordinates of the given face.  Note
        // that the cross product *must* be computed in the original (x,y,z)
        // coordinate system because RobustCrossProd (unlike the mathematical cross
        // product) can produce different results in different coordinate systems
        // when one argument is a linear multiple of the other, due to the use of
        // symbolic perturbations.
        var n: S2PointUVW = S2Coords.faceXyzToUvw(face, S2PointUtil.robustCrossProd(aXyz, bXyz))
        val a: S2PointUVW = S2Coords.faceXyzToUvw(face, aXyz)
        val b: S2PointUVW = S2Coords.faceXyzToUvw(face, bXyz)

        // Padding is handled by scaling the u- and v-components of the normal.
        // Letting R=1+padding, this means that when we compute the dot product of
        // the normal with a cube face vertex (such as (-1,-1,1)), we will actually
        // compute the dot product with the scaled vertex (-R,-R,1).  This allows
        // methods such as IntersectsFace(), GetExitAxis(), etc, to handle padding
        // with no further modifications.
        val scaleUv = 1 + padding;
        val scaledN: S2PointUVW = S2Point(scaleUv * n[0], scaleUv * n[1], n[2])
        if (!intersectsFace(scaledN)) return null

        // TODO(ericv): This is a temporary hack until I rewrite S2::RobustCrossProd;
        // it avoids loss of precision in normalize() when the vector is so small
        // that it underflows.
        if (max(abs(n[0]), max(abs(n[1]), abs(n[2]))) < ldexp(1.0, -511)) {
            logger.warn { "N = $n is so small that it underflows." }
            n = n * ldexp(1.0, 563);
        }  // END OF HACK

        n.normalize()
        val aTangent: S2PointUVW = n.crossProd(a)
        val bTangent: S2Point = b.crossProd(n)
        // As described above, if the sum of the scores from clipping the two
        // endpoints is 3 or more, then the segment does not intersect this face.
        val (aScore, aClippedUv) = clipDestination(b, a, -scaledN, bTangent, aTangent, scaleUv)
        val (bScore, bClippedUv) = clipDestination(a, b, scaledN, aTangent, bTangent, scaleUv)
        return if (aScore + bScore < 3) Pair(aClippedUv, bClippedUv) else null
    }

    /**
     * Indicates if the edge AB intersects the rectangle rect.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param rect A rectangle.
     *
     * @return true if the edge AB intersects the given (closed) rectangle to within the error bound below.
     */
    fun intersectsRect(a: R2Point, b: R2Point, rect: R2Rect): Boolean {
        // First check whether the bound of AB intersects "rect".
        val bound = R2Rect.fromPointPair(a, b)
        if (!rect.intersects(bound)) return false

        // Otherwise AB intersects "rect" if and only if all four vertices of "rect"
        // do not lie on the same side of the extended line AB.  We test this by
        // finding the two vertices of "rect" with minimum and maximum projections
        // onto the normal of AB, and computing their dot products with the edge
        // normal.
        val n = (b - a).ortho()
        val i = if (n[0] >= 0) 1 else 0
        val j = if (n[1] >= 0) 1 else 0
        val max = n.dotProd(rect.getVertex(i, j) - a)
        val min = n.dotProd(rect.getVertex(1 - i, 1 - j) - a)
        return (max >= 0) && (min <= 0)
    }

    /**
     * Given an edge AB, returns the portion of AB that is contained by the given rectangle "clip".
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param rect A rectangle.
     *
     * @return The clipped edge or null if the edge does not intersect the rectangle.
     */
    fun clipEdge(a: R2Point, b: R2Point, clip: R2Rect): Pair<R2Point, R2Point>? {
        // Compute the bounding rectangle of AB, clip it, and then extract the new
        // endpoints from the clipped bound.
        val bound = R2Rect.fromPointPair(a, b)
        if (clipEdgeBound(a, b, clip, bound)) {
            val ai = if (a[0] > b[0]) 1 else 0
            val aj = if (a[1] > b[1]) 1 else 0
            val aClipped = bound.getVertex(ai, aj)
            val bClipped = bound.getVertex(1 - ai, 1 - aj)
            return Pair(aClipped, bClipped)
        }
        return null
    }

    /**
     * Given an edge AB and a rectangle "clip", returns the bounding rectangle of the portion of AB intersected by
     * "clip". The resulting bound may be empty. This is a convenience function built on top of clipEdgeBound.
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param clip The target rectangle.
     *
     * @return the bounding rectangle of the portion of AB intersected by "clip"
     */
    fun getClippedEdgeBound(a: R2Point, b: R2Point, clip: R2Rect): R2Rect {
        val bound = R2Rect.fromPointPair(a, b)
        if (clipEdgeBound(a, b, clip, bound)) return bound
        return R2Rect.empty()
    }

    /**
     * This function can be used to clip an edge AB to sequence of rectangles efficiently. It represents the clipped
     * edges by their bounding boxes rather than as a pair of endpoints. Specifically, let A'B' be some portion of an
     * edge AB, and let "bound" be a tight bound of A'B'. This function updates "bound" (in place) to be a tight bound
     * of A'B' intersected with a given rectangle "clip". If A'B' does not intersect "clip", returns false and does not
     * necessarily update "bound".
     *
     * REQUIRES: "bound" is a tight bounding rectangle for some portion of AB. (This condition is automatically
     * satisfied if you start with the bounding box of AB and clip to a sequence of rectangles, stopping when the method
     * returns false.)
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param clip The target rectangle.
     * @param bound The computed bound.
     * @return false if the edge does not intersect the "clip" rectangle.
     */
    fun clipEdgeBound(a: R2Point, b: R2Point, clip: R2Rect, bound: R2Rect): Boolean {
        // "diag" indicates which diagonal of the bounding box is spanned by AB: it
        // is 0 if AB has positive slope, and 1 if AB has negative slope.  This is
        // used to determine which interval endpoints need to be updated each time
        // the edge is clipped.
        val diag = if ((a[0] > b[0]) != (a[1] > b[1])) 1 else 0
        return (clipBoundAxis(a[0], b[0], bound[0], a[1], b[1], bound[1], diag, clip[0]) &&
                clipBoundAxis(a[1], b[1], bound[1], a[0], b[0], bound[0], diag, clip[1]))
    }

    /**
     * Given a value x that is some linear combination of a and b, returns the value x1 that is the same linear
     * combination of a1 and b1. This function makes the following guarantees:
     *  - If x == a, then x1 = a1 (exactly).
     *  - If x == b, then x1 = b1 (exactly).
     *  - If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
     *
     *  REQUIRES: a != b
     *
     * @param a The first endpoint of the edge AB.
     * @param b The second endpoint of the edge AB.
     * @param a1 The first endpoint of the edge A'B'.
     * @param b1 The second endpoint of the edge A'B'.
     */
    fun interpolateDouble(x: Double, a: Double, b: Double, a1: Double, b1: Double): Double {
        requireNE(a, b)
        // To get results that are accurate near both A and B, we interpolate
        // starting from the closer of the two points.
        if (abs(a - x) <= abs(b - x)) {
            return a1 + (b1 - a1) * (x - a) / (b - a)
        } else {
            return b1 + (a1 - b1) * (x - b) / (a - b)
        }
    }

    //////////////////////// Internal methods //////////////////////////////////

    // The three functions below all compare a sum (u + v) to a third value w.
    // They are implemented in such a way that they produce an exact result even
    // though all calculations are done with ordinary floating-point operations.
    // Here are the principles on which these functions are based:
    //
    // A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
    //
    // B. If u + v < w in exact arithmetic, then at least one of the following
    //    expressions is true in floating-point:
    //       u + v < w
    //       u < w - v
    //       v < w - u
    //
    //    Proof: By rearranging terms and substituting ">" for "<", we can assume
    //    that all values are non-negative.  Now clearly "w" is not the smallest
    //    value, so assume WLOG that "u" is the smallest.  We want to show that
    //    u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
    //    exact since the result is smaller in magnitude than either input value,
    //    so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
    //    (even in floating point), so the result also holds.

    // Return true if u + v == w exactly.
    private fun sumEquals(u: Double, v: Double, w: Double): Boolean = (u + v == w) && (u == w - v) && (v == w - u)

    // Return true if a given directed line L intersects the cube face F.  The
    // line L is defined by its normal N in the (u,v,w) coordinates of F.
    private fun intersectsFace(n: S2PointUVW): Boolean {
        // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
        // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
        // and (-1,1,1) do not all have the same sign.  This is true exactly when
        // |Nu| + |Nv| >= |Nw|.  The code below evaluates this expression exactly
        // (see comments above).
        val u = abs(n[0])
        val v = abs(n[1])
        val w = abs(n[2])
        // We only need to consider the cases where u or v is the smallest value,
        // since if w is the smallest then both expressions below will have a
        // positive LHS and a negative RHS.
        return (v >= w - u) && (u >= w - v)
    }

    // Given a directed line L intersecting a cube face F, return true if L
    // intersects two opposite edges of F (including the case where L passes
    // exactly through a corner vertex of F).  The line L is defined by its
    // normal N in the (u,v,w) coordinates of F.
    private fun intersectsOppositeEdges(n: S2PointUVW): Boolean {
        // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
        // and only exactly two of the corner vertices lie on each side of L.  This
        // is true exactly when ||Nu| - |Nv|| >= |Nw|.  The code below evaluates this
        // expression exactly (see comments above).
        val u = abs(n[0])
        val v = abs(n[1])
        val w = abs(n[2])
        // If w is the smallest, the following line returns an exact result.
        if (abs(u - v) != w) return abs(u - v) >= w
        // Otherwise u - v = w exactly, or w is not the smallest value.  In either
        // case the following line returns the correct result.
        return if (u >= v) (u - w >= v) else (v - w >= u)
    }

    // Given cube face F and a directed line L (represented by its CCW normal N in
    // the (u,v,w) coordinates of F), compute the axis of the cube face edge where
    // L exits the face: return 0 if L exits through the u=-1 or u=+1 edge, and 1
    // if L exits through the v=-1 or v=+1 edge.  Either result is acceptable if L
    // exits exactly through a corner vertex of the cube face.
    private fun getExitAxis(n: S2PointUVW): Int {
        requireArgument { intersectsFace(n) }
        if (intersectsOppositeEdges(n)) {
            // The line passes through through opposite edges of the face.
            // It exits through the v=+1 or v=-1 edge if the u-component of N has a
            // larger absolute magnitude than the v-component.
            return if (abs(n[0]) >= abs(n[1])) 1 else 0
        } else {
            // The line passes through through two adjacent edges of the face.
            // It exits the v=+1 or v=-1 edge if an even number of the components of N
            // are negative.  We test this using signbit() rather than multiplication
            // to avoid the possibility of underflow.
            checkState { n[0] != 0.0 && n[1] != 0.0 && n[2] != 0.0 }
            return if (!(signbit(n[0]) xor signbit(n[1]) xor signbit(n[2]))) 1 else 0
        }
    }

    private fun signbit(v: Double): Boolean = signum(v) < 0.0

    // Given a cube face F, a directed line L (represented by its CCW normal N in
    // the (u,v,w) coordinates of F), and result of GetExitAxis(N), return the
    // (u,v) coordinates of the point where L exits the cube face.
    private fun getExitPoint(n: S2PointUVW, axis: Int): R2Point {
        if (axis == 0) {
            val u = if (n[1] > 0) 1.0 else -1.0
            return R2Point(u, (-u * n[0] - n[2]) / n[1])
        } else {
            val v = if (n[0] < 0) 1.0 else -1.0
            return R2Point((-v * n[1] - n[2]) / n[0], v)
        }
    }

    // Given a line segment AB whose origin A has been projected onto a given cube
    // face, determine whether it is necessary to project A onto a different face
    // instead.  This can happen because the normal of the line AB is not computed
    // exactly, so that the line AB (defined as the set of points perpendicular to
    // the normal) may not intersect the cube face containing A.  Even if it does
    // intersect the face, the "exit point" of the line from that face may be on
    // the wrong side of A (i.e., in the direction away from B).  If this happens,
    // we reproject A onto the adjacent face where the line AB approaches A most
    // closely.  This moves the origin by a small amount, but never more than the
    // error tolerances documented in the header file.
    private fun moveOriginToValidFace(face: Int, a: S2Point, ab: S2Point, aUv: R2Point): Int {
        var currentFace = face
        // Fast path: if the origin is sufficiently far inside the face, it is
        // always safe to use it.
        val kMaxSafeUVCoord = 1 - kFaceClipErrorUVCoord
        if (max(abs(aUv[0]), abs(aUv[1])) <= kMaxSafeUVCoord) {
            return currentFace;
        }
        // Otherwise check whether the normal AB even intersects this face.
        val n = S2Coords.faceXyzToUvw(currentFace, ab)
        if (intersectsFace(n)) {
            // Check whether the point where the line AB exits this face is on the
            // wrong side of A (by more than the acceptable error tolerance).
            val exit = S2Coords.faceUvToXyz(currentFace, getExitPoint(n, getExitAxis(n)))
            val aTangent = ab.normalized().crossProd(a)
            if ((exit - a).dotProd(aTangent) >= -kFaceClipErrorRadians) {
                return currentFace;  // We can use the given face.
            }
        }
        // Otherwise we reproject A to the nearest adjacent face.  (If line AB does
        // not pass through a given face, it must pass through all adjacent faces.)
        if (abs(aUv[0]) >= abs(aUv[1])) {
            currentFace = S2Coords.uvwFace(currentFace, 0 /*U axis*/, if (aUv[0] > 0.0) 1 else 0)
        } else {
            currentFace = S2Coords.uvwFace(currentFace, 1 /*V axis*/, if (aUv[1] > 0.0) 1 else 0)
        }
        checkState { intersectsFace(S2Coords.faceXyzToUvw(currentFace, ab)) }
        S2Coords.validFaceXyzToUv(currentFace, a, aUv)
        aUv[0] = max(-1.0, min(1.0, aUv[0]))
        aUv[1] = max(-1.0, min(1.0, aUv[1]))
        return currentFace
    }

    // Return the next face that should be visited by GetFaceSegments, given that
    // we have just visited "face" and we are following the line AB (represented
    // by its normal N in the (u,v,w) coordinates of that face).  The other
    // arguments include the point where AB exits "face", the corresponding
    // exit axis, and the "target face" containing the destination point B.
    private fun getNextFace(face: Int, exit: R2Point, axis: Int, n: S2PointUVW, targetFace: Int): Int {
        // We return the face that is adjacent to the exit point along the given
        // axis.  If line AB exits *exactly* through a corner of the face, there are
        // two possible next faces.  If one is the "target face" containing B, then
        // we guarantee that we advance to that face directly.
        //
        // The three conditions below check that (1) AB exits approximately through
        // a corner, (2) the adjacent face along the non-exit axis is the target
        // face, and (3) AB exits *exactly* through the corner.  (The SumEquals()
        // code checks whether the dot product of (u,v,1) and "n" is exactly zero.)
        if (abs(exit[1 - axis]) == 1.0 &&
            S2Coords.uvwFace(face, 1 - axis, if (exit[1 - axis] > 0.0) 1 else 0) == targetFace &&
            sumEquals(exit[0] * n[0], exit[1] * n[1], -n[2])
        ) {
            return targetFace
        }
        // Otherwise return the face that is adjacent to the exit point in the
        // direction of the exit axis.
        return S2Coords.uvwFace(face, axis, if (exit[axis] > 0.0) 1 else 0)
    }

    // This helper function does two things.  First, it clips the line segment AB
    // to find the clipped destination B' on a given face.  (The face is specified
    // implicitly by expressing *all arguments* in the (u,v,w) coordinates of that
    // face.)  Second, it partially computes whether the segment AB intersects
    // this face at all.  The actual condition is fairly complicated, but it turns
    // out that it can be expressed as a "score" that can be computed
    // independently when clipping the two endpoints A and B.  This function
    // returns the score for the given endpoint, which is an integer ranging from
    // 0 to 3.  If the sum of the two scores is 3 or more, then AB does not
    // intersect this face.  See the calling function for the meaning of the
    // various parameters.
    private fun clipDestination(
        a: S2PointUVW, b: S2PointUVW, scaledN: S2PointUVW,
        aTangent: S2PointUVW, bTangent: S2PointUVW,
        scaleUv: Double
    ): Pair<Int, R2Point> {
        requireArgument { intersectsFace(scaledN) }

        var uv: R2Point
        // Optimization: if B is within the safe region of the face, use it.
        val kMaxSafeUVCoord = 1 - kFaceClipErrorUVCoord
        if (b[2] > 0) {
            uv = R2Point(b[0] / b[2], b[1] / b[2])
            if (max(abs(uv[0]), abs(uv[1])) <= kMaxSafeUVCoord)
                return 0 to uv
        }
        // Otherwise find the point B' where the line AB exits the face.
        uv = scaleUv * getExitPoint(scaledN, getExitAxis(scaledN))
        val p = S2PointUVW(uv[0], uv[1], 1.0)

        // Determine if the exit point B' is contained within the segment.  We do this
        // by computing the dot products with two inward-facing tangent vectors at A
        // and B.  If either dot product is negative, we say that B' is on the "wrong
        // side" of that point.  As the point B' moves around the great circle AB past
        // the segment endpoint B, it is initially on the wrong side of B only; as it
        // moves further it is on the wrong side of both endpoints; and then it is on
        // the wrong side of A only.  If the exit point B' is on the wrong side of
        // either endpoint, we can't use it; instead the segment is clipped at the
        // original endpoint B.
        //
        // We reject the segment if the sum of the scores of the two endpoints is 3
        // or more.  Here is what that rule encodes:
        //  - If B' is on the wrong side of A, then the other clipped endpoint A'
        //    must be in the interior of AB (otherwise AB' would go the wrong way
        //    around the circle).  There is a similar rule for A'.
        //  - If B' is on the wrong side of either endpoint (and therefore we must
        //    use the original endpoint B instead), then it must be possible to
        //    project B onto this face (i.e., its w-coordinate must be positive).
        //    This rule is only necessary to handle certain zero-length edges (A=B).
        var score = 0;
        if ((p - a).dotProd(aTangent) < 0) {
            score = 2;  // B' is on wrong side of A.
        } else if ((p - b).dotProd(bTangent) < 0) {
            score = 1;  // B' is on wrong side of B.
        }
        if (score > 0) {  // B' is not in the interior of AB.
            if (b[2] <= 0) {
                score = 3;    // B cannot be projected onto this face.
            } else {
                uv = R2Point(b[0] / b[2], b[1] / b[2]);
            }
        }
        return score to uv
    }

    private fun updateEndpoint(bound: R1Interval, end: Int, value: Double): Boolean {
        if (end == 0) {
            if (bound.hi < value) return false
            if (bound.lo < value) bound.lo = value
        } else {
            if (bound.lo > value) return false
            if (bound.hi > value) bound.hi = value
        }
        return true
    }

    // Given a line segment from (a0,a1) to (b0,b1) and a bounding interval for
    // each axis, clip the segment further if necessary so that "bound0" does not
    // extend outside the given interval "clip".  "diag" is a a precomputed helper
    // variable that indicates which diagonal of the bounding box is spanned by AB:
    // it is 0 if AB has positive slope, and 1 if AB has negative slope.
    private fun clipBoundAxis(
        a0: Double, b0: Double, bound0: R1Interval,
        a1: Double, b1: Double, bound1: R1Interval,
        diag: Int, clip0: R1Interval
    ): Boolean {
        if (bound0.lo < clip0.lo) {
            if (bound0.hi < clip0.lo) return false
            bound0[0] = clip0.lo
            if (!updateEndpoint(bound1, diag, interpolateDouble(clip0.lo, a0, b0, a1, b1))) return false
        }
        if (bound0.hi > clip0.hi) {
            if (bound0.lo > clip0.hi) return false
            bound0[1] = clip0.hi
            if (!updateEndpoint(bound1, 1 - diag, interpolateDouble(clip0.hi, a0, b0, a1, b1))) return false
        }
        return true
    }
}




