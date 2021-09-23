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
package dilivia.s2.region

import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkLE
import dilivia.PreConditions.checkState
import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S2Centroids
import dilivia.s2.S2Measures
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min

// Author: ericv@google.com (Eric Veach)
//
// Defines various angle and area measures for loops on the sphere.  These are
// low-level methods that work directly with arrays of S2Points.  They are
// used to implement the methods in s2shapeindex_measures.h,
// s2shape_measures.h, s2loop.h, and s2polygon.h.
//
// See s2polyline_measures.h, s2edge_distances.h, and s2measures.h for
// additional low-level methods.

object S2LoopMeasures {

    private val logger = KotlinLogging.logger { }

    // Returns the perimeter of the loop.
    fun getPerimeter(loop: S2PointLoopSpan): S1Angle {
        val perimeter = S1Angle.zero()
        if (loop.size <= 1) return perimeter;
        for (i in 0 until loop.size) {
            perimeter.radians = perimeter.radians + S1Angle(loop[i], loop[i + 1]).radians
        }
        return perimeter;
    }

    // Returns the area of the loop interior, i.e. the region on the left side of
    // the loop.  The result is between 0 and 4*Pi steradians.  The implementation
    // ensures that nearly-degenerate clockwise loops have areas close to zero,
    // while nearly-degenerate counter-clockwise loops have areas close to 4*Pi.
    fun getArea(loop: S2PointLoopSpan): Double {
        var area = getSignedArea(loop)
        checkLE(abs(area), 2 * M_PI)
        if (area < 0.0) area += 4 * M_PI
        return area;
    }

    // Like GetArea(), except that this method is faster and has more error.  The
    // result is between 0 and 4*Pi steradians.  The maximum error is 2.22e-15
    // steradians per loop vertex, which works out to about 0.09 square meters per
    // vertex on the Earth's surface.  For example, a loop with 100 vertices has a
    // maximum error of about 9 square meters.  (The actual error is typically
    // much smaller than this.)  The error bound can be computed using
    // GetCurvatureMaxError(), which returns the maximum error in steradians.
    fun getApproxArea(loop: S2PointLoopSpan): Double = 2 * M_PI - getCurvature(loop)

    // Returns either the positive area of the region on the left side of the
    // loop, or the negative area of the region on the right side of the loop,
    // whichever is smaller in magnitude.  The result is between -2*Pi and 2*Pi
    // steradians.  This method is used to accurately compute the area of polygons
    // consisting of multiple loops.
    //
    // The following cases are handled specially:
    //
    //  - Counter-clockwise loops are guaranteed to have positive area, and
    //    clockwise loops are guaranteed to have negative area.
    //
    //  - Degenerate loops (consisting of an isolated vertex or composed entirely
    //    of sibling edge pairs) have an area of exactly zero.
    //
    //  - The full loop (containing all points, and represented as a loop with no
    //    vertices) has a negative area with the minimum possible magnitude.
    //    (This is the "signed equivalent" of having an area of 4*Pi.)
    fun getSignedArea(loop: S2PointLoopSpan): Double {
        // It is suprisingly difficult to compute the area of a loop robustly.  The
        // main issues are (1) whether degenerate loops are considered to be CCW or
        // not (i.e., whether their area is close to 0 or 4*Pi), and (2) computing
        // the areas of small loops with good relative accuracy.
        //
        // With respect to degeneracies, we would like GetArea() to be consistent
        // with S2Loop::Contains(S2Point) in that loops that contain many points
        // should have large areas, and loops that contain few points should have
        // small areas.  For example, if a degenerate triangle is considered CCW
        // according to s2pred::Sign(), then it will contain very few points and
        // its area should be approximately zero.  On the other hand if it is
        // considered clockwise, then it will contain virtually all points and so
        // its area should be approximately 4*Pi.
        //
        // More precisely, let U be the set of S2Points for which S2::IsUnitLength()
        // is true, let P(U) be the projection of those points onto the mathematical
        // unit sphere, and let V(P(U)) be the Voronoi diagram of the projected
        // points.  Then for every loop x, we would like GetArea() to approximately
        // equal the sum of the areas of the Voronoi regions of the points p for
        // which x.Contains(p) is true.
        //
        // The second issue is that we want to compute the area of small loops
        // accurately.  This requires having good relative precision rather than
        // good absolute precision.  For example, if the area of a loop is 1e-12 and
        // the error is 1e-15, then the area only has 3 digits of accuracy.  (For
        // reference, 1e-12 is about 40 square meters on the surface of the earth.)
        // We would like to have good relative accuracy even for small loops.
        //
        // To achieve these goals, we combine two different methods of computing the
        // area.  This first method is based on the Gauss-Bonnet theorem, which says
        // that the area enclosed by the loop equals 2*Pi minus the total geodesic
        // curvature of the loop (i.e., the sum of the "turning angles" at all the
        // loop vertices).  The big advantage of this method is that as long as we
        // use s2pred::Sign() to compute the turning angle at each vertex, then
        // degeneracies are always handled correctly.  In other words, if a
        // degenerate loop is CCW according to the symbolic perturbations used by
        // s2pred::Sign(), then its turning angle will be approximately 2*Pi.
        //
        // The disadvantage of the Gauss-Bonnet method is that its absolute error is
        // about 2e-15 times the number of vertices (see GetCurvatureMaxError).
        // So, it cannot compute the area of small loops accurately.
        //
        // The second method is based on splitting the loop into triangles and
        // summing the area of each triangle.  To avoid the difficulty and expense
        // of decomposing the loop into a union of non-overlapping triangles,
        // instead we compute a signed sum over triangles that may overlap (see the
        // comments for S2Loop::GetSurfaceIntegral).  The advantage of this method
        // is that the area of each triangle can be computed with much better
        // relative accuracy (using l'Huilier's theorem).  The disadvantage is that
        // the result is a signed area: CCW loops may yield a small positive value,
        // while CW loops may yield a small negative value (which is converted to a
        // positive area by adding 4*Pi).  This means that small errors in computing
        // the signed area may translate into a very large error in the result (if
        // the sign of the sum is incorrect).
        //
        // So, our strategy is to combine these two methods as follows.  First we
        // compute the area using the "signed sum over triangles" approach (since it
        // is generally more accurate).  We also estimate the maximum error in this
        // result.  If the signed area is too close to zero (i.e., zero is within
        // the error bounds), then we double-check the sign of the result using the
        // Gauss-Bonnet method.  If the two methods disagree, we return the smallest
        // possible positive or negative area based on the result of GetCurvature().
        // Otherwise we return the area that we computed originally.

        // The signed area should be between approximately -4*Pi and 4*Pi.
        // Normalize it to be in the range [-2*Pi, 2*Pi].
        var area = getSurfaceIntegral(loop, SummableDouble(0.0)) { a: S2Point, b: S2Point, c: S2Point ->
            SummableDouble(S2Measures.signedArea(a, b, c))
        }
        val max_error = getCurvatureMaxError(loop)
        checkLE(abs(area), 4 * M_PI + max_error)
        area =  IEEEremainder(area, 4 * M_PI)

        // If the area is a small negative or positive number, verify that the sign
        // of the result is consistent with the loop orientation.
        if (abs(area) <= max_error) {
            val curvature = getCurvature(loop)
            // Zero-area loops should have a curvature of approximately +/- 2*Pi.
            checkState { !(area == 0.0 && curvature == 0.0) }
            if (curvature == 2 * M_PI) return 0.0;  // Degenerate
            if (area <= 0 && curvature > 0) {
                return Double.MIN_VALUE
            }
            // Full loops are handled by the case below.
            if (area >= 0 && curvature < 0) {
                return -Double.MIN_VALUE
            }
        }
        return area;
    }

    // Returns a new loop obtained by removing all degeneracies from "loop".  In
    // particular, the result will not contain any adjacent duplicate vertices or
    // sibling edge pairs, i.e. vertex sequences of the form (A, A) or (A, B, A).
    //
    // "new_vertices" represents storage where new loop vertices may be written.
    // Note that the S2PointLoopSpan result may be a subsequence of either "loop"
    // or "new_vertices", and therefore "new_vertices" must persist until the
    // result of this method is no longer needed.
    internal fun pruneDegeneracies(loop: S2PointLoopSpan): S2PointLoopSpan {
        logger.trace { "Prune Degeneracies of $loop" }
        val vertices = mutableListOf<S2Point>()
        logger.trace { "Remove edge pairs of the form ABA" }
        loop.forEachIndexed { index, v ->
            logger.trace { "Before iteration $index: $vertices" }
            // Remove duplicate vertices.
            if (vertices.isEmpty() || v != vertices.last()) {
                // Remove edge pairs of the form ABA.
                if (vertices.size >= 2 && v == vertices[vertices.size - 2]) {
                    logger.trace { "Vertex[$index] = $v == Vertex[${vertices.size - 2}] = ${vertices[vertices.size - 2]} => Remove last vertex" }
                    vertices.removeLast();
                } else {
                    logger.trace { "Add vertex" }
                    vertices.add(v)
                }
            }
            logger.trace { "After iteration $index: $vertices" }
        }
        // Check whether the loop was completely degenerate.
        if (vertices.size < 3) {
            logger.trace { "Loop contains less than 3 points, return empty loop" }
            return S2PointLoopSpan(emptyList())
        }

        // Otherwise some portion of the loop is guaranteed to be non-degenerate.
        // However there may still be some degenerate portions to remove.
        if (vertices[0] == vertices.last()) {
            logger.trace { "First point == last point, remove last point" }
            vertices.removeLast()
        }

        // If the loop begins with BA and ends with A, then there is an edge pair of
        // the form ABA at the end/start of the loop.  Remove all such pairs.  As
        // noted above, this is guaranteed to leave a non-degenerate loop.
        var k = 0
        while (k + 1 < vertices.size && vertices[k + 1] == vertices[vertices.size - (k + 1)]) {
            logger.trace { "k = $k, vertices[k + 1] = ${vertices[k + 1]}, vertices[vertices.size-(k + 1) = ${vertices.size - (k + 1)}] = ${vertices[vertices.size - (k + 1)]}" }
            ++k
        }
        val spanLength = vertices.size - 2 * k
        return if (spanLength < 3) S2PointLoopSpan(emptyList()) else S2PointLoopSpan(
            vertices.subList(
                k,
                k + spanLength
            )
        )
    }

    private val kMaxCurvature = 2 * M_PI - 4 * DoubleType.epsilon

    // Returns the geodesic curvature of the loop, defined as the sum of the turn
    // angles at each vertex (see S2::TurnAngle).  The result is positive if the
    // loop is counter-clockwise, negative if the loop is clockwise, and zero if
    // the loop is a great circle.  The geodesic curvature is equal to 2*Pi minus
    // the area of the loop.
    //
    // The following cases are handled specially:
    //
    //  - Degenerate loops (consisting of an isolated vertex or composed entirely
    //    of sibling edge pairs) have a curvature of 2*Pi exactly.
    //
    //  - The full loop (containing all points, and represented as a loop with no
    //    vertices) has a curvature of -2*Pi exactly.
    //
    //  - All other loops have a non-zero curvature in the range (-2*Pi, 2*Pi).
    //    For any such loop, reversing the order of the vertices is guaranteed to
    //    negate the curvature.  This property can be used to define a unique
    //    normalized orientation for every loop.
    fun getCurvature(loop: S2PointLoopSpan): Double {
        // By convention, a loop with no vertices contains all points on the sphere.
        if (loop.isEmpty()) return -2 * M_PI

        // Remove any degeneracies from the loop.
        val prunedLoop = pruneDegeneracies(loop)

        // If the entire loop was degenerate, it's turning angle is defined as 2*Pi.
        if (prunedLoop.isEmpty()) return 2 * M_PI

        // To ensure that we get the same result when the vertex order is rotated,
        // and that the result is negated when the vertex order is reversed, we need
        // to add up the individual turn angles in a consistent order.  (In general,
        // adding up a set of numbers in a different order can change the sum due to
        // rounding errors.)
        //
        // Furthermore, if we just accumulate an ordinary sum then the worst-case
        // error is quadratic in the number of vertices.  (This can happen with
        // spiral shapes, where the partial sum of the turning angles can be linear
        // in the number of vertices.)  To avoid this we use the Kahan summation
        // algorithm (http://en.wikipedia.org/wiki/Kahan_summation_algorithm).
        val order = getCanonicalLoopOrder(prunedLoop)
        var i = order.first
        val dir = order.dir
        var n = loop.size
        var sum = S2Measures.turnAngle(loop[(i + n - dir) % n], loop[i], loop[(i + dir) % n])
        var compensation = 0.0;  // Kahan summation algorithm
        while (--n > 0) {
            i += dir;
            var angle = S2Measures.turnAngle(loop[i - dir], loop[i], loop[i + dir])
            val old_sum = sum
            angle += compensation;
            sum += angle;
            compensation = (old_sum - sum) + angle;
        }
        sum += compensation;
        return max(-kMaxCurvature, min(kMaxCurvature, dir * sum))
    }

    private val kMaxErrorPerVertex = 9.73 * DoubleType.epsilon

    // Returns the maximum error in GetCurvature() for the given loop.  This value
    // is also an upper bound on the error in GetArea(), GetSignedArea(), and
    // GetApproxArea().
    fun getCurvatureMaxError(loop: S2PointLoopSpan): Double {
        // The maximum error can be bounded as follows:
        //   2.24 * DBL_EPSILON    for RobustCrossProd(b, a)
        //   2.24 * DBL_EPSILON    for RobustCrossProd(c, b)
        //   3.25 * DBL_EPSILON    for Angle()
        //   2.00 * DBL_EPSILON    for each addition in the Kahan summation
        //   ------------------
        //   9.73 * DBL_EPSILON
        //
        // TODO(ericv): This error estimate is approximate.  There are two issues:
        // (1) SignedArea needs some improvements to ensure that its error is
        // actually never higher than GirardArea, and (2) although the number of
        // triangles in the sum is typically N-2, in theory it could be as high as
        // 2*N for pathological inputs.  But in other respects this error bound is
        // very conservative since it assumes that the maximum error is achieved on
        // every triangle.
        return kMaxErrorPerVertex * loop.size;
    }

    // Returns the true centroid of the loop multiplied by the area of the loop
    // (see s2centroids.h for details on centroids).  The result is not unit
    // length, so you may want to normalize it.  Also note that in general, the
    // centroid may not be contained by the loop.
    //
    // The result is scaled by the loop area for two reasons: (1) it is cheaper to
    // compute this way, and (2) it makes it easier to compute the centroid of
    // more complicated shapes (by splitting them into disjoint regions and adding
    // their centroids).
    fun getCentroid(loop: S2PointLoopSpan): S2Point {
        // GetSurfaceIntegral() returns either the integral of position over loop
        // interior, or the negative of the integral of position over the loop
        // exterior.  But these two values are the same (!), because the integral of
        // position over the entire sphere is (0, 0, 0).
        return getSurfaceIntegral(
            loop,
            SummableS2Point(S2Point())
        ) { a: S2Point, b: S2Point, c: S2Point -> SummableS2Point(S2Centroids.trueCentroid(a, b, c)) }
    }

    // Returns true if the loop area is at most 2*Pi.  (A small amount of error is
    // allowed in order to ensure that loops representing an entire hemisphere are
    // always considered normalized.)
    //
    // Degenerate loops are handled consistently with s2pred::Sign(), i.e., if a
    // loop can be expressed as the union of degenerate or nearly-degenerate
    // counter-clockwise triangles then this method will return true.
    fun isNormalized(loop: S2PointLoopSpan): Boolean {
        // We allow some error so that hemispheres are always considered normalized.
        //
        // TODO(ericv): This is no longer required by the S2Polygon implementation,
        // so alternatively we could create the invariant that a loop is normalized
        // if and only if its complement is not normalized.
        return getCurvature(loop) >= -getCurvatureMaxError(loop);
    }

    // LoopOrder represents a cyclic ordering of the loop vertices, starting at
    // the index "first" and proceeding in direction "dir" (either +1 or -1).
    // "first" and "dir" must be chosen such that (first, ..., first + n * dir)
    // are all in the range [0, 2*n-1] as required by S2PointLoopSpan::operator[].
    data class LoopOrder(val first: Int, val dir: Int)

    private fun isOrderLess(order1: LoopOrder, order2: LoopOrder, loop: S2PointLoopSpan): Boolean {
        if (order1 == order2) return false;

        var i1 = order1.first
        var i2 = order2.first
        val dir1 = order1.dir
        val dir2 = order2.dir
        checkEQ(loop[i1], loop[i2])
        var n = loop.size
        while (--n > 0) {
            i1 += dir1;
            i2 += dir2;
            if (loop[i1] < loop[i2]) return true;
            if (loop[i1] > loop[i2]) return false;
        }
        return false;
    }

    // Returns an index "first" and a direction "dir" such that the vertex
    // sequence (first, first + dir, ..., first + (n - 1) * dir) does not change
    // when the loop vertex order is rotated or reversed.  This allows the loop
    // vertices to be traversed in a canonical order.
    fun getCanonicalLoopOrder(loop: S2PointLoopSpan): LoopOrder {
        // In order to handle loops with duplicate vertices and/or degeneracies, we
        // return the LoopOrder that minimizes the entire corresponding vertex
        // *sequence*.  For example, suppose that vertices are sorted
        // alphabetically, and consider the loop CADBAB.  The canonical loop order
        // would be (4, 1), corresponding to the vertex sequence ABCADB.  (For
        // comparison, loop order (4, -1) yields the sequence ABDACB.)
        //
        // If two or more loop orders yield identical minimal vertex sequences, then
        // it doesn't matter which one we return (since they yield the same result).

        // For efficiency, we divide the process into two steps.  First we find the
        // smallest vertex, and the set of vertex indices where that vertex occurs
        // (noting that the loop may contain duplicate vertices).  Then we consider
        // both possible directions starting from each such vertex index, and return
        // the LoopOrder corresponding to the smallest vertex sequence.
        val n = loop.size
        if (n == 0) return LoopOrder(0, 1)

        val min_indices = mutableListOf<Int>()
        min_indices.add(0)
        for (i in 1 until n) {
            if (loop[i] <= loop[min_indices[0]]) {
                if (loop[i] < loop[min_indices[0]]) min_indices.clear()
                min_indices.add(i)
            }
        }
        var min_order = LoopOrder(min_indices[0], 1)
        for (min_index in min_indices) {
            val order1 = LoopOrder(min_index, 1)
            val order2 = LoopOrder(min_index + n, -1)
            if (isOrderLess(order1, min_order, loop)) min_order = order1
            if (isOrderLess(order2, min_order, loop)) min_order = order2
        }
        return min_order;
    }

    // The maximum length of an edge for it to be considered numerically stable.
    // The exact value is fairly arbitrary since it depends on the stability of
    // the "f_tri" function.  The value below is quite conservative but could be
    // reduced further if desired.
    private val kMaxLength = M_PI - 1e-5

    // Returns the oriented surface integral of some quantity f(x) over the loop
    // interior, given a function f_tri(A,B,C) that returns the corresponding
    // integral over the spherical triangle ABC.  Here "oriented surface integral"
    // means:
    //
    // (1) f_tri(A,B,C) must be the integral of f if ABC is counterclockwise,
    //     and the integral of -f if ABC is clockwise.
    //
    // (2) The result of this function is *either* the integral of f over the
    //     loop interior, or the integral of (-f) over the loop exterior.
    //
    // Note that there are at least two common situations where property (2) above
    // is not a limitation:
    //
    //  - If the integral of f over the entire sphere is zero, then it doesn't
    //    matter which case is returned because they are always equal.
    //
    //  - If f is non-negative, then it is easy to detect when the integral over
    //    the loop exterior has been returned, and the integral over the loop
    //    interior can be obtained by adding the integral of f over the entire
    //    unit sphere (a constant) to the result.
    //
    // REQUIRES: The default constructor for T must initialize the value to zero.
    //           (This is true for built-in types such as "double".)
    fun <T, V> getSurfaceIntegral(
        loop: S2PointLoopSpan,
        zero: T,
        f_tri: (S2Point, S2Point, S2Point) -> T
    ): V where T : Summable<V> {
        // We sum "f_tri" over a collection T of oriented triangles, possibly
        // overlapping.  Let the sign of a triangle be +1 if it is CCW and -1
        // otherwise, and let the sign of a point "x" be the sum of the signs of the
        // triangles containing "x".  Then the collection of triangles T is chosen
        // such that either:
        //
        //  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
        //  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
        //
        // The triangles basically consist of a "fan" from vertex 0 to every loop
        // edge that does not include vertex 0.  These triangles will always satisfy
        // either (1) or (2).  However, what makes this a bit tricky is that
        // spherical edges become numerically unstable as their length approaches
        // 180 degrees.  Of course there is not much we can do if the loop itself
        // contains such edges, but we would like to make sure that all the triangle
        // edges under our control (i.e., the non-loop edges) are stable.  For
        // example, consider a loop around the equator consisting of four equally
        // spaced points.  This is a well-defined loop, but we cannot just split it
        // into two triangles by connecting vertex 0 to vertex 2.
        //
        // We handle this type of situation by moving the origin of the triangle fan
        // whenever we are about to create an unstable edge.  We choose a new
        // location for the origin such that all relevant edges are stable.  We also
        // create extra triangles with the appropriate orientation so that the sum
        // of the triangle signs is still correct at every point.


        // The default constructor for T must initialize the value to zero.
        // (This is true for built-in types such as "double".)
        var sum: T = zero
        if (loop.size < 3) return sum.value()

        var origin = loop[0]
        var i = 1
        while (i + 1 < loop.size) {
            // Let V_i be loop[i], let O be the current origin, and let length(A,B)
            // be the length of edge (A,B).  At the start of each loop iteration, the
            // "leading edge" of the triangle fan is (O,V_i), and we want to extend
            // the triangle fan so that the leading edge is (O,V_i+1).
            //
            // Invariants:
            //  1. length(O,V_i) < kMaxLength for all (i > 1).
            //  2. Either O == V_0, or O is approximately perpendicular to V_0.
            //  3. "sum" is the oriented integral of f over the area defined by
            //     (O, V_0, V_1, ..., V_i).
            checkState { i == 1 || origin.angle(loop[i]) < kMaxLength }
            checkState { origin == loop[0] || abs(origin.dotProd(loop[0])) < 1e-15 }

            if (loop[i + 1].angle(origin) > kMaxLength) {
                // We are about to create an unstable edge, so choose a new origin O'
                // for the triangle fan.
                val old_origin = origin
                if (origin == loop[0]) {
                    // The following point is well-separated from V_i and V_0 (and
                    // therefore V_i+1 as well).
                    origin = S2PointUtil.robustCrossProd(loop[0], loop[i]).normalize()
                } else if (loop[i].angle(loop[0]) < kMaxLength) {
                    // All edges of the triangle (O, V_0, V_i) are stable, so we can
                    // revert to using V_0 as the origin.
                    origin = loop[0]
                } else {
                    // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
                    // perpendicular.  Therefore V_0.CrossProd(O) is approximately
                    // perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
                    // this point O' as the new origin.
                    origin = loop[0].crossProd(old_origin)

                    // Advance the edge (V_0,O) to (V_0,O').
                    sum += f_tri(loop[0], old_origin, origin);
                }
                // Advance the edge (O,V_i) to (O',V_i).
                sum += f_tri(old_origin, loop[i], origin);
            }
            // Advance the edge (O,V_i) to (O,V_i+1).
            sum += f_tri(origin, loop[i], loop[i + 1]);
            ++i
        }
        // If the origin is not V_0, we need to sum one more triangle.
        if (origin != loop[0]) {
            // Advance the edge (O,V_n-1) to (O,V_0).
            sum += f_tri(origin, loop[loop.size - 1], loop[0]);
        }
        return sum.value();
    }

}


interface Summable<V> {

    fun value(): V

    operator fun plusAssign(other: Summable<V>): Unit

}

class SummableDouble(private var value: Double) : Summable<Double> {
    override fun value(): Double = value
    override fun plusAssign(other: Summable<Double>) {
        this.value += other.value()
    }
}

class SummableS2Point(private var value: S2Point) : Summable<S2Point> {
    override fun value(): S2Point = value
    override fun plusAssign(other: Summable<S2Point>) {
        this.value = this.value + other.value()
    }
}

