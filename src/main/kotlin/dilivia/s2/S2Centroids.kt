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
import dilivia.math.vectors.times
import dilivia.s2.S2PointUtil.isUnitLength
import org.apache.commons.math3.util.FastMath.sin
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * Utilities to compute "centroid" of a triangle.
 *
 * There are several notions of the "centroid" of a triangle.
 *  - First, there is the planar centroid, which is simply the centroid of the ordinary (non-spherical) triangle defined
 *    by the three vertices.
 *  - Second, there is the surface centroid, which is defined as the intersection of the three medians of the spherical
 *    triangle. It is possible to show that this point is simply the planar centroid projected to the surface of the
 *    sphere.
 *  - Finally, there is the true centroid (mass centroid), which is defined as the surface integral over the spherical
 *    triangle of (x,y,z) divided by the triangle area. This is the point that the triangle would rotate around if it
 *    was spinning in empty space.
 *
 * The best centroid for most purposes is the true centroid. Unlike the planar and surface centroids, the true centroid
 * behaves linearly as regions are added or subtracted. That is, if you split a triangle into pieces and compute the
 * average of their centroids (weighted by triangle area), the result equals the centroid of the original triangle.
 * This is not true of the other centroids.
 *
 * Also note that the surface centroid may be nowhere near the intuitive "center" of a spherical triangle.
 * For example, consider the triangle with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere). The surface
 * centroid of this triangle is at S=(0, 2*eps, 1), which is within a distance of 2*eps of the vertex B. Note that the
 * median from A (the segment connecting A to the midpoint of BC) passes through S, since this is the shortest path
 * connecting the two endpoints. On the other hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
 * the surface is a much more reasonable interpretation of the "center" of this triangle.
 *
 * This class is a port of the s2centroids methods of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2Centroids {

    /**
     * Computes the centroid of the planar triangle ABC. This can be normalized to unit length to obtain the
     * "surface centroid" of the corresponding spherical triangle, i.e. the intersection of the three medians.
     * However, note that for large spherical triangles the surface centroid may be nowhere near the intuitive "center"
     * (see example above).
     *
     * @param a Point A of the ABC triangle.
     * @param b Point B of the ABC triangle.
     * @param c Point C of the ABC triangle.
     * @return The planar centroid of the triangle ABC.
     */
    @JvmStatic
    fun planarCentroid(a: S2Point, b: S2Point, c: S2Point): S2Point = (1.0 / 3.0) * (a + b + c)

    /**
     * Computes the true centroid of the spherical triangle ABC multiplied by the signed area of spherical triangle ABC.
     * The reasons for multiplying by the signed area are
     * (1) this is the quantity that needs to be summed to compute the centroid of a union or difference of triangles,
     * (2) it's actually easier to calculate this way.
     *
     * All points must have unit length.
     *
     * Note that the result of this function is defined to be S2Point(0, 0, 0) if the triangle is degenerate (and that
     * this is intended behavior).
     *
     * @param a Point A of the ABC triangle.
     * @param b Point B of the ABC triangle.
     * @param c Point C of the ABC triangle.
     * @return The true centroid of the triangle ABC.
     */
    @JvmStatic
    fun trueCentroid(a: S2Point, b: S2Point, c: S2Point): S2Point {
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        requireArgument { isUnitLength(c) }

        // I couldn't find any references for computing the true centroid of a
        // spherical triangle...  I have a truly marvellous demonstration of this
        // formula which this margin is too narrow to contain :)

        // Use Angle() in order to get accurate results for small triangles.
        val angleA = b.angle(c)
        val angleB = c.angle(a)
        val angleC = a.angle(b)
        val ra = if (angleA == 0.0) 1.0 else (angleA / sin(angleA))
        val rb = if (angleB == 0.0) 1.0 else (angleB / sin(angleB))
        val rc = if (angleC == 0.0) 1.0 else (angleC / sin(angleC))

        // Now compute a point M such that:
        //
        //  [Ax Ay Az] [Mx]                       [ra]
        //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
        //  [Cx Cy Cz] [Mz]                       [rc]
        //
        // To improve the numerical stability we subtract the first row (A) from the
        // other two rows; this reduces the cancellation error when A, B, and C are
        // very close together.  Then we solve it using Cramer's rule.
        //
        // The result is the true centroid of the triangle multiplied by the
        // triangle's area.
        //
        // TODO(ericv): This code still isn't as numerically stable as it could be.
        // The biggest potential improvement is to compute B-A and C-A more
        // accurately so that (B-A)x(C-A) is always inside triangle ABC.
        val x = S2Point(a.x, b.x - a.x, c.x - a.x)
        val y = S2Point(a.y, b.y - a.y, c.y - a.y)
        val z = S2Point(a.z, b.z - a.z, c.z - a.z)
        val r = S2Point(ra, rb - ra, rc - ra)
        return 0.5 * S2Point(y.crossProd(z).dotProd(r), z.crossProd(x).dotProd(r), x.crossProd(y).dotProd(r))
    }

    /**
     * Conputes the true centroid of the spherical geodesic edge AB multiplied by the length of the edge AB.
     * As with triangles, the true centroid of a collection of edges may be computed simply by summing the result of
     * this method for each edge.
     *
     * Note that the planar centroid of a geodesic edge simply 0.5 * (a + b), while the surface centroid is
     * (a + b).normalize(). However neither of these values is appropriate for computing the centroid of a collection of
     * edges (such as a polyline).
     *
     * Also note that the result of this function is defined to be S2Point(0, 0, 0)
     * if the edge is degenerate (and that this is intended behavior).
     *
     * @param a Point A of the spherical geodesic edge AB.
     * @param b Point B of the spherical geodesic edge AB.
     * @return The true centroid of the spherical geodesic edge AB.
     */
    fun trueCentroid(a: S2Point, b: S2Point): S2Point {
        // The centroid (multiplied by length) is a vector toward the midpoint
        // of the edge, whose length is twice the sine of half the angle between
        // the two vertices.  Defining theta to be this angle, we have:
        val vdiff = a - b  // Length == 2*sin(theta)
        val vsum = a + b   // Length == 2*cos(theta)
        val sin2 = vdiff.norm2()
        val cos2 = vsum.norm2()
        if (cos2 == 0.0) return S2Point()  // Ignore antipodal edges.
        return sqrt(sin2 / cos2) * vsum  // Length == 2*sin(theta)
    }

}
