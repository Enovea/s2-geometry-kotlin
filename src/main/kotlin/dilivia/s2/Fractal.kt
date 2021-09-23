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

import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLT
import dilivia.math.matrix.Matrix3x3Double
import dilivia.math.vectors.R2Point
import dilivia.math.vectors.R2VectorDouble
import dilivia.math.vectors.times
import dilivia.s2.region.S2Loop
import org.apache.commons.math3.util.FastMath.log10
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.pow
import org.apache.commons.math3.util.FastMath.round
import org.apache.commons.math3.util.FastMath.sqrt
import kotlin.random.Random
import org.apache.commons.math3.util.FastMath.log as log2

// A simple class that generates "Koch snowflake" fractals (see Wikipedia
// for an introduction).  There is an option to control the fractal
// dimension (between 1.0 and 2.0); values between 1.02 and 1.50 are
// reasonable simulations of various coastlines.  The default dimension
// (about 1.26) corresponds to the standard Koch snowflake.  (The west coast
// of Britain has a fractal dimension of approximately 1.25.)
//
// The fractal is obtained by starting with an equilateral triangle and
// recursively subdividing each edge into four segments of equal length.
// Therefore the shape at level "n" consists of 3*(4**n) edges.  Multi-level
// fractals are also supported: if you set min_level() to a non-negative
// value, then the recursive subdivision has an equal probability of
// stopping at any of the levels between the given min and max (inclusive).
// This yields a fractal where the perimeter of the original triangle is
// approximately equally divided between fractals at the various possible
// levels.  If there are k distinct levels {min,..,max}, the expected number
// of edges at each level "i" is approximately 3*(4**i)/k.
class Fractal(

        // Set the maximum subdivision level for the fractal (see above).
        // REQUIRES: max_level >= 0
        private var maxLevel: Int = -1,

        // Set the minimum subdivision level for the fractal (see above).  The
        // default value of -1 causes the min and max levels to be the same.  A
        // min_level of 0 should be avoided since this creates a significant
        // chance that none of the three original edges will be subdivided at all.
        //
        // DEFAULT: max_level()
        private var minLevel: Int = maxLevel,

        // Set the fractal dimension.  The default value of approximately 1.26
        // corresponds to the stardard Koch curve.  The value must lie in the
        // range [1.0, 2.0).
        //
        // DEFAULT: log(4) / log(3) ~= 1.26
        private var dimension: Double = log10(4.0) / log10(3.0)

) {

    init {
        computeOffsets()
    }

    private var minLevelArg: Int = -1  // Value set by user

    // The ratio of the sub-edge length to the original edge length at each
    // subdivision step.
    private var edgeFraction: Double = 0.0

    // The distance from the original edge to the middle vertex at each
    // subdivision step, as a fraction of the original edge length.
    private var offsetFraction: Double = 0.0

    fun setMaxLevel(level: Int) {
        requireGE(level, 0)
        this.maxLevel = level
        computeMinLevel()
    }
    fun getMaxLevel(): Int = maxLevel

    fun setMinLevel(level: Int) {
        requireGE(level, -1)
        minLevelArg = level
        computeMinLevel()
    }
    fun getMinLevel(): Int = minLevel

    fun setDimension(dimension: Double) {
        requireGE(dimension, 1.0)
        requireLT(dimension, 2.0)
        this.dimension = dimension
        computeOffsets()
    }

    // Set the min and/or max level to produce approximately the given number
    // of edges.  (The values are rounded to a nearby value of 3*(4**n).)
    fun setLevelForApproxMinEdges(minEdges: Int) {
        // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
        setMinLevel(round(0.5 * log2(minEdges / 3.0)).toInt())
    }

    fun setLevelForApproxMaxEdges(maxEdges: Int) {
        // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
        setMaxLevel(round(0.5 * log2(maxEdges / 3.0)).toInt())
    }


    // Return a lower bound on ratio (Rmin / R), where "R" is the radius
    // passed to MakeLoop() and "Rmin" is the minimum distance from the
    // fractal boundary to its center, where all distances are measured in the
    // tangent plane at the fractal's center.  This can be used to inscribe
    // another geometric figure within the fractal without intersection.
    fun minRadiusFactor(): Double {
        // The minimum radius is attained at one of the vertices created by the
        // first subdivision step as long as the dimension is not too small (at
        // least kMinDimensionForMinRadiusAtLevel1, see below).  Otherwise we fall
        // back on the incircle radius of the original triangle, which is always a
        // lower bound (and is attained when dimension = 1).
        //
        // The value below was obtained by letting AE be an original triangle edge,
        // letting ABCDE be the corresponding polyline after one subdivision step,
        // and then letting BC be tangent to the inscribed circle at the center of
        // the fractal O.  This gives rise to a pair of similar triangles whose edge
        // length ratios can be used to solve for the corresponding "edge fraction".
        // This method is slightly conservative because it is computed using planar
        // rather than spherical geometry.  The value below is equal to
        // -log(4)/log((2 + cbrt(2) - cbrt(4))/6).
        val kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407
        if (dimension >= kMinDimensionForMinRadiusAtLevel1) {
            return sqrt(1 + 3 * edgeFraction * (edgeFraction - 1))
        }
        return 0.5
    }

    // Return the ratio (Rmax / R), where "R" is the radius passed to
    // MakeLoop() and "Rmax" is the maximum distance from the fractal boundary
    // to its center, where all distances are measured in the tangent plane at
    // the fractal's center.  This can be used to inscribe the fractal within
    // some other geometric figure without intersection.
    fun maxRadiusFactor(): Double {
        // The maximum radius is always attained at either an original triangle
        // vertex or at a middle vertex from the first subdivision step.
        return max(1.0, offsetFraction * sqrt(3.0) + 0.5)
    }

    // Return a fractal loop centered around the z-axis of the given
    // coordinate frame, with the first vertex in the direction of the
    // positive x-axis.  In order to avoid self-intersections, the fractal is
    // generated by first drawing it in a 2D tangent plane to the unit sphere
    // (touching at the fractal's center point) and then projecting the edges
    // onto the sphere.  This has the side effect of shrinking the fractal
    // slightly compared to its nominal radius.
    fun makeLoop(frame: Matrix3x3Double, nominalRadius: S1Angle): S2Loop {
        val r2vertices = mutableListOf<R2Point>()
        getR2Vertices(r2vertices)
        val vertices = mutableListOf<S2Point>()
        val r = nominalRadius.radians
        for (v in r2vertices) {
            val p = S2Point(v[0] * r, v[1] * r, 1.0)
            vertices.add(S2PointUtil.fromFrame(frame, p).normalize())
        }
        return S2Loop(vertices)
    }

    private fun computeMinLevel() {
        if (minLevelArg in 0..maxLevel) {
            minLevel = minLevelArg
        } else {
            minLevel = maxLevel
        }
    }

    private fun computeOffsets() {
        edgeFraction = pow(4.0, -1.0 / dimension)
        offsetFraction = sqrt(edgeFraction - 0.25)
    }

    private fun getR2Vertices(vertices: MutableList<R2Point>) {
        // The Koch "snowflake" consists of three Koch curves whose initial edges
        // form an equilateral triangle.
        val v0 = R2VectorDouble(1.0, 0.0)
        val v1 = R2VectorDouble(-0.5, sqrt(3.0)/2)
        val v2 = R2VectorDouble(-0.5, -sqrt(3.0)/2)
        getR2VerticesHelper(v0, v1, 0, vertices)
        getR2VerticesHelper(v1, v2, 0, vertices)
        getR2VerticesHelper(v2, v0, 0, vertices)
    }

    private fun getR2VerticesHelper(v0: R2Point, v4: R2Point, level: Int, vertices: MutableList<R2Point>) {
        if (level >= minLevel && Random.Default.nextInt(maxLevel - level + 1) == 0) {
            // Stop subdivision at this level.
            vertices.add(v0)
            return
        }
        // Otherwise compute the intermediate vertices v1, v2, and v3.
        val dir = v4 - v0
        val v1 = v0 + edgeFraction * dir
        val v2 = 0.5 * (v0 + v4) - offsetFraction * dir.ortho()
        val v3 = v4 - edgeFraction * dir

        // And recurse on the four sub-edges.
        getR2VerticesHelper(v0, v1, level+1, vertices)
        getR2VerticesHelper(v1, v2, level+1, vertices)
        getR2VerticesHelper(v2, v3, level+1, vertices)
        getR2VerticesHelper(v3, v4, level+1, vertices)
    }

    override fun toString(): String {
        return "Fractal(maxLevel=$maxLevel, minLevel=$minLevel, dimension=$dimension, minLevelArg=$minLevelArg, edgeFraction=$edgeFraction, offsetFraction=$offsetFraction)"
    }


}
