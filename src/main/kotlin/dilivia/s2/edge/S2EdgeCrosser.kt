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

import dilivia.PreConditions.requireArgument
import dilivia.math.DoubleType
import dilivia.math.M_SQRT3
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil.isUnitLength
import dilivia.s2.S2PointUtil.robustCrossProd
import dilivia.s2.S2Predicates
import mu.KotlinLogging

/**
 * This class allows edges to be efficiently tested for intersection with a given fixed edge AB.
 * It is especially efficient when testing for intersection with an edge chain connecting vertices v0, v1, v2, ...
 *
 * Example usages:
 *
 *   fun countIntersections(a: S2Point, b: S2Point, edges: List<Pair<S2Point, S2Point>>) {
 *     val count: Int = 0
 *     val crosser = S2EdgeCrosser(a, b)
 *     for (edge : edges) {
 *       if (crosser.crossingSign(edge.first, edge.second) >= 0) {
 *         ++count
 *       }
 *     }
 *     return count
 *   }
 *
 * This class, store last edge of tested point in order to compute tangents twince when processing a chain. So, it is
 * also possible to do:
 *
 *   val chain: List<S2Point>
 *   crosser.restartAt(chain[0])
 *   for (i in 1..chain.lastIndex) {
 *     if (crosser.edgeOrVertexCrossing(chain[i])) { ++count }
 *   }
 *
 * This class is a port of the S2EdgeCrosser class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @constructor Default constructor. must be followed by a call to init().
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2EdgeCrosser() {

    /** The point A of the edge AB to test. */
    private lateinit var a: S2Point

    /** The point B of the edge AB to test. */
    private lateinit var b: S2Point

    /** The cross product a ^ b. */
    private lateinit var aCrossB: S2Point

    // To reduce the number of calls to s2pred::ExpensiveSign(), we compute an
    // outward-facing tangent at A and B if necessary.  If the plane
    // perpendicular to one of these tangents separates AB from CD (i.e., one
    // edge on each side) then there is no intersection.

    /** Indicates of the rangents have been computed. */
    private var tangentsComputed: Boolean = false

    /** Outward-facing tangent at A. */
    private lateinit var aTangent: S2Point

    /** Outward-facing tangent at B. */
    private lateinit var bTangent: S2Point

    // The fields below are updated for each vertex in the chain.
    /** Previous vertex in the vertex chain. */
    private var c: S2Point? = null

    /** The orientation of triangle ACB. */
    private var acb: Int = 0

    // The field below is a temporary used by crossingSignInternal().
    /** The orientation of triangle BDA. */
    private var bda: Int = 0

    /**
     * Convenience constructor that calls init() with the given fixed edge AB.
     *
     * @param a The first point of the edge AB.
     * @param b The second point of the edge AB.
     */
    constructor(a: S2Point, b: S2Point) : this() {
        init(a, b)
    }

    /**
     * Convenience constructor that uses AB as the fixed edge, and C as the first vertex of the vertex chain
     * (equivalent to calling restartAt(c)).
     *
     * @param a The first point of the edge AB.
     * @param b The second point of the edge AB.
     * @param c The first point of the vertex chain to test.
     */
    constructor(a: S2Point, b: S2Point, c: S2Point) : this() {
        init(a, b)
        restartAt(c)
    }

    /**
     * Initialize the S2EdgeCrosser with the given fixed edge AB.
     *
     * @param a The first point of the edge AB.
     * @param b The second point of the edge AB.
     * @param aCrossB A pre-computed a^b.
     */
    fun init(a: S2Point, b: S2Point, aCrossB: S2Point = a.crossProd(b)) {
        requireArgument { isUnitLength(a) }
        requireArgument { isUnitLength(b) }
        this.a = a
        this.b = b
        this.aCrossB = aCrossB
        tangentsComputed = false
        c = null
    }

    /**
     * Get the current edge A point.
     *
     * @return The A point.
     */
    fun a(): S2Point = a

    /**
     * Get the current edge B point.
     *
     * @return The B point.
     */
    fun b(): S2Point = b

    /**
     * Get the current C point.
     *
     * @return The C point.
     */
    fun c(): S2Point? = c

    ///////////////////////// Single Edge Methods ///////////////////////////

    /**
     * This function determines whether the edge AB intersects the edge CD.
     *
     * Note that if an edge is degenerate (A == B or C == D), the return value is 0, if two vertices from different
     * edges are the same and -1 otherwise.
     *
     * Properties of crossingSign:
     *
     * (1) crossingSign(b,a,c,d) == crossingSign(a,b,c,d)
     * (2) crossingSign(c,d,a,b) == crossingSign(a,b,c,d)
     * (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
     * (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
     *
     * This function implements an exact, consistent perturbation model such that no three points are ever considered to
     * be collinear. This means that even if you have 4 points A, B, C, D that lie exactly in a line
     * (say, around the equator), C and D will be treated as being slightly to one side or the other of AB.
     * This is done in a way such that the results are always consistent (see S2Predicates.sign).
     *
     * Note that if you want to check an edge against a chain of other edges, it is slightly more efficient to use the
     * single-argument version of crossingSign below.
     *
     * The arguments must point to values that persist until the next call.
     *
     * @param c The first point of the edge CD
     * @param d The second point of the edge CD
     *
     * @return +1 if AB crosses CD at a point that is interior to both edges, 0 if any two vertices from different edges
     * are the same and -1 otherwise.
     */
    fun crossingSign(c: S2Point, d: S2Point): Int {
        if (this.c != c) restartAt(c)
        return crossingSign(d)
    }

    /**
     * This method extends the concept of a "crossing" to the case where AB and CD have a vertex in common. The two
     * edges may or may not cross, according to the rules defined in vertexCrossing() below. The rules are designed so
     * that point containment tests can be implemented simply by counting edge crossings. Similarly, determining whether
     * one edge chain crosses another edge chain can be implemented by counting.
     *
     * @param c The first point of the edge CD
     * @param d The second point of the edge CD
     *
     * @return true if CrossingSign(c, d) > 0, or AB and CD share a vertex and vertexCrossing(a, b, c, d) returns true.
     */
    fun edgeOrVertexCrossing(c: S2Point, d: S2Point): Boolean {
        if (this.c != c) restartAt(c)
        return edgeOrVertexCrossing(d)
    }

    ///////////////////////// Edge Chain Methods ///////////////////////////

    /**
     * Call this method when your chain 'jumps' to a new place.
     *
     * @param c The first point of the next edge to test.
     */
    fun restartAt(c: S2Point) {
        requireArgument { isUnitLength(c) }
        this.c = c
        this.acb = -S2Predicates.triageSign(a, b, c, aCrossB)
    }

    /**
     * Like crossingSign above, but uses the last vertex passed to one of the crossing methods (or restartAt) as the
     * first vertex of the current edge.
     *
     * @param d The second point of the edge CD
     *
     * @return +1 if AB crosses CD at a point that is interior to both edges, 0 if any two vertices from different edges
     * are the same and -1 otherwise.
     */
    fun crossingSign(d: S2Point): Int {
        requireArgument({ isUnitLength(d) }, { "Point $d is not unit length: norm = ${d.norm()}" })
        // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
        // all be oriented the same way (CW or CCW).  We keep the orientation of ACB
        // as part of our state.  When each new point D arrives, we compute the
        // orientation of BDA and check whether it matches ACB.  This checks whether
        // the points C and D are on opposite sides of the great circle through AB.

        // Recall that TriageSign is invariant with respect to rotating its
        // arguments, i.e. ABD has the same orientation as BDA.
        val bda = S2Predicates.triageSign(a, b, d, aCrossB)
        if (acb == -bda && bda != 0) {
            // The most common case -- triangles have opposite orientations.  Save the
            // current vertex D as the next vertex C, and also save the orientation of
            // the new triangle ACB (which is opposite to the current triangle BDA).
            this.c = d
            this.acb = -bda
            return -1;
        }
        this.bda = bda
        return crossingSignInternal(d)
    }

    /**
     * Like EdgeOrVertexCrossing above, but uses the last vertex passed to one of the crossing methods (or restartAt)
     * as the first vertex of the current edge.
     *
     * @param d The second point of the edge CD
     *
     * @return true if crossingSign(c, d) > 0, or AB and CD share a vertex and vertexCrossing(a, b, c, d) returns true.
     */
    fun edgeOrVertexCrossing(d: S2Point): Boolean {
        // We need to copy c_ since it is clobbered by CrossingSign().
        val c = this.c!!
        val crossing = crossingSign(d)
        logger.trace { "edgeOrVertexCrossing(d = $d): crossing sign = $crossing" }
        if (crossing < 0) return false;
        if (crossing > 0) return true;
        val vertexCrossing = S2EdgeCrossings.vertexCrossing(a, b, c, d)
        logger.trace { "edgeOrVertexCrossing(d = $d): vertex crossing = $vertexCrossing" }
        return vertexCrossing
    }

    ///////////////////////// Internal Methods ///////////////////////////

    // These functions handle the "slow path" of CrossingSign().
    private fun crossingSignInternal(d: S2Point): Int {
        // Compute the actual result, and then save the current vertex D as the next
        // vertex C, and save the orientation of the next triangle ACB (which is
        // opposite to the current triangle BDA).
        val result = crossingSignInternal2(d)
        this.c = d
        this.acb = -this.bda
        return result;
    }

    private fun crossingSignInternal2(d: S2Point): Int {
        val c = this.c!!

        // At this point, a very common situation is that A,B,C,D are four points on
        // a line such that AB does not overlap CD.  (For example, this happens when
        // a line or curve is sampled finely, or when geometry is constructed by
        // computing the union of S2CellIds.)  Most of the time, we can determine
        // that AB and CD do not intersect by computing the two outward-facing
        // tangents at A and B (parallel to AB) and testing whether AB and CD are on
        // opposite sides of the plane perpendicular to one of these tangents.  This
        // is moderately expensive but still much cheaper than s2pred::ExpensiveSign.
        if (!tangentsComputed) {
            val norm = robustCrossProd(a, b).normalize()
            aTangent = a.crossProd(norm)
            bTangent = norm.crossProd(b)
            tangentsComputed = true;
        }
        // The error in RobustCrossProd() is insignificant.  The maximum error in
        // the call to CrossProd() (i.e., the maximum norm of the error vector) is
        // (0.5 + 1/sqrt(3)) * DBL_EPSILON.  The maximum error in each call to
        // DotProd() below is DBL_EPSILON.  (There is also a small relative error
        // term that is insignificant because we are comparing the result against a
        // constant that is very close to zero.)
        val kError = (1.5 + 1 / M_SQRT3) * DoubleType.epsilon
        if ((c.dotProd(aTangent) > kError && d.dotProd(aTangent) > kError)
            || (c.dotProd(bTangent) > kError && d.dotProd(bTangent) > kError)
        ) {
            return -1
        }

        // Otherwise, eliminate the cases where two vertices from different edges
        // are equal.  (These cases could be handled in the code below, but we would
        // rather avoid calling ExpensiveSign whenever possible.)
        if (a == c || a == d || b == c || b == d) return 0;

        // Eliminate cases where an input edge is degenerate.  (Note that in most
        // cases, if CD is degenerate then this method is not even called because
        // acb_ and bda have different signs.)
        if (a == b || c == d) return -1;

        // Otherwise it's time to break out the big guns.
        if (acb == 0) acb = -S2Predicates.expensiveSign(a, b, c)
        assert(acb != 0)
        if (bda == 0) bda = S2Predicates.expensiveSign(a, b, d)
        assert(bda != 0)
        if (bda != acb) return -1;

        val cCrossD = c.crossProd(d)
        val cbd = -S2Predicates.sign(c, d, b, cCrossD)
        assert(cbd != 0)
        if (cbd != acb) return -1;
        val dac = S2Predicates.sign(c, d, a, cCrossD)
        assert(dac != 0);
        return if (dac != acb) -1 else 1
    }

    companion object {
        val logger = KotlinLogging.logger(S2EdgeCrosser::class.java.name)
    }
}
