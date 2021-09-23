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

import dilivia.PreConditions
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLT
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.R1Interval
import dilivia.math.matrix.Matrix3x3Double
import dilivia.s2.S1Angle
import dilivia.s2.S1Interval
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Predicates
import dilivia.s2.S2WedgeRelations
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2XYZFaceSiTi
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.edge.S2EdgeClipping.kFaceClipErrorUVCoord
import dilivia.s2.edge.S2EdgeClipping.kIntersectsRectErrorUVDist
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.shape.AbstractMutableS2ShapeIndex
import dilivia.s2.index.shape.CellRelation
import dilivia.s2.index.shape.CellRelation.INDEXED
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.index.shape.S2CrossingEdgePairsScanner
import dilivia.s2.index.shape.S2ShapeIndex.RangeIterator
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.TypeTag
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.sin
import java.util.concurrent.atomic.AtomicInteger

/**
 * An S2Loop represents a simple spherical polygon. It consists of a single chain of vertices where the first vertex is
 * implicitly connected to the last. All loops are defined to have a CCW orientation, i.e. the interior of the loop is
 * on the left side of the edges. This implies that a clockwise loop enclosing a small area is interpreted to be a CCW
 * loop enclosing a very large area.
 *
 * Loops are not allowed to have any duplicate vertices (whether adjacent or not). Non-adjacent edges are not allowed
 * to intersect, and furthermore edges of length 180 degrees are not allowed (i.e., adjacent vertices cannot be
 * antipodal).  Loops must have at least 3 vertices (except for the empty and full loops discussed below).
 * Although these restrictions are not enforced in optimized code, you may get unexpected results if they are violated.
 *
 * There are two special loops: the "empty loop" contains no points, while the "full loop" contains all points.
 * These loops do not have any edges, but to preserve the invariant that every loop can be represented as a vertex
 * chain, they are defined as having exactly one vertex each (see kEmpty and kFull).
 *
 * Point containment of loops is defined such that if the sphere is subdivided into faces (loops), every point is
 * contained by exactly one face. This implies that loops do not necessarily contain their vertices.
 *
 * Note: The reason that duplicate vertices and intersecting edges are not allowed is that they make it harder to
 * define and implement loop relationships, e.g. whether one loop contains another.  If your data does not satisfy
 * these restrictions, you can use S2Builder to normalize it.
 *
 * This class is a port of the S2Loop class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2Loop : S2Region {

    /////////////////////////// Fields ///////////////////////////

    /**
     * Allows overriding the automatic validity checks controlled by the -ea flag. If this flag is true, then loops are
     * automatically checked for validity as they are initialized. The main reason to disable this flag is if you
     * intend to call isValid() explicitly, like this:
     *
     *   val loop = S2Loop()
     *   loop.debugOverride = S2Debug.DISABLE
     *   loop.init(...)
     *   if (!loop.isValid()) { ... }
     *
     * Without the call to loop.debugOverride = S2Debug.DISABLE, invalid data would cause a fatal error in init()
     * whenever the -ea flag is enabled.
     */
    var debugOverride: S2Debug = S2Debug.ALLOW

    /**
     * The nesting depth, if this field belongs to an S2Polygon.
     *
     * The depth of a loop is defined as its nesting level within its containing polygon. "Outer shell" loops have
     * depth 0, holes within those loops have depth 1, shells within those holes have depth 2, etc.  This field is only
     * used by the S2Polygon implementation.
     */
    var depth = 0

    /** Vertices of the polyline. */
    private val vertices: ArrayList<S2Point> = ArrayList()
    fun vertices(): List<S2Point> = vertices

    /** The number of vertices of the loop. */
    val numVertices: Int
        get() = vertices.size

    /** Indicates if the loop contains S2PointUtil.origin(). */
    private var originInside = false

    /**
     * In general we build the index the first time it is needed, but we make an exception for contains(S2Point)
     * because this method has a simple brute force implementation that is also relatively cheap. For this one method
     * we keep track of the number of calls made and only build the index once enough calls have been made that we
     * think an index would be worthwhile.
     */
    private val unindexedContainsCalls: AtomicInteger = AtomicInteger(0)

    /**
     * "bound" is a conservative bound on all points contained by this loop:
     *      if A.contains(P), then A.bound.contains(S2LatLng(P)).
     */
    val bound: S2LatLngRect = S2LatLngRect.empty()

    /**
     * Since "bound" is not exact, it is possible that a loop A contains another loop B whose bounds are slightly
     * larger.  "subregionBound" has been expanded sufficiently to account for this error, i.e.
     * if A.contains(B), then A.subregionBound.contains(B.bound).
     */
    val subregionBound: S2LatLngRect = S2LatLngRect.empty()

    // Spatial index for this loop.
    val index: MutableS2ShapeIndex = MutableS2ShapeIndex()

    /////////////////////////// Init /////////////////////////////

    /**
     * Convenience constructor to disable the automatic validity checking.
     *
     * The main reason to use this constructor is if you intend to call isValid() explicitly.
     */
    constructor(vertices: List<S2Point>, depth: Int = 0, debugOverride: S2Debug = S2Debug.ALLOW) {
        this.depth = depth
        this.debugOverride = debugOverride
        init(vertices)
    }

    /**
     * Default constructor. The loop must be initialized by calling init() before it is used.
     */
    constructor() : this(emptyList())

    /**
     * Construct a loop corresponding to the given cell.
     *
     * Note that the loop and cell *do not* contain exactly the same set of points, because S2Loop and S2Cell have
     * slightly different definitions of point containment. For example, an S2Cell vertex is contained by all four
     * neighboring S2Cells, but it is contained by exactly one of four S2Loops constructed from those cells. As another
     * example, the S2Cell coverings of "cell" and "S2Loop(cell)" will be different, because the loop contains points
     * on its boundary that actually belong to other cells (i.e., the covering will include a layer of neighboring
     * cells).
     *
     * @param cell A s2 cell.
     */
    constructor(cell: S2Cell) {
        logger.trace { "Create loop from cell: ${cell.id()}" }
        depth = 0
        vertices.ensureCapacity(4)
        debugOverride = S2Debug.ALLOW
        unindexedContainsCalls.set(0)
        for (i in 0..3) {
            vertices.add(cell.getVertex(i))
        }

        logger.trace { "Loop vertices: $vertices" }
        // We recompute the bounding rectangle ourselves, since S2Cell uses a
        // different method and we need all the bounds to be consistent.
        initOriginAndBound()
    }

    // Internal copy constructor used only by clone() that makes a deep copy of
    // its argument.
    private constructor(src: S2Loop) {
        this.depth = src.depth
        this.vertices.clear()
        this.vertices.addAll(src.vertices)
        this.debugOverride = src.debugOverride
        this.originInside = src.originInside
        this.unindexedContainsCalls.set(0)
        this.bound.init(src.bound)
        this.subregionBound.init(src.subregionBound)
        this.initIndex()
    }

    /**
     * Initialize a loop with given vertices. The last vertex is implicitly connected to the first. All points should
     * be unit length. Loops must have at least 3 vertices (except for the empty and full loops, see kEmpty and kFull).
     * This method may be called multiple times.
     *
     * @param vertices The vertices of the loop.
     */
    fun init(vertices: List<S2Point>) {
        clearIndex()
        this.vertices.clear()
        this.vertices.addAll(vertices)
        initOriginAndBound()
    }

    /**
     * Indicates if the loop is valid.
     *
     * Note that validity is checked automatically during initialization when -ea is enabled.
     * @return true if this is a valid loop.
     */
    fun isValid(): Boolean {
        val error = S2Error()
        var valid = true
        if (findValidationError(error)) {
            logger.error { error.text }
            valid = false
        }
        return valid
    }

    /**
     * Validates the loop.
     *
     * @return true if this is *not* a valid loop and sets "error" appropriately. Otherwise returns false and leaves
     * "error" unchanged.
     */
    fun findValidationError(error: S2Error): Boolean =
        (findValidationErrorNoIndex(error) || S2CrossingEdgePairsScanner.findSelfIntersection(index, error))


    fun toDebugString(): String {
        if (this.isEmpty()) {
            return "empty"
        } else if (this.isFull()) {
            return "full"
        }
        return this.vertices().joinToString(", ", prefix = "depth=$depth,") { p -> S2PointUtil.toDegreesString(p) }
    }

    fun containsOrigin(): Boolean = originInside

    /**
     * Get the vertex i of this loop.
     *
     * For convenience, we make two entire copies of the vertex list available: vertex(n..2*n-1) is mapped to
     * vertex(0..n-1), where n == numVertices.
     *
     * REQUIRES: 0 <= i < 2 * numVertices
     *
     * @param i The vertex index.
     * @return the vertex at index i.
     */
    fun vertex(i: Int): S2Point {
        requireGE(i, 0)
        requireLT(i, 2 * numVertices)
        val j = i - numVertices
        return vertices[if (j < 0) i else j]
    }

    /*
     * Like vertex(), but this method returns vertices in reverse order if the loop represents a polygon hole.
     * For example, arguments 0, 1, 2 are mapped to vertices n-1, n-2, n-3, where n == numVertices. This ensures that
     * the interior of the polygon is always to the left of the vertex chain.
     *
     * REQUIRES: 0 <= i < 2 * num_vertices()
     *
     * @param i The vertex index.
     * @return the vertex at index i.
     */
    fun orientedVertex(i: Int): S2Point {
        requireGE(i, 0)
        requireLT(i, 2 * numVertices)
        var j = i - numVertices
        if (j < 0) j = i
        if (isHole()) j = numVertices - 1 - j
        return vertices[j]
    }

    /**
     * Indicates if this loop is empty.
     *
     * @return true if this is the special empty loop that contains no points.
     */
    fun isEmpty(): Boolean = isEmptyOrFull() && !originInside

    /**
     * Indicates if this loop is full.
     *
     * @return true if this is the special full loop that contains all points.
     */
    fun isFull(): Boolean = isEmptyOrFull() && originInside

    /**
     *  Indicates if this loop is empty or full.
     *
     * @return true if this loop is either empty or full.
     */
    fun isEmptyOrFull(): Boolean = numVertices == 1

    /**
     * Indicates is this loop is a hole.
     * @return true if this loop represents a hole in its containing polygon.
     */
    fun isHole(): Boolean = (depth and 1) != 0

    /**
     * The sign of a loop is -1 if the loop represents a hole in its containing polygon, and +1 otherwise.
     * @return The sign of the loop.
     */
    fun sign(): Int = if (isHole()) -1 else 1

    /**
     * Indicates if the loop is normalized.
     *
     * Degenerate loops are handled consistently with S2Predicates.sign(), i.e., if a loop can be expressed as the
     * union of degenerate or nearly-degenerate CCW triangles, then it will always be considered normalized.
     * @return true if the loop area is at most 2*Pi.
     */
    fun isNormalized(): Boolean {
        // Optimization: if the longitude span is less than 180 degrees, then the
        // loop covers less than half the sphere and is therefore normalized.
        if (bound.lng.length < M_PI) return true

        return S2LoopMeasures.isNormalized(verticesSpan())
    }

    /**
     * Invert the loop if necessary so that the area enclosed by the loop is at most 2*Pi.
     */
    fun normalize() {
        if (!isNormalized()) invert()
        checkState { isNormalized() }
    }

    /**
     * Reverse the order of the loop vertices, effectively complementing the region represented by the loop.
     * For example, the loop ABCD (with edges AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
     * Notice that the last edge is the same in both cases except that its direction has been reversed.
     */
    fun invert() {
        clearIndex()
        if (isEmptyOrFull()) {
            vertices[0] = if (isFull()) kEmptyVertex() else kFullVertex()
        } else {
            vertices.reverse()
        }
        // origin_inside_ must be set correctly before building the S2ShapeIndex.
        originInside = originInside xor true
        if (bound.lat.lo > -M_PI_2 && bound.lat.hi < M_PI_2) {
            // The complement of this loop contains both poles.
            subregionBound.init(S2LatLngRect.full())
            bound.init(S2LatLngRect.full())
        } else {
            initBound()
        }
        initIndex()
    }

    /**
     * The area of the loop interior, i.e. the region on the left side of the loop. The return value is between 0 and
     * 4*Pi.  (Note that the return value is not affected by whether this loop is a "hole" or a "shell".)
     */
    val area: Double
        get() {
            // S2Loop has its own convention for empty and full loops.
            if (isEmptyOrFull()) {
                return if (originInside) (4 * M_PI) else 0.0
            }
            return S2LoopMeasures.getArea(verticesSpan())
        }

    /**
     * The true centroid of the loop multiplied by the area of the loop (see S2Centroids for details on centroids).
     * The result is not unit length, so you may want to normalize it.  Also note that in general, the centroid may not
     * be contained by the loop.
     *
     * We prescale by the loop area for two reasons:
     *   (1) it is cheaper to compute this way, and
     *   (2) it makes it easier to compute the centroid of more complicated shapes (by splitting them into disjoint
     *       regions and adding their centroids).
     *
     * Note that the return value is not affected by whether this loop is a "hole" or a "shell".
     */
    val centroid: S2Point
        get() = S2LoopMeasures.getCentroid(verticesSpan())

    /**
     * The geodesic curvature of the loop, defined as the sum of the turn angles at each vertex
     * (see S2Measures.turnAngle). The result is positive if the loop is counter-clockwise, negative if the loop is
     * clockwise, and zero if the loop is a great circle. The geodesic curvature is equal to 2*Pi minus the area of the
     * loop.
     *
     * Degenerate and nearly-degenerate loops are handled consistently with S2Predicates.sign(). So for example, if a
     * loop has zero area (i.e., it is a very small CCW loop) then its geodesic curvature will always be positive.
     */
    val curvature: Double
        get() {
            // S2Loop has its own convention for empty and full loops.  For such loops,
            // we return the limit value as the area approaches 0 or 4*Pi respectively.
            if (isEmptyOrFull()) {
                return if (originInside) (-2 * M_PI) else (2 * M_PI)
            }
            return S2LoopMeasures.getCurvature(verticesSpan())
        }

    /**
     * Gets the maximum error in curvature. The return value is not constant; it depends on the loop.
     *
     * @return The curvature max error.
     */
    fun getCurvatureMaxError(): Double = S2LoopMeasures.getCurvatureMaxError(verticesSpan())

    /**
     * Compute the minimum distance from a point "x" to this loop.
     * "x" should be unit length.
     *
     * @return the distance from the given point to the loop interior. If the loop is empty, return S1Angle.infinity().
     */
    fun distance(x: S2Point): S1Angle {
        // Note that S2Loop::Contains(S2Point) is slightly more efficient than the
        // generic version used by S2ClosestEdgeQuery.
        if (contains(x)) return S1Angle.zero()
        return distanceToBoundary(x);
    }

    /**
     * Compute the minimum distance from a point "x" and the loop boundary.
     * "x" should be unit length.
     *
     * @param x A point.
     * @return the distance from the given point to the loop boundary. If the loop is empty or full, return
     * S1Angle.infinity() (since the loop has no boundary).
     */
    fun distanceToBoundary(x: S2Point): S1Angle {
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = false
        val t = S2ClosestEdgeQuery.PointTarget(x)
        return S2ClosestEdgeQuery(index, options).getDistance(t).toAngle()
    }

    /**
     * Projects a point "x" to the loop boundary. If the given point is contained by the loop, return it. Otherwise
     * return the closest point on the loop boundary.
     * If the loop is empty, return the input argument.
     *
     * Note that the result may or may not be contained by the loop.  "x" should be unit length.
     *
     * @param x The point to project.
     * @return The projection of "x" to the loop boundary.
     */
    fun project(x: S2Point): S2Point {
        if (contains(x)) return x
        return projectToBoundary(x)
    }

    /**
     * Project a point "x" to the loop boundary. Returns the closest point on the loop boundary to the given point.
     * If the loop is empty or full, return the input argument (since the loop has no boundary).
     * "x" should be unit length.
     *
     * @param x The point to project.
     * @return The projection of "x" to the loop boundary.
     */
    fun projectToBoundary(x: S2Point): S2Point {
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = false
        val q = S2ClosestEdgeQuery(index, options)
        val target = S2ClosestEdgeQuery.PointTarget(x)
        val edge = q.findClosestEdge(target)
        return q.project(x, edge)
    }

    /**
     * Indicates if this loop contains an other loop.
     *
     * @param b A loop.
     * @return true if the region contained by this loop is a superset of the region contained by the given other loop.
     */
    fun contains(b: S2Loop): Boolean {
        // For this loop A to contains the given loop B, all of the following must
        // be true:
        //
        //  (1) There are no edge crossings between A and B except at vertices.
        //
        //  (2) At every vertex that is shared between A and B, the local edge
        //      ordering implies that A contains B.
        //
        //  (3) If there are no shared vertices, then A must contain a vertex of B
        //      and B must not contain a vertex of A.  (An arbitrary vertex may be
        //      chosen in each case.)
        //
        // The second part of (3) is necessary to detect the case of two loops whose
        // union is the entire sphere, i.e. two loops that contains each other's
        // boundaries but not each other's interiors.
        if (!subregionBound.contains(b.bound)) return false

        // Special cases to handle either loop being empty or full.
        if (isEmptyOrFull() || b.isEmptyOrFull()) {
            return isFull() || b.isEmpty()
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = ContainsRelation()
        if (hasCrossingRelation(this, b, relation)) return false

        // There are no crossings, and if there are any shared vertices then A
        // contains B locally at each shared vertex.
        if (relation.foundSharedVertex()) return true

        // Since there are no edge intersections or shared vertices, we just need to
        // test condition (3) above.  We can skip this test if we discovered that A
        // contains at least one point of B while checking for edge crossings.
        if (!contains(b.vertex(0))) return false

        // We still need to check whether (A union B) is the entire sphere.
        // Normally this check is very cheap due to the bounding box precondition.
        if ((b.subregionBound.contains(bound) || b.bound.union(bound).isFull) && b.contains(vertex(0))) {
            return false
        }
        return true
    }

    // Returns true if the region contained by this loop intersects the region
    // contained by the given other loop.
    fun intersects(b: S2Loop): Boolean {
        // a.intersects(b) if and only if !a.complement().contains(b).
        // This code is similar to contains(), but is optimized for the case
        // where both loops enclose less than half of the sphere.
        if (!bound.intersects(b.bound)) return false

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = IntersectsRelation()
        if (hasCrossingRelation(this, b, relation)) return true
        if (relation.foundSharedVertex()) return false

        // Since there are no edge intersections or shared vertices, the loops
        // intersect only if A contains B, B contains A, or the two loops contain
        // each other's boundaries.  These checks are usually cheap because of the
        // bounding box preconditions.  Note that neither loop is empty (because of
        // the bounding box check above), so it is safe to access vertex(0).

        // Check whether A contains B, or A and B contain each other's boundaries.
        // (Note that A contains all the vertices of B in either case.)
        if (subregionBound.contains(b.bound) || bound.union(b.bound).isFull) {
            if (contains(b.vertex(0))) return true
        }
        // Check whether B contains A.
        if (b.subregionBound.contains(bound) && b.contains(vertex(0))) return true
        return false
    }

    // Returns true if two loops have the same vertices in the same linear order
    // (i.e., cyclic rotations are not allowed).
    override fun equals(other: Any?): Boolean {
        if (other !is S2Loop) return false
        if (numVertices != other.numVertices) return false
        for (i in 0 until numVertices) {
            if (vertex(i) != other.vertex(i)) return false
        }
        return true
    }

    /**
     * Indicates if the boundaries of this loop and loop B are equals. This is true if and only if the loops have the
     * same vertices in the same cyclic order (i.e., the vertices may be cyclically rotated). The empty and full loops
     * are considered to have different boundaries.
     *
     * @param b A loop
     * @return true if two loops have the same boundary.
     */
    fun boundaryEquals(b: S2Loop): Boolean {
        if (numVertices != b.numVertices) return false

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (isEmptyOrFull()) return isEmpty() == b.isEmpty()

        for (offset in 0 until numVertices) {
            if (vertex(offset) == b.vertex(0)) {
                // There is at most one starting offset since loop vertices are unique.
                for (i in 0 until numVertices) {
                    if (vertex(i + offset) != b.vertex(i)) return false
                }
                return true
            }
        }
        return false
    }

    /**
     * Indicates if the boundaries of this loop and loop B are approx equals. More precisely, the vertices in the two
     * loops must be in the same cyclic order, and corresponding vertex pairs must be separated by no more than
     * "maxError".
     *
     * @param b A loop.
     * @param maxError Max authorized error.
     * @return true if two loops have the same boundary except for vertex perturbations.
     */
    fun boundaryApproxEquals(b: S2Loop, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        if (numVertices != b.numVertices) return false

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (isEmptyOrFull()) return isEmpty() == b.isEmpty()

        for (offset in 0 until numVertices) {
            if (S2PointUtil.approxEquals(vertex(offset), b.vertex(0), maxError)) {
                var success = true
                for (i in 0 until numVertices) {
                    if (!S2PointUtil.approxEquals(vertex(i + offset), b.vertex(i), maxError)) {
                        success = false
                        break
                    }
                }
                if (success) return true
                // Otherwise continue looping.  There may be more than one candidate
                // starting offset since vertices are only matched approximately.
            }
        }
        return false
    }

    /**
     * Indicates if the boundary of this loop is near the loop b. The two loops may have different numbers of vertices.
     * More precisely, this method returns true if the two loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2
     * such that distance(a(t), b(t)) <= max_error for all t.  You can think of this as testing whether it is possible
     * to drive two cars all the way around the two loops such that no car ever goes backward and the cars are always
     * within "max_error" of each other.
     *
     * @param b A loop.
     * @param maxError The max loop distance.
     * @return true if the two loop boundaries are within "maxError" of each other along their entire lengths.
     */
    fun boundaryNear(b: S2Loop, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        // Special case to handle empty or full loops.
        if (isEmptyOrFull() || b.isEmptyOrFull()) {
            return (isEmpty() && b.isEmpty()) || (isFull() && b.isFull())
        }

        for (aOffset in 0 until numVertices) {
            if (matchBoundaries(this, b, aOffset, maxError)) return true
        }
        return false
    }

    /////////////////////////// S2Region ///////////////////////////

    override fun clone(): S2Loop = S2Loop(this)

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.

    override val capBound: S2Cap
        get() = bound.capBound

    override val rectBound: S2LatLngRect
        get() = bound

    override fun contains(cell: S2Cell): Boolean {
        val it = AbstractMutableS2ShapeIndex.Iterator(index)
        val relation = it.locate(cell.id())

        // If "target" is disjoint from all index cells, it is not contained.
        // Similarly, if "target" is subdivided into one or more index cells then it
        // is not contained, since index cells are subdivided only if they (nearly)
        // intersect a sufficient number of edges.  (But note that if "target" itself
        // is an index cell then it may be contained, since it could be a cell with
        // no edges in the loop interior.)
        if (relation != INDEXED) return false

        // Otherwise check if any edges intersect "target".
        if (boundaryApproxIntersects(it, cell)) return false

        // Otherwise check if the loop contains the center of "target".
        return contains(it, cell.getCenter())
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        val it = AbstractMutableS2ShapeIndex.Iterator(index)
        val relation = it.locate(cell.id())

        // If "target" does not overlap any index cell, there is no intersection.
        if (relation == CellRelation.DISJOINT) return false

        // If "target" is subdivided into one or more index cells, there is an
        // intersection to within the S2ShapeIndex error bound (see Contains).
        if (relation == CellRelation.SUBDIVIDED) return true

        // If "target" is an index cell, there is an intersection because index cells
        // are created only if they have at least one edge or they are entirely
        // contained by the loop.
        if (it.id() == cell.id()) return true

        // Otherwise check if any edges intersect "target".
        if (boundaryApproxIntersects(it, cell)) return true

        // Otherwise check if the loop contains the center of "target".
        return contains(it, cell.getCenter())
    }

    override fun contains(p: S2Point): Boolean {
        // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
        // worthwhile only when it might allow us to delay building the index.
        if (!index.isFresh() && !bound.contains(p)) return false

        // For small loops it is faster to just check all the crossings.  We also
        // use this method during loop initialization because InitOriginAndBound()
        // calls Contains() before InitIndex().  Otherwise, we keep track of the
        // number of calls to Contains() and only build the index when enough calls
        // have been made so that we think it is worth the effort.  Note that the
        // code below is structured so that if many calls are made in parallel only
        // one thread builds the index, while the rest continue using brute force
        // until the index is actually available.
        //
        // The constants below were tuned using the benchmarks.  It turns out that
        // building the index costs roughly 50x as much as Contains().  (The ratio
        // increases slowly from 46x with 64 points to 61x with 256k points.)  The
        // textbook approach to this problem would be to wait until the cumulative
        // time we would have saved with an index approximately equals the cost of
        // building the index, and then build it.  (This gives the optimal
        // competitive ratio of 2; look up "competitive algorithms" for details.)
        // We set the limit somewhat lower than this (20 rather than 50) because
        // building the index may be forced anyway by other API calls, and so we
        // want to err on the side of building it too early.
        if (index.nextNewShapeId() == 0 ||  // InitIndex() not called yet
            numVertices <= kMaxBruteForceVertices ||
            (!index.isFresh() && unindexedContainsCalls.incrementAndGet() != kMaxUnindexedContainsCalls)
        ) {
            return bruteForceContains(p)
        }
        // Otherwise we look up the S2ShapeIndex cell containing this point.  Note
        // the index is built automatically the first time an iterator is created.
        val it = AbstractMutableS2ShapeIndex.Iterator(index)
        if (!it.locate(p)) return false
        return contains(it, p)
    }

    ////////////////////////////////////////////////////////////////////////
    // Methods intended primarily for use by the S2Polygon implementation:

    /**
     * Given two loops of a polygon, return true if A contains B. This version of Contains() is cheap because it does
     * not test for edge intersections. The loops must meet all the S2Polygon requirements; for example this implies
     * that their boundaries may not cross or have any shared edges (although they may have shared vertices).
     *
     * @param b A loop of the same polygon.
     * @return true if A contains B.
     */
    fun containsNested(b: S2Loop): Boolean {
        if (!subregionBound.contains(b.bound)) return false

        // Special cases to handle either loop being empty or full.  Also bail out
        // when B has no vertices to avoid heap overflow on the vertex(1) call
        // below.  (This method is called during polygon initialization before the
        // client has an opportunity to call IsValid().)
        if (isEmptyOrFull() || b.numVertices < 2) {
            return isFull() || b.isEmpty()
        }

        // We are given that A and B do not share any edges, and that either one
        // loop contains the other or they do not intersect.
        val m = findVertex(b.vertex(1))
        if (m < 0) {
            // Since b->vertex(1) is not shared, we can check whether A contains it.
            return contains(b.vertex(1))
        }
        // Check whether the edge order around b->vertex(1) is compatible with
        // A containing B.
        return S2WedgeRelations.wedgeContains(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(0), b.vertex(2))
    }

    /**
     * Compares the boundaries of this loop and an other one.
     * Returns +1 if A contains the boundary of B, -1 if A excludes the boundary of B, and 0 if the boundaries of A and
     * B cross. Shared edges are handled as follows:
     * If XY is a shared edge, define Reversed(XY) to be true if XY appears in opposite directions in A and B. Then A
     * contains XY if and only if Reversed(XY) == B.isHole(). (Intuitively, this checks whether A contains a
     * vanishingly small region extending from the boundary of B toward the interior of the polygon to which loop B
     * belongs.)
     *
     * This method is used for testing containment and intersection of multi-loop polygons. Note that this method is
     * not symmetric, since the result depends on the direction of loop A but not on the direction of loop B
     * (in the absence of shared edges).
     *
     * REQUIRES: neither loop is empty.
     * REQUIRES: if b.isFull, then !b.isHole().
     *
     * @param b A loop
     * @return The comparison result (+1 if A contains the boundary of B, -1 if A excludes the boundary of B and 0
     * if boundaries of A and B crosses)
     */
    fun compareBoundary(b: S2Loop): Int {
        requireArgument { !isEmpty() && !b.isEmpty() }
        requireArgument { !b.isFull() || !b.isHole() }

        // The bounds must intersect for containment or crossing.
        if (!bound.intersects(b.bound)) return -1

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (isFull()) return 1
        if (b.isFull()) return -1

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = CompareBoundaryRelation(b.isHole())
        if (hasCrossingRelation(this, b, relation)) return 0
        if (relation.foundSharedVertex()) {
            return if (relation.containsEdge()) 1 else -1
        }

        // There are no edge intersections or shared vertices, so we can check
        // whether A contains an arbitrary vertex of B.
        return if (contains(b.vertex(0))) 1 else -1
    }

    /**
     * Given two loops whose boundaries do not cross (see CompareBoundary), return true if A contains the boundary of B.
     * If "reverse_b" is true, the boundary of B is reversed first (which only affects the result when there are shared
     * edges). This method is cheaper than compareBoundary() because it does not test for edge intersections.
     *
     * REQUIRES: neither loop is empty.
     * REQUIRES: if b.isFull(), then reverseB == false.
     *
     * @param b A loop
     * @param reverseB Indicates if B must be reversed first.
     * @return true if A contains the boundary of B and false otherwise.
     */
    fun containsNonCrossingBoundary(b: S2Loop, reverseB: Boolean): Boolean {
        requireArgument { !isEmpty() && !b.isEmpty() }
        requireArgument { b.isFull() || !reverseB }

        // The bounds must intersect for containment.
        if (!bound.intersects(b.bound)) return false

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (isFull()) return true
        if (b.isFull()) return false

        val m = findVertex(b.vertex(0))
        if (m < 0) {
            // Since vertex b0 is not shared, we can check whether A contains it.
            return contains(b.vertex(0))
        }
        // Otherwise check whether the edge (b0, b1) is contained by A.
        return wedgeContainsSemiwedge(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(1), reverseB)
    }

    /////////////////////////// Internal ///////////////////////////

    // Returns an S2PointLoopSpan containing the loop vertices, for use with the
    // functions defined in s2loop_measures.h.
    fun verticesSpan(): S2PointLoopSpan = S2PointLoopSpan(vertices)

    private fun initOriginAndBound() {
        if (numVertices < 3) {
            // Check for the special empty and full loops (which have one vertex).
            if (!isEmptyOrFull()) {
                originInside = false
                return  // Bail out without trying to access non-existent vertices.
            }
            // If the vertex is in the southern hemisphere then the loop is full,
            // otherwise it is empty.
            originInside = (vertex(0).z < 0)
        } else {
            // Point containment testing is done by counting edge crossings starting
            // at a fixed point on the sphere (S2::Origin()).  Historically this was
            // important, but it is now no longer necessary, and it may be worthwhile
            // experimenting with using a loop vertex as the reference point.  In any
            // case, we need to know whether the reference point (S2::Origin) is
            // inside or outside the loop before we can construct the S2ShapeIndex.
            // We do this by first guessing that it is outside, and then seeing
            // whether we get the correct containment result for vertex 1.  If the
            // result is incorrect, the origin must be inside the loop.
            //
            // A loop with consecutive vertices A,B,C contains vertex B if and only if
            // the fixed vector R = S2::Ortho(B) is contained by the wedge ABC.  The
            // wedge is closed at A and open at C, i.e. the point B is inside the loop
            // if A=R but not if C=R.  This convention is required for compatibility
            // with S2::VertexCrossing.  (Note that we can't use S2::Origin()
            // as the fixed vector because of the possibility that B == S2::Origin().)
            //
            // TODO(ericv): Investigate using vertex(0) as the reference point.

            originInside = false  // Initialize before calling Contains().
            val v1Inside = S2Predicates.orderedCCW(S2PointUtil.ortho(vertex(1)), vertex(0), vertex(2), vertex(1))
            // Note that Contains(S2Point) only does a bounds check once InitIndex()
            // has been called, so it doesn't matter that bound_ is undefined here.
            if (v1Inside != contains(vertex(1))) {
                originInside = true
            }
            logger.trace { "Origin inside = $originInside" }
        }

        // We *must* call InitBound() before InitIndex(), because InitBound() calls
        // Contains(S2Point), and Contains(S2Point) does a bounds check whenever the
        // index is not fresh (i.e., the loop has been added to the index but the
        // index has not been updated yet).
        //
        // TODO(ericv): When fewer S2Loop methods depend on internal bounds checks,
        // consider computing the bound on demand as well.
        initBound()
        initIndex()
    }

    private fun initBound() {
        // Check for the special empty and full loops.
        if (isEmpty()) {
            subregionBound.setEmpty()
            bound.setEmpty()
            logger.trace { "Loop is empty: subregionBound = $subregionBound, bound = $bound" }
            return
        } else if (isFull()) {
            subregionBound.setFull()
            bound.setFull()
            logger.trace { "Loop is full: subregionBound = $subregionBound, bound = $bound" }
            return
        }

        // The bounding rectangle of a loop is not necessarily the same as the
        // bounding rectangle of its vertices.  First, the maximal latitude may be
        // attained along the interior of an edge.  Second, the loop may wrap
        // entirely around the sphere (e.g. a loop that defines two revolutions of a
        // candy-cane stripe).  Third, the loop may include one or both poles.
        // Note that a small clockwise loop near the equator contains both poles.

        val bounder = S2LatLngRectBounder()
        for (i in 0..numVertices) {
            bounder.addPoint(vertex(i))
        }
        var b = bounder.getBound()
        if (contains(S2Point(0, 0, 1))) {
            b = S2LatLngRect(R1Interval(b.lat.lo, M_PI_2), S1Interval(S1Interval.full()))
        }
        // If a loop contains the south pole, then either it wraps entirely
        // around the sphere (full longitude range), or it also contains the
        // north pole in which case b.lng().is_full() due to the test above.
        // Either way, we only need to do the south pole containment test if
        // b.lng().is_full().
        if (b.lng.isFull && contains(S2Point(0, 0, -1))) {
            b.lat.lo = -M_PI_2
        }
        bound.init(b)
        subregionBound.init(S2LatLngRectBounder.expandForSubregions(bound))

        logger.trace { "Loop bounds: subregionBound = $subregionBound, bound = $bound" }
    }

    private fun initIndex() {
        index.add(Shape(loop = this))
        if (!lazyIndexing) {
            index.forceBuild()
        }
        if (PreConditions.enabled && debugOverride == S2Debug.ALLOW) {
            // Note that FLAGS_s2debug is false in optimized builds (by default).
            checkState { isValid() }
        }
    }

    // A version of Contains(S2Point) that does not use the S2ShapeIndex.
    // Used by the S2Polygon implementation.
    internal fun bruteForceContains(p: S2Point): Boolean {
        logger.trace { "bruteForceContains(p = $p)" }
        // Empty and full loops don't need a special case, but invalid loops with
        // zero vertices do, so we might as well handle them all at once.
        if (numVertices < 3) return originInside

        val origin = S2PointUtil.origin()
        val crosser = S2EdgeCrosser(origin, p, vertex(0))
        var inside = originInside
        for (i in 1..numVertices) {
            val edgeOrVertexCrossing = crosser.edgeOrVertexCrossing(vertex(i))
            logger.trace { "EdgeOrVertexCrossign $i - ${vertex(i)}: $edgeOrVertexCrossing" }
            inside = inside xor edgeOrVertexCrossing
        }

        logger.trace { "bruteForceContains(p = $p) = $inside" }
        return inside
    }

    /**
     * Like findValidationError(), but skips any checks that would require building the S2ShapeIndex
     * (i.e., self-intersection tests). This is used by the S2Polygon implementation, which uses its own index to
     * check for loop self-intersections.
     *
     * @param error The error to update according the validation state.
     * @return true is the loop is valid and false otherwise.
     */
    fun findValidationErrorNoIndex(error: S2Error): Boolean {
        // subregionBound must be at least as large as Bound.  (This is an internal consistency check rather than a
        // test of client data.)
        checkState { subregionBound.contains(bound) }

        // All vertices must be unit length.  (Unfortunately this check happens too
        // late in debug mode, because S2Loop construction calls s2pred::Sign which
        // expects vertices to be unit length.  But it is still a useful check in
        // optimized builds.)
        for (i in 0 until numVertices) {
            if (!S2PointUtil.isUnitLength(vertex(i))) {
                error.init(S2Error.NOT_UNIT_LENGTH, "Vertex $i is not unit length")
                return true
            }
        }
        // Loops must have at least 3 vertices (except for the empty and full loops).
        if (numVertices < 3) {
            if (isEmptyOrFull()) {
                return false  // Skip remaining tests.
            }
            error.init(S2Error.LOOP_NOT_ENOUGH_VERTICES, "Non-empty, non-full loops must have at least 3 vertices")
            return true
        }
        // Loops are not allowed to have any duplicate vertices or edge crossings.
        // We split this check into two parts.  First we check that no edge is
        // degenerate (identical endpoints).  Then we check that there are no
        // intersections between non-adjacent edges (including at vertices).  The
        // second part needs the S2ShapeIndex, so it does not fall within the scope
        // of this method.
        for (i in 0 until numVertices) {
            if (vertex(i) == vertex(i + 1)) {
                error.init(S2Error.DUPLICATE_VERTICES, "Edge $i is degenerate (duplicate vertex)")
                return true
            }
            if (vertex(i) == -vertex(i + 1)) {
                error.init(S2Error.ANTIPODAL_VERTICES, "Vertices $i and ${(i + 1) % numVertices} are antipodal")
                return true
            }
        }
        return false
    }

    // Converts the loop vertices to the S2XYZFaceSiTi format and store the result
    // in the given array, which must be large enough to store all the vertices.
    private fun getXYZFaceSiTiVertices(vertices: MutableList<S2XYZFaceSiTi>): List<S2XYZFaceSiTi> {
        vertices.addAll((0 until numVertices).map { i ->
            val xyz = vertex(i)
            val (cellLevel, faceSiTi) = S2Coords.xyzToFaceSiTi(xyz)
            S2XYZFaceSiTi(xyz = xyz, face = faceSiTi.face, si = faceSiTi.si, ti = faceSiTi.ti, cellLevel = cellLevel)
        })
        return vertices
    }

    // Given an iterator that is already positioned at the S2ShapeIndexCell
    // containing "p", returns Contains(p).
    private fun contains(it: AbstractMutableS2ShapeIndex.Iterator, p: S2Point): Boolean {
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        val aClipped = it.cell().clipped(0)
        var inside = aClipped.containsCenter
        val aNumEdges = aClipped.numEdges
        if (aNumEdges > 0) {
            val center = it.center()
            val crosser = S2EdgeCrosser(center, p)
            var aiPrev = -2
            for (i in 0 until aNumEdges) {
                val ai = aClipped.edge(i)
                if (ai != aiPrev + 1) crosser.restartAt(vertex(ai))
                aiPrev = ai
                inside = inside xor crosser.edgeOrVertexCrossing(vertex(ai + 1))
            }
        }
        return inside
    }

    val kMaxError = (kFaceClipErrorUVCoord + kIntersectsRectErrorUVDist)

    // Returns true if the loop boundary intersects "target".  It may also
    // return true when the loop boundary does not intersect "target" but
    // some edge comes within the worst-case error tolerance.
    //
    // REQUIRES: it.id().contains(target.id())
    // [This condition is true whenever it.Locate(target) returns INDEXED.]
    private fun boundaryApproxIntersects(it: AbstractMutableS2ShapeIndex.Iterator, target: S2Cell): Boolean {
        requireArgument { it.id().contains(target.id()) }
        val a_clipped = it.cell().clipped(0)
        val a_num_edges = a_clipped.numEdges

        // If there are no edges, there is no intersection.
        if (a_num_edges == 0) return false

        // We can save some work if "target" is the index cell itself.
        if (it.id() == target.id()) return true

        // Otherwise check whether any of the edges intersect "target".
        val bound = target.boundUV().expanded(kMaxError)
        for (i in 0 until a_num_edges) {
            val ai = a_clipped.edge(i)
            val clippedEdge = S2EdgeClipping.clipToPaddedFace(vertex(ai), vertex(ai + 1), target.face(), kMaxError)
            if (clippedEdge != null && S2EdgeClipping.intersectsRect(clippedEdge.first, clippedEdge.second, bound)) {
                return true
            }
        }

        return false
    }

    // Returns an index "first" and a direction "dir" such that the vertex
    // sequence (first, first + dir, ..., first + (n - 1) * dir) does not change
    // when the loop vertex order is rotated or reversed.  This allows the loop
    // vertices to be traversed in a canonical order.
    internal fun getCanonicalLoopOrder(): S2LoopMeasures.LoopOrder = S2LoopMeasures.getCanonicalLoopOrder(verticesSpan())

    // Returns the index of a vertex at point "p", or -1 if not found.
    // The return value is in the range 1..num_vertices_ if found.
    private fun findVertex(p: S2Point): Int {
        if (numVertices < 10) {
            // Exhaustive search.  Return value must be in the range [1..N].
            for (i in 1..numVertices) {
                if (vertex(i) == p) return i
            }
            return -1
        }
        val it = AbstractMutableS2ShapeIndex.Iterator(index)
        if (!it.locate(p)) return -1

        val aClipped = it.cell().clipped(0)
        var i = aClipped.numEdges - 1
        while (i >= 0) {
            val ai = aClipped.edge(i)
            // Return value must be in the range [1..N].
            if (vertex(ai) == p) return if (ai == 0) numVertices else ai
            if (vertex(ai + 1) == p) return ai + 1
            --i
        }
        return -1
    }

    /**
     * When the loop is modified (Invert(), or Init() called again) then the indexing structures need to be cleared
     * since they become invalid.
     */
    private fun clearIndex() {
        unindexedContainsCalls.set(0)
        index.removeAll()
    }

    companion object {

        private val logger = KotlinLogging.logger(S2Loop::class.java.name)

        /**
         * A special vertex chain of length 1 that creates an empty loop (i.e., a loop with no edges that contains no
         * points).  Example usage:
         *
         *    val empty = S2Loop(S2Loop.kEmpty)
         */
        val kEmpty: List<S2Point> = listOf(kEmptyVertex())

        /**
         * A special vertex chain of length 1 that creates a full loop (i.e., a loop with no edges that contains all
         * points).  See kEmpty for details.
         */
        val kFull: List<S2Point> = listOf(kFullVertex())

        // The single vertex in the "empty loop" vertex chain.
        fun kEmptyVertex(): S2Point = S2Point(0, 0, 1)

        // The single vertex in the "full loop" vertex chain.
        fun kFullVertex(): S2Point = S2Point(0, 0, -1)

        private val kMaxBruteForceVertices = 32
        private val kMaxUnindexedContainsCalls = 20

        /**
         * Build the S2ShapeIndex only when it is first needed. This can save significant amounts of memory and time
         * when geometry is constructed but never queried, for example when loops are passed directly to S2Polygon, or
         * when geometry is being converted from one format to another.
         */
        private var lazyIndexing: Boolean = true

        /**
         * Constructs a regular polygon with the given number of vertices, all located on a circle of the specified
         * radius around "center".  The radius is the actual distance from "center" to each vertex.
         *
         * @param center The circle center.
         * @param radius The circle radius
         * @param numVertices The number of verticles.
         *
         * @return The built loop.
         */
        fun makeRegularLoop(center: S2Point, radius: S1Angle, numVertices: Int): S2Loop {
            val m = S2PointUtil.getFrame(center)
            return makeRegularLoop(m, radius, numVertices)
        }

        // Like the function above, but this version constructs a loop centered
        // around the z-axis of the given coordinate frame, with the first vertex in
        // the direction of the positive x-axis.  (This allows the loop to be
        // rotated for testing purposes.)
        fun makeRegularLoop(frame: Matrix3x3Double, radius: S1Angle, numVertices: Int): S2Loop {
            // We construct the loop in the given frame coordinates, with the center at
            // (0, 0, 1).  For a loop of radius "r", the loop vertices have the form
            // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r).  The distance on the
            // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
            val z = cos(radius.radians)
            val r = sin(radius.radians)
            val radian_step = 2 * M_PI / numVertices
            val vertices = mutableListOf<S2Point>()
            for (i in 0 until numVertices) {
                val angle = i * radian_step
                val p = S2Point(r * cos(angle), r * sin(angle), z)
                vertices.add(S2PointUtil.fromFrame(frame, p).normalize())
            }
            return S2Loop(vertices)
        }

        // This method checks all edges of loop A for intersection against all edges
        // of loop B.  If there is any shared vertex, the wedges centered at this
        // vertex are sent to "relation".
        //
        // If the two loop boundaries cross, this method is guaranteed to return
        // true.  It also returns true in certain cases if the loop relationship is
        // equivalent to crossing.  For example, if the relation is Contains() and a
        // point P is found such that B contains P but A does not contain P, this
        // method will return true to indicate that the result is the same as though
        // a pair of crossing edges were found (since Contains() returns false in
        // both cases).
        //
        // See Contains(), Intersects() and CompareBoundary() for the three uses of
        // this function.
        private fun hasCrossingRelation(a: S2Loop, b: S2Loop, relation: LoopRelation): Boolean {
            // We look for S2CellId ranges where the indexes of A and B overlap, and
            // then test those edges for crossings.
            val ai = RangeIterator(a.index)
            val bi = RangeIterator(b.index)
            val ab = LoopCrosser(a, b, relation, false)  // Tests edges of A against B
            val ba = LoopCrosser(b, a, relation, true)   // Tests edges of B against A
            while (!ai.done() || !bi.done()) {
                if (ai.rangeMax() < bi.rangeMin()) {
                    // The A and B cells don't overlap, and A precedes B.
                    ai.seekTo(bi)
                } else if (bi.rangeMax() < ai.rangeMin()) {
                    // The A and B cells don't overlap, and B precedes A.
                    bi.seekTo(ai)
                } else {
                    // One cell contains the other.  Determine which cell is larger.
                    val aiLsb = ai.id().lsb()
                    val biLsb = bi.id().lsb()
                    if (aiLsb > biLsb) {
                        // A's index cell is larger.
                        if (ab.hasCrossingRelation(ai, bi)) return true
                    } else if (aiLsb < biLsb) {
                        // B's index cell is larger.
                        if (ba.hasCrossingRelation(bi, ai)) return true
                    } else {
                        // The A and B cells are the same.  Since the two cells have the same
                        // center point P, check whether P satisfies the crossing targets.
                        if (ai.containsCenter(0) && ab.aCrossingTarget == 1 &&
                            bi.containsCenter(0) && ab.bCrossingTarget == 1
                        ) {
                            return true;
                        }
                        // Otherwise test all the edge crossings directly.
                        if (ai.numEdges(0) > 0 && bi.numEdges(0) > 0 &&
                            ab.cellCrossesCell(ai.clipped(0), bi.clipped(0))
                        ) {
                            return true;
                        }
                        ai.next()
                        bi.next()
                    }
                }
            }
            return false
        }

        private fun matchBoundaries(a: S2Loop, b: S2Loop, a_offset: Int, max_error: S1Angle): Boolean {
            // The state consists of a pair (i,j).  A state transition consists of
            // incrementing either "i" or "j".  "i" can be incremented only if
            // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
            // applies to "j".  The function returns true iff we can proceed all the way
            // around both loops in this way.
            //
            // Note that when "i" and "j" can both be incremented, sometimes only one
            // choice leads to a solution.  We handle this using a stack and
            // backtracking.  We also keep track of which states have already been
            // explored to avoid duplicating work.

            val pending = mutableListOf<Pair<Int, Int>>()
            val done = mutableSetOf<Pair<Int, Int>>()
            pending.add(0 to 0)
            while (pending.isNotEmpty()) {
                val (i, j) = pending.removeLast()
                if (i == a.numVertices && j == b.numVertices) {
                    return true
                }
                done.add(i to j)

                // If (i == na && offset == na-1) where na == a->num_vertices(), then
                // then (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
                // So we reduce the range if necessary.
                var io = i + a_offset
                if (io >= a.numVertices) io -= a.numVertices

                if (i < a.numVertices && !done.contains((i + 1) to j) &&
                    S2EdgeDistances.getDistance(a.vertex(io + 1), b.vertex(j), b.vertex(j + 1)) <= max_error) {
                    pending.add((i + 1) to j)
                }
                if (j < b.numVertices && !done.contains(i to (j + 1)) &&
                    S2EdgeDistances.getDistance(b.vertex(j + 1), a.vertex(io), a.vertex(io + 1)) <= max_error) {
                    pending.add(i to (j + 1))
                }
            }
            return false
        }
    }


    // Wrapper class for indexing a loop (see S2ShapeIndex).  Once this object
    // is inserted into an S2ShapeIndex it is owned by that index, and will be
    // automatically deleted when no longer needed by the index.  Note that this
    // class does not take ownership of the loop itself (see OwningShape below).
    // You can also subtype this class to store additional data (see S2Shape for
    // details).
    class Shape(id: Int = 0, val loop: S2Loop) : S2Shape(id) {

        // S2Shape interface:

        override val numEdges: Int
            get() = if (loop.isEmptyOrFull()) 0 else loop.numVertices

        override fun edge(edgeId: Int): Edge = Edge(loop.vertex(edgeId), loop.vertex(edgeId + 1))

        override val dimension: Int = 2

        override fun getReferencePoint(): ReferencePoint = ReferencePoint(S2PointUtil.origin(), loop.originInside)

        override val numChains: Int
            get() = if (loop.isEmpty()) 0 else 1

        override fun chain(chainId: Int): Chain {
            requireEQ(chainId, 0)
            return Chain(0, numEdges)
        }

        override fun chainEdge(chainId: Int, offset: Int): Edge {
            requireEQ(chainId, 0)
            return Edge(loop.vertex(offset), loop.vertex(offset + 1))
        }

        override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

        override val typeTag: TypeTag = kNoTypeTag

        override fun toString(): String {
            return "S2Loop.Shape(id=$id)"
        }

    }

}

// List of vertices with operator[] maps index values in the range
// [n, 2*n-1] to the range [0, n-1] by subtracting n (where n == size()).
// In other words, two full copies of the vertex array are available.  (This
// is a compromise between convenience and efficiency, since computing the
// index modulo "n" is surprisingly expensive.)
//
// This property is useful for implementing algorithms where the elements of
// the span represent the vertices of a loop.
class S2PointLoopSpan(val points: List<S2Point> = emptyList()) : List<S2Point> by points {

    // Like operator[], but allows index values in the range [0, 2*size()-1]
    // where each index i >= size() is mapped to i - size().
    override operator fun get(index: Int): S2Point {
        requireGE(index, 0)
        requireLT(index, 2 * points.size)
        val j = index - points.size
        return points[if (j < 0) index else j]
    }

    fun rotated(startIndex: Int): S2PointLoopSpan {
        requireGE(startIndex, 0)
        requireLT(startIndex, points.size)
        val vertices = mutableListOf<S2Point>()
        for (i in startIndex until (startIndex + points.size)) {
            vertices.add(this[i])
        }
        return S2PointLoopSpan(vertices)
    }

    override fun toString(): String = points.toString()
}
