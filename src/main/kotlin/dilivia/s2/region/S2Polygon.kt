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
import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireEQ
import dilivia.PreConditions.requireLT
import dilivia.collections.upperBound
import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2TextParser
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.layers.S2PolygonLayer
import dilivia.s2.builder.layers.S2PolylineLayer
import dilivia.s2.builder.layers.S2PolylineVectorLayer
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.S2CellIdSnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2BooleanOperation
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.index.shape.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.index.shape.S2CrossingEdgePairsScanner
import dilivia.s2.region.S2ShapeIndexRegion.Companion.makeS2ShapeIndexRegion
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.TypeTag
import dilivia.s2.shape.TypeTags
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.sqrt
import java.util.*
import java.util.concurrent.atomic.AtomicInteger


// A map from each loop to its immediate children with respect to nesting.
// This map is built during initialization of multi-loop polygons to
// determine which are shells and which are holes, and then discarded.
typealias LoopMap = MutableMap<S2Loop?, MutableList<S2Loop>>

// Effectively create a priority queue of polygons in order of number of
// vertices.  Repeatedly union the two smallest polygons and add the result
// to the queue until we have a single polygon to return.
typealias QueueType = TreeMap<Int, MutableList<S2Polygon>>

fun QueueType.polygonCount(): Int = this.entries.map { entry -> entry.value.size }.sum()

// An S2Polygon is an S2Region object that represents a polygon.  A polygon is
// defined by zero or more loops; recall that the interior of a loop is
// defined to be its left-hand side (see S2Loop).  There are two different
// conventions for creating an S2Polygon:
//
//   - InitNested() expects the input loops to be nested hierarchically.  The
//     polygon interior then consists of the set of points contained by an odd
//     number of loops.  So for example, a circular region with a hole in it
//     would be defined as two CCW loops, with one loop containing the other.
//     The loops can be provided in any order.
//
//     When the orientation of the input loops is unknown, the nesting
//     requirement is typically met by calling S2Loop::Normalize() on each
//     loop (which inverts the loop if necessary so that it encloses at most
//     half the sphere).  But in fact any set of loops can be used as long as
//     (1) there is no pair of loops that cross, and (2) there is no pair of
//     loops whose union is the entire sphere.
//
//   - InitOriented() expects the input loops to be oriented such that the
//     polygon interior is on the left-hand side of every loop.  So for
//     example, a circular region with a hole in it would be defined using a
//     CCW outer loop and a CW inner loop.  The loop orientations must all be
//     consistent; for example, it is not valid to have one CCW loop nested
//     inside another CCW loop, because the region between the two loops is on
//     the left-hand side of one loop and the right-hand side of the other.
//
// Most clients will not call these methods directly; instead they should use
// S2Builder, which has better support for dealing with imperfect data.
//
// When the polygon is initialized, the given loops are automatically
// converted into a canonical form consisting of "shells" and "holes".  Shells
// and holes are both oriented CCW, and are nested hierarchically.  The loops
// are reordered to correspond to a preorder traversal of the nesting
// hierarchy; InitOriented may also invert some loops. The set of input S2Loop
// pointers is always preserved; the caller can use this to determine how the
// loops were reordered if desired.
//
// Polygons may represent any region of the sphere with a polygonal boundary,
// including the entire sphere (known as the "full" polygon).  The full
// polygon consists of a single full loop (see S2Loop), whereas the empty
// polygon has no loops at all.
//
// Polygons have the following restrictions:
//
//  - Loops may not cross, i.e. the boundary of a loop may not intersect
//    both the interior and exterior of any other loop.
//
//  - Loops may not share edges, i.e. if a loop contains an edge AB, then
//    no other loop may contain AB or BA.
//
//  - Loops may share vertices, however no vertex may appear twice in a
//    single loop (see S2Loop).
//
//  - No loop may be empty.  The full loop may appear only in the full polygon.

// The default constructor creates an empty polygon.  It can be made
// non-empty by calling Init(), Decode(), etc.
class S2Polygon() : S2Region {

    /////////////////////////// Fields ///////////////////////////

    /**
     * Allows overriding the automatic validity checks controlled by the -ea flag. If this flag is true, then loops are
     * automatically checked for validity as they are initialized.
     * Without the call to loop.debugOverride = S2Debug.DISABLE, invalid data would cause a fatal error in init()
     * whenever the -ea flag is enabled.
     */
    var debugOverride: S2Debug = S2Debug.ALLOW

    private val loops: MutableList<S2Loop> = mutableListOf()

    // Cache for num_vertices().
    private var numVertices: Int = 0

    // "bound_" is a conservative bound on all points contained by this polygon:
    // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
    private var bound: S2LatLngRect = S2LatLngRect.empty()

    // Since "bound_" is not exact, it is possible that a polygon A contains
    // another polygon B whose bounds are slightly larger.  "subregion_bound_"
    // has been expanded sufficiently to account for this error, i.e.
    // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
    private var subregionBound: S2LatLngRect = S2LatLngRect.empty()

    // True if InitOriented() was called and the given loops had inconsistent
    // orientations (i.e., it is not possible to construct a polygon such that
    // the interior is on the left-hand side of all loops).  We need to remember
    // this error so that it can be returned later by FindValidationError(),
    // since it is not possible to detect this error once the polygon has been
    // initialized.  This field is not preserved by Encode/Decode.
    private var errorInconsistentLoopOrientations: Boolean = false

    // Spatial index containing this polygon.
    // Returns the built-in S2ShapeIndex associated with every S2Polygon.  This
    // can be used in conjunction with the various S2ShapeIndex query classes
    // (S2ClosestEdgeQuery, S2BooleanOperation, etc) to do things beyond what is
    // possible with S2Polygon built-in convenience methods.
    //
    // For example, to measure the distance from one S2Polygon to another, you
    // can write:
    //   S2ClosestEdgeQuery query(&polygon1.index());
    //   S2ClosestEdgeQuery::ShapeIndexTarget target(&polygon2.index());
    //   S1ChordAngle distance = query.GetDistance(&target);
    //
    // The index contains a single S2Polygon::Shape object.
    val index: MutableS2ShapeIndex = MutableS2ShapeIndex()

    // In general we build the index the first time it is needed, but we make an
    // exception for Contains(S2Point) because this method has a simple brute
    // force implementation that is also relatively cheap.  For this one method
    // we keep track of the number of calls made and only build the index once
    // enough calls have been made that we think an index would be worthwhile.
    val unindexedContainsCalls: AtomicInteger = AtomicInteger(0)

    /////////////////////////// Init /////////////////////////////

    // Convenience constructor that calls InitNested() with the given loops.
    //
    // When called with override == S2Debug::ALLOW, the automatic validity
    // checking is controlled by --s2debug (which is true by default in
    // non-optimized builds).  When this flag is enabled, a fatal error is
    // generated whenever an invalid polygon is constructed.
    //
    // With override == S2Debug::DISABLE, the automatic validity checking
    // is disabled.  The main reason to do this is if you intend to call
    // IsValid() explicitly.  (See set_s2debug_override() for details.)
    // Example:
    //
    //   std::vector<std::unique_ptr<S2Loop>> loops;
    //   // ... set up loops ...
    //   S2Polygon* polygon = new S2Polygon(std::move(loops), S2Debug::DISABLE);
    //
    // This is equivalent to:
    //
    //   S2Polygon* polygon = new S2Polygon;
    //   polygon->set_s2debug_override(S2Debug::DISABLE);
    //   polygon->InitNested(std::move(loops));
    constructor(loops: List<S2Loop>, debugOverride: S2Debug = S2Debug.ALLOW): this() {
        this.debugOverride = debugOverride
        initNested(loops)
    }

    // Convenience constructor that creates a polygon with a single loop
    // corresponding to the given cell.
    constructor(cell: S2Cell, debugOverride: S2Debug = S2Debug.ALLOW): this() {
        this.debugOverride = debugOverride
        init(S2Loop(cell))
    }

    // Convenience constructor that calls Init(S2Loop*).  Note that this method
    // automatically converts the special empty loop (see S2Loop) into an empty
    // polygon, unlike the vector-of-loops constructor which does not allow
    // empty loops at all.
    constructor(loop: S2Loop, debugOverride: S2Debug = S2Debug.ALLOW): this() {
        this.debugOverride = debugOverride
        init(loop)
    }

    fun loops(): List<S2Loop> = loops.toList()

    // Create a polygon from a set of hierarchically nested loops.  The polygon
    // interior consists of the points contained by an odd number of loops.
    // (Recall that a loop contains the set of points on its left-hand side.)
    //
    // This method figures out the loop nesting hierarchy and assigns every
    // loop a depth.  Shells have even depths, and holes have odd depths.  Note
    // that the loops are reordered so the hierarchy can be traversed more
    // easily (see GetParent(), GetLastDescendant(), and S2Loop::depth()).
    //
    // This method may be called more than once, in which case any existing
    // loops are deleted before being replaced by the input loops.
    fun initNested(loops: List<S2Loop>) {
        clearLoops()
        this.loops.addAll(loops)

        logger.trace { """Init nested loops:
            |${loops.mapIndexed { idx, l -> "$idx: ${l.toDebugString()}" }.joinToString("\n")}
        """.trimMargin() }

        if (numLoops() == 1) {
            initOneLoop()
            return
        }
        val loopMap: LoopMap = mutableMapOf()
        for (i in 0 until numLoops()) {
            insertLoop(loop(i), loopMap)
        }

        logger.trace { """Nested loops:
            |${loopMap.entries.joinToString("\n") { entry -> "Parent: ${entry.key?.toDebugString()}\nChildren: [\n${entry.value.joinToString(",\n") { c -> "   " + c.toDebugString() }} \n]" }}
        """.trimMargin() }

        // Reorder the loops in depth-first traversal order.
        // Loops are now owned by loop_map, don't let them be
        // deleted by clear().
        //for (auto& loop : loops_) loop.release();
        this.loops.clear()
        initLoops(loopMap)

        // Compute num_vertices_, bound_, subregion_bound_.
        initLoopProperties()
    }

    // Like InitNested(), but expects loops to be oriented such that the polygon
    // interior is on the left-hand side of all loops.  This implies that shells
    // and holes should have opposite orientations in the input to this method.
    // (During initialization, loops representing holes will automatically be
    // inverted.)
    fun initOriented(loops: List<S2Loop>) {
        // Here is the algorithm:
        //
        // 1. Remember which of the given loops contain S2::Origin().
        //
        // 2. Invert loops as necessary to ensure that they are nestable (i.e., no
        //    loop contains the complement of any other loop).  This may result in a
        //    set of loops corresponding to the complement of the given polygon, but
        //    we will fix that problem later.
        //
        //    We make the loops nestable by first normalizing all the loops (i.e.,
        //    inverting any loops whose curvature is negative).  This handles
        //    all loops except those whose curvature is very close to zero
        //    (within the maximum error tolerance).  Any such loops are inverted if
        //    and only if they contain S2::Origin().  (In theory this step is only
        //    necessary if there are at least two such loops.)  The resulting set of
        //    loops is guaranteed to be nestable.
        //
        // 3. Build the polygon.  This yields either the desired polygon or its
        //    complement.
        //
        // 4. If there is at least one loop, we find a loop L that is adjacent to
        //    S2::Origin() (where "adjacent" means that there exists a path
        //    connecting S2::Origin() to some vertex of L such that the path does
        //    not cross any loop).  There may be a single such adjacent loop, or
        //    there may be several (in which case they should all have the same
        //    contains_origin() value).  We choose L to be the loop containing the
        //    origin whose depth is greatest, or loop(0) (a top-level shell) if no
        //    such loop exists.
        //
        // 5. If (L originally contained origin) != (polygon contains origin), we
        //    invert the polygon.  This is done by inverting a top-level shell whose
        //    curvature is minimal and then fixing the nesting hierarchy.  Note
        //    that because we normalized all the loops initially, this step is only
        //    necessary if the polygon requires at least one non-normalized loop to
        //    represent it.

        val containedOrigin = mutableSetOf<S2Loop>()
        loops.forEach { loop ->
            if (loop.containsOrigin()) {
                containedOrigin.add(loop)
            }
            val angle = loop.curvature
            if (abs(angle) > loop.getCurvatureMaxError()) {
                // Normalize the loop.
                if (angle < 0) loop.invert()
            } else {
                // Ensure that the loop does not contain the origin.
                if (loop.containsOrigin()) loop.invert()
            }
        }
        initNested(loops)
        if (numLoops() > 0) {
            var originLoop = loop(0)
            var polygonContainsOrigin = false
            for (i in 0 until numLoops()) {
                if (loop(i).containsOrigin()) {
                    polygonContainsOrigin = polygonContainsOrigin xor true
                    originLoop = loop(i)
                }
            }
            if (containedOrigin.contains(originLoop) != polygonContainsOrigin) {
                invert()
            }
        }
        // Verify that the original loops had consistent shell/hole orientations.
        // Each original loop L should have been inverted if and only if it now
        // represents a hole.
        this.loops.indices.forEach { i ->
            val containsOrigin = loop(i).containsOrigin();
            if ((containedOrigin.contains(loop(i)) != containsOrigin) != loop(i).isHole()) {
                // There is no point in saving the loop index, because the error is a
                // property of the entire set of loops.  In general there is no way to
                // determine which ones are incorrect.
                errorInconsistentLoopOrientations = true
                if (PreConditions.enabled && debugOverride == S2Debug.ALLOW) {
                    // Note that FLAGS_s2debug is false in optimized builds (by default).
                    checkState { isValid() }
                }
            }
        }

        logger.trace { "initOriented | errorInconsistentLoopOrientations = $errorInconsistentLoopOrientations" }
    }

    // Initialize a polygon from a single loop.  Note that this method
    // automatically converts the special empty loop (see S2Loop) into an empty
    // polygon, unlike the vector-of-loops InitNested() method which does not
    // allow empty loops at all.
    fun init(loop: S2Loop) {
        logger.trace { "init | loop = ${S2TextParser.toString(loop)}" }
        // We don't allow empty loops in the other Init() methods because deleting
        // them changes the number of loops, which is awkward to handle.
        clearLoops()
        if (loop.isEmpty()) {
            initLoopProperties()
        } else {
            loops.add(loop)
            initOneLoop()
        }
    }

    // Makes a deep copy of the given source polygon.  The destination polygon
    // will be cleared if necessary.
    fun copy(src: S2Polygon) {
        clearLoops()
        for (i in 0 until src.numLoops()) {
            loops.add(src.loop(i).clone())
        }
        // Don't copy error_inconsistent_loop_orientations_, since this is not a
        // property of the polygon but only of the way the polygon was constructed.
        numVertices = src.numVertices()
        unindexedContainsCalls.set(0)
        bound = src.bound
        subregionBound = src.subregionBound
        initIndex();  // TODO(ericv): Copy the index efficiently.
    }

    // Returns true if this is a valid polygon (including checking whether all
    // the loops are themselves valid).  Note that validity is checked
    // automatically during initialization when --s2debug is enabled (true by
    // default in debug binaries).
    fun isValid(): Boolean {
        val error = findValidationError()
        if (!error.isOk()) {
            logger.error { error }
            return false;
        }
        return true;
    }

    // Returns true if this is *not* a valid polygon and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // Note that in error messages, loops that represent holes have their edges
    // numbered in reverse order, starting from the last vertex of the loop.
    //

    fun findValidationError(): S2Error {
        val error = S2Error()
        findValidationError(error)
        return error
    }

    fun findValidationError(error: S2Error): Boolean {
        for (i in 0 until numLoops()) {
            val loop = loop(i)
            logger.trace { "FindValidationError: Validate loop $i = ${S2TextParser.toString(loop)}" }
            // Check for loop errors that don't require building an S2ShapeIndex.
            if (loop.findValidationErrorNoIndex(error)) {
                error.text = "Loop $i: ${error.text}"
                return true
            }

            // Check that no loop is empty, and that the full loop only appears in the
            // full polygon.
            if (loop.isEmpty()) {
                error.init(code = S2Error.POLYGON_EMPTY_LOOP, text = "Loop $i: empty loops are not allowed")
                return true
            }
            if (loop.isFull() && numLoops() > 1) {
                error.init(
                        code = S2Error.POLYGON_EXCESS_FULL_LOOP,
                        text = "Loop $i: full loop appears in non-full polygon"
                );
                return true
            }
        }

        // Check for loop self-intersections and loop pairs that cross
        // (including duplicate edges and vertices).
        if (S2CrossingEdgePairsScanner.findSelfIntersection(index, error)) {
            return true
        }

        // Check whether InitOriented detected inconsistent loop orientations.
        if (errorInconsistentLoopOrientations) {
            error.init(code = S2Error.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS, text = "Inconsistent loop orientations detected")
            return true
        }

        // Finally, verify the loop nesting hierarchy.
        return !error.init(findLoopNestingError()).isOk()
    }

    // Return true if this is the empty polygon (consisting of no loops).
    fun isEmpty(): Boolean = loops.isEmpty()

    // Return true if this is the full polygon (consisting of a single loop that
    // encompasses the entire sphere).
    fun isFull(): Boolean = numLoops() == 1 && loop(0).isFull()

    // Return the number of loops in this polygon.
    fun numLoops(): Int = loops.size

    // Total number of vertices in all loops.
    fun numVertices(): Int = numVertices

    // Return the loop at the given index.  Note that during initialization, the
    // given loops are reordered according to a preorder traversal of the loop
    // nesting hierarchy.  This implies that every loop is immediately followed
    // by its descendants.  This hierarchy can be traversed using the methods
    // GetParent(), GetLastDescendant(), and S2Loop::depth().
    fun loop(k: Int): S2Loop = loops[k]

    // Return the index of the parent of loop k, or -1 if it has no parent.
    fun getParent(k: Int): Int {
        val depth = loop(k).depth
        if (depth == 0) return -1;  // Optimization.
        var i = k
        while (--i >= 0 && loop(i).depth >= depth) continue
        return i
    }

    // Return the index of the last loop that is contained within loop k.
    // Returns num_loops() - 1 if k < 0.  Note that loops are indexed according
    // to a preorder traversal of the nesting hierarchy, so the immediate
    // children of loop k can be found by iterating over loops
    // (k+1)..GetLastDescendant(k) and selecting those whose depth is equal to
    // (loop(k)->depth() + 1).
    fun getLastDescendant(k: Int): Int {
        if (k < 0) return numLoops() - 1
        val depth = loop(k).depth
        var i = k
        while (++i < numLoops() && loop(i).depth > depth) continue
        return i - 1
    }

    // Return the area of the polygon interior, i.e. the region on the left side
    // of an odd number of loops.  The return value is between 0 and 4*Pi.
    fun getArea(): Double = loops.map { l -> l.sign() * l.area }.sum()

    // Return the true centroid of the polygon multiplied by the area of the
    // polygon (see s2centroids.h for details on centroids).  The result is not
    // unit length, so you may want to normalize it.  Also note that in general,
    // the centroid may not be contained by the polygon.
    //
    // We prescale by the polygon area for two reasons: (1) it is cheaper to
    // compute this way, and (2) it makes it easier to compute the centroid of
    // more complicated shapes (by splitting them into disjoint regions and
    // adding their centroids).
    fun getCentroid(): S2Point =
            loops.map { l -> l.centroid * l.sign().toDouble() }.fold(S2Point()) { acc, point -> acc + point }

    // If all of the polygon's vertices happen to be the centers of S2Cells at
    // some level, then return that level, otherwise return -1.  See also
    // InitToSnapped() and s2builderutil::S2CellIdSnapFunction.
    // Returns -1 if the polygon has no vertices.
    fun getSnapLevel(): Int {
        var snapLevel = -1
        loops.forEach { child ->
            for (j in 0 until child.numVertices) {
                val (level, _) = S2Coords.xyzToFaceSiTi(child.vertex(j))
                if (level < 0) return level;  // Vertex is not a cell center.
                if (level != snapLevel) {
                    if (snapLevel < 0) {
                        snapLevel = level  // First vertex.
                    } else {
                        return -1;  // Vertices at more than one cell level.
                    }
                }
            }
        }
        return snapLevel
    }

    // Return the distance from the given point to the polygon interior.  If the
    // polygon is empty, return S1Angle::Infinity().  "x" should be unit length.
    fun getDistance(x: S2Point): S1Angle {
        // Note that S2Polygon::Contains(S2Point) is slightly more efficient than
        // the generic version used by S2ClosestEdgeQuery.
        if (contains(x)) return S1Angle.zero()
        return getDistanceToBoundary(x)
    }

    // Return the distance from the given point to the polygon boundary.  If the
    // polygon is empty or full, return S1Angle::Infinity() (since the polygon
    // has no boundary).  "x" should be unit length.
    fun getDistanceToBoundary(x: S2Point): S1Angle {
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = false
        val t = S2ClosestEdgeQuery.PointTarget(x)
        return S2ClosestEdgeQuery(index, options).getDistance(t).toAngle()
    }

    // If the given point is contained by the polygon, return it.  Otherwise
    // return the closest point on the polygon boundary.  If the polygon is
    // empty, return the input argument.  Note that the result may or may not be
    // contained by the polygon.  "x" should be unit length.
    fun project(x: S2Point): S2Point {
        if (contains(x)) return x
        return projectToBoundary(x)
    }

    // Return the closest point on the polygon boundary to the given point.  If
    // the polygon is empty or full, return the input argument (since the
    // polygon has no boundary).  "x" should be unit length.
    fun projectToBoundary(x: S2Point): S2Point {
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = false
        val q = S2ClosestEdgeQuery(index, options)
        val target = S2ClosestEdgeQuery.PointTarget(x)
        val edge = q.findClosestEdge(target)
        return q.project(x, edge)
    }

    // Return true if this polygon contains the given other polygon, i.e.
    // if polygon A contains all points contained by polygon B.
    fun contains(b: S2Polygon): Boolean {
        // It's worth checking bounding rectangles, since they are precomputed.
        // Note that the first bound has been expanded to account for possible
        // numerical errors in the second bound.
        if (!subregionBound.contains(b.bound)) {
            // It is possible that A contains B even though Bound(A) does not contain
            // Bound(B).  This can only happen when polygon B has at least two outer
            // shells and the union of the two bounds spans all longitudes.  For
            // example, suppose that B consists of two shells with a longitude gap
            // between them, while A consists of one shell that surrounds both shells
            // of B but goes the other way around the sphere (so that it does not
            // intersect the longitude gap).
            //
            // For convenience we just check whether B has at least two loops rather
            // than two outer shells.
            if (b.numLoops() == 1 || !bound.lng.union(b.bound.lng).isFull) {
                return false;
            }
        }

        // The following case is not handled by S2BooleanOperation because it only
        // determines whether the boundary of the result is empty (which does not
        // distinguish between the full and empty polygons).
        if (isEmpty() && b.isFull()) return false

        return S2BooleanOperation.contains(index, b.index)
    }

    // Returns true if this polgyon (A) approximately contains the given other
    // polygon (B). This is true if it is possible to move the vertices of B
    // no further than "tolerance" such that A contains the modified B.
    //
    // For example, the empty polygon will contain any polygon whose maximum
    // width is no more than "tolerance".
    fun approxContains(b: S2Polygon, tolerance: S1Angle): Boolean {
        val difference = S2Polygon()
        difference.initToApproxDifference(b, this, tolerance)
        return difference.isEmpty()
    }

    // Return true if this polygon intersects the given other polygon, i.e.
    // if there is a point that is contained by both polygons.
    fun intersects(b: S2Polygon): Boolean {
        // It's worth checking bounding rectangles, since they are precomputed.
        if (!bound.intersects(b.bound)) return false

        // The following case is not handled by S2BooleanOperation because it only
        // determines whether the boundary of the result is empty (which does not
        // distinguish between the full and empty polygons).
        if (isFull() && b.isFull()) return true

        return S2BooleanOperation.intersects(index, b.index)
    }

    // Returns true if this polgyon (A) and the given polygon (B) are
    // approximately disjoint.  This is true if it is possible to ensure that A
    // and B do not intersect by moving their vertices no further than
    // "tolerance".  This implies that in borderline cases where A and B overlap
    // slightly, this method returns true (A and B are approximately disjoint).
    //
    // For example, any polygon is approximately disjoint from a polygon whose
    // maximum width is no more than "tolerance".
    fun approxDisjoint(b: S2Polygon, tolerance: S1Angle): Boolean {
        val intersection = S2Polygon()
        intersection.initToApproxIntersection(b, this, tolerance)
        return intersection.isEmpty()
    }

    // Invert the polygon (replace it by its complement).
    fun invert(): Unit {
        // Inverting any one loop will invert the polygon.  The best loop to invert
        // is the one whose area is largest, since this yields the smallest area
        // after inversion.  The loop with the largest area is always at depth 0.
        // The descendents of this loop all have their depth reduced by 1, while the
        // former siblings of this loop all have their depth increased by 1.

        // The empty and full polygons are handled specially.
        if (isEmpty()) {
            loops.add(S2Loop(S2Loop.kFull))
        } else if (isFull()) {
            clearLoops()
        } else {
            // Find the loop whose area is largest (i.e., whose curvature is
            // smallest), minimizing calls to GetCurvature().  In particular, for
            // polygons with a single shell at level 0 there is not need to call
            // GetCurvature() at all.  (This method is relatively expensive.)
            var best = 0
            val kNone = 10.0;  // Flag that means "not computed yet"
            var bestAngle = kNone
            for (i in 1 until numLoops()) {
                if (loop(i).depth == 0) {
                    // We defer computing the curvature of loop 0 until we discover
                    // that the polygon has another top-level shell.
                    if (bestAngle == kNone) bestAngle = loop(best).curvature
                    val angle = loop(i).curvature
                    // We break ties deterministically in order to avoid having the output
                    // depend on the input order of the loops.
                    if (angle < bestAngle ||
                            (angle == bestAngle && compareLoops(loop(i), loop(best)) < 0)) {
                        best = i
                        bestAngle = angle
                    }
                }
            }
            // Build the new loops vector, starting with the inverted loop.
            loop(best).invert()
            val newLoops = mutableListOf<S2Loop>()
            // Add the former siblings of this loop as descendants.
            val lastBest = getLastDescendant(best)
            newLoops.add(loops[best])
            for (i in 0 until numLoops()) {
                if (i < best || i > lastBest) {
                    loop(i).depth = loop(i).depth + 1
                    newLoops.add(loops[i])
                }
            }
            // Add the former children of this loop as siblings.
            for (i in 0 until numLoops()) {
                if (i in (best + 1)..lastBest) {
                    loop(i).depth = loop(i).depth - 1
                    newLoops.add(loops[i])
                }
            }
            loops.clear()
            loops.addAll(newLoops)
            checkEQ(newLoops.size, numLoops())
        }
        clearIndex()
        initLoopProperties()
    }

    // Initialize this polygon to the complement of the given polygon.
    fun initToComplement(a: S2Polygon) {
        copy(a)
        invert()
    }

    // Initialize this polygon to the intersection, union, difference (A - B),
    // or symmetric difference (XOR) of the given two polygons.
    //
    // "snap_function" allows you to specify a minimum spacing between output
    // vertices, and/or that the vertices should be snapped to a discrete set of
    // points (e.g. S2CellId centers or E7 lat/lng coordinates).  Any snap
    // function can be used, including the IdentitySnapFunction with a
    // snap_radius of zero (which preserves the input vertices exactly).
    //
    // The boundary of the output polygon before snapping is guaranteed to be
    // accurate to within S2::kIntersectionError of the exact result.
    // Snapping can move the boundary by an additional distance that depends on
    // the snap function.  Finally, any degenerate portions of the output
    // polygon are automatically removed (i.e., regions that do not contain any
    // points) since S2Polygon does not allow such regions.
    //
    // See S2Builder and s2builderutil for more details on snap functions.  For
    // example, you can snap to E7 coordinates by setting "snap_function" to
    // s2builderutil::IntLatLngSnapFunction(7).
    //
    // The default snap function is the IdentitySnapFunction with a snap radius
    // of S2::kIntersectionMergeRadius (equal to about 1.8e-15 radians
    // or 11 nanometers on the Earth's surface).  This means that vertices may
    // be positioned arbitrarily, but vertices that are extremely close together
    // can be merged together.  The reason for a non-zero default snap radius is
    // that it helps to eliminate narrow cracks and slivers when T-vertices are
    // present.  For example, adjacent S2Cells at different levels do not share
    // exactly the same boundary, so there can be a narrow crack between them.
    // If a polygon is intersected with those cells and the pieces are unioned
    // together, the result would have a narrow crack unless the snap radius is
    // set to a non-zero value.
    //
    // Note that if you want to encode the vertices in a lower-precision
    // representation (such as S2CellIds or E7), it is much better to use a
    // suitable SnapFunction rather than rounding the vertices yourself, because
    // this will create self-intersections unless you ensure that the vertices
    // and edges are sufficiently well-separated first.  In particular you need
    // to use a snap function whose min_edge_vertex_separation() is at least
    // twice the maximum distance that a vertex can move when rounded.
    //
    // The versions of these functions with an S2Error argument return true on
    // success and set "error" appropriately otherwise.  However note that these
    // functions should never return an error provided that both input polygons
    // are valid (i.e., IsValid() returns true).
    fun initToIntersection(a: S2Polygon, b: S2Polygon) {
        initToApproxIntersection(a, b, S2EdgeCrossings.kIntersectionMergeRadius)
    }

    fun initToIntersection(a: S2Polygon, b: S2Polygon, snapFunction: SnapFunction, error: S2Error = S2Error()): Boolean {
        if (!a.bound.intersects(b.bound)) return true
        return initToOperation(S2BooleanOperation.OpType.INTERSECTION, snapFunction, a, b, error)
    }

    fun initToUnion(a: S2Polygon, b: S2Polygon) {
        initToApproxUnion(a, b, S2EdgeCrossings.kIntersectionMergeRadius)
    }

    fun initToUnion(a: S2Polygon, b: S2Polygon, snapFunction: SnapFunction, error: S2Error = S2Error()): Boolean {
        return initToOperation(S2BooleanOperation.OpType.UNION, snapFunction, a, b, error)
    }

    fun initToDifference(a: S2Polygon, b: S2Polygon) {
        initToApproxDifference(a, b, S2EdgeCrossings.kIntersectionMergeRadius)
    }

    fun initToDifference(a: S2Polygon, b: S2Polygon, snapFunction: SnapFunction, error: S2Error = S2Error()): Boolean {
        return initToOperation(S2BooleanOperation.OpType.DIFFERENCE, snapFunction, a, b, error)
    }

    fun initToSymmetricDifference(a: S2Polygon, b: S2Polygon) {
        initToApproxSymmetricDifference(a, b, S2EdgeCrossings.kIntersectionMergeRadius)
    }

    fun initToSymmetricDifference(a: S2Polygon, b: S2Polygon, snapFunction: SnapFunction, error: S2Error = S2Error()): Boolean {
        return initToOperation(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snapFunction, a, b, error)
    }

    // Convenience functions that use the IdentitySnapFunction with the given
    // snap radius.  TODO(ericv): Consider deprecating these and require the
    // snap function to be specified explcitly?
    fun initToApproxIntersection(a: S2Polygon, b: S2Polygon, snapRadius: S1Angle) {
        initToIntersection(a, b, IdentitySnapFunction(snapRadius))
    }

    fun initToApproxUnion(a: S2Polygon, b: S2Polygon, snapRadius: S1Angle) {
        initToUnion(a, b, IdentitySnapFunction(snapRadius))
    }

    fun initToApproxDifference(a: S2Polygon, b: S2Polygon, snapRadius: S1Angle) {
        initToDifference(a, b, IdentitySnapFunction(snapRadius))
    }

    fun initToApproxSymmetricDifference(a: S2Polygon, b: S2Polygon, snapRadius: S1Angle) {
        initToSymmetricDifference(a, b, IdentitySnapFunction(snapRadius))
    }

    /**
     * Initializes the polygon from input polygon "a" using the given S2Builder. If the result has an empty boundary
     * (no loops), also decides whether the result should be the full polygon rather than the empty one based on the
     * area of the input polygon.
     *
     * @param a A polygon
     * @param builder A builder.
     *
     * @see initToApproxIntersection
     */
    fun initFromBuilder(a: S2Polygon, builder: S2Builder) {
        builder.startLayer(S2PolygonLayer(this))
        builder.addPolygon(a)
        val error = S2Error()
        if (!builder.build(error)) {
            logger.error { "Could not build polygon: $error" }
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  (See InitToApproxIntersection.)
        if (numLoops() == 0) {
            if (a.bound.area > 2 * M_PI && a.getArea() > 2 * M_PI) invert()
        }
    }

    fun operationWithPolyline(opType: S2BooleanOperation.OpType, snapFunction: SnapFunction, a: S2Polyline): List<S2Polyline> {
        val options = S2BooleanOperation.Options()
        options.snapFunction = snapFunction
        val result = mutableListOf<S2Polyline>()
        val layerOptions = S2PolylineVectorLayer.Options()
        layerOptions.polylineType = Graph.PolylineType.WALK
        val op = S2BooleanOperation(opType, S2PolylineVectorLayer(result, options = layerOptions), options)
        val aIndex = MutableS2ShapeIndex()
        aIndex.add(S2Polyline.Shape(polyline = a))
        val error = S2Error()
        if (!op.build(aIndex, index, error = error)) {
            logger.error { "Polyline $opType operation failed: $error" }
        }
        return result
    }

    // Snaps the vertices of the given polygon using the given SnapFunction
    // (e.g., s2builderutil::IntLatLngSnapFunction(6) snaps to E6 coordinates).
    // This can change the polygon topology (merging loops, for example), but
    // the resulting polygon is guaranteed to be valid, and no vertex will move
    // by more than snap_function.snap_radius().  See S2Builder for other
    // guarantees (e.g., minimum edge-vertex separation).
    //
    // Note that this method is a thin wrapper over S2Builder, so if you are
    // starting with data that is not in S2Polygon format (e.g., integer E7
    // coordinates) then it is faster to just use S2Builder directly.
    fun initToSnapped(polygon: S2Polygon, snapFunction: SnapFunction) {
        val builder = S2Builder(S2Builder.Options(snapFunction = snapFunction))
        initFromBuilder(polygon, builder)
    }

    // Convenience function that snaps the vertices to S2CellId centers at the
    // given level (default level 30, which has S2CellId centers spaced about 1
    // centimeter apart).  Polygons can be efficiently encoded by Encode() after
    // they have been snapped.
    fun initToSnapped(polygon: S2Polygon, snapLevel: Int = S2CellId.kMaxLevel) {
        val builder = S2Builder(S2Builder.Options(snapFunction = S2CellIdSnapFunction(snapLevel), verbose = false))
        initFromBuilder(polygon, builder)
    }

    // Snaps the input polygon according to the given "snap_function" and
    // reduces the number of vertices if possible, while ensuring that no vertex
    // moves further than snap_function.snap_radius().
    //
    // Simplification works by replacing nearly straight chains of short edges
    // with longer edges, in a way that preserves the topology of the input
    // polygon up to the creation of degeneracies.  This means that loops or
    // portions of loops may become degenerate, in which case they are removed.
    // For example, if there is a very small island in the original polygon, it
    // may disappear completely.  (Even if there are dense islands, they could
    // all be removed rather than being replaced by a larger simplified island
    // if more area is covered by water than land.)
    fun initToSimplified(a: S2Polygon, snap_function: SnapFunction) {
        val builder = S2Builder(S2Builder.Options(snapFunction = snap_function, simplifyEdgeChains = true))
        initFromBuilder(a, builder)
    }

    fun toDebugString(loop_separator: String = " ; "): String {
        if (this.isEmpty()) {
            return "empty"
        } else if (this.isFull()) {
            return "full"
        }
        var out = ""
        for (i in 0 until this.numLoops()) {
            if (i > 0) out += loop_separator
            val loop = this.loop(i)
            out += loop.vertices().joinToString(", ") { p -> S2PointUtil.toDegreesString(p) }
        }
        return out
    }

    // Given a point "p" inside an S2Cell or on its boundary, return a mask
    // indicating which of the S2Cell edges the point lies on.  All boundary
    // comparisons are to within a maximum "u" or "v" error of "tolerance_uv".
    // Bit "i" in the result is set if and only "p" is incident to the edge
    // corresponding to S2Cell::edge(i).
    private fun getCellEdgeIncidenceMask(cell: S2Cell, p: S2Point, tolerance_uv: Double, debugOverride: S2Debug = S2Debug.ALLOW): UByte {
        var mask: UByte = 0.toUByte()
        val uv = S2Coords.faceXyztoUv(cell.face(), p)
        if (uv != null) {
            val bound = cell.boundUV()
            if (PreConditions.enabled && debugOverride == S2Debug.ALLOW) { checkState { bound.expanded(tolerance_uv).contains(uv) } }
            if (abs(uv[1] - bound[1][0]) <= tolerance_uv) mask = mask or 1.toUByte()
            if (abs(uv[0] - bound[0][1]) <= tolerance_uv) mask = mask or 2.toUByte()
            if (abs(uv[1] - bound[1][1]) <= tolerance_uv) mask = mask or 4.toUByte()
            if (abs(uv[0] - bound[0][0]) <= tolerance_uv) mask = mask or 8.toUByte()
        }
        return mask;
    }

    // Like InitToSimplified, except that any vertices or edges on the boundary
    // of the given S2Cell are preserved if possible.  This method requires that
    // the polygon has already been clipped so that it does not extend outside
    // the cell by more than "boundary_tolerance".  In other words, it operates
    // on polygons that have already been intersected with a cell.
    //
    // Typically this method is used in geometry-processing pipelines that
    // intersect polygons with a collection of S2Cells and then process those
    // cells in parallel, where each cell generates some geometry that needs to
    // be simplified.  In contrast, if you just need to simplify the *input*
    // geometry then it is easier and faster to do the simplification before
    // computing the intersection with any S2Cells.
    //
    // "boundary_tolerance" specifies how close a vertex must be to the cell
    // boundary to be kept.  The default tolerance is large enough to handle any
    // reasonable way of interpolating points along the cell boundary, such as
    // S2::GetIntersection(), S2::Interpolate(), or direct (u,v)
    // interpolation using S2::FaceUVtoXYZ().  However, if the vertices have
    // been snapped to a lower-precision representation (e.g., S2CellId centers
    // or E7 coordinates) then you will need to set this tolerance explicitly.
    // For example, if the vertices were snapped to E7 coordinates then
    // "boundary_tolerance" should be set to
    //
    //   s2builderutil::IntLatLngSnapFunction::MinSnapRadiusForExponent(7)
    //
    // Degenerate portions of loops are always removed, so if a vertex on the
    // cell boundary belongs only to degenerate regions then it will not be
    // kept.  For example, if the input polygon is a narrow strip of width less
    // than "snap_radius" along one side of the cell, then the entire loop may
    // become degenerate and be removed.
    //
    // REQUIRES: all vertices of "a" are within "boundary_tolerance" of "cell".
    fun initToSimplifiedInCell(a: S2Polygon, cell: S2Cell, snap_radius: S1Angle, boundary_tolerance: S1Angle = S1Angle.radians(1e-15)) {
        // The polygon to be simplified consists of "boundary edges" that follow the
        // cell boundary and "interior edges" that do not.  We want to simplify the
        // interior edges while leaving the boundary edges unchanged.  It's not
        // sufficient to call S2Builder::ForceVertex() on all boundary vertices.
        // For example, suppose the polygon includes a triangle ABC where all three
        // vertices are on the cell boundary and B is a cell corner.  Then if
        // interior edge AC snaps to vertex B, this loop would become degenerate and
        // be removed.  Similarly, we don't want boundary edges to snap to interior
        // vertices, since this also would cause portions of the polygon along the
        // boundary to be removed.
        //
        // Instead we use a two-pass algorithm.  In the first pass, we simplify
        // *only* the interior edges, using ForceVertex() to ensure that any edge
        // endpoints on the cell boundary do not move.  In the second pass, we add
        // the boundary edges (which are guaranteed to still form loops with the
        // interior edges) and build the output polygon.
        //
        // Note that in theory, simplifying the interior edges could create an
        // intersection with one of the boundary edges, since if two interior edges
        // intersect very near the boundary then the intersection point could be
        // slightly outside the cell (by at most S2::kIntersectionError).
        // This is the *only* way that a self-intersection can be created, and it is
        // expected to be extremely rare.  Nevertheless we use a small snap radius
        // in the second pass in order to eliminate any such self-intersections.
        //
        // We also want to preserve the cyclic vertex order of loops, so that the
        // original polygon can be reconstructed when no simplification is possible
        // (i.e., idempotency).  In order to do this, we group the input edges into
        // a sequence of polylines.  Each polyline contains only one type of edge
        // (interior or boundary).  We use S2Builder to simplify the interior
        // polylines, while the boundary polylines are passed through unchanged.
        // Each interior polyline is in its own S2Builder layer in order to keep the
        // edges in sequence.  This lets us ensure that in the second pass, the
        // edges are added in their original order so that S2PolygonLayer can
        // reconstruct the original loops.

        // We want an upper bound on how much "u" or "v" can change when a point on
        // the boundary of the S2Cell is moved away by up to "boundary_tolerance".
        // Inverting this, instead we could compute a lower bound on how far a point
        // can move away from an S2Cell edge when "u" or "v" is changed by a given
        // amount.  The latter quantity is simply (S2::kMinWidth.deriv() / 2)
        // under the S2_LINEAR_PROJECTION model, where we divide by 2 because we
        // want the bound in terms of (u = 2 * s - 1) rather than "s" itself.
        // Consulting s2metrics.cc, this value is sqrt(2/3)/2 = sqrt(1/6).
        // Going back to the original problem, this gives:
        val boundary_tolerance_uv = boundary_tolerance.radians * sqrt(6.0)

        // The first pass yields a collection of simplified polylines that preserve
        // the original cyclic vertex order.
        val polylines = simplifyEdgesInCell(a, cell, boundary_tolerance_uv, snap_radius)

        // The second pass eliminates any intersections between interior edges and
        // boundary edges, and then assembles the edges into a polygon.
        val builder = S2Builder(S2Builder.Options(
                snapFunction = IdentitySnapFunction(S2EdgeCrossings.kIntersectionError),
                idempotent = false  // Force snapping up to the given radius
        ))
        builder.startLayer(S2PolygonLayer(this))
        for (polyline in polylines) {
            builder.addPolyline(polyline)
        }
        val error = S2Error()
        if (!builder.build(error)) {
            logger.error { "Could not build polygon: $error" }
            return;
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  (See InitToApproxIntersection.)
        if (numLoops() == 0) {
            if (a.bound.area > 2 * M_PI && a.getArea() > 2 * M_PI) invert()
        }
    }

    // Return true if this polygon contains the given polyline.  This method
    // returns an exact result, according to the following model:
    //
    //  - All edges are geodesics (of course).
    //
    //  - Vertices are ignored for the purposes of defining containment.
    //    (This is because polygons often do not contain their vertices, in
    //    order to that when a set of polygons tiles the sphere then every point
    //    is contained by exactly one polygon.)
    //
    //  - Points that lie exactly on geodesic edges are resolved using symbolic
    //    perturbations (i.e., they are considered to be infinitesmally offset
    //    from the edge).
    //
    //  - If the polygon and polyline share an edge, it is handled as follows.
    //    First, the polygon edges are oriented so that the interior is always
    //    on the left.  Then the shared polyline edge is contained if and only
    //    if it is in the same direction as the corresponding polygon edge.
    //    (This model ensures that when a polyline is intersected with a polygon
    //    and its complement, the edge only appears in one of the two results.)
    //
    // TODO(ericv): Update the implementation to correspond to the model above.
    fun contains(b: S2Polyline): Boolean = approxContains(b, S2EdgeCrossings.kIntersectionMergeRadius)

    // Returns true if this polgyon approximately contains the given polyline
    // This is true if it is possible to move the polyline vertices no further
    // than "tolerance" such that the polyline is now contained.
    fun approxContains(b: S2Polyline, tolerance: S1Angle): Boolean {
        val difference = approxSubtractFromPolyline(b, tolerance);
        return difference.isEmpty()
    }

    // Return true if this polygon intersects the given polyline.  This method
    // returns an exact result; see Contains(S2Polyline) for details.
    fun intersects(b: S2Polyline): Boolean = !approxDisjoint(b, S2EdgeCrossings.kIntersectionMergeRadius)

    // Returns true if this polgyon is approximately disjoint from the given
    // polyline.  This is true if it is possible to avoid intersection by moving
    // their vertices no further than "tolerance".
    //
    // This implies that in borderline cases where there is a small overlap,
    // this method returns true (i.e., they are approximately disjoint).
    fun approxDisjoint(b: S2Polyline, tolerance: S1Angle): Boolean {
        val intersection = approxIntersectWithPolyline(b, tolerance)
        return intersection.isEmpty()
    }

    // Intersect this polygon with the polyline "in" and return the resulting
    // zero or more polylines.  The polylines are returned in the order they
    // would be encountered by traversing "in" from beginning to end.
    // Note that the output may include polylines with only one vertex,
    // but there will not be any zero-vertex polylines.
    //
    // This is equivalent to calling ApproxIntersectWithPolyline() with the
    // "snap_radius" set to S2::kIntersectionMergeRadius.
    fun intersectWithPolyline(line: S2Polyline): List<S2Polyline> =
            approxIntersectWithPolyline(line, S2EdgeCrossings.kIntersectionMergeRadius)

    // Similar to IntersectWithPolyline(), except that vertices will be
    // dropped as necessary to ensure that all adjacent vertices in the
    // sequence obtained by concatenating the output polylines will be
    // farther than "snap_radius" apart.  Note that this can change
    // the number of output polylines and/or yield single-vertex polylines.
    fun approxIntersectWithPolyline(line: S2Polyline, snapRadius: S1Angle): List<S2Polyline> {
        return intersectWithPolyline(line, IdentitySnapFunction(snapRadius))
    }

    // TODO(ericv): Update documentation.
    fun intersectWithPolyline(line: S2Polyline, snapFunction: SnapFunction): List<S2Polyline> =
            operationWithPolyline(S2BooleanOperation.OpType.INTERSECTION, snapFunction, line)

    // Same as IntersectWithPolyline, but subtracts this polygon from
    // the given polyline.
    fun subtractFromPolyline(line: S2Polyline): List<S2Polyline> =
            approxSubtractFromPolyline(line, S2EdgeCrossings.kIntersectionMergeRadius)

    // Same as ApproxIntersectWithPolyline, but subtracts this polygon
    // from the given polyline.
    fun approxSubtractFromPolyline(line: S2Polyline, snapRadius: S1Angle): List<S2Polyline> =
            subtractFromPolyline(line, IdentitySnapFunction(snapRadius))

    fun subtractFromPolyline(line: S2Polyline, snapFunction: SnapFunction): List<S2Polyline> =
            operationWithPolyline(S2BooleanOperation.OpType.DIFFERENCE, snapFunction, line)

    // Initialize this polygon to the outline of the given cell union.
    // In principle this polygon should exactly contain the cell union and
    // this polygon's inverse should not intersect the cell union, but rounding
    // issues may cause this not to be the case.
    fun initToCellUnionBorder(cells: S2CellUnion) {
        /// We use S2Builder to compute the union.  Due to rounding errors, we can't
        // compute an exact union - when a small cell is adjacent to a larger cell,
        // the shared edges can fail to line up exactly.  Two cell edges cannot come
        // closer then kMinWidth, so if we have S2Builder snap edges within half
        // that distance, then we should always merge shared edges without merging
        // different edges.
        val snapRadius = 0.5 * S2Coords.projection.kMinWidth.getValue(S2CellId.kMaxLevel)
        val builder = S2Builder(options = S2Builder.Options(snapFunction = IdentitySnapFunction(S1Angle.radians(snapRadius))))
        builder.startLayer(S2PolygonLayer(this))
        for (id in cells) {
            builder.addLoop(S2Loop(S2Cell(id)))
        }
        val error = S2Error()
        if (!builder.build(error)) {
            logger.error { "InitToCellUnionBorder failed: $error" }
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  There are only two ways that this can
        // happen: either the cell union is empty, or it consists of all six faces.
        if (numLoops() == 0) {
            if (cells.isEmpty()) return
            checkState { 6UL shl (2 * S2CellId.kMaxLevel) == cells.leafCellsCovered() }
            invert()
        }
    }

    // Return true if every loop of this polygon shares at most one vertex with
    // its parent loop.  Every polygon has a unique normalized form.  A polygon
    // can be normalized by passing it through S2Builder (with no snapping) in
    // order to reconstruct the polygon from its edges.
    //
    // Generally there is no reason to convert polygons to normalized form.  It
    // is mainly useful for testing in order to compare whether two polygons
    // have exactly the same interior, even when they have a different loop
    // structure.  For example, a diamond nested within a square (touching at
    // four points) could be represented as a square with a diamond-shaped hole,
    // or as four triangles.  Methods such as BoundaryApproxEquals() will report
    // these polygons as being different (because they have different
    // boundaries) even though they contain the same points.  However if they
    // are both converted to normalized form (the "four triangles" version) then
    // they can be compared more easily.
    //
    // Also see ApproxEquals(), which can determine whether two polygons contain
    // approximately the same set of points without any need for normalization.
    fun isNormalized(): Boolean {
        // TODO(ericv): The condition tested here is insufficient.  The correct
        // condition is that each *connected component* of child loops can share at
        // most one vertex with their parent loop.  Example: suppose loop A has
        // children B, C, D, and the following pairs are connected: AB, BC, CD, DA.
        // Then the polygon is not normalized.
        val vertices = mutableSetOf<S2Point>()
        var last_parent: S2Loop? = null
        for (i in 0 until numLoops()) {
            val child = loop(i);
            if (child.depth == 0) continue
            val parent = loop(getParent(i))
            if (parent != last_parent) {
                vertices.clear();
                for (j in 0 until parent.numVertices) {
                    vertices.add(parent.vertex(j))
                }
                last_parent = parent;
            }
            var count = 0;
            for (j in 0 until child.numVertices) {
                if (vertices.contains(child.vertex(j))) ++count
            }
            if (count > 1) return false
        }
        return true
    }

    // Return true if two polygons have exactly the same loops.  The loops must
    // appear in the same order, and corresponding loops must have the same
    // linear vertex ordering (i.e., cyclic rotations are not allowed).
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as S2Polygon

        if (numLoops() != other.numLoops()) return false
        for (i in 0 until numLoops()) {
            val a_loop = loop(i)
            val b_loop = other.loop(i)
            if ((b_loop.depth != a_loop.depth) || !b_loop.equals(a_loop)) {
                return false;
            }
        }
        return true;
    }

    override fun hashCode(): Int {
        return loops.hashCode()
    }

    override fun toString(): String {
        return S2TextParser.toString(this)
    }


    // Return true if two polygons are approximately equal to within the given
    // tolerance.  This is true if it is possible to move the vertices of the
    // two polygons so that they contain the same set of points.
    //
    // Note that according to this model, small regions less than "tolerance" in
    // width do not need to be considered, since these regions can be collapsed
    // into degenerate loops (which contain no points) by moving their vertices.
    //
    // This model is not as strict as using the Hausdorff distance would be, and
    // it is also not as strict as BoundaryNear (defined below).  However, it is
    // a good choice for comparing polygons that have been snapped, simplified,
    // unioned, etc, since these operations use a model similar to this one
    // (i.e., degenerate loops or portions of loops are automatically removed).
    fun approxEquals(b: S2Polygon, tolerance: S1Angle): Boolean {
        // TODO(ericv): This can be implemented more cheaply with S2Builder, by
        // simply adding all the edges from one polygon, adding the reversed edge
        // from the other polygon, and turning on the options to split edges and
        // discard sibling pairs.  Then the polygons are approximately equal if the
        // output graph has no edges.
        val symmetric_difference = S2Polygon()
        symmetric_difference.initToApproxSymmetricDifference(b, this, tolerance)
        return symmetric_difference.isEmpty()
    }

    // Returns true if two polygons have the same boundary.  More precisely,
    // this method requires that both polygons have loops with the same cyclic
    // vertex order and the same nesting hierarchy.  (This implies that vertices
    // may be cyclically rotated between corresponding loops, and the loop
    // ordering may be different between the two polygons as long as the nesting
    // hierarchy is the same.)
    fun boundaryEquals(b: S2Polygon): Boolean {
        if (numLoops() != b.numLoops()) return false

        for (i in 0 until numLoops()) {
            val a_loop = loop(i)
            var success = false
            for (j in 0 until numLoops()) {
                val b_loop = b.loop(j)
                if ((b_loop.depth == a_loop.depth) && b_loop.boundaryEquals(a_loop)) {
                    success = true
                    break
                }
            }
            if (!success) return false
        }
        return true
    }

    // Return true if two polygons have the same boundary except for vertex
    // perturbations.  Both polygons must have loops with the same cyclic vertex
    // order and the same nesting hierarchy, but the vertex locations are
    // allowed to differ by up to "max_error".
    fun boundaryApproxEquals(b: S2Polygon, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        if (numLoops() != b.numLoops()) return false

        // For now, we assume that there is at most one candidate match for each
        // loop.  (So far this method is just used for testing.)

        for (i in 0 until numLoops()) {
            val a_loop = loop(i)
            var success = false
            for (j in 0 until numLoops()) {
                val b_loop = b.loop(j)
                if (b_loop.depth == a_loop.depth && b_loop.boundaryApproxEquals(a_loop, maxError)) {
                    success = true
                    break
                }
            }
            if (!success) return false
        }
        return true;
    }

    // Return true if two polygons have boundaries that are within "max_error"
    // of each other along their entire lengths.  More precisely, there must be
    // a bijection between the two sets of loops such that for each pair of
    // loops, "a_loop->BoundaryNear(b_loop)" is true.
    fun boundaryNear(b: S2Polygon, maxError: S1Angle = S1Angle.radians(1e-15)): Boolean {
        if (numLoops() != b.numLoops()) return false

        // For now, we assume that there is at most one candidate match for each
        // loop.  (So far this method is just used for testing.)

        for (i in 0 until numLoops()) {
            val a_loop = loop(i)
            var success = false
            for (j in 0 until numLoops()) {
                val b_loop = b.loop(j)
                if (b_loop.depth == a_loop.depth && b_loop.boundaryNear(a_loop, maxError)) {
                    success = true
                    break
                }
            }
            if (!success) return false
        }
        return true
    }

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.

    public override fun clone(): S2Polygon {
        val result = S2Polygon()
        result.copy(this)
        return result
    }

    override val capBound: S2Cap
        get() = bound.capBound

    override val rectBound: S2LatLngRect
        get() = bound

    override fun getCellUnionBound(cellIds: MutableList<S2CellId>) {
        return makeS2ShapeIndexRegion(index).getCellUnionBound(cellIds)
    }

    override fun contains(cell: S2Cell): Boolean = makeS2ShapeIndexRegion(index).contains(cell)

    override fun mayIntersect(cell: S2Cell): Boolean = makeS2ShapeIndexRegion(index).mayIntersect(cell)

    // The point 'p' does not need to be normalized.
    override fun contains(p: S2Point): Boolean {
        // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
        // worthwhile only when it might allow us to delay building the index.
        if (!index.isFresh() && !bound.contains(p)) return false

        // For small polygons it is faster to just check all the crossings.
        // Otherwise we keep track of the number of calls to Contains() and only
        // build the index once enough calls have been made so that we think it is
        // worth the effort.  See S2Loop::Contains(S2Point) for detailed comments.
        if (numVertices <= kMaxBruteForceVertices ||
                (!index.isFresh() && unindexedContainsCalls.incrementAndGet() != kMaxUnindexedContainsCalls)) {
            var inside = false
            loops.forEach { l ->
                inside = inside xor l.bruteForceContains(p)
            }
            return inside
        }
        // Otherwise we look up the S2ShapeIndex cell containing this point.
        return makeS2ContainsPointQuery(index).contains(p)
    }

    // Wrapper class for indexing a polygon (see S2ShapeIndex).  Once this
    // object is inserted into an S2ShapeIndex it is owned by that index, and
    // will be automatically deleted when no longer needed by the index.  Note
    // that this class does not take ownership of the polygon itself (see
    // OwningShape below).  You can also subtype this class to store additional
    // data (see S2Shape for details).
    //
    // Note that unlike S2Polygon, the edges of S2Polygon::Shape are directed
    // such that the polygon interior is always on the left.
    class Shape(id: Int = -1, val polygon: S2Polygon) : S2Shape(id) {
        // The total number of edges in the polygon.  This is the same as
        // polygon_->num_vertices() except in one case (polygon_->is_full()).  On
        // the other hand this field doesn't take up any extra space due to field
        // packing with S2Shape::id_.
        //
        // TODO(ericv): Consider using this field instead as an atomic<int> hint to
        // speed up edge location when there are a large number of loops.  Also
        // consider changing S2Polygon::num_vertices to num_edges instead.
        override val numEdges: Int

        // An array where element "i" is the total number of edges in loops 0..i-1.
        // This field is only used for polygons that have a large number of loops.
        private val cumulativeEdges: IntArray

        init {
            var cumulativeNumEdges = 0
            if (!polygon.isFull()) {
                val kMaxLinearSearchLoops = 12  // From benchmarks.
                val numLoops = polygon.numLoops()
                cumulativeEdges = if (numLoops > kMaxLinearSearchLoops) {
                    IntArray(numLoops)
                } else {
                    IntArray(0)
                }
                for (i in 0 until numLoops) {
                    if (cumulativeEdges.size > i) cumulativeEdges[i] = cumulativeNumEdges
                cumulativeNumEdges += polygon.loop(i).numVertices
            }
            }
            else {
                cumulativeEdges = IntArray(0)
            }
            numEdges = cumulativeNumEdges
        }

        // S2Shape interface:

        override fun edge(edgeId: Int): Edge {
            requireLT(edgeId, numEdges)
            var i = 0
            var e = edgeId
            if (cumulativeEdges.isNotEmpty()) {
                // "upper_bound" finds the loop just beyond the one we want.
                i = cumulativeEdges.upperBound( 0, polygon.numLoops(), edgeId) - 1
                e -= cumulativeEdges[i]
            } else {
                // When the number of loops is small, linear search is faster.  Most often
                // there is exactly one loop and the code below executes zero times.
                while (e >= polygon.loop(i).numVertices) {
                    e -= polygon.loop(i).numVertices
                    ++i
                }
            }
            return Edge(polygon.loop(i).orientedVertex(e), polygon.loop(i).orientedVertex(e + 1))
        }

        override val dimension: Int = 2

        override fun getReferencePoint(): ReferencePoint {
            var containsOrigin = false
            for (i in 0 until polygon.numLoops()) {
                containsOrigin = containsOrigin xor polygon.loop(i).containsOrigin()
            }
            return ReferencePoint(S2PointUtil.origin(), containsOrigin)
        }

        override val numChains: Int
            get() = polygon.numLoops()

        override fun chain(chainId: Int): Chain {
            requireLT(chainId, numChains)
            if (cumulativeEdges.isNotEmpty()) {
                return Chain(cumulativeEdges[chainId], polygon.loop(chainId).numVertices)
            } else {
                var e = 0
                for (j in 0 until chainId) e += polygon.loop(j).numVertices
                // S2Polygon represents a full loop as a loop with one vertex, while
                // S2Shape represents a full loop as a chain with no vertices.
                val numVertices = polygon.loop(chainId).numVertices
                return Chain(e, if(numVertices == 1) 0 else numVertices)
            }
        }

        override fun chainEdge(chainId: Int, offset: Int): Edge {
            logger.trace { "Shape.chain_edge | i = $chainId, j = $offset\nnum_chains = $numChains\nnum_vertices=${polygon.loop(chainId).numVertices}" }
            requireLT(chainId, numChains)
            requireLT(offset, polygon.loop(chainId).numVertices)
            return Edge(polygon.loop(chainId).orientedVertex(offset), polygon.loop(chainId).orientedVertex(offset + 1))
        }

        override fun chainPosition(edgeId: Int): ChainPosition {
            requireLT(edgeId, numEdges)
            var i = 0
            var e = edgeId
            if (cumulativeEdges.isNotEmpty()) {
                // "upper_bound" finds the loop just beyond the one we want.
                val start = cumulativeEdges.upperBound(0, polygon.numLoops(), edgeId) - 1
                i = start
                e -= cumulativeEdges[start]
            } else {
                // When the number of loops is small, linear search is faster.  Most often
                // there is exactly one loop and the code below executes zero times.
                while (e >= polygon.loop(i).numVertices) {
                    e -= polygon.loop(i).numVertices
                    ++i
                }
            }
            return ChainPosition(i, e)
        }

        override val typeTag: TypeTag = TypeTags.kPolygonTypeTag

    }

    companion object {

        private val logger = KotlinLogging.logger(S2Polygon::class.java.name)

        private val kMaxBruteForceVertices = 32
        private val kMaxUnindexedContainsCalls = 20

        // Return the overlap fractions between two polygons, i.e. the ratios of the
        // area of intersection to the area of each polygon.
        fun getOverlapFractions(a: S2Polygon, b: S2Polygon): Pair<Double, Double> {
            val intersection = S2Polygon()
            intersection.initToIntersection(a, b)
            val intersection_area = intersection.getArea()
            val a_area = a.getArea()
            val b_area = b.getArea()
            return Pair(
                    if (intersection_area >= a_area) 1.0 else intersection_area / a_area,
                    if (intersection_area >= b_area) 1.0 else intersection_area / b_area)
        }

        // Return a polygon which is the union of the given polygons.
        fun destructiveUnion(polygons: List<S2Polygon>): S2Polygon {
            return destructiveApproxUnion(polygons, S2EdgeCrossings.kIntersectionMergeRadius)
        }

        fun destructiveApproxUnion(polygons: List<S2Polygon>, snap_radius: S1Angle): S2Polygon {
            val queue: QueueType = TreeMap();  // Map from # of vertices to polygon.
            for (polygon in polygons) {
                queue.getOrPut(polygon.numVertices, { mutableListOf() }).add(polygon)
            }

            while (queue.polygonCount() > 1) {
                //QueueType::iterator smallest_it = queue.begin();
                val a_size = queue.firstKey()
                val a_polygonList = queue.getValue(a_size)
                val a_polygon = a_polygonList.removeAt(0)
                if (a_polygonList.isEmpty()) {
                    queue.remove(a_size)
                }
                val b_size = queue.firstKey()
                val b_polygonList = queue.getValue(b_size)
                val b_polygon = b_polygonList.removeAt(0)
                if (b_polygonList.isEmpty()) {
                    queue.remove(b_size)
                }

                // Union and add result back to queue.
                val union_polygon = S2Polygon()
                union_polygon.initToApproxUnion(a_polygon, b_polygon, snap_radius)
                queue.getOrPut(a_size + b_size, { mutableListOf() }).add(union_polygon)
                // We assume that the number of vertices in the union polygon is the
                // sum of the number of vertices in the original polygons, which is not
                // always true, but will almost always be a decent approximation, and
                // faster than recomputing.
            }

            if (queue.isEmpty())
                return S2Polygon()
            else
                return queue.firstEntry().value.first()
        }

    }

    // Given that loops_ contains a single loop, initialize all other fields.
    // This is an internal method that expects that loops_ has already been
    // initialized with a single non-empty loop.
    private fun initOneLoop() {
        logger.trace { "initOneLoop | loops = \n${loops.joinToString(separator = "\n") { l -> S2TextParser.toString(l) }}" }
        requireEQ(1, numLoops())
        val loop = loops[0]
        loop.depth = 0
        errorInconsistentLoopOrientations = false
        numVertices = loop.numVertices
        bound = loop.rectBound
        subregionBound = S2LatLngRectBounder.expandForSubregions(bound)
        initIndex()

        logger.trace { "initOneLoop | numVertices = $numVertices" }
        logger.trace { "initOneLoop | bound = $bound" }
        logger.trace { "initOneLoop | subregionBound = $subregionBound" }
    }

    // Compute num_vertices_, bound_, subregion_bound_.
    private fun initLoopProperties() {
        numVertices = 0
        bound = S2LatLngRect.empty()
        loops.indices.forEach { i ->
            if (loop(i).depth == 0) {
                bound = bound.union(loop(i).rectBound)
            }
            numVertices += loop(i).numVertices
        }
        subregionBound = S2LatLngRectBounder.expandForSubregions(bound)
        initIndex()
    }

    // Deletes the contents of the loops_ vector and resets the polygon state.
    private fun clearLoops() {
        logger.trace { "clearLoops | " }
        clearIndex()
        loops.clear()
        errorInconsistentLoopOrientations = false
    }

    // Return true if there is an error in the loop nesting hierarchy.
    private fun findLoopNestingError(): S2Error {
        // First check that the loop depths make sense.
        var last_depth = -1
        for (i in 0 until numLoops()) {
            val depth = loop(i).depth
            if (depth < 0 || depth > last_depth + 1) {
                return S2Error(code = S2Error.POLYGON_INVALID_LOOP_DEPTH, text = "Loop $i: invalid loop depth ($depth)")
            }
            last_depth = depth;
        }
        // Then check that they correspond to the actual loop nesting.  This test
        // is quadratic in the number of loops but the cost per iteration is small.
        for (i in 0 until numLoops()) {
            val last = getLastDescendant(i)
            for (j in 0 until numLoops()) {
                if (i == j) continue
                val nested = (j >= i + 1) && (j <= last)
                if (loop(i).containsNonCrossingBoundary(loop(j), false) != nested) {
                    return S2Error(code = S2Error.POLYGON_INVALID_LOOP_NESTING, text = "Invalid nesting: loop $i should %scontain loop $j".format(
                            if (nested) "" else "not "
                    ))
                }
            }
        }
        return S2Error(code = S2Error.OK)
    }


    private fun insertLoop(new_loop: S2Loop, loop_map: LoopMap) {
        var parent: S2Loop? = null
        var children: MutableList<S2Loop>
        var done: Boolean
        do {
            if (loop_map.containsKey(parent)) {
                children = loop_map.getValue(parent)
            } else {
                children = mutableListOf()
                loop_map[parent] = children
            }
            done = true
            for (child in children) {
                if (child.containsNested(new_loop)) {
                    parent = child
                    done = false
                    break
                }
            }
        } while (!done)

        // Some of the children of the parent loop may now be children of
        // the new loop.
        var new_children = loop_map[new_loop]
        if (new_children == null) {
            new_children = mutableListOf()
            loop_map[new_loop] = new_children
        }
        var i = 0
        while (i < children.size) {
            val child = children[i]
            if (new_loop.containsNested(child)) {
                new_children.add(child)
                children.removeAt(i)
            } else {
                ++i
            }
        }
        children.add(new_loop)
    }

    fun initLoops(loop_map: LoopMap) {
        val loopStack: Deque<S2Loop?> = LinkedList<S2Loop?>(listOf(null))
        var depth = -1
        while (!loopStack.isEmpty()) {
            val loop = loopStack.pop()
            if (loop != null) {
                depth = loop.depth
                loops.add(loop)
            }
            val children = loop_map.getOrDefault(loop, mutableListOf())
            for (i in (children.size - 1) downTo 0) {
                val child = children[i]
                child.depth = depth + 1
                loopStack.push(child);
            }
        }
    }

    // Add the polygon's loops to the S2ShapeIndex.  (The actual work of
    // building the index only happens when the index is first used.)
    private fun initIndex() {
        requireEQ(0, index.nextNewShapeId())
        index.add(Shape(polygon = this));
        //if (!FLAGS_s2polygon_lazy_indexing) {
        //    index_.ForceBuild();
        //}

        logger.trace { """initIndex | Loops:
            |${loops.joinToString("\n") { l -> S2TextParser.toString(l)} }
        """.trimMargin() }
        logger.trace { """initIndex | Index:
            |${index.toDebugString()}
        """.trimMargin() }

        if (PreConditions.enabled && debugOverride == S2Debug.ALLOW) {
            // Note that FLAGS_s2debug is false in optimized builds (by default).
            checkState { isValid() }
        }

    }

    // When the loop is modified (Invert(), or Init() called again) then the
    // indexing structures need to be cleared since they become invalid.
    private fun clearIndex() {
        logger.trace { "clearIndex |" }
        unindexedContainsCalls.set(0)
        index.removeAll()
    }

    // Initializes the polygon to the result of the given boolean operation,
    // returning an error on failure.
    fun initToOperation(opType: S2BooleanOperation.OpType, snap_function: SnapFunction, a: S2Polygon, b: S2Polygon, error: S2Error): Boolean {
        val options = S2BooleanOperation.Options()
        options.snapFunction = snap_function
        logger.trace { "initToOperation | opType = $opType\na = ${S2TextParser.toString(a)}\nb = ${S2TextParser.toString(b)}" }
        val op = S2BooleanOperation(opType, S2PolygonLayer(this), options)
        val ok = op.build(a.index, b.index, error)
        logger.trace { "initToOperation | ok = $ok: $error" }
        return ok
    }

    fun simplifyEdgesInCell(a: S2Polygon, cell: S2Cell, tolerance_uv: Double, snap_radius: S1Angle): List<S2Polyline> {
        val builder = S2Builder(S2Builder.Options(snapFunction = IdentitySnapFunction(snap_radius), simplifyEdgeChains = true))
        // The output consists of a sequence of polylines.  Polylines consisting of
        // interior edges are simplified using S2Builder, while polylines consisting
        // of boundary edges are returned unchanged.
        val polylines = mutableListOf<S2Polyline>()
        for (i in 0 until a.numLoops()) {
            val a_loop = a.loop(i)
            var v0 = a_loop.orientedVertex(0)
            var mask0 = getCellEdgeIncidenceMask(cell, v0, tolerance_uv)
            var in_interior = false  // Was the last edge an interior edge?
            for (j in 1..a_loop.numVertices) {
                val v1 = a_loop.orientedVertex(j)
                val mask1 = getCellEdgeIncidenceMask(cell, v1, tolerance_uv)
                if ((mask0 and mask1) != 0.toUByte()) {
                    // This is an edge along the cell boundary.  Such edges do not get
                    // simplified; we add them directly to the output.  (We create a
                    // separate polyline for each edge to keep things simple.)  We call
                    // ForceVertex on all boundary vertices to ensure that they don't
                    // move, and so that nearby interior edges are snapped to them.
                    checkState { !in_interior }
                    builder.forceVertex(v1)
                    polylines.add(S2Polyline(listOf(v0, v1)))
                } else {
                    // This is an interior edge.  If this is the first edge of an interior
                    // chain, then start a new S2Builder layer.  Also ensure that any
                    // polyline vertices on the boundary do not move, so that they will
                    // still connect with any boundary edge(s) there.
                    if (!in_interior) {
                        val polyline = S2Polyline()
                        builder.startLayer(S2PolylineLayer(polyline))
                        polylines.add(polyline)
                        in_interior = true
                    }
                    builder.addEdge(v0, v1)
                    if (mask1 != 0.toUByte()) {
                        builder.forceVertex(v1)
                        in_interior = false;  // Terminate this polyline.
                    }
                }
                v0 = v1;
                mask0 = mask1;
            }
        }
        val error = S2Error()
        if (!builder.build(error)) {
            logger.error { "InitToSimplifiedInCell failed: $error" }
        }
        return polylines
    }

    // Defines a total ordering on S2Loops that does not depend on the cyclic
    // order of loop vertices.  This function is used to choose which loop to
    // invert in the case where several loops have exactly the same area.
    // TODO(ericv): Consider adding this to the S2Loop API.  (May also want an
    // undirected version (CompareDirected vs CompareUndirected); should they
    // return a sign, or have separate "<" and "==" methods?)
    fun compareLoops(a: S2Loop, b: S2Loop): Int {
        if (a.numVertices != b.numVertices) {
            return a.numVertices - b.numVertices
        }
        val ao = a.getCanonicalLoopOrder()
        val bo = b.getCanonicalLoopOrder()
        if (ao.dir != bo.dir) return ao.dir - bo.dir
        var n = a.numVertices
        var ai = ao.first
        var bi = bo.first
        while (--n >= 0) {
            if (a.vertex(ai) < b.vertex(bi)) return -1
            if (a.vertex(ai) > b.vertex(bi)) return 1
            ai += ao.dir
            bi += bo.dir
        }
        return 0;

    }

}

