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
package dilivia.s2.shape

import dilivia.PreConditions.checkLE
import dilivia.PreConditions.requireEQ
import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.region.S2LoopMeasures
import dilivia.s2.region.S2LoopMeasures.getApproxArea
import dilivia.s2.region.S2LoopMeasures.getSignedArea
import dilivia.s2.region.S2PointLoopSpan
import dilivia.s2.region.S2Polyline
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs

// A 32-bit tag that can be used to identify the type of an encoded S2Shape.
// All encodable types have a non-zero type tag.  The tag associated with a
// given shape type can be accessed as Shape::kTypeTag, while the tag
// associated with a given object can be accessed as shape.type_tag().
//
// Type tags in the range 0..8191 are reserved for use by the S2 library.
typealias TypeTag = UInt

object TypeTags {
    val kNoTypeTag: TypeTag = 0U
    val kPolygonTypeTag: TypeTag = 1U
    val kPolylineTypeTag: TypeTag = 2U
    val kPointVectorTypeTag: TypeTag = 3U
    val kLaxPolylineTypeTag: TypeTag = 4U
    val kLaxPolygonTypeTag: TypeTag = 5U
    val kPointTypeTag: TypeTag = 6U
}

// An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
// allowed, and can be used to represent points.
data class Edge(val v0: S2Point = S2Point(), val v1: S2Point = S2Point()): Comparable<Edge> {

    override fun compareTo(other: Edge): Int {
        val v0Comparison = v0.compareTo(other.v0)
        return if (v0Comparison == 0) {
            v1.compareTo(other.v1)
        } else v0Comparison
    }

    fun reversed(): Edge = Edge(v0 = v1, v1 = v0)

    override fun toString(): String {
        return "Edge(v0=$v0, v1=$v1)"
    }


};
// The purpose of S2Shape is to represent polygonal geometry in a flexible
// way.  It is organized as a collection of edges that optionally defines an
// interior.  All geometry represented by an S2Shape must have the same
// dimension, which means that an S2Shape can represent either a set of
// points, a set of polylines, or a set of polygons.
//
// S2Shape is defined as an abstract base class in order to give clients
// control over the underlying data representation.  Sometimes an S2Shape does
// not have any data of its own, but instead "wraps" some other class.  There
// are various useful subtypes defined in *_shape.h, and some S2 classes also
// have a nested "Shape" class (e.g., S2Polygon::Shape).  It is easy for
// clients to implement their own subtypes, since the interface is minimal.
//
// S2Shape operations are typically defined on S2ShapeIndex objects rather
// than individual shapes.  An S2ShapeIndex is simply a collection of
// S2Shapes, possibly of different dimensions (e.g. 10 points and 3 polygons),
// organized into a data structure for efficient edge access.
//
// The edges of an S2Shape are identified by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, an S2Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an S2Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BooleanOperation).
abstract class S2Shape(var id: Int = -1) {

    // A range of edge ids corresponding to a chain of zero or more connected
    // edges, specified as a (start, length) pair.  The chain is defined to
    // consist of edge ids {start, start + 1, ..., start + length - 1}.
    data class Chain(val start: Int, val length: Int)

    // The position of an edge within a given edge chain, specified as a
    // (chain_id, offset) pair.  Chains are numbered sequentially starting from
    // zero, and offsets are measured from the start of each chain.
    data class ChainPosition(val chainId: Int, val offset: Int)

    // A ReferencePoint consists of a point P and a boolean indicating whether P
    // is contained by a particular shape.
    data class ReferencePoint(val point: S2Point = S2PointUtil.origin(), val contained: Boolean)

    // Returns the number of edges in this shape.  Edges have ids ranging from 0
    // to num_edges() - 1.
    abstract val numEdges: Int

    // Returns the endpoints of the given edge id.
    //
    // REQUIRES: 0 <= id < num_edges()
    abstract fun edge(edgeId: Int): Edge

    // Returns the dimension of the geometry represented by this shape.
    //
    //  0 - Point geometry.  Each point is represented as a degenerate edge.
    //
    //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
    //      represent any number of polylines.  Polylines edges may intersect.
    //
    //  2 - Polygon geometry.  Edges should be oriented such that the polygon
    //      interior is always on the left.  In theory the edges may be returned
    //      in any order, but typically the edges are organized as a collection
    //      of edge chains where each chain represents one polygon loop.
    //      Polygons may have degeneracies (e.g., degenerate edges or sibling
    //      pairs consisting of an edge and its corresponding reversed edge).
    //      A polygon loop may also be full (containing all points on the
    //      sphere); by convention this is represented as a chain with no edges.
    //      (See S2LaxPolygonShape for details.)
    //
    // Note that this method allows degenerate geometry of different dimensions
    // to be distinguished, e.g. it allows a point to be distinguished from a
    // polyline or polygon that has been simplified to a single point.
    abstract val dimension: Int

    // Returns true if the shape contains no points.  (Note that the full
    // polygon is represented as a chain with zero edges.)
    fun isEmpty(): Boolean {
        return numEdges == 0 && (dimension < 2 || numChains == 0);
    }
    // Returns true if the shape contains all points on the sphere.
    fun isFull(): Boolean {
        return numEdges == 0 && dimension == 2 && numChains > 0;
    }

    // Returns an arbitrary point P along with a boolean indicating whether P is
    // contained by the shape.  (The boolean value must be false for shapes that
    // do not have an interior.)
    //
    // This ReferencePoint may then be used to compute the containment of other
    // points by counting edge crossings.
    abstract fun getReferencePoint(): ReferencePoint

    // Returns the number of contiguous edge chains in the shape.  For example,
    // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
    // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
    // sequentially starting from zero.
    //
    // Note that it is always acceptable to implement this method by returning
    // num_edges() (i.e. every chain consists of a single edge), but this may
    // reduce the efficiency of some algorithms.
    abstract val numChains: Int

    // Returns the range of edge ids corresponding to the given edge chain.  The
    // edge chains must form contiguous, non-overlapping ranges that cover the
    // entire range of edge ids.  This is spelled out more formally below:
    //
    // REQUIRES: 0 <= i < num_chains()
    // REQUIRES: chain(i).length >= 0, for all i
    // REQUIRES: chain(0).start == 0
    // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
    //           for i < num_chains() - 1
    // REQUIRES: chain(i).start + chain(i).length == num_edges(),
    //           for i == num_chains() - 1
    abstract fun chain(chainId: Int): Chain

    // Returns the edge at offset "offset" within edge chain "chain_id".
    // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
    // but may be more efficient.
    abstract fun chainEdge(chainId: Int, offset: Int): Edge

    // Finds the chain containing the given edge, and returns the position of
    // that edge as a (chain_id, offset) pair.
    //
    // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
    // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
    //
    // where     pos == shape.chain_position(edge_id).
    abstract fun chainPosition(edgeId: Int): ChainPosition

    // Returns an integer that can be used to identify the type of an encoded
    // S2Shape (see TypeTag above).
    open val typeTag: TypeTag = kNoTypeTag

    // Virtual methods that return pointers of your choice.  These methods are
    // intended to help with the problem of attaching additional data to S2Shape
    // objects.  For example, you could return a pointer to a source object, or
    // a pointer to a bundle of additional data allocated with the S2Shape.
    // Because this method exists in all S2Shapes, you can override it in each
    // type of shape you have and call it without knowing the concrete subtype.
    // For example, if you have polyline and polygon shapes, you can do this:
    //
    //   class MyPolyline : public S2Polyline::Shape {
    //    public:
    //     virtual void* mutable_user_data() { return &my_data_; }
    //    private:
    //     MyData my_data_;
    //   };
    //   class MyPolygon : public S2Polygon::Shape {
    //    public:
    //     virtual void* mutable_user_data() { return &my_data_; }
    //    private:
    //     MyData my_data_;
    //   };
    //   ...
    //   S2Shape* shape = index.shape(id);
    //   const MyData* data = static_cast<const MyData*>(shape->user_data());
    //
    // This is not the only way to map from an S2Shape back to your source
    // data.  Other reasonable techniques include:
    //
    //  - Every shape has an id() assigned by S2ShapeIndex.  Ids are assigned
    //    sequentially starting from 0 in the order the shapes are added to the
    //    index.  You can use this id to look up arbitrary data stored in your
    //    own vector.
    //
    //  - If all of your shapes are the same type, then you can create your own
    //    subclass of some existing S2Shape type (such as S2Polyline::Shape) and
    //    add your own methods and fields.  You can access this data by
    //    downcasting the S2Shape pointers returned by S2ShapeIndex methods.
    open fun userData(): Any? = null

    companion object {

        private val logger = KotlinLogging.logger(S2Shape::class.java.name)

        // Indicates that a given S2Shape type cannot be encoded.
        val kNoTypeTag: TypeTag = 0U

        // The minimum allowable tag for user-defined S2Shape types.
        val kMinUserTypeTag: TypeTag = 8192.toUInt()
        // This is a helper function for implementing S2Shape::GetReferencePoint().
        //
        // Given a shape consisting of closed polygonal loops, the interior of the
        // shape is defined as the region to the left of all edges (which must be
        // oriented consistently).  This function then chooses an arbitrary point and
        // returns true if that point is contained by the shape.
        //
        // Unlike S2Loop and S2Polygon, this method allows duplicate vertices and
        // edges, which requires some extra care with definitions.  The rule that we
        // apply is that an edge and its reverse edge "cancel" each other: the result
        // is the same as if that edge pair were not present.  Therefore shapes that
        // consist only of degenerate loop(s) are either empty or full; by convention,
        // the shape is considered full if and only if it contains an empty loop (see
        // S2LaxPolygonShape for details).
        //
        // Determining whether a loop on the sphere contains a point is harder than
        // the corresponding problem in 2D plane geometry.  It cannot be implemented
        // just by counting edge crossings because there is no such thing as a "point
        // at infinity" that is guaranteed to be outside the loop.
        fun getReferencePoint(shape: S2Shape): ReferencePoint {
            requireEQ(shape.dimension, 2)
            if (shape.numEdges == 0) {
                // A shape with no edges is defined to be full if and only if it
                // contains at least one chain.
                return ReferencePoint(contained = shape.numChains > 0)
            }
            // Define a "matched" edge as one that can be paired with a corresponding
            // reversed edge.  Define a vertex as "balanced" if all of its edges are
            // matched. In order to determine containment, we must find an unbalanced
            // vertex.  Often every vertex is unbalanced, so we start by trying an
            // arbitrary vertex.
            val edge = shape.edge(0)
            var result = getReferencePointAtVertex(shape, edge.v0)
            if (result != null) {
                return result
            }
            // That didn't work, so now we do some extra work to find an unbalanced
            // vertex (if any).  Essentially we gather a list of edges and a list of
            // reversed edges, and then sort them.  The first edge that appears in one
            // list but not the other is guaranteed to be unmatched.
            val n = shape.numEdges
            val edges = ArrayList<Edge>(n)
            val revEdges = ArrayList<Edge>(n)
            for (i in 0 until n) {
                val e = shape.edge(i)
                edges.add(e)
                revEdges.add(Edge(e.v1, e.v0))
            }
            edges.sort()
            revEdges.sort()
            for (i in 0 until n) {
                if (edges[i] < revEdges[i]) {  // edges[i] is unmatched
                    result = getReferencePointAtVertex(shape, edges[i].v0)
                    check(result != null)
                    return result;
                }
                if (revEdges[i] < edges[i]) {  // rev_edges[i] is unmatched
                    result = getReferencePointAtVertex(shape, revEdges[i].v0)
                    check(result != null)
                    return result;
                }
            }
            // All vertices are balanced, so this polygon is either empty or full except
            // for degeneracies.  By convention it is defined to be full if it contains
            // any chain with no edges.
            for (i in 0 until shape.numChains) {
                if (shape.chain(i).length == 0) return ReferencePoint(contained = true)
            }
            return ReferencePoint(contained = false)
        }

        // This is a helper function for GetReferencePoint() above.
        //
        // If the given vertex "vtest" is unbalanced (see definition below), sets
        // "result" to a ReferencePoint indicating whther "vtest" is contained and
        // returns true.  Otherwise returns false.
        private fun getReferencePointAtVertex(shape: S2Shape, vtest: S2Point): ReferencePoint? {
            // Let P be an unbalanced vertex.  Vertex P is defined to be inside the
            // region if the region contains a particular direction vector starting from
            // P, namely the direction S2::Ortho(P).  This can be calculated using
            // S2ContainsVertexQuery.
            val contains_query = S2ContainsVertexQuery(vtest)
            val n = shape.numEdges
            for (e in 0 until n) {
                val edge = shape.edge(e)
                if (edge.v0 == vtest) contains_query.addEdge(edge.v1, 1)
                if (edge.v1 == vtest) contains_query.addEdge(edge.v0, -1)
            }
            val contains_sign = contains_query.containsSign()
            if (contains_sign == 0) {
                return null;  // There are no unmatched edges incident to this vertex.
            }
            return ReferencePoint(point = vtest, contained = contains_sign > 0)
        }

        // Returns true if the given shape contains the given point.  Most clients
        // should not use this method, since its running time is linear in the number
        // of shape edges.  Instead clients should create an S2ShapeIndex and use
        // S2ContainsPointQuery, since this strategy is much more efficient when many
        // points need to be tested.
        //
        // Polygon boundaries are treated as being semi-open (see S2ContainsPointQuery
        // and S2VertexModel for other options).
        //
        // CAVEAT: Typically this method is only used internally.  Its running time is
        //         linear in the number of shape edges.
        fun containsBruteForce(shape: S2Shape, focus: S2Point): Boolean {
            if (shape.dimension < 2) return false

            val refPoint = shape.getReferencePoint()
            if (refPoint.point == focus) {
                logger.trace { "Focus point $focus = ref point of shape $shape: Focus point is contained." }
                return refPoint.contained
            }

            val crosser = S2EdgeCrosser(refPoint.point, focus)
            var inside = refPoint.contained;
            for (e in 0 until  shape.numEdges) {
                val edge = shape.edge(e)
                val crossing = crosser.edgeOrVertexCrossing(edge.v0, edge.v1)
                inside = inside xor crossing

                //logger.trace { "Process edge $e: crossing = $crossing => inside = $inside" }
            }

            logger.trace { "Focus point $focus is inside ? $inside" }
            return inside
        }

        // For shapes of dimension 1, returns the sum of all polyline lengths on the
        // unit sphere.  Otherwise returns zero.  (See GetPerimeter for shapes of
        // dimension 2.)
        //
        // All edges are modeled as spherical geodesics.  The result can be converted
        // to a distance on the Earth's surface (with a worst-case error of 0.562%
        // near the equator) using the functions in s2earth.h.
        fun getLength(shape: S2Shape): S1Angle {
            if (shape.dimension != 1) return S1Angle.zero()
            val length = S1Angle()
            val vertices = mutableListOf<S2Point>()
            val numChains = shape.numChains
            for (chain_id in 0 until numChains) {
                getChainVertices(shape, chain_id, vertices)
                length += S2Polyline.getLength(vertices)
            }
            return length
        }

        // For shapes of dimension 2, returns the sum of all loop perimeters on the
        // unit sphere.  Otherwise returns zero.  (See GetLength for shapes of
        // dimension 1.)
        //
        // All edges are modeled as spherical geodesics.  The result can be converted
        // to a distance on the Earth's surface (with a worst-case error of 0.562%
        // near the equator) using the functions in s2earth.h.
        fun getPerimeter(shape: S2Shape): S1Angle {
            if (shape.dimension != 2) return S1Angle.zero()
            val perimeter = S1Angle()
            val vertices = mutableListOf<S2Point>()
            val numChains = shape.numChains
            for (chain_id in 0 until numChains) {
                getChainVertices(shape, chain_id, vertices)
                perimeter += S2LoopMeasures.getPerimeter(S2PointLoopSpan(vertices))
            }
            return perimeter
        }

        // For shapes of dimension 2, returns the area of the shape on the unit
        // sphere.  The result is between 0 and 4*Pi steradians.  Otherwise returns
        // zero.  This method has good relative accuracy for both very large and very
        // small regions.
        //
        // All edges are modeled as spherical geodesics.  The result can be converted
        // to an area on the Earth's surface (with a worst-case error of 0.900% near
        // the poles) using the functions in s2earth.h.
        fun getArea(shape: S2Shape): Double {
            if (shape.dimension != 2) return 0.0

            // Since S2Shape uses the convention that the interior of the shape is to
            // the left of all edges, in theory we could compute the area of the polygon
            // by simply adding up all the loop areas modulo 4*Pi.  The problem with
            // this approach is that polygons holes typically have areas near 4*Pi,
            // which can create large cancellation errors when computing the area of
            // small polygons with holes.  For example, a shell with an area of 4 square
            // meters (1e-13 steradians) surrounding a hole with an area of 3 square
            // meters (7.5e-14 sterians) would lose almost all of its accuracy if the
            // area of the hole was computed as 12.566370614359098.
            //
            // So instead we use S2::GetSignedArea() to ensure that all loops have areas
            // in the range [-2*Pi, 2*Pi].
            var area = 0.0
            val vertices = mutableListOf<S2Point>()
            val numChains = shape.numChains
            for (chain_id in 0 until numChains) {
                getChainVertices(shape, chain_id, vertices)
                area += getSignedArea(S2PointLoopSpan(vertices))
            }
            // Note that S2::GetSignedArea() guarantees that the full loop (containing
            // all points on the sphere) has a very small negative area.
            checkLE(abs(area), 4 * M_PI)
            if (area < 0.0) area += 4 * M_PI
            return area
        }

        // Like GetArea(), except that this method is faster and has more error.  The
        // additional error is at most 2.22e-15 steradians per vertex, which works out
        // to about 0.09 square meters per vertex on the Earth's surface.  For
        // example, a loop with 100 vertices has a maximum error of about 9 square
        // meters.  (The actual error is typically much smaller than this.)
        fun getApproxArea(shape: S2Shape): Double {
            if (shape.dimension != 2) return 0.0

            var area = 0.0
            val vertices = mutableListOf<S2Point>()
            val numChains = shape.numChains
            for (chain_id in 0 until numChains) {
                getChainVertices(shape, chain_id, vertices)
                area += getApproxArea(S2PointLoopSpan(vertices))
            }
            // Special case to ensure that full polygons are handled correctly.
            if (area <= 4 * M_PI) return area
            return area % (4 * M_PI)
        }

        // Returns the centroid of the shape multiplied by the measure of the shape,
        // which is defined as follows:
        //
        //  - For dimension 0 shapes, the measure is shape->num_edges().
        //  - For dimension 1 shapes, the measure is GetLength(shape).
        //  - For dimension 2 shapes, the measure is GetArea(shape).
        //
        // Note that the result is not unit length, so you may need to call
        // Normalize() before passing it to other S2 functions.
        //
        // The result is scaled by the measure defined above for two reasons: (1) it
        // is cheaper to compute this way, and (2) this makes it easier to compute the
        // centroid of a collection of shapes.  (This requires simply summing the
        // centroids of all shapes in the collection whose dimension is maximal.)
        fun getCentroid(shape: S2Shape): S2Point {
            var centroid = S2Point()
            val vertices = mutableListOf<S2Point>()
            val dimension = shape.dimension
            val numChains = shape.numChains
            for (chain_id in 0 until numChains) {
                when (dimension) {
                    0 -> centroid = centroid + shape.edge(chain_id).v0
                    1 -> {
                        getChainVertices(shape, chain_id, vertices)
                        centroid = centroid + S2Polyline.getCentroid(vertices)
                    }
                    else -> {
                        getChainVertices(shape, chain_id, vertices)
                        centroid = centroid + S2LoopMeasures.getCentroid(S2PointLoopSpan(vertices))
                    }
                }
            }
            return centroid
        }

        // Overwrites "vertices" with the vertices of the given edge chain of "shape".
        // If dimension == 1, the chain will have (chain.length + 1) vertices, and
        // otherwise it will have (chain.length) vertices.
        //
        // This is a low-level helper method used in the implementations of some of
        // the methods above.
        fun getChainVertices(shape: S2Shape, chainId: Int, vertices: MutableList<S2Point>) {
            val chain = shape.chain(chainId)
            val numVertices = chain.length + (if(shape.dimension == 1) 1 else 0)
            vertices.clear()
            var e = 0
            if ((numVertices and 1) != 0) {
                vertices.add(shape.chainEdge(chainId, e++).v0)
            }
            while (e < numVertices) {
                val edge = shape.chainEdge(chainId, e)
                vertices.add(edge.v0)
                vertices.add(edge.v1)
                e += 2
            }
        }
    }

}
