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
package dilivia.s2.index.shape

import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkLE
import dilivia.PreConditions.requireArgument
import dilivia.PreConditions.requireEQ
import dilivia.collections.assign
import dilivia.collections.lowerBound
import dilivia.collections.resize
import dilivia.collections.resizeWith
import dilivia.s2.S1Angle
import dilivia.s2.S2Error
import dilivia.s2.S2Measures.turnAngle
import dilivia.s2.S2Predicates
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.InputEdgeIdSetId
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.graph.VertexId
import dilivia.s2.builder.graph.VertexInMap
import dilivia.s2.builder.graph.VertexOutMap
import dilivia.s2.builder.layers.Layer
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.index.lexicon.ValueLexicon
import dilivia.s2.index.shape.S2BooleanOperation.OpType.DIFFERENCE
import dilivia.s2.index.shape.S2BooleanOperation.OpType.INTERSECTION
import dilivia.s2.index.shape.S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE
import dilivia.s2.index.shape.S2BooleanOperation.SourceId
import dilivia.s2.index.shape.S2BooleanOperationImpl.Companion.getInputEdgeChainOrder
import dilivia.s2.index.shape.S2BooleanOperationImpl.Companion.kSentinel
import dilivia.s2.index.shape.S2BooleanOperationImpl.Companion.kSetInside
import dilivia.s2.index.shape.S2BooleanOperationImpl.Companion.kSetInvertB
import dilivia.s2.index.shape.S2BooleanOperationImpl.Companion.kSetReverseA
import dilivia.s2.index.shape.S2BooleanOperationImpl.CrossingGraphEdge
import dilivia.s2.index.shape.S2BooleanOperationImpl.CrossingInputEdge
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import java.util.*
import dilivia.s2.builder.Edge as GraphEdge

/**
 * This class implements boolean operations (intersection, union, difference, and symmetric difference) for regions
 * whose boundaries are defined by geodesic edges.
 *
 * S2BooleanOperation operates on exactly two input regions at a time. Each region is represented as an S2ShapeIndex
 * and may contain any number of points, polylines, and polygons. The region is essentially the union of these objects,
 * except that polygon interiors must be disjoint from all other geometry (including other polygon interiors). If the
 * input geometry for a region does not meet this condition, it can be normalized by computing its union first.
 * Note that points or polylines are allowed to coincide with the boundaries of polygons.
 *
 * Degeneracies are supported. A polygon loop or polyline may consist of a single edge from a vertex to itself, and
 * polygons may contain "sibling pairs" consisting of an edge and its corresponding reverse edge. Polygons must not have
 * any duplicate edges (due to the requirement that polygon interiors are disjoint), but polylines may have duplicate
 * edges or can even be self-intersecting.
 *
 * Points and polyline edges are treated as multisets: if the same point or polyline edge appears multiple times in the
 * input, it will appear multiple times in the output. For example, the union of a point with an identical point
 * consists of two points. This feature is useful for modeling large sets of points or polylines as a single region
 * while maintaining their distinct identities, even when the points or polylines intersect each other. It is also
 * useful for reconstructing polylines that loop back on themselves. If duplicate geometry is not desired, it can be
 * merged by GraphOptions.DuplicateEdges.MERGE in the S2Builder output layer.
 *
 * Polylines are always considered to be directed. Polyline edges between the same pair of vertices are defined to
 * intersect even if the two edges are in opposite directions.  (Undirected polylines can be modeled by specifying
 * GraphOptions.EdgeType.UNDIRECTED in the S2Builder output layer.)
 *
 * The output of each operation is sent to an S2Builder.Layer provided by the client. This allows clients to build any
 * representation of the geometry they choose.  It also allows the client to do additional postprocessing of
 * the output before building data structures; for example, the client can easily discard degeneracies or convert them
 * to another data type. The boundaries of polygons and polylines can be modeled as open, semi-open, or closed.
 * Polyline boundaries are controlled by the PolylineModel class, whose options are as follows:
 *
 *  - In the OPEN model, polylines do not contain their first or last vertex except for one special case: namely, if
 *    the polyline forms a loop and the polyline_loops_have_boundaries() option is set to false, then the first/last
 *    vertex is contained.
 *
 *  - In the SEMI_OPEN model, polylines contain all vertices except the last. Therefore if one polyline starts where
 *    another polyline stops, the two polylines do not intersect.
 *
 *  - In the CLOSED model, polylines contain all of their vertices.
 *
 * When multiple polylines are present, they are processed independently and have no effect on each other. For example,
 * in the OPEN boundary model the polyline ABC contains the vertex B, while set of polylines {AB, BC} does not.
 * (If you want to treat the polylines as a union instead, with boundaries merged according to the "mod 2" rule,
 * this can be achieved by reassembling the edges into maximal polylines using S2PolylineVectorLayer with
 * EdgeType.UNDIRECTED, DuplicateEdges.MERGE, and PolylineType.WALK.)
 *
 * Polygon boundaries are controlled by the PolygonModel class, which has the following options:
 *
 *  - In the OPEN model, polygons do not contain their vertices or edges. This implies that a polyline that follows the
 *    boundary of a polygon will not intersect it.
 *
 *  - In the SEMI_OPEN model, polygon point containment is defined such that if several polygons tile the region around
 *    a vertex, then exactly one of those polygons contains that vertex.  Similarly polygons contain all of their edges,
 *    but none of their reversed edges.  This implies that a polyline and polygon edge with the same endpoints
 *    intersect if and only if they are in the same direction.  (This rule ensures that if a polyline is intersected
 *    with a polygon and its complement, the two resulting polylines do not have any edges in common.)
 *
 *  - In the CLOSED model, polygons contain all their vertices, edges, and reversed edges.  This implies that a
 *    polyline that shares an edge (in either direction) with a polygon is defined to intersect it.  Similarly, this
 *    is the only model where polygons that touch at a vertex or along an edge intersect.
 *
 * Note that PolylineModel and PolygonModel are defined as separate classes in order to allow for possible future
 * extensions.
 *
 * Operations between geometry of different dimensions are defined as follows:
 *
 *  - For UNION, the higher-dimensional shape always wins.  For example the union of a closed polygon A with a polyline
 *    B that coincides with the boundary of A consists only of the polygon A.
 *
 *  - For INTERSECTION, the lower-dimensional shape always wins.  For example, the intersection of a closed polygon A
 *    with a point B that coincides with a vertex of A consists only of the point B.
 *
 *  - For DIFFERENCE, higher-dimensional shapes are not affected by subtracting lower-dimensional shapes. For example,
 *    subtracting a point or polyline from a polygon A yields the original polygon A. This rule exists because in
 *    general, it is impossible to represent the output using the specified boundary model(s). (Consider subtracting
 *    one vertex from a PolylineModel.CLOSED polyline, or subtracting one edge from a PolygonModel.CLOSED polygon.)
 *    If you want to perform operations like this, consider representing all boundaries explicitly (topological
 *    boundaries) using OPEN boundary models.  Another option for polygons is to subtract a degenerate loop, which
 *    yields a polygon with a degenerate hole (see S2LaxPolygonShape).
 *
 * Note that in the case of Precision::EXACT operations, the above remarks only apply to the output before snapping.
 * Snapping may cause nearby distinct edges to become coincident, e.g. a polyline may become coincident with a polygon
 * boundary.  However also note that S2BooleanOperation is perfectly happy to accept such geometry as input.
 *
 * Note the following differences between S2BooleanOperation and the similar S2MultiBooleanOperation class:
 *
 *  - S2BooleanOperation operates on exactly two regions at a time, whereas S2MultiBooleanOperation operates on any
 *    number of regions.
 *
 *  - S2BooleanOperation is potentially much faster when the input is already represented as S2ShapeIndexes.
 *    The algorithm is output sensitive and is often sublinear in the input size.  This can be a big advantage if, say,
 *
 *  - S2BooleanOperation supports exact predicates and the corresponding exact operations (i.e., operations that are
 *    equivalent to computing the exact result and then snap rounding it).
 *
 *  - S2MultiBooleanOperation has better error guarantees when there are many regions, since it requires only one
 *    snapping operation for any number of input regions.
 *
 * Example usage:
 *   val a = S2ShapeIndex()
 *   val b = S2ShapeIndex()     // Input geometry, e.g. containing polygons.
 *   val polygon = S2Polygon()  // Output geometry.
 *   val options = S2BooleanOperation.Options()
 *   options.snapFunction = snapFunction
 *   val op = S2BooleanOperation(OpType.INTERSECTION, S2PolygonLayer(polygon), options)
 *   val error = S2Error()
 *   if (!op.build(a, b, error)) {
 *     logger.error { "Error: $error" }
 *     ...
 *   }
 *
 * If the output includes objects of different dimensions, they can be assembled into different layers with code like
 * this:
 *
 *   val points = mutableListOf<S2Point>()
 *   val polylines = S2Polyline()
 *   val polygon = S2Polygon()
 *   val op = S2BooleanOperation(
 *       OpType.UNION,
 *       PointVectorLayer(points),
 *       S2PolylineVectorLayer(polylines),
 *       S2PolygonLayer(polygon))
*/
class S2BooleanOperation
/**
 * Primary constructor.
 * Specifies that "result_empty" should be set to indicate whether the exact result of the operation is empty.
 * This constructor is used to efficiently test boolean relationships (see isEmpty above).
 * The following field is set if and only if there are no output layers.
 *
 * @param opType The operation type.
 * @param resultEmpty The boolean result destination or null if none.
 * @param options The operation options.
 */
private constructor(
    /** Operation type. */
    val opType: OpType,
    /** The destination result for boolean operation. */
    val resultEmpty: BooleanResult? = null,
    /** The operation options. */
    val options: Options
) {

    ////////////////////////////////// Inner Types /////////////////////////////////

    /** The supported operation types. */
    enum class OpType {
        /** Contained by either region. */
        UNION,
        /** Contained by both regions. */
        INTERSECTION,
        /** Contained by the first region but not the second. */
        DIFFERENCE,
        /** // Contained by one region but not the other. */
        SYMMETRIC_DIFFERENCE
    }

    /** Defines whether polygons are considered to contain their vertices and/or edges (see definitions above). */
    enum class PolygonModel { OPEN, SEMI_OPEN, CLOSED }

    /** Defines whether polylines are considered to contain their endpoints (see definitions above). */
    enum class PolylineModel { OPEN, SEMI_OPEN, CLOSED }

    /**
     * With Precision.EXACT, the operation is evaluated using the exact input geometry. Predicates that use this option
     * will produce exact results; for example, they can distinguish between a polyline that barely intersects a
     * polygon from one that barely misses it.  Constructive operations (ones that yield new geometry, as opposed to
     * predicates) are implemented by computing the exact result and then snap rounding it according to the given
     * snapFunction (see below). This is as close as it is possible to get to the exact result while requiring that
     * vertex coordinates have type "double".
     *
     * With Precision.SNAPPED, the input regions are snapped together *before* the operation is evaluated. So for
     * example, two polygons that overlap slightly will be treated as though they share a common boundary, and similarly
     * two polygons that are slightly separated from each other will be treated as though they share a common boundary.
     * Snapped results are useful for dealing with points, since in S2 the only points that lie exactly on a polyline
     * or polygon edge are the endpoints of that edge.
     *
     * Conceptually, the difference between these two options is that with Precision.SNAPPED, the inputs are snap
     * rounded (together), whereas with Precision.EXACT only the result is snap rounded.
     */
    enum class Precision { EXACT, SNAPPED }

    /**
     * SourceId identifies an edge from one of the two input S2ShapeIndexes. It consists of a region id (0 or 1),
     * a shape id within that region's S2ShapeIndex, and an edge id within that shape.
     *
     * @property regionId The identifier of the region (0: index a, 1: index b)
     * @property shapeId The identifier of the shape in its index.
     * @property edgeId The identifier of the edge in the shape.
     */
    data class SourceId(val regionId: Int = 0, val shapeId: Int = 0, val edgeId: Int = -1) : Comparable<SourceId> {

        override fun compareTo(other: SourceId): Int = when {
            regionId < other.regionId -> -1
            regionId > other.regionId -> 1
            shapeId < other.shapeId -> -1
            shapeId > other.shapeId -> 1
            else -> edgeId.compareTo(other.edgeId)
        }

    }

    /**
     * Defines options of the S2BooleanOperation.
     */
    data class Options(

        /**
         * Specifies the function to be used for snap rounding.
         *
         * DEFAULT: IdentitySnapFunction(S1Angle.zero())
         *  - This does no snapping and preserves all input vertices exactly unless there are crossing edges, in which
         *    case the snap radius is increased to the maximum intersection point error (kIntersectionError).
         */
        var snapFunction: SnapFunction = IdentitySnapFunction(S1Angle.zero()),

        /**
         * Defines whether polygons are considered to contain their vertices and/or edges (see comments above).
         *
         * DEFAULT: PolygonModel.SEMI_OPEN
         */
        var polygonModel: PolygonModel = PolygonModel.SEMI_OPEN,

        /**
         * Defines whether polylines are considered to contain their vertices (see comments above).
         *
         * DEFAULT: PolylineModel.CLOSED
         */
        var polylineModel: PolylineModel = PolylineModel.CLOSED,

        /**
         * Specifies whether a polyline loop is considered to have a non-empty boundary. By default this option is true,
         * meaning that even if the first and last vertices of a polyline are the same, the polyline is considered to
         * have a well-defined "start" and "end".  For example, if the polyline boundary model is OPEN then the polyline
         * loop would not include the start/end vertices. These are the best semantics for most applications, such as
         * GPS tracks or road network segments.
         *
         * If the polyline forms a loop and this option is set to false, then instead the first and last vertices are
         * considered to represent a single vertex in the interior of the polyline.  In this case the boundary of the
         * polyline is empty, meaning that the first/last vertex will be contained by the polyline even if the boundary
         * model is OPEN.
         * (Note that this option also has a small effect on the CLOSED boundary model, because the first/last vertices
         * of a polyline loop are considered to represent one vertex rather than two.)
         *
         * The main reason for this option is to implement the "mod 2 union" boundary semantics of the OpenGIS Simple
         * Features spec.  This can be achieved by making sure that all polylines are constructed using
         * Graph.PolylineType.WALK (which ensures that all polylines are as long as possible), and then setting this
         * option to false.
         *
         * DEFAULT: true
         */
        var polylineLoopsHaveBoundaries: Boolean = true,

        /**
         * Specifies whether the operation should use the exact input geometry (Precision.EXACT), or whether the two
         * input regions should be snapped together first (Precision.SNAPPED).
         *
         * DEFAULT: Precision.EXACT
         */
        val precision: Precision = Precision.EXACT,

        /**
         * If true, the input geometry is interpreted as representing nearby geometry that has been snapped or
         * simplified.  It then outputs a conservative result based on the value of polygon_model() and polylineModel.
         * For the most part, this only affects the handling of degeneracies.
         *
         *  - If the model is OPEN, the result is as open as possible.  For example, the intersection of two identical
         *    degenerate shells is empty under PolygonModel.OPEN because they could have been disjoint before snapping.
         *    Similarly, two identical degenerate polylines have an empty intersection under PolylineModel.OPEN.
         *
         *  - If the model is CLOSED, the result is as closed as possible.  In the case of the DIFFERENCE operation,
         *    this is equivalent to evaluating A - B as Closure(A) - Interior(B).  For other operations, it affects
         *    only the handling of degeneracies.  For example, the union of two identical degenerate holes is empty
         *    under PolygonModel.CLOSED (i.e., the hole disappears) because the holes could have been disjoint before
         *    snapping.
         *
         *  - If the model is SEMI_OPEN, the result is as degenerate as possible. New degeneracies will not be created,
         *    but all degeneracies that coincide with the opposite region's boundary are retained unless this would
         *    cause a duplicate polygon edge to be created.  This model is is very useful for working with input data
         *    that has both positive and negative degeneracies (i.e., degenerate shells and holes).
         *
         * DEFAULT: false
         */
        val conservative_output: Boolean = false,

        /**
         * If specified, then each output edge will be labelled with one or more SourceIds indicating which input
         * edge(s) it corresponds to.  This can be useful if your input geometry has additional data that needs to
         * be propagated from the input to the output (e.g., elevations).
         *
         * You can access the labels by using an Layer type that supports labels, such as S2PolygonLayer. The layer
         * outputs a "label_set_lexicon" and an "label_set_id" for each edge.  You can then look up the source
         * information for each edge like this:
         *
         * for (label inn label_set_lexicon.id_set(label_set_id)) {
         *   val src: SourceId = source_id_lexicon.value(label)
         *   // region_id() specifies which S2ShapeIndex the edge is from (0 or 1).
         *   DoSomething(src.regionId, src.shapeId, src.edgeId)
         * }
         *
         * DEFAULT: nullptr
         */
        val sourceIdLexicon: ValueLexicon<SourceId>? = null,
    )

    /**
     * Boolean result.
     * @property value The result value.
     */
    data class BooleanResult(var value: Boolean = false)

    /////////////////////////////////// Fields /////////////////////////////////

    /** The input region a. */
    internal lateinit var a: S2ShapeIndex
    /** The input region b. */
    internal lateinit var b: S2ShapeIndex
    /** The output consists either of zero layers, one layer, or three layers. */
    internal lateinit var layers: MutableList<Layer>

    /////////////////////////////////// Init ///////////////////////////////////

    /**
     *
     * @param opType The operation type.
     * @param layer The destination layer.
     * @param options The operation options.
     */
    constructor(opType: OpType, layer: Layer, options: Options = Options()) : this(opType, null, options) {
        this.layers = mutableListOf(layer)
    }

    /**
     * Specifies that the output boundary edges should be sent to three different layers according to their dimension.
     * Points (represented by degenerate edges) are sent to layer 0, polyline edges are sent to layer 1, and polygon
     * edges are sent to layer 2.
     *
     * The dimension of an edge is defined as the minimum dimension of the two input edges that produced it. For
     * example, the intersection of two crossing polyline edges is a considered to be a degenerate polyline rather than
     * a point, so it is sent to layer 1.  Clients can easily reclassify such polylines as points if desired, but
     * this rule makes it easier for clients that want to process point, polyline, and polygon inputs differently.
     *
     * The layers are always built in the order 0, 1, 2, and all arguments to the build() calls are guaranteed to be
     * valid until the last call returns.
     * All Graph objects have the same set of vertices and the same lexicon objects, in order to make it easier to write
     * classes that process all the edges in parallel.
     *
     * @param opType The operation type.
     * @param layers The destination layers.
     * @param options The operation options.
     */
    constructor(opType: OpType, layers: MutableList<Layer>, options: Options = Options()) : this(opType, null, options) {
        this.layers = layers
    }

    /**
     * Executes the given operation.
     *
     * @return true on success, and otherwise sets "error" appropriately.  (This class does not generate any errors
     * itself, but the Layer might.)
     */
    fun build(a: S2ShapeIndex, b: S2ShapeIndex, error: S2Error): Boolean {
        this.a = a
        this.b = b
        return S2BooleanOperationImpl(this).build(error)
    }

    ///////////////////////////////// Methods ///////////////////////////////

    /**
     * Gets the index corresponding to the given regionId (0 -> a, 1 -> b).
     *
     * @param regionId The region identifier.
     * @return The shape index of the region.
     */
    fun getShapeIndex(regionId: Int): S2ShapeIndex = when (regionId) {
        0 -> a
        1 -> b
        else -> throw IllegalArgumentException("Wrong region id $regionId")
    }

    companion object {

        private val logger = KotlinLogging.logger(S2BooleanOperation::class.java.name)

        /**
         * Convenience method that returns true if the result of the given operation is empty.
         *
         * @param opType The operation type.
         * @param a The region index a.
         * @param b The region index b.
         * @param options The operation options.
         * @return true if the operation is empty and false otherwise.
         */
        fun isEmpty(opType: OpType, a: S2ShapeIndex, b: S2ShapeIndex, options: Options = Options()): Boolean {
            val result = BooleanResult()
            val op = S2BooleanOperation(opType, result, options)
            val error = S2Error()
            op.build(a, b, error)
            check(error.isOk())
            return result.value
        }

        /**
         * Convenience method that returns true if A intersects B.
         *
         * @param a The region index a.
         * @param b The region index b.
         * @param options The operation options.
         * @return true if A intersects B and false otherwise.
         */
        fun intersects(a: S2ShapeIndex, b: S2ShapeIndex, options: Options = Options()): Boolean =
                !isEmpty(INTERSECTION, b, a, options)

        /**
         * Convenience method that returns true if A contains B, i.e., if the difference (B - A) is empty.
         *
         * @param a The region index a.
         * @param b The region index b.
         * @param options The operation options.
         * @return true if A contains B and false otherwise.
         */
        fun contains(a: S2ShapeIndex, b: S2ShapeIndex, options: Options = Options()): Boolean = isEmpty(DIFFERENCE, b, a, options)

        // Convenience method that returns true if the symmetric difference of A and
        // B is empty.  (Note that A and B may still not be identical, e.g. A may
        // contain two copies of a polyline while B contains one.)
        /**
         * Convenience method that returns true if the symmetric difference of A and B is empty. (Note that A and B
         * may still not be identical, e.g. A may contain two copies of a polyline while B contains one.)
         *
         * @param a The region index a.
         * @param b The region index b.
         * @param options The operation options.
         * @return true if symmetric difference of A and B is empty.
         */
        fun equals(a: S2ShapeIndex, b: S2ShapeIndex, options: Options = Options()): Boolean = isEmpty(SYMMETRIC_DIFFERENCE, b, a, options)
    }

}

//
// Boolean operations are implemented by constructing the boundary of the
// result and then using S2Builder to assemble the edges.  The boundary is
// obtained by clipping each of the two input regions to the interior or
// exterior of the other region.  For example, to compute the union of A and
// B, we clip the boundary of A to the exterior of B and the boundary of B to
// the exterior of A; the resulting set of edges defines the union of the two
// regions.
//
// We use exact predicates, but inexact constructions (e.g. computing the
// intersection point of two edges).  Nevertheless, the following algorithm is
// guaranteed to be 100% robust, in that the computed boundary stays within a
// small tolerance (snap_radius + S2::kIntersectionError) of the exact
// result, and also preserves the correct topology (i.e., no crossing edges).
//
// Unfortunately this robustness cannot quite be achieved using the strategy
// outlined above (clipping the two input regions and assembling the
// resulting edges).  Since computed intersection points are not exact, the
// input geometry passed to S2Builder might contain self-intersections, and
// these self-intersections cannot be eliminated reliably by snap rounding.
//
// So instead, we pass S2Builder the entire set of input edges where at least
// some portion of each edge belongs to the output boundary.  We allow
// S2Builder to compute the intersection points and snap round the edges
// (which it does in a way that is guaranteed to preserve the input topology).
// Then once this is finished, we remove the portions of each edge that would
// have been clipped if we had done the clipping first.  This step only
// involves deciding whether to keep or discard each edge in the output, since
// all intersection points have already been resolved, and therefore there is
// no risk of creating new self-intersections.
//
// This is implemented using the following classes:
//
//  - S2BooleanOperation::Impl: the top-level class that clips each of
//                              the two regions to the other region.
//
//  - CrossingProcessor: a class that processes edge crossings and maintains
//                       the necessary state in order to clip the boundary
//                       of one region to the interior or exterior of the
//                       other region.
//
//  - EdgeClippingLayer: an S2Builder::Layer that removes graph edges that
//                       correspond to clipped portions of input edges, and
//                       passes the result to another layer for assembly.
//
//  - GraphEdgeClipper: a helper class that does the actual work of the
//                      EdgeClippingLayer.


// Given a set of clipping instructions encoded as a set of InputEdgeCrossings,
// GraphEdgeClipper determines which graph edges correspond to clipped
// portions of input edges and removes them.
//
// The clipping model is as follows.  The input consists of edge chains.  The
// clipper maintains an "inside" boolean state as it clips each chain, and
// toggles this state whenever an input edge is crossed.  Any edges that are
// deemed to be "outside" after clipping are removed.
//
// The "inside" state can be reset when necessary (e.g., when jumping to the
// start of a new chain) by adding a special crossing marked kSetInside.
// There are also two other special "crossings" that modify the clipping
// parameters: kSetInvertB specifies that edges should be clipped to the
// exterior of the other region, and kSetReverseA specifies that edges should
// be reversed before emitting them (which is needed to implement difference
// operations).
// "input_dimensions" is a vector specifying the dimension of each input
// edge (0, 1, or 2).  "input_crossings" is the set of all crossings to be
// used when clipping the edges of "g", sorted in lexicographic order.
//
// The clipped set of edges and their corresponding set of input edge ids
// are returned in "new_edges" and "new_input_edge_ids".  (These can be used
// to construct a new S2Builder::Graph.)
class GraphEdgeClipper(
        val g: Graph,
        val inputDimensions: List<Byte>,
        val inputCrossings: InputEdgeCrossings,
        val newEdges: MutableList<dilivia.s2.builder.Edge>,
        val newInputEdgeIds: MutableList<InputEdgeIdSetId>
) {
    private val logger = KotlinLogging.logger {  }

    private val inMap: VertexInMap = VertexInMap(g)
    private val outMap: VertexOutMap = VertexOutMap(g)

    // Every graph edge is associated with exactly one input edge in our case,
    // which means that we can declare g_.input_edge_id_set_ids() as a vector of
    // InputEdgeIds rather than a vector of InputEdgeIdSetIds.  (This also takes
    // advantage of the fact that IdSetLexicon represents a singleton set as the
    // value of its single element.)
    private val inputIds: List<InputEdgeId> = g.inputEdgeIdSetIds

    private val order: List<EdgeId> = getInputEdgeChainOrder(g, inputIds)  // Graph edges sorted in input edge id order.
    private val rank: ArrayList<Int> = ArrayList(order.size)      // The rank of each graph edge within order_.

    init {
        rank.assign(order.size, 0)
        for (i in order.indices) {
            rank[order[i]] = i
        }

        logger.trace { "constructor | g = ${g.toDebugString()}" }
        logger.trace { "constructor | inputIds = $inputIds" }
        logger.trace { "constructor | order = $order" }
    }

    fun run() {
        // Declare vectors here and reuse them to avoid reallocation.
        val a_vertices = mutableListOf<VertexId>()
        val a_num_crossings = mutableListOf<Int>()
        val a_isolated = mutableListOf<Boolean>()
        val b_input_edges = mutableListOf<CrossingInputEdge>()
        val b_edges = mutableListOf<CrossingGraphEdgeVector>()

        var inside = false
        var invert_b = false
        var reverse_a = false
        var nextIdx = 0
        var i = 0
        while (i < order.size) {
            // For each input edge (the "A" input edge), gather all the input edges
            // that cross it (the "B" input edges).
            val aInputId = inputIds[order[i]]
            val edge0 = g.edge(order[i])
            logger.trace { "run | iteration i = $i, aInputId = $aInputId, edge0 = $edge0" }

            b_input_edges.clear()
            while (nextIdx < inputCrossings.size) {
                val next = inputCrossings[nextIdx]
                if (next.first != aInputId) {
                    break
                }
                if (next.second.inputId >= 0) {
                    b_input_edges.add(next.second)
                } else if (next.second.inputId == kSetInside) {
                    inside = next.second.leftToRight
                } else if (next.second.inputId == kSetInvertB) {
                    invert_b = next.second.leftToRight
                } else {
                    checkEQ(next.second.inputId, kSetReverseA)
                    reverse_a = next.second.leftToRight
                }
                ++nextIdx
            }
            // Optimization for degenerate edges.
            // TODO(ericv): If the output layer for this edge dimension specifies
            // DegenerateEdges::DISCARD, then remove the edge here.
            if (edge0.first == edge0.second) {
                inside = inside xor ((b_input_edges.size and 1) != 0)
                addEdge(edge0, aInputId)
                ++i
                continue
            }
            // Optimization for the case where there are no crossings.
            if (b_input_edges.isEmpty()) {
                // In general the caller only passes edges that are part of the output
                // (i.e., we could S2_DCHECK(inside) here).  The one exception is for
                // polyline/polygon operations, where the polygon edges are needed to
                // compute the polyline output but are not emitted themselves.
                if (inside) {
                    addEdge(if(reverse_a) Graph.reverse(edge0) else edge0, aInputId)
                }
                ++i
                continue
            }
            // Walk along the chain of snapped edges for input edge A, and at each
            // vertex collect all the incident edges that belong to one of the
            // crossing edge chains (the "B" input edges).
            a_vertices.clear()
            a_vertices.add(edge0.first)
            b_edges.clear()
            b_edges.resizeWith(b_input_edges.size) { ArrayList() }
            gatherIncidentEdges(a_vertices, 0, b_input_edges, b_edges)
            while (i < order.size && inputIds[order[i]] == aInputId) {
                a_vertices.add(g.edge(order[i]).second)
                gatherIncidentEdges(a_vertices, a_vertices.size - 1, b_input_edges, b_edges)
                ++i
            }
            --i

            // Now for each B edge chain, decide which vertex of the A chain it
            // crosses, and keep track of the number of signed crossings at each A
            // vertex.  The sign of a crossing depends on whether the other edge
            // crosses from left to right or right to left.
            //
            // This would not be necessary if all calculations were done in exact
            // arithmetic, because crossings would have strictly alternating signs.
            // But because we have already snapped the result, some crossing locations
            // are ambiguous, and GetCrossedVertexIndex() handles this by choosing a
            // candidate vertex arbitrarily.  The end result is that rarely, we may
            // see two crossings in a row with the same sign.  We correct for this by
            // adding extra output edges that essentially link up the crossings in the
            // correct (alternating sign) order.  Compared to the "correct" behavior,
            // the only difference is that we have added some extra sibling pairs
            // (consisting of an edge and its corresponding reverse edge) which do not
            // affect the result.
            a_num_crossings.clear();
            a_num_crossings.resize(a_vertices.size)
            a_isolated.clear();
            a_isolated.resize(a_vertices.size)
            for (bi in 0 until b_input_edges.size) {
                val leftToRight = b_input_edges[bi].leftToRight
                val aIndex = getCrossedVertexIndex (a_vertices, b_edges[bi], leftToRight)
                if (aIndex >= 0) {
                    // Keep track of the number of signed crossings (see above).
                    val isLine = inputDimensions[b_input_edges[bi].inputId] == 1.toByte()
                    val sign = if(isLine) 0 else if(leftToRight == invert_b) -1 else 1
                    a_num_crossings[aIndex] += sign

                    // Any polyline or polygon vertex that has at least one crossing but no
                    // adjacent emitted edge may be emitted as an isolated vertex.
                    a_isolated[aIndex] = true
                } else {
                    // TODO(b/112043775): fix this condition.
                    logger.error { "Failed to get crossed vertex index." }
                }
            }

            // Finally, we iterate through the A edge chain, keeping track of the
            // number of signed crossings as we go along.  The "multiplicity" is
            // defined as the cumulative number of signed crossings, and indicates how
            // many edges should be output (and in which direction) in order to link
            // up the edge crossings in the correct order.  (The multiplicity is
            // almost always either 0 or 1 except in very rare cases.)
            var multiplicity = (if(inside) 1 else 0) + a_num_crossings[0]
            for (ai in 1 until a_vertices.size) {
                if (multiplicity != 0) {
                    a_isolated[ai - 1] = false
                    a_isolated[ai] = false
                }
                val edgeCount = if(reverse_a) -multiplicity else multiplicity
                // Output any forward edges required.
                if (edgeCount > 0) repeat(edgeCount) {
                    addEdge(GraphEdge(a_vertices[ai - 1], a_vertices[ai]), aInputId)
                }
                // Output any reverse edges required.
                else repeat (-edgeCount) {
                    addEdge(GraphEdge(a_vertices[ai], a_vertices[ai - 1]), aInputId)
                }
                multiplicity += a_num_crossings[ai]
            }
            // Multiplicities other than 0 or 1 can only occur in the edge interior.
            assert(multiplicity == 0 || multiplicity == 1)
            inside = (multiplicity != 0)

            // Output any isolated polyline vertices.
            // TODO(ericv): Only do this if an output layer wants degenerate edges.
            if (inputDimensions[aInputId] != 0.toByte()) {
                for (ai in 0 until a_vertices.size) {
                    if (a_isolated[ai]) {
                        addEdge(GraphEdge(a_vertices[ai], a_vertices[ai]), aInputId)
                    }
                }
            }

            ++i
        }
    }

    private fun addEdge(edge: dilivia.s2.builder.Edge, inputEdgeId: InputEdgeId) {
        newEdges.add(edge)
        newInputEdgeIds.add(inputEdgeId)
    }

    // Given the vertices of the snapped edge chain for an input edge A and the
    // set of input edges B that cross input edge A, this method gathers all of
    // the snapped edges of B that are incident to a given snapped vertex of A.
    // The incident edges for each input edge of B are appended to a separate
    // output vector.  (A and B can refer to either the input edge or the
    // corresponding snapped edge chain.)
    private fun gatherIncidentEdges(a: List<VertexId>, ai: Int, bInputEdges: List<CrossingInputEdge>, bEdges: List<CrossingGraphEdgeVector>) {
        // Examine all of the edges incident to the given vertex of A.  If any edge
        // comes from a B input edge, append it to the appropriate vector.
        requireEQ(bInputEdges.size, bEdges.size)

        logger.trace { "gatherIncidentEdges | a = $a\nai = $ai\nbInputEdges = $bInputEdges" }

        for (e:EdgeId in inMap.edgeIds(a[ai])) {
            val id = inputIds[e]
            val idx = bInputEdges.lowerBound(value = id)
            logger.trace { "gatherIncidentEdges | process in edge $id. edges index = $idx" }
            if (idx != bInputEdges.size && bInputEdges[idx].inputId == id) {
                val edges = bEdges[idx]
                edges.add(CrossingGraphEdge(e, ai, false, g.edge(e).first))
            }
        }
        for (e: EdgeId in outMap.edgeIds(a[ai])) {
            val id = inputIds[e]
            val idx = bInputEdges.lowerBound(value = id)
            logger.trace { "gatherIncidentEdges | process out edge $id. edges index = $idx" }
            if (idx != bInputEdges.size && bInputEdges[idx].inputId == id) {
                val edges = bEdges[idx]
                edges.add(CrossingGraphEdge(e, ai, true, g.edge(e).second))
            }
        }

        logger.trace { "gatherIncidentEdges | bEdges = ${bEdges.joinToString(separator = "\n  - ", prefix = "[\n  - ", postfix = "\n]")}" }
    }

    // Given an edge chain A that is crossed by another edge chain B (where
    // "left_to_right" indicates whether B crosses A from left to right), this
    // method decides which vertex of A the crossing takes place at.  The
    // parameters are the vertices of the A chain ("a") and the set of edges in
    // the B chain ("b") that are incident to vertices of A.  The B chain edges
    // are sorted in increasing order of (a_index, outgoing) tuple.
    private fun getCrossedVertexIndex(a: List<VertexId>, b: CrossingGraphEdgeVector, leftToRight: Boolean): Int {
        logger.trace { "getCrossedVertexIndex | a = $a,\nb = $b,\nleftToRight = $leftToRight" }
        requireArgument { a.isNotEmpty() }
        requireArgument { b.isNotEmpty() }

        // The reason this calculation is tricky is that after snapping, the A and B
        // chains may meet and separate several times.  For example, if B crosses A
        // from left to right, then B may touch A, make an excursion to the left of
        // A, come back to A, then make an excursion to the right of A and come back
        // to A again, like this:
        //
        //  *--B--*-\             /-*-\
        //           B-\       /-B     B-\      6     7     8     9
        //  *--A--*--A--*-A,B-*--A--*--A--*-A,B-*--A--*--A--*-A,B-*
        //  0     1     2     3     4     5      \-B     B-/
        //                                          \-*-/
        //
        // (where "*" is a vertex, and "A" and "B" are edge labels).  Note that B
        // may also follow A for one or more edges whenever they touch (e.g. between
        // vertices 2 and 3).  In this case the only vertices of A where the
        // crossing could take place are 5 and 6, i.e. after all excursions of B to
        // the left of A, and before all excursions of B to the right of A.
        //
        // Other factors to consider are that the portion of B before and/or after
        // the crossing may be degenerate, and some or all of the B edges may be
        // reversed relative to the A edges.

        // First, check whether edge A is degenerate.
        val n = a.size
        if (n == 1) return 0

        // If edge chain B is incident to only one vertex of A, we're done.
        if (b[0].aIndex == b.last().aIndex) return b[0].aIndex

        // Determine whether the B chain visits the first and last vertices that it
        // shares with the A chain in the same order or the reverse order.  This is
        // only needed to implement one special case (see below).
        val b_reversed = getVertexRank(b[0]) > getVertexRank(b.last())

        // Examine each incident B edge and use it to narrow the range of positions
        // where the crossing could occur in the B chain.  Vertex positions are
        // represented as a range [lo, hi] of vertex ranks in the B chain (see
        // GetVertexRank).
        //
        // Note that if an edge of B is incident to the first or last vertex of A,
        // we can't test which side of the A chain it is on.  (An s2pred::Sign test
        // doesn't work; e.g. if the B edge is XY and the first edge of A is YZ,
        // then snapping can change the sign of XYZ while maintaining topological
        // guarantees.)  There can be up to 4 such edges (one incoming and one
        // outgoing edge at each endpoint of A).  Two of these edges logically
        // extend past the end of the A chain and place no restrictions on the
        // crossing vertex.  The other two edges define the ends of the subchain
        // where B shares vertices with A.  We save these edges in order to handle a
        // special case (see below).
        var lo = -1
        var hi = order.size   // Vertex ranks of acceptable crossings
        var b_first: Int = -1
        var b_last: Int = -1  // "b" subchain connecting "a" endpoints
        for (e in b) {
            val ai = e.aIndex
            if (ai == 0) {
                if (e.outgoing != b_reversed && e.dst != a[1]) b_first = e.id
            } else if (ai == n - 1) {
                if (e.outgoing == b_reversed && e.dst != a[n - 2]) b_last = e.id
            } else {
                // This B edge is incident to an interior vertex of the A chain.  First
                // check whether this edge is identical (or reversed) to an edge in the
                // A chain, in which case it does not create any restrictions.
                if (e.dst == a[ai - 1] || e.dst == a[ai + 1]) continue

                // Otherwise we can test which side of the A chain the edge lies on.
                val on_left = S2Predicates.orderedCCW(g.vertex(a[ai + 1]), g.vertex(e.dst), g.vertex(a[ai - 1]), g.vertex(a[ai]))

                // Every B edge that is incident to an interior vertex of the A chain
                // places some restriction on where the crossing vertex could be.
                if (leftToRight == on_left) {
                    // This is a pre-crossing edge, so the crossing cannot be before the
                    // destination vertex of this edge.  (For example, the input B edge
                    // crosses the input A edge from left to right and this edge of the B
                    // chain is to the left of the A chain.)
                    lo = max(lo, rank[e.id] + 1)
                } else {
                    // This is a post-crossing edge, so the crossing cannot be after the
                    // source vertex of this edge.
                    hi = min(hi, rank[e.id])
                }
            }
        }

        // There is one special case.  If a subchain of B connects the first and
        // last vertices of A, then together with the edges of A this forms a loop
        // whose orientation can be tested to determine whether B is on the left or
        // right side of A.  This is only possible (and only necessary) if the B
        // subchain does not include any interior vertices of A, since otherwise the
        // B chain might cross from one side of A to the other.
        //
        // Note that it would be possible to avoid this test in some situations by
        // checking whether either endpoint of the A chain has two incident B edges,
        // in which case we could check which side of the B chain the A edge is on
        // and use this to limit the possible crossing locations.
        if (b_first >= 0 && b_last >= 0) {
            // The B subchain connects the first and last vertices of A.  Test whether
            // the chain includes any interior vertices of A.  We do this indirectly
            // by testing whether any edge of B has restricted the range of allowable
            // crossing vertices (since any interior edge of the B subchain incident
            // to any interior edge of A is guaranteed to do so).
            var min_rank = order.size
            var max_rank = -1
            for (e in b) {
                min_rank = min(min_rank, getVertexRank(e))
                max_rank = max(max_rank, getVertexRank(e))
            }
            if (lo <= min_rank && hi >= max_rank) {
                // The B subchain is not incident to any interior vertex of A.
                // Swap the edges if necessary so that they are in B chain order.
                if (b_reversed) {
                    val tmp = b_first
                    b_first = b_last
                    b_last = tmp
                }
                val on_left = edgeChainOnLeft(a, b_first, b_last)
                if (leftToRight == on_left) {
                    lo = max(lo, rank[b_last] + 1)
                } else {
                    hi = min(hi, rank[b_first])
                }
            }
        }

        // Otherwise we choose the smallest shared VertexId in the acceptable range,
        // in order to ensure that both chains choose the same crossing vertex.
        var best = -1
        checkLE(lo, hi)
        for (e in b) {
            val ai = e.aIndex
            val vrank = getVertexRank(e)
            if (vrank in lo..hi && (best < 0 || a[ai] < a[best])) {
                best = ai
            }
        }
        return best
    }

    // Returns the "vertex rank" of the shared vertex associated with the given
    // CrossingGraphEdge.  Recall that graph edges are sorted in input edge order,
    // and that the rank of an edge is its position in this order (rank_[e]).
    // VertexRank(e) is defined such that VertexRank(e.src) == rank_[e] and
    // VertexRank(e.dst) == rank_[e] + 1.  Note that the concept of "vertex rank"
    // is only defined within a single edge chain (since different edge chains can
    // have overlapping vertex ranks).
    private fun getVertexRank(e: CrossingGraphEdge): Int = rank[e.id] + (if (!e.outgoing) 1 else 0)

    // Given edge chains A and B that form a loop (after possibly reversing the
    // direction of chain B), returns true if chain B is to the left of chain A.
    // Chain A is given as a sequence of vertices, while chain B is specified as
    // the first and last edges of the chain.
    private fun edgeChainOnLeft(a: List<VertexId>, b_first: EdgeId, b_last: EdgeId): Boolean {
        // Gather all the interior vertices of the B subchain.
        val loop = mutableListOf<VertexId>()
        for (i in rank[b_first] until rank[b_last]) {
            loop.add(g.edge(order[i]).second)
        }
        // Possibly reverse the chain so that it forms a loop when "a" is appended.
        if (g.edge(b_last).second != a[0]) loop.reverse()
        loop.addAll(a)
        // Duplicate the first two vertices to simplify vertex indexing.
        for (j in 0..1) {
            loop.add(loop[j])
        }
        // Now B is to the left of A if and only if the loop is counterclockwise.
        var sum = 0.0
        for (i in 2 until loop.size) {
            sum += turnAngle(g.vertex(loop[i - 2]), g.vertex(loop[i - 1]), g.vertex(loop[i]))
        }
        return sum > 0.0
    }

}

// Given a set of clipping instructions encoded as a set of intersections
// between input edges, EdgeClippingLayer determines which graph edges
// correspond to clipped portions of input edges and removes them.  It
// assembles the remaining edges into a new S2Builder::Graph and passes the
// result to the given output layer for assembly.
class EdgeClippingLayer(
        private val layers: List<Layer>,
        private val inputDimensions: List<Byte>,
        private val inputCrossings: InputEdgeCrossings
) : Layer() {

    // Layer interface:

    // We keep all edges, including degenerate ones, so that we can figure out
    // the correspondence between input edge crossings and output edge
    // crossings.
    override fun graphOptions(): GraphOptions {
        return GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.KEEP, SiblingPairs.KEEP)
    }

    override fun build(g: Graph, error: S2Error) {
        // The bulk of the work is handled by GraphEdgeClipper.
        val new_edges = ArrayList<GraphEdge>()
        val new_input_edge_ids = ArrayList<InputEdgeIdSetId>()
        // Destroy the GraphEdgeClipper immediately to save memory.
        GraphEdgeClipper(g, inputDimensions, inputCrossings, new_edges, new_input_edge_ids).run()
        logger.trace { """Edges after clipping:
            |${new_edges.mapIndexed { index, edge -> "${new_input_edge_ids[index]} (${edge.first}, ${edge.second})" }.joinToString("\n")}
            |""".trimMargin()
        }
        // Construct one or more graphs from the clipped edges and pass them to the
        // given output layer(s).
        val new_input_edge_id_set_lexicon = IdSetLexicon()
        if (layers.size == 1) {
            val options = layers[0].graphOptions()
            val new_graph = makeGraph (g, options, new_edges, new_input_edge_ids, new_input_edge_id_set_lexicon, error)
            layers[0].build(new_graph, error)
        } else {
            // The Graph objects must be valid until the last Build() call completes,
            // so we store all of the graph data in arrays with 3 elements.
            checkEQ(3, layers.size)
            val layer_edges = Array<ArrayList<GraphEdge>>(3) { ArrayList() }
            val layer_input_edge_ids = Array<ArrayList<InputEdgeIdSetId>>(3) { ArrayList() }
            val layer_options = ArrayList<GraphOptions>(3)
            val layer_graphs = ArrayList<Graph>(3)
            // Separate the edges according to their dimension.
            for (i in 0 until new_edges.size) {
                val d = inputDimensions[new_input_edge_ids[i]]
                layer_edges[d.toInt()].add(new_edges[i])
                layer_input_edge_ids[d.toInt()].add(new_input_edge_ids[i])
            }
            // Clear variables to save space.
            new_edges.clear()
            new_input_edge_ids.clear()
            for (d in 0..2) {
                layer_options.add(d, layers[d].graphOptions())
                layer_graphs.add(makeGraph(g, layer_options[d], layer_edges[d], layer_input_edge_ids[d], new_input_edge_id_set_lexicon, error));
                layers[d].build(layer_graphs[d], error)
            }
        }
    }

    companion object {

        private val logger = KotlinLogging.logger(EdgeClippingLayer::class.java.name)

        // Helper function (in anonymous namespace) to create an S2Builder::Graph from
        // a vector of edges.
        fun makeGraph(g: Graph, options: GraphOptions,
                      new_edges: ArrayList<GraphEdge>, new_input_edge_ids: ArrayList<InputEdgeIdSetId>,
                      new_input_edge_id_set_lexicon: IdSetLexicon,
                      error: S2Error
        ): Graph {
            if (options.edgeType == EdgeType.UNDIRECTED) {
                // Create a reversed edge for every edge.
                val n = new_edges.size
                new_edges.ensureCapacity(2 * n)
                new_input_edge_ids.ensureCapacity(2 * n)
                for (i in 0 until n) {
                    new_edges.add(Graph.reverse(new_edges[i]))
                    new_input_edge_ids.add(IdSetLexicon.emptySetId())
                }
            }
            Graph.processEdges(options, new_edges, new_input_edge_ids, new_input_edge_id_set_lexicon, error)
            return Graph(options, g.vertices,
                    new_edges, new_input_edge_ids, new_input_edge_id_set_lexicon,
                    g.labelSetIds, g.labelSetLexicon, g.isFullPolygonPredicate
            )
        }
    }
}

// An IndexCrossing represents a pair of intersecting S2ShapeIndex edges
// ("a_edge" and "b_edge").  We store all such intersections because the
// algorithm needs them twice, once when processing the boundary of region A
// and once when processing the boundary of region B.
data class IndexCrossing(
        var a: ShapeEdgeId,
        var b: ShapeEdgeId,

        // True if S2::CrossingSign(a_edge, b_edge) > 0.
        var isInteriorCrossing: Boolean = false,

        // True if "a_edge" crosses "b_edge" from left to right.  Undefined if
        // is_interior_crossing is false.
        var leftToRight: Boolean = false,

        // Equal to S2::VertexCrossing(a_edge, b_edge).  Undefined if "a_edge" and
        // "b_edge" do not share exactly one vertex or either edge is degenerate.
        var isVertexCrossing: Boolean = false,
) : Comparable<IndexCrossing> {

    override fun compareTo(y: IndexCrossing): Int {
        // The compiler (2017) doesn't optimize the following as well:
        // return x.a < y.a || (x.a == y.a && x.b < y.b);
        if (a.shapeId < y.a.shapeId) return -1
        if (y.a.shapeId < a.shapeId) return 1
        if (a.edgeId < y.a.edgeId) return -1
        if (y.a.edgeId < a.edgeId) return 1
        if (b.shapeId < y.b.shapeId) return -1
        if (y.b.shapeId < b.shapeId) return 1
        return b.edgeId.compareTo(y.b.edgeId)
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IndexCrossing) return false

        if (a != other.a) return false
        if (b != other.b) return false

        return true
    }

    override fun hashCode(): Int {
        var result = a.hashCode()
        result = 31 * result + b.hashCode()
        return result
    }
}

typealias IndexCrossings = MutableList<IndexCrossing>

typealias CrossingGraphEdgeVector = ArrayList<CrossingGraphEdge>

// InputEdgeCrossings represents all pairs of intersecting input edges.
// It is sorted in lexicographic order.
typealias InputEdgeCrossings = ArrayList<Pair<InputEdgeId, CrossingInputEdge>>

// SourceEdgeCrossing represents an input edge that crosses some other
// edge; it crosses the edge from left to right iff the second parameter
// is "true".
typealias SourceEdgeCrossing = Pair<SourceId, Boolean>
typealias SourceEdgeCrossings = MutableList<Pair<InputEdgeId, SourceEdgeCrossing>>
typealias SourceIdMap = TreeMap<SourceId, InputEdgeId>


// A helper class for iterating through the edges from region B that cross a
// particular edge from region A.  It caches information from the current
// shape, chain, and edge so that it doesn't need to be looked up repeatedly.
// Typical usage:
//
//  void SomeFunction(ShapeEdgeId a_id, CrossingIterator *it) {
//    // Iterate through the edges that cross edge "a_id".
//    for (; !it->Done(a_id); it->Next()) {
//      ... use it->b_shape(), it->b_edge(), etc ...
//    }
// Creates an iterator over crossing edge pairs (a, b) where "b" is an edge
// from "b_index".  "crossings_complete" indicates that "crossings" contains
// all edge crossings between the two regions (rather than a subset).
class CrossingIterator(
        val bIndex: S2ShapeIndex,
        private val crossings: IndexCrossings,
        // True if all edge crossings are available (see above).
        val crossingsComplete: Boolean
) {
    private var currentIdx = 0;
    private var current: IndexCrossing = crossings[currentIdx]

    private lateinit var b_shape: S2Shape
    private var b_shape_id: Int = -1
    private var b_dimension: Int = -1
    private val b_info: ChainInfo = ChainInfo()  // Computed on demand.

    init {
        update()
    }

    fun next() {
        ++currentIdx
        if (currentIdx < crossings.size) current = crossings[currentIdx]
        update()
    }

    fun done(id: ShapeEdgeId): Boolean = a_id() != id

    // True if this crossing occurs at a point interior to both edges.
    fun is_interior_crossing(): Boolean = current.isInteriorCrossing

    // Equal to S2::VertexCrossing(a_edge, b_edge), provided that a_edge and
    // b_edge have exactly one vertex in common and neither edge is degenerate.
    fun is_vertex_crossing(): Boolean = current.isVertexCrossing

    // True if a_edge crosses b_edge from left to right (for interior crossings).
    fun left_to_right(): Boolean = current.leftToRight

    fun a_id(): ShapeEdgeId = current.a

    fun b_id(): ShapeEdgeId = current.b

    fun b_shape(): S2Shape = b_shape

    fun b_dimension(): Int = b_dimension

    fun b_shape_id(): Int = b_shape_id

    fun b_edge_id(): Int = b_id().edgeId

    fun b_edge(): Edge = b_shape.edge(b_edge_id()) // Opportunity to cache this.

    // Information about the chain to which an edge belongs.
    data class ChainInfo(
        var chain_id: Int = -1,  // chain id
        var start: Int = -1,     // starting edge id
        var limit: Int = -1     // limit edge id
    )

    // Returns a description of the chain to which the current B edge belongs.
    fun b_chain_info(): ChainInfo {
        if (b_info.chain_id < 0) {
            b_info.chain_id = b_shape().chainPosition(b_edge_id()).chainId
            val chain = b_shape ().chain(b_info.chain_id)
            b_info.start = chain.start
            b_info.limit = chain.start + chain.length
        }
        return b_info
    }

    // Updates information about the B shape whenever it changes.
    private fun update() {
        if (a_id() != kSentinel && b_id().shapeId != b_shape_id) {
            b_shape_id = b_id().shapeId
            b_shape = bIndex.shape(b_shape_id)!!
            b_dimension = b_shape.dimension
            b_info.chain_id = -1  // Computed on demand.
        }
    }

}

// PointCrossingResult describes the relationship between a point from region A
// and a set of crossing edges from region B.  For example, "matches_polygon"
// indicates whether a polygon vertex from region B matches the given point.
data class PointCrossingResult(
        // Note that "matches_polyline" is true only if the point matches a polyline
        // vertex of B *and* the polyline contains that vertex, whereas
        // "matches_polygon" is true if the point matches any polygon vertex.
        var matchesPoint: Boolean = false,     // Matches point.
        var matchesPolyline: Boolean = false,  // Matches contained polyline vertex.
        var matchesPolygon: Boolean = false   // Matches polygon vertex.
)

// EdgeCrossingResult describes the relationship between an edge from region A
// ("a_edge") and a set of crossing edges from region B.  For example,
// "matches_polygon" indicates whether "a_edge" matches a polygon edge from
// region B.
data class EdgeCrossingResult(
        // These fields indicate that "a_edge" exactly matches an edge of B.
        var matchesPolyline: Boolean = false,     // Matches polyline edge (either direction).
        var matchesPolygon: Boolean = false,       // Matches polygon edge (same direction).
        var matchesSibling: Boolean = false,       // Matches polygon edge (reverse direction).

        // These fields indicate that a vertex of "a_edge" matches a polyline vertex
        // of B *and* the polyline contains that vertex.
        var a0MatchesPolyline: Boolean = false,    // Start vertex matches contained polyline vertex.
        var a1MatchesPolyline: Boolean = false,    // End vertex matches contained polyline vertex.

        // These fields indicate that a vertex of "a_edge" matches a polygon vertex
        // of B.  (Unlike with polylines, the polygon may not contain that vertex.)
        var a0MatchesPolygon: Boolean = false,     // Start vertex matches polygon vertex.
        var a1MatchesPolygon: Boolean = false,     // End vertex matches polygon vertex.

        // These fields count the number of edge crossings at the start vertex, end
        // vertex, and interior of "a_edge".
        var a0Crossings: Int = 0,          // Count of polygon crossings at start vertex.
        var a1Crossings: Int = 0,          // Count of polygon crossings at end vertex.
        var interiorCrossings: Int = 0,    // Count of polygon crossings in edge interior.
)
