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
package dilivia.s2.builder.graph

import dilivia.PreConditions
import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkLT
import dilivia.PreConditions.checkState
import dilivia.collections.assign
import dilivia.collections.erase
import dilivia.collections.get
import dilivia.collections.remove
import dilivia.collections.resize
import dilivia.collections.reverse
import dilivia.collections.rotate
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2TextParser
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.InputEdgeIdSetId
import dilivia.s2.builder.IsFullPolygonPredicate
import dilivia.s2.builder.IsFullPolygonUnspecified
import dilivia.s2.builder.Label
import dilivia.s2.builder.LabelSetId
import mu.KotlinLogging

/** Identifies a vertex in the graph.  Vertices are numbered sequentially starting from zero. */
typealias VertexId = Int

/** A loop consisting of a sequence of edges. */
typealias EdgeLoop = List<EdgeId>

typealias DirectedComponent = MutableList<EdgeLoop>

typealias UndirectedComponent = Pair<MutableList<EdgeLoop>, MutableList<EdgeLoop>>

typealias EdgePolyline = MutableList<EdgeId>

/**
 * A struct for sorting the incoming and outgoing edges around a vertex "v0".
 *
 * @property incoming Is this an incoming edge to "v0"?
 * @property index Index of this edge in "edges" or "inEdgeIds"
 * @property endpoint The other (not "v0") endpoint of this edge
 * @property rank Secondary key for edges with the same endpoint.
 */
data class VertexEdge(
        val incoming: Boolean,
        val index: EdgeId,
        val endpoint: VertexId,
        val rank: Int
)


// represents the outgoing edges from a given vertex.
typealias VertexOutEdges = List<Edge>

// represents the incoming edge *ids*  to a given vertex.
typealias VertexInEdgeIds = List<EdgeId>

typealias VertexOutEdgeIds = IntRange
fun VertexOutEdgeIds.size() = if (isEmpty()) 0 else (this.last - this.first + 1)

/**
 * An Graph represents a collection of snapped edges that is passed to a Layer for assembly. (Example layers include
 * polygons, polylines, and polygon meshes.) The Graph object does not own any of its underlying data it is simply a
 * view of data that is stored elsewhere. You will only need this interface if you want to implement a new Layer
 * subtype.
 *
 * The graph consists of vertices and directed edges. Vertices are numbered sequentially starting from zero. An edge is
 * represented as a pair of vertex ids. The edges are sorted in lexicographic order, therefore all of the outgoing
 * edges from a particular vertex form a contiguous range.
 *
// TODO(ericv): Consider pulling out the methods that are helper functions for Layer implementations (such as GetDirectedLoops) into s2builderutil_graph.h.
 */
class Graph(
        // "options":
        //    - the GraphOptions used to build the Graph.  In some cases these
        //      can be different than the options provided by the Layer.
    val options: GraphOptions = GraphOptions(),

        // "vertices":
        //   - a vector of S2Points indexed by VertexId.
    val vertices: List<S2Point> = emptyList(),

        // "edges":
        //   - a vector of VertexId pairs (sorted in lexicographic order)
        //     indexed by EdgeId.
    val edges: List<Edge> = emptyList(),

        // "input_edge_id_set_ids":
        //   - a vector indexed by EdgeId that allows access to the set of
        //     InputEdgeIds that were mapped to the given edge, by looking up the
        //     returned value (an InputEdgeIdSetId) in "input_edge_id_set_lexicon".
    val inputEdgeIdSetIds: List<InputEdgeIdSetId> = emptyList(),

        // "input_edge_id_set_lexicon":
        //   - a class that maps an InputEdgeIdSetId to a set of InputEdgeIds.
    val inputEdgeIdSetLexicon: IdSetLexicon = IdSetLexicon(),

        // "label_set_ids":
        //   - a vector indexed by InputEdgeId that allows access to the set of
        //     labels that were attached to the given input edge, by looking up the
        //     returned value (a LabelSetId) in the "label_set_lexicon".
    val labelSetIds: List<LabelSetId> = emptyList(),

        // "label_set_lexicon":
        //   - a class that maps a LabelSetId to a set of S2Builder::Labels.
    val labelSetLexicon: IdSetLexicon = IdSetLexicon(),

        // "is_full_polygon_predicate":
        //   - a predicate called to determine whether a graph consisting only of
        //     polygon degeneracies represents the empty polygon or the full polygon
        //     (see s2builder.h for details).
    val isFullPolygonPredicate: IsFullPolygonPredicate = IsFullPolygonUnspecified(),

    ) {

    constructor(g: Graph): this(
        options = g.options.copy(),
        vertices = g.vertices.map { it.clone() },
        edges = g.edges.map { it.copy() },
        inputEdgeIdSetIds = g.inputEdgeIdSetIds.toList(),
        inputEdgeIdSetLexicon = g.inputEdgeIdSetLexicon.clone(),
        labelSetIds = g.labelSetIds.toList(),
        labelSetLexicon = g.labelSetLexicon.clone(),
        isFullPolygonPredicate = g.isFullPolygonPredicate
    )

    // the number of vertices in the graph.
    val numVertices: Int
        get() = vertices.size

    //  the total number of edges in the graph.
    val numEdges: Int
        get() = edges.size

    // Returns the vertex at the given index.
    fun vertex(v: VertexId): S2Point = vertices[v]

    // Returns the endpoints of the given edge (as vertex indices).
    fun edge(e: EdgeId): Edge = edges[e]

    // Returns a vector of edge ids sorted in lexicographic order by
    // (destination, origin).  All of the incoming edges to each vertex form a
    // contiguous subrange of this ordering.
    fun getInEdgeIds(): MutableList<EdgeId> {
        val inEdgeIds = ArrayList<EdgeId>(numEdges)
        repeat(numEdges) { i -> inEdgeIds.add(i) }
        inEdgeIds.sortWith { ai, bi -> if (stableLessThan(reverse(edge(ai)), reverse(edge(bi)), ai, bi)) -1 else 1 }
        return inEdgeIds
    }

    // Given a graph such that every directed edge has a sibling, returns a map
    // from EdgeId to the sibling EdgeId.  This method is identical to
    // GetInEdgeIds() except that (1) it requires edges to have siblings, and
    // (2) undirected degenerate edges are grouped together in pairs such that
    // one edge is the sibling of the other.  Handles duplicate edges correctly
    // and is also consistent with GetLeftTurnMap().
    //
    // REQUIRES: An option is chosen that guarantees sibling pairs:
    //     (options.sibling_pairs() == { REQUIRE, CREATE } ||
    //      options.edge_type() == UNDIRECTED)
    fun getSiblingMap(): List<EdgeId> {
        val inEdgeIds = getInEdgeIds()
        makeSiblingMap(inEdgeIds)
        return inEdgeIds
    }

    // Like GetSiblingMap(), but constructs the map starting from the vector of
    // incoming edge ids returned by GetInEdgeIds().  (This operation is a no-op
    // except unless undirected degenerate edges are present, in which case such
    // edges are grouped together in pairs to satisfy the requirement that every
    // edge must have a sibling edge.)
    fun makeSiblingMap(inEdgeIds: MutableList<EdgeId>) {
        checkState {
            (options.siblingPairs == SiblingPairs.REQUIRE ||
                    options.siblingPairs == SiblingPairs.CREATE ||
                    options.edgeType == EdgeType.UNDIRECTED)
        }
        if (PreConditions.enabled) {
            repeat(numEdges) { e -> checkState { (edge(e) == reverse(edge(inEdgeIds[e]))) } }
        }
        if (options.edgeType == EdgeType.DIRECTED) return
        if (options.degenerateEdges == DegenerateEdges.DISCARD) return

        var e = 0
        while (e < numEdges) {
            val v = edge(e).first
            if (edge(e).second == v) {
                checkLT(e + 1, numEdges)
                checkEQ(edge(e + 1).first, v)
                checkEQ(edge(e + 1).second, v)
                checkEQ(inEdgeIds[e], e)
                checkEQ(inEdgeIds[e + 1], e + 1)
                inEdgeIds[e] = e + 1
                inEdgeIds[e + 1] = e
                ++e
            }
            ++e
        }
    }

    // Returns the set of input edge ids that were snapped to the given
    // edge.  ("Input edge ids" are assigned to input edges sequentially in
    // the order they are added to the builder.)  For example, if input
    // edges 2 and 17 were snapped to edge 12, then input_edge_ids(12)
    // returns a set containing the numbers 2 and 17.  Example usage:
    //
    //   for (InputEdgeId input_edge_id : g.input_edge_ids(e)) { ... }
    //
    // Please note the following:
    //
    //  - When edge chains are simplified, the simplified edge is assigned all
    //    the input edge ids associated with edges of the chain.
    //
    //  - Edges can also have multiple input edge ids due to edge merging
    //    (if DuplicateEdges::MERGE is specified).
    //
    //  - Siblings edges automatically created by EdgeType::UNDIRECTED or
    //    SiblingPairs::CREATE have an empty set of input edge ids.  (However
    //    you can use a LabelFetcher to retrieve the set of labels associated
    //    with both edges of a given sibling pair.)
    fun inputEdgeIds(e: EdgeId): IdSetLexicon.IdSet = inputEdgeIdSetLexicon.idSet(inputEdgeIdSetIds[e])

    // Low-level method that returns an integer representing the entire set of
    // input edge ids that were snapped to the given edge.  The elements of the
    // IdSet can be accessed using input_edge_id_set_lexicon().
    fun inputEdgeIdSetId(e: EdgeId): InputEdgeIdSetId = inputEdgeIdSetIds[e]

    // Returns the minimum input edge id that was snapped to this edge, or -1 if
    // no input edges were snapped (see SiblingPairs.CREATE).  This is
    // useful for layers that wish to preserve the input edge ordering as much
    // as possible (e.g., to ensure idempotency).
    fun minInputEdgeId(e: EdgeId): InputEdgeId {
        val idSet = inputEdgeIds(e)
        val minEdgeId = if (idSet.size() == 0) kNoInputEdgeId else idSet.values.first()
        logger.trace { "minInputEdgeId(e = $e) => $minEdgeId" }
        return minEdgeId
    }

    // Returns a vector containing the minimum input edge id for every edge.
    // If an edge has no input ids, kNoInputEdgeId is used.
    fun getMinInputEdgeIds(): List<InputEdgeId> {
        val minInputIds = ArrayList<InputEdgeId>(numEdges)
        repeat(numEdges) { e -> minInputIds.add(minInputEdgeId(e)) }
        return minInputIds
    }

    // Returns a vector of EdgeIds sorted by minimum input edge id.  This is an
    // approximation of the input edge ordering.x
    fun getInputEdgeOrder(minInputEdgeIds: List<InputEdgeId>): MutableList<EdgeId> {
        val order = mutableListOf<EdgeId>()
        minInputEdgeIds.indices.forEach { i -> order.add(i) }
        order.sortWith { a: EdgeId, b: EdgeId ->
            if (minInputEdgeIds[a] == minInputEdgeIds[b]) a.compareTo(b)
            else minInputEdgeIds[a].compareTo(minInputEdgeIds[b])
        }
        return order
    }

    // Returns the set of labels associated with a given input edge.  Example:
    //   for (Label label : g.labels(input_edge_id)) { ... }
    fun labels(e: InputEdgeId): IdSetLexicon.IdSet = labelSetLexicon.idSet(labelSetIds[e])

    // Low-level method that returns an integer representing the set of
    // labels associated with a given input edge.  The elements of
    // the IdSet can be accessed using label_set_lexicon().
    fun labelSetId(e: InputEdgeId): LabelSetId = labelSetIds[e]

    // Convenience method that calls is_full_polygon_predicate() to determine
    // whether a graph that consists only of polygon degeneracies represents the
    // empty polygon or the full polygon (see s2builder.h for details).
    fun isFullPolygon(error: S2Error): Boolean = isFullPolygonPredicate.test(this, error)

    // Returns a map "m" that maps each edge e=(v0,v1) to the following outgoing edge around "v1" in clockwise order.
    // (This corresponds to making a "left turn" at the vertex.)  By starting at a given edge and making only left
    // turns, you can construct a loop whose interior does not contain any edges in the same connected component.
    //
    // If the incoming and outgoing edges around a vertex do not alternate perfectly (e.g., there are two incoming
    // edges in a row), then adjacent (incoming, outgoing) pairs are repeatedly matched and removed.  This is
    // similar to finding matching parentheses in a string such as "(()())()".
    //
    // For sibling edge pairs, the incoming edge is assumed to immediately
    // follow the outgoing edge in clockwise order.  Thus a left turn is made
    // from an edge to its sibling only if there are no other outgoing edges.
    // With respect to the parentheses analogy, a sibling pair is ")(".
    // Similarly, if there are multiple copies of a sibling edge pair then the
    // duplicate incoming and outgoing edges are sorted in alternating order
    // (e.g., ")()(").
    //
    // Degenerate edges (edges from a vertex to itself) are treated as loops
    // consisting of a single edge.  This avoids the problem of deciding the
    // connectivity and ordering of such edges when they share a vertex with
    // other edges (possibly including other degenerate edges).
    //
    // If it is not possible to make a left turn from every input edge, this
    // method returns false and sets "error" appropriately.  In this situation
    // the left turn map is still valid except that any incoming edge where it
    // is not possible to make a left turn will have its entry set to -1.
    //
    // "in_edge_ids" should be equal to GetInEdgeIds() or GetSiblingMap().
    fun getLeftTurnMap(inEdgeIds: List<EdgeId>, leftTurnMap: MutableList<EdgeId>, error: S2Error): Boolean {
        logger.trace { "getLeftTurnMap: $inEdgeIds" }
        leftTurnMap.assign(numEdges, -1)
        if (numEdges == 0) return true

        // Declare vectors outside the loop to avoid reallocating them each time.
        val v0Edges = mutableListOf<VertexEdge>()
        val eIn = mutableListOf<EdgeId>()
        val eOut = mutableListOf<EdgeId>()

        // Walk through the two sorted arrays of edges (outgoing and incoming) and
        // gather all the edges incident to each vertex.  Then we sort those edges
        // and add an entry to the left turn map from each incoming edge to the
        // immediately following outgoing edge in clockwise order.
        var outIdx = 0
        var inIdx = 0
        var outEdge = edge(outIdx)
        var inEdge = edge(inEdgeIds[inIdx])
        val sentinel = Edge(numVertices, numVertices)
        var minEdge = minOf(outEdge, reverse(inEdge))
        while (minEdge != sentinel) {

            // Gather all incoming and outgoing edges around vertex "v0".
            val v0: VertexId = minEdge.first
            while (minEdge.first == v0) {
                val v1: VertexId = minEdge.second
                // Count the number of copies of "min_edge" in each direction.
                val outBegin = outIdx
                var inBegin = inIdx

                while (outEdge == minEdge) {
                    outEdge = if (++outIdx == numEdges) sentinel else edge(outIdx)
                }
                while (reverse(inEdge) == minEdge) {
                    inEdge = if (++inIdx == numEdges) sentinel else edge(inEdgeIds[inIdx])
                }

                if (v0 != v1) {
                    addVertexEdges(outBegin, outIdx, inBegin, inIdx, v1, v0Edges)
                } else {
                    // Each degenerate edge becomes its own loop.
                    while (inBegin < inIdx) {
                        leftTurnMap[inBegin] = inBegin
                        ++inBegin
                    }
                }

                minEdge = minOf(outEdge, reverse(inEdge))
            }
            if (v0Edges.isEmpty()) continue

            // Sort the edges in clockwise order around "v0".
            val minEndpoint: VertexId = v0Edges.first().endpoint
            v0Edges.subList(1, v0Edges.size).sortWith { a: VertexEdge, b: VertexEdge ->
                when {
                    a.endpoint == b.endpoint -> if (a.rank < b.rank) -1 else 1
                    a.endpoint == minEndpoint -> -1
                    b.endpoint == minEndpoint -> 1
                    !S2Predicates.orderedCCW(vertex(a.endpoint), vertex(b.endpoint), vertex(minEndpoint), vertex(v0)) -> -1
                    else -> 1
                }
            }

            // Match incoming with outgoing edges.  We do this by keeping a stack of
            // unmatched incoming edges.  We also keep a stack of outgoing edges with
            // no previous incoming edge, and match these at the end by wrapping
            // around circularly to the start of the edge ordering.
            for (e in v0Edges) {
                when {
                    e.incoming -> eIn.add(inEdgeIds[e.index])
                    eIn.isNotEmpty() -> {
                        leftTurnMap[eIn.last()] = e.index
                        eIn.removeLast()
                    }
                    else -> eOut.add(e.index)  // Matched below.
                }
            }
            // Pair up additional edges using the fact that the ordering is circular.
            eOut.reverse()
            while (eOut.isNotEmpty() && eIn.isNotEmpty()) {
                leftTurnMap[eIn.last()] = eOut.last()
                eOut.removeLast()
                eIn.removeLast()
            }
            // We only need to process unmatched incoming edges, since we are only
            // responsible for creating left turn map entries for those edges.
            if (eIn.isNotEmpty() && error.isOk()) {
                error.code = S2Error.BUILDER_EDGES_DO_NOT_FORM_LOOPS
                error.text = "Given edges do not form loops (indegree != outdegree)"
            }
            eIn.clear()
            eOut.clear()
            v0Edges.clear()
        }

        logger.trace { "getLeftTurnMap: result = $leftTurnMap, error = $error" }
        return error.isOk()
    }

    // Builds loops from a set of directed edges, turning left at each vertex
    // until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
    // (for LoopType::CIRCUIT) is found.  (Use LoopType::SIMPLE if you intend to
    // construct an S2Loop.)
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // If any degenerate edges are present, then each such edge is treated as a
    // separate loop.  This is mainly useful in conjunction with
    // options.degenerate_edges() == DISCARD_EXCESS, in order to build polygons
    // that preserve degenerate geometry.
    //
    // REQUIRES: options.degenerate_edges() == {DISCARD, DISCARD_EXCESS}
    // REQUIRES: options.edge_type() == DIRECTED
    fun getDirectedLoops(loopType: LoopType, loops: MutableList<EdgeLoop>, error: S2Error): Boolean {
        logger.trace { "getDirectedLoops($loopType: LoopType, loops: MutableList<EdgeLoop>, error: S2Error)" }
        checkState { options.degenerateEdges == DegenerateEdges.DISCARD || options.degenerateEdges == DegenerateEdges.DISCARD_EXCESS }
        checkState { options.edgeType == EdgeType.DIRECTED }

        val leftTurnMap = mutableListOf<EdgeId>()
        if (!getLeftTurnMap(getInEdgeIds(), leftTurnMap, error)) return false
        val minInputIds = getMinInputEdgeIds()
        logger.trace { "numEdges = $numEdges, leftTurnMap = $leftTurnMap, minInputIds = $minInputIds" }

        // If we are breaking loops at repeated vertices, we maintain a map from
        // VertexId to its position in "path".
        val pathIndex = mutableListOf<Int>()
        if (loopType == LoopType.SIMPLE) repeat(numVertices) { pathIndex.add(-1) }

        // Visit edges in arbitrary order, and try to build a loop from each edge.
        val path = mutableListOf<EdgeId>()
        for (start in 0 until numEdges) {
            logger.trace { "Visit edge $start" }
            if (leftTurnMap[start] < 0) {
                logger.trace { "Left turn of edge $start is negative: skip" }
                continue
            }

            // Build a loop by making left turns at each vertex until we return to
            // "start".  We use "left_turn_map" to keep track of which edges have
            // already been visited by setting its entries to -1 as we go along.  If
            // we are building vertex cycles, then whenever we encounter a vertex that
            // is already part of the path, we "peel off" a loop by removing those
            // edges from the path so far.
            var e = start
            var next: EdgeId
            while (leftTurnMap[e] >= 0) {
                path.add(e)
                next = leftTurnMap[e]
                leftTurnMap[e] = -1
                if (loopType == LoopType.SIMPLE) {
                    pathIndex[edge(e).first] = path.size - 1
                    val loopStart = pathIndex[edge(e).second]
                    if (loopStart < 0) {
                        e = next
                        continue
                    }
                    // Peel off a loop from the path.
                    val loop = path.subList(loopStart, path.size).toMutableList()
                    path.remove(loopStart, path.size)
                    for (e2 in loop) pathIndex[edge(e2).first] = -1
                    canonicalizeLoopOrder(minInputIds, loop)
                    loops.add(loop)
                }
                e = next
            }
            if (loopType == LoopType.SIMPLE) {
                check(path.isEmpty())  // Invariant.
            } else {
                canonicalizeLoopOrder(minInputIds, path)
                loops.add(path.toList())
                path.clear()
            }
        }
        canonicalizeVectorOrder(minInputIds, loops)
        logger.trace { "loops = $loops" }
        return true
    }

    // Builds loops from a set of directed edges, turning left at each vertex
    // until a repeated edge is found (i.e., LoopType::CIRCUIT).  The loops are
    // further grouped into connected components, where each component consists
    // of one or more loops connected by shared vertices.
    //
    // This method is used to build polygon meshes from directed or undirected
    // input edges.  To convert the output of this method into a mesh, the
    // client must determine how the loops in different components are related
    // to each other: for example, several loops from different components may
    // bound the same region on the sphere, in which case all of those loops are
    // combined into a single polygon.  (See s2shapeutil::BuildPolygonBoundaries
    // and s2builderutil::LaxPolygonVectorLayer for details.)
    //
    // Note that loops may include both edges of a sibling pair.  When several
    // such edges are connected in a chain or a spanning tree, they form a
    // zero-area "filament".  The entire loop may be a filament (i.e., a
    // degenerate loop with an empty interior), or the loop may have have
    // non-empty interior with several filaments that extend inside it, or the
    // loop may consist of several "holes" connected by filaments.  These
    // filaments do not change the interior of any loop, so if you are only
    // interested in point containment then they can safely be removed by
    // setting the "degenerate_boundaries" parameter to DISCARD.  (They can't be
    // removed by setting (options.sibling_pairs() == DISCARD) because the two
    // siblings might belong to different polygons of the mesh.)  Note that you
    // can prevent multiple copies of sibling pairs by specifying
    // options.duplicate_edges() == MERGE.
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
    //           (but requires DISCARD if degenerate_boundaries == DISCARD)
    // REQUIRES: options.sibling_pairs() == { REQUIRE, CREATE }
    //           [i.e., every edge must have a sibling edge]
    fun getDirectedComponents(degenerate_boundaries: DegenerateBoundaries, components: MutableList<DirectedComponent>, error: S2Error): Boolean {
        assert(options.degenerateEdges == DegenerateEdges.DISCARD ||
                (options.degenerateEdges == DegenerateEdges.DISCARD_EXCESS && degenerate_boundaries == DegenerateBoundaries.KEEP))
        assert(options.siblingPairs == SiblingPairs.REQUIRE || options.siblingPairs == SiblingPairs.CREATE)
        assert(options.edgeType == EdgeType.DIRECTED)  // Implied by above.

        val siblingMap = getInEdgeIds()
        val leftTurnMap = mutableListOf<EdgeId>()
        if (!getLeftTurnMap(siblingMap, leftTurnMap, error)) return false
        makeSiblingMap(siblingMap)
        val minInputIds = getMinInputEdgeIds()
        val frontier = mutableListOf<EdgeId>()  // Unexplored sibling edges.

        // A map from EdgeId to the position of that edge in "path".  Only needed if
        // degenerate boundaries are being discarded.
        val pathIndex = mutableListOf<Int>()
        if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
            pathIndex.assign(numEdges, -1)
        }
        for (min_start: EdgeId in 0 until numEdges) {
            if (leftTurnMap[min_start] < 0) continue  // Already used.

            // Build a connected component by keeping a stack of unexplored siblings
            // of the edges used so far.
            val component: DirectedComponent = mutableListOf()
            frontier.add(min_start)
            while (frontier.isNotEmpty()) {
                val start: Int = frontier.removeLast()
                if (leftTurnMap[start] < 0) continue  // Already used.

                // Build a path by making left turns at each vertex until we return to
                // "start".  Whenever we encounter an edge that is a sibling of an edge
                // that is already on the path, we "peel off" a loop consisting of any
                // edges that were between these two edges.
                val path = mutableListOf<EdgeId>()
                var e: EdgeId = start
                var next: EdgeId
                while (leftTurnMap[e] >= 0) {
                    path.add(e)
                    next = leftTurnMap[e]
                    leftTurnMap[e] = -1
                    // If the sibling hasn't been visited yet, add it to the frontier.
                    val sibling: EdgeId = siblingMap[e]
                    if (leftTurnMap[sibling] >= 0) {
                        frontier.add(sibling)
                    }
                    if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
                        pathIndex[e] = path.size - 1
                        val siblingIndex = pathIndex[sibling]
                        if (siblingIndex < 0) {
                            e = next
                            continue
                        }

                        // Common special case: the edge and its sibling are adjacent, in
                        // which case we can simply remove them from the path and continue.
                        if (siblingIndex == path.size - 2) {
                            path.resize(siblingIndex)
                            // We don't need to update "path_index" for these two edges
                            // because both edges of the sibling pair have now been used.
                            e = next
                            continue
                        }
                        // Peel off a loop from the path.
                        val loop = path.subList(siblingIndex + 1, path.size - 1).toMutableList() // TODO check
                        path.erase(siblingIndex, path.size)
                        // Mark the edges that are no longer part of the path.
                        for (e2: EdgeId in loop) pathIndex[e2] = -1
                        canonicalizeLoopOrder(minInputIds, loop)
                        component.add(loop)
                    }

                    e = next
                }
                // Mark the edges that are no longer part of the path.
                if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
                    for (e2: EdgeId in path) pathIndex[e2] = -1
                }
                canonicalizeLoopOrder(minInputIds, path)
                component.add(path)
            }
            canonicalizeVectorOrder(minInputIds, component)
            components.add(component)
        }
        // Sort the components to correspond to the input edge ordering.
        components.sortWith { a, b -> minInputIds[a[0][0]].compareTo(minInputIds[b[0][0]]) }

        return true
    }

    // Builds loops from a set of undirected edges, turning left at each vertex
    // until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
    // (for LoopType::CIRCUIT) is found.  The loops are further grouped into
    // "components" such that all the loops in a component are connected by
    // shared vertices.  Finally, the loops in each component are divided into
    // two "complements" such that every edge in one complement is the sibling
    // of an edge in the other complement.  This corresponds to the fact that
    // given any set of non-crossing undirected loops, there are exactly two
    // possible interpretations of the region that those loops represent (where
    // one possibility is the complement of the other).  This method does not
    // attempt to resolve this ambiguity, but instead returns both possibilities
    // for each connected component and lets the client choose among them.
    //
    // This method is used to build single polygons.  (Use GetDirectedComponents
    // to build polygon meshes, even when the input edges are undirected.)  To
    // convert the output of this method into a polygon, the client must choose
    // one complement from each component such that the entire set of loops is
    // oriented consistently (i.e., they define a region such that the interior
    // of the region is always on the left).  The non-chosen complements form
    // another set of loops that are also oriented consistently but represent
    // the complementary region on the sphere.  Finally, the client needs to
    // choose one of these two sets of loops based on heuristics (e.g., the area
    // of each region), since both sets of loops are equally valid
    // interpretations of the input.
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
    // REQUIRES: options.edge_type() == UNDIRECTED
    // REQUIRES: options.siblings_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
    //           [since REQUIRE, CREATE convert the edge_type() to DIRECTED]
    fun getUndirectedComponents(loopType: LoopType, components: MutableList<UndirectedComponent>, error: S2Error): Boolean {
        assert(options.degenerateEdges == DegenerateEdges.DISCARD || options.degenerateEdges == DegenerateEdges.DISCARD_EXCESS)
        assert(options.edgeType == EdgeType.UNDIRECTED)

        val siblingMap = getInEdgeIds()
        val leftTurnMap = mutableListOf<EdgeId>()
        if (!getLeftTurnMap(siblingMap, leftTurnMap, error)) return false
        makeSiblingMap(siblingMap)
        val minInputIds = getMinInputEdgeIds()

        // A stack of unexplored sibling edges.  Each sibling edge has a "slot"
        // (0 or 1) that indicates which of the two complements it belongs to.
        val frontier = mutableListOf<Pair<EdgeId, Int>>()

        // If we are breaking loops at repeated vertices, we maintain a map from
        // VertexId to its position in "path".
        val pathIndex = mutableListOf<Int>()
        if (loopType == LoopType.SIMPLE) pathIndex.assign(numVertices, -1)

        for (min_start: EdgeId in 0 until numEdges) {
            if (leftTurnMap[min_start] < 0) continue  // Already used.

            // Build a connected component by keeping a stack of unexplored siblings
            // of the edges used so far.
            var component = UndirectedComponent(mutableListOf(), mutableListOf())
            frontier.add(Pair(min_start, 0))
            while (frontier.isNotEmpty()) {
                val start: Int = frontier.last().first
                val slot = frontier.last().second
                frontier.removeLast()
                if (leftTurnMap[start] < 0) continue  // Already used.

                // Build a path by making left turns at each vertex until we return to
                // "start".  We use "left_turn_map" to keep track of which edges have
                // already been visited, and which complement they were assigned to, by
                // setting its entries to negative values as we go along.
                val path = mutableListOf<EdgeId>()
                var next: EdgeId
                var e: EdgeId = start
                while (leftTurnMap[e] >= 0) {
                    path.add(e)
                    next = leftTurnMap[e]
                    leftTurnMap[e] = markEdgeUsed(slot)
                    // If the sibling hasn't been visited yet, add it to the frontier.
                    val sibling: Int = siblingMap[e]
                    if (leftTurnMap[sibling] >= 0) {
                        frontier.add(Pair(sibling, 1 - slot))
                    } else if (leftTurnMap[sibling] != markEdgeUsed(1 - slot)) {
                        // Two siblings edges can only belong the same complement if the
                        // given undirected edges do not form loops.
                        error.init(S2Error.BUILDER_EDGES_DO_NOT_FORM_LOOPS, "Given undirected edges do not form loops")
                        return false
                    }
                    if (loopType == LoopType.SIMPLE) {
                        // Whenever we encounter a vertex that is already part of the path,
                        // we "peel off" a loop by removing those edges from the path.
                        pathIndex[edge(e).first] = path.size - 1
                        val loopStart = pathIndex[edge(e).second]
                        if (loopStart < 0) {
                            e = next
                            continue
                        }
                        val loop = path.subList(loopStart, path.size).toMutableList()
                        path.erase(loopStart, path.size)
                        // Mark the vertices that are no longer part of the path.
                        for (e2: EdgeId in loop) pathIndex[edge(e2).first] = -1
                        canonicalizeLoopOrder(minInputIds, loop)
                        component[slot].add(loop)
                    }

                    e = next
                }
                if (loopType == LoopType.SIMPLE) {
                    assert(path.isEmpty())  // Invariant.
                } else {
                    canonicalizeLoopOrder(minInputIds, path)
                    component[slot].add(path)
                }
            }
            canonicalizeVectorOrder(minInputIds, component[0])
            canonicalizeVectorOrder(minInputIds, component[1])
            // To save some work in S2PolygonLayer, we swap the two loop sets of the
            // component so that the loop set whose first loop most closely follows
            // the input edge ordering is first.  (If the input was a valid S2Polygon,
            // then this component will contain normalized loops.)
            if (minInputIds[component[0][0][0]] > minInputIds[component[1][0][0]]) {
                component = component.reverse()
            }
            components.add(component)
        }
        // Sort the components to correspond to the input edge ordering.
        components.sortWith { a: UndirectedComponent, b: UndirectedComponent -> minInputIds[a[0][0][0]].compareTo(minInputIds[b[0][0][0]]) }

        return true
    }

    /**
     * Builds polylines from a set of edges. If "polylineType" is PATH, then only vertices of indegree and outdegree 1
     * (or degree 2 in the case of undirected edges) will appear in the interior of polylines. This essentially
     * generates one polyline for each edge chain in the graph.  If "polyline_type" is WALK, then polylines may pass
     * through the same vertex or even the same edge multiple times (if duplicate edges are present), and each polyline
     * will be as long as possible.  This option is useful for reconstructing a polyline that has been snapped to a
     * lower resolution, since snapping can cause edges to become identical.
     *
     * This method attempts to preserve the input edge ordering in order to implement idempotency, even when there are
     * repeated edges or loops.  This is true whether directed or undirected edges are used.  Degenerate edges are also
     * handled appropriately.
     *
     * REQUIRES: options.sibling_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
     *
     * @param polylineType The type of polyline (PATH or WALK)
     * @return All the polylines built from the set of edges.
     */
    fun getPolylines(polylineType: PolylineType): List<EdgePolyline> {
        logger.trace { "getPolylines | polylineType = $polylineType" }
        checkState { options.siblingPairs in listOf(SiblingPairs.DISCARD, SiblingPairs.DISCARD_EXCESS, SiblingPairs.KEEP) }
        val builder = PolylineBuilder(this)
        return if (polylineType == PolylineType.PATH) {
            builder.buildPaths()
        } else {
            builder.buildWalks()
        }
    }


    // Convenience class to return the set of labels associated with a given
    // graph edge.  Note that due to snapping, one graph edge may correspond to
    // several different input edges and will have all of their labels.
    // This class is the preferred way to retrieve edge labels.
    //
    // The reason this is a class rather than a graph method is because for
    // undirected edges, we need to fetch the labels associated with both
    // siblings.  This is because only the original edge of the sibling pair has
    // labels; the automatically generated sibling edge does not.
    class LabelFetcher() {

        private lateinit var g: Graph
        private lateinit var edgeType: EdgeType
        private lateinit var siblingMap: List<EdgeId>

        constructor(g: Graph, edgeType: EdgeType) : this() {
            init(g, edgeType)
        }

        // Prepares to fetch labels associated with the given edge type.  For
        // EdgeType::UNDIRECTED, labels associated with both edges of the sibling
        // pair will be returned.  "edge_type" is a parameter (rather than using
        // g.options().edge_type()) so that clients can explicitly control whether
        // labels from one or both siblings are returned.
        fun init(g: Graph, edgeType: EdgeType) {
            this.g = g
            this.edgeType = edgeType
            if (edgeType == EdgeType.UNDIRECTED) siblingMap = g.getSiblingMap()
        }

        // Returns the set of labels associated with edge "e" (and also the labels
        // associated with the sibling of "e" if edge_type() is UNDIRECTED).
        // Labels are sorted and duplicate labels are automatically removed.
        //
        // This method uses an output parameter rather than returning by value in
        // order to avoid allocating a new vector on every call to this method.
        fun fetch(e: EdgeId, labels: MutableList<Label>) {
            labels.clear()
            for (input_edge_id in g.inputEdgeIds(e)) {
                for (label in g.labels(input_edge_id)) {
                    labels.add(label)
                }
            }
            if (edgeType == EdgeType.UNDIRECTED) {
                for (input_edge_id in g.inputEdgeIds(siblingMap[e])) {
                    for (label in g.labels(input_edge_id)) {
                        labels.add(label)
                    }
                }
            }
            if (labels.size > 1) {
                labels.sortAndRemoveDuplicates()
            }
        }

    }

    // Indicates whether loops should be simple cycles (no repeated vertices) or
    // circuits (which allow repeated vertices but not repeated edges).  In
    // terms of how the loops are built, this corresponds to closing off a loop
    // at the first repeated vertex vs. the first repeated edge.
    enum class LoopType { SIMPLE, CIRCUIT }

    enum class DegenerateBoundaries { DISCARD, KEEP }

    // Indicates whether polylines should be "paths" (which don't allow
    // duplicate vertices, except possibly the first and last vertex) or
    // "walks" (which allow duplicate vertices and edges).
    enum class PolylineType { PATH, WALK }

    fun toDebugString(): String {
        return """Graph:
            |vertices: ${vertices.map { v -> S2TextParser.toString(v) }}
            |edges: $edges
        """.trimMargin()
    }

    companion object {

        private val logger = KotlinLogging.logger(Graph::class.java.name)

        // Defines a value larger than any valid InputEdgeId.
        val kMaxInputEdgeId: InputEdgeId = Integer.MAX_VALUE

        // The following value of InputEdgeId means that an edge does not
        // corresponds to any input edge.
        val kNoInputEdgeId: InputEdgeId = kMaxInputEdgeId - 1

        // Given an edge (src, dst), returns the reverse edge (dst, src).
        fun reverse(e: Edge): Edge = Edge(e.second, e.first)

        // Rotates the edges of "loop" if necessary so that the edge(s) with the
        // largest input edge ids are last.  This ensures that when an output loop
        // is equivalent to an input loop, their cyclic edge orders are the same.
        // "min_input_ids" is the output of GetMinInputEdgeIds().
        fun canonicalizeLoopOrder(minInputIds: List<InputEdgeId>, loop: MutableList<EdgeId>) {
            if (loop.isEmpty()) return
            // Find the position of the element with the highest input edge id.  If
            // there are multiple such elements together (i.e., the edge was split
            // into several pieces by snapping it to several vertices), then we choose
            // the last such position in cyclic order (this attempts to preserve the
            // original loop order even when new vertices are added).  For example, if
            // the input edge id sequence is (7, 7, 4, 5, 6, 7) then we would rotate
            // it to obtain (4, 5, 6, 7, 7, 7).

            // The reason that we put the highest-numbered edge last, rather than the
            // lowest-numbered edge first, is that S2Loop::Invert() reverses the loop
            // edge order *except* for the last edge.  For example, the loop ABCD (with
            // edges AB, BC, CD, DA) becomes DCBA (with edges DC, CB, BA, AD).  Note
            // that the last edge is the same except for its direction (DA vs. AD).
            // This has the advantage that if an undirected loop is assembled with the
            // wrong orientation and later inverted (e.g. by S2Polygon::InitOriented),
            // we still end up preserving the original cyclic vertex order.
            var pos = 0
            var sawGap = false
            for (i in loop.indices) {
                val cmp = minInputIds[loop[i]] - minInputIds[loop[pos]]
                if (cmp < 0) {
                    sawGap = true
                } else if (cmp > 0 || !sawGap) {
                    pos = i
                    sawGap = false
                }
            }
            if (++pos == loop.size) pos = 0  // Convert loop end to loop start.
            loop.rotate(0, pos, loop.size)
        }

        // Sorts the given edge chains (i.e., loops or polylines) by the minimum
        // input edge id of each chains's first edge.  This ensures that when the
        // output consists of multiple loops or polylines, they are sorted in the
        // same order as they were provided in the input.
        fun canonicalizeVectorOrder(minInputIds: List<InputEdgeId>, chains: MutableList<List<EdgeId>>) {
            chains.sortWith { a: List<EdgeId>, b: List<EdgeId> -> minInputIds[a[0]].compareTo(minInputIds[b[0]]) }
        }

        ////////////////////////////////////////////////////////////////////////
        //////////////// Helper Functions for Creating Graphs //////////////////

        /**
         * Given an unsorted collection of edges, transform them according to the given set of GraphOptions.
         * This includes actions such as discarding degenerate edges; merging duplicate edges; and canonicalizing
         * sibling edge pairs in several possible ways (e.g. discarding or creating them).
         * The output is suitable for passing to the Graph constructor.
         *
         * If options.edge_type) == EdgeType.UNDIRECTED, then all input edges should already have been transformed
         * into a pair of directed edges.
         *
         * Note that "options" may be modified by this method: in particular, the edge_type) can be changed if
         * sibling_pairs is CREATE or REQUIRE (see the description of GraphOptions)
         *
         * @param input_ids vector of the same length as "edges" that indicates which input edges were snapped to each
         * edge. This vector is also updated appropriately as edges are discarded, merged, etc.
         */
        fun processEdges(options: GraphOptions, edges: ArrayList<Edge>, input_ids: ArrayList<InputEdgeIdSetId>, id_set_lexicon: IdSetLexicon, error: S2Error) {
            val processor = EdgeProcessor(options, edges, input_ids, id_set_lexicon)
            processor.run(error)
            // Certain values of sibling_pairs() discard half of the edges and change
            // the edge_type() to DIRECTED (see the description of GraphOptions).
            if (options.siblingPairs == SiblingPairs.REQUIRE || options.siblingPairs == SiblingPairs.CREATE) {
                options.edgeType = EdgeType.DIRECTED
            }
        }

        // Given a set of vertices and edges, removes all vertices that do not have
        // any edges and returned the new, minimal set of vertices.  Also updates
        // each edge in "edges" to correspond to the new vertex numbering.  (Note
        // that this method does *not* merge duplicate vertices, it simply removes
        // vertices of degree zero.)
        //
        // The new vertex ordering is a subsequence of the original ordering,
        // therefore if the edges were lexicographically sorted before calling this
        // method then they will still be sorted after calling this method.
        //
        // The extra argument "tmp" points to temporary storage used by this method.
        // All calls to this method from a single thread can reuse the same
        // temporary storage.  It should initially point to an empty vector.  This
        // can make a big difference to efficiency when this method is called many
        // times (e.g. to extract the vertices for different layers), since the
        // incremental running time for each layer becomes O(edges.size()) rather
        // than O(vertices.size() + edges.size()).
        fun filterVertices(vertices: List<S2Point>, edges: List<Edge>, vmap: MutableList<VertexId>): List<S2Point> {
            // Gather the vertices that are actually used.
            val used = ArrayList<VertexId>(2 * edges.size)
            for (e: Edge in edges) {
                used.add(e.first)
                used.add(e.second)
            }
            // Sort the vertices and find the distinct ones.
            used.sortAndRemoveDuplicates()

            // Build the list of new vertices, and generate a map from old vertex id to
            // new vertex id.
            vmap.resize(vertices.size)
            val new_vertices = ArrayList<S2Point>(used.size)
            for (i in 0 until used.size) {
                new_vertices[i] = vertices[used[i]]
                vmap[used[i]] = i
            }
            // Update the edges.
            for (e: Edge in edges) {
                e.first = vmap[e.first]
                e.second = vmap[e.second]
            }
            return new_vertices
        }

        // A comparison function that allows stable sorting with std::sort (which is
        // fast but not stable).  It breaks ties between equal edges by comparing
        // their edge ids.
        fun stableLessThan(a: Edge, b: Edge, ai: EdgeId, bi: EdgeId): Boolean {
            // The following is simpler but the compiler (2016) doesn't optimize it as
            // well as it should:
            //   return make_pair(a, ai) < make_pair(b, bi)
            if (a.first < b.first) return true
            if (b.first < a.first) return false
            if (a.second < b.second) return true
            if (b.second < a.second) return false
            return ai < bi  // Stable sort.
        }

        // Given a set of duplicate outgoing edges (v0, v1) and a set of duplicate
        // incoming edges (v1, v0), this method assigns each edge an integer "rank" so
        // that the edges are sorted in a consistent order with respect to their
        // orderings around "v0" and "v1".  Usually there is just one edge, in which
        // case this is easy.  Sometimes there is one edge in each direction, in which
        // case the outgoing edge is always ordered before the incoming edge.
        //
        // In general, we allow any number of duplicate edges in each direction, in
        // which case outgoing edges are interleaved with incoming edges so as to
        // create as many degenerate (two-edge) loops as possible.  In order to get a
        // consistent ordering around "v0" and "v1", we move forwards through the list
        // of outgoing edges and backwards through the list of incoming edges.  If
        // there are more incoming edges, they go at the beginning of the ordering,
        // while if there are more outgoing edges then they go at the end.
        //
        // For example, suppose there are 2 edges "a,b" from "v0" to "v1", and 4 edges
        // "w,x,y,z" from "v1" to "v0".  Using lower/upper case letters to represent
        // incoming/outgoing edges, the clockwise ordering around v0 would be zyAxBw,
        // and the clockwise ordering around v1 would be WbXaYZ.  (Try making a
        // diagram with each edge as a separate arc.)
        private fun addVertexEdges(outBegin: EdgeId, outEnd: EdgeId, inBegin: EdgeId, inEnd: EdgeId, v1: VertexId, v0Edges: MutableList<VertexEdge>) {
            var rank = 0
            var inIdx = inEnd
            var outIdx = outBegin
            // Any extra incoming edges go at the beginning of the ordering.
            while (inIdx - inBegin > outEnd - outIdx) {
                v0Edges.add(VertexEdge(true, --inIdx, v1, rank++))
            }
            // Next we interleave as many outgoing and incoming edges as possible.
            while (inIdx > inBegin) {
                v0Edges.add(VertexEdge(false, outIdx++, v1, rank++))
                v0Edges.add(VertexEdge(true, --inIdx, v1, rank++))
            }
            // Any extra outgoing edges to at the end of the ordering.
            while (outEnd > outIdx) {
                v0Edges.add(VertexEdge(false, outIdx++, v1, rank++))
            }

            logger.trace {
                """AddVertexEdges: outBegin=$outBegin, outEnd=$outEnd, inBegin=$inBegin, inEnd=$inEnd, v1=$v1
                |v0Edges: $v0Edges 
            """.trimMargin()
            }
        }

        // Encodes the index of one of the two complements of each component
        // (a.k.a. the "slot", either 0 or 1) as a negative EdgeId.
        private fun markEdgeUsed(slot: Int): EdgeId = -1 - slot
    }
}
