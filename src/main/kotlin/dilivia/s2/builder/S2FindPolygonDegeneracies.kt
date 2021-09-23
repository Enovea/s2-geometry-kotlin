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
package dilivia.s2.builder

import dilivia.ComparisonChain
import dilivia.PreConditions.checkNE
import dilivia.collections.resize
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Predicates
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.graph.VertexId
import dilivia.s2.builder.graph.VertexInMap
import dilivia.s2.builder.graph.VertexOutMap
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2CrossingEdgeQuery
import dilivia.s2.shape.S2ContainsVertexQuery
import dilivia.s2.shape.S2GraphShape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging

object S2FindPolygonDegeneracies {
    // A polygon degeneracy is either a degenerate edge (an edge from a vertex to
// itself) or a sibling edge pair (consisting of an edge and its corresponding
// reverse edge).  "is_hole" indicates whether the degeneracy corresponds to a
// polygon hole (as opposed to a polygon shell).
//
// Degeneracies are not allowed to coincide with any non-degenerate portions
// of the polygon's boundary (since that would make it impossible to classify
// the degeneracy as a shell or hole).  Specifically, degenerate edges must
// coincide only with other degenerate edges, and sibling pairs must coincide
// only with other sibling pairs.  (Below we require a slightly stronger
// condition, namely that sibling pairs cannot coincide with any other edges.)
    data class PolygonDegeneracy(val edgeId: Int = 0, val isHole: Boolean = false) : Comparable<PolygonDegeneracy> {

        override fun compareTo(other: PolygonDegeneracy): Int {
            return ComparisonChain.start().compare(edgeId, other.edgeId).compareFalseFirst(isHole, other.isHole).result()
        }

    }

    // Given a graph representing a polygon, finds all degenerate edges and
// sibling pairs and classifies them as being either shells or holes.  The
// result vector is sorted by edge id.
//
// REQUIRES: g.options().edge_type() == DIRECTED
// REQUIRES: g.options().sibling_pairs() == DISCARD_EXCESS (or DISCARD)
// REQUIRES: g.options().degenerate_edges() == DISCARD_EXCESS (or DISCARD)
//
// Usually callers will want to specify SiblingPairs::DISCARD_EXCESS and
// DegenerateEdges::DISCARD_EXCESS in order to remove all redundant
// degeneracies.  DISCARD is also allowed in case you want to keep only one
// type of degeneracy (i.e., degenerate edges or sibling pairs).
//
// If the graph edges cannot be assembled into loops, the result is undefined.
// (An error may or may not be returned.)
    fun findPolygonDegeneracies(g: Graph, error: S2Error): List<PolygonDegeneracy> {
        checkGraphOptions(g)
        if (g.options.degenerateEdges == DegenerateEdges.DISCARD &&
                g.options.siblingPairs == SiblingPairs.DISCARD) {
            return emptyList()  // All degeneracies have already been discarded.
        }
        return DegeneracyFinder(g).run(error)
    }

    // Given a graph representing a polygon, returns true the graph consists
// entirely of degenerate edges and/or sibling pairs.  Such a graph represents
// either the empty polygon together with a collection of degenerate shells,
// or the full polygon together with a collection of degenerate holes.
//
// REQUIRES: g.options().edge_type() == DIRECTED
// REQUIRES: g.options().sibling_pairs() == DISCARD_EXCESS (or DISCARD)
// REQUIRES: g.options().degenerate_edges() == DISCARD_EXCESS (or DISCARD)
    fun isFullyDegenerate(g: Graph): Boolean {
        checkGraphOptions(g)
        val edges = g.edges
        for (e in 0 until g.numEdges) {
            val edge = edges[e]
            if (edge.first == edge.second) continue
            if (edges.binarySearch(Graph.reverse(edge)) < 0) {
                return false
            }
        }
        return true
    }


    private fun checkGraphOptions(g: Graph) {
        check(g.options.edgeType == EdgeType.DIRECTED)
        check(g.options.degenerateEdges == DegenerateEdges.DISCARD ||
                g.options.degenerateEdges == DegenerateEdges.DISCARD_EXCESS)
        check(g.options.siblingPairs == SiblingPairs.DISCARD ||
                g.options.siblingPairs == SiblingPairs.DISCARD_EXCESS)
    }

    // The algorithm builds a set of connected components containing all edges
// that form degeneracies.  The shell/hole status of each degeneracy is
// initially unknown, and is expressed relative to the root vertex: "is_hole"
// means that the degeneracy is a hole if and only if the root vertex turns
// out to be inside the polygon.
    private data class Component(
            // The root vertex from which this component was built.
        var root: VertexId = -1,

            // +1 if "root" inside the polygon, -1 if outside, and 0 if unknown.
        var rootSign: Int = 0,

            // The degeneracies found in this component.  "is_hole" is expressed
            // relative to the root vertex: the degeneracy is a hole iff the root vertex
            // turns out to be inside the polygon (i.e., root_sign > 0).
        val degeneracies: MutableList<PolygonDegeneracy> = mutableListOf()
    )

    // The actual implementation of FindPolygonDegeneracies.
    private class DegeneracyFinder(val g: Graph) {

        private val inMap: VertexInMap = VertexInMap(g)
        private val outMap: VertexOutMap = VertexOutMap(g)
        private val is_vertex_used = mutableListOf<Boolean>()        // Has vertex been visited?
        private val is_edge_degeneracy = mutableListOf<Boolean>()    // Belongs to a degeneracy?
        private val is_vertex_unbalanced = mutableListOf<Boolean>()  // Has unbalanced sibling pairs?

        fun run(error: S2Error): List<PolygonDegeneracy> {
            // Mark all degenerate edges and sibling pairs in the "is_edge_degeneracy_"
            // vector, and mark any vertices with unbalanced edges in the
            // "is_vertex_unbalanced_" vector.
            val numDegeneracies = computeDegeneracies()
            if (numDegeneracies == 0) return emptyList()

            // If all edges are degenerate, then use IsFullPolygon() to classify the
            // degeneracies (they are necessarily all the same type).
            if (numDegeneracies == g.numEdges) {
                val isHole = g.isFullPolygon(error)
                return (0 until g.numEdges).map { edgeId -> PolygonDegeneracy(edgeId, isHole) }
            }

            // Otherwise repeatedly build components starting from an unvisited
            // degeneracy.  (This avoids building components that don't contain any
            // degeneracies.)  Each component records the "is_hole" status of each
            // degeneracy relative to the root vertex of that component.  If the
            // component contains any non-degenerate portions, then we also determine
            // whether the root vertex is contained by the component (root_sign).
            // In addition we keep track of the number of components that were
            // completely degenerate (to help us decide whether to build an index).
            val components = mutableListOf<Component>()
            var knownVertex: VertexId = -1
            var knownVertexSign = 0
            var numUnknownSigns = 0;
            is_vertex_used.resize(g.numVertices)
            for (e in 0 until g.numEdges) {
                if (is_edge_degeneracy[e]) {
                    val root: VertexId = g.edge(e).first
                    if (is_vertex_used[root]) continue
                    val component: Component = buildComponent(root)
                    if (component.rootSign == 0) {
                        ++numUnknownSigns
                    } else {
                        knownVertex = root
                        knownVertexSign = component.rootSign
                    }
                    components.add(component)
                }
            }

            // If some components have an unknown root_sign (i.e., it is unknown whether
            // the root vertex is contained by the polygon or not), we determine the
            // sign of those root vertices by counting crossings starting from a vertex
            // whose sign is known.  Depending on how many components we need to do this
            // for, it may be worthwhile to build an index first.
            if (numUnknownSigns > 0) {
                if (knownVertexSign == 0) {
                    knownVertex = findUnbalancedVertex()
                    knownVertexSign = containsVertexSign(knownVertex);
                }
                val kMaxUnindexedContainsCalls = 20;  // Tuned using benchmarks.
                if (numUnknownSigns <= kMaxUnindexedContainsCalls) {
                    computeUnknownSignsBruteForce(knownVertex, knownVertexSign, components)
                } else {
                    computeUnknownSignsIndexed(knownVertex, knownVertexSign, components)
                }
            }
            // Finally we convert the "is_hole" status of each degeneracy from a
            // relative value (compared to the component's root vertex) to an absolute
            // one, and sort all the degeneracies by EdgeId.
            return mergeDegeneracies(components)
        }

        //private:

        private fun computeDegeneracies(): Int {
            is_edge_degeneracy.resize(g.numEdges)
            is_vertex_unbalanced.resize(g.numVertices)
            var numDegeneracies = 0
            val inEdgeIds = inMap.inEdgeIds
            val n = g.numEdges
            var inId = 0
            for (outId in 0 until n) {
                val out_edge = g.edge(outId);
                if (out_edge.first == out_edge.second) {
                    is_edge_degeneracy[outId] = true
                    ++numDegeneracies
                } else {
                    while (inId < n && Graph.reverse(g.edge(inEdgeIds[inId])) < out_edge) {
                        ++inId
                    }
                    if (inId < n && Graph.reverse(g.edge(inEdgeIds[inId])) == out_edge) {
                        is_edge_degeneracy[outId] = true
                        ++numDegeneracies;
                    } else {
                        // This edge does not have a sibling, which mean that we can determine
                        // whether either vertex is contained by the polygon (using semi-open
                        // boundaries) by examining only the edges incident to that vertex.
                        // We only mark the first vertex since there is no advantage to
                        // finding more than one unbalanced vertex per connected component.
                        is_vertex_unbalanced[out_edge.first] = true
                    }
                }
            }
            return numDegeneracies
        }


        // Build a connected component starting at the given root vertex.  The
        // information returned includes: the root vertex, whether the containment
        // status of the root vertex could be determined using only the edges in this
        // component, and a vector of the edges that belong to degeneracies along with
        // the shell/hole status of each such edge relative to the root vertex.
        private fun buildComponent(root: VertexId): Component {
            val result: Component = Component()
            result.root = root;
            // We keep track of the frontier of unexplored vertices, and whether each
            // vertex is on the same side of the polygon boundary as the root vertex.
            val frontier = mutableListOf<Pair<VertexId, Boolean>>()
            frontier.add(Pair(root, true))
            is_vertex_used[root] = true
            while (frontier.isNotEmpty()) {
                val v0: VertexId = frontier.last().first
                val v0_same_inside = frontier.last().second;  // Same as root vertex?
                frontier.removeLast();
                if (result.rootSign == 0 && is_vertex_unbalanced[v0]) {
                    val v0_sign = containsVertexSign(v0);
                    checkNE(v0_sign, 0)
                    result.rootSign = if (v0_same_inside) v0_sign else -v0_sign
                }
                for (e: EdgeId in outMap.edgeIds(v0)) {
                    val v1: VertexId = g.edge(e).second
                    var same_inside = v0_same_inside xor crossingParity(v0, v1, false)
                    if (is_edge_degeneracy[e]) {
                        result.degeneracies.add(PolygonDegeneracy(e, same_inside))
                    }
                    if (is_vertex_used[v1]) continue
                    same_inside = same_inside xor crossingParity(v1, v0, true)
                    frontier.add(Pair(v1, same_inside))
                    is_vertex_used[v1] = true
                }
            }
            return result
        }


        // Counts the number of times that (v0, v1) crosses the edges incident to v0,
        // and returns the result modulo 2.  This is equivalent to calling
        // S2::VertexCrossing for the edges incident to v0, except that this
        // implementation is more efficient (since it doesn't need to determine which
        // two edge vertices are the same).
        //
        // If "include_same" is false, then the edge (v0, v1) and its sibling (v1, v0)
        // (if any) are excluded from the parity calculation.
        private fun crossingParity(v0: VertexId, v1: VertexId, includeSame: Boolean): Boolean {
            var crossings = 0
            val p0: S2Point = g.vertex(v0)
            val p1: S2Point = g.vertex(v1)
            val p0Ref: S2Point = S2PointUtil.ortho(p0)
            for (edge: Edge in outMap.edges(v0)) {
                if (edge.second == v1) {
                    if (includeSame) ++crossings;
                } else if (S2Predicates.orderedCCW(p0Ref, g.vertex(edge.second), p1, p0)) {
                    ++crossings
                }
            }
            for (e: EdgeId in inMap.edgeIds(v0)) {
                val edge: Edge = g.edge(e)
                if (edge.first == v1) {
                    if (includeSame) ++crossings;
                } else if (S2Predicates.orderedCCW(p0Ref, g.vertex(edge.first), p1, p0)) {
                    ++crossings
                }
            }
            return (crossings and 1) != 0
        }

        private fun findUnbalancedVertex(): VertexId {
            for (v: VertexId in 0 until g.numVertices) {
                if (is_vertex_unbalanced[v]) return v
            }
            logger.error { "Could not find previously marked unbalanced vertex" }
            return -1
        }

        private fun containsVertexSign(v0: VertexId): Int {
            val query = S2ContainsVertexQuery(g.vertex(v0))
            for (edge: Edge in outMap.edges(v0)) {
                query.addEdge(g.vertex(edge.second), 1)
            }
            for (e: EdgeId in inMap.edgeIds(v0)) {
                query.addEdge(g.vertex(g.edge(e).first), -1)
            }
            return query.containsSign()
        }

        // Determines any unknown signs of component root vertices by counting
        // crossings starting from a vertex whose sign is known.  This version simply
        // tests all edges for crossings.
        private fun computeUnknownSignsBruteForce(knownVertex: VertexId, knownVertexSign: Int, components: List<Component>) {
            val crosser = S2EdgeCrosser()
            for (component: Component in components) {
                if (component.rootSign != 0) continue
                var inside = knownVertexSign > 0
                crosser.init(g.vertex(knownVertex), g.vertex(component.root))
                for (e: EdgeId in 0 until g.numEdges) {
                    if (is_edge_degeneracy[e]) continue
                    val edge: Edge = g.edge(e)
                    inside = inside xor crosser.edgeOrVertexCrossing(g.vertex(edge.first), g.vertex(edge.second))
                }
                component.rootSign = if (inside) 1 else -1
            }
        }


        // Like ComputeUnknownSignsBruteForce, except that this method uses an index
        // to find the set of edges that cross a given edge.
        private fun computeUnknownSignsIndexed(known_vertex: VertexId, known_vertex_sign: Int, components: List<Component>) {
            val index = MutableS2ShapeIndex()
            index.add(S2GraphShape(g))
            val query = S2CrossingEdgeQuery(index)
            val crossing_edges = mutableListOf<ShapeEdgeId>()
            val crosser = S2EdgeCrosser()
            for (component: Component in components) {
                if (component.rootSign != 0) continue
                var inside = known_vertex_sign > 0
                crosser.init(g.vertex(known_vertex), g.vertex(component.root))
                query.getCandidates(g.vertex(known_vertex), g.vertex(component.root), index.shape(0)!!, crossing_edges)
                for (id: ShapeEdgeId in crossing_edges) {
                    var e = id.edgeId
                    if (is_edge_degeneracy[e]) continue
                    inside = inside xor crosser.edgeOrVertexCrossing(g.vertex(g.edge(e).first), g.vertex(g.edge(e).second))
                }
                component.rootSign = if (inside) 1 else -1
            }
        }

        // Merges the degeneracies from all components together, and computes the
        // final "is_hole" status of each edge (since up to this point, the "is_hole"
        // value has been expressed relative to the root vertex of each component).
        private fun mergeDegeneracies(components: List<Component>): List<PolygonDegeneracy> {
            val result = mutableListOf<PolygonDegeneracy>()
            for (component: Component in components) {
                checkNE(component.rootSign, 0)
                val invert = component.rootSign < 0
                for (d in component.degeneracies) {
                    result.add(PolygonDegeneracy(d.edgeId, d.isHole xor invert))
                }
            }
            result.sort()
            return result
        }

        companion object {
            private val logger = KotlinLogging.logger(DegeneracyFinder::class.java.name)
        }
    }

}
