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
import dilivia.collections.assign
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.InputEdgeId
import mu.KotlinLogging
import java.util.*

class PolylineBuilder(val g: Graph) {

    private val inMap: VertexInMap = VertexInMap(g)
    private val outMap: VertexOutMap = VertexOutMap(g)
    private val siblingMap: MutableList<EdgeId>

    private val minInputIds: List<InputEdgeId> = g.getMinInputEdgeIds()
    private val directed: Boolean = g.options.edgeType == EdgeType.DIRECTED
    private var edgesLeft: Int = g.numEdges / (if (directed) 1 else 2)
    private val used: MutableList<Boolean>

    // A map of (outdegree(v) - indegree(v)) considering used edges only.
    private val excessUsed: TreeMap<VertexId, Int> = TreeMap()

    init {

        this.used = ArrayList(g.numEdges)
        this.used.assign(g.numEdges, false)
        if (!directed) {
            this.siblingMap = inMap.inEdgeIds.toMutableList()
            g.makeSiblingMap(siblingMap)
        } else this.siblingMap = mutableListOf()

        logger.trace { "constructor | g.edges = ${g.edges}" }
        logger.trace { "constructor | out = $outMap" }
        logger.trace { "constructor | in = $inMap, inEdges = ${inMap.inEdgeIds}" }
    }

    fun buildPaths(): List<EdgePolyline> {
        // First build polylines starting at all the vertices that cannot be in the
        // polyline interior (i.e., indegree != 1 or outdegree != 1 for directed
        // edges, or degree != 2 for undirected edges).  We consider the possible
        // starting edges in input edge id order so that we preserve the input path
        // direction even when undirected edges are used.  (Undirected edges are
        // represented by sibling pairs where only the edge in the input direction
        // is labeled with an input edge id.)
        val polylines = mutableListOf<EdgePolyline>()
        val edges = g.getInputEdgeOrder(minInputIds)
        for (e: EdgeId in edges) {
            if (!used[e] && !is_interior(g.edge(e).first)) {
                polylines.add(buildPath(e))
            }
        }
        // If there are any edges left, they form non-intersecting loops.  We build
        // each loop and then canonicalize its edge order.  We consider candidate
        // starting edges in input edge id order in order to preserve the input
        // direction of undirected loops.  Even so, we still need to canonicalize
        // the edge order to ensure that when an input edge is split into an edge
        // chain, the loop does not start in the middle of such a chain.
        for (e: EdgeId in edges) {
            if (edgesLeft == 0) break
            if (used[e]) continue
            val polyline: EdgePolyline = buildPath(e)
            Graph.canonicalizeLoopOrder(minInputIds, polyline)
            polylines.add(polyline)
        }
        PreConditions.checkEQ(0, edgesLeft)

        // Sort the polylines to correspond to the input order (if possible).
        @Suppress("UNCHECKED_CAST")
        Graph.canonicalizeVectorOrder(minInputIds, polylines as MutableList<List<EdgeId>>)
        return polylines
    }

    fun buildWalks(): List<EdgePolyline> {
        logger.trace { "buildWalks |" }
        // Note that some of this code is worst-case quadratic in the maximum vertex
        // degree.  This could be fixed with a few extra arrays, but it should not
        // be a problem in practice.

        // First, build polylines from all vertices where outdegree > indegree (or
        // for undirected edges, vertices whose degree is odd).  We consider the
        // possible starting edges in input edge id order, for idempotency in the
        // case where multiple input polylines share vertices or edges.
        val polylines = mutableListOf<EdgePolyline>()
        val edges = g.getInputEdgeOrder(minInputIds)
        for (e in edges) {
            logger.trace { "buildWalks | Process edge $e" }
            if (used[e]) continue
            val v = g.edge(e).first
            var excess = excess_degree(v)
            if (excess <= 0) continue
            excess -= excessUsed.getOrDefault(v, 0)
            if ((directed && excess <= 0) || (!directed && (excess % 2 == 0))) continue
            excessUsed.compute(v) { _, value -> value?.let { v -> v + 1 } ?: 1 }
            polylines.add(buildWalk(v))
            excessUsed.compute(g.edge(polylines.last().last()).second) { _, value -> value?.let { v -> v - 1 } ?: -1 }
        }
        // Now all vertices have outdegree == indegree (or even degree if undirected
        // edges are being used).  Therefore all remaining edges can be assembled
        // into loops.  We first try to expand the existing polylines if possible by
        // adding loops to them.
        if (edgesLeft > 0) {
            for (polyline in polylines) {
                maximizeWalk(polyline)
            }
        }
        // Finally, if there are still unused edges then we build loops.  If the
        // input is a polyline that forms a loop, then for idempotency we need to
        // start from the edge with minimum input edge id.  If the minimal input
        // edge was split into several edges, then we start from the first edge of
        // the chain.
        var i = 0
        while (i < edges.size && edgesLeft > 0) {
            val e = edges[i]
            if (used[e]) {
                ++i
                continue
            }

            // Determine whether the origin of this edge is the start of an edge
            // chain.  To do this, we test whether (outdegree - indegree == 1) for the
            // origin, considering only unused edges with the same minimum input edge
            // id.  (Undirected edges have input edge ids in one direction only.)
            val v = g.edge(e).first
            val id = minInputIds[e]
            var excess = 0
            var j = i
            while (j < edges.size && minInputIds[edges[j]] == id) {
                val e2 = edges[j]
                if (used[e2]) {
                    ++j
                    continue
                }
                if (g.edge(e2).first == v) ++excess
                if (g.edge(e2).second == v) --excess
                ++j
            }
            // It is also acceptable to start a polyline from any degenerate edge.
            if (excess == 1 || g.edge(e).second == v) {
                val polyline = buildWalk(v)
                maximizeWalk(polyline)
                polylines.add(polyline)
            }
            ++i
        }
        checkEQ(0, edgesLeft)

        // Sort the polylines to correspond to the input order (if possible).
        @Suppress("UNCHECKED_CAST")
        Graph.canonicalizeVectorOrder(minInputIds, polylines as MutableList<List<EdgeId>>)
        return polylines;
    }

    private fun is_interior(v: VertexId): Boolean = if (directed) {
        inMap.degree(v) == 1 && outMap.degree(v) == 1
    } else {
        outMap.degree(v) == 2
    }

    private fun excess_degree(v: VertexId): Int =
        if (directed) outMap.degree(v) - inMap.degree(v) else outMap.degree(v) % 2

    private fun buildPath(e: EdgeId): EdgePolyline {
        // We simply follow edges until either we reach a vertex where there is a
        // choice about which way to go (where is_interior(v) is false), or we
        // return to the starting vertex (if the polyline is actually a loop).
        var currentEdgeId = e
        val polyline: EdgePolyline = mutableListOf()
        val start: Int = g.edge(currentEdgeId).first
        while (true) {
            polyline.add(currentEdgeId)
            assert(!used[currentEdgeId])
            used[currentEdgeId] = true
            if (!directed) used[siblingMap[currentEdgeId]] = true
            --edgesLeft
            val v: VertexId = g.edge(currentEdgeId).second
            if (!is_interior(v) || v == start) break
            if (directed) {
                PreConditions.checkEQ(1, outMap.degree(v))
                currentEdgeId = outMap.edgeIds(v).first
            } else {
                PreConditions.checkEQ(2, outMap.degree(v))
                for (e2 in outMap.edgeIds(v)) {
                    if (!used[e2]) currentEdgeId = e2
                }
            }
        }
        return polyline
    }

    private fun buildWalk(v: VertexId): EdgePolyline {
        logger.trace { "buildWalk | v = $v" }
        val polyline = mutableListOf<EdgeId>()
        var currentVertex = v
        while (true) {
            logger.trace { "buildWalk | used = $used" }
            // Follow the edge with the smallest input edge id.
            var bestEdge: EdgeId = -1
            var best_out_id: InputEdgeId = Int.MAX_VALUE
            logger.trace { "buildWalk | outEdges($currentVertex) = ${outMap.edgeIds(currentVertex)}" }
            for (e: EdgeId in outMap.edgeIds(currentVertex)) {
                logger.trace { "buildWalk | e = $e, bestEdge = $bestEdge, best_out_id = $best_out_id, minInputIds[e] = ${minInputIds[e]}" }
                if (used[e] || minInputIds[e] >= best_out_id) continue
                best_out_id = minInputIds[e]
                bestEdge = e
            }
            logger.trace { "buildWalk | bestEdge = $bestEdge" }
            if (bestEdge < 0) return polyline
            // For idempotency when there are multiple input polylines, we stop the
            // walk early if "best_edge" might be a continuation of a different
            // incoming edge.
            val excess = excess_degree(currentVertex) - excessUsed.getOrDefault(currentVertex, 0)
            logger.trace { "buildWalk | excess = $excess" }
            if ((directed && excess < 0) || (!directed && (excess % 2 == 1))) {
                logger.trace { "buildWalk | process in edges = ${inMap.edgeIds(currentVertex)}" }
                for (e: EdgeId in inMap.edgeIds(currentVertex)) {
                    if (!used[e] && minInputIds[e] <= best_out_id) {
                        return polyline
                    }
                }
            }
            polyline.add(bestEdge)
            used[bestEdge] = true
            if (!directed) used[siblingMap[bestEdge]] = true
            --edgesLeft
            currentVertex = g.edge(bestEdge).second
        }

    }

    private fun maximizeWalk(polyline: EdgePolyline) {
        // Examine all vertices of the polyline and check whether there are any
        // unused outgoing edges.  If so, then build a loop starting at that vertex
        // and insert it into the polyline.  (The walk is guaranteed to be a loop
        // because this method is only called when all vertices have equal numbers
        // of unused incoming and outgoing edges.)
        for (i in 0..polyline.size) {
            val v: VertexId = if (i == 0) g.edge(polyline[i]).first else g.edge(polyline[i - 1]).second
            for (e: EdgeId in outMap.edgeIds(v)) {
                if (!used[e]) {
                    val loop: EdgePolyline = buildWalk(v)
                    PreConditions.checkEQ(v, g.edge(loop.last()).second)
                    polyline.addAll(i, loop)
                    assert(used[e])  // All outgoing edges from "v" are now used.
                    break
                }
            }
        }
    }


    companion object {
        private val logger = KotlinLogging.logger(PolylineBuilder::class.java.name)
    }
}
