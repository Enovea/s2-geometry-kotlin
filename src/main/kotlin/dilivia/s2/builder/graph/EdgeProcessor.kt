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

import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkState
import dilivia.s2.S2Error
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.InputEdgeIdSetId
import org.apache.commons.math3.util.FastMath.max


class EdgeProcessor(
    val options: GraphOptions,
    val edges: ArrayList<Edge>,
    val inputIds: ArrayList<InputEdgeIdSetId>,
    val idSetLexicon: IdSetLexicon
) {

    private val outEdges = mutableListOf<EdgeId>()
    private val inEdges = mutableListOf<EdgeId>()

    private val new_edges_ = ArrayList<Edge>()
    private val new_input_ids_ = ArrayList<InputEdgeIdSetId>()

    private val tmp_ids_ = mutableListOf<InputEdgeId>()

    init {
        // Sort the outgoing and incoming edges in lexigraphic order.  We use a
        // stable sort to ensure that each undirected edge becomes a sibling pair,
        // even if there are multiple identical input edges.
        repeat(edges.size) { i ->
            outEdges.add(i)
            inEdges.add(i)
        }
        outEdges.sortWith { a, b -> if (Graph.stableLessThan(edges[a], edges[b], a, b)) -1 else 1 }
        inEdges.sortWith { a, b -> if (Graph.stableLessThan(edges[a].reversed(), edges[b].reversed(), a, b)) -1 else 1 }

        new_edges_.ensureCapacity(edges.size)
        new_input_ids_.ensureCapacity(edges.size)
    }

    fun run(error: S2Error) {
        val numEdges = edges.size
        if (numEdges == 0) return

        // Walk through the two sorted arrays performing a merge join.  For each
        // edge, gather all the duplicate copies of the edge in both directions
        // (outgoing and incoming).  Then decide what to do based on "options_" and
        // how many copies of the edge there are in each direction.
        var outIdx = 0
        var inIdx = 0
        var outEdge = edges[outEdges[outIdx]]
        var inEdge = edges[inEdges[inIdx]]
        val sentinel = Edge(VertexId.MAX_VALUE, VertexId.MAX_VALUE)
        loop@ while (true) {
            val edge = minOf(outEdge, inEdge.reversed())
            if (edge == sentinel) break

            val outBegin = outIdx
            val inBegin = inIdx
            while (outEdge == edge) {
                outEdge = if (++outIdx == numEdges) sentinel else edges[outEdges[outIdx]]
            }
            while (inEdge.reversed() == edge) {
                inEdge = if (++inIdx == numEdges) sentinel else edges[inEdges[inIdx]]
            }
            val nOut = outIdx - outBegin
            val nIn = inIdx - inBegin

            when {
                edge.first == edge.second -> {
                    checkEQ(nOut, nIn)
                    if (options.degenerateEdges == DegenerateEdges.DISCARD) {
                        continue@loop
                    }
                    if (options.degenerateEdges == DegenerateEdges.DISCARD_EXCESS &&
                            ((outBegin > 0 && edges[outEdges[outBegin - 1]].first == edge.first) ||
                                    (outIdx < numEdges && edges[outEdges[outIdx]].first == edge.first) ||
                                    (inBegin > 0 && edges[inEdges[inBegin - 1]].second == edge.first) ||
                                    (inIdx < numEdges && edges[inEdges[inIdx]].second == edge.first))) {
                        continue@loop  // There were non-degenerate incident edges, so discard.
                    }
                    if (options.edgeType == EdgeType.UNDIRECTED &&
                            (options.siblingPairs == SiblingPairs.REQUIRE || options.siblingPairs == SiblingPairs.CREATE)) {
                        // When we have undirected edges and are guaranteed to have siblings,
                        // we cut the number of edges in half (see s2builder.h).
                        checkEQ(0, nOut and 1)  // Number of edges is always even.
                        addEdges(if (options.duplicateEdges == DuplicateEdges.MERGE) 1 else (nOut / 2), edge, mergeInputIds(outBegin, outIdx))
                    } else if (options.duplicateEdges == DuplicateEdges.MERGE) {
                        addEdges(if (options.edgeType == EdgeType.UNDIRECTED) 2 else 1, edge, mergeInputIds(outBegin, outIdx))
                    } else if (options.siblingPairs == SiblingPairs.DISCARD || options.siblingPairs == SiblingPairs.DISCARD_EXCESS) {
                        // Any SiblingPair option that discards edges causes the labels of all
                        // duplicate edges to be merged together (see s2builder.h).
                        addEdges(nOut, edge, mergeInputIds(outBegin, outIdx))
                    } else {
                        copyEdges(outBegin, outIdx)
                    }
                }
                options.siblingPairs == SiblingPairs.KEEP -> {
                    if (nOut > 1 && options.duplicateEdges == DuplicateEdges.MERGE) {
                        addEdge(edge, mergeInputIds(outBegin, outIdx))
                    } else {
                        copyEdges(outBegin, outIdx)
                    }
                }
                options.siblingPairs == SiblingPairs.DISCARD -> {
                    if (options.edgeType == EdgeType.DIRECTED) {
                        // If n_out == n_in: balanced sibling pairs
                        // If n_out < n_in:  unbalanced siblings, in the form AB, BA, BA
                        // If n_out > n_in:  unbalanced siblings, in the form AB, AB, BA
                        if (nOut <= nIn) continue@loop
                        // Any option that discards edges causes the labels of all duplicate
                        // edges to be merged together (see s2builder.h).
                        addEdges(if (options.duplicateEdges == DuplicateEdges.MERGE) 1 else (nOut - nIn), edge, mergeInputIds(outBegin, outIdx))
                    } else {
                        if ((nOut and 1) == 0) continue@loop
                        addEdge(edge, mergeInputIds(outBegin, outIdx))
                    }
                }
                options.siblingPairs == SiblingPairs.DISCARD_EXCESS -> {
                    if (options.edgeType == EdgeType.DIRECTED) {
                        // See comments above.  The only difference is that if there are
                        // balanced sibling pairs, we want to keep one such pair.
                        if (nOut < nIn) continue@loop
                        addEdges(if (options.duplicateEdges == DuplicateEdges.MERGE) 1 else max(1, nOut - nIn), edge, mergeInputIds(outBegin, outIdx))
                    } else {
                        addEdges(if ((nOut and 1) != 0) 1 else 2, edge, mergeInputIds(outBegin, outIdx))
                    }
                }
                else -> {
                    checkState { options.siblingPairs == SiblingPairs.REQUIRE || options.siblingPairs == SiblingPairs.CREATE }
                    if (error.isOk() && options.siblingPairs == SiblingPairs.REQUIRE &&
                            if (options.edgeType == EdgeType.DIRECTED) (nOut != nIn) else ((nOut and 1) != 0)
                    ) {
                        error.init(S2Error.BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
                                "Expected all input edges to have siblings,  but some were missing")
                    }
                    if (options.duplicateEdges == DuplicateEdges.MERGE) {
                        addEdge(edge, mergeInputIds(outBegin, outIdx))
                    } else if (options.edgeType == EdgeType.UNDIRECTED) {
                        // Convert graph to use directed edges instead (see documentation of
                        // REQUIRE/CREATE for undirected edges).
                        addEdges((nOut + 1) / 2, edge, mergeInputIds(outBegin, outIdx))
                    } else {
                        copyEdges(outBegin, outIdx)
                        if (nIn > nOut) {
                            // Automatically created edges have no input edge ids or labels.
                            addEdges(nIn - nOut, edge, IdSetLexicon.emptySetId())
                        }
                    }
                }
            }
        }

        edges.clear()
        edges.addAll(new_edges_)
        edges.trimToSize()
        inputIds.clear()
        inputIds.addAll(new_input_ids_)
        inputIds.trimToSize()
    }


    private fun addEdge(edge: Edge, inputEdgeIdSetId: InputEdgeIdSetId) {
        new_edges_.add(edge)
        new_input_ids_.add(inputEdgeIdSetId)
    }

    private fun addEdges(numEdges: Int, edge: Edge, inputEdgeIdSetId: InputEdgeIdSetId) {
        for (i in 0 until numEdges) {
            addEdge(edge, inputEdgeIdSetId)
        }
    }

    private fun copyEdges(outBegin: Int, outEnd: Int) {
        for (i in outBegin until outEnd) {
            addEdge(edges[outEdges[i]], inputIds[outEdges[i]])
        }
    }

    private fun mergeInputIds(outBegin: Int, outEnd: Int): InputEdgeIdSetId {
        if (outEnd - outBegin == 1) {
            return inputIds[outEdges[outBegin]]
        }
        tmp_ids_.clear()
        for (i in outBegin until outEnd) {
            for (id in idSetLexicon.idSet(inputIds[outEdges[i]])) {
                tmp_ids_.add(id)
            }
        }
        return idSetLexicon.add(tmp_ids_)
    }

}
