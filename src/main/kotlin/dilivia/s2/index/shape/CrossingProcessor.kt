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

import dilivia.PreConditions
import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.requireEQ
import dilivia.s2.S2Point
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.S2Builder
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import java.util.*

// CrossingProcessor is a helper class that processes all the edges from one
// region that cross a specific edge of the other region.  It outputs the
// appropriate edges to an S2Builder, and outputs other information required
// by GraphEdgeClipper to the given vectors.
// The output of the CrossingProcessor consists of a subset of the input
// edges that are emitted to "builder_", and some auxiliary information
// that allows GraphEdgeClipper to determine which segments of those input
// edges belong to the output.  The auxiliary information consists of the
// dimension of each input edge, and set of input edges from the other
// region that cross each input input edge.
// Prepares to build output for the given polygon and polyline boundary
// models.  Edges are emitted to "builder", while other auxiliary data is
// appended to the given vectors.
//
// If a predicate is being evaluated (i.e., we do not need to construct the
// actual result), then "builder" and the various output vectors should all
// be nullptr.
class CrossingProcessor(
        val polygonModel: S2BooleanOperation.PolygonModel,
        val polylineModel: S2BooleanOperation.PolylineModel,
        val polylineLoopsHaveBoundaries: Boolean,
        val builder: S2Builder?,
        val inputDimensions: MutableList<Byte>,
        val inputCrossings: InputEdgeCrossings
) {
    // Fields set by StartBoundary:

    private var aRegionId: Int = -1
    private var bRegionId: Int = -1
    private var invertA: Boolean = false
    private var invertB: Boolean = false
    private var invertResult: Boolean = false
    private var isUnion: Boolean = false  // True if this is a UNION operation.

    // Fields set by StartShape:

    private lateinit var aShape: S2Shape
    private var aDimension: Int = -1

    // Fields set by StartChain:

    private var chainId: Int = -1
    private var chainStart: Int = -1
    private var chainLimit: Int = -1

    // Fields updated by ProcessEdge:

    // A temporary representation of input_crossings_ that is used internally
    // until all necessary edges from *both* polygons have been emitted to the
    // S2Builder.  This field is then converted by DoneBoundaryPair() into
    // the InputEdgeCrossings format expected by GraphEdgeClipper.
    //
    // The reason that we can't construct input_crossings_ directly is that it
    // uses InputEdgeIds to identify the edges from both polygons, and when we
    // are processing edges from the first polygon, InputEdgeIds have not yet
    // been assigned to the second polygon.  So instead this field identifies
    // edges from the first polygon using an InputEdgeId, and edges from the
    // second polygon using a (region_id, shape_id, edge_id) tuple (i.e., a
    // SourceId).
    //
    // All crossings are represented twice, once to indicate that an edge from
    // polygon 0 is crossed by an edge from polygon 1, and once to indicate that
    // an edge from polygon 1 is crossed by an edge from polygon 0.
    private val sourceEdgeCrossings: SourceEdgeCrossings = mutableListOf()

    // A map that translates from SourceId (the (region_id, shape_id,
    // edge_id) triple that identifies an S2ShapeIndex edge) to InputEdgeId (the
    // sequentially increasing numbers assigned to input edges by S2Builder).
    private val sourceIdMap: SourceIdMap = TreeMap()

    // Indicates whether the point being processed along the current edge chain
    // is in the polygonal interior of the opposite region, using semi-open
    // boundaries.  If "invert_b_" is true then this field is inverted.
    //
    // Equal to: b_index_.Contains(current point) ^ invert_b_
    private var inside: Boolean = false

    // The value of that "inside_" would have just before the end of the
    // previous edge added to S2Builder.  This value is used to determine
    // whether the GraphEdgeClipper state needs to be updated when jumping from
    // one edge chain to another.
    private var prevInside: Boolean = false

    // The maximum edge id of any edge in the current chain whose v0 vertex has
    // already been emitted.  This is used to determine when an isolated vertex
    // needs to be emitted, e.g. when two closed polygons share only a vertex.
    private var v0EmittedMaxEdgeId: Int = -1

    // True if the first vertex of the current chain has been emitted.  This is
    // used when processing loops in order to determine whether the first/last
    // vertex of the loop should be emitted as an isolated vertex.
    private var chainV0Emitted: Boolean = false

    // Starts processing edges from the given region.  "invert_a", "invert_b",
    // and "invert_result" indicate whether region A, region B, and/or the
    // result should be inverted, which allows operations such as union and
    // difference to be implemented.  For example, union is ~(~A & ~B).
    //
    // This method should be called in pairs, once to process the edges from
    // region A and once to process the edges from region B.
    fun startBoundary(aRegionId: Int, invertA: Boolean, invertB: Boolean, invertResult: Boolean) {
        logger.trace { "startBoundary(aRegionId: $aRegionId, invertA: $invertA, invertB: $invertB, invertResult: $invertResult)" }

        this.aRegionId = aRegionId
        this.bRegionId = 1 - aRegionId
        this.invertA = invertA
        this.invertB = invertB
        this.invertResult = invertResult
        isUnion = invertB && invertResult

        // Specify to GraphEdgeClipper how these edges should be clipped.
        setClippingState(S2BooleanOperationImpl.kSetReverseA, invertA != invertResult)
        setClippingState(S2BooleanOperationImpl.kSetInvertB, invertB);
    }

    // Starts processing edges from the given shape.
    fun startShape(aShape: S2Shape) {
        logger.trace { "startShape(aShape: ${aShape.id})" }
        this.aShape = aShape
        this.aDimension = aShape.dimension
    }

    // Starts processing edges from the given chain.
    fun startChain(chainId: Int, chain: S2Shape.Chain, inside: Boolean) {
        logger.trace { "startChain(chainId: $chainId, chain: $chain, inside: $inside)" }
        this.chainId = chainId
        this.chainStart = chain.start
        this.chainLimit = chain.start + chain.length
        this.inside = inside
        this.v0EmittedMaxEdgeId = chain.start - 1  // No edges emitted yet.
        this.chainV0Emitted = false
    }

    // Processes the given edge "a_id".  "it" should be positioned to the set of
    // edges from the other region that cross "a_id" (if any).
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    fun processEdge(aId: ShapeEdgeId, iter: CrossingIterator): Boolean {
        logger.trace { "processEdge(aId: $aId, iter)" }
        // chain_edge() is faster than edge() when there are multiple chains.
        val a = aShape.chainEdge(chainId, aId.edgeId - chainStart)
        val result = when (aDimension) {
            0 -> processEdge0(aId, a, iter)
            1 -> processEdge1(aId, a, iter)
            else -> {
                checkEQ(2, aDimension)
                processEdge2(aId, a, iter)
            }
        }
        logger.trace { "processEdge: $result" }
        return result
    }

    // This method should be called after each pair of calls to StartBoundary.
    // (The only operation that processes more than one pair of boundaries is
    // SYMMETRIC_DIFFERENCE, which computes the union of A-B and B-A.)
    //
    // Resets the state of the CrossingProcessor.
    // Translates the temporary representation of crossing edges (SourceId) into
    // the format expected by EdgeClippingLayer (InputEdgeId).
    fun doneBoundaryPair() {
        logger.trace { "doneBoundaryPair()" }
        // Add entries that translate the "special" crossings.
        sourceIdMap[S2BooleanOperation.SourceId(edgeId = S2BooleanOperationImpl.kSetInside)] = S2BooleanOperationImpl.kSetInside
        sourceIdMap[S2BooleanOperation.SourceId(edgeId = S2BooleanOperationImpl.kSetInvertB)] = S2BooleanOperationImpl.kSetInvertB
        sourceIdMap[S2BooleanOperation.SourceId(edgeId = S2BooleanOperationImpl.kSetReverseA)] = S2BooleanOperationImpl.kSetReverseA
        inputCrossings.ensureCapacity(inputCrossings.size + sourceEdgeCrossings.size)

        logger.trace { "sourceEdgeCrossings: $sourceEdgeCrossings" }
        logger.trace { "sourceIdMap: $sourceIdMap" }

        for (tmp in sourceEdgeCrossings) {
            val it = sourceIdMap[tmp.second.first]
            check(it != null)
            inputCrossings.add(Pair(tmp.first, S2BooleanOperationImpl.CrossingInputEdge(it, tmp.second.second)))
        }
        sourceEdgeCrossings.clear()
        sourceIdMap.clear()
    }

    // Indicates whether the point being processed along the current edge chain
    // is in the polygonal interior of the opposite region, using semi-open
    // boundaries.  If "invert_b_" is true then this field is inverted.
    //
    // This value along with the set of incident edges can be used to compute
    // whether the opposite region contains this point under any of the
    // supported boundary models (PolylineModel::CLOSED, etc).
    fun inside(): Boolean = inside

    private fun inputEdgeId(): InputEdgeId = inputDimensions.size

    // Returns true if the edges on either side of the first vertex of the
    // current edge have not been emitted.
    //
    // REQUIRES: This method is called just after updating "inside_" for "v0".
    private fun isV0Isolated(aId: ShapeEdgeId): Boolean = !inside && v0EmittedMaxEdgeId < aId.edgeId

    // Returns true if "a_id" is the last edge of the current chain, and the
    // edges on either side of the last vertex have not been emitted (including
    // the possibility that the chain forms a loop).
    private fun isChainLastVertexIsolated(aId: ShapeEdgeId): Boolean =
            (aId.edgeId == chainLimit - 1 && !chainV0Emitted && v0EmittedMaxEdgeId <= aId.edgeId)

    // Returns true if the given polyline edge contains "v0", taking into
    // account the specified PolylineModel.
    private fun polylineContainsV0(edgeId: Int, chainStart: Int): Boolean =
            (polylineModel != S2BooleanOperation.PolylineModel.OPEN || edgeId > chainStart)

    private fun addCrossing(crossing: SourceEdgeCrossing) {
        sourceEdgeCrossings.add(Pair(inputEdgeId(), crossing))
    }

    private fun setClippingState(parameter: InputEdgeId, state: Boolean) {
        addCrossing(SourceEdgeCrossing(S2BooleanOperation.SourceId(edgeId = parameter), state))
    }

    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun addEdge(aId: ShapeEdgeId, a: Edge, dimension: Int, interiorCrossings: Int): Boolean {
        if (builder == null) return false  // Boolean output.
        if (interiorCrossings > 0) {
            // Build a map that translates temporary edge ids (SourceId) to
            // the representation used by EdgeClippingLayer (InputEdgeId).
            val src_id = S2BooleanOperation.SourceId(aRegionId, aId.shapeId, aId.edgeId)
            sourceIdMap[src_id] = inputEdgeId()
        }
        // Set the GraphEdgeClipper's "inside" state to match ours.
        if (inside != prevInside) setClippingState(S2BooleanOperationImpl.kSetInside, inside)
        inputDimensions.add(dimension.toByte())
        builder.addEdge(a.v0, a.v1)
        inside = inside xor ((interiorCrossings and 1) != 0)
        prevInside = inside
        return true
    }

    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun addPointEdge(p: S2Point, dimension: Int): Boolean {
        if (builder == null) return false  // Boolean output.
        if (!prevInside) setClippingState(S2BooleanOperationImpl.kSetInside, true)
        inputDimensions.add(dimension.toByte())
        builder.addEdge(p, p)
        prevInside = true
        return true;
    }

    // Processes an edge of dimension 0 (i.e., a point) from region A.
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun processEdge0(a_id: ShapeEdgeId, a: Edge, iter: CrossingIterator): Boolean {
        requireEQ(a.v0, a.v1)
        // When a region is inverted, all points and polylines are discarded.
        if (invertA != invertResult) {
            skipCrossings(a_id, iter)
            return true
        }
        val r = processPointCrossings(a_id, a.v0, iter)

        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region, using semi-open boundaries.
        var contained = inside xor invertB
        if (r.matchesPolygon && polygonModel != S2BooleanOperation.PolygonModel.SEMI_OPEN) {
            contained = (polygonModel == S2BooleanOperation.PolygonModel.CLOSED)
        }
        if (r.matchesPolyline) contained = true

        // The output of UNION includes duplicate values, so ensure that points are
        // not suppressed by other points.
        if (r.matchesPoint && !isUnion) contained = true

        // Test whether the point is contained after region B is inverted.
        if (contained == invertB) return true  // Don't exit early.
        return addPointEdge(a.v0, 0)
    }

    // Processes an edge of dimension 1 (i.e., a polyline edge) from region A.
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun processEdge1(a_id: ShapeEdgeId, a: Edge, iter: CrossingIterator): Boolean {
        // When a region is inverted, all points and polylines are discarded.
        if (invertA != invertResult) {
            skipCrossings(a_id, iter)
            return true
        }
        // Evaluate whether the start vertex should belong to the output, in case it
        // needs to be emitted as an isolated vertex.
        val r = processEdgeCrossings(a_id, a, iter)
        val a0_inside = isPolylineVertexInside(r.a0MatchesPolyline, r.a0MatchesPolygon)

        // Test whether the entire polyline edge should be emitted (or not emitted)
        // because it matches a polyline or polygon edge.
        inside = inside xor ((r.a0Crossings and 1) != 0)
        if (inside != isPolylineEdgeInside(r)) {
            inside = inside xor true   // Invert the inside_ state.
            ++r.a1Crossings  // Restore the correct (semi-open) state later.
        }

        // If neither edge adjacent to v0 was emitted, and this polyline contains
        // v0, and the other region contains v0, then emit an isolated vertex.
        if (!polylineLoopsHaveBoundaries &&
                a_id.edgeId == chainStart &&
                a.v0 == aShape.chainEdge(chainId, chainLimit - chainStart - 1).v1) {
            // This is the first vertex of a polyline loop, so we can't decide if it
            // is isolated until we process the last polyline edge.
            chainV0Emitted = inside
        } else if (isV0Isolated(a_id) && polylineContainsV0(a_id.edgeId, chainStart) && a0_inside) {
            if (!addPointEdge(a.v0, 1)) return false
        }

        // Test whether the entire edge or any part of it belongs to the output.
        if (inside || r.interiorCrossings > 0) {
            // Note: updates "inside_" to correspond to the state just before a1.
            if (!addEdge(a_id, a, 1 /*dimension*/, r.interiorCrossings)) {
                return false
            }
        }
        // Remember whether the edge portion just before "a1" was emitted, so that
        // we can decide whether "a1" need to be emitted as an isolated vertex.
        if (inside) v0EmittedMaxEdgeId = a_id.edgeId + 1

        // Verify that edge crossings are being counted correctly.
        inside = inside xor ((r.a1Crossings and 1) != 0)
        if (iter.crossingsComplete) {
            if (PreConditions.enabled) {
                checkEQ(S2ContainsPointQuery.makeS2ContainsPointQuery(iter.bIndex).contains(a.v1), inside xor invertB)
            }
        }

        // Special case to test whether the last vertex of a polyline should be
        // emitted as an isolated vertex.
        if (iter.crossingsComplete && isChainLastVertexIsolated(a_id) &&
                (polylineModel == S2BooleanOperation.PolylineModel.CLOSED ||
                        (!polylineLoopsHaveBoundaries && a.v1 == aShape.chainEdge(chainId, chainStart).v0)) &&
                isPolylineVertexInside(r.a1MatchesPolyline, r.a1MatchesPolygon)
        ) {
            if (!addPointEdge(a.v1, 1)) return false
        }
        return true
    }

    // Processes an edge of dimension 2 (i.e., a polygon edge) from region A.
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun processEdge2(aId: ShapeEdgeId, a: Edge, iter: CrossingIterator): Boolean {
        // In order to keep only one copy of any shared polygon edges, we only
        // output shared edges when processing the second region.
        val emitShared = (aRegionId == 1)

        // Degeneracies such as isolated vertices and sibling pairs can only be
        // created by intersecting CLOSED polygons or unioning OPEN polygons.
        val emitDegenerate = (polygonModel == S2BooleanOperation.PolygonModel.CLOSED && !invertA && !invertB) ||
                (polygonModel == S2BooleanOperation.PolygonModel.OPEN && invertA && invertB)

        val r = processEdgeCrossings(aId, a, iter)
        check(!r.matchesPolyline)
        inside = inside xor ((r.a0Crossings and 1) != 0)

        // If only one region is inverted, matching/sibling relations are reversed.
        // TODO(ericv): Update the following code to handle degenerate loops.
        check(!r.matchesPolygon || !r.matchesSibling)
        if (invertA != invertB) {
            val tmp = r.matchesPolygon
            r.matchesPolygon = r.matchesSibling
            r.matchesSibling = tmp
        }

        // Test whether the entire polygon edge should be emitted (or not emitted)
        // because it matches a polygon edge or its sibling.
        var newInside = inside

        // Shared edge are emitted only while processing the second region.
        if (r.matchesPolygon) newInside = emitShared

        // Sibling pairs are emitted only when degeneracies are desired.
        if (r.matchesSibling) newInside = emitDegenerate
        if (inside != newInside) {
            inside = inside xor true   // Invert the inside_ state.
            ++r.a1Crossings  // Restore the correct (semi-open) state later.
        }

        // Test whether the first vertex of this edge should be emitted as an
        // isolated degenerate vertex.
        if (aId.edgeId == chainStart) {
            chainV0Emitted = inside
        } else if (emitShared && emitDegenerate && r.a0MatchesPolygon && isV0Isolated(aId)) {
            if (!addPointEdge(a.v0, 2)) return false
        }

        // Test whether the entire edge or any part of it belongs to the output.
        if (inside || r.interiorCrossings > 0) {
            // Note: updates "inside_" to correspond to the state just before a1.
            if (!addEdge(aId, a, 2 /*dimension*/, r.interiorCrossings)) {
                return false
            }
        }

        // Remember whether the edge portion just before "a1" was emitted, so that
        // we can decide whether "a1" need to be emitted as an isolated vertex.
        if (inside) v0EmittedMaxEdgeId = aId.edgeId + 1
        inside = inside xor ((r.a1Crossings and 1) != 0)

        // Verify that edge crossings are being counted correctly.
        if (iter.crossingsComplete) {
            if (PreConditions.enabled) {
                checkEQ(S2ContainsPointQuery.makeS2ContainsPointQuery(iter.bIndex).contains(a.v1), inside xor invertB)
            }
        }

        // Special case to test whether the last vertex of a loop should be emitted
        // as an isolated degenerate vertex.
        if (emitShared && emitDegenerate && r.a1MatchesPolygon && iter.crossingsComplete && isChainLastVertexIsolated(aId)) {
            if (!addPointEdge(a.v1, 2)) return false;
        }
        return true;
    }

    // Skip any crossings that were not needed to determine the result.
    private fun skipCrossings(a_id: ShapeEdgeId, iter: CrossingIterator) {
        while (!iter.done(a_id)) iter.next()
    }

    // Returns a summary of the relationship between a point from region A and
    // a set of crossing edges from region B (see PointCrossingResult).
    private fun processPointCrossings(a_id: ShapeEdgeId, a0: S2Point, iter: CrossingIterator): PointCrossingResult {
        val r = PointCrossingResult()
        while (!iter.done(a_id)) {
            if (iter.b_dimension() == 0) {
                r.matchesPoint = true
            } else if (iter.b_dimension() == 1) {
                if (polylineEdgeContainsVertex(a0, iter)) {
                    r.matchesPolyline = true
                }
            } else {
                r.matchesPolygon = true
            }
            iter.next()
        }
        return r
    }

    // Returns a summary of the relationship between a test edge from region A and
    // a set of crossing edges from region B (see EdgeCrossingResult).
    //
    // NOTE(ericv): We could save a bit of work when matching polygon vertices by
    // passing in a flag saying whether this information is needed.  For example
    // if is only needed in ProcessEdge2 when (emit_shared && emit_degenerate).
    private fun processEdgeCrossings(a_id: ShapeEdgeId, a: Edge, iter: CrossingIterator): EdgeCrossingResult {
        val r = EdgeCrossingResult()
        if (iter.done(a_id)) return r

        // TODO(ericv): bool a_degenerate = (a.v0 == a.v1);
        while (!iter.done(a_id)) {
            // Polylines and polygons are not affected by point geometry.
            if (iter.b_dimension() == 0) {
                iter.next()
                continue
            }
            val b: Edge = iter.b_edge()
            if (iter.is_interior_crossing()) {
                // The crossing occurs in the edge interior.  The condition below says
                // that (1) polyline crossings don't affect polygon output, and (2)
                // subtracting a crossing polyline from a polyline has no effect.
                if (aDimension <= iter.b_dimension() && !(invertB != invertResult && iter.b_dimension() == 1)) {
                    val src_id = S2BooleanOperation.SourceId(bRegionId, iter.b_shape_id(), iter.b_edge_id())
                    addCrossing(Pair(src_id, iter.left_to_right()))
                }
                r.interiorCrossings += if (iter.b_dimension() == 1) 2 else 1
            } else if (iter.b_dimension() == 1) {
                // Polygons are not affected by polyline geometry.
                if (aDimension == 2) {
                    iter.next()
                    continue
                }
                if ((a.v0 == b.v0 && a.v1 == b.v1) || (a.v0 == b.v1 && a.v1 == b.v0)) {
                    r.matchesPolyline = true
                }
                if ((a.v0 == b.v0 || a.v0 == b.v1) && polylineEdgeContainsVertex(a.v0, iter)) {
                    r.a0MatchesPolyline = true
                }
                if ((a.v1 == b.v0 || a.v1 == b.v1) && polylineEdgeContainsVertex(a.v1, iter)) {
                    r.a1MatchesPolyline = true
                }
            } else {
                checkEQ(2, iter.b_dimension())
                if (a.v0 == b.v0 && a.v1 == b.v1) {
                    ++r.a0Crossings
                    r.matchesPolygon = true
                } else if (a.v0 == b.v1 && a.v1 == b.v0) {
                    ++r.a0Crossings
                    r.matchesSibling = true
                } else if (iter.is_vertex_crossing()) {
                    if (a.v0 == b.v0 || a.v0 == b.v1) {
                        ++r.a0Crossings
                    } else {
                        ++r.a1Crossings
                    }
                }
                if (a.v0 == b.v0 || a.v0 == b.v1) {
                    r.a0MatchesPolygon = true
                }
                if (a.v1 == b.v0 || a.v1 == b.v1) {
                    r.a1MatchesPolygon = true
                }
            }
            iter.next()
        }
        return r
    }

    // Returns true if the current point being processed (which must be a polyline
    // vertex) is contained by the opposite region (after inversion if "invert_b_"
    // is true).  "matches_polyline" and "matches_polygon" indicate whether the
    // vertex matches a polyline/polygon vertex of the opposite region.
    private fun isPolylineVertexInside(matches_polyline: Boolean, matches_polygon: Boolean): Boolean {
        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region using semi-open boundaries.
        var contained = inside xor invertB

        // For UNION the output includes duplicate polylines.  The test below
        // ensures that isolated polyline vertices are not suppressed by other
        // polyline vertices in the output.
        if (matches_polyline && !isUnion) {
            contained = true
        } else if (matches_polygon && polygonModel != S2BooleanOperation.PolygonModel.SEMI_OPEN) {
            contained = (polygonModel == S2BooleanOperation.PolygonModel.CLOSED)
        }
        // Finally, invert the result if the opposite region should be inverted.
        return contained xor invertB
    }

    // Returns true if the current polyline edge is contained by the opposite
    // region (after inversion if "invert_b_" is true).
    private fun isPolylineEdgeInside(r: EdgeCrossingResult): Boolean {
        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region using semi-open boundaries.
        var contained = inside xor invertB
        if (r.matchesPolyline && !isUnion) {
            contained = true
        } else if (r.matchesPolygon) {
            // In the SEMI_OPEN model, polygon sibling pairs cancel each other and
            // have no effect on point or edge containment.
            if (!(r.matchesSibling && polygonModel == S2BooleanOperation.PolygonModel.SEMI_OPEN)) {
                contained = (polygonModel != S2BooleanOperation.PolygonModel.OPEN);
            }
        } else if (r.matchesSibling) {
            contained = (polygonModel == S2BooleanOperation.PolygonModel.CLOSED)
        }
        // Finally, invert the result if the opposite region should be inverted.
        return contained xor invertB
    }

    // Returns true if the vertex "v" is contained by the polyline edge referred
    // to by the CrossingIterator "it", taking into account the PolylineModel.
    //
    // REQUIRES: it.b_dimension() == 1
    // REQUIRES: "v" is an endpoint of it.b_edge()
    private fun polylineEdgeContainsVertex(v: S2Point, iter: CrossingIterator): Boolean {
        requireEQ(1, iter.b_dimension())
        require(iter.b_edge().v0 == v || iter.b_edge().v1 == v)

        // Closed polylines contain all their vertices.
        if (polylineModel == S2BooleanOperation.PolylineModel.CLOSED) return true

        // Note that the code below is structured so that it.b_edge() is not usually
        // needed (since accessing the edge can be relatively expensive).
        val b_chain = iter.b_chain_info()
        val b_edge_id = iter.b_edge_id()

        // The last polyline vertex is never contained.  (For polyline loops, it is
        // sufficient to treat the first vertex as begin contained.)  This case also
        // handles degenerate polylines (polylines with one edge where v0 == v1),
        // which do not contain any points.
        if (b_edge_id == b_chain.limit - 1 && v == iter.b_edge().v1) return false

        // Otherwise all interior vertices are contained.  The first polyline
        // vertex is contained if either the polyline model is not OPEN, or the
        // polyline forms a loop and polyline_loops_have_boundaries_ is false.
        if (polylineContainsV0(b_edge_id, b_chain.start)) return true
        if (v != iter.b_edge().v0) return true
        if (polylineLoopsHaveBoundaries) return false
        return v == iter.b_shape().chainEdge(b_chain.chain_id, b_chain.limit - b_chain.start - 1).v1
    }

    companion object {

        private val logger = KotlinLogging.logger(CrossingProcessor::class.java.name)

    }

}
