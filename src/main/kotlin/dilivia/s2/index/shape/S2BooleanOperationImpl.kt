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

import com.google.common.collect.ComparisonChain
import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkLE
import dilivia.PreConditions.checkLT
import dilivia.PreConditions.requireEQ
import dilivia.collections.assign
import dilivia.collections.lowerBound
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.s2.S2CellId
import dilivia.s2.S2Error
import dilivia.s2.S2Predicates
import dilivia.s2.S2TextParser
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.graph.VertexId
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.index.CrossingType
import dilivia.s2.shape.ShapeEdge
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min

class S2BooleanOperationImpl(val op: S2BooleanOperation) {

    // The S2Builder used to construct the output.
    private var builder: S2Builder? = null

    // A vector specifying the dimension of each edge added to S2Builder.
    private val inputDimensions: MutableList<Byte> = mutableListOf()

    // The set of all input edge crossings, which is used by EdgeClippingLayer
    // to construct the clipped output polygon.
    private val inputCrossings: InputEdgeCrossings = ArrayList<Pair<InputEdgeId, CrossingInputEdge>>()

    // A vector containing all pairs of crossing edges from the two input
    // regions (including edge pairs that share a common vertex).  The first
    // element of each pair is an edge from "index_crossings_first_region_id_",
    // while the second element of each pair is an edge from the other region.
    private val indexCrossings: IndexCrossings = mutableListOf()

    // Indicates that the first element of each crossing edge pair in
    // "index_crossings_" corresponds to an edge from the given region.
    // This field is negative if index_crossings_ has not been computed yet.
    private var indexCrossingsFirstRegionId: Int = -1

    // Temporary storage used in GetChainStarts(), declared here to avoid
    // repeatedly allocating memory.
    private val tmpCrossings: IndexCrossings = mutableListOf()

    fun build(error: S2Error): Boolean {
        logger.trace {
            "build | opType = ${op.opType}, " +
                    "\noptions = ${op.options}," +
                    "\na       = ${S2TextParser.toString(op.a)}" +
                    "\nb       = ${S2TextParser.toString(op.b)}"
        }
        error.clear()
        if (isBooleanOutput()) {
            // BuildOpType() returns true if and only if the result has no edges.
            op.resultEmpty!!.value = buildOpType(op.opType) && !isFullPolygonResult()
            return true
        }
        // TODO(ericv): Rather than having S2Builder split the edges, it would be
        // faster to call AddVertex() in this class and have a new S2Builder
        // option that increases the edge_snap_radius_ to account for errors in
        // the intersection point (the way that split_crossing_edges does).
        val options = S2Builder.Options(op.options.snapFunction)
        options.splitCrossingEdges = true

        // TODO(ericv): Ideally idempotent() should be true, but existing clients
        // expect vertices closer than the full "snap_radius" to be snapped.
        options.idempotent = false
        builder = S2Builder(options).also { b ->
            b.startLayer(EdgeClippingLayer(op.layers, inputDimensions, inputCrossings))
            // Add a predicate that decides whether a result with no polygon edges should
            // be interpreted as the empty polygon or the full polygon.
            b.addIsFullPolygonPredicate { _: Graph, _: S2Error -> isFullPolygonResult() }
        }

        buildOpType(op.opType)
        return builder?.build(error) ?: throw IllegalStateException()
    }

    private fun isBooleanOutput(): Boolean {
        return op.resultEmpty != null
    }

    // All of the methods below support "early exit" in the case of boolean
    // results by returning "false" as soon as the result is known to be
    // non-empty.

    // Clips the boundary of A to the interior of the opposite region B and adds
    // the resulting edges to the output.  Optionally, any combination of region
    // A, region B, and the result may be inverted, which allows operations such
    // as union and difference to be implemented.
    //
    // Note that when an input region is inverted with respect to the output
    // (e.g., invert_a != invert_result), all polygon edges are reversed and all
    // points and polylines are discarded, since the complement of such objects
    // cannot be represented.  (If you want to compute the complement of points
    // or polylines, you can use S2LaxPolygonShape to represent your geometry as
    // degenerate polygons instead.)
    //
    // This method must be called an even number of times (first to clip A to B
    // and then to clip B to A), calling DoneBoundaryPair() after each pair.
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun addBoundary(
        aRegionId: Int,
        invertA: Boolean,
        invertB: Boolean,
        invertResult: Boolean,
        aChainStarts: List<ShapeEdgeId>,
        cp: CrossingProcessor
    ): Boolean {
        logger.trace { "addBoundary | aRegionId = $aRegionId, invertA = $invertA, invertB = $invertB, invertResult = $invertResult" }

        val aIndex: S2ShapeIndex = op.getShapeIndex(aRegionId)
        val bIndex: S2ShapeIndex = op.getShapeIndex(1 - aRegionId)
        if (!getIndexCrossings(aRegionId)) {
            return false
        }
        cp.startBoundary(aRegionId, invertA, invertB, invertResult)

        // Walk the boundary of region A and build a list of all edge crossings.
        // We also keep track of whether the current vertex is inside region B.
        var nextStartIdx = 0
        val nextCrossing = CrossingIterator(bIndex, indexCrossings, true /*crossings_complete*/)
        var nextId = minOf(aChainStarts[nextStartIdx], nextCrossing.a_id())
        while (nextId != kSentinel) {
            logger.trace { "addBoundary | nextId = $nextId" }
            val aShapeId = nextId.shapeId
            val aShape = aIndex.shape(aShapeId)!!
            cp.startShape(aShape)
            while (nextId.shapeId == aShapeId) {
                logger.trace { "addBoundary | process edge = $nextId" }
                // TODO(ericv): Special handling of dimension 0?  Can omit most of this
                // code, including the loop, since all chains are of length 1.
                var edgeId = nextId.edgeId
                val chainPosition = aShape.chainPosition(edgeId)
                val chainId = chainPosition.chainId
                logger.trace { "chain id = $chainPosition" }

                val chain = aShape.chain(chainId)
                val startInside = (nextId == aChainStarts[nextStartIdx])
                if (startInside) ++nextStartIdx
                cp.startChain(chainId, chain, startInside)
                val chainLimit = chain.start + chain.length
                while (edgeId < chainLimit) {
                    val aId = ShapeEdgeId(aShapeId, edgeId)
                    check(cp.inside() || nextCrossing.a_id() == aId)
                    if (!cp.processEdge(aId, nextCrossing)) {
                        return false
                    }
                    if (cp.inside()) {
                        ++edgeId
                    } else if (nextCrossing.a_id().shapeId == aShapeId && nextCrossing.a_id().edgeId < chainLimit) {
                        edgeId = nextCrossing.a_id().edgeId
                    } else {
                        break;
                    }
                }
                nextId = minOf(aChainStarts[nextStartIdx], nextCrossing.a_id())
            }
        }
        return true
    }

    // Returns the first edge of each edge chain from "a_region_id" whose first
    // vertex is contained by opposite region's polygons (using the semi-open
    // boundary model).  Each input region and the result region are inverted as
    // specified (invert_a, invert_b, and invert_result) before testing for
    // containment.  The algorithm uses these "chain starts" in order to clip the
    // boundary of A to the interior of B in an output-senstive way.
    //
    // This method supports "early exit" in the case where a boolean predicate is
    // being evaluated and the algorithm discovers that the result region will be
    // non-empty.
    private fun getChainStarts(
        aRegionId: Int,
        invertA: Boolean,
        invertB: Boolean,
        invertResult: Boolean,
        cp: CrossingProcessor,
        chainStarts: MutableList<ShapeEdgeId>
    ): Boolean {
        val aIndex: S2ShapeIndex = op.getShapeIndex(aRegionId)
        val bIndex: S2ShapeIndex = op.getShapeIndex(1 - aRegionId)

        logger.trace { "getChainStarts | aRegionId = $aRegionId, invertA = $invertA, invertB = $invertB, invertResult = $invertResult" }

        if (isBooleanOutput()) {
            // If boolean output is requested, then we use the CrossingProcessor to
            // determine whether the first edge of each chain will be emitted to the
            // output region.  This lets us terminate the operation early in many
            // cases.
            cp.startBoundary(aRegionId, invertA, invertB, invertResult)
        }

        // If region B has no two-dimensional shapes and is not inverted, then by
        // definition no chain starts are contained.  However if boolean output is
        // requested then we check for containment anyway, since as a side effect we
        // may discover that the result region is non-empty and terminate the entire
        // operation early.
        val bHasInterior = hasInterior(bIndex);
        if (bHasInterior || invertB || isBooleanOutput()) {
            val query = S2ContainsPointQuery.makeS2ContainsPointQuery(bIndex)
            val numShapeIds = aIndex.nextNewShapeId()
            for (shapeId in 0 until numShapeIds) {
                val aShape = aIndex.shape(shapeId) ?: continue

                // If region A is being subtracted from region B, points and polylines
                // in region A can be ignored since these shapes never contribute to the
                // output (they can only remove edges from region B).
                if (invertA != invertResult && aShape.dimension < 2) continue

                if (isBooleanOutput()) cp.startShape(aShape)
                val numChains = aShape.numChains
                for (chainId in 0 until numChains) {
                    val chain = aShape.chain(chainId)
                    if (chain.length == 0) continue
                    val a = ShapeEdge(shapeId, chain.start, aShape.chainEdge(chainId, 0))
                    val inside = (bHasInterior && query.contains(a.v0)) != invertB
                    if (inside) {
                        chainStarts.add(ShapeEdgeId(shapeId, chain.start))
                    }
                    if (isBooleanOutput()) {
                        cp.startChain(chainId, chain, inside)
                        if (!processIncidentEdges(a, query, cp)) {
                            logger.trace { "getChainStarts | false" }
                            return false
                        }
                    }
                }
            }
        }
        chainStarts.add(kSentinel)

        logger.trace { "getChainStarts | chainStarts = $chainStarts" }
        return true
    }

    private fun processIncidentEdges(
        a: ShapeEdge,
        query: S2ContainsPointQuery<S2ShapeIndex>,
        cp: CrossingProcessor
    ): Boolean {
        tmpCrossings.clear()
        query.visitIncidentEdges(a.v0) { b: ShapeEdge -> addIndexCrossing(a, b, false /*is_interior*/, tmpCrossings) }
        // Fast path for the common case where there are no incident edges.  We
        // return false (terminating early) if the first chain edge will be emitted.
        if (tmpCrossings.isEmpty()) {
            return !cp.inside()
        }
        // Otherwise we invoke the full CrossingProcessor logic to determine whether
        // the first chain edge will be emitted.
        if (tmpCrossings.size > 1) {
            tmpCrossings.sort()
            // VisitIncidentEdges() should not generate any duplicate values.
            //S2_DCHECK(std::adjacent_find(tmp_crossings_.begin(), tmp_crossings_.end()) == tmp_crossings_.end());
        }
        tmpCrossings.add(IndexCrossing(kSentinel, kSentinel))
        val nextCrossing = CrossingIterator(query.index(), tmpCrossings, false /*crossings_complete*/)
        return cp.processEdge(a.id, nextCrossing)
    }

    // Initialize index_crossings_ to the set of crossing edge pairs such that the
    // first element of each pair is an edge from "region_id".
    //
    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun getIndexCrossings(regionId: Int): Boolean {
        logger.trace { "getIndexCrossings | regionId = $regionId, indexCrossingsFirstRegionId = $indexCrossingsFirstRegionId\na = ${S2TextParser.toString(op.a)}\nb = ${S2TextParser.toString(op.b)}" }
        if (regionId == indexCrossingsFirstRegionId) return true
        if (indexCrossingsFirstRegionId < 0) {
            requireEQ(regionId, 0)  // For efficiency, not correctness.
            if (!S2CrossingEdgePairsScanner.visitCrossingEdgePairs(
                    op.a,
                    op.b,
                    CrossingType.ALL
                ) { a: ShapeEdge, b: ShapeEdge, isInterior: Boolean ->
                    // For all supported operations (union, intersection, and
                    // difference), if the input edges have an interior crossing
                    // then the output is guaranteed to have at least one edge.
                    if (isInterior && isBooleanOutput()) false else addIndexCrossing(a, b, isInterior, indexCrossings)
                }
            ) {
                return false
            }
            if (indexCrossings.size > 1) {
                indexCrossings.sortAndRemoveDuplicates()
            }
            // Add a sentinel value to simplify the loop logic.
            indexCrossings.add(IndexCrossing(kSentinel, kSentinel))
            indexCrossingsFirstRegionId = 0
        }
        if (regionId != indexCrossingsFirstRegionId) {
            for (crossing in indexCrossings) {
                val tmp = crossing.a
                crossing.a = crossing.b
                crossing.b = tmp
                // The following predicates get inverted when the edges are swapped.
                crossing.leftToRight = crossing.leftToRight xor true
                crossing.isVertexCrossing = crossing.isVertexCrossing xor true
            }
            indexCrossings.sort()
            indexCrossingsFirstRegionId = regionId
        }
        return true
    }

    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun addBoundaryPair(
        invertA: Boolean,
        invertB: Boolean,
        invertResult: Boolean,
        cp: CrossingProcessor
    ): Boolean {
        logger.trace { "addBoundaryPair | invertA = $invertA, invertB = $invertB, invertResult = $invertResult" }

        // Optimization: if the operation is DIFFERENCE or SYMMETRIC_DIFFERENCE,
        // it is worthwhile checking whether the two regions are identical (in which
        // case the output is empty).
        //
        // TODO(ericv): When boolean output is requested there are other quick
        // checks that could be done here, such as checking whether a full cell from
        // one S2ShapeIndex intersects a non-empty cell of the other S2ShapeIndex.
        if (op.opType == S2BooleanOperation.OpType.DIFFERENCE || op.opType == S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE) {
            if (areRegionsIdentical()) {
                return true
            }
        } else if (!isBooleanOutput()) {
        }

        val aStarts: MutableList<ShapeEdgeId> = mutableListOf()
        val bStarts: MutableList<ShapeEdgeId> = mutableListOf()
        if (!getChainStarts(0, invertA, invertB, invertResult, cp, aStarts) ||
            !getChainStarts(1, invertB, invertA, invertResult, cp, bStarts) ||
            !addBoundary(0, invertA, invertB, invertResult, aStarts, cp) ||
            !addBoundary(1, invertB, invertA, invertResult, bStarts, cp)
        ) {
            return false
        }
        if (!isBooleanOutput()) cp.doneBoundaryPair()
        return true
    }

    private fun areRegionsIdentical(): Boolean {
        val a = op.a
        val b = op.b
        if (a == b) return true
        val num_shape_ids = a.nextNewShapeId()
        if (num_shape_ids != b.nextNewShapeId()) return false
        for (s in 0 until num_shape_ids) {
            val a_shape = a.shape(s)!!
            val b_shape = b.shape(s)!!
            if (a_shape.dimension != b_shape.dimension) return false
            if (a_shape.dimension == 2) {
                val a_ref = a_shape.getReferencePoint()
                val b_ref = b_shape.getReferencePoint()
                if (a_ref.point != b_ref.point) return false
                if (a_ref.contained != b_ref.contained) return false
            }
            val num_chains = a_shape.numChains
            if (num_chains != b_shape.numChains) return false
            for (c in 0 until num_chains) {
                val a_chain = a_shape.chain(c)
                val b_chain = b_shape.chain(c)
                checkEQ(a_chain.start, b_chain.start)
                if (a_chain.length != b_chain.length) return false
                for (i in 0 until a_chain.length) {
                    val a_edge = a_shape.chainEdge(c, i)
                    val b_edge = b_shape.chainEdge(c, i)
                    if (a_edge.v0 != b_edge.v0) return false
                    if (a_edge.v1 != b_edge.v1) return false
                }
            }
        }
        return true
    }

    // Supports "early exit" in the case of boolean results by returning false
    // as soon as the result is known to be non-empty.
    private fun buildOpType(opType: S2BooleanOperation.OpType): Boolean {
        // CrossingProcessor does the real work of emitting the output edges.
        val cp = CrossingProcessor(
            op.options.polygonModel,
            op.options.polylineModel,
            op.options.polylineLoopsHaveBoundaries,
            builder,
            inputDimensions,
            inputCrossings
        )
        return when (opType) {
            // A | B == ~(~A & ~B)
            S2BooleanOperation.OpType.UNION -> addBoundaryPair(true, true, true, cp)

            // A & B
            S2BooleanOperation.OpType.INTERSECTION -> addBoundaryPair(false, false, false, cp)

            // A - B = A & ~B
            S2BooleanOperation.OpType.DIFFERENCE -> addBoundaryPair(false, true, false, cp)

            // Compute the union of (A - B) and (B - A).
            S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE ->
                (addBoundaryPair(false, true, false, cp) &&
                        addBoundaryPair(true, false, false, cp))
        }
    }

    // Given a polygon edge graph containing only degenerate edges and sibling edge
    // pairs, the purpose of this function is to decide whether the polygon is empty
    // or full except for the degeneracies, i.e. whether the degeneracies represent
    // shells or holes.
    private fun isFullPolygonResult(): Boolean {
        // If there are no edges of dimension 2, the result could be either the
        // empty polygon or the full polygon.  Note that this is harder to determine
        // than you might think due to snapping.  For example, the union of two
        // non-empty polygons can be empty, because both polygons consist of tiny
        // loops that are eliminated by snapping.  Similarly, even if two polygons
        // both contain a common point their intersection can still be empty.
        //
        // We distinguish empty from full results using two heuristics:
        //
        //  1. We compute a bit mask representing the subset of the six S2 cube faces
        //     intersected by each input geometry, and use this to determine if only
        //     one of the two results is possible.  (This test is very fast.)  Note
        //     that snapping will never cause the result to cover an entire extra
        //     cube face because the maximum allowed snap radius is too small.
        checkLE(SnapFunction.kMaxSnapRadius().degrees(), 70.0)
        //
        //  2. We compute the area of each input geometry, and use this to bound the
        //     minimum and maximum area of the result.  If only one of {0, 4*Pi} is
        //     possible then we are done.  If neither is possible then we choose the
        //     one that is closest to being possible (since snapping can change the
        //     result area).  Both results are possible only when computing the
        //     symmetric difference of two regions of area 2*Pi each, in which case we
        //     must resort to additional heuristics (see below).
        //
        // TODO(ericv): Implement a predicate that uses the results of edge snapping
        // directly, rather than computing areas.  This would not only be much faster
        // but would also allows all cases to be handled 100% robustly.
        val a: S2ShapeIndex = op.a
        val b: S2ShapeIndex = op.b
        return when (op.opType) {
            S2BooleanOperation.OpType.UNION -> isFullPolygonUnion(a, b)
            S2BooleanOperation.OpType.INTERSECTION -> isFullPolygonIntersection(a, b)
            S2BooleanOperation.OpType.DIFFERENCE -> isFullPolygonDifference(a, b)
            S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE -> isFullPolygonSymmetricDifference(a, b)
        }
    }

    private fun isFullPolygonUnion(a: S2ShapeIndex, b: S2ShapeIndex): Boolean {
        // See comments in IsFullPolygonResult().  The most common case is that
        // neither input polygon is empty but the result is empty due to snapping.

        // The result can be full only if the union of the two input geometries
        // intersects all six faces of the S2 cube.  This test is fast.
        if ((getFaceMask(a) or getFaceMask(b)) != kAllFacesMask) return false

        // The union area satisfies:
        //
        //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
        //
        // where A, B can refer to a polygon or its area.  We then choose the result
        // that assumes the smallest amount of error.
        val aArea = a.getArea()
        val bArea = b.getArea()
        val minArea = max(aArea, bArea)
        val maxArea = min(4 * M_PI, aArea + bArea)
        return minArea > 4 * M_PI - maxArea
    }

    private fun isFullPolygonIntersection(a: S2ShapeIndex, b: S2ShapeIndex): Boolean {
        // See comments in IsFullPolygonResult().  By far the most common case is
        // that the result is empty.

        // The result can be full only if each of the two input geometries
        // intersects all six faces of the S2 cube.  This test is fast.
        if ((getFaceMask(a) and getFaceMask(b)) != kAllFacesMask) return false

        // The intersection area satisfies:
        //
        //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
        //
        // where A, B can refer to a polygon or its area.  We then choose the result
        // that assumes the smallest amount of error.
        val a_area = a.getArea()
        val b_area = b.getArea()
        val min_area = max(0.0, a_area + b_area - 4 * M_PI)
        val max_area = min(a_area, b_area)
        return min_area > 4 * M_PI - max_area
    }

    private fun isFullPolygonDifference(a: S2ShapeIndex, b: S2ShapeIndex): Boolean {
        // See comments in IsFullPolygonResult().  By far the most common case is
        // that the result is empty.

        // The result can be full only if each cube face is intersected by the first
        // geometry.  (The second geometry is irrelevant, since for example it could
        // consist of a tiny loop on each S2 cube face.)  This test is fast.
        if (getFaceMask(a) != kAllFacesMask) return false

        // The difference area satisfies:
        //
        //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
        //
        // where A, B can refer to a polygon or its area.  We then choose the result
        // that assumes the smallest amount of error.
        val aArea = a.getArea()
        val bArea = b.getArea()
        val minArea = max(0.0, aArea - bArea)
        val maxArea = min(aArea, 4 * M_PI - bArea)
        return minArea > 4 * M_PI - maxArea
    }

    private fun isFullPolygonSymmetricDifference(a: S2ShapeIndex, b: S2ShapeIndex): Boolean {
        // See comments in IsFullPolygonResult().  By far the most common case is
        // that the result is empty.

        // The result can be full only if the union of the two input geometries
        // intersects all six faces of the S2 cube.  This test is fast.
        val aMask = getFaceMask(a)
        val bMask = getFaceMask(b)
        if ((aMask or bMask) != kAllFacesMask) return false

        // The symmetric difference area satisfies:
        //
        //   |A - B| <= SymmetricDifference(A, B) <= 4*Pi - |4*Pi - (A + B)|
        //
        // where A, B can refer to a polygon or its area.
        val aArea = a.getArea()
        val bArea = b.getArea()
        val minArea = abs(aArea - bArea)
        val maxArea = 4 * M_PI - abs(4 * M_PI - (aArea + bArea))

        // Now we choose the result that assumes the smallest amount of error
        // (min_area in the empty case, and (4*Pi - max_area) in the full case).
        // However in the case of symmetric difference these two errors may be equal,
        // meaning that the result is ambiguous.  This happens when both polygons have
        // area 2*Pi.  Furthermore, this can happen even when the areas are not
        // exactly 2*Pi due to snapping and area calculation errors.
        //
        // To determine whether the result is ambiguous, we compute a rough estimate
        // of the maximum expected area error (including errors due to snapping),
        // using the worst-case error bound for a hemisphere defined by 4 vertices.
        val edgeSnapRadius =
            op.options.snapFunction.snapRadius() + S2EdgeCrossings.kIntersectionError  // split_crossing_edges
        val hemisphereAreaError = 2 * M_PI * edgeSnapRadius.radians + 40 * DoubleType.epsilon  // GetCurvatureMaxError

        // The following sign is the difference between the error needed for an empty
        // result and the error needed for a full result.  It is negative if an
        // empty result is possible, positive if a full result is possible, and zero
        // if both results are possible.
        val errorSign = minArea - (4 * M_PI - maxArea)
        if (abs(errorSign) <= hemisphereAreaError) {
            // Handling the ambiguous case correctly requires a more sophisticated
            // algorithm (see below), but we can at least handle the simple cases by
            // testing whether both input geometries intersect all 6 cube faces.  If
            // not, then the result is definitely full.
            if ((aMask and bMask) != kAllFacesMask) return true

            // Otherwise both regions have area 2*Pi and intersect all 6 cube faces.
            // We choose "empty" in this case under the assumption that it is more
            // likely that the user is computing the difference between two nearly
            // identical polygons.
            //
            // TODO(ericv): Implement a robust algorithm based on examining the edge
            // snapping results directly, or alternatively add another heuristic (such
            // as testing containment of random points, or using a larger bit mask in
            // the tests above, e.g. a 24-bit mask representing all level 1 cells).
            return false
        }
        return errorSign > 0
    }

    // CrossingInputEdge represents an input edge B that crosses some other input
    // edge A.  It stores the input edge id of edge B and also whether it crosses
    // edge A from left to right (or vice versa).
    // Indicates that input edge "input_id" crosses another edge (from left to
    // right if "left_to_right" is true).
    data class CrossingInputEdge(val inputId: InputEdgeId, val leftToRight: Boolean) : Comparable<InputEdgeId> {

        operator fun compareTo(other: CrossingInputEdge): Int {
            return inputId.compareTo(other.inputId)
        }

        override operator fun compareTo(other: InputEdgeId): Int = inputId.compareTo(other)

    }

    // Given two input edges A and B that intersect, suppose that A maps to a
    // chain of snapped edges A_0, A_1, ..., A_m and B maps to a chain of snapped
    // edges B_0, B_1, ..., B_n.  CrossingGraphEdge represents an edge from chain
    // B that shares a vertex with chain A.  It is used as a temporary data
    // representation while processing chain A.  The arguments are:
    //
    //   "id" - the Graph::EdgeId of an edge from chain B.
    //   "a_index" - the index of the vertex (A_i) that is shared with chain A.
    //   "outgoing" - true if the shared vertex is the first vertex of the B edge.
    //   "dst" - the Graph::VertexId of the vertex that is not shared with chain A.
    //
    // Note that if an edge from the B chain shares both vertices with the A
    // chain, there will be two entries: an outgoing edge that treats its first
    // vertex as being shared, and an incoming edge that treats its second vertex
    // as being shared.
    data class CrossingGraphEdge(val id: EdgeId, val aIndex: Int, val outgoing: Boolean, val dst: VertexId)

    companion object {

        private val logger = KotlinLogging.logger(S2BooleanOperationImpl::class.java.name)

        // A collection of special InputEdgeIds that allow the GraphEdgeClipper state
        // modifications to be inserted into the list of edge crossings.
        const val kSetInside: InputEdgeId = -1
        const val kSetInvertB: InputEdgeId = -2
        const val kSetReverseA: InputEdgeId = -3

        // A bit mask representing all six faces of the S2 cube.
        const val kAllFacesMask: Int = 0x3f

        // kSentinel is a sentinel value used to mark the end of vectors.
        val kSentinel: ShapeEdgeId = ShapeEdgeId(Int.MAX_VALUE, 0)

        fun hasInterior(index: S2ShapeIndex): Boolean {
            for (s in index.nextNewShapeId() - 1 downTo 0) {
                val shape = index.shape(s)
                if (shape != null && shape.dimension == 2) return true
            }
            return false
        }

        fun addIndexCrossing(a: ShapeEdge, b: ShapeEdge, isInterior: Boolean, crossings: IndexCrossings): Boolean {
            crossings.add(IndexCrossing(a.id, b.id))
            val crossing = crossings.last()
            if (isInterior) {
                crossing.isInteriorCrossing = true
                if (S2Predicates.sign(a.v0, a.v1, b.v0) > 0) {
                    crossing.leftToRight = true
                }
            } else {
                // TODO(ericv): This field isn't used unless one shape is a polygon and
                // the other is a polyline or polygon, but we don't have the shape
                // dimension information readily available here.
                if (S2EdgeCrossings.vertexCrossing(a.v0, a.v1, b.v0, b.v1)) {
                    crossing.isVertexCrossing = true
                }
            }
            return true;  // Continue visiting.
        }

        // Returns a bit mask indicating which of the 6 S2 cube faces intersect the
        // index contents.
        fun getFaceMask(index: S2ShapeIndex): Int {
            var mask: Int = 0
            val it = index.newIterator(InitialPosition.BEGIN)
            while (!it.done()) {
                val face = it.id().face()
                mask = mask or (1 shl face)
                it.seek(S2CellId.fromFace(face + 1).rangeMin())
            }
            return mask
        }

        // Returns a vector of EdgeIds sorted by input edge id.  When more than one
        // output edge has the same input edge id (i.e., the input edge snapped to a
        // chain of edges), the edges are sorted so that they form a directed edge
        // chain.
        //
        // This function could possibily be moved to S2Builder::Graph, but note that
        // it has special requirements.  Namely, duplicate edges and sibling pairs
        // must be kept in order to ensure that every output edge corresponds to
        // exactly one input edge.  (See also S2Builder::Graph::GetInputEdgeOrder.)
        fun getInputEdgeChainOrder(g: Graph, input_ids: List<InputEdgeId>): List<EdgeId> {
            require(g.options.edgeType == EdgeType.DIRECTED);
            require(g.options.duplicateEdges == DuplicateEdges.KEEP)
            require(g.options.siblingPairs == SiblingPairs.KEEP)

            // First, sort the edges so that the edges corresponding to each input edge
            // are consecutive.  (Each input edge was snapped to a chain of output
            // edges, or two chains in the case of undirected input edges.)
            val order = g.getInputEdgeOrder(input_ids)

            // Now sort the group of edges corresponding to each input edge in edge
            // chain order (e.g.  AB, BC, CD).
            val vmap = mutableListOf<Pair<VertexId, EdgeId>>()     // Map from source vertex to edge id.
            val indegree = mutableListOf<Int>()
            indegree.assign(g.numVertices, 0)  // Restricted to current input edge.
            var end = 0
            var begin = 0
            while (begin < order.size) {
                // Gather the edges that came from a single input edge.
                val input_id = input_ids[order[begin]]
                end = begin
                while (end < order.size) {
                    if (input_ids[order[end]] != input_id) break
                    ++end
                }
                if (end - begin == 1) {
                    begin = end
                    continue
                }

                // Build a map from the source vertex of each edge to its edge id,
                // and also compute the indegree at each vertex considering only the edges
                // that came from the current input edge.
                for (i in begin until end) {
                    val e = order[i]
                    vmap.add(Pair(g.edge(e).first, e))
                    indegree[g.edge(e).second] += 1
                }
                val vmapComparator = Comparator<Pair<VertexId, EdgeId>> { o1, o2 ->
                    ComparisonChain.start().compare(o1.first, o2.first).compare(o1.second, o2.second).result()
                }
                vmap.sortWith(vmapComparator)

                // Find the starting edge for building the edge chain.
                var next = g.numEdges
                for (i in begin until end) {
                    val e = order[i]
                    if (indegree[g.edge(e).first] == 0) next = e
                }
                // Build the edge chain.
                var i = begin
                while (true) {
                    order[i] = next
                    val v = g.edge(next).second
                    indegree[v] = 0  // Clear as we go along.
                    if (++i == end) break
                    val outIdx = vmap.lowerBound(0, vmap.size, Pair(v, 0), vmapComparator)
                    checkLT(outIdx, vmap.size)
                    val out = vmap[outIdx]
                    checkEQ(v, out.first)
                    next = out.second
                }
                vmap.clear()
                begin = end
            }
            return order
        }

    }

}
