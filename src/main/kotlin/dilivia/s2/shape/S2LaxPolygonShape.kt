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

import dilivia.PreConditions.requireLT
import dilivia.collections.lowerBound
import dilivia.s2.S2Point
import dilivia.s2.region.S2PointLoopSpan
import dilivia.s2.region.S2Polygon

typealias Loop = List<S2Point>

// S2LaxPolygonShape represents a region defined by a collection of zero or
// more closed loops.  The interior is the region to the left of all loops.
// This is similar to S2Polygon::Shape except that this class supports
// polygons with degeneracies.  Degeneracies are of two types: degenerate
// edges (from a vertex to itself) and sibling edge pairs (consisting of two
// oppositely oriented edges).  Degeneracies can represent either "shells" or
// "holes" depending on the loop they are contained by.  For example, a
// degenerate edge or sibling pair contained by a "shell" would be interpreted
// as a degenerate hole.  Such edges form part of the boundary of the polygon.
//
// Loops with fewer than three vertices are interpreted as follows:
//  - A loop with two vertices defines two edges (in opposite directions).
//  - A loop with one vertex defines a single degenerate edge.
//  - A loop with no vertices is interpreted as the "full loop" containing
//    all points on the sphere.  If this loop is present, then all other loops
//    must form degeneracies (i.e., degenerate edges or sibling pairs).  For
//    example, two loops {} and {X} would be interpreted as the full polygon
//    with a degenerate single-point hole at X.
//
// S2LaxPolygonShape does not have any error checking, and it is perfectly
// fine to create S2LaxPolygonShape objects that do not meet the requirements
// below (e.g., in order to analyze or fix those problems).  However,
// S2LaxPolygonShapes must satisfy some additional conditions in order to
// perform certain operations:
//
//  - In order to be valid for point containment tests, the polygon must
//    satisfy the "interior is on the left" rule.  This means that there must
//    not be any crossing edges, and if there are duplicate edges then all but
//    at most one of thm must belong to a sibling pair (i.e., the number of
//    edges in opposite directions must differ by at most one).
//
//  - To be valid for boolean operations (S2BooleanOperation), degenerate
//    edges and sibling pairs cannot coincide with any other edges.  For
//    example, the following situations are not allowed:
//
//      {AA, AA}      // degenerate edge coincides with another edge
//      {AA, AB}      // degenerate edge coincides with another edge
//      {AB, BA, AB}  // sibling pair coincides with another edge
//
// Note that S2LaxPolygonShape is must faster to initialize and is more
// compact than S2Polygon, but unlike S2Polygon it does not have any built-in
// operations.  Instead you should use S2ShapeIndex operations
// (S2BooleanOperation, S2ClosestEdgeQuery, etc).
class S2LaxPolygonShape : S2Shape {

    private var numLoops: Int = 0
    private var vertices: Array<S2Point> = emptyArray()

    // array of size (num_loops + 1) where element "i"
    // is the total number of vertices in loops 0..i-1.
    private var cumulativeVertices: IntArray = IntArray(1)

    // Constructs an empty polygon.
    constructor()

    // Constructs an S2LaxPolygonShape from the given vertex loops.
    constructor(loops: List<Loop>) {
        init(loops)
    }

    // Constructs an S2LaxPolygonShape from an S2Polygon, by copying its data.
    // Full and empty S2Polygons are supported.
    constructor(polygon: S2Polygon) {
        init(polygon)
    }

    // Initializes an S2LaxPolygonShape from the given vertex loops.
    fun init(loops: List<Loop>) {
        initSpans(loops.map { S2PointLoopSpan(it) })
    }

    // Initializes an S2LaxPolygonShape from an S2Polygon, by copying its data.
    // Full and empty S2Polygons are supported.
    fun init(polygon: S2Polygon) {
        val spans = mutableListOf<S2PointLoopSpan>()
        for (i in 0 until polygon.numLoops()) {
            val loop = polygon.loop(i)
            if (loop.isFull()) {
                spans.add(S2PointLoopSpan(emptyList()));  // Empty span.
            } else {
                spans.add(loop.verticesSpan())
            }
        }
        initSpans(spans)

        // S2Polygon and S2LaxPolygonShape holes are oriented oppositely, so we need
        // to reverse the orientation of any loops representing holes.
        for (i in 0 until polygon.numLoops()) {
            if (polygon.loop(i).isHole()) {
                val loopStartIdx = cumulativeVertices[i]
                vertices.reverse(loopStartIdx, loopStartIdx + numLoopVertices(i))
            }
        }
    }

    // Returns the number of loops.
    fun numLoops(): Int = numLoops

    // Returns the total number of vertices in all loops.
    fun numVertices(): Int = cumulativeVertices[numLoops]

    // Returns the number of vertices in the given loop.
    fun numLoopVertices(i: Int): Int {
        requireLT(i, numLoops)
        return cumulativeVertices[i + 1] - cumulativeVertices[i]
    }

    // Returns the vertex from loop "i" at index "j".
    // REQUIRES: 0 <= i < num_loops()
    // REQUIRES: 0 <= j < num_loop_vertices(i)
    fun loopVertex(i: Int, j: Int): S2Point {
        requireLT(i, numLoops)
        requireLT(j, numLoopVertices(i))
        return vertices[cumulativeVertices[i] + j]
    }

    // S2Shape interface:

    override val numEdges: Int
        get() = numVertices()

    override fun edge(edgeId: Int): Edge {
        requireLT(edgeId, numEdges)
        var e1 = edgeId + 1
        if (numLoops == 1) {
            if (e1 == numVertices()) { e1 = 0 }
        } else {
            // Find the index of the first vertex of the loop following this one.
            val kMaxLinearSearchLoops = 12;  // From benchmarks.

            var nextIdx = 1
            if (numLoops <= kMaxLinearSearchLoops) {
                while (cumulativeVertices[nextIdx] <= edgeId) ++nextIdx
            } else {
                nextIdx = cumulativeVertices.lowerBound(nextIdx, nextIdx + numLoops, e1)
            }
            // Wrap around to the first vertex of the loop if necessary.
            if (e1 == cumulativeVertices[nextIdx]) { e1 = cumulativeVertices[nextIdx-1] }
        }
        return Edge(vertices[edgeId], vertices[e1])
    }

    override val dimension: Int = 2

    override fun getReferencePoint(): ReferencePoint = S2Shape.getReferencePoint(this)

    override val numChains: Int
        get() = numLoops()

    override fun chain(chainId: Int): Chain {
        requireLT(chainId, numLoops)
        return if (numLoops == 1) {
            Chain(0, numVertices())
        } else {
            val start = cumulativeVertices[chainId]
            Chain(start, cumulativeVertices[chainId + 1] - start)
        }
    }

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        requireLT(chainId, numLoops)
        requireLT(offset, numLoopVertices(chainId))
        val n = numLoopVertices(chainId)
        val k = if(offset + 1 == n) 0 else offset + 1
        return if (numLoops() == 1) {
            Edge(vertices[offset], vertices[k])
        } else {
            val base = cumulativeVertices[chainId]
            Edge(vertices[base + offset], vertices[base + k])
        }
    }

    override fun chainPosition(edgeId: Int): ChainPosition {
        requireLT(edgeId, numEdges)
        val kMaxLinearSearchLoops = 12;  // From benchmarks.
        if (numLoops == 1) {
            return ChainPosition(0, edgeId);
        } else {
            // Find the index of the first vertex of the loop following this one.
            var nextIdx = 1
            if (numLoops <= kMaxLinearSearchLoops) {
                while (cumulativeVertices[nextIdx] <= edgeId) ++nextIdx
            } else {
                nextIdx =  cumulativeVertices.lowerBound(nextIdx, nextIdx + numLoops, edgeId + 1)
            }
            return ChainPosition(nextIdx - 1, edgeId - cumulativeVertices[nextIdx-1])
        }
    }

    override val typeTag: TypeTag = TypeTags.kLaxPolygonTypeTag

    private fun initSpans(loops: List<S2PointLoopSpan>) {
        numLoops = loops.size
        if (numLoops == 0) {
            vertices = emptyArray()
            cumulativeVertices = intArrayOf(0)
        } else if (numLoops == 1) {
            vertices = Array(loops[0].size) { i -> loops[0][i] }
            cumulativeVertices = intArrayOf(0, vertices.size)
        } else {
            cumulativeVertices = IntArray(numLoops + 1)
            var numVertices = 0
            for (i in 0 until numLoops) {
                cumulativeVertices[i] = numVertices
                numVertices += loops[i].size
            }
            cumulativeVertices[numLoops] = numVertices
            vertices = loops.asSequence()
                    .map { loop -> loop.points.toMutableList() }
                    .reduce { acc, list -> acc.addAll(list); acc }
                    .toTypedArray()
        }
    }


}

/*
// Exactly like S2LaxPolygonShape, except that the vertices are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
class EncodedS2LaxPolygonShape : public S2Shape {
 public:
  // Constructs an uninitialized object; requires Init() to be called.
  EncodedS2LaxPolygonShape() {}

  // Initializes an EncodedS2LaxPolygonShape.
  //
  // REQUIRES: The Decoder data buffer must outlive this object.
  bool Init(Decoder* decoder);

  int num_loops() const { return num_loops_; }
  int num_vertices() const;
  int num_loop_vertices(int i) const;
  S2Point loop_vertex(int i, int j) const;

  // S2Shape interface:
  int num_edges() const final { return num_vertices(); }
  Edge edge(int e) const final;
  int dimension() const final { return 2; }
  ReferencePoint GetReferencePoint() const final;
  int num_chains() const final { return num_loops(); }
  Chain chain(int i) const final;
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final;

 private:
  int32 num_loops_;
  s2coding::EncodedS2PointVector vertices_;
  s2coding::EncodedUintVector<uint32> cumulative_vertices_;
};
*/
