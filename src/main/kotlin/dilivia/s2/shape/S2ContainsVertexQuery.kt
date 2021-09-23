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
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Predicates
import org.apache.commons.math3.util.FastMath.abs
import java.util.*

// This class determines whether a polygon contains one of its vertices given
// the edges incident to that vertex.  The result is +1 if the vertex is
// contained, -1 if it is not contained, and 0 if the incident edges consist
// of matched sibling pairs (in which case the result cannot be determined
// locally).
//
// Point containment is defined according to the "semi-open" boundary model
// (see S2VertexModel), which means that if several polygons tile the region
// around a vertex, then exactly one of those polygons contains that vertex.
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
// "target" is the vertex whose containment will be determined.
class S2ContainsVertexQuery(val target: S2Point) {

    private val edge_map = TreeMap<S2Point, Int>()

    // Indicates that the polygon has an edge between "target" and "v" in the
    // given direction (+1 = outgoing, -1 = incoming, 0 = degenerate).
    fun addEdge(v: S2Point, direction: Int) {
        edge_map.compute(v) { _, d -> (d ?: 0) + direction }
    }

    // Returns +1 if the vertex is contained, -1 if it is not contained, and 0
    // if the incident edges consisted of matched sibling pairs.
    fun containsSign(): Int {
        // Find the unmatched edge that is immediately clockwise from S2::Ortho(P).
        val referenceDir = S2PointUtil.ortho(target)
        var best = Pair(referenceDir, 0)
        for (e in edge_map) {
            checkLE(abs(e.value), 1)
            if (e.value == 0) continue;  // This is a "matched" edge.
            if (S2Predicates.orderedCCW(referenceDir, best.first, e.key, target)) {
                best = e.toPair()
            }
        }
        return best.second;
    }

}

