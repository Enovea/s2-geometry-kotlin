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

import dilivia.ComparisonChain
import dilivia.s2.index.Distance

/**
 * Each "Result" object represents a closest edge. Note the following special cases:
 *
 *  - (shapeId >= 0) && (edgeId < 0) represents the interior of a shape.
 *    Such results may be returned when options.includeInteriors is true.
 *    Such results can be identified using the isInterior() method.
 *
 *  - (shapeId < 0) && (edgeId < 0) is returned by `findClosestEdge`  to indicate that no edge satisfies the given query
 *    options.  Such results can be identified using isEmpty() method.
 *
 * This class is a port of the S2ClosestEdgeQueryResult class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T The distance type.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
data class S2ClosestEdgeQueryResult<T : Distance<T>>(
    /** The distance from the target to this edge. */
    val distance: T,
    /** Identifies an indexed shape. */
    val shapeId: Int = -1,
    /** Identifies an edge within the shape. */
    val edgeId: Int = -1
) : Comparable<S2ClosestEdgeQueryResult<T>> {

    /**
     * Indicates if the result objec represents the interior of a shape.
     * (Such results may be returned when options.includeInteriors is true.)
     *
     * @return true if this Result object represents the interior of a shape.
     */
    fun isInterior(): Boolean = shapeId >= 0 && edgeId < 0

    /**
     * Indicates if the result is empty.
     * (This result is only returned in one special case, namely when findClosestEdge() does not find any suitable edges.
     * It is never returned by methods that return a list of results.)
     *
     * @return true if this Result object indicates that no edge satisfies the given query options.
     */
    fun isEmpty(): Boolean = shapeId < 0

    // Compares edges first by distance, then by (shape_id, edge_id).
    override fun compareTo(other: S2ClosestEdgeQueryResult<T>): Int = ComparisonChain.start()
        .compare(distance, other.distance)
        .compare(shapeId, other.shapeId)
        .compare(edgeId, other.edgeId)
        .result()

}
