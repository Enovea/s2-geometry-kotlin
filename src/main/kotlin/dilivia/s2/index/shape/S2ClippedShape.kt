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

/**
 * S2ClippedShape represents the part of a shape that intersects an S2Cell. It consists of the set of edge ids that
 * intersect that cell, and a boolean indicating whether the center of the cell is inside the shape (for shapes that
 * have an interior).
 *
 * Note that the edges themselves are not clipped; we always use the original edges for intersection tests so that the
 * results will be the same as the original shape.
 *
 * @property shapeId The shape id of the clipped shape.
 * @property edges The ids of the edges that intersects the cell.
 * @property containsCenter Indicates if the center of the cell is inside the shape (for shape that have an interior).
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
data class S2ClippedShape(val shapeId: Int, val edges: List<Int>, val containsCenter: Boolean) {

    /** The number of edges that intersect the S2CellId. */
    val numEdges: Int
        get() = edges.size

    /**
     * Gets the edge id of the given edge in this clipped shape. Edges are sorted in increasing order of edge id.
     *
     * REQUIRES: 0 <= i < num_edges()
     * @return the id of the edge at index i.
     */
    fun edge(i: Int): Int = edges[i]

    /**
     * Checks if an edge is contains in the cell.
     *
     * @param id The identifier of the edge to check.
     * @return true if the clipped shape contains the given edge id.
     */
    fun containsEdge(id: Int): Boolean {
        // Linear search is fast because the number of edges per shape is typically
        // very small (less than 10).
        for (e in 0 until numEdges) {
            if (edge(e) == id) return true
        }
        return false
    }

}
