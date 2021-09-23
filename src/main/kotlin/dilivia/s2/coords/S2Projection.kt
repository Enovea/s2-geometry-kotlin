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
package dilivia.s2.coords

/**
 * Defines a projection of sphere cells to the cube-space.
 *
 * This interface defines various constants that describe the shapes and sizes of S2Cells.
 * They are useful for deciding which cell level to use in order to satisfy a given condition (e.g. that cell vertices
 * must be no further than "x" apart).  All of the raw constants are differential quantities; you can use the
 * getValue(level) method to compute the corresponding length or area on the unit sphere for cells at a given level.
 * The minimum and maximum bounds are valid for cells at all levels, but they may be somewhat conservative for very
 * large cells (e.g. face cells).
 *
 * Three different projections from cell-space (s,t) to cube-space (u,v) are implemented: linear, quadratic, and
 * tangent. They have the following tradeoffs:
 *
 *  Linear - This is the fastest transformation, but also produces the least
 *    uniform cell sizes.  Cell areas vary by a factor of about 5.2, with the
 *    largest cells at the center of each face and the smallest cells in
 *    the corners.
 *
 *  Tangent - Transforming the coordinates via atan() makes the cell sizes
 *    more uniform.  The areas vary by a maximum ratio of 1.4 as opposed to a
 *    maximum ratio of 5.2.  However, each call to atan() is about as expensive
 *    as all of the other calculations combined when converting from points to
 *    cell ids, i.e. it reduces performance by a factor of 3.
 *
 *  Quadratic - This is an approximation of the tangent projection that
 *    is much faster and produces cells that are almost as uniform in size.
 *    It is about 3 times faster than the tangent projection for converting
 *    cell ids to points or vice versa.  Cell areas vary by a maximum ratio of
 *    about 2.1.
 *
 * Here is a table comparing the cell uniformity using each projection.  "Area ratio" is the maximum ratio over all
 * subdivision levels of the largest cell area to the smallest cell area at that level, "edge ratio" is the maximum
 * ratio of the longest edge of any cell to the shortest edge of any cell at the same level, and "diag ratio" is the
 * ratio of the longest diagonal of any cell to the shortest diagonal of any cell at the same level. "ToPoint" and
 * "FromPoint" are the times in microseconds required to convert cell ids to and from points (unit vectors)
 * respectively.  "ToPointRaw" is the time to convert to a non-unit-length vector, which is all that is needed for
 * some purposes.
 *
 *  Area    Edge    Diag   ToPointRaw  ToPoint  FromPoint
 *  Ratio   Ratio   Ratio             (microseconds)
 *  -------------------------------------------------------------------
 *  Linear:      5.200   2.117   2.959      0.020     0.087     0.085
 *  Tangent:     1.414   1.414   1.704      0.237     0.299     0.258
 *  Quadratic:   2.082   1.802   1.932      0.033     0.096     0.108
 *
 * The worst-case cell aspect ratios are about the same with all three projections. The maximum ratio of the longest
 * edge to the shortest edge within the same cell is about 1.4 and the maximum ratio of the diagonals within the same
 * cell is about 1.7.
 *
 * This data was produced using s2cell_test and s2cell_id_test.
 *
 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
interface S2Projection {

    /**
     * Convert an s- or t-value to the corresponding u- or v-value. This is a non-linear transformation from [-1,1] to
     * [-1,1] that attempts to make the cell sizes more uniform.
     *
     * @param s s- or t-value (cell-space).
     * @return u- or v-value
     */
    fun stToUv(s: Double): Double

    /**
     * The inverse of the stToUV transformation.  Note that it is not always true that UVtoST(STtoUV(x)) == x due to
     * numerical errors.
     *
     * @param u u- or v-value (cube-space).
     * @return s- or t-value
    */
    fun uvToSt(u: Double): Double

    /**
     * Cell minimum angle span.
     *
     * Each cell is bounded by four planes passing through its four edges and the center of the sphere. These metrics
     * relate to the angle between each pair of opposite bounding planes, or equivalently, between the planes
     * corresponding to two different s-values or two different t-values. For example, the maximum angle between
     * opposite bounding planes for a cell at level k is kMaxAngleSpan.getValue(k), and the average angle span for all
     * cells at level k is approximately kAvgAngleSpan.getValue(k).
     */
    val kMinAngleSpan: LengthMetric

    /**
     * Cell maximum angle span.
     * @see kMinAngleSpan
     */
    val kMaxAngleSpan: LengthMetric

    /**
     * Cell average angle span.
     * @see kMinAngleSpan
     */
    val kAvgAngleSpan: LengthMetric

    /**
     * Cell min width
     *
     * The width of geometric figure is defined as the distance between two parallel bounding lines in a given
     * direction.  For cells, the minimum width is always attained between two opposite edges, and the maximum width is
     * attained between two opposite vertices.  However, for our purposes we redefine the width of a cell as the
     * perpendicular distance between a pair of opposite edges.  A cell therefore has two widths, one in each direction.
     * The minimum width according to this definition agrees with the classic geometric one, but the maximum width is
     * different.  (The maximum geometric width corresponds to kMaxDiag defined below.)
     *
     * For a cell at level k, the distance between opposite edges is at least kMinWidth.getValue(k) and at most
     * kMaxWidth.getValue(k). The average width in both directions for all cells at level k is approximately
     * kAvgWidth.getValue(k).
     *
     * The width is useful for bounding the minimum or maximum distance from a point on one edge of a cell to the
     * closest point on the opposite edge. For example, this is useful when "growing" regions by a fixed distance.
     *
     * Note that because S2Cells are not usually rectangles, the minimum width of a cell is generally smaller than its
     * minimum edge length.  (The interior angles of an S2Cell range from 60 to 120 degrees.)
     */
    val kMinWidth: LengthMetric

    /**
     * Cell max width
     *
     * @see kMinWidth
     */
    val kMaxWidth: LengthMetric

    /**
     * Cell average width
     *
     * @see kMinWidth
     */
    val kAvgWidth: LengthMetric

    /**
     * Cell min edge.
     *
     * The minimum edge length of any cell at level k is at least kMinEdge.getValue(k), and the maximum is at most
     * kMaxEdge.getValue(k). The average edge length is approximately kAvgEdge.getValue(k).
     *
     * The edge length metrics can also be used to bound the minimum, maximum, or average distance from the center of
     * one cell to the center of one of its edge neighbors.  In particular, it can be used to bound the distance
     * between adjacent cell centers along the space-filling Hilbert curve for cells at any given level.
     */
    val kMinEdge: LengthMetric

    /**
     * Cell max edge.
     * @see kMinEdge
     */
    val kMaxEdge: LengthMetric

    /**
     * Cell average edge.
     * @see kMinEdge
     */
    val kAvgEdge: LengthMetric

    /**
     * Cell minimum diagonal.
     *
     * The minimum diagonal length of any cell at level k is at least kMinDiag.getValue(k), and the maximum is at most
     * kMaxDiag.getValue(k). The average diagonal length is approximately kAvgDiag.getValue(k).
     *
     * The maximum diagonal also happens to be the maximum diameter of any cell, and also the maximum geometric width
     * (see the discussion above).  So for example, the distance from an arbitrary point to the closest cell center
     * at a given level is at most half the maximum diagonal length.
     */
    val kMinDiag: LengthMetric

    /**
     * Cell maximum diagonal
     * @see kMinDiag
     */
    val kMaxDiag: LengthMetric

    /**
     * Cell average diagonal
     * @see kMinDiag
     */
    val kAvgDiag: LengthMetric

    /**
     * Minimum cell area.
     *
     * The minimum area of any cell at level k is at least kMinArea.getValue(k), and the maximum is at most
     * kMaxArea.getValue(k).  The average area of all cells at level k is exactly kAvgArea.getValue(k).
     */
    val kMinArea: AreaMetric

    /**
     * Maximum cell area.
     * @see kMinArea
     */
    val kMaxArea: AreaMetric

    /**
     * Maximum average area.
     * @see kMinArea
     */
    val kAvgArea: AreaMetric

    /**
     * This is the maximum edge aspect ratio over all cells at any level, where the edge aspect ratio of a cell is
     * defined as the ratio of its longest edge length to its shortest edge length.
     */
    val kMaxEdgeAspect: Double

    /**
     * This is the maximum diagonal aspect ratio over all cells at any level, where the diagonal aspect ratio of a cell
     * is defined as the ratio of its longest diagonal length to its shortest diagonal length.
     */
    val kMaxDiagAspect: Double

}
