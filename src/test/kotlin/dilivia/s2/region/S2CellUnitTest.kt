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
package dilivia.s2.region

import dilivia.math.DoubleType
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.vectors.R2Point
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random.oneIn
import dilivia.s2.S2Random.randomCellId
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomInt
import dilivia.s2.S2Random.randomPoint
import dilivia.s2.S2Random.samplePoint
import dilivia.s2.coords.S2CellMetric
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2Coords.kMaxCellLevel
import dilivia.s2.coords.S2Coords.kSwapMask
import dilivia.s2.edge.S2EdgeCrossings
import dilivia.s2.edge.S2EdgeDistances
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.abs
import kotlin.math.ln
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow
import kotlin.math.sqrt

class S2CellUnitTest {

    private val logger = KotlinLogging.logger { }

    @Test
    fun faces() {
        val edgeCounts = mutableMapOf<S2Point, Int>()
        val vertexCounts = mutableMapOf<S2Point, Int>()
        for (face in 0..5) {
            val id = S2CellId.fromFace(face)
            val cell = S2Cell(id)
            assertThat(cell.id()).isEqualTo(id)
            assertThat(cell.face()).isEqualTo(face)
            assertThat(cell.level()).isEqualTo(0)
            // Top-level faces have alternating orientations to get RHS coordinates.
            assertThat(cell.orientation()).isEqualTo(face and kSwapMask)
            assertThat(cell.isLeaf()).isFalse()

            for (k in 0..3) {
                edgeCounts.compute(cell.getEdgeRaw(k)) { _, i -> (i ?: 0) + 1 }
                vertexCounts.compute(cell.getVertexRaw(k)) { _, i -> (i ?: 0) + 1 }
                assertThat(cell.getVertexRaw(k).dotProd(cell.getEdgeRaw(k))).isCloseTo(0.0, Offset.offset(1e-15))
                assertThat(cell.getVertexRaw(k + 1).dotProd(cell.getEdgeRaw(k))).isCloseTo(0.0, Offset.offset(1e-15))
                assertThat(
                    cell.getVertexRaw(k).crossProd(cell.getVertexRaw(k + 1)).normalize().dotProd(cell.getEdge(k))
                ).isCloseTo(1.0, Offset.offset(1e-15))
            }
        }
        // Check that edges have multiplicity 2 and vertices have multiplicity 3.
        for (p in edgeCounts) {
            assertThat(p.value).isEqualTo(2)
        }
        for (p in vertexCounts) {
            assertThat(p.value).isEqualTo(3)
        }
    }

    data class LevelStats(
        var count: Double = 0.0,
        var min_area: Double = 100.0,
        var max_area: Double = 0.0,
        var avg_area: Double = 0.0,
        var min_width: Double = 100.0,
        var max_width: Double = 0.0,
        var avg_width: Double = 0.0,
        var min_edge: Double = 100.0,
        var max_edge: Double = 0.0,
        var avg_edge: Double = 0.0,
        var max_edge_aspect: Double = 0.0,
        var min_diag: Double = 100.0,
        var max_diag: Double = 0.0,
        var avg_diag: Double = 0.0,
        var max_diag_aspect: Double = 0.0,
        var min_angle_span: Double = 100.0,
        var max_angle_span: Double = 0.0,
        var avg_angle_span: Double = 0.0,
        var min_approx_ratio: Double = 100.0,
        var max_approx_ratio: Double = 0.0
    )

    private fun gatherStats(cell: S2Cell, level_stats: MutableList<LevelStats>) {
        val exactArea = cell.exactArea()
        val approxArea = cell.approxArea()
        var minEdge = 100.0
        var maxEdge = 0.0
        var avgEdge = 0.0
        var minDiag = 100.0
        var maxDiag = 0.0
        var minWidth = 100.0
        var maxWidth = 0.0
        var minAngleSpan = 100.0
        var maxAngleSpan = 0.0
        for (i in 0..3) {
            val edge = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 1))
            minEdge = min(edge, minEdge)
            maxEdge = max(edge, maxEdge)
            avgEdge += 0.25 * edge
            val mid = cell.getVertexRaw(i) + cell.getVertexRaw(i + 1)
            val width = M_PI_2 - mid.angle(cell.getEdgeRaw(i + 2))
            minWidth = min(width, minWidth)
            maxWidth = max(width, maxWidth)
            if (i < 2) {
                val diag = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 2))
                minDiag = min(diag, minDiag)
                maxDiag = max(diag, maxDiag)
                val angleSpan = cell.getEdgeRaw(i).angle(-cell.getEdgeRaw(i + 2))
                minAngleSpan = min(angleSpan, minAngleSpan)
                maxAngleSpan = max(angleSpan, maxAngleSpan)
            }
        }

        val s = level_stats[cell.level()]
        s.count += 1
        s.min_area = min(exactArea, s.min_area)
        s.max_area = max(exactArea, s.max_area)
        s.avg_area += exactArea
        s.min_width = min(minWidth, s.min_width)
        s.max_width = max(maxWidth, s.max_width)
        s.avg_width += 0.5 * (minWidth + maxWidth)
        s.min_edge = min(minEdge, s.min_edge)
        s.max_edge = max(maxEdge, s.max_edge)
        s.avg_edge += avgEdge
        s.max_edge_aspect = max(maxEdge / minEdge, s.max_edge_aspect)
        s.min_diag = min(minDiag, s.min_diag)
        s.max_diag = max(maxDiag, s.max_diag)
        s.avg_diag += 0.5 * (minDiag + maxDiag)
        s.max_diag_aspect = max(maxDiag / minDiag, s.max_diag_aspect)
        s.min_angle_span = min(minAngleSpan, s.min_angle_span)
        s.max_angle_span = max(maxAngleSpan, s.max_angle_span)
        s.avg_angle_span += 0.5 * (minAngleSpan + maxAngleSpan)
        val apprapproxratioxRatio = approxArea / exactArea
        s.min_approx_ratio = min(apprapproxratioxRatio, s.min_approx_ratio)
        s.max_approx_ratio = max(apprapproxratioxRatio, s.max_approx_ratio)
    }

    private fun testSubdivide(cell: S2Cell, level_stats: MutableList<LevelStats>) {
        logger.trace { "Test subdivide: $cell" }
        gatherStats(cell, level_stats)
        if (cell.isLeaf()) return

        val children = Array<S2Cell>(4) { S2Cell() }
        assertThat(cell.subdivide(children)).isTrue()
        var childId = cell.id().childBegin()
        var exactArea = 0.0
        var approxArea = 0.0
        var averageArea = 0.0
        for (i in 0..3) {
            logger.trace { "Check children $i: ${children[i]}" }
            exactArea += children[i].exactArea()
            approxArea += children[i].approxArea()
            averageArea += children[i].averageArea()

            // Check that the child geometry is consistent with its cell ID.
            assertThat(children[i].id()).isEqualTo(childId)
            assertThat(S2PointUtil.approxEquals(children[i].getCenter(), childId.toPoint())).isTrue()
            val direct = S2Cell(childId)
            assertThat(children[i].face()).isEqualTo(direct.face())
            assertThat(children[i].level()).isEqualTo(direct.level())
            assertThat(children[i].orientation()).isEqualTo(direct.orientation())
            assertThat(children[i].getCenterRaw()).isEqualTo(direct.getCenterRaw())
            for (k in 0..3) {
                assertThat(children[i].getVertexRaw(k)).isEqualTo(direct.getVertexRaw(k))
                assertThat(children[i].getEdgeRaw(k)).isEqualTo(direct.getEdgeRaw(k))
            }

            // Test Contains() and MayIntersect().
            assertThat(cell.contains(children[i])).isTrue()
            assertThat(cell.mayIntersect(children[i])).isTrue()
            assertThat(children[i].contains(cell)).isFalse()
            assertThat(cell.contains(children[i].getCenterRaw())).isTrue()
            for (j in 0..3) {
                assertThat(cell.contains(children[i].getVertexRaw(j))).isTrue()
                if (j != i) {
                    assertThat(children[i].contains(children[j].getCenterRaw())).isFalse()
                    assertThat(children[i].mayIntersect(children[j])).isFalse()
                }
            }

            // Test GetCapBound and GetRectBound.
            val parentCap = cell.capBound
            val parentRect = cell.rectBound
            if (cell.contains(S2Point(0, 0, 1)) || cell.contains(S2Point(0, 0, -1))) {
                assertThat(parentRect.lng.isFull).isTrue()
            }
            val child_cap = children[i].capBound
            val child_rect = children[i].rectBound
            assertThat(child_cap.contains(children[i].getCenter())).isTrue()
            assertThat(child_rect.contains(children[i].getCenterRaw())).isTrue()
            assertThat(parentCap.contains(children[i].getCenter())).isTrue()
            assertThat(parentRect.contains(children[i].getCenterRaw())).isTrue()
            for (j in 0..3) {
                val vertex = children[i].getVertex(j)
                val vertexRaw = children[i].getVertexRaw(j)
                logger.trace { "Check vertex $j = $vertex, raw = $vertexRaw" }
                assertThat(child_cap.contains(vertex)).isTrue()
                assertThat(child_rect.contains(vertex)).isTrue()
                assertThat(child_rect.contains(vertexRaw)).isTrue()
                assertThat(parentCap.contains(vertex)).isTrue()
                assertThat(parentRect.contains(vertex)).isTrue()
                assertThat(parentRect.contains(vertexRaw)).isTrue()
                if (j != i) {
                    // The bounding caps and rectangles should be tight enough so that
                    // they exclude at least two vertices of each adjacent cell.
                    var cap_count = 0
                    var rect_count = 0
                    for (k in 0..3) {
                        if (child_cap.contains(children[j].getVertex(k))) {
                            logger.trace { "Child cap $child_cap contains vertex $k = ${children[j].getVertex(k)} of child $j = ${children[j]}" }
                            ++cap_count
                        }
                        if (child_rect.contains(children[j].getVertexRaw(k))) {
                            ++rect_count
                        }
                    }
                    assertThat(cap_count <= 2).isTrue()
                    if (child_rect.latLo().radians > -M_PI_2 && child_rect.latHi().radians < M_PI_2) {
                        // Bounding rectangles may be too large at the poles because the
                        // pole itself has an arbitrary fixed longitude.
                        assertThat(rect_count <= 2).isTrue()
                    }
                }
            }

            // Check all children for the first few levels, and then sample randomly.
            // We also always subdivide the cells containing a few chosen points so
            // that we have a better chance of sampling the minimum and maximum metric
            // values.  kMaxSizeUV is the absolute value of the u- and v-coordinate
            // where the cell size at a given level is maximal.
            val kMaxSizeUV = 0.3964182625366691
            val special_uv = arrayOf(
                R2Point(DoubleType.epsilon, DoubleType.epsilon),  // Face center
                R2Point(DoubleType.epsilon, 1.0),       // Edge midpoint
                R2Point(1, 1),                // Face corner
                R2Point(kMaxSizeUV, kMaxSizeUV),    // Largest cell area
                R2Point(DoubleType.epsilon, kMaxSizeUV),   // Longest edge/diagonal
            )
            var force_subdivide = false
            for (uv in special_uv) {
                if (children[i].boundUV().contains(uv))
                    force_subdivide = true
            }
            if (force_subdivide || cell.level() < 6 || oneIn(5)) {
                testSubdivide(children[i], level_stats)
            }
            childId = childId.next()
        }

        // Check sum of child areas equals parent area.
        //
        // For ExactArea(), the best relative error we can expect is about 1e-6
        // because the precision of the unit vector coordinates is only about 1e-15
        // and the edge length of a leaf cell is about 1e-9.
        //
        // For ApproxArea(), the areas are accurate to within a few percent.
        //
        // For AverageArea(), the areas themselves are not very accurate, but
        // the average area of a parent is exactly 4 times the area of a child.

        assertThat(abs(ln(exactArea / cell.exactArea())) <= abs(ln(1 + 1e-6))).isTrue()
        assertThat(abs(ln(approxArea / cell.approxArea())) <= abs(ln(1.03))).isTrue()
        assertThat(abs(ln(averageArea / cell.averageArea())) <= abs(ln(1 + 1e-15))).isTrue()
    }

    private fun checkMinMaxAvg(
        dim: Int, label: String, level: Int, count: Double, abs_error: Double,
        min_value: Double, max_value: Double, avg_value: Double,
        min_metric: S2CellMetric, max_metric: S2CellMetric, avg_metric: S2CellMetric
    ) {

        // All metrics are minimums, maximums, or averages of differential
        // quantities, and therefore will not be exact for cells at any finite
        // level.  The differential minimum is always a lower bound, and the maximum
        // is always an upper bound, but these minimums and maximums may not be
        // achieved for two different reasons.  First, the cells at each level are
        // sampled and we may miss the most extreme examples.  Second, the actual
        // metric for a cell is obtained by integrating the differential quantity,
        // which is not constant across the cell.  Therefore cells at low levels
        // (bigger cells) have smaller variations.
        //
        // The "tolerance" below is an attempt to model both of these effects.
        // At low levels, error is dominated by the variation of differential
        // quantities across the cells, while at high levels error is dominated by
        // the effects of random sampling.
        var tolerance =
            (max_metric.getValue(level) - min_metric.getValue(level)) / sqrt(min(count, 0.5 * (1 shl level).toDouble()))
        if (tolerance == 0.0) tolerance = abs_error

        val minError = min_value - min_metric.getValue(level)
        val maxError = max_metric.getValue(level) - max_value
        val avgError = abs(avg_metric.getValue(level) - avg_value)
        logger.debug(
            String.format(
                "%-10s (%6.0f samples, tolerance %8.3g) - min %9.4g (%9.3g : %9.3g) max %9.4g (%9.3g : %9.3g), avg %9.4g (%9.3g : %9.3g)",
                label, count, tolerance,
                min_value, minError / min_value, minError / tolerance,
                max_value, maxError / max_value, maxError / tolerance,
                avg_value, avgError / avg_value, avgError / tolerance
            )
        )

        assertThat(min_metric.getValue(level) <= min_value + abs_error).isTrue()
        assertThat(min_metric.getValue(level) >= min_value - tolerance).isTrue()
        assertThat(max_metric.getValue(level) <= max_value + tolerance).isTrue()
        assertThat(max_metric.getValue(level) >= max_value - abs_error).isTrue()
        assertThat(avg_metric.getValue(level)).isCloseTo(avg_value, Offset.offset(10 * tolerance))
    }

    @Test
    fun subdivide() {
        val levelStats = (0..kMaxCellLevel).map { LevelStats() }.toMutableList()
        // Only test a sample of faces to reduce the runtime.
        testSubdivide(S2Cell.fromFace(0), levelStats)
        testSubdivide(S2Cell.fromFace(3), levelStats)
        testSubdivide(S2Cell.fromFace(5), levelStats)

        // This table is useful in evaluating the quality of the various S2
        // projections.
        //
        // The maximum edge *ratio* is the ratio of the longest edge of any cell to
        // the shortest edge of any cell at the same level (and similarly for the
        // maximum diagonal ratio).
        //
        // The maximum edge *aspect* is the maximum ratio of the longest edge of a
        // cell to the shortest edge of that same cell (and similarly for the
        // maximum diagonal aspect).
        println(
            """
                | Ratio:  (Max value for any cell) / (Min value for any cell)
                | Aspect: (Max value / min value) for any cell
                | Edge          Diag       Approx Area/    Avg Area/
                | Area     Length        Length       Exact Area    Exact Area
                | Level   Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max
                | --------------------------------------------------------------------
                | """.trimMargin()
        )
        for (i in 0..S2CellId.kMaxLevel) {
            val s = levelStats[i]
            if (s.count > 0) {
                s.avg_area /= s.count
                s.avg_width /= s.count
                s.avg_edge /= s.count
                s.avg_diag /= s.count
                s.avg_angle_span /= s.count
            }
            logger.debug(
                String.format(
                    "%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f",
                    i, s.max_area / s.min_area,
                    s.max_edge / s.min_edge, s.max_edge_aspect,
                    s.max_diag / s.min_diag, s.max_diag_aspect,
                    s.min_approx_ratio, s.max_approx_ratio,
                    S2Cell.averageArea(i) / s.max_area,
                    S2Cell.averageArea(i) / s.min_area
                )
            )
        }

        // Now check the validity of the S2 length and area metrics.
        for (i in 0..S2CellId.kMaxLevel) {
            val s = levelStats[i]
            if (s.count == 0.0) continue

            logger.debug(String.format("Level %2d - metric value (error/actual : error/tolerance)", i))

            // The various length calculations are only accurate to 1e-15 or so,
            // so we need to allow for this amount of discrepancy with the theoretical
            // minimums and maximums.  The area calculation is accurate to about 1e-15
            // times the cell width.
            checkMinMaxAvg(
                2,
                "area",
                i,
                s.count,
                1e-15 * s.min_width,
                s.min_area,
                s.max_area,
                s.avg_area,
                S2Coords.projection.kMinArea,
                S2Coords.projection.kMaxArea,
                S2Coords.projection.kAvgArea
            )
            checkMinMaxAvg(
                1,
                "width",
                i,
                s.count,
                1e-15,
                s.min_width,
                s.max_width,
                s.avg_width,
                S2Coords.projection.kMinWidth,
                S2Coords.projection.kMaxWidth,
                S2Coords.projection.kAvgWidth
            )
            checkMinMaxAvg(
                1,
                "edge",
                i,
                s.count,
                1e-15,
                s.min_edge,
                s.max_edge,
                s.avg_edge,
                S2Coords.projection.kMinEdge,
                S2Coords.projection.kMaxEdge,
                S2Coords.projection.kAvgEdge
            )
            checkMinMaxAvg(
                1,
                "diagonal",
                i,
                s.count,
                1e-15,
                s.min_diag,
                s.max_diag,
                s.avg_diag,
                S2Coords.projection.kMinDiag,
                S2Coords.projection.kMaxDiag,
                S2Coords.projection.kAvgDiag
            )
            checkMinMaxAvg(
                1,
                "angle span",
                i,
                s.count,
                1e-15,
                s.min_angle_span,
                s.max_angle_span,
                s.avg_angle_span,
                S2Coords.projection.kMinAngleSpan,
                S2Coords.projection.kMaxAngleSpan,
                S2Coords.projection.kAvgAngleSpan
            )

            // The aspect ratio calculations are ratios of lengths and are therefore
            // less accurate at higher subdivision levels.
            assertThat(s.max_edge_aspect <= S2Coords.projection.kMaxEdgeAspect + 1e-15 * (1 shl i)).isTrue()
            assertThat(s.max_diag_aspect <= S2Coords.projection.kMaxDiagAspect + 1e-15 * (1 shl i)).isTrue()
        }
    }

    @Test
    fun cellVsLoopRectBound() {
        // This test verifies that the S2Cell and S2Loop bounds contain each other
        // to within their maximum errors.
        //
        // The S2Cell and S2Loop calculations for the latitude of a vertex can differ
        // by up to 2 * DoubleType.epsilon, therefore the S2Cell bound should never exceed
        // the S2Loop bound by more than this (the reverse is not true, because the
        // S2Loop code sometimes thinks that the maximum occurs along an edge).
        // Similarly, the longitude bounds can differ by up to 4 * DoubleType.epsilon since
        // the S2Cell bound has an error of 2 * DoubleType.epsilon and then expands by this
        // amount, while the S2Loop bound does no expansion at all.

        // Possible additional S2Cell error compared to S2Loop error:
        val kCellError = S2LatLng.fromRadians(2 * DoubleType.epsilon, 4 * DoubleType.epsilon)
        // Possible additional S2Loop error compared to S2Cell error:
        val kLoopError = S2LatLngRectBounder.maxErrorForTests()

        for (iter in 0 until 1000) {
            val cell = S2Cell(randomCellId())
            val loop = S2Loop(cell)
            val cell_bound = cell.rectBound
            val loop_bound = loop.rectBound
            assertThat(loop_bound.expanded(kCellError).contains(cell_bound)).isTrue()
            assertThat(cell_bound.expanded(kLoopError).contains(loop_bound)).isTrue()
        }
    }

    @Test
    fun rectBoundIsLargeEnough() {
        // Construct many points that are nearly on an S2Cell edge, and verify that
        // whenever the cell contains a point P then its bound contains S2LatLng(P).
        var iter = 0
        while (iter < 1000 /* advanced in loop below */) {
            val cell = S2Cell(randomCellId())
            val i = randomInt(4)
            val v1 = cell.getVertex(i)
            val v2 = samplePoint(S2Cap.fromCenterAngle(cell.getVertex(i + 1), S1Angle.radians(1e-15)))
            val fraction = randomDouble()
            val p = S2EdgeDistances.interpolate(fraction, v1, v2)

            logger.debug {
                """cell = ${cell.id()}
                |i = $i
                |v1 = $v1
                |v2 = $v2
                |fraction = $fraction
                |p = $p 
                |contains = ${S2Loop(cell).contains(p)}
            """.trimMargin()
            }

            if (S2Loop(cell).contains(p)) {
                assertThat(cell.rectBound.contains(S2LatLng.fromPoint(p)))
                    .withFailMessage("cell ${cell.id()} contains point p $p but bound ${cell.rectBound} does not.")
                    .isTrue();
                ++iter;
            }
        }
    }

    @Test
    fun consistentWithS2CellIdFromPoint() {
        // Construct many points that are nearly on an S2Cell edge, and verify that
        // S2Cell(S2CellId(p)).contains(p) is always true.
        for (iter in 0 until 1000) {
            val cell = S2Cell(randomCellId())
            val i = randomInt(4)
            val v1 = cell.getVertex(i)
            val v2 = samplePoint(S2Cap.fromCenterAngle(cell.getVertex(i + 1), S1Angle.radians(1e-15)))
            val p = S2EdgeDistances.interpolate(randomDouble(), v1, v2)
            assertThat(S2Cell(S2CellId.fromPoint(p)).contains(p)).isTrue();
        }
    }

    @Test
    fun ambiguousContainsPoint() {
        // This tests a case where S2CellId returns the "wrong" cell for a point
        // that is very close to the cell edge. (ConsistentWithS2CellIdFromPoint
        // generates more examples like this.)
        //
        // The S2Point below should have x = 0, but conversion from latlng to
        // (x,y,z) gives x = 6.1e-17.  When xyz is converted to uv, this gives u =
        // -6.1e-17.  However when converting to st, which is centered at 0.5 rather
        // than 0, the low precision bits of u are lost and we wind up with s = 0.5.
        // S2CellId(const S2Point&) then chooses an arbitrary neighboring cell.
        //
        // This tests that S2Cell::Contains() expands the cell bounds sufficiently
        // so that the returned cell is still considered to contain "p".
        val p = S2LatLng.fromDegrees(-2, 90).toPoint()
        val cell_id = S2CellId.fromPoint(p).parent(1)
        val cell = S2Cell(cell_id)
        assertThat(cell.contains(p)).isTrue()
    }

    private fun getDistanceToPointBruteForce(cell: S2Cell, target: S2Point): S1ChordAngle {
        val minDistance = S1ChordAngle.infinity().clone()
        for (i in 0..3) {
            S2EdgeDistances.updateMinDistance(target, cell.getVertex(i), cell.getVertex(i + 1), minDistance)
        }
        return minDistance
    }

    private fun getMaxDistanceToPointBruteForce(cell: S2Cell, target: S2Point): S1ChordAngle {
        if (cell.contains(-target)) {
            return S1ChordAngle.straight()
        }
        val maxDistance = S1ChordAngle.negative().clone()
        for (i in 0..3) {
            S2EdgeDistances.updateMaxDistance(target, cell.getVertex(i), cell.getVertex(i + 1), maxDistance)
        }
        return maxDistance
    }

    @Test
    fun getDistanceToPoint() {
        for (iter in 0 until 1000) {
            val cell = S2Cell(randomCellId())
            val target = randomPoint()
            val expected_to_boundary = getDistanceToPointBruteForce(cell, target).toAngle()
            val expected_to_interior = if (cell.contains(target)) S1Angle.zero() else expected_to_boundary
            val expected_max = getMaxDistanceToPointBruteForce(cell, target).toAngle()
            val actual_to_boundary = cell.getBoundaryDistance(target).toAngle()
            val actual_to_interior = cell.getDistance(target).toAngle()
            val actual_max = cell.getMaxDistance(target).toAngle()
            // The error has a peak near Pi/2 for edge distance, and another peak near
            // Pi for vertex distance.
            assertThat(actual_to_boundary.radians).isCloseTo(expected_to_boundary.radians, Offset.offset(1e-12))
            assertThat(actual_to_interior.radians).isCloseTo(expected_to_interior.radians, Offset.offset(1e-12))
            assertThat(actual_max.radians).isCloseTo(expected_max.radians, Offset.offset(1e-12))
            if (expected_to_boundary.radians <= M_PI / 3) {
                assertThat(actual_to_boundary.radians).isCloseTo(expected_to_boundary.radians, Offset.offset(1e-15))
                assertThat(actual_to_interior.radians).isCloseTo(expected_to_interior.radians, Offset.offset(1e-15))
            }
            if (expected_max.radians <= M_PI / 3) {
                assertThat(actual_max.radians).isCloseTo(expected_max.radians, Offset.offset(1e-15))
            }
        }
    }

    fun chooseEdgeNearCell(cell: S2Cell): Pair<S2Point, S2Point> {
        val cap = cell.capBound
        var a = if (oneIn(5)) {
            // Choose a point anywhere on the sphere.
            randomPoint()
        } else {
            // Choose a point inside or somewhere near the cell.
            samplePoint(S2Cap.fromCenterAngle(cap.center, 1.5 * cap.radius()))
        }
        // Now choose a maximum edge length ranging from very short to very long
        // relative to the cell size, and choose the other endpoint.
        val max_length = min(100 * 1e-4.pow(randomDouble()) * cap.radius().radians, M_PI_2)
        var b = samplePoint(S2Cap.fromCenterAngle(a, S1Angle.radians(max_length)))

        if (oneIn(20)) {
            // Occasionally replace edge with antipodal edge.
            a = -a
            b = -b
        }

        return a to b
    }

    fun getDistanceToEdgeBruteForce(cell: S2Cell, a: S2Point, b: S2Point): S1ChordAngle {
        if (cell.contains(a) || cell.contains(b)) {
            return S1ChordAngle.zero()
        }

        val minDist = S1ChordAngle.infinity().clone()
        for (i in 0..3) {
            val v0 = cell.getVertex(i)
            val v1 = cell.getVertex(i + 1)
            // If the edge crosses through the cell, max distance is 0.
            if (S2EdgeCrossings.crossingSign(a, b, v0, v1) >= 0) {
                return S1ChordAngle.zero()
            }
            S2EdgeDistances.updateMinDistance(a, v0, v1, minDist)
            S2EdgeDistances.updateMinDistance(b, v0, v1, minDist)
            S2EdgeDistances.updateMinDistance(v0, a, b, minDist)
        }
        return minDist
    }

    fun getMaxDistanceToEdgeBruteForce(cell: S2Cell, a: S2Point, b: S2Point): S1ChordAngle {
        // If any antipodal endpoint is within the cell, the max distance is Pi.
        if (cell.contains(-a) || cell.contains(-b)) {
            return S1ChordAngle.straight()
        }

        val maxDist = S1ChordAngle.negative().clone()
        for (i in 0..3) {
            val v0 = cell.getVertex(i)
            val v1 = cell.getVertex(i + 1)
            // If the antipodal edge crosses through the cell, max distance is Pi.
            if (S2EdgeCrossings.crossingSign(-a, -b, v0, v1) >= 0) {
                return S1ChordAngle.straight()
            }
            S2EdgeDistances.updateMaxDistance(a, v0, v1, maxDist)
            S2EdgeDistances.updateMaxDistance(b, v0, v1, maxDist)
            S2EdgeDistances.updateMaxDistance(v0, a, b, maxDist)
        }
        return maxDist
    }

    @Test
    fun getDistanceToEdge() {
        for (iter in 0 until 1000) {
            //S2Random.reset(iter + 1)

            val cell = S2Cell(randomCellId())
            val (a, b) = chooseEdgeNearCell(cell)
            val expectedMin = getDistanceToEdgeBruteForce(cell, a, b).toAngle()
            val expectedMax = getMaxDistanceToEdgeBruteForce(cell, a, b).toAngle()
            val actualMin = cell.getDistance(a, b).toAngle()
            val actualMax = cell.getMaxDistance(a, b).toAngle()

            logger.trace { """
                | cell = ${cell.id()}
                | edge: a = $a ; b = $b
                | expectedMin = $expectedMin
                | expectedMax = $expectedMax
                | actualMin = $actualMin
                | actualMax = $actualMax
            """.trimMargin() }
            // The error has a peak near Pi/2 for edge distance, and another peak near
            // Pi for vertex distance.
            if (expectedMin.radians > M_PI / 2) {
                // Max error for S1ChordAngle as it approaches Pi is about 2e-8.
                assertThat(actualMin.radians)
                    .withFailMessage("Error = ${abs(expectedMin.radians - actualMin.radians)} > 3e-8")
                    .isCloseTo(expectedMin.radians, Offset.offset(3e-8))
            } else if (expectedMin.radians <= M_PI / 3) {
                assertThat(actualMin.radians).isCloseTo(expectedMin.radians, Offset.offset(1e-15))
            } else {
                assertThat(actualMin.radians).isCloseTo(expectedMin.radians, Offset.offset(1e-12))
            }

            assertThat(actualMax.radians).isCloseTo(expectedMax.radians, Offset.offset(1e-12))
            if (expectedMax.radians <= M_PI / 3) {
                assertThat(actualMax.radians).isCloseTo(expectedMax.radians, Offset.offset(1e-15))
            }
        }
    }

    @Test
    fun getMaxDistanceToEdge() {
        // Test an edge for which its antipode crosses the cell. Validates both the
        // standard and brute force implementations for this case.
        val cell = S2Cell.fromFacePosLevel(0, 0, 20)
        val a = -S2EdgeDistances.interpolate(2.0, cell.getCenter(), cell.getVertex(0))
        val b = -S2EdgeDistances.interpolate(2.0, cell.getCenter(), cell.getVertex(2))

        val actual = cell.getMaxDistance(a, b)
        val expected = getMaxDistanceToEdgeBruteForce(cell, a, b)

        assertThat(S1ChordAngle.straight().radians()).isCloseTo(expected.radians(), Offset.offset(1e-15))
        assertThat(S1ChordAngle.straight().radians()).isCloseTo(actual.radians(), Offset.offset(1e-15))
    }

    @Test
    fun getMaxDistanceToCellAntipodal() {
        val p = S2LatLng(0, 0).toPoint()
        val cell = S2Cell(p)
        val antipodal_cell = S2Cell(-p)
        val dist = cell.getMaxDistance(antipodal_cell)
        assertThat(dist).isEqualTo(S1ChordAngle.straight())
    }

    @Test
    fun getMaxDistanceToCell() {
        for (i in 0 until 1000) {
            val cell = S2Cell(randomCellId())
            val test_cell = S2Cell(randomCellId());
            val antipodal_leaf_id = S2CellId.fromPoint(-test_cell.getCenter())
            val antipodal_test_cell = S2Cell(antipodal_leaf_id.parent(test_cell.level()))

            val dist_from_min = S1ChordAngle.straight() - cell.getDistance(antipodal_test_cell)
            val dist_from_max = cell.getMaxDistance(test_cell)
            assertThat(dist_from_max.radians()).isCloseTo(dist_from_min.radians(), Offset.offset(1e-8));
        }
    }
}

