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
package dilivia.s2.builder.snap

import dilivia.ComparisonChain
import dilivia.PreConditions
import dilivia.collections.get
import dilivia.collections.sortAndRemoveDuplicates
import dilivia.collections.sortAndRemoveDuplicatesWith
import dilivia.math.M_PI
import dilivia.math.M_SQRT1_2
import dilivia.math.toExactFloat
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.max
import dilivia.s2.S1Angle.Companion.min
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng
import dilivia.s2.S2Measures
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.region.S2Cell
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.IEEErem
import kotlin.math.abs
import kotlin.math.min
import kotlin.math.roundToLong

// S2CellIdMinEdgeSeparationFunction defines an objective function that will
// be optimized by GetS2CellIdMinEdgeSeparation() by finding worst-case
// configurations of S2CellIds.  We use this to find the worst cases under
// various conditions (e.g., when the minimum snap radius at a given level is
// being used).  The objective function is called for a specific configuration
// of vertices that are snapped at the given S2CellId level.  "edge_sep" is
// the edge-vertex distance that is achieved by this configuration, and
// "min_snap_radius" and "max_snap_radius" are the minimum and maximum snap
// radii for which this configuration is valid (i.e., where the desired
// snapping will take place).
interface S2CellIdMinEdgeSeparationFunction {
    fun apply(level: Int, edge_sep: S1Angle, min_snap_radius: S1Angle, max_snap_radius: S1Angle): Double
}

// A scaled S2LatLng with integer coordinates, similar to E7 coordinates,
// except that the scale is variable (see LatLngConfig below).
typealias IntLatLng = Pair<Long, Long>

fun IntLatLng.isValid(scale: Long): Boolean {
    // A coordinate value of "scale" corresponds to 180 degrees.
    return (abs(this[0]) <= scale / 2 && abs(this[1]) <= scale)
}

fun IntLatLng.hasValidVertices(scale: Long): Boolean {
    // Like IsValid, but excludes latitudes of 90 and longitudes of 180.
    // A coordinate value of "scale" corresponds to 180 degrees.
    return (abs(this[0]) < scale / 2 && abs(this[1]) < scale)
}

fun IntLatLng.rescale(scale_factor: Double): IntLatLng {
    return IntLatLng(
        (scale_factor * this[0]).roundToLong(),
        (scale_factor * this[1]).roundToLong()
    )
}

fun IntLatLng.toPoint(scale: Long): S2Point {
    return S2LatLng.fromRadians(
        this[0] * (M_PI / scale),
        this[1] * (M_PI / scale)
    ).toPoint()
}

fun IntLatLng.getVertex(scale: Long, i: Int): S2Point {
    // Return the points in CCW order starting from the lower left.
    val dlat = if (i == 0 || i == 3) -1 else 1
    val dlng = if (i == 0 || i == 1) -1 else 1
    return (2L * this + IntLatLng(dlat.toLong(), dlng.toLong())).toPoint(2 * scale)
}

operator fun IntLatLng.plus(other: IntLatLng): IntLatLng =
    IntLatLng(this.first + other.first, this.second + other.second)

operator fun IntLatLng.times(other: Long): IntLatLng = IntLatLng(this.first * other, this.second * other)
operator fun Long.times(other: IntLatLng): IntLatLng = other * this

operator fun IntLatLng.compareTo(other: IntLatLng): Int = IntLatLngComparator.compare(this, other)

object IntLatLngComparator : Comparator<IntLatLng> {
    override fun compare(o1: Pair<Long, Long>, o2: Pair<Long, Long>): Int = ComparisonChain.start()
        .compare(o1.first, o2.first)
        .compare(o1.second, o2.second)
        .result()
}


// A triple of scaled S2LatLng coordinates.  The coordinates are multiplied by
// (M_PI / scale) to convert them to radians.
data class LatLngConfig(val scale: Long, var ll0: IntLatLng, var ll1: IntLatLng, var ll2: IntLatLng) :
    Comparable<LatLngConfig> {

    override fun compareTo(other: LatLngConfig): Int {
        PreConditions.requireEQ(scale, other.scale)
        return ComparisonChain.start()
            .compare(ll0, other.ll0, IntLatLngComparator)
            .compare(ll1, other.ll1, IntLatLngComparator)
            .compare(ll2, other.ll2, IntLatLngComparator)
            .result()
    }

}

interface LatLngMinEdgeSeparationFunction {
    fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double
}

class SnapFunctionsUnitTest {

    private val logger = KotlinLogging.logger { }

    @Test
    fun levelToFromSnapRadius() {
        for (level in 0..S2CellId.kMaxLevel) {
            val radius = S2CellIdSnapFunction.minSnapRadiusForLevel(level)
            assertThat(S2CellIdSnapFunction.levelForMaxSnapRadius(radius)).isEqualTo(level)
            assertThat(S2CellIdSnapFunction.levelForMaxSnapRadius(0.999.times(radius))).isEqualTo(
                min(
                    level + 1,
                    S2CellId.kMaxLevel
                )
            )
        }
        assertThat(S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(5))).isEqualTo(0)
        assertThat(S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(1e-30))).isEqualTo(S2CellId.kMaxLevel)
    }

    @Test
    fun snapPoint() {
        for (iter in 0 until 1000) {
            for (level in 0..S2CellId.kMaxLevel) {
                // This checks that points are snapped to the correct level, since
                // S2CellId centers at different levels are always different.
                val f = S2CellIdSnapFunction(level)
                val p = S2Random.randomCellId(level).toPoint()
                assertThat(f.snapPoint(p)).isEqualTo(p)
            }
        }
    }

    @Test
    fun intLatLngSnapFunctionExponentToFromSnapRadius() {
        var exponent = IntLatLngSnapFunction.kMinExponent

        while (exponent <= IntLatLngSnapFunction.kMaxExponent) {
            val radius = IntLatLngSnapFunction.minSnapRadiusForExponent(exponent)
            assertThat(IntLatLngSnapFunction.exponentForMaxSnapRadius(radius)).isEqualTo(exponent)
            assertThat(IntLatLngSnapFunction.exponentForMaxSnapRadius(0.999 * radius)).isEqualTo(
                min(
                    exponent + 1,
                    IntLatLngSnapFunction.kMaxExponent
                )
            )
            ++exponent
        }
        assertThat(IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(5))).isEqualTo(IntLatLngSnapFunction.kMinExponent)
        assertThat(IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(1e-30))).isEqualTo(
            IntLatLngSnapFunction.kMaxExponent
        )
    }

    @Test
    fun intLatLngSnapFunctionSnapPoint() {
        repeat(1000) {
            // Test that IntLatLngSnapFunction does not modify points that were
            // generated using the S2LatLng::From{E5,E6,E7} methods.  This ensures
            // that both functions are using bitwise-compatible conversion methods.
            val p = S2Random.randomPoint()
            val ll = S2LatLng.fromPoint(p)
            val p5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5()).toPoint()
            assertThat(IntLatLngSnapFunction(5).snapPoint(p5)).isEqualTo(p5)
            val p6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6()).toPoint()
            assertThat(IntLatLngSnapFunction(6).snapPoint(p6)).isEqualTo(p6)
            val p7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7()).toPoint()
            assertThat(IntLatLngSnapFunction(7).snapPoint(p7)).isEqualTo(p7)

            // Make sure that we're not snapping using some lower exponent.
            val p7not6 = S2LatLng.fromE7(10 * ll.lat().e6() + 1, 10 * ll.lng().e6() + 1).toPoint()
            assertThat(IntLatLngSnapFunction(6).snapPoint(p7not6)).isNotEqualTo(p7not6)
        }
    }

    val kSearchRootId = S2CellId.fromFace(0)
    val kSearchFocusId = S2CellId.fromFace(0).child(3)

    fun getMaxVertexDistance(p: S2Point, id: S2CellId): S1Angle {
        val cell = S2Cell(id)
        return max(
            max(
                S1Angle(p, cell.getVertex(0)),
                S1Angle(p, cell.getVertex(1))
            ),
            max(
                S1Angle(p, cell.getVertex(2)),
                S1Angle(p, cell.getVertex(3))
            )
        )
    }

    // Helper function that computes the vertex separation between "id0" and its
// neighbors.
    fun updateS2CellIdMinVertexSeparation(id0: S2CellId, scores: MutableList<Pair<Double, S2CellId>>) {
        val site0 = id0.toPoint()
        val nbrs = mutableListOf<S2CellId>()
        id0.appendAllNeighbors(id0.level(), nbrs)
        for (id1 in nbrs) {
            val site1 = id1.toPoint()
            val vertex_sep = S1Angle(site0, site1)
            val max_snap_radius = getMaxVertexDistance(site0, id1)
            assertThat(max_snap_radius).isGreaterThanOrEqualTo(S2CellIdSnapFunction.minSnapRadiusForLevel(id0.level()))
            val r = vertex_sep / max_snap_radius
            scores.add(Pair(r.radians, id0))
        }
    }

    fun getS2CellIdMinVertexSeparation(level: Int, best_cells: MutableSet<S2CellId>): Double {
        // The worst-case separation ratios always occur when the snap_radius is not
        // much larger than the minimum, since this allows the site spacing to be
        // reduced by as large a fraction as possible.
        //
        // For the minimum vertex separation ratio, we choose a site and one of its
        // 8-way neighbors, then look at the ratio of the distance to the center of
        // that neighbor to the distance to the furthest corner of that neighbor
        // (which is the largest possible snap radius for this configuration).
        val scores = mutableListOf<Pair<Double, S2CellId>>()
        if (level == 0) {
            updateS2CellIdMinVertexSeparation(kSearchRootId, scores)
        } else {
            for (parent in best_cells) {
                var id0 = parent.childBegin()
                while (id0 != parent.childEnd()) {
                    updateS2CellIdMinVertexSeparation(id0, scores)
                    id0 = id0.next()
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and keep
        // the best "num_to_keep" of them to seed the next round.
        scores.sortAndRemoveDuplicatesWith { o1, o2 ->
            ComparisonChain.start()
                .compare(o1.first, o2.first)
                .compare(o1.second, o2.second)
                .result()
        }

        best_cells.clear()
        var num_to_keep = 300
        var num_to_print = 1
        for (entry in scores) {
            val id = entry.second
            if (--num_to_print >= 0) {
                val uv = id.centerUV
                logger.debug(
                    String.format(
                        "Level %2d: min_vertex_sep_ratio = %.15f u=%.6f v=%.6f %s\n",
                        level, entry.first, uv[0], uv[1], id.toToken()
                    )
                )
            }
            if (kSearchFocusId.contains(id) || id.contains(kSearchFocusId)) {
                if (best_cells.add(id) && --num_to_keep <= 0) break
            }
        }
        return scores[0].first
    }

    @Test
    fun minVertexSeparationSnapRadiusRatio() {
        // The purpose of this "test" is to compute a lower bound to the fraction
        // (min_vertex_separation() / snap_radius()).  Essentially this involves
        // searching for two adjacent cells A and B such when one of the corner
        // vertices of B is snapped to the center of B, the distance to the center
        // of A decreases as much as possible.  In other words, we want the ratio
        //
        //   distance(center(A), center(B)) / distance(center(A), vertex(B))
        //
        // to be as small as possible.  We do this by considering one cell level at
        // a time, and remembering the cells that had the lowest ratios.  When we
        // proceed from one level to the next, we consider all the children of those
        // cells and keep the best ones.
        //
        // The reason we can restrict the search to children of cells at the
        // previous level is that the ratio above is essentially a function of the
        // local distortions created by projecting the S2 cube space onto the
        // sphere.  These distortions change smoothly over the sphere, so by keeping
        // a fairly large number of candidates ("num_to_keep"), we are essentially
        // keeping all the neighbors of the optimal cell as well.
        var best_score = 1e10
        val best_cells = mutableSetOf<S2CellId>()
        for (level in 0..S2CellId.kMaxLevel) {
            val score = getS2CellIdMinVertexSeparation(level, best_cells)
            best_score = min(best_score, score)
        }
        println(String.format("min_vertex_sep / snap_radius ratio: %.15f", best_score))
    }

    fun getCircumRadius(a: S2Point, b: S2Point, c: S2Point): S1Angle {
        // We return this value is the circumradius is very large.
        val kTooBig = S1Angle.radians(M_PI)
        val turn_angle = S2Measures.turnAngle(a, b, c)
        if (abs(turn_angle.IEEErem(M_PI)) < 1e-2) return kTooBig

        val a2 = (b - c).norm2().toExactFloat() //.toBigDecimal(MathContext.UNLIMITED)
        val b2 = (c - a).norm2().toExactFloat() //.toBigDecimal(MathContext.UNLIMITED)
        val c2 = (a - b).norm2().toExactFloat() //.toBigDecimal(MathContext.UNLIMITED)
        if (a2 > 2.0.toExactFloat() || b2 > 2.0.toExactFloat() || c2 > 2.0.toExactFloat()) return kTooBig
        val ma = a2 * (b2 + c2 - a2);
        val mb = b2 * (c2 + a2 - b2);
        val mc = c2 * (a2 + b2 - c2);
        val p = (ma * a.toExactFloat() + mb * b.toExactFloat() + mc * c.toExactFloat()).toDecimal128() / (ma + mb + mc)
        return S1Angle(p.toDouble(), a)
    }


    fun getNeighbors(id: S2CellId): List<S2CellId> {
        val kNumLayers = 2;
        val nbrs = mutableListOf<S2CellId>()
        nbrs.add(id);
        for (layer in 0 until kNumLayers) {
            val new_nbrs = mutableListOf<S2CellId>()
            for (nbr in nbrs) {
                nbr.appendAllNeighbors(id.level(), new_nbrs)
            }
            nbrs.addAll(new_nbrs)
            nbrs.removeAll { cellId -> cellId == id }
            nbrs.sortAndRemoveDuplicates()
        }
        return nbrs;
    }

    // Returns the minimum value of the given objective function over sets of
    // nearby vertices that are designed to minimize the edge-vertex separation
    // when an edge is snapped.
    private fun getS2CellIdMinEdgeSeparation(
        label: String,
        objective: S2CellIdMinEdgeSeparationFunction,
        level: Int,
        best_cells: MutableSet<S2CellId>
    ): Double {
        // To find minimum edge separations, we choose a cell ("id0") and two nearby
        // cells ("id1" and "id2"), where "nearby" is defined by GetNeighbors().
        // Let "site0", "site1", and "site2" be the centers of these cells.  The
        // idea is to consider an input edge E that intersects the Voronoi regions
        // of "site1" and "site2" (and therefore snaps to an edge E' between these
        // sites) but does not not intersect the Voronoi region of "site0" (and
        // therefore can't be snapped to site0).  The goal is to search for snapped
        // edges E' that approach site0 as closely as possible.
        //
        // To do this, we first compute the circumradius of the three cell centers
        // ("site0", "site1", and "site2"); this is the minimum snap radius in order
        // for it to be possible to construct an edge E that snaps to "site1" and
        // "site2" but not to "site0".  We also compute the distance from "site0" to
        // the snapped edge.  Next we find the corner vertex of "id1" and "id2" that
        // is furthest from "site0"; the smaller of these two distances is the
        // maximum snap radius such that "site1" and "site2" can be chosen as
        // sites after choosing "site0".  If the maximum is less than the minimum,
        // then this configuration is rejected; otherwise we evaluate the given
        // objective function and keep the configurations that result in the
        // smallest values.
        //
        // The optimization process works by keeping track of the set of S2CellIds
        // that yielded the best results at the previous level, and exploring all
        // the nearby neighbor combinations of the children of those cells at the
        // next level.  In order to get better coverage, we keep track of the best
        // score and configuration (i.e. the two neighboring cells "id1" and "id2")
        // for each initial cell "id0".
        val bestScores = mutableMapOf<S2CellId, Double>()
        val bestConfigs = mutableMapOf<S2CellId, Pair<S2CellId, S2CellId>>()
        for (parent: S2CellId in best_cells) {

            var id0: S2CellId = parent.childBegin(level)
            while (id0 != parent.childEnd(level)) {
                val site0 = id0.toPoint()
                val nbrs = getNeighbors(id0)
                for (id1: S2CellId in nbrs) {
                    val site1: S2Point = id1.toPoint()
                    val max_v1 = getMaxVertexDistance(site0, id1)
                    for (id2: S2CellId in nbrs) {
                        if (id2 <= id1) continue
                        val site2 = id2.toPoint()
                        val minSnapRadius = getCircumRadius(site0, site1, site2)
                        if (minSnapRadius > SnapFunction.kMaxSnapRadius()) continue
                        // Note that it is only the original points *before* snapping that
                        // need to be at least "snap_radius" away from "site0".  The points
                        // after snapping ("site1" and "site2") may be closer.
                        val maxV2 = getMaxVertexDistance(site0, id2)
                        val maxSnapRadius = min(max_v1, maxV2)
                        if (minSnapRadius > maxSnapRadius) continue
                        assertThat(maxSnapRadius).isGreaterThanOrEqualTo(S2CellIdSnapFunction.minSnapRadiusForLevel(level))

                        // This is a valid configuration, so evaluate it.
                        val edgeSep = S2EdgeDistances.getDistance(site0, site1, site2)
                        val score = objective.apply(level, edgeSep, minSnapRadius, maxSnapRadius)
                        val bestScore = bestScores[id0] ?: 0.0
                        if (bestScore == 0.0 || bestScore > score) {
                            bestScores[id0] = score
                            bestConfigs[id0] = Pair(id1, id2)
                        }
                    }
                }
                id0 = id0.next()
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and
        // generate a set of candidates for the next round by generating all the
        // 8-way neighbors of the best candidates, and keeping up to"num_to_keep" of
        // them.  The results vary slightly according to how many candidates we
        // keep, but the variations are much smaller than the conservative
        // assumptions made by the S2CellIdSnapFunction implementation.
        var numToKeep = 20 // google::DEBUG_MODE ? 20 : 100;
        var numToPrint = 3
        val sorted = mutableListOf<Pair<Double, S2CellId>>()
        for (entry in bestScores) {
            sorted.add(Pair(entry.value, entry.key))
        }
        sorted.sortWith { o1, o2 ->
            ComparisonChain.start().compare(o1.first, o2.first).compare(o1.second, o2.second).result()
        }
        best_cells.clear()
        logger.debug(String.format("Level %d:", level))
        for (entry in sorted) {
            val id = entry.second
            if (--numToPrint >= 0) {
                val uv = id.centerUV
                val nbrs = bestConfigs.getValue(id)
                logger.debug(
                    String.format(
                        "  %s = %.15f u=%7.4f v=%7.4f %s %s %s",
                        label, entry.first, uv[0], uv[1], id.toToken(), nbrs.first.toToken(), nbrs.second.toToken()
                    )
                )
            }
            val nbrs = mutableListOf(id)
            id.appendAllNeighbors(id.level(), nbrs)
            for (nbr in nbrs) {
                // The S2Cell hierarchy has many regions that are symmetrical.  We can
                // eliminate most of the "duplicates" by restricting the search to cells
                // in kS2CellIdFocus.
                if (kSearchFocusId.contains(nbr) || nbr.contains(kSearchFocusId)) {
                    if (best_cells.add(nbr) && --numToKeep <= 0) {
                        return sorted[0].first
                    }
                }
            }
        }
        return sorted[0].first
    }

    fun getS2CellIdMinEdgeSeparation(label: String, objective: S2CellIdMinEdgeSeparationFunction): Double {
        var bestScore = 1e10
        val bestCells = mutableSetOf<S2CellId>()
        bestCells.add(kSearchRootId)
        for (level in 0..S2CellId.kMaxLevel) {
            val score = getS2CellIdMinEdgeSeparation(label, objective, level, bestCells)
            bestScore = min(bestScore, score);
        }
        return bestScore
    }

    @Test
    fun minEdgeVertexSeparationForLevel() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        val score = getS2CellIdMinEdgeSeparation("min_sep_for_level", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(
                level: Int,
                edge_sep: S1Angle,
                min_snap_radius: S1Angle,
                max_snap_radius: S1Angle
            ): Double {
                return edge_sep.radians / S2Coords.projection.kMinDiag.getValue(level)
            }
        })
        logger.debug(String.format("min_edge_vertex_sep / kMinDiag ratio: %.15f", score))
    }

    @Test
    fun minEdgeVertexSeparationAtMinSnapRadius() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for the
        // special case where the minimum snap radius is being used.
        val score = getS2CellIdMinEdgeSeparation("min_sep_at_min_radius", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(
                level: Int,
                edge_sep: S1Angle,
                min_snap_radius: S1Angle,
                max_snap_radius: S1Angle
            ): Double {
                val min_radius_at_level = S2Coords.projection.kMaxDiag.getValue(level) / 2
                return if (min_snap_radius.radians <= (1 + 1e-10) * min_radius_at_level) (edge_sep.radians / S2Coords.projection.kMinDiag.getValue(
                    level
                )) else 100.0
            }
        })
        logger.debug(String.format("min_edge_vertex_sep / kMinDiag at MinSnapRadiusForLevel: %.15f", score))
    }

    @Test
    fun minEdgeVertexSeparationSnapRadiusRatio() {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        val score = getS2CellIdMinEdgeSeparation("min_sep_at_min_radius", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(
                level: Int,
                edge_sep: S1Angle,
                min_snap_radius: S1Angle,
                max_snap_radius: S1Angle
            ): Double {
                return edge_sep.radians / max_snap_radius.radians
            }
        })
        logger.debug(String.format("min_edge_vertex_sep / snap_radius ratio: %.15f", score))
    }


    fun getMaxVertexDistance(p: S2Point, ll: IntLatLng, scale: Long): S1Angle {
        return max(
            max(S1Angle(p, ll.getVertex(scale, 0)), S1Angle(p, ll.getVertex(scale, 1))),
            max(S1Angle(p, ll.getVertex(scale, 2)), S1Angle(p, ll.getVertex(scale, 3)))
        )
    }


    fun getLatLngMinVertexSeparation(old_scale: Long, scale: Long, best_configs: MutableSet<IntLatLng>): Double {
        // The worst-case separation ratios always occur when the snap_radius is not
        // much larger than the minimum, since this allows the site spacing to be
        // reduced by as large a fraction as possible.
        //
        // For the minimum vertex separation ratio, we choose a site and one of its
        // 8-way neighbors, then look at the ratio of the distance to the center of
        // that neighbor to the distance to the furthest corner of that neighbor
        // (which is the largest possible snap radius for this configuration).
        val min_snap_radius_at_scale = S1Angle.radians(M_SQRT1_2 * M_PI / scale)
        val scores = mutableListOf<Pair<Double, IntLatLng>>()
        val scale_factor = scale.toDouble() / old_scale
        for (parent in best_configs) {
            val new_parent = parent.rescale(scale_factor)
            for (dlat0 in -7L..7L) {
                val ll0 = new_parent + IntLatLng(dlat0, 0);
                if (!ll0.isValid(scale) || ll0[0] < 0) continue
                val site0 = ll0.toPoint(scale)
                for (dlat1 in 0L..2L) {
                    for (dlng1 in 0L..5L) {
                        val ll1 = ll0 + IntLatLng(dlat1, dlng1)
                        if (ll1 == ll0 || !ll1.hasValidVertices(scale)) continue
                        val max_snap_radius = getMaxVertexDistance(site0, ll1, scale)
                        if (max_snap_radius < min_snap_radius_at_scale) continue;
                        val site1 = ll1.toPoint(scale)
                        val vertex_sep = S1Angle(site0, site1)
                        val r = vertex_sep / max_snap_radius
                        scores.add(Pair(r.radians, ll0))
                    }
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and keep
        // the best "num_to_keep" of them to seed the next round.
        scores.sortAndRemoveDuplicatesWith { o1, o2 ->
            ComparisonChain.start()
                .compare(o1.first, o1.first)
                .compare(o1.second.first, o2.second.first)
                .compare(o1.second.second, o1.second.second)
                .result()
        }
        best_configs.clear()
        var num_to_keep = 100
        var num_to_print = 1
        for (entry in scores) {
            if (--num_to_print >= 0) {
                logger.debug(
                    String.format(
                        "Scale %14d: min_vertex_sep_ratio = %.15f, %s",
                        scale,
                        entry.first,
                        S2TextParser.toString(entry.second.toPoint(scale))
                    )
                )
            }
            if (best_configs.add(entry.second) && --num_to_keep <= 0) break
        }
        return scores[0].first;
    }

    @Test
    fun intLatLngSnapFunctionMinVertexSeparationSnapRadiusRatio() {
        var best_score = 1e10;
        val best_configs = mutableSetOf<IntLatLng>()
        var scale = 18L
        for (lat0 in 0L..9L) {
            best_configs.add(IntLatLng(lat0, 0L))
        }
        for (exp in 0L..10L) {
            val score = getLatLngMinVertexSeparation(scale, 10 * scale, best_configs)
            best_score = min(best_score, score)
            scale *= 10
        }
        logger.debug(String.format("min_vertex_sep / snap_radius ratio: %.15f", best_score))
    }

    fun getLatLngMinEdgeSeparation(
        label: String,
        objective: LatLngMinEdgeSeparationFunction,
        scale: Long,
        best_configs: MutableList<LatLngConfig>
    ): Double {
        val min_snap_radius_at_scale = S1Angle.radians(M_SQRT1_2 * M_PI / scale)
        val scores = mutableListOf<Pair<Double, LatLngConfig>>()
        for (parent in best_configs) {
            // To reduce duplicates, we require that site0 always has longitude 0.
            PreConditions.checkEQ(0, parent.ll0[1])
            val scale_factor = scale.toDouble() / parent.scale
            parent.ll0 = parent.ll0.rescale(scale_factor)
            parent.ll1 = parent.ll1.rescale(scale_factor)
            parent.ll2 = parent.ll2.rescale(scale_factor)
            for (dlat0 in -1L..1L) {
                val ll0 = parent.ll0 + IntLatLng(dlat0, 0)
                // To reduce duplicates, we require that site0.latitude >= 0.
                if (!ll0.isValid(scale) || ll0[0] < 0) continue
                val site0 = ll0.toPoint(scale)
                for (dlat1 in -1L..1L) {
                    for (dlng1 in -2L..2L) {
                        val ll1 = parent.ll1 + IntLatLng(dlat0 + dlat1, dlng1)
                        if (ll1 == ll0 || !ll1.hasValidVertices(scale)) continue
                        // Only consider neighbors within 2 latitude units of site0.
                        if (abs(ll1[0] - ll0[0]) > 2) continue

                        val site1 = ll1.toPoint(scale)
                        val max_v1 = getMaxVertexDistance(site0, ll1, scale)
                        for (dlat2 in -1L..1L) {
                            for (dlng2 in -2L..2L) {
                                val ll2 = parent.ll2 + IntLatLng(dlat0 + dlat2, dlng2)
                                if (!ll2.hasValidVertices(scale)) continue
                                // Only consider neighbors within 2 latitude units of site0.
                                if (abs(ll2[0] - ll0[0]) > 2) continue;
                                // To reduce duplicates, we require ll1 < ll2 lexicographically
                                // and site2.longitude >= 0.  (It's *not* okay to
                                // require site1.longitude >= 0, because then some configurations
                                // with site1.latitude == site2.latitude would be missed.)
                                if (ll2 <= ll1 || ll2[1] < 0) continue;

                                val site2 = ll2.toPoint(scale)
                                val min_snap_radius = getCircumRadius(site0, site1, site2)
                                if (min_snap_radius > SnapFunction.kMaxSnapRadius()) {
                                    continue
                                }
                                // Only the original points *before* snapping that need to be at
                                // least "snap_radius" away from "site0".  The points after
                                // snapping ("site1" and "site2") may be closer.
                                val max_v2 = getMaxVertexDistance(site0, ll2, scale)
                                val max_snap_radius = min(max_v1, max_v2)
                                if (min_snap_radius > max_snap_radius) continue;
                                if (max_snap_radius < min_snap_radius_at_scale) continue;

                                // This is a valid configuration, so evaluate it.
                                val edge_sep = S2EdgeDistances.getDistance(site0, site1, site2)
                                val score = objective.apply(scale, edge_sep, max_snap_radius);
                                val config = LatLngConfig(scale, ll0, ll1, ll2);
                                scores.add(Pair(score, config))
                            }
                        }
                    }
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and keep
        // the best "num_to_keep" of them to seed the next round.
        scores.sortAndRemoveDuplicatesWith { o1, o2 ->
            ComparisonChain.start().compare(o1.first, o2.first).compare(o1.second, o2.second).result()
        }

        best_configs.clear()
        var num_to_keep = 50 // google ::DEBUG_MODE ? 50 : 200;
        var num_to_print = 3
        logger.debug(String.format("Scale %d :", scale))
        for (entry in scores) {
            val config = entry.second
            val scale = config.scale
            if (--num_to_print >= 0) {
                logger.debug(
                    String.format(
                        "  %s = %.15f %s %s %s",
                        label, entry.first,
                        S2TextParser.toString(config.ll0.toPoint(scale)),
                        S2TextParser.toString(config.ll1.toPoint(scale)),
                        S2TextParser.toString(config.ll2.toPoint(scale))
                    )
                )
            }
            // Optional: filter the candidates to concentrate on a specific region
            // (e.g., the north pole).
            best_configs.add(config)
            if (--num_to_keep <= 0) break
        }
        return scores[0].first
    }

    fun getLatLngMinEdgeSeparation(label: String, objective: LatLngMinEdgeSeparationFunction): Double {
        var best_score = 1e10
        val best_configs = mutableListOf<LatLngConfig>()
        var scale = 6L  // Initially points are 30 degrees apart.
        val maxLng = scale
        val maxLat = scale / 2
        for (lat0 in 0..maxLat) {
            for (lat1 in (lat0 - 2)..min(maxLat, lat0 + 2)) {
                for (lng1 in 0..maxLng) {
                    for (lat2 in lat1..min(maxLat, lat0 + 2)) {
                        for (lng2 in 0..maxLng) {
                            val ll0 = IntLatLng(lat0, 0)
                            val ll1 = IntLatLng(lat1, lng1)
                            val ll2 = IntLatLng(lat2, lng2)
                            if (ll2 <= ll1) continue;
                            best_configs.add(LatLngConfig(scale, ll0, ll1, ll2));
                        }
                    }
                }
            }
        }
        logger.info { "Starting with ${best_configs.size} configurations" }
        var target_scale = 180L
        for (exp in 0..10) {
            while (scale < target_scale) {
                scale = min((1.8 * scale).toLong(), target_scale);
                val score = getLatLngMinEdgeSeparation(label, objective, scale, best_configs)
                if (scale == target_scale) {
                    best_score = min(best_score, score);
                }
                target_scale *= 10
            }
        }
        return best_score
    }

    @Test
    fun intLatLngSnapFunctionMinEdgeVertexSeparationForLevel() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        val score = getLatLngMinEdgeSeparation("min_sep_for_level", object : LatLngMinEdgeSeparationFunction {
            override fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double {
                val e_unit = M_PI / scale
                return edge_sep.radians / e_unit
            }
        })
        logger.debug(String.format("min_edge_vertex_sep / e_unit ratio: %.15f", score))
    }

    @Test
    fun intLatLngSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio() {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        val score = getLatLngMinEdgeSeparation("min_sep_snap_radius_ratio", object : LatLngMinEdgeSeparationFunction {
            override fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double {
                return edge_sep.radians / max_snap_radius.radians
            }
        })
        logger.debug(String.format("min_edge_vertex_sep / snap_radius ratio: %.15f", score))
    }

}
