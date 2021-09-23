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
package dilivia.s2.index.cell

import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Earth
import dilivia.s2.S2GeometryTestCase.checkDistanceResults
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.coords.S2Coords
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.index.shape.FractalLoopShapeIndexFactory
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2LatLngRect
import dilivia.s2.region.S2RegionCoverer
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import java.lang.Integer.min
import java.lang.Math.pow


typealias TestingResult = Pair<S2MinDistance, ValuedCell<Int>>

class S2ClosestCellQueryUnitTest {

    val kNumIndexes = 20
    val kNumCells =   100
    val kNumQueries = 100

    @Test
    fun noCells() {
        val index = S2CellIndex()
        index.build()
        val query = S2ClosestCellQuery(index)
        val target = S2ClosestCellQuery.PointTarget(S2Point(1, 0, 0))
        val result = query.findClosestCell(target)
        assertThat(result.distance.value).isEqualTo(S1ChordAngle.infinity())
        assertThat(result.cellId).isEqualTo(S2CellId.none)
        assertThat(result.value).isEqualTo(-1)
        assertThat(result.isEmpty()).isTrue()
        assertThat(query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())
    }

    @Test
    fun optionsNotModified() {
        // Tests that findClosestCell(), GetDistance(), and IsDistanceLess() do not
        // modify query.options, even though all of these methods have their own
        // specific options requirements.
        val options = S2ClosestCellQuery.Options()
        options.maxResults = 3
        options.setMaxDistance(S1ChordAngle.degrees(3))
        options.maxError = S1ChordAngle.degrees(0.001)
        val index = S2CellIndex()
        index.add(S2CellId.fromPoint(S2TextParser.makePoint("1:1")), 1)
        index.add(S2CellId.fromPoint(S2TextParser.makePoint("1:2")), 2)
        index.add(S2CellId.fromPoint(S2TextParser.makePoint("1:3")), 3)
        index.build()
        val query = S2ClosestCellQuery(index, options)
        val target = S2ClosestCellQuery.PointTarget(S2TextParser.makePoint("2:2"))
        assertThat(query.findClosestCell(target).value).isEqualTo(2)
        assertThat(query.getDistance(target).degrees()).isCloseTo(1.0, Offset.offset(1e-7))
        assertThat(query.isDistanceLess(target, S1ChordAngle.degrees(1.5))).isTrue()

        // Verify that none of the options above were modified.
        assertThat(query.options.maxResults).isEqualTo(options.maxResults)
        assertThat(query.options.maxDistance).isEqualTo(options.maxDistance)
        assertThat(query.options.maxError).isEqualTo(options.maxError)
    }

    @Test
    fun distanceEqualToLimit() {
        // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
        // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
        // the distance to the target exactly equals the chosen limit.
        val id0 = S2CellId.fromPoint(S2TextParser.makePoint("23:12"))
        val id1 = S2CellId.fromPoint(S2TextParser.makePoint("47:11"))
        val index = S2CellIndex()
        index.add(id0, 0)
        index.build()
        val query = S2ClosestCellQuery(index)

        // Start with two identical cells and a zero distance.
        val target0 = S2ClosestCellQuery.CellTarget(S2Cell(id0))
        val dist0 = S1ChordAngle.zero()
        assertThat(query.isDistanceLess(target0, dist0)).isFalse()
        assertThat(query.isDistanceLessOrEqual(target0, dist0)).isTrue()
        assertThat(query.isConservativeDistanceLessOrEqual(target0, dist0)).isTrue()

        // Now try two cells separated by a non-zero distance.
        val target1 = S2ClosestCellQuery.CellTarget(S2Cell(id1))
        val dist1 = S2Cell(id0).getDistance(S2Cell(id1))
        assertThat(query.isDistanceLess(target1, dist1)).isFalse()
        assertThat(query.isDistanceLessOrEqual(target1, dist1)).isTrue()
        assertThat(query.isConservativeDistanceLessOrEqual(target1, dist1)).isTrue()
    }

    @Test
    fun targetPointInsideIndexedCell() {
        // Tests a target point in the interior of an indexed cell.
        val cell_id = S2TextParser.makeCellId("4/012")
        val index = S2CellIndex()
        index.add(cell_id, 1)
        index.build()
        val query = S2ClosestCellQuery(index)
        val target = S2ClosestCellQuery.PointTarget(cell_id.toPoint())
        val result = query.findClosestCell(target)
        assertThat(result.distance.value).isEqualTo(S1ChordAngle.zero())
        assertThat(result.cellId).isEqualTo(cell_id)
        assertThat(result.value).isEqualTo(1)
    }

    @Test
    fun emptyTargetOptimized() {
        // Ensure that the optimized algorithm handles empty targets when a distance
        // limit is specified.
        val index = S2CellIndex()
        repeat(1000) { i ->
            index.add(S2Random.randomCellId(), i)
        }
        index.build()
        val query = S2ClosestCellQuery(index)
        query.options.setMaxDistance(S1Angle.radians(1e-5))
        val targetIndex = MutableS2ShapeIndex()
        val target = S2ClosestCellQuery.ShapeIndexTarget(targetIndex)
        assertThat(query.findClosestCells(target).size).isEqualTo(0)
    }

    @Test
    fun emptyCellUnionTarget() {
        // Verifies that distances are measured correctly to empty S2CellUnion
        // targets.
        val target = S2ClosestCellQuery.CellUnionTarget(S2CellUnion())

        val emptyIndex = S2CellIndex()
        emptyIndex.build()
        val emptyQuery = S2ClosestCellQuery(emptyIndex)
        assertThat(emptyQuery.getDistance(target)).isEqualTo(S1ChordAngle.infinity())

        val oneCellIndex = S2CellIndex()
        oneCellIndex.add(S2TextParser.makeCellId("1/123123"), 1)
        oneCellIndex.build()
        val oneCellQuery = S2ClosestCellQuery(oneCellIndex)
        assertThat(oneCellQuery.getDistance(target)).isEqualTo(S1ChordAngle.infinity())
    }

    // An abstract class that adds cells to an S2CellIndex for benchmarking.
    interface CellIndexFactory {

        // Requests that approximately "num_cells" cells located within the given
        // S2Cap bound should be added to "index".
        fun addCells(index_cap: S2Cap, num_cells: Int, index: S2CellIndex)

    }

    // Generates a cloud of points that approximately fills the given S2Cap, and
    // adds a leaf S2CellId for each one.
    class PointCloudCellIndexFactory : CellIndexFactory {

        override fun addCells(index_cap: S2Cap, num_cells: Int, index: S2CellIndex) {
            for (i in 0 until num_cells) {
                index.add(S2CellId.fromPoint(S2Random.samplePoint(index_cap)), i)
            }
        }

    }

    // Generates a collection of S2Caps that are approximately within the given
    // "index_cap", generates a covering with "max_cells_per_cap" for each one,
    // and adds the coverings to the index.  The radius of each cap is chosen
    // randomly such that the total area of the coverings is approximately
    // "cap_density" times the area of "index_cap".  In other words, a random
    // point inside "index_cap" is likely to intersect about "cap_density"
    // coverings (within a factor of 2 or so).
    class CapsCellIndexFactory(val max_cells_per_cap: Int, val cap_density: Double) : CellIndexFactory {

        override fun addCells(index_cap: S2Cap, num_cells: Int, index: S2CellIndex) {
            // All of this math is fairly approximate, since the coverings don't have
            // exactly the given number of cells, etc.
            val num_caps = (num_cells - 1) / max_cells_per_cap + 1
            val max_area = index_cap.area * cap_density / num_caps
            for (i in 0 until num_caps) {
                // The coverings are bigger than the caps, so we compensate for this by
                // choosing the cap area randomly up to the limit value.
                val cap = S2Cap.fromCenterArea(S2Random.samplePoint(index_cap), S2Random.randomDouble() * max_area)
                val coverer = S2RegionCoverer()
                coverer.setMaxCells(max_cells_per_cap)
                index.add(coverer.getCovering(cap), i)
            }
        }
    }

    // Converts to the format required by CheckDistanceResults() in s2testing.h.
    fun convertResults(results: List<S2ClosestCellQueryResult<S2MinDistance>>): List<TestingResult> {
        val testing_results = mutableListOf<TestingResult>()
        for (result in results) {
            testing_results.add(TestingResult(result.distance, ValuedCell(result.cellId, result.value)))
        }
        return testing_results
    }

    companion object {

        private val logger = KotlinLogging.logger(S2ClosestCellQueryUnitTest::class.qualifiedName!!)

        // The approximate radius of S2Cap from which query cells are chosen.
        val kTestCapRadius = S1Angle.radians(S2Earth.kmToRadians(10.0))

        // Use "query" to find the closest cell(s) to the given target, and extract
        // the query results into the given vector.  Also verify that the results
        // satisfy the search criteria.
        fun getClosestCells(
            target: S2MinDistanceTarget,
            query: S2ClosestCellQuery,
            results: MutableList<TestingResult>
        ) {
            val query_results = query.findClosestCells(target)
            assertThat(query_results.size).isLessThanOrEqualTo(query.options.maxResults)
            val region = query.options.region
            if (region == null && query.options.maxDistance.value == S1ChordAngle.infinity()) {
                // We can predict exactly how many cells should be returned.
                assertThat(query_results.size).isEqualTo(min(query.options.maxResults, query.index.numCells))
            }
            for (result in query_results) {
                // Check that the cell satisfies the region() condition.
                if (region != null) assertThat(region.mayIntersect(S2Cell(result.cellId))).isTrue()

                // Check that it satisfies the max_distance() condition.
                assertThat(result.distance.value).isLessThan(query.options.maxDistance.value)
                results.add(TestingResult(result.distance, ValuedCell(result.cellId, result.value)))
            }
        }

        fun testFindClosestCells(target: S2MinDistanceTarget, query: S2ClosestCellQuery) {
            val expected = mutableListOf<TestingResult>()
            val actual = mutableListOf<TestingResult>()
            query.options.useBruteForce = true
            getClosestCells(target, query, expected)
            query.options.useBruteForce = false
            getClosestCells(target, query, actual)

            logger.debug { """
                |
                |Test Find Closest Cells
                |=================================================================
                |target: $target
                |options: ${query.options}
                |Cells: 
                |-------
                |
                |-------
                |Bruteforce:
                |-------
                |${expected.joinToString("\n")}
                |
                |-------
                |Optimized
                |-------
                |${actual.joinToString(("\n"))}
                |----------
            """.trimMargin() }


            assertThat(
                checkDistanceResults(
                    S2MinDistanceFactory,
                    expected,
                    actual,
                    query.options.maxResults,
                    query.options.maxDistance,
                    query.options.maxError
                )
            )
                .withFailMessage("""max_results=${query.options.maxResults}, max_distance=${query.options.maxDistance}, max_error=${query.options.maxError}
                    |result=$actual
                """.trimMargin())
                .isTrue()

            if (expected.isEmpty()) return

            // Note that when options.maxError > 0, expected[0].distance() may not be
            // the minimum distance.  It is never larger by more than max_error(), but
            // the actual value also depends on max_results().
            //
            // Here we verify that GetDistance() and IsDistanceLess() return results
            // that are consistent with the max_error() setting.
            val max_error = query.options.maxError
            val min_distance = expected[0].first
            assertThat(query.getDistance(target)).isLessThanOrEqualTo(min_distance.value + max_error)

            // Test IsDistanceLess().
            assertThat(query.isDistanceLess(target, min_distance.value - max_error))
                .withFailMessage("Distance ${query.getDistance(target)} is less that ${min_distance.value - max_error} (min_distance.value - max_error)")
                .isFalse()
            assertThat(query.isConservativeDistanceLessOrEqual(target, min_distance.value)).isTrue()
        }

        // The running time of this test is proportional to
        //    (num_indexes + num_queries) * num_cells.
        // (Note that every query is checked using the brute force algorithm.)
        fun testWithIndexFactory(factory: CellIndexFactory, num_indexes: Int, num_cells: Int, num_queries: Int, seed: Int = 0) {
            // Build a set of S2CellIndexes containing the desired geometry.
            val index_caps = mutableListOf<S2Cap>()
            val indexes = mutableListOf<S2CellIndex>()
            repeat(num_indexes) { i ->
                S2Random.reset(seed + i)
                index_caps.add(S2Cap.fromCenterAngle(S2Random.randomPoint(), kTestCapRadius))
                val index = S2CellIndex()
                factory.addCells(index_caps.last(), num_cells, index)
                index.build()
                indexes.add(index)
            }
            repeat (num_queries) { i ->
                S2Random.reset(seed + i)
                val i_index = S2Random.randomInt(num_indexes)
                val index_cap = index_caps [i_index]

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                val query_radius = 2 * index_cap.radius()
                val query_cap = S2Cap.fromCenterAngle (index_cap.center, query_radius)
                val query = S2ClosestCellQuery(indexes[i_index])

                // Occasionally we don't set any limit on the number of result cells.
                // (This may return all cells if we also don't set a distance limit.)
                if (!S2Random.oneIn(10)) {
                    query.options.maxResults = 1 + S2Random.randomInt(10)
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Random.oneIn(3)) {
                    query.options.setMaxDistance(S2Random.randomDouble() * query_radius)
                }
                if (S2Random.oneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.options.setMaxError(S1Angle.radians(pow(1e-4, S2Random.randomDouble()) * query_radius.radians))
                }
                val filter_rect = S2LatLngRect.fromCenterSize(
                    S2LatLng.fromPoint(S2Random.samplePoint(query_cap)),
                    S2LatLng.fromLatLng(
                        S2Random.randomDouble() * kTestCapRadius,
                        S2Random.randomDouble() * kTestCapRadius
                    )
                )
                if (S2Random.oneIn(5)) {
                    query.options.region = filter_rect
                }


                val target_type = S2Random.randomInt(4) // TODO(fmeurisse) add S2ShapeIndex target

                if (target_type == 0) {
                    // Find the cells closest to a given point.
                    val point = S2Random.samplePoint(query_cap)
                    val target = S2ClosestCellQuery.PointTarget(point)
                    testFindClosestCells(target, query)
                }
                else if (target_type == 1) {
                    // Find the cells closest to a given edge.
                    val a = S2Random.samplePoint(query_cap)
                    val b = S2Random.samplePoint(
                        S2Cap.fromCenterAngle(a, pow(1e-4, S2Random.randomDouble()) * query_radius)
                    )
                    val target = S2ClosestCellQuery.EdgeTarget(a, b)
                    testFindClosestCells(target, query)
                } else if (target_type == 2) {
                    // Find the cells closest to a given cell.
                    val min_level = S2Coords.projection.kMaxDiag.getLevelForMaxValue(query_radius.radians)
                    val level = min_level + S2Random.randomInt(S2CellId.kMaxLevel - min_level + 1)
                    val a = S2Random.samplePoint(query_cap)
                    val cell = S2Cell(S2CellId.fromPoint(a).parent(level))
                    val target = S2ClosestCellQuery.CellTarget(cell)
                    testFindClosestCells(target, query)
                } else if (target_type == 3) {
                    // Find the cells closest to an S2Cap covering.
                    val cap = S2Cap.fromCenterAngle(S2Random.samplePoint(query_cap), 0.1 * pow(1e-4, S2Random.randomDouble()) * query_radius)
                    val coverer = S2RegionCoverer()
                    coverer.setMaxCells(4)
                    val target = S2ClosestCellQuery.CellUnionTarget(coverer.getCovering(cap))
                    testFindClosestCells(target, query)
                } else {
                    assertThat(target_type).isEqualTo(4)
                    val target_index = MutableS2ShapeIndex()
                    FractalLoopShapeIndexFactory().addEdges(index_cap, 100, target_index)
                    val target = S2ClosestCellQuery.ShapeIndexTarget(target_index)
                    target.setIncludeInteriors(S2Random.oneIn(2))
                    testFindClosestCells(target, query)
                }
            }
        }
    }



    @Test
    fun pointCloudCells() {
        testWithIndexFactory(
            PointCloudCellIndexFactory(),
            kNumIndexes, kNumCells, kNumQueries
        );
    }

    @Test
    fun capsCells() {
        testWithIndexFactory(
            CapsCellIndexFactory(
                16 /* max_cells_per_cap*/,
                0.1 /*density*/
            ),
            kNumIndexes, kNumCells, kNumQueries
        );
    }

    @Test
    fun conservativeCellDistanceIsUsed() {
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2ShapeIndex cells.
        for (seed in listOf(32, 109, 253, 342, 948, 1535, 1884, 1887, 2133)) {
            testWithIndexFactory(PointCloudCellIndexFactory(), 5, 100, 10, seed)
        }
    }


    @Test
    fun test() {
        val index = S2CellIndex()
        index.add(S2CellId.fromDebugString("1/303002302200222"), 1)
        index.add(S2CellId.fromDebugString("1/303002302201223"), 1)
        index.add(S2CellId.fromDebugString("1/30300230220123"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022013"), 1)
        index.add(S2CellId.fromDebugString("1/303002302202"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022030"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022031"), 1)
        index.add(S2CellId.fromDebugString("1/30300230220320"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022033"), 1)
        index.add(S2CellId.fromDebugString("1/30300230221003"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022101"), 1)
        index.add(S2CellId.fromDebugString("1/303002302210211"), 1)
        index.add(S2CellId.fromDebugString("1/303002302213"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022200"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022311"), 1)
        index.add(S2CellId.fromDebugString("1/3030023022312"), 1)
        index.add(S2CellId.fromDebugString("1/30300233111230012"), 0)
        index.add(S2CellId.fromDebugString("1/303002331112300130"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123002"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123003"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123010"), 0)
        index.add(S2CellId.fromDebugString("1/30300233111230113"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123012"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123013"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123020"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123021"), 0)
        index.add(S2CellId.fromDebugString("1/30300233111230220"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123023"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123030"), 0)
        index.add(S2CellId.fromDebugString("1/3030023311123031"), 0)
        index.add(S2CellId.fromDebugString("1/30300233111230320"), 0)
        index.add(S2CellId.fromDebugString("1/30300233111230321"), 0)
        index.build()

        val options = S2ClosestCellQuery.Options()
        options.maxResults = 3
        options.maxDistance = S2MinDistance(S1ChordAngle.infinity())
        options.maxError = S1ChordAngle.zero()
        options.useBruteForce = true
        options.region = null

        val query = S2ClosestCellQuery(index, options)

        val target = S2ClosestCellQuery.PointTarget(
            S2LatLng.fromLatLng(S1Angle.degrees(24.45708934580098), S1Angle.degrees(81.89577772365111)).toPoint()
        )
        testFindClosestCells(target, query)

//        var bruteForceResults = query.findClosestCells(target)
//
//
//        options.useBruteForce = false
//        var optimizedResults = query.findClosestCells(target)
//
//
//        println("Bruteforce:")
//        println("-----------------------------")
//        bruteForceResults.forEach { println(it) }
//        println("Optimized:")
//        println("-----------------------------")
//        optimizedResults.forEach { println(it) }
    }

}
