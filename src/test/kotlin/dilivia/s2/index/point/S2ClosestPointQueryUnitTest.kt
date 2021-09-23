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
package dilivia.s2.index.point

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Earth
import dilivia.s2.S2GeometryTestCase.checkDistanceResults
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.coords.S2Coords
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryCellTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryEdgeTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryPointTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryShapeIndexTarget
import dilivia.s2.index.shape.FractalLoopShapeIndexFactory
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2LatLngRect
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import tech.units.indriya.quantity.Quantities
import tech.units.indriya.unit.Units.METRE
import java.lang.Math.pow
import kotlin.math.min
import kotlin.math.pow
import dilivia.s2.index.point.S2ClosestPointQuery as TestQuery
import dilivia.s2.index.point.S2PointIndex as TestIndex


class S2ClosestPointQueryUnitTest  {

    @Test
    fun noPoints() {
        val index = TestIndex<Int>()
        val query = TestQuery(index)
        val target = S2ClosestPointQueryPointTarget(S2Point(1, 0, 0))
        val results = query.findClosestPoints(target)
        assertThat(results.size).isEqualTo(0)
    }

    @Test
    fun manyDuplicatePoints() {
        val kNumPoints = 10000
        val kTestPoint = S2Point(1, 0, 0)
        val index = TestIndex<Int>()
        repeat(kNumPoints) { i ->
            index.add(kTestPoint, i)
        }
        val query = TestQuery(index)
        val target = S2ClosestPointQueryPointTarget(kTestPoint)
        val results = query.findClosestPoints(target)
        assertThat(results.size).isEqualTo(kNumPoints)
    }

    @Test
    fun emptyTargetOptimized() {
        // Ensure that the optimized algorithm handles empty targets when a distance
        // limit is specified.
        val index = TestIndex<Int>()
        repeat(1000) { i ->
            index.add(S2Random.randomPoint(), i)
        }
        val query = TestQuery(index, TestQuery.Options())
        query.options().setMaxDistance(S1Angle.radians(1e-5))
        val targetIndex = MutableS2ShapeIndex()
        val target = S2ClosestPointQueryShapeIndexTarget(targetIndex)
        assertThat(query.findClosestPoints(target).size).isEqualTo(0)
    }

    @Test
    fun circlePoints() {
        testWithIndexFactory(CirclePointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    @Test
    fun fractalPoints() {
        testWithIndexFactory(FractalPointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    @Test
    fun gridPoints() {
        testWithIndexFactory(GridPointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    @Test
    fun conservativeCellDistanceIsUsed() {
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2PointIndex cells.  They
        // all involve S2ShapeIndexTarget, which takes advantage of max_error() to
        // optimize its distance calculation.
        for (seed in listOf(16, 586, 589, 822, 1959, 2298, 3155, 3490, 3723, 4953)) {
            testWithIndexFactory(FractalPointIndexFactory(), 5, 100, 10, seed)
        }
    }

    companion object {

        val logger = KotlinLogging.logger(S2ClosestPointQueryUnitTest::class.java.name)

        // The approximate radius of S2Cap from which query points are chosen.
        val kTestCapRadius: S1Angle = S2Earth.toAngle(Quantities.getQuantity(10.0, METRE))

        const val kNumIndexes = 1 //10
        const val kNumPoints = 1000 //1000
        const val kNumQueries = 50

        // Use "query" to find the closest point(s) to the given target, and extract
        // the query results into the given vector.  Also verify that the results
        // satisfy the search criteria.
        fun getClosestPoints(target: S2MinDistanceTarget, query: TestQuery<Int>, results: MutableList<Pair<S2MinDistance, Int>>) {
            val queryResults = query.findClosestPoints(target)
            assertThat(queryResults.size <= query.options().getMaxResult()).isTrue()
            val region = query.options().region
            if (region == null && query.options().maxDistance.value == S1ChordAngle.infinity()) {
                // We can predict exactly how many points should be returned.
                assertThat(queryResults.size).isEqualTo(min(query.options().getMaxResult(), query.index().numPoints()))
            }
            for (result in queryResults) {
                // Check that the point satisfies the region() condition.
                if (region != null) assertThat(region.contains(result.point())).isTrue()

                // Check that it satisfies the max_distance() condition.
                assertThat(result.distance < query.options().maxDistance).isTrue()
                results.add(Pair(result.distance, result.data()))
            }
        }

        fun testFindClosestPoints(target: S2MinDistanceTarget, query: TestQuery<Int>) {
            val expected = mutableListOf<Pair<S2MinDistance, Int>>()
            val actual = mutableListOf<Pair<S2MinDistance, Int>>()
            query.options().useBruteForce = true
            getClosestPoints(target, query, expected)
            query.options().useBruteForce = false
            getClosestPoints(target, query, actual)
            assertThat(
                checkDistanceResults(S2MinDistanceFactory, expected, actual, query.options().getMaxResult(), query.options().maxDistance, query.options().maxError))
                .withFailMessage("""
               max_results  = ${query.options().getMaxResult()},
               max_distance = ${query.options().maxDistance},
               max_error    = ${query.options().maxError}
            """.trimIndent(),
                    ).isTrue()


            if (expected.isEmpty()) return

            // Note that when options.max_error() > 0, expected[0].distance may not be
            // the minimum distance.  It is never larger by more than max_error(), but
            // the actual value also depends on max_results().
            //
            // Here we verify that GetDistance() and IsDistanceLess() return results
            // that are consistent with the max_error() setting.
            val maxError = query.options().maxError
            val minDistance = expected[0].first
            assertThat(query.getDistance(target)).isLessThanOrEqualTo(minDistance.value + maxError)

            // Test IsDistanceLess().
            assertThat(query.isDistanceLess(target, minDistance.value - maxError)).isFalse()
            assertThat(query.isConservativeDistanceLessOrEqual(target, minDistance.value)).isTrue()
        }

        // (Note that every query is checked using the brute force algorithm.)
        fun testWithIndexFactory(factory: PointIndexFactory, num_indexes: Int, num_points: Int, num_queries: Int, randomSeed: Int = 1) {
            // Build a set of S2PointIndexes containing the desired geometry.
            val index_caps = mutableListOf<S2Cap>()
            val indexes = mutableListOf<TestIndex<Int>>()
            for (i in 0 until num_indexes) {
                S2Random.reset(randomSeed + i)
                index_caps.add(S2Cap.fromCenterAngle(S2Random.randomPoint(), kTestCapRadius))
                indexes.add(TestIndex())
                factory.addPoints(index_caps.last(), num_points, indexes.last())
            }
            for (i in 0 until num_queries) {

                S2Random.reset(randomSeed + i)
                val iIndex = S2Random.randomInt(num_indexes)
                val indexCap = index_caps[iIndex]

                logger.trace { """
                    |
                    |=============================================================================================
                    | Query $i - index: $iIndex
                    |=============================================================================================
                    |Index:
                    |${indexes[iIndex]}
                """.trimMargin() }

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                val queryRadius = indexCap.radius() * 2.0
                val queryCap = S2Cap.fromCenterAngle(indexCap.center, queryRadius)
                val query = TestQuery(indexes[iIndex])

                // Occasionally we don't set any limit on the number of result points.
                // (This may return all points if we also don't set a distance limit.)
                if (!S2Random.oneIn(5)) {
                    query.options().setMaxResult(1 + S2Random.randomInt(10))
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Random.oneIn(3)) {
                    query.options().setMaxDistance(queryRadius * S2Random.randomDouble())
                }
                if (S2Random.oneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.options().setMaxError(S1Angle.radians(pow(1e-4, S2Random.randomDouble()) * queryRadius.radians))
                }
                val filter_rect = S2LatLngRect.fromCenterSize(
                        S2LatLng.fromPoint(S2Random.samplePoint(queryCap)),
                        S2LatLng.fromLatLng(kTestCapRadius * S2Random.randomDouble(), kTestCapRadius * S2Random.randomDouble()))

                if (S2Random.oneIn(5)) {
                    query.options().region = filter_rect
                }
                val randomTargetType = S2Random.randomInt(3) // TODO(fmeurisse): test MutableS2ShapeIndex
                when (val targetType = randomTargetType) {
                    0 -> {
                        // Find the points closest to a given point.
                        val point = S2Random.samplePoint(queryCap)
                        val target = S2ClosestPointQueryPointTarget(point)
                        testFindClosestPoints(target, query)
                    }
                    1 -> {
                        // Find the points closest to a given edge.
                        val a = S2Random.samplePoint(queryCap)
                        val b = S2Random.samplePoint(S2Cap.fromCenterAngle(a, queryRadius * 1e-4.pow(S2Random.randomDouble())))
                        val target = S2ClosestPointQueryEdgeTarget(a, b)
                        testFindClosestPoints(target, query)
                    }
                    2 -> {
                        // Find the points closest to a given cell.
                        val minLevel = S2Coords.projection.kMaxDiag.getLevelForMaxValue(queryRadius.radians)
                        val level = minLevel + S2Random.randomInt(S2CellId.kMaxLevel - minLevel + 1)
                        val a = S2Random.samplePoint(queryCap)
                        val cell = S2Cell(S2CellId.fromPoint(a).parent(level))
                        val target = S2ClosestPointQueryCellTarget(cell)
                        testFindClosestPoints(target, query)
                    }
                    else -> {
                        assert(3 == targetType)
                        val targetIndex = MutableS2ShapeIndex()
                        FractalLoopShapeIndexFactory().addEdges(indexCap, 100, targetIndex)
                        val target = S2ClosestPointQueryShapeIndexTarget(targetIndex)
                        target.setIncludeInteriors(S2Random.oneIn(2))
                        testFindClosestPoints(target, query)
                    }
                }
            }
        }

    }


}
