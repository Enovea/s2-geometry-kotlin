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

import dilivia.PreConditions.checkEQ
import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Earth
import dilivia.s2.S2GeometryTestCase.checkDistanceResults
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser.makeIndex
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.S2TextParser.makePolygon
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.shape.S2PointVectorShape
import dilivia.s2.shape.ShapeEdgeId
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import tech.units.indriya.quantity.Quantities
import tech.units.indriya.unit.Units
import java.lang.Math.min
import java.lang.Math.pow
import javax.measure.MetricPrefix

typealias TestingResult = Pair<S2MinDistance, ShapeEdgeId>

class S2ClosestEdgeQueryUnitTest {

    @Test
    fun noEdges() {
        val index = MutableS2ShapeIndex()
        val query = S2ClosestEdgeQuery(index)
        val target = S2ClosestEdgeQuery.PointTarget(S2Point(1, 0, 0))
        val edge = query.findClosestEdge(target)
        assertThat(edge.distance.value).isEqualTo(S1ChordAngle.infinity())
        assertThat(edge.shapeId).isEqualTo(-1)
        assertThat(edge.edgeId).isEqualTo(-1)
        assertThat(edge.isInterior()).isFalse()
        assertThat(edge.isEmpty()).isTrue()
        assertThat(query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())
    }

    @Test
    fun optionsNotModified() {
        // Tests that FindClosestEdge(), GetDistance(), and IsDistanceLess() do not
        // modify query.options, even though all of these methods have their own
        // specific options requirements.
        val options = S2ClosestEdgeQuery.Options()
        options.maxResults = 3
        options.setMaxDistance(S1ChordAngle.degrees(3))
        options.maxError = S1ChordAngle.degrees(0.001)
        val index = makeIndex("1:1 | 1:2 | 1:3 # #")
        val query = S2ClosestEdgeQuery(index, options)
        val target = S2ClosestEdgeQuery.PointTarget(makePoint("2:2"))
        assertThat(query.findClosestEdge(target).edgeId).isEqualTo(1)
        assertThat(query.getDistance(target).degrees()).isCloseTo(1.0, Offset.offset(1e-15))
        assertThat(query.isDistanceLess(target, S1ChordAngle.degrees(1.5))).isTrue()

        // Verify that none of the options above were modified.
        assertThat(query.options.maxResults).isEqualTo(options.maxResults)
        assertThat(query.options.maxDistance).isEqualTo(options.maxDistance)
        assertThat(query.options.maxError).isEqualTo(options.maxError)
    }

    @Test
    fun DistanceEqualToLimit() {
        // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
        // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
        // the distance to the target exactly equals the chosen limit.
        val p0 = makePoint("23:12")
        val p1 = makePoint("47:11")
        val index_points = listOf(p0)
        val index = MutableS2ShapeIndex()
        index.add(S2PointVectorShape(points = index_points))
        val query = S2ClosestEdgeQuery(index)

        // Start with two identical points and a zero distance.
        val target0 = S2ClosestEdgeQuery.PointTarget(p0)
        val dist0 = S1ChordAngle.zero()
        assertThat(query.isDistanceLess(target0, dist0)).isFalse()
        assertThat(query.isDistanceLessOrEqual(target0, dist0)).isTrue()
        assertThat(query.isConservativeDistanceLessOrEqual(target0, dist0)).isTrue()

        // Now try two points separated by a non-zero distance.
        val target1 = S2ClosestEdgeQuery.PointTarget(p1)
        val dist1 = S1ChordAngle(p0, p1)
        assertThat(query.isDistanceLess(target1, dist1)).isFalse()
        assertThat(query.isDistanceLessOrEqual(target1, dist1)).isTrue()
        assertThat(query.isConservativeDistanceLessOrEqual(target1, dist1)).isTrue()
    }

    @Test
    fun TrueDistanceLessThanS1ChordAngleDistance() {
        // Tests that IsConservativeDistanceLessOrEqual returns points where the
        // true distance is slightly less than the one computed by S1ChordAngle.
        //
        // The points below had the worst error from among 100,000 random pairs.
        val p0 = S2Point(0.78516762584829192, -0.50200400690845970, -0.36263449417782678)
        val p1 = S2Point(0.78563011732429433, -0.50187655940493503, -0.36180828883938054)

        // The S1ChordAngle distance is ~4 ulps greater than the true distance.
        val dist1 = S1ChordAngle(p0, p1)
        val limit = dist1.predecessor().predecessor().predecessor().predecessor()
        assertThat(S2Predicates.compareDistance(p0, p1, limit)).isLessThan(0)

        // Verify that IsConservativeDistanceLessOrEqual() still returns "p1".
        val index_points = listOf(p0)
        val index = MutableS2ShapeIndex()
        index.add(S2PointVectorShape(points = index_points))
        val query = S2ClosestEdgeQuery(index)
        val target1 = S2ClosestEdgeQuery.PointTarget(p1)
        assertThat(query.isDistanceLess(target1, limit)).isFalse()
        assertThat(query.isDistanceLessOrEqual(target1, limit)).isFalse()
        assertThat(query.isConservativeDistanceLessOrEqual(target1, limit)).isTrue()
    }

    @Test
    fun TestReuseOfQuery() {
        // Tests that between queries, the internal mechanism for de-duplicating
        // results is re-set.  See b/71646017.
        val index = makeIndex("2:2 # #")
        val query = S2ClosestEdgeQuery(index)
        query.options.setMaxError(S1Angle.degrees(1))
        val target_index = makeIndex("## 0:0, 0:5, 5:5, 5:0")
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(target_index)
        val results1 = query.findClosestEdges(target)
        val results2 = query.findClosestEdges(target)
        assertThat(results2.size).isEqualTo(results1.size)
    }

    @Test
    fun TargetPointInsideIndexedPolygon() {
        // Tests a target point in the interior of an indexed polygon.
        // (The index also includes a polyline loop with no interior.)
        val index = makeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10")
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = true
        options.setMaxDistance(S1Angle.degrees(1))
        val query = S2ClosestEdgeQuery(index, options)
        val target = S2ClosestEdgeQuery.PointTarget(makePoint("2:12"));
        val results = query.findClosestEdges(target)
        assertThat(results.size).isEqualTo(1)
        assertThat(results[0].distance.value).isEqualTo(S1ChordAngle.zero())
        assertThat(results[0].shapeId).isEqualTo(1)
        assertThat(results[0].edgeId).isEqualTo(-1)
        assertThat(results[0].isInterior()).isTrue()
        assertThat(results[0].isEmpty()).isFalse()
    }

    @Test
    fun TargetPointOutsideIndexedPolygon() {
        // Tests a target point in the interior of a polyline loop with no
        // interior.  (The index also includes a nearby polygon.)
        val index = makeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10")
        val options = S2ClosestEdgeQuery.Options()
        options.includeInteriors = true
        options.setMaxDistance(S1Angle.degrees(1));
        val query = S2ClosestEdgeQuery(index, options)
        val target = S2ClosestEdgeQuery.PointTarget(makePoint("2:2"));
        val results = query.findClosestEdges(target)
        assertThat(results.size).isEqualTo(0)
    }

    @Test
    fun TargetPolygonContainingIndexedPoints() {
        // Two points are contained within a polyline loop (no interior) and two
        // points are contained within a polygon.
        val index = makeIndex("2:2 | 3:3 | 1:11 | 3:13 # #")
        val query = S2ClosestEdgeQuery(index)
        query.options.setMaxDistance(S1Angle.degrees(1))
        val target_index = makeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10")
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(target_index)
        target.setIncludeInteriors(true)
        val results = query.findClosestEdges(target)
        assertThat(results.size).isEqualTo(2)
        assertThat(results[0].distance.value).isEqualTo(S1ChordAngle.zero())
        assertThat(results[0].shapeId).isEqualTo(0)
        assertThat(results[0].edgeId).isEqualTo(2)  // 1:11
        assertThat(results[1].distance.value).isEqualTo(S1ChordAngle.zero())
        assertThat(results[1].shapeId).isEqualTo(0)
        assertThat(results[1].edgeId).isEqualTo(3)  // 3:13
    }

    @Test
    fun EmptyTargetOptimized() {
        // Ensure that the optimized algorithm handles empty targets when a distance
        // limit is specified.
        val index = MutableS2ShapeIndex()
        index.add(S2Polygon.Shape(polygon = S2Polygon(S2Loop.makeRegularLoop(S2Point(1, 0, 0), S1Angle.radians(0.1), 1000))));
        val query = S2ClosestEdgeQuery(index)
        query.options.setMaxDistance(S1Angle.radians(1e-5))
        val target_index = MutableS2ShapeIndex()
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(target_index)
        assertThat(query.findClosestEdges(target).size).isEqualTo(0)
    }

    @Test
    fun EmptyPolygonTarget() {
        // Verifies that distances are measured correctly to empty polygon targets.
        val empty_polygon_index = makeIndex("# # empty")
        val point_index = makeIndex("1:1 # #")
        val full_polygon_index = makeIndex("# # full")
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(empty_polygon_index)
        target.setIncludeInteriors(true)

        val empty_query = S2ClosestEdgeQuery(empty_polygon_index)
        empty_query.options.includeInteriors = true
        assertThat(empty_query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())

        val point_query = S2ClosestEdgeQuery(point_index)
        point_query.options.includeInteriors = true;
        assertThat(point_query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())

        val full_query = S2ClosestEdgeQuery(full_polygon_index)
        full_query.options.includeInteriors = true;
        assertThat(full_query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())
    }

    @Test
    fun FullLaxPolygonTarget() {
        // Verifies that distances are measured correctly to full LaxPolygon targets.
        val empty_polygon_index = makeIndex("# # empty")
        val point_index = makeIndex("1:1 # #")
        val full_polygon_index = makeIndex("# # full")
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(full_polygon_index)
        target.setIncludeInteriors(true)

        val empty_query = S2ClosestEdgeQuery(empty_polygon_index)
        empty_query.options.includeInteriors = true
        assertThat(empty_query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())

        val point_query = S2ClosestEdgeQuery(point_index)
        point_query.options.includeInteriors = true;
        assertThat(point_query.getDistance(target)).isEqualTo(S1ChordAngle.zero())

        println()
        println()
        println("====================================================")
        println(full_polygon_index)
        println("====================================================")
        println()
        println()

        val full_query = S2ClosestEdgeQuery(full_polygon_index)
        full_query.options.includeInteriors = true;
        assertThat(full_query.getDistance(target)).isEqualTo(S1ChordAngle.zero())
    }

    @Test
    fun fullS2PolygonTarget() {
        // Verifies that distances are measured correctly to full S2Polygon targets
        // (which use a different representation of "full" than LaxPolygon does).
        val empty_polygon_index = makeIndex("# # empty");
        val point_index = makeIndex("1:1 # #");
        val full_polygon_index = makeIndex("# #");
        full_polygon_index.add(S2Polygon.Shape(polygon = makePolygon("full")))

        val target = S2ClosestEdgeQuery.ShapeIndexTarget(full_polygon_index);
        target.setIncludeInteriors(true)

        val empty_query = S2ClosestEdgeQuery(empty_polygon_index)
        empty_query.options.includeInteriors = true;
        assertThat(empty_query.getDistance(target)).isEqualTo(S1ChordAngle.infinity())

        val point_query = S2ClosestEdgeQuery(point_index)
        point_query.options.includeInteriors = true;
        assertThat(point_query.getDistance(target)).isEqualTo(S1ChordAngle.zero())

        val full_query = S2ClosestEdgeQuery(full_polygon_index);
        full_query.options.includeInteriors = true;
        assertThat(full_query.getDistance(target)).isEqualTo(S1ChordAngle.zero())
    }

    @Test
    fun isConservativeDistanceLessOrEqual() {
        // Test
        var num_tested = 0
        var num_conservative_needed = 0
        repeat(1000) { iter ->
            S2Random.reset(iter + 1)  // Easier to reproduce a specific case.
            val x = S2Random.randomPoint()
            val dir = S2Random.randomPoint()
            val r = S1Angle.radians(M_PI * pow(1e-30, S2Random.randomDouble()))
            val y = S2EdgeDistances.interpolateAtDistance(r, x, dir)
            val limit = S1ChordAngle(r)
            if (S2Predicates.compareDistance(x, y, limit) <= 0) {
                val index = MutableS2ShapeIndex()
                index.add(S2PointVectorShape(points = listOf<S2Point>(x)))
                val query = S2ClosestEdgeQuery(index)
                val target = S2ClosestEdgeQuery.PointTarget(y);
                assertThat(query.isConservativeDistanceLessOrEqual(target, limit)).isTrue()
                ++num_tested
                if (!query.isDistanceLess(target, limit)) ++num_conservative_needed;
            }
        }

        println("Num tested = $num_tested, Num conservative = $num_conservative_needed")

        // Verify that in most test cases, the distance between the target points
        // was close to the desired value.  Also verify that at least in some test
        // cases, the conservative distance test was actually necessary.
        assertThat(num_tested).isGreaterThanOrEqualTo(300)
        assertThat(num_tested).isLessThanOrEqualTo(700)
        assertThat(num_conservative_needed).isGreaterThanOrEqualTo(20)
    }

    @Test
    fun circleEdges() {
        testWithIndexFactory(
            RegularLoopShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    @Test
    fun fractalEdges() {
        testWithIndexFactory(
            FractalLoopShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    @Test
    fun pointCloudEdges() {
        testWithIndexFactory(
            PointCloudShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    @Test
    fun conservativeCellDistanceIsUsed() {
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2ShapeIndex cells.
        for (seed in listOf(42, 681, 894, 1018, 1750, 1759, 2401)) {
            testWithIndexFactory(
                FractalLoopShapeIndexFactory(),
                5, 100, 10
            );
        }
    }

    companion object {

        val kNumIndexes = 50
        val kNumEdges = 100
        val kNumQueries = 200

        // The approximate radius of S2Cap from which query edges are chosen.
        val kTestCapRadius: S1Angle = S2Earth.toAngle(Quantities.getQuantity(10, MetricPrefix.KILO(Units.METRE)))

        // An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
        val kTestChordAngleError = 1e-15


        // Converts to the format required by CheckDistanceResults() in s2testing.h.
        fun convertResults(results: List<S2ClosestEdgeQueryResult<S2MinDistance>>): List<TestingResult> {
            val testing_results = mutableListOf<TestingResult>()
            for (result in results) {
                testing_results.add(
                    TestingResult(
                        result.distance,
                        ShapeEdgeId(result.shapeId, result.edgeId)
                    )
                )
            }
            return testing_results
        }

        // Use "query" to find the closest edge(s) to the given target.  Verify that
// the results satisfy the search criteria.
        fun getClosestEdges(
            target: S2MinDistanceTarget,
            query: S2ClosestEdgeQuery,
            edges: MutableList<S2ClosestEdgeQueryResult<S2MinDistance>>
        ) {
            query.findClosestEdges(target, edges)
            assertThat(edges.size).isLessThanOrEqualTo(query.options.maxResults)
            if (query.options.maxDistance == S2MinDistanceFactory.infinity()) {
                val min_expected = min(query.options.maxResults, S2CountEdges.countEdges(query.index()))
                if (!query.options.includeInteriors) {
                    // We can predict exactly how many edges should be returned.
                    assertThat(edges.size).isEqualTo(min_expected)
                } else {
                    // All edges should be returned, and possibly some shape interiors.
                    assertThat(edges.size).isGreaterThanOrEqualTo(min_expected)
                }
            }
            for (edge in edges) {
                // Check that the edge satisfies the max_distance() condition.
                assertThat(edge.distance).isLessThan(query.options.maxDistance)
            }
        }

        fun testFindClosestEdges(
            target: S2MinDistanceTarget,
            query: S2ClosestEdgeQuery
        ): S2ClosestEdgeQueryResult<S2MinDistance> {
            val expected = mutableListOf<S2ClosestEdgeQueryResult<S2MinDistance>>()
            val actual = mutableListOf<S2ClosestEdgeQueryResult<S2MinDistance>>()
            query.options.useBruteForce = true
            getClosestEdges(target, query, expected)
            query.options.useBruteForce = false
            getClosestEdges(target, query, actual)
            assertThat(
                checkDistanceResults(
                    S2MinDistanceFactory,
                    convertResults(expected),
                    convertResults(actual),
                    query.options.maxResults,
                    query.options.maxDistance,
                    query.options.maxError
                )
            )
                .withFailMessage(
                    """
                | max_results=${query.options.maxResults}
                | max_distance=${query.options.maxDistance}
                | max_error=${query.options.maxError}
              """.trimMargin()
                )
                .isTrue()

            if (expected.isEmpty()) return S2ClosestEdgeQueryResult(S2MinDistanceFactory.infinity())

            // Note that when options.maxError > 0, expected[0].distance() may not
            // be the minimum distance.  It is never larger by more than max_error(),
            // but the actual value also depends on max_results().
            //
            // Here we verify that GetDistance() and IsDistanceLess() return results
            // that are consistent with the max_error() setting.
            val max_error = query.options.maxError
            val min_distance = expected[0].distance
            assertThat(query.getDistance(target)).isLessThanOrEqualTo(min_distance.value + max_error)

            // Test IsDistanceLess().
            assertThat(query.isDistanceLess(target, min_distance.value - max_error)).isFalse()
            assertThat(query.isConservativeDistanceLessOrEqual(target, min_distance.value)).isTrue()

            // Return the closest edge result so that we can also test Project.
            return expected[0]
        }

        // The running time of this test is proportional to
//    (num_indexes + num_queries) * num_edges.
// (Note that every query is checked using the brute force algorithm.)
        fun testWithIndexFactory(
            factory: ShapeIndexFactory,
            num_indexes: Int,
            num_edges: Int,
            num_queries: Int,
            seed: Int = 1
        ) {
            // Build a set of MutableS2ShapeIndexes containing the desired geometry.
            val index_caps = mutableListOf<S2Cap>()
            val indexes = mutableListOf<MutableS2ShapeIndex>()
            repeat(num_indexes) { i ->
                S2Random.reset(seed + i)
                index_caps.add(S2Cap.fromCenterAngle(S2Random.randomPoint(), kTestCapRadius))
                indexes.add(MutableS2ShapeIndex())
                factory.addEdges(index_caps.last(), num_edges, indexes.last())
            }
            repeat(num_queries) { i ->
                S2Random.reset(seed + i)
                val i_index = S2Random.randomInt(num_indexes)
                val index_cap = index_caps[i_index]

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                val query_radius = 2 * index_cap.radius()
                val query_cap = S2Cap.fromCenterAngle(index_cap.center, query_radius)
                val query = S2ClosestEdgeQuery(indexes[i_index])

                // Occasionally we don't set any limit on the number of result edges.
                // (This may return all edges if we also don't set a distance limit.)
                if (!S2Random.oneIn(5)) {
                    query.options.maxResults = 1 + S2Random.randomInt(10)
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Random.oneIn(3)) {
                    query.options.setMaxDistance(S2Random.randomDouble() * query_radius)
                }
                if (S2Random.oneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.options.setMaxError(
                        S1Angle.radians(
                            pow(
                                1e-4,
                                S2Random.randomDouble()
                            ) * query_radius.radians
                        )
                    )
                }
                query.options.includeInteriors = S2Random.oneIn(2)
                val target_type = S2Random.randomInt(4)

                if (target_type == 0) {
                    // Find the edges closest to a given point.
                    val point = S2Random.samplePoint(query_cap)
                    val target = S2ClosestEdgeQuery.PointTarget(point)
                    val closest = testFindClosestEdges(target, query)
                    if (!closest.distance.value.isInfinity()) {
                        // Also test the Project method.
                        assertThat(
                            closest.distance.value.toAngle().radians
                        )
                            .isCloseTo(
                                S1Angle(point, query.project(point, closest)).radians,
                                Offset.offset(kTestChordAngleError)
                            )
                    }
                } else if (target_type == 1) {
                    // Find the edges closest to a given edge.
                    val a = S2Random.samplePoint(query_cap)
                    val b = S2Random.samplePoint(
                        S2Cap.fromCenterAngle(
                            a,
                            pow(1e-4, S2Random.randomDouble()) * query_radius
                        )
                    )
                    val target = S2ClosestEdgeQuery.EdgeTarget(a, b)
                    testFindClosestEdges(target, query)
                } else if (target_type == 2) {
                    // Find the edges closest to a given cell.
                    val min_level = S2Coords.projection.kMaxDiag.getLevelForMaxValue(query_radius.radians)
                    val level = min_level + S2Random.randomInt(S2CellId.kMaxLevel - min_level + 1)
                    val a = S2Random.samplePoint(query_cap)
                    val cell = S2Cell(S2CellId.fromPoint(a).parent(level))
                    val target = S2ClosestEdgeQuery.CellTarget(cell)
                    testFindClosestEdges(target, query)
                } else {
                    checkEQ(3, target_type)
                    // Use another one of the pre-built indexes as the target.
                    val j_index = S2Random.randomInt(num_indexes)
                    val target = S2ClosestEdgeQuery.ShapeIndexTarget(indexes[j_index]);
                    target.setIncludeInteriors(S2Random.oneIn(2))
                    testFindClosestEdges(target, query)
                }
            }
        }

    }

//    @Test
//    fun test() {
//        val factory = PointCloudShapeIndexFactory()
//        S2Random.reset(0)
//        val cap = S2Cap.fromCenterAngle(S2Random.randomPoint(), kTestCapRadius)
//        val index = MutableS2ShapeIndex()
//        factory.addEdges(cap, 30, index)
//        val options = S2ClosestEdgeQuery.Options()
//        options.maxResults = 3
//        options.maxDistance = S2MinDistance(S1ChordAngle.infinity())
//        options.maxError = S1ChordAngle.zero()
//        options.includeInteriors = true
//        options.useBruteForce = false
//    }

}
