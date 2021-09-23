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
package dilivia.s2

import dilivia.collections.isSorted
import dilivia.s2.index.Delta
import dilivia.s2.index.Distance
import dilivia.s2.index.DistanceFactory
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Region
import org.assertj.core.api.Assertions.assertThat

object S2GeometryTestCase {

    fun checkCovering(region: S2Region, covering: S2CellUnion, check_tight: Boolean, id: S2CellId = S2CellId()) {
        if (!id.isValid) {
            for (face in 0..5) {
                checkCovering(region, covering, check_tight, S2CellId.fromFace(face))
            }
            return
        }

        if (!region.mayIntersect(S2Cell(id))) {
            // If region does not intersect id, then neither should the covering.
            if (check_tight) assertThat(!covering.intersects(id)).isTrue()

        } else if (!covering.contains(id)) {
            // The region may intersect id, but we can't assert that the covering
            // intersects id because we may discover that the region does not actually
            // intersect upon further subdivision.  (MayIntersect is not exact.)
            assertThat(!region.contains(S2Cell(id))).isTrue()
            assertThat(!id.isLeaf).isTrue()
            val end = id.childEnd()
            var child = id.childBegin()
            while (child != end) {
                checkCovering(region, covering, check_tight, child)
                child = child.next()
            }
        }
    }

    // Compare two sets of "closest" items, where "expected" is computed via brute
    // force (i.e., considering every possible candidate) and "actual" is computed
    // using a spatial data structure.  Here "max_size" is a bound on the maximum
    // number of items, "max_distance" is a limit on the distance to any item, and
    // "max_error" is the maximum error allowed when selecting which items are
    // closest (see S2ClosestEdgeQuery::Options::max_error).
    fun <T : Distance<T>, D : Comparable<D>> checkDistanceResults(
        distanceFactory: DistanceFactory<T>,
        expected: List<Pair<T, D>>,
        actual: List<Pair<T, D>>,
        max_size: Int,
        max_distance: T,
        max_error: Delta
    ): Boolean {
        // This is a conservative bound on the error in computing the distance from
        // the target geometry to an S2Cell.  Such errors can cause candidates to be
        // pruned from the result set even though they may be slightly closer.
        val kMaxPruningError: Delta = S1ChordAngle.radians(1e-15)
        return (checkResultSet(
            distanceFactory,
            actual,
            expected,
            max_size,
            max_distance,
            max_error,
            kMaxPruningError,
            "Missing"
        ) and /*not &&*/
                checkResultSet(
                    distanceFactory,
                    expected,
                    actual,
                    max_size,
                    max_distance,
                    max_error,
                    Delta.zero(),
                    "Extra"
                ))
    }


    // Check that result set "x" contains all the expected results from "y", and
    // does not include any duplicate results.
    private fun <T : Distance<T>, D : Comparable<D>> checkResultSet(
        distanceFactory: DistanceFactory<T>,
        x: List<Pair<T, D>>,
        y: List<Pair<T, D>>,
        max_size: Int,
        max_distance: T,
        max_error: Delta,
        max_pruning_error: Delta,
        label: String
    ): Boolean {
        // Results should be sorted by distance, but not necessarily then by Id.
        assertThat(x.isSorted { p1, p2 -> p1.first.compareTo(p2.first) }).isTrue()

        // Result set X should contain all the items from Y whose distance is less
        // than "limit" computed below.
        var limit = distanceFactory.zero()
        if (x.size < max_size) {
            // Result set X was not limited by "max_size", so it should contain all
            // the items up to "max_distance", except that a few items right near the
            // distance limit may be missed because the distance measurements used for
            // pruning S2Cells are not conservative.
            if (max_distance == distanceFactory.infinity()) {
                limit = max_distance;
            } else {
                limit = max_distance - max_pruning_error;
            }
        } else if (!x.isEmpty()) {
            // Result set X contains only the closest "max_size" items, to within a
            // tolerance of "max_error + max_pruning_error".
            limit = (x.last().first - max_error) - max_pruning_error;
        }

        var result = true
        for (yp in y) {
            // Note that this test also catches duplicate values.
            val count = x.count { xp -> xp.second == yp.second }
            if (yp.first < limit && count != 1) {
                result = false
                println((if (count > 1) "Duplicate" else label) + " distance = ${yp.first}, id = ${yp.second}\n")
            }
        }

        return result
    }


    fun KmToAngle(km: Double): S1Angle = S1Angle.radians(km / S2Earth.radiusKm)

}
