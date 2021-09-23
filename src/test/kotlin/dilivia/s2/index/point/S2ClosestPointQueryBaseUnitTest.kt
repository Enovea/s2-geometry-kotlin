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

import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.index.S2MaxDistance
import dilivia.s2.index.S2MaxDistanceFactory
import dilivia.s2.index.S2MaxDistancePointTarget
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test

//
// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_point_query_test.cc.


// This is a proof-of-concept prototype of a possible S2FurthestPointQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
typealias FurthestPointQuery = S2ClosestPointQueryBase<S2MaxDistance, Int>

class FurthestPointTarget(point: S2Point) : S2MaxDistancePointTarget(point) {

    override fun maxBruteForceIndexSize(): Int = 10

}

class S2ClosestPointQueryBaseUnitTest {

    @Test
    fun closestPointQueryBaseMaxDistance() {
        val index = S2PointIndex<Int>()
      val points = S2TextParser.parsePoints("0:0, 1:0, 2:0, 3:0")
        for (i in points.indices) {
            index.add(points[i], i);
        }
        val query = FurthestPointQuery(S2MaxDistanceFactory, index)
        val options = S2ClosestPointQueryBase.Options(S2MaxDistanceFactory)
        options.setMaxResult(1)
        val target = FurthestPointTarget(S2TextParser.makePoint("4:0"))
        val results = query . findClosestPoints (target, options)
        assertThat(results.size).isEqualTo(1)
        assertThat(results[0].point()).isEqualTo(points[0])
        assertThat(results[0].data()).isEqualTo(0)
        assertThat(S1ChordAngle(results[0].distance.value).toAngle().degrees()).isCloseTo(4.0, Offset.offset(1e-13))
    }

}
