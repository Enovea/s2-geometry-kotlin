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

import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng
import dilivia.s2.S2TextParser
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test


class S2RegionUnionUnitTest  {

    @Test
    fun basic() {
        val ruEmpty = S2RegionUnion()
        assertThat(ruEmpty.numRegions()).isEqualTo(0)
        assertThat(ruEmpty.capBound).isEqualTo(S2Cap.empty)
        assertThat(ruEmpty.rectBound).isEqualTo(S2LatLngRect.empty())
        val emptyClone = ruEmpty.clone()

        val twoPointRegion = listOf<S2Region>(
                S2PointRegion(S2LatLng.fromDegrees(35, 40).toPoint()),
                S2PointRegion(S2LatLng.fromDegrees(-35, -40).toPoint())
        )

        val twoPointsOrig = S2RegionUnion(twoPointRegion)
        // two_point_region is in a valid, but unspecified, state.

        val twoPoints = twoPointsOrig.clone()
        // The bounds below may not be exactly equal because the S2PointRegion
        // version converts each S2LatLng value to an S2Point and back.
        assertThat(S2TextParser.makeLatLngRect("-35:-40,35:40").approxEquals(twoPoints.rectBound)).isTrue()

        val face0 = S2Cell.fromFace(0)
        assertThat(twoPoints.mayIntersect(face0)).isTrue()
        assertThat(twoPoints.contains(face0)).isFalse()

        assertThat(twoPoints.contains(S2LatLng.fromDegrees(35, 40).toPoint())).isTrue()
        assertThat(twoPoints.contains(S2LatLng.fromDegrees(-35, -40).toPoint())).isTrue()
        assertThat(twoPoints.contains(S2LatLng.fromDegrees(0, 0).toPoint())).isFalse()

        // Check that we can Add() another region.
        val threePoints = twoPoints.clone()
        assertThat(threePoints.contains(S2LatLng.fromDegrees(10, 10).toPoint())).isFalse()
        threePoints.add(S2PointRegion(S2LatLng.fromDegrees(10, 10).toPoint()))
        assertThat(threePoints.contains(S2LatLng.fromDegrees(10, 10).toPoint())).isTrue()

        val coverer = S2RegionCoverer(maxCells = 1)
        val covering = mutableListOf<S2CellId>()
        coverer.getCovering(twoPoints, covering)
        assertThat(covering.size).isEqualTo(1)
        assertThat(covering[0]).isEqualTo(face0.id())
    }

}  
