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

import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2PointRegionUnitTest {

    @Test
    fun basic() {
        val p = S2Point(1, 0, 0)
        val r0 = S2PointRegion(p)
        assertThat(r0.point).isEqualTo(p)
        assertThat(r0.contains(p)).isTrue()
        assertThat(r0.contains(r0.point)).isTrue()
        assertThat(r0.contains(S2Point(1, 0, 1))).isFalse()
        val r0_clone = r0.clone()
        assertThat(r0_clone.point).isEqualTo(r0.point)
        assertThat(r0.capBound).isEqualTo(S2Cap.fromPoint(p))
        val ll = S2LatLng.fromPoint(p)
        assertThat(r0.rectBound).isEqualTo(S2LatLngRect(ll, ll))

        // The leaf cell containing a point is still much larger than the point.
        val cell = S2Cell(p)
        assertThat(r0.contains(cell)).isFalse()
        assertThat(r0.mayIntersect(cell)).isTrue()
    }


}
