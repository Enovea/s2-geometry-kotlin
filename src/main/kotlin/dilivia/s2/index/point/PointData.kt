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

import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil

// PointData is essentially std::pair with named fields.  It stores an
// S2Point and its associated data, taking advantage of the "empty base
// optimization" to ensure that no extra space is used when Data is empty.
data class PointData<T : Comparable<T>>(val point: S2Point, val data: T) : Comparable<PointData<T>> {

    override fun compareTo(other: PointData<T>): Int {
        val pointComparison = point.compareTo(other.point)
        return if (pointComparison == 0) data.compareTo(other.data) else pointComparison
    }

    override fun toString(): String {
        return "(p=${S2PointUtil.toDegreesString(point)}, d=$data)"
    }


}
