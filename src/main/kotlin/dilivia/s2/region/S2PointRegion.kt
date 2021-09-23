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

import dilivia.PreConditions.requireArgument
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil

/**
// An S2PointRegion is a region that contains a single point.  It is more
// expensive than the raw S2Point type and is useful mainly for completeness.
// Create a region containing the given point, which must be unit length.
 */
class S2PointRegion(val point: S2Point) : S2Region {

    init {
        requireArgument { S2PointUtil.isUnitLength(point) }
    }

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public override fun clone(): S2PointRegion = this

    override val capBound: S2Cap = S2Cap.fromPoint(point)

    override val rectBound: S2LatLngRect = S2LatLngRect.fromPoint(S2LatLng.fromPoint(point))

    override fun contains(cell: S2Cell): Boolean = false

    override fun mayIntersect(cell: S2Cell): Boolean = cell.contains(point)

    override fun contains(p: S2Point): Boolean = point == p

}
