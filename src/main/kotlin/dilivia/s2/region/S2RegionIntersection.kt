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

import dilivia.s2.S2Point

// An S2RegionIntersection represents the intersection of a set of regions.
// It is convenient for computing a covering of the intersection of a set of
// regions.
// Creates an empty intersection that should be initialized by calling Init().
// Note: an intersection of no regions covers the entire sphere.
class S2RegionIntersection() : S2Region {

    private val regions: MutableList<S2Region> = mutableListOf()

    // Create a region representing the intersection of the given regions.
    constructor(regions: List<S2Region>): this() {
        init(regions)
    }

    // Internal copy constructor used only by Clone() that makes a deep copy of
    // its argument.
    constructor(intersection: S2RegionIntersection): this() {
        init(intersection.regions)
    }

    // Initialize region by taking ownership of the given regions.
    fun init(regions: List<S2Region>) {
        assert(regions.isEmpty())
        this.regions.addAll(regions)
    }

    // Releases ownership of the regions of this intersection and returns them,
    // leaving this region empty.
    fun clear() {
        regions.clear()
    }

    // Accessor methods.
    fun numRegions(): Int = regions.size
    fun region(i: Int): S2Region = regions[i]

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    override fun clone(): S2RegionIntersection {
        return S2RegionIntersection(this)
    }

    override val capBound: S2Cap
        get() {
            // TODO(ericv): This could be optimized to return a tighter bound, but
            // doesn't seem worth it unless profiling shows otherwise.
            return rectBound.capBound
        }

    override val rectBound: S2LatLngRect
        get() {
            var result = S2LatLngRect.full()
            for (i in 0 until numRegions()) {
                result = result.intersection(region(i).rectBound)
            }
            return result
        }

    override fun contains(p: S2Point): Boolean {
        for (i in 0 until numRegions()) {
            if (!region(i).contains(p)) return false
        }
        return true
    }

    override fun contains(cell: S2Cell): Boolean {
        for (i in 0 until numRegions()) {
            if (!region(i).contains(cell)) return false
        }
        return true
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        for (i in 0 until numRegions()) {
            if (!region(i).mayIntersect(cell)) return false
        }
        return true
    }

}
