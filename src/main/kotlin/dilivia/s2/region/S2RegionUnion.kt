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

// An S2RegionUnion represents a union of possibly overlapping regions.
// It is convenient for computing a covering of a set of regions.
// Create an empty region.  Can be made non-empty by calling Init() or Add().
class S2RegionUnion() : S2Region {

    private val regions: MutableList<S2Region> = mutableListOf()

    // Create a region representing the union of the given regions.
    constructor(regions: List<S2Region>) : this() {
        init(regions)
    }

    constructor(region: S2RegionUnion) : this(region.regions)


    // Initialize region by taking ownership of the given regions.
    fun init(regions: List<S2Region>) {
        assert(this.regions.isEmpty())
        this.regions.addAll(regions)
    }

    // Releases ownership of the regions of this union and returns them,
    // leaving this region empty.
    fun clear() {
        regions.clear()
    }

    // Add the given region to the union.  This method can be called repeatedly
    // as an alternative to Init().
    fun add(region: S2Region) {
        this.regions.add(region)
    }

    // Accessor methods.
    fun numRegions(): Int = regions.size
    fun region(i: Int): S2Region = regions[i]

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    // TODO(ericv): This could be optimized to return a tighter bound,
    // but doesn't seem worth it unless profiling shows otherwise.
    override val capBound: S2Cap
        get() = rectBound.capBound

    override val rectBound: S2LatLngRect
        get() {
            var result = S2LatLngRect.empty()
            for (i in 0 until numRegions()) {
                result = result.union(region(i).rectBound)
            }
            return result
        }

    public override fun clone(): S2RegionUnion = S2RegionUnion(this)

    // Note that this method is allowed to return false even if the cell
    // is contained by the region.
    override fun contains(cell: S2Cell): Boolean = regions.any { r -> r.contains(cell) }

    override fun mayIntersect(cell: S2Cell): Boolean = regions.any { r -> r.mayIntersect(cell) }

    override fun contains(p: S2Point): Boolean = regions.any { r -> r.contains(p) }

}

