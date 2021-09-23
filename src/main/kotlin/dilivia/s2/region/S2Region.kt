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
import dilivia.s2.S2Point

/**
 * A S2Region represents a two-dimensional region over the unit sphere. It is an abstract interface with various
 * concrete subtypes.
 *
 * The main purpose of this interface is to allow complex regions to be approximated as simpler regions. So rather than
 * having a wide variety of virtual methods that are implemented by all subtypes, the interface is restricted to
 * methods that are useful for computing approximations.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
interface S2Region : Cloneable {

    /**
     * Bounding spherical cap that contains the region. The bound may not be tight.
     */
    val capBound: S2Cap

    /**
     * Bounding latitude-longitude rectangle that contains the region. The bound may not be tight.
     */
    val rectBound: S2LatLngRect

    /**
     * Returns a small collection of S2CellIds whose union covers the region.
     * The cells are not sorted, may have redundancies (such as cells that contain other cells), and may cover much more
     * area than necessary.
     *
     * This method is not intended for direct use by client code.  Clients should typically use
     * S2RegionCoverer.getCovering, which has options to control the size and accuracy of the covering. Alternatively,
     * if you want a fast covering and don't care about accuracy, consider calling S2RegionCoverer.getFastCovering
     * (which returns a cleaned-up version of the covering computed by this method).
     *
     * getCellUnionBound() implementations should attempt to return a small covering (ideally 4 cells or fewer) that
     * covers the region and can be computed quickly. The result is used by S2RegionCoverer as a starting point for
     * further refinement.
     *
     * Default implementation returns capBound.getCellUnionBound(cellIds)
     *
     * @param cellIds The destination list of the cell union bounds.
     */
    @JvmDefault
    fun getCellUnionBound(cellIds: MutableList<S2CellId>) {
        capBound.getCellUnionBound(cellIds)
    }

    /**
     * Check if this region contains a given S2Cell
     *
     * @param cell A cell.
     * @return true, if the region completely contains the given cell. Otherwise, either the region does not contain
     * the cell or the containment relationship could not be determined.
     */
    operator fun contains(cell: S2Cell): Boolean

    /**
     * Check if this region may intersects a cell.
     *
     * Note that there is currently exactly one implementation of this method (S2LatLngRect.mayIntersect) that takes
     * advantage of the semantics above to be more efficient. For all other S2Region subtypes, this method returns true
     * if the region intersect the cell and false otherwise.
     *
     * @param cell A cell.
     * @return false, if the region does not intersect the given cell. Otherwise, either region intersects the cell, or
     * the intersection relationship could not be determined.
     */
    fun mayIntersect(cell: S2Cell): Boolean

    /**
     * Check if this region contains a given point.
     * The point 'p' is generally required to be unit length, although some subtypes may relax this restriction.
     *
     * @param p A point
     * @return true if and only if the given point is contained by the region.
     */
    @JvmDefault
    fun contains(p: S2Point): Boolean {
        throw NotImplementedError()
    }

    public override fun clone(): S2Region

}
