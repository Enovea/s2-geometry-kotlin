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

import dilivia.s2.S2CellId
import dilivia.s2.shape.S2Shape
import java.util.*

typealias CellMap = TreeMap<S2CellId, S2ShapeIndexCell>

typealias ShapeIdSet = List<Int>

/**
 * MutableS2ShapeIndex is a class for in-memory indexing of polygonal geometry. The objects in the index are known as
 * "shapes", and may consist of points, polylines, and/or polygons, possibly overlapping. The index makes it very
 * fast to answer queries such as finding nearby shapes, measuring distances, testing for intersection and
 * containment, etc.
 *
 * Internally, MutableS2ShapeIndex is essentially a map from S2CellIds to the set of shapes that intersect each
 * S2CellId.  It is adaptively refined to ensure that no cell contains more than a small number of edges.
 *
 * For efficiency, updates are batched together and applied lazily on the first subsequent query. Locking is used to
 * ensure that MutableS2ShapeIndex is not completely thread-safe. This means that if one thread updates the index,
 * you must ensure that no other thread is reading or updating the index at the same time.
 *
 * @property options The index options.
 * @constructor Create a MutableS2ShapeIndex with the given options. Option values may be changed by calling Init().
 * @param options The options supplied for this index.
 *
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
class MutableS2ShapeIndex(options: Options = Options()) : AbstractMutableS2ShapeIndex(options) {

    override val state: S2ShapeIndexState = DefaultS2ShapeIndexState()

    /** The shapes in the index, accessed by their shape id. Removed shapes are replaced by null value. */
    private val shapes: MutableList<S2Shape?> = mutableListOf()

    /** A map from S2CellId to the set of clipped shapes that intersect that cell. The cell ids cover a set of
     * non-overlapping regions on the sphere. Note that this field is updated lazily (see below). */
    private val cellMap = CellMap()

    val numCells: Int
        get() = cellMap.size

    /**
     * Initialize a MutableS2ShapeIndex with the given options.  This method may only be called when the index is empty
     * (i.e. newly created or clear() has just been called).
     *
     * @param options The index options.
     */
    fun init(options: Options) {
        check(shapes.isEmpty())
        this.options = options
    }

    fun init(options: Options, shapes: List<S2Shape?>, cells: List<Pair<S2CellId, S2ShapeIndexCell>>, state: S2ShapeIndexState) {
        check(this.shapes.isEmpty())
        this.options = options
        this.shapes.addAll(shapes)
        this.cellMap.putAll(cells)
        this.state.pendingAdditionsBegin(state.pendingAdditionsBegin())
        this.state.status(state.status())
        this.state.clearRemovedShapes()
        state.pendingRemovals().forEach { this.state.addRemovedShape(it) }
    }

    /**
     * Gets the number of distinct shape ids that have been assigned. This equals the number of shapes in the index
     * provided that no shapes have ever been removed. (Shape ids are not reused.)
     *
     * @return The number of already assigned shape id.
     */
    override fun nextNewShapeId(): Int = (shapes.lastOrNull()?.id ?: -1) + 1

    /**
     * Gets the shape with the given id, or null if the shape has been removed from the index.
     *
     * @param id A shape id.
     * @return The shape with the given id or null if it has been removed.
     */
    override fun shape(id: Int): S2Shape? = shapes[id]

    override fun addShapeToStore(shape: S2Shape) {
        shape.id = nextNewShapeId()
        shapes.add(shape)
    }

    override fun removeShapeFromStore(shapeId: Int): S2Shape? {
        val shape = shape(shapeId)
        shapes[shapeId] = null
        return shape
    }

    override fun removeAllShapesFromStore() {
        shapes.clear()
        cellMap.clear()
    }

    override fun getShapeIndexCell(id: S2CellId): S2ShapeIndexCell = cellMap.getValue(id)

    override fun setShapeIndexCell(id: S2CellId, cell: S2ShapeIndexCell) {
        cellMap[id] = cell
    }

    override fun removeEmptyShapeIndexCells() {
        cellMap.entries.removeIf { entry -> entry.value.numClipped == 0 }
    }

    override fun removeShapeIndexCell(id: S2CellId) {
        cellMap.remove(id)
    }

    override fun isEmpty(): Boolean = cellMap.isEmpty()

    override fun firstCellId(): S2CellId? = cellMap.firstKey()

    override fun lowerCellId(id: S2CellId): S2CellId? = cellMap.lowerKey(id)

    override fun higherCellId(id: S2CellId): S2CellId? = cellMap.higherKey(id)

    override fun ceilingCellId(id: S2CellId): S2CellId? = cellMap.ceilingKey(id)
    

}
