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

import dilivia.PreConditions.checkState
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.shape.Edge
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import org.apache.commons.math3.util.FastMath.max
import java.util.concurrent.atomic.AtomicReference
import kotlin.reflect.full.createInstance

/**
 * S2ShapeIndex is an abstract base class for indexing polygonal geometry. The objects in the index are known as
 * "shapes", and may consist of points, polylines, and/or polygons, possibly overlapping.  The index makes it very fast
 * to answer queries such as finding nearby shapes, measuring distances, testing for intersection and containment, etc.
 *
 * Each object in the index implements the S2Shape interface. An S2Shape is a collection of edges that optionally
 * defines an interior. The edges do not need to be connected, so for example an S2Shape can represent a polygon
 * with multiple shells and/or holes, or a set of polylines, or a set of points.
 *
 * All geometry within a single S2Shape must have the same dimension, so for example if you want to create an
 * S2ShapeIndex containing a polyline and 10 points, then you will need at least two different S2Shape objects.
 *
 * The most important type of S2ShapeIndex is MutableS2ShapeIndex, which allows you to build an index incrementally by
 * adding or removing shapes.
 *
 * There are a number of built-in classes that work with S2ShapeIndex objects. Generally these classes accept any
 * collection of geometry that can be represented by an S2ShapeIndex, i.e. any combination of points, polylines,
 * and polygons.  Such classes include:
 *
 * - S2ContainsPointQuery: returns the shape(s) that contain a given point.
 * - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge, S2CellId, or S2ShapeIndex.
 * - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
 * - S2BooleanOperation: computes boolean operations such as union, and boolean predicates such as containment.
 * - S2ShapeIndexRegion: computes approximations for a collection of geometry.
 * - S2ShapeIndexBufferedRegion: computes approximations that have been expanded by a given radius.
 *
 * Here is an example showing how to index a set of polygons and then
 * determine which polygon(s) contain each of a set of query points:
 *
 * <pre>
 *   fun testContainment(points: List<S2Point>, polygons: List<S2Polygon>) {
 *     val index = MutableS2ShapeIndex()
 *     for (polygon in polygons) {
 *       index.add(polygon)
 *     }
 *     val query = makeS2ContainsPointQuery(index)
 *     for (point in points) {
 *       for (shape : query.getContainingShapes(point)) {
 *         val polygon = polygons[shape.id()]
 *         ... do something with (point, polygon) ...
 *       }
 *     }
 *   }
 * </pre>
 *
 * This example uses S2Polygon.Shape, which is one example of an S2Shape object. S2Polyline and S2Loop also have
 * nested Shape classes, and there are additional S2Shape types defined in *Shape.kt
 *
 * Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set of shapes that intersect each S2CellId.
 * It is adaptively refined to ensure that no cell contains more than a small number of edges.
 *
 * In addition to implementing a set of abstract methods, all S2ShapeIndex subtypes define an Iterator type with the
 * same API. This makes it easy to convert code that uses a particular S2ShapeIndex subtype to instead use the abstract
 * base class (or vice versa).
 *
 * This class is a port of the S2ShapeIndex class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
abstract class S2ShapeIndex : Iterable<S2Shape?> {

    /**
     * Gets the number of distinct shape ids in the index. This is the same as the number of shapes provided that no
     * shapes have ever been removed. (Shape ids are never reused.)
     *
     * @return The num of used shape ids.
     */
    abstract fun nextNewShapeId(): Int

    /**
     * Gets the shape with the given id, or null if the shape has been removed from the index.
     *
     * @param id A shape id.
     * @return The shape with the given id or null if it has been removed.
     */
    abstract fun shape(id: Int): S2Shape?

    /**
     * Gets an iterator over the indexed shape. Example usage:
     * <pre>
     *     for(shape in index.shapeIterator()) {
     *          ...
     *     }
     * </pre>
     *
     * @return an iterator on the indexed shapes.
     */
    override fun iterator(): Iterator<S2Shape?> = object : Iterator<S2Shape?> {
        private var shapeId: Int = 0
        override fun hasNext(): Boolean = shapeId < nextNewShapeId()
        override fun next(): S2Shape? = shape(shapeId++)
    }

    /**
     * Gets an iterator over the cells of this index. The cells are sorted by ids.
     *
     * @param pos The initial position of the iterator.
     * @return A cell iterator at the desired position.
     */
    fun cellIterator(pos: InitialPosition = InitialPosition.UNPOSITIONED): CellIterator = CellIterator(index = this, pos = pos)

    /**
     * Gets a new iterator positioned as specified.
     *
     * @param pos The initial position of the iterator.
     * @return An iterator instance.
     */
    abstract fun newIterator(pos: InitialPosition): IteratorBase

    /**
     * Prints the content of the index for debug purpose.
     *
     * @return A debug string of the index content.
     */
    fun toDebugString(): String {
        var cellMap = ""
        val iterator = newIterator(pos = InitialPosition.BEGIN)
        while (!iterator.done()) {
            val cell = iterator.cell()
            val id = iterator.id()
            cellMap += "\n$id => ${cell.toDebugString(separator = "", prefix = "\n - ")}"
            iterator.next()
        }
        return """
            |${this.javaClass.simpleName} content
            |---------------------------------------------
            |shapes: ${this.iterator().asSequence().joinToString(separator = ",\n", prefix = "[\n", postfix = "\n]")}
            |CellMap:$cellMap
            |---------------------------------------------
        """.trimMargin()
    }

    /**
     * Gets the maximum dimension of any shape in the index.  Returns -1 if the index does not contain any shapes.
     *
     * Note that the dimension does *not* depend on whether the shapes in the index contain any points; for example,
     * the dimension of an empty point set is 0, and the dimension of an empty polygon is 2.
     *
     * @return The dimension of the index.
     */
    fun getDimension(): Int {
        var dim = -1;
        for (shape in this) {
            if (shape != null) dim = max(dim, shape.dimension)
        }
        return dim
    }

    /**
     * Gets the number of points (objects of dimension zero) in the index.
     * Note that polyline and polygon vertices are *not* included in this count.
     *
     * @return The number of point in the index.
     */
    fun getNumPoints(): Int {
        var count = 0
        for (shape in this) {
            if (shape != null && shape.dimension == 0) {
                count += shape.numEdges
            }
        }
        return count
    }

    /**
     * Computes the total length of all polylines in the index. Returns zero if no polylines are present.
     *
     * All edges are modeled as spherical geodesics.  The result can be converted to a distance on the Earth's
     * surface (with a worst-case error of 0.562% near the equator) using the functions in S2Earth
     *
     * @return The total length of all polylines.
     */
    fun getLength(): S1Angle {
        val length = S1Angle()
        for (shape in this) {
            if (shape != null) length += S2Shape.getLength(shape)
        }
        return length
    }

    /**
     * Computes the total perimeter of all polygons in the index (including both "shells" and "holes"). Returns zero
     * if no polygons are present.
     *
     * All edges are modeled as spherical geodesics. The result can be converted to a distance on the Earth's
     * surface (with a worst-case error of 0.562% near the equator) using the functions in S2Earth.
     *
     * @return The total perimeter of all polygons.
     */
    fun getPerimeter(): S1Angle {
        val perimeter = S1Angle()
        for (shape in this) {
            if (shape != null) perimeter += S2Shape.getPerimeter(shape)
        }
        return perimeter
    }

    /**
     * Computes the total area of all polygons in the index.  Returns zero if no polygons are present.
     * This method has good relative accuracy for both very large and very small regions.
     *
     * Note that the result may exceed 4*Pi if the index contains overlapping polygons.
     *
     * All edges are modeled as spherical geodesics.  The result can be converted to an area on the Earth's surface
     * (with a worst-case error of 0.900% near the poles) using the functions in S2Earth.
     *
     * @return The total area of all polygons.
     */
    fun getArea(): Double {
        var area = 0.0
        for (shape in this) {
            if (shape != null) area += S2Shape.getArea(shape)
        }
        return area
    }

    /**
     * Like getArea(), except that this method is faster and has more error. The additional error is at most
     * 2.22e-15 steradians per vertex, which works out to about 0.09 square meters per vertex on the Earth's
     * surface.  For example, a loop with 100 vertices has a maximum error of about 9 square meters.
     * (The actual error is typically much smaller than this.)
     *
     * @return The total approx area of all polygons.
     */
    fun getApproxArea(): Double {
        var area = 0.0
        for (shape in this) {
            if (shape != null) area += S2Shape.getApproxArea(shape)
        }
        return area
    }

    /**
     * Computes the centroid of all shapes whose dimension is maximal within the index, multiplied by the measure
     * of those shapes.  For example, if the index contains points and polylines, then the result is defined as the
     * centroid of the polylines multiplied by the total length of those polylines.
     * The points would be ignored when computing the centroid.
     *
     * The measure of a given shape is defined as follows:
     *
     *  - For dimension 0 shapes, the measure is shape.numEdges().
     *  - For dimension 1 shapes, the measure is getLength(shape).
     *  - For dimension 2 shapes, the measure is getArea(shape).
     *
     * Note that the centroid is not unit length, so you may need to call normalize() before passing it to other
     * S2 functions.  Note that this function returns (0, 0, 0) if the index contains no geometry.
     *
     * The centroid is scaled by the total measure of the shapes for two reasons:
     *   (1) it is cheaper to compute this way, and
     *   (2) this makes it easier to compute the centroid of a collection of shapes (since the individual centroids
     *       can simply be summed).
     *
     * @return The centroid of the index.
     */
    fun getCentroid(): S2Point {
        val dim = getDimension()
        val centroid = S2Point()
        for (shape in this) {
            if (shape != null && shape.dimension == dim) {
                centroid += S2Shape.getCentroid(shape)
            }
        }
        return centroid
    }


    /**
     * Base class for ShapeIndex dilivia.iterators.
     * Each subtype of S2ShapeIndex should define an Iterator type derived from the following base class.
     *
     * @author Google S2Geometry Project
     * @author Fabien Meurisse (fabien.meurisse@enovea.net)
     * @since 1.0
     */
    abstract class IteratorBase(): Cloneable {

        /** The current sphere cell identifier. */
        protected var id: S2CellId = S2CellId.sentinel
        /** The contents of the current index cell. */
        private val cell: AtomicReference<S2ShapeIndexCell?> = AtomicReference()

        /**
         * Makes a copy of the given source iterator.
         */
        constructor(iter: IteratorBase): this() {
            this.id = iter.id
            this.cell.set(iter.cell())
        }

        /**
         * Gets the S2CellId of the current index cell. If done() is true, returns a value larger than any valid S2CellId
         * (S2CellId.sentinel()).
         *
         * @return The current cell.
         */
        fun id(): S2CellId = id

        /**
         * Gets the center point of the current sphere cell.
         * REQUIRES: !done()
         *
         * @return The center of the current cell.
         */
        fun center(): S2Point {
            checkState { !done() }
            return id().toPoint()
        }

        /**
         * Gets the contents of the current index cell.
         * REQUIRES: !done()
         *
         * @return The current index cell.
         */
        open fun cell(): S2ShapeIndexCell {
            // Like other const methods, this method is thread-safe provided that it
            // does not overlap with calls to non-const methods.
            checkState { !done() }
            var currentCell = rawCell()
            if (currentCell == null) {
                currentCell = getCell()
                this.cell.set(currentCell)
            }
            return currentCell!!
        }


        /**
         * Gets the current contents of the "cell" field, which may be null if the cell contents have not been decoded yet.
         *
         * @return The current cell.
         */
        protected fun rawCell(): S2ShapeIndexCell? = cell.get()

        /**
         * Indicates if the iterator is positioned past the last index cell.
         *
         * @return true if the iterator is on the last index cell.
         */
        fun done(): Boolean = id == S2CellId.sentinel

        /** Positions the iterator at the first index cell (if any). */
        abstract fun begin()

        /** Positions the iterator past the last index cell. */
        abstract fun finish()

        /**
         * Positions the iterator at the next index cell.
         * REQUIRES: !done()
         */
        abstract fun next()

        /**
         * If the iterator is already positioned at the beginning, returns false. Otherwise positions the iterator at the
         * previous entry and returns true.
         */
        abstract fun prev(): Boolean

        /**
         * Positions the iterator at the first cell with id() >= target, or at the end of the index if no such cell exists.
         *
         * @param target The target sphere cell.
         */
        abstract fun seek(target: S2CellId)

        /**
         * Positions the iterator at the cell containing "target".  If no such cell exists, returns false and leaves the
         * iterator positioned arbitrarily.
         * The returned index cell is guaranteed to contain all edges that might intersect the line segment between
         * "target" and the cell center.
         *
         * @param targetPoint The target point.
         * @return true if a cell containing the "target" exists and false otherwise.
         */
        open fun locate(targetPoint: S2Point): Boolean {
            // Let I = cellMap.lowerBound(T), where T is the leaf cell containing "target_point".  Then if T is
            // contained by an index cell, then the containing cell is either I or I'.  We test for containment by comparing
            // the ranges of leaf cells spanned by T, I, and I'.

            val targetCellId = S2CellId.fromPoint(targetPoint)
            seek(targetCellId)
            if (!done() && id().rangeMin() <= targetCellId) return true
            if (prev() && id().rangeMax() >= targetCellId) return true
            return false
        }

        /**
         * Let T be the target S2CellId. If T is contained by some index cell I (including equality), this method positions
         * the iterator at I and returns INDEXED.  Otherwise if T contains one or more (smaller) index cells, it positions
         * the iterator at the first such cell I and returns SUBDIVIDED. Otherwise it returns DISJOINT and leaves the
         * iterator positioned arbitrarily.
         *
         * @param target The target sphere cell.
         */
        open fun locate(target: S2CellId): CellRelation {
            // Let T be the target, let I = cellMap.lowerBound(T.rangeMin()), and  let I' be the predecessor of I.
            // If T contains any index cells, then T contains I. Similarly, if T is contained by an index cell, then the
            // containing cell is either I or I'.  We test for containment by comparing the ranges of leaf cells spanned
            // by T, I, and I'.
            seek(target.rangeMin())
            if (!done()) {
                if (id() >= target && id().rangeMin() <= target) return CellRelation.INDEXED
                if (id() <= target.rangeMax()) return CellRelation.SUBDIVIDED
            }
            if (prev() && id().rangeMax() >= target) return CellRelation.INDEXED
            return CellRelation.DISJOINT;
        }

        /**
         * Sets the iterator state. "cell" typically points to the cell contents, but may also be given as "null" in order
         * to implement decoding on demand.  In that situation, the first that the client attempts to access the cell
         * contents, the GetCell() method is called and "cell_" is updated in a thread-safe way.
         *
         * @param id The current sphere cell.
         * @param cell The current index cell.
         */
        protected fun setState(id: S2CellId, cell: S2ShapeIndexCell?) {
            this.id = id
            this.cell.set(cell)
        }

        /**
         * This method is called to decode the contents of the current cell, if setState() was previously called with a
         * null "cell" argument.  This allows decoding on demand for subtypes that keep the cell contents in an encoded
         * state. It does not need to be implemented at all if setState() is always called with (cell != nullptr).
         *
         * REQUIRES: This method is thread-safe.
         * REQUIRES: Multiple calls to this method return the same value.
         */
        protected open fun getCell(): S2ShapeIndexCell? = null

        /** Sets the iterator state so that done() is true. */
        protected fun setFinished() {
            this.id = S2CellId.sentinel
            this.cell.set(null)
        }

        public override fun clone(): IteratorBase {
            val clone = this::class.createInstance()
            clone.setState(id, cell())
            return clone
        }

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (other !is IteratorBase) return false

            if (id != other.id) return false
            if (cell.get() != other.cell.get()) return false

            return true
        }

        override fun hashCode(): Int {
            var result = id.hashCode()
            result = 31 * result + cell.hashCode()
            return result
        }


    }

    /**
     * A random access iterator that provides low-level access to the cells of the index.
     * Cells are sorted in increasing order of S2CellId.
     *
     * This shape index iterator delegates the iteration action to a concrete iterator base instance associated with the
     * shape index type we are scanning.
     *
     * @constructor Default constructor; must be followed by a call to init().
     *
     * @author Fabien Meurisse (fabien.meurisse@enovea.net)
     * @since 1.0
     */
    class CellIterator(): Cloneable {

        /** The concrete iterator. */
        private lateinit var iter: IteratorBase

        /**
         * Constructs an iterator positioned as specified. By default dilivia.iterators are unpositioned, since this avoids an
         * extra seek in this situation where one of the seek methods (such as locate) is immediately called.
         *
         * If you want to position the iterator at the beginning, e.g. in order to loop through the entire index,
         * do this:
         *
         * <pre>
         *     val iter = S2ShapeIndexCellIterator(index, BEGIN)
         *     while(!iter.done() {
         *          ...
         *          iter.next()
         *     }
         * </pre>
         *
         * @param index The shape index.
         * @param pos The initial position of the iterator (Default = UNPOSITIONED).
         */
        constructor(index: S2ShapeIndex, pos: InitialPosition = InitialPosition.UNPOSITIONED): this() {
            this.iter = index.newIterator(pos)
        }

        /**
         * Makes a copy of the given source iterator.
         *
         * @param cellIter The iterator to copy.
         */
        constructor(cellIter: CellIterator): this() {
            this.iter = cellIter.iter.clone()
        }

        /**
         * Initializes an iterator for the given S2ShapeIndex.  This method may also be called in order to restore an
         * iterator to a valid state after the underlying index has been updated (although it is usually easier just to
         * declare a new iterator whenever required, since iterator construction is cheap).
         *
         * @param index The shape index.
         * @param pos The initial position of the iterator (Default = UNPOSITIONED).
         */
        fun init(index: S2ShapeIndex, pos: InitialPosition = InitialPosition.UNPOSITIONED) {
            this.iter = index.newIterator(pos)
        }

        /**
         * Gets the S2CellId of the current index cell. If done() is true, returns a value larger than any valid S2CellId
         * (S2CellId.sentinel()).
         *
         * @return The current sphere cell id.
         */
        fun id(): S2CellId = iter.id()

        /**
         * Gets the center point of the sphere cell.
         *
         * REQUIRES: !done()
         * @return The center point.
         */
        fun center(): S2Point {
            checkState { !done() }
            return id().toPoint()
        }

        /**
         * Gets the contents of the current index cell.
         *
         * REQUIRES: !done()
         * @return the contents of the current cell.
         */
        fun cell(): S2ShapeIndexCell {
            checkState { !done() }
            return iter.cell()
        }

        /**
         * Indicates if the iterator is positioned past the last index cell.
         *
         * @return true if the iterator is at the end and false otherwise.
         */
        fun done(): Boolean { return iter.done() }

        /** Positions the iterator at the first index cell (if any). */
        fun begin() { iter.begin() }

        /** Positions the iterator past the last index cell. */
        fun finish() { iter.finish() }

        /**
         * Positions the iterator at the next index cell.
         * REQUIRES: !done()
         */
        fun next() { iter.next() }

        /**
         * If the iterator is already positioned at the beginning, returns false.
         * Otherwise positions the iterator at the previous entry and returns true.
         *
         * @return true if the iterator have been moved.
         */
        fun prev(): Boolean = iter.prev()

        /**
         * Positions the iterator at the first cell with id() >= target, or at the end of the index if no such cell exists.
         *
         * @param target The target cell.
         */
        fun seek(target: S2CellId) { iter.seek(target) }

        /**
         * Positions the iterator at the cell containing "target". If no such cell exists, returns false and leaves the
         * iterator positioned arbitrarily. The returned index cell is guaranteed to contain all edges that might intersect
         * the line segment between "target" and the cell center.
         *
         * @param target The target point
         * @return true if the iterator has been moved to target point, false otherwise.
         */
        fun locate(target: S2Point): Boolean = iter.locate(target)

        /**
         * Moves the iterator to a target cell.
         *
         * Let T be the target S2CellId.
         * - If T is contained by some index cell I (including equality), this method positions the iterator at I and
         *   returns INDEXED.
         * - Otherwise if T contains one or more (smaller) index cells, it positions the iterator at the first such cell I
         *   and returns SUBDIVIDED.
         * - Otherwise it returns DISJOINT and leaves the iterator positioned arbitrarily.
         *
         * @param target The target cell.
         * @return The cell relation between the target and the underlying index.
         */
        fun locate(target: S2CellId): CellRelation = iter.locate(target)

        public override fun clone(): CellIterator {
            val cellIter = CellIterator()
            cellIter.iter = iter.clone()
            return cellIter
        }
    }


    /**
     * RangeIterator is a wrapper over S2ShapeIndex.Iterator with extra methods that are useful for merging the contents
     * of two or more S2ShapeIndexes.
     *
     * @author Google S2Geometry Project
     * @author Fabien Meurisse (fabien.meurisse@enovea.net)
     * @since 1.0
     */
    class RangeIterator(index: S2ShapeIndex) {

        private val iterator: CellIterator = index.cellIterator(InitialPosition.BEGIN)

        // The min and max leaf cell ids covered by the current cell.  If done() is
        // true, these methods return a value larger than any valid cell id.
        private var rangeMin: S2CellId = S2CellId.sentinel
        private var rangeMax: S2CellId = S2CellId.sentinel

        init {
            refresh()
        }

        // The current S2CellId and cell contents.
        fun id(): S2CellId = iterator.id()
        fun cell(): S2ShapeIndexCell = iterator.cell()

        // The min and max leaf cell ids covered by the current cell.  If done() is
        // true, these methods return a value larger than any valid cell id.
        fun rangeMin(): S2CellId = rangeMin
        fun rangeMax(): S2CellId = rangeMax

        fun next() {
            iterator.next()
            refresh()
        }

        fun done(): Boolean = iterator.done()

        // Position the iterator at the first cell that overlaps or follows
        // "target", i.e. such that range_max() >= target.range_min().
        fun seekTo(target: RangeIterator) {
            iterator.seek(target.rangeMin)
            // If the current cell does not overlap "target", it is possible that the
            // previous cell is the one we are looking for.  This can only happen when
            // the previous cell contains "target" but has a smaller S2CellId.
            if ((iterator.done() || iterator.id().rangeMin() > target.rangeMax) &&
                iterator.prev() &&
                iterator.id().rangeMax() < target.id()
            ) iterator.next()
            refresh()
        }

        // Position the iterator at the first cell that follows "target", i.e. the
        // first cell such that range_min() > target.range_max().
        fun seekBeyond(target: RangeIterator) {
            iterator.seek(target.rangeMax.next())
            if (!iterator.done() && iterator.id().rangeMin() <= target.rangeMax) {
                iterator.next()
            }
            refresh()
        }

        // Various other convenience methods for the current cell.
        fun clipped(i: Int): S2ClippedShape = cell().clipped(i)
        fun numEdges(i: Int): Int = clipped(i).numEdges
        fun containsCenter(i: Int): Boolean = clipped(i).containsCenter

        // Updates internal state after the iterator has been repositioned.
        private fun refresh() {
            rangeMin = id().rangeMin()
            rangeMax = id().rangeMax()
        }
    }

    // An iterator that advances through all edges in an S2ShapeIndex.
//
// Example usage:
//
// for (EdgeIterator it(index); !it.Done(); it.Next()) {
//   auto edge = it.edge();
//   //...
// }
    class EdgeIterator : Cloneable {

        val index: S2ShapeIndex
        private var shape_id = -1
        private var num_edges = 0
        private var edge_id = -1

        constructor(index: S2ShapeIndex) {
            this.index = index
            next()
        }

        constructor(iterator: EdgeIterator) {
            this.index = iterator.index
            this.shape_id = iterator.shape_id
            this.num_edges = iterator.num_edges
            this.edge_id = iterator.edge_id
        }

        // Returns the current shape id.
        fun shapeId() = shape_id

        // Returns the current edge id.
        fun edgeId() = edge_id

        // Returns the current (shape_id, edge_id).
        fun shapeEdgeId(): ShapeEdgeId = ShapeEdgeId(shape_id, edge_id)

        // Returns the current edge.
        fun edge(): Edge {
            checkState { !done() }
            return index.shape(shape_id)!!.edge(edge_id)
        }

        // Returns true if there are no more edges in the index.
        fun done() = shapeId() >= index.nextNewShapeId()

        // Advances to the next edge.
        fun next(): Unit {
            while (++edge_id >= num_edges) {
                if (++shape_id >= index.nextNewShapeId()) break
                val shape = index.shape(shape_id)
                num_edges = shape?.numEdges ?: 0
                edge_id = -1
            }
        }

        fun debugString(): String = "(shape=$shape_id, edge=$edge_id)"

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (other !is EdgeIterator) return false

            if (index != other.index) return false
            if (shape_id != other.shape_id) return false
            if (edge_id != other.edge_id) return false

            return true
        }

        override fun hashCode(): Int {
            var result = index.hashCode()
            result = 31 * result + shape_id
            result = 31 * result + edge_id
            return result
        }

        public override fun clone(): EdgeIterator {
            return EdgeIterator(this)
        }
    }

}
