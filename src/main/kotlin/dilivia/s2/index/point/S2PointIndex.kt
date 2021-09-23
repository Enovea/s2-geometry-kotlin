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

import com.google.common.collect.SortedMultiset
import com.google.common.collect.TreeMultiset
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import mu.KotlinLogging
import java.util.*

/**
 * S2PointIndex maintains an index of points sorted by leaf S2CellId.
 *
 * Each point can optionally store auxiliary data such as an integer or object. This can be used to map results back
 * to client data structures.
 *
 * The class supports adding or removing points dynamically, and provides a seekable iterator interface for navigating
 * the index.
 *
 *  You can use this class in conjuction with S2ClosestPointQuery to find the closest index points to a given query
 *  point.  For example:
 *
 *  <pre>
 * fun test(indexPoints: List<S2Point>, target_points: List<S2Point) {
 *   // The template argument allows auxiliary data to be attached to each
 *   // point (in this case, the array index).
 *   val index = S2PointIndex<Int>()
 *   indexPoints.forEachIndexed { i, p -> index.add(p, i) }
 *
 *   val query = S2ClosestPointQuery<Int>(index)
 *   query.mutable_options()->set_max_results(5);
 *   for (const S2Point& target_point : target_points) {
 *     S2ClosestPointQueryPointTarget target(target_point);
 *     for (const auto& result : query.FindClosestPoints(&target)) {
 *       // The Result class contains the following methods:
 *       //   distance() is the distance to the target.
 *       //   point() is the indexed point.
 *       //   data() is the auxiliary data.
 *       DoSomething(target_point, result);
 *     }
 *   }
 * }
 * </pre>
 *
 * Points can be added or removed from the index at any time by calling add() or remove(). However when the index is
 * modified, you must call init() on each iterator before using it again (or simply create a new iterator).
 *
 * <pre>
 *   index.add(newPoint, 123456)
 *   it.init(index)
 *   it.seek(target.rangeMin())
 * </pre>
 *
 * You can also access the index directly using the iterator interface. For example, here is how to iterate through
 * all the points in a given S2CellId "targetId":
 *
 * <pre>
 *   val it = S2PointIndexIterator<Int>(index)
 *   it.seek(targetId.rangeMin());
 *   while (!it.done() && it.id() <= targetId.rangeMax()) {
 *     doSomething(it.id(), it.point(), it.data())
 *     it.next()
 *   }
 * </pre>
 *
 * TODO(ericv): Consider adding an S2PointIndexRegion class, which could be used to efficiently compute coverings of a collection of S2Points.
 *
 * This class is a port of the S2CellIndex class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @param T The type of indexed data.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2PointIndex<T : Comparable<T>> {

    internal val map: TreeMap<S2CellId, SortedMultiset<PointData<T>>> = TreeMap()

    // Returns the number of points in the index.
    fun numPoints(): Int = map.map { entry -> entry.value.size }.sum()

    // Adds the given point to the index.  Invalidates all dilivia.iterators.
    fun add(point: S2Point, data: T) = add(PointData(point, data))
    fun add(pointData: PointData<T>) {
        val id = S2CellId.fromPoint(pointData.point)
        map.getOrPut(id, { TreeMultiset.create() }).add(pointData)
    }

    // Removes the given point from the index.  Both the "point" and "data"
    // fields must match the point to be removed.  Returns false if the given
    // point was not present.  Invalidates all dilivia.iterators.
    fun remove(point: S2Point, data: T): Boolean = remove(PointData(point, data))
    fun remove(pointData: PointData<T>): Boolean {
        val id = S2CellId.fromPoint(pointData.point)
        val dataSet = map[id] ?: return false
        val removed = dataSet.remove(pointData)
        if (removed && dataSet.isEmpty()) {
            map.remove(id)
        }
        return removed
    }

    // Resets the index to its original empty state.  Invalidates all dilivia.iterators.
    fun clear(): Unit = map.clear()

    fun iterator(): Iterator<T> = Iterator(this)

    override fun toString(): String {
        return map.entries.joinToString("\n") { entry -> "- ${entry.key}: ${entry.value.joinToString(";")}" }
    }

    companion object {
        private val logger = KotlinLogging.logger(S2PointIndex::class.java.name)
    }


    class Iterator<T : Comparable<T>>(index: S2PointIndex<T>) {

        private lateinit var index: S2PointIndex<T>
        private var currentCellId: S2CellId = S2CellId.none
        private var currentPointData: PointData<T>? = null
        private var currentOccurence: Int = 0

        init {
            init(index)
        }

        // Initializes an iterator for the given S2PointIndex.  If the index is
        // non-empty, the iterator is positioned at the first cell.
        //
        // This method may be called multiple times, e.g. to make an iterator
        // valid again after the index is modified.
        fun init(index: S2PointIndex<T>) {
            this.index = index
            begin()
        }

        // The S2CellId for the current index entry.
        // REQUIRES: !done()
        fun id(): S2CellId = currentCellId

        // The point associated with the current index entry.
        // REQUIRES: !done()
        fun point(): S2Point = currentPointData!!.point

        // The client-supplied data associated with the current index entry.
        // REQUIRES: !done()
        fun data(): T = currentPointData!!.data

        // The (S2Point, data) pair associated with the current index entry.
        fun pointData(): PointData<T> = currentPointData!!

        // Returns true if the iterator is positioned past the last index entry.
        fun done(): Boolean = currentPointData == null

        // Positions the iterator at the first index entry (if any).
        fun begin() {
            currentPointData = null
            currentCellId = if (index.map.isNotEmpty()) index.map.firstKey() else S2CellId.sentinel
            while (currentPointData == null && currentCellId != S2CellId.sentinel) {
                currentPointData = index.map[currentCellId]?.firstOrNull()
                if (currentPointData == null) {
                    currentCellId = index.map.higherKey(currentCellId) ?: S2CellId.sentinel
                }
            }
            if (currentPointData != null) currentOccurence = 1

            logger.trace { """
                |Iterator.begin()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
        }

        // Positions the iterator so that done() is true.
        fun finish() {
            currentCellId = S2CellId.sentinel
            currentPointData = null
            currentOccurence = 0

            logger.trace { """
                |Iterator.finish()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
        }

        // Advances the iterator to the next index entry.
        // REQUIRES: !done()
        fun next() {
            var nextPointData: PointData<T>? = null
            var cellId = currentCellId
            var nextOccurence = currentOccurence + 1
            while (nextPointData == null && cellId != S2CellId.sentinel) {
                val pointMultiset = index.map[cellId]
                if (pointMultiset != null) {
                    if (nextOccurence > pointMultiset.count(currentPointData)) {
                        nextPointData = pointMultiset.elementSet()?.higher(currentPointData)
                        nextOccurence = 1
                    } else {
                        nextPointData = currentPointData
                    }
                }
                if (nextPointData == null) {
                    cellId = index.map.higherKey(cellId) ?: S2CellId.sentinel
                    nextPointData = cellId.let { index.map[it]?.firstOrNull() }
                    nextOccurence = if (nextPointData == null) 0 else 1
                }
            }
            currentCellId = cellId
            currentPointData = nextPointData
            currentOccurence = nextOccurence

            logger.trace { """
                |Iterator.next()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }

        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        fun prev(): Boolean {
            var previousPointData: PointData<T>? = null
            var cellId = currentCellId
            var previousOccurence = currentOccurence - 1
            while (previousPointData == null && cellId != S2CellId.none) {
                var pointMultiset = index.map[cellId]
                if (pointMultiset != null) {
                    if (previousOccurence <= 0) {
                        previousPointData = index.map[cellId]?.elementSet()?.lower(currentPointData)
                        previousOccurence = pointMultiset.count(previousPointData)
                    } else {
                        previousPointData = currentPointData
                    }
                }
                if (previousPointData == null) {
                    cellId = index.map.lowerKey(cellId) ?: S2CellId.none
                    pointMultiset = cellId.let { index.map[it] }
                    previousPointData = pointMultiset?.lastOrNull()
                    previousOccurence = pointMultiset?.size ?: 0
                }
            }
            val result = if (previousPointData != null) {
                currentCellId = cellId
                currentPointData = previousPointData
                currentOccurence = previousOccurence
                true
            } else false

            logger.trace { """
                |Iterator.prev()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
                | Has moved: $result
            """.trimMargin() }

            return result
        }

        // Positions the iterator at the first entry with id() >= target, or at the
        // end of the index if no such entry exists.
        fun seek(target: S2CellId) {
            currentCellId = index.map.ceilingKey(target) ?: S2CellId.sentinel
            currentPointData = null
            currentOccurence = 0
            while (currentPointData == null && currentCellId != S2CellId.sentinel) {
                val pointList = index.map.getValue(currentCellId)
                if (pointList.isNotEmpty()) {
                    currentPointData = pointList.first()
                    currentOccurence = 1
                }
                else {
                    currentCellId = index.map.higherKey(currentCellId) ?: S2CellId.sentinel
                }
            }


            logger.trace { """
                |Iterator.seek($target)
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
        }

        companion object {
            private val logger = KotlinLogging.logger(Iterator::class.java.name)
        }

    }

}
