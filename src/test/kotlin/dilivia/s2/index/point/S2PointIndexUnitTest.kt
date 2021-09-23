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
import dilivia.s2.S2Random
import dilivia.s2.region.S2CellUnion
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test


class S2PointIndexUnitTest {

    private val logger = KotlinLogging.logger { }

    private val index: S2PointIndex<Int> = S2PointIndex()
    private val contents: SortedMultiset<PointData<Int>> = TreeMultiset.create()


    fun add(point: S2Point, data: Int) {
        index.add(point, data)
        contents.add(PointData(point, data))
    }

    fun remove(point: S2Point, data: Int): Boolean {
        // If there are multiple copies, remove only one.
        contents.remove(PointData(point, data))
        return index.remove(point, data)  // Invalidates "point".
    }

    fun verify() {
        verifyContents()
        verifyIteratorMethods()
    }

    private fun verifyContents() {
        logger.debug { "Contents: ${contents.size} elements" }
        logger.debug { "Index: ${index.numPoints()} elements" }
        val remaining = TreeMultiset.create(contents)
        val iter = S2PointIndex.Iterator(index)
        var i = 0
        while (!iter.done()) {
            val currentPointData = iter.pointData()
            logger.debug { "Iteration $i" }
            logger.debug { "Current point : $currentPointData" }
            logger.debug { "Point remaining occurences: ${remaining.count(currentPointData)}" }
            val removed = remaining.remove(currentPointData)
            logger.debug { "Removed from remaining: $removed" }
            logger.debug { "Point remaining occurences after removal: ${remaining.count(currentPointData)}" }
            logger.debug { "Remaining: ${remaining.size}" }
            assertThat(removed).isTrue
            iter.next()
            ++i
        }
        assertThat(remaining.isEmpty()).isTrue
    }

    private fun verifyIteratorMethods() {
        val iter = index.iterator()
        assertThat(iter.prev()).isFalse
        iter.finish()
        assertThat(iter.done()).isTrue

        // Iterate through all the cells in the index.
        var prev_cellid = S2CellId.none
        var min_cellid = S2CellId.begin(S2CellId.kMaxLevel)
        iter.begin()
        while (!iter.done()) {
            val cellid = iter.id()
            assertThat(S2CellId.fromPoint(iter.point())).isEqualTo(cellid)
            assertThat(cellid >= prev_cellid).isTrue()

            val iter2 = index.iterator()
            if (cellid == prev_cellid) {
                iter2.seek(cellid)
            }

            // Generate a cellunion that covers the range of empty leaf cells between
            // the last cell and this one.  Then make sure that seeking to any of
            // those cells takes us to the immediately following cell.
            if (cellid > prev_cellid) {
                for (skipped in S2CellUnion.fromBeginEnd(min_cellid, cellid)) {
                    iter2.seek(skipped)
                    assertThat(iter2.id()).isEqualTo(cellid)
                }

                // Test Prev(), Next(), and Seek().
                if (prev_cellid.isValid) {
                    iter2.seek(iter.id())
                    assertThat(iter2.prev()).isTrue()
                    assertThat(iter2.id()).isEqualTo(prev_cellid);
                    iter2.next();
                    assertThat(iter2.id()).isEqualTo(cellid);
                    iter2.seek(prev_cellid);
                    assertThat(iter2.id()).isEqualTo(prev_cellid);
                }
            }
            prev_cellid = cellid;
            min_cellid = cellid.next()
            iter.next()
        }
    }

  @Test
    fun noPoints() {
        verify()
    }

  @Test
    fun duplicatePoints() {
        repeat(10) {
            add(S2Point(1, 0, 0), 123);  // All points have same Data argument.
        }
        assertThat(index.numPoints()).isEqualTo(10)
        assertThat(contents.size).isEqualTo(10)
        verify()
        // Now remove half of the points.
        repeat(5) {
            assertThat(remove(S2Point(1, 0, 0), 123)).isTrue()
        }
        verify()
        assertThat(index.numPoints()).isEqualTo(5)
    }

  @Test
    fun randomPoints() {
        repeat(100) {
            add(S2Random.randomPoint(), S2Random.randomInt(100))
        }
        verify();
        // Now remove some of the points.
        repeat(10) {
            val iter = index.iterator()
            do {
                iter.seek(S2Random.randomCellId(S2CellId.kMaxLevel))
            } while (iter.done())
            remove(iter.point(), iter.data())
            verify();
        }
    }

  @Test
    fun emptyData() {
        // Verify that when Data is an empty class, no space is used.
        //assertEquals(sizeof(S2Point), sizeof(S2PointIndex<>::PointData));

        // Verify that points can be added and removed with an empty Data class.
        val index = S2PointIndex<Int>()
        index.add(S2Point(1, 0, 0), 1)
        index.remove(S2Point(1, 0, 0), 1)
        assertThat(index.numPoints()).isEqualTo(0)
    }

}
