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
package dilivia.s2.index.cell

import dilivia.s2.S2CellId
import dilivia.s2.S2Random
import dilivia.s2.index.cell.S2CellIndex.CellIterator
import dilivia.s2.index.cell.S2CellIndex.ContentsIterator
import dilivia.s2.index.cell.S2CellIndex.NonEmptyRangeIterator
import dilivia.s2.index.lexicon.Label
import dilivia.s2.region.S2CellUnion
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.util.*


class S2CellIndexUnitTest {

    private val logger = KotlinLogging.logger { }

    private val index = S2CellIndex()
    private val contents = mutableListOf<ValuedCell<Int>>()


    @Test
    fun test() {
        add("1/01302", 1)
        add("1/0130231", 2)
        add("1/0130231", 2)


        build()
        println(index.toDebugString())

        println("Cells:")
        println("-------------")
        val cellIterator = CellIterator(index)
        while (!cellIterator.done()) {
            val cell = cellIterator.valuedCell()
            println("  - id: ${cell.cellId}, label = ${cell.label}")
            cellIterator.next()
        }


        val rangeIterator = S2CellIndex.RangeIterator(index)
        val contentIterator = ContentsIterator(index)
        rangeIterator.begin()
        while (!rangeIterator.done()) {
            println("Current : startId = ${rangeIterator.startId()}, limit = ${rangeIterator.limitId()}, content = ${rangeIterator.contents()}")
            println("----------------------------------------------------------")

            contentIterator.startUnion(rangeIterator)
            while(!contentIterator.done()) {
                println("- ${contentIterator.cellId()} : ${contentIterator.value()}")
                contentIterator.next()
            }

            println("----------------------------------------------------------")
            rangeIterator.next()
        }
    }

    @Test
    fun empty() {
        quadraticValidate()
        // Deltas: [(cellId = Invalid: 0000000000000000, label-1) ; (cellId = Invalid: 0000000000000000, label-1) ; ]
    }

    @Test
    fun oneFaceCell() {
        add("0/", 0)
        quadraticValidate()

        println(index)
    }

    @Test
    fun oneLeafCell() {
        add("1/012301230123012301230123012301", 12)
        quadraticValidate()

        println(index)
    }

    @Test
    fun duplicateValues() {
        add("0/", 0)
        add("0/", 0)
        add("0/", 1)
        add("0/", 17)
        quadraticValidate()

        println(index)
    }

    @Test
    fun disjointCells() {
        add("0/", 0)
        add("3/", 0)
        quadraticValidate()

        println(index)
    }

    @Test
    fun nestedCells() {
        // Tests nested cells, including cases where several cells have the same
        // range_min() or range_max() and with randomly ordered labels.
        add("1/", 3)
        add("1/0", 15)
        add("1/000", 9)
        add("1/00000", 11)
        add("1/012", 6)
        add("1/01212", 5)
        add("1/312", 17)
        add("1/31200", 4)
        add("1/3120000", 10)
        add("1/333", 20)
        add("1/333333", 18)
        add("5/", 3)
        add("5/3", 31)
        add("5/3333", 27)
        quadraticValidate()

        println(index)

        val rangeIterator = S2CellIndex.RangeIterator(index)
        val contentsIterator = ContentsIterator(index)
        rangeIterator.begin()
        while (!rangeIterator.done()) {
            println("Range: start = ${rangeIterator.startId()}, limit = ${rangeIterator.limitId()}, content = ${rangeIterator.contents()}")
            contentsIterator.startUnion(rangeIterator)
            println("Contents")
            while (!contentsIterator.done()) {
                println(" - ${contentsIterator.valuedCell()}")
                contentsIterator.next()
            }
            rangeIterator.next()
        }
    }

    @Test
    fun randomCellUnions() {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        repeat(100) { i ->
            add(getRandomCellUnion(), i)
        }
        quadraticValidate()
    }

    @Test
    fun contentsIteratorSuppressesDuplicates() {
        // Checks that ContentsIterator stops reporting values once it reaches a
        // node of the cell tree that was visited by the previous call to Begin().
        add("2/1", 1)
        add("2/1", 2)
        add("2/10", 3)
        add("2/100", 4)
        add("2/102", 5)
        add("2/1023", 6)
        add("2/31", 7)
        add("2/313", 8)
        add("2/3132", 9)
        add("3/1", 10)
        add("3/12", 11)
        add("3/13", 12)
        quadraticValidate()

        val contents = ContentsIterator(index)
        expectContents("1/123", contents, emptyList())
        expectContents("2/100123", contents, listOf(Pair("2/1", 1), Pair("2/1", 2), Pair("2/10", 3), Pair("2/100", 4)))
        // Check that a second call with the same key yields no additional results.
        expectContents("2/100123", contents, emptyList())
        // Check that seeking to a different branch yields only the new values.
        expectContents("2/10232", contents, listOf(Pair("2/102", 5), Pair("2/1023", 6)))
        // Seek to a node with a different root.
        expectContents("2/313", contents, listOf(Pair("2/31", 7), Pair("2/313", 8)))
        // Seek to a descendant of the previous node.
        expectContents("2/3132333", contents, listOf(Pair("2/3132", 9)))
        // Seek to an ancestor of the previous node.
        expectContents("2/213", contents, emptyList())
        // A few more tests of incremental reporting.
        expectContents("3/1232", contents, listOf(Pair("3/1", 10), Pair("3/12", 11)))
        expectContents("3/133210", contents, listOf(Pair("3/13", 12)))
        expectContents("3/133210", contents, emptyList())
        expectContents("5/0", contents, emptyList())

        // Now try moving backwards, which is expected to yield values that were
        // already reported above.
        expectContents("3/13221", contents, listOf(Pair("3/1", 10), Pair("3/13", 12)))
        expectContents("2/31112", contents, listOf(Pair("2/31", 7)))
    }

    @Test
    fun intersectionOptimization() {
        // Tests various corner cases for the binary search optimization in
        // VisitIntersectingCells.

        add("1/001", 1)
        add("1/333", 2)
        add("2/00", 3)
        add("2/0232", 4)
        build()
        testIntersection(makeCellUnion(listOf("1/010", "1/3")))
        testIntersection(makeCellUnion(listOf("2/010", "2/011", "2/02")))
    }

    @Test
    fun intersectionRandomUnions() {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        repeat(100) { i ->
            add(getRandomCellUnion(), i)
        }
        build()
        // Now repeatedly query a cell union constructed in the same way.
        repeat(200) {
            testIntersection(getRandomCellUnion())
        }
    }

    @Test
    fun intersectionSemiRandomUnions() {
        // This test also uses random S2CellUnions, but the unions are specially
        // constructed so that interesting cases are more likely to arise.
        repeat(200) {
            S2Random.reset(it)
            index.clear()
            var id = S2CellId.fromDebugString("1/0123012301230123")
            val target = mutableListOf<S2CellId>()
            repeat(100) { i ->
                if (S2Random.oneIn(10)) add(id, i)
                if (S2Random.oneIn(4)) target.add(id)
                if (S2Random.oneIn(2)) id = id.nextWrap()
                if (S2Random.oneIn(6) && !id.isFace) id = id.parent()
                if (S2Random.oneIn(6) && !id.isLeaf) id = id.childBegin()
            }
            build()
            testIntersection(S2CellUnion(target))
        }
    }

    @Test
    fun intersection() {
        add(cellStr = "1/012301230123020021", label = 18)
        add(cellStr = "1/012301230123020110", label = 40)
        add(cellStr = "1/012301230123020112", label = 43)
        add(cellStr = "1/01230123012302012", label = 56)
        add(cellStr = "1/0123012301230211", label = 75)
        add(cellStr = "1/0123012301230212", label = 76)
        add(cellStr = "1/012301230123022", label = 85)
        add(cellStr = "1/01230123012310", label = 93)
        build()

        val union = S2CellUnion(
            S2CellId.fromDebugString("1/0123012301230123"),
            S2CellId.fromDebugString("1/012301230123013"),
            S2CellId.fromDebugString("1/0123012301230200"),
            S2CellId.fromDebugString("1/012301230123020100"),
            S2CellId.fromDebugString("1/012301230123020101"),
            S2CellId.fromDebugString("1/012301230123020102"),
            S2CellId.fromDebugString("1/01230123012302011"),
            S2CellId.fromDebugString("1/01230123012302012"),
            S2CellId.fromDebugString("1/0123012301230210"),
            S2CellId.fromDebugString("1/0123012301230211"),
            S2CellId.fromDebugString("1/01230123012310"),
            S2CellId.fromDebugString("1/01230123012311"),
        )

        testIntersection(union)
    }

    // Adds the (cell_id, label) pair to index_ and also contents_ (which is
    // used for independent validation).
    private fun add(cellId: S2CellId, label: Label) {
        index.add(cellId, label)
        contents.add(ValuedCell(cellId, label))
    }

    private fun add(cellStr: String, label: Label) {
        add(S2CellId.fromDebugString(cellStr), label)
    }

    private fun add(cellUnion: S2CellUnion, label: Label) {
        index.add(cellUnion, label)
        for (cell_id in cellUnion) {
            contents.add(ValuedCell(cell_id, label))
        }
    }

    private fun build() {
        index.build()
    }

    // Verifies that the index computes the correct set of (cell_id, label) pairs
    // for every possible leaf cell.  The running time of this function is
    // quadratic in the size of the index.
    private fun quadraticValidate() {
        build()
        verifyCellIterator()
        verifyIndexContents()
        verifyRangeIterators()
    }

    // Verifies that S2CellIndex::CellIterator visits each (cell_id, label) pair
    // exactly once.
    private fun verifyCellIterator() {
        val actual = mutableListOf<ValuedCell<Int>>()
        val iterator = CellIterator(index)
        while (!iterator.done()) {
            actual.add(ValuedCell(iterator.cellId(), iterator.value()))
            iterator.next()
        }
        expectEqual(contents, actual)
    }

    private fun verifyRangeIterators() {
        // Test Finish(), which is not otherwise tested below.
        val iterator = S2CellIndex.RangeIterator(index)
        iterator.begin()
        iterator.finish()
        assertThat(iterator.done()).isTrue()

        // And also for non-empty ranges.
        val nonEmpty = NonEmptyRangeIterator(index)
        nonEmpty.begin()
        nonEmpty.finish()
        assertThat(nonEmpty.done()).isTrue()

        // Iterate through all the ranges in the index.  We simultaneously iterate
        // through the non-empty ranges and check that the correct ranges are found.
        var prevStart = S2CellId.none
        var nonEmptyPrevStart = S2CellId.none
        iterator.begin()
        nonEmpty.begin()
        while (!iterator.done()) {
            // Check that seeking in the current range takes us to this range.
            val rangeIterator = S2CellIndex.RangeIterator(index)
            val start = iterator.startId()
            rangeIterator.seek(iterator.startId())
            assertThat(rangeIterator.startId()).isEqualTo(start)
            rangeIterator.seek(iterator.limitId().previous())
            assertThat(rangeIterator.startId()).isEqualTo(start)

            // And also for non-empty ranges.
            val nonEmpty2 = NonEmptyRangeIterator(index)
            val nonEmptyStart = nonEmpty.startId()
            nonEmpty2.seek(iterator.startId())
            assertThat(nonEmpty2.startId()).isEqualTo(nonEmptyStart)
            nonEmpty2.seek(iterator.limitId().previous())
            assertThat(nonEmpty2.startId()).isEqualTo(nonEmptyStart)

            // Test Prev() and Next().
            if (rangeIterator.prev()) {
                assertThat(rangeIterator.startId()).isEqualTo(prevStart)
                rangeIterator.next()
                assertThat(rangeIterator.startId()).isEqualTo(start)
            } else {
                assertThat(rangeIterator.startId()).isEqualTo(start)
                assertThat(prevStart).isEqualTo(S2CellId.none)
            }

            // And also for non-empty ranges.
            if (nonEmpty2.prev()) {
                assertThat(nonEmpty2.startId()).isEqualTo(nonEmptyPrevStart)
                nonEmpty2.next()
                assertThat(nonEmpty2.startId()).isEqualTo(nonEmptyStart)
            } else {
                assertThat(nonEmpty2.startId()).isEqualTo(nonEmptyStart)
                assertThat(nonEmptyPrevStart).isEqualTo(S2CellId.none)
            }

            // Keep the non-empty iterator synchronized with the regular one.
            if (!iterator.isEmpty()) {
                assertThat(nonEmpty.startId()).isEqualTo(iterator.startId())
                assertThat(nonEmpty.limitId()).isEqualTo(iterator.limitId())
                assertThat(nonEmpty.done()).isFalse()
                nonEmptyPrevStart = nonEmptyStart
                nonEmpty.next()
            }
            prevStart = start
            iterator.next()
        }
        // Verify that the NonEmptyRangeIterator is also finished.
        assertThat(nonEmpty.done()).isTrue()
    }

    // Verifies that RangeIterator and ContentsIterator can be used to determine
    // the exact set of (s2cell_id, label) pairs that contain any leaf cell.
    private fun verifyIndexContents() {
        // "min_cellid" is the first S2CellId that has not been validated yet.
        var minCellId = S2CellId.begin(S2CellId.kMaxLevel)
        val range = S2CellIndex.RangeIterator(index)
        range.begin()
        while (!range.done()) {
            logger.trace { "VerifyIndexContents: process range (startId = ${range.startId()}, contents = ${range.contents()}, limitId = ${range.limitId()} )" }
            assertThat(range.startId()).isEqualTo(minCellId)
            assertThat(range.limitId()).isGreaterThan(minCellId)
            assertThat(range.limitId().isLeaf).isTrue()
            minCellId = range.limitId()

            logger.trace { "minCellId: $minCellId" }

            // Build a list of expected (cell_id, label) pairs for this range.
            val expected = mutableListOf<ValuedCell<Int>>()
            for (x in contents) {
                if (x.cellId.rangeMin() <= range.startId() && x.cellId.rangeMax().next() >= range.limitId()) {
                    // The cell contains the entire range.
                    expected.add(x)
                } else {
                    // Verify that the cell does not intersect the range.
                    assertThat(
                        x.cellId.rangeMin() <= range.limitId().previous() && x.cellId.rangeMax() >= range.startId()
                    ).isFalse()
                }
            }
            val actual = mutableListOf<ValuedCell<Int>>()
            val contents = ContentsIterator(index)
            contents.startUnion(range)
            while (!contents.done()) {
                actual.add(ValuedCell(contents.cellId(), contents.value()))
                contents.next()
            }
            expectEqual(expected, actual)
            range.next()
        }
        assertThat(minCellId).isEqualTo(S2CellId.end(S2CellId.kMaxLevel))
    }

    // Tests that VisitIntersectingCells() and GetIntersectingLabels() return
    // correct results for the given target.
    private fun testIntersection(target: S2CellUnion) {
        val expected = mutableListOf<ValuedCell<Int>>()
        val actual = mutableListOf<ValuedCell<Int>>()
        val expectedLabels = TreeSet<Label>()
        val it = CellIterator(index)
        while (!it.done()) {
            if (target.intersects(it.cellId())) {
                expected.add(ValuedCell(it.cellId(), it.value()))
                expectedLabels.add(it.value())
            }
            it.next()
        }
        index.visitIntersectingCells(target, object : CellVisitor {
            override fun visit(cellId: S2CellId, label: Label): Boolean {
                actual.add(ValuedCell(cellId, label))
                return true
            }
        })
        expectEqual(expected, actual)
        var actualLabels = index.getIntersectingValues(target)
        actualLabels = actualLabels.distinct().sorted()
        assertThat(actualLabels).isEqualTo(expectedLabels.toList())
    }

    // Given an S2CellId "target_str" in human-readable form, expects that the
    // first leaf cell contained by this target will intersect the exact set of
    // (cell_id, label) pairs given by "expected_strs".
    private fun expectContents(
        targetStr: String,
        contents: ContentsIterator,
        expectedStrs: List<Pair<String, Label>>
    ) {
        val range = S2CellIndex.RangeIterator(index)
        range.seek(S2CellId.fromDebugString(targetStr).rangeMin())
        val expected = mutableListOf<ValuedCell<Int>>()
        val actual = mutableListOf<ValuedCell<Int>>()
        for (p in expectedStrs) {
            expected.add(ValuedCell(S2CellId.fromDebugString(p.first), p.second))
        }
        contents.startUnion(range)
        while (!contents.done()) {
            actual.add(ValuedCell(contents.cellId(), contents.value()))
            contents.next()
        }
        expectEqual(expected, actual)
    }

    // Verifies that "expected" and "actual" have the same contents.  Note that
    // duplicate values are allowed.
    private fun expectEqual(expected: List<ValuedCell<Int>>, actual: List<ValuedCell<Int>>) {
        assertThat(actual.sorted()).isEqualTo(expected.sorted())
    }

    // Creates a cell union from a small number of random cells at random levels.
    private fun getRandomCellUnion(): S2CellUnion {
        val ids = mutableListOf<S2CellId>()
        repeat(10) {
            ids.add(S2Random.randomCellId())
        }
        return S2CellUnion(ids)
    }

    private fun makeCellUnion(strs: List<String>): S2CellUnion {
        val ids = mutableListOf<S2CellId>()
        for (str in strs) {
            ids.add(S2CellId.fromDebugString(str))
        }
        return S2CellUnion(ids)
    }

}


