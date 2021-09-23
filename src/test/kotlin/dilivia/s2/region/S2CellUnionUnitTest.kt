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

import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2Random.oneIn
import dilivia.s2.S2Random.randomCap
import dilivia.s2.S2Random.randomCellId
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomInt
import dilivia.s2.S2Random.skewed
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.region.S2CellUnion.S2CellUnionTestPeer.fromVerbatimNoChecks
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.lang.Math.pow
import kotlin.math.max
import kotlin.math.min

class S2CellUnionUnitTest {

    private val logger = KotlinLogging.logger {  }

    /*
    class S2CellUnionTestPeer {
        public:
        // Creates a possibly invalid S2CellUnion without any checks.
        static S2CellUnion FromVerbatimNoChecks(vector<S2CellId> cell_ids) {
            return S2CellUnion(std::move(cell_ids), S2CellUnion::VERBATIM)
        }
    };*/

    @Test
    fun defaultConstructor() {
        val empty = S2CellUnion(emptyList())
        assertThat(empty.isEmpty()).isTrue()
    }

    @Test
    fun cellIdConstructor() {
        val face1Id = S2CellId.fromFace(1)
        val face1Union = S2CellUnion(listOf(face1Id))
        assertThat(face1Union.numCells()).isEqualTo(1)
        assertThat(face1Union.cellId(0)).isEqualTo(face1Id)
    }

    @Test
    fun wholeSphere() {
        val wholeSphere = S2CellUnion.wholeSphere()
        assertThat(6UL * (1UL shl 60)).isEqualTo(wholeSphere.leafCellsCovered())
        wholeSphere.expand(0)
        assertThat(S2CellUnion.wholeSphere()).isEqualTo(wholeSphere)
    }

    @Test
    fun duplicateCellsNotValid() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0))
        val cellUnion = fromVerbatimNoChecks(listOf(id, id))
        assertThat(cellUnion.isValid()).isFalse()
    }

    @Test
    fun unsortedCellsNotValid() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0)).parent(10)
        val cellUnion = fromVerbatimNoChecks(listOf(id, id.previous()))
        assertThat(cellUnion.isValid()).isFalse()
    }

    @Test
    fun invalidCellIdNotValid() {
        assertThat(S2CellId.none.isValid).isFalse()
        val cellUnion = fromVerbatimNoChecks(listOf(S2CellId.none))
        assertThat(cellUnion.isValid()).isFalse()
    }

    @Test
    fun invalidCellIdNotValidWithDebugFlag() {
        // Manually save and restore flag, to preserve test state in opensource
        // without gflags.
        assertThat(S2CellId.none.isValid).isFalse()
        val cellUnion = fromVerbatimNoChecks(listOf(S2CellId.none))
        assertThat(cellUnion.isValid()).isFalse()
    }

    @Test
    fun isNormalized() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0)).parent(10)
        val cellUnion = S2CellUnion.fromVerbatim((0..3).map { i -> id.child(i) })
        assertThat(cellUnion.isValid()).isTrue()
        assertThat(cellUnion.isNormalized()).isFalse()
    }

    fun addCells(id: S2CellId, selected: Boolean, input: MutableList<S2CellId>, expected: MutableList<S2CellId>) {
        var selected = selected
        // Decides whether to add "id" and/or some of its descendants to the
        // test case.  If "selected" is true, then the region covered by "id"
        // *must* be added to the test case (either by adding "id" itself, or
        // some combination of its descendants, or both).  If cell ids are to
        // the test case "input", then the corresponding expected result after
        // simplification is added to "expected".

        if (id == S2CellId.none) {
            // Initial call: decide whether to add cell(s) from each face.
            for (face in 0..5) {
                addCells(S2CellId.fromFace(face), false, input, expected)
            }
            return
        }
        if (id.isLeaf) {
            // The rnd.OneIn() call below ensures that the parent of a leaf cell
            // will always be selected (if we make it that far down the hierarchy).
            assert(selected)
            input.add(id)
            return
        }
        // The following code ensures that the probability of selecting a cell
        // at each level is approximately the same, i.e. we test normalization
        // of cells at all levels.
        if (!selected && oneIn(S2CellId.kMaxLevel - id.level())) {
            // Once a cell has been selected, the expected output is predetermined.
            // We then make sure that cells are selected that will normalize to
            // the desired output.
            expected.add(id)
            selected = true
        }

        // With the rnd.OneIn() constants below, this function adds an average
        // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
        // level at which the cell was first selected (level 15 on average).
        // Therefore the average number of input cells in a test case is about
        // (5/6 * 15 * 6) = 75.  The average number of output cells is about 6.

        // If a cell is selected, we add it to "input" with probability 5/6.
        var added = false
        if (selected && !oneIn(6)) {
            input.add(id)
            added = true
        }
        var numChildren = 0
        var child = id.childBegin()
        for (pos in 0..3) {
            // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
            // This intentionally may result in a cell and some of its children
            // being included in the test case.
            //
            // If the cell is not selected, on average we recurse on one child.
            // We also make sure that we do not recurse on all 4 children, since
            // then we might include all 4 children in the input case by accident
            // (in which case the expected output would not be correct).
            if (oneIn(if (selected) 12 else 4) && numChildren < 3) {
                addCells(child, selected, input, expected)
                ++numChildren
            }
            // If this cell was selected but the cell itself was not added, we
            // must ensure that all 4 children (or some combination of their
            // descendants) are added.
            if (selected && !added) addCells(child, selected, input, expected)

            child = child.next()
        }
    }

    @Test
    fun normalize() {
        // Try a bunch of random test cases, and keep track of average
        // statistics for normalization (to see if they agree with the
        // analysis above).
        var inSum = 0
        var outSum = 0
        val kIters = 2000
        for (i in 0 until kIters) {
            S2Random.reset(i)
            val input = mutableListOf<S2CellId>()
            val expected = mutableListOf<S2CellId>()
            addCells(S2CellId.none, false, input, expected)
            inSum += input.size
            outSum += expected.size
            val cellUnion = S2CellUnion(input)

            logger.trace { """
                | Test Normalize: Iteration $i
                | ----------------------------------
                | union: $cellUnion
                | expected: $expected
            """.trimMargin() }

            assertThat(cellUnion.numCells()).isEqualTo(expected.size)
            for (t in 0 until expected.size) {
                assertThat(cellUnion[t]).isEqualTo(expected[t])
            }

            // Test GetCapBound().
            val cap = cellUnion.capBound
            for (id in cellUnion) {
                assertThat(cap.contains(S2Cell(id)))
                    .withFailMessage("Cell $id is not in cap bound $cap")
                    .isTrue()
            }

            // Test Contains(S2CellId) and Intersects(S2CellId).
            for (inputId in input) {
                assertThat(cellUnion.contains(inputId)).isTrue()
                assertThat(cellUnion.contains(inputId.toPoint())).isTrue()
                assertThat(cellUnion.intersects(inputId)).isTrue()
                if (!inputId.isFace) {
                    assertThat(cellUnion.intersects(inputId.parent())).isTrue()
                    if (inputId.level() > 1) {
                        assertThat(cellUnion.intersects(inputId.parent().parent())).isTrue()
                        assertThat(cellUnion.intersects(inputId.parent(0))).isTrue()
                    }
                }
                if (!inputId.isLeaf) {
                    assertThat(cellUnion.contains(inputId.childBegin())).isTrue()
                    assertThat(cellUnion.intersects(inputId.childBegin())).isTrue()
                    assertThat(cellUnion.contains(inputId.childEnd().previous())).isTrue()
                    assertThat(cellUnion.intersects(inputId.childEnd().previous())).isTrue()
                    assertThat(cellUnion.contains(inputId.childBegin(S2CellId.kMaxLevel))).isTrue()
                    assertThat(cellUnion.intersects(inputId.childBegin(S2CellId.kMaxLevel))).isTrue()
                }
            }
            for (expected_id in expected) {
                if (!expected_id.isFace) {
                    assertThat(!cellUnion.contains(expected_id.parent())).isTrue()
                    assertThat(!cellUnion.contains(expected_id.parent(0))).isTrue()
                }
            }

            // Test Contains(S2CellUnion), Intersects(S2CellUnion), Union(),
            // Intersection(), and Difference().
            val x = mutableListOf<S2CellId>()
            val y = mutableListOf<S2CellId>()
            val xOrY = mutableListOf<S2CellId>()
            val xAndY = mutableListOf<S2CellId>()
            for (input_id in input) {
                val inX = oneIn(2)
                val inY = oneIn(2)
                if (inX) x.add(input_id)
                if (inY) y.add(input_id)
                if (inX || inY) xOrY.add(input_id)
            }
            val xcells = S2CellUnion(x)
            val ycells = S2CellUnion(y)
            val xOrYExpected = S2CellUnion(xOrY)
            val xOrYCells = xcells.union(ycells)
            assertThat(xOrYCells == xOrYExpected).isTrue()

            // Compute the intersection of "x" with each cell of "y",
            // check that this intersection is correct, and append the
            // results to x_and_y_expected.
            for (yid in ycells) {
                val ucells = xcells.intersection(yid)
                for (xid in xcells) {
                    if (xid.contains(yid)) {
                        assertThat(ucells.numCells() == 1 && ucells[0] == yid).isTrue()
                    } else if (yid.contains(xid)) {
                        assertThat(ucells.contains(xid)).isTrue()
                    }
                }
                for (uid in ucells) {
                    assertThat(xcells.contains(uid)).isTrue()
                    assertThat(yid.contains(uid)).isTrue()
                }
                xAndY.addAll(ucells)
            }
            val xAndYExpected = S2CellUnion(xAndY)
            val xAndYCells = xcells.intersection(ycells)
            assertThat(xAndYCells == xAndYExpected).isTrue()

            val xMinusYCells = xcells.difference(ycells)
            val yMinusXCells = ycells.difference(xcells)
            assertThat(xcells.contains(xMinusYCells)).isTrue()
            assertThat(!xMinusYCells.intersects(ycells)).isTrue()
            assertThat(ycells.contains(yMinusXCells)).isTrue()
            assertThat(!yMinusXCells.intersects(xcells)).isTrue()
            assertThat(!xMinusYCells.intersects(yMinusXCells)).isTrue()

            val diffIntersectionUnion = xMinusYCells.union(yMinusXCells).union(xAndYCells)
            assertThat(diffIntersectionUnion == xOrYCells).isTrue()

            val test = mutableListOf<S2CellId>()
            val dummy = mutableListOf<S2CellId>()
            addCells(S2CellId.none, false, test, dummy)
            for (testId in test) {
                var contains = false
                var intersects = false
                for (expected_id in expected) {
                    if (expected_id.contains(testId)) contains = true
                    if (expected_id.intersects(testId)) intersects = true
                }
                assertThat(cellUnion.contains(testId)).isEqualTo(contains)
                assertThat(cellUnion.intersects(testId)).isEqualTo(intersects)
            }
        }
        println(String.format("avg in %.2f, avg out %.2f\n", inSum.toDouble() / kIters, outSum.toDouble() / kIters))
    }

    // Return the maximum geodesic distance from "axis" to any point of
// "covering".
    fun getRadius(covering: S2CellUnion, axis: S2Point): Double {
        var maxDist = 0.0
        for (id in covering) {
            val cell = S2Cell(id)
            for (j in 0..3) {
                val a = cell.getVertex(j)
                val b = cell.getVertex(j + 1)
                var dist = 0.0
                // The maximum distance is not always attained at a cell vertex: if at
                // least one vertex is in the opposite hemisphere from "axis" then the
                // maximum may be attained along an edge.  We solve this by computing
                // the minimum distance from the edge to (-axis) instead.  We can't
                // simply do this all the time because S2::GetDistance() has
                // poor accuracy when the result is close to Pi.
                //
                // TODO(ericv): Improve S2::GetDistance() accuracy near Pi.
                if (a.angle(axis) > M_PI_2 || b.angle(axis) > M_PI_2) {
                    dist = M_PI - S2EdgeDistances.getDistance(-axis, a, b).radians
                } else {
                    dist = a.angle(axis)
                }
                maxDist = max(maxDist, dist)
            }
        }
        return maxDist
    }

    @Test
    fun expand() {
        // This test generates coverings for caps of random sizes, expands
        // the coverings by a random radius, and then make sure that the new
        // covering covers the expanded cap.  It also makes sure that the
        // new covering is not too much larger than expected.

        val coverer = S2RegionCoverer()
        for (i in 0 until 1000) {
            val cap = randomCap(S2Cell.averageArea(S2CellId.kMaxLevel), 4 * M_PI)

            // Expand the cap area by a random factor whose log is uniformly
            // distributed between 0 and log(1e2).
            val expandedCap = S2Cap.fromCenterHeight(cap.center, min(2.0, pow(1e2, randomDouble()) * cap.height))

            val radius = (expandedCap.radius() - cap.radius()).radians
            val maxLevelDiff = randomInt(8)

            // Generate a covering for the original cap, and measure the maximum
            // distance from the cap center to any point in the covering.
            coverer.setMaxCells(1 + skewed(10))
            val covering = coverer.getCovering(cap)
            checkCovering(cap, covering, true)
            val coveringRadius = getRadius(fromVerbatimNoChecks(covering.cellIds()), cap.center)

            // This code duplicates the logic in Expand(min_radius, max_level_diff)
            // that figures out an appropriate cell level to use for the expansion.
            var minLevel = S2CellId.kMaxLevel
            for (id in covering) {
                minLevel = min(minLevel, id.level())
            }
            val expandLevel = min(minLevel + maxLevelDiff, S2Coords.projection.kMinWidth.getLevelForMinValue(radius))

            // Generate a covering for the expanded cap, and measure the new maximum
            // distance from the cap center to any point in the covering.
            covering.expand(S1Angle.radians(radius), maxLevelDiff)
            checkCovering(expandedCap, covering, false)
            val expandedCoveringRadius = getRadius(fromVerbatimNoChecks(covering.cellIds()), cap.center)

            // If the covering includes a tiny cell along the boundary, in theory the
            // maximum angle of the covering from the cap center can increase by up to
            // twice the maximum length of a cell diagonal.
            assertThat(expandedCoveringRadius - coveringRadius <= 2 * S2Coords.projection.kMaxDiag.getValue(expandLevel)).isTrue()
        }
    }

    fun testFromMinMax(min_id: S2CellId, max_id: S2CellId) {
        val cellUnion = S2CellUnion.fromMinMax(min_id, max_id)
        val cellIds = cellUnion.cellIds()

        assertThat(cellIds.size > 0).isTrue()
        assertThat(cellIds.first().rangeMin()).isEqualTo(min_id)
        assertThat(cellIds.last().rangeMax()).isEqualTo(max_id)
        for (i in cellIds.indices) {
            if (i > 0) {
                assertThat(cellIds[i - 1].rangeMax().next()).isEqualTo(cellIds[i].rangeMin())
            }
        }
        assertThat(cellUnion.isNormalized()).isTrue()
    }

    @Test
    fun fromMinMax() {
        // Check the very first leaf cell and face cell.
        val face1Id = S2CellId.fromFace(0)
        testFromMinMax(face1Id.rangeMin(), face1Id.rangeMin())
        testFromMinMax(face1Id.rangeMin(), face1Id.rangeMax())

        // Check the very last leaf cell and face cell.
        val face5Id = S2CellId.fromFace(5)
        testFromMinMax(face5Id.rangeMin(), face5Id.rangeMax())
        testFromMinMax(face5Id.rangeMax(), face5Id.rangeMax())

        // Check random ranges of leaf cells.
        for (iter in 0 until 100) {
            var x = randomCellId(S2CellId.kMaxLevel)
            var y = randomCellId(S2CellId.kMaxLevel)
            if (x > y) {
                val temp = x
                x = y
                y = temp
            }
            testFromMinMax(x, y)
        }
    }

    @Test
    fun fromBeginEnd() {
        // Since FromMinMax() is implemented in terms of FromBeginEnd(), we
        // focus on test cases that generate an empty range.
        val idBegin = S2CellId.begin(S2CellId.kMaxLevel)
        var cellUnion = S2CellUnion.fromBeginEnd(idBegin, idBegin)
        assertThat(cellUnion.isEmpty()).isTrue()

        // Test the full sphere.
        val idEnd = S2CellId.end(S2CellId.kMaxLevel)
        cellUnion = S2CellUnion.fromBeginEnd(idBegin, idEnd)
        assertThat(cellUnion.numCells()).isEqualTo(6)
        for (id in cellUnion) {
            assertThat(id.isFace).isTrue()
        }
    }

    @Test
    fun empty() {
        val emptyCellUnion = S2CellUnion()
        val face1Id = S2CellId.fromFace(1)

        // Normalize()
        emptyCellUnion.normalize()
        assertThat(emptyCellUnion.isEmpty()).isTrue()

        // Denormalize(...)
        val output = emptyCellUnion.denormalize(0, 2)
        assertThat(emptyCellUnion.isEmpty()).isTrue()

        // Contains(...)
        assertThat(emptyCellUnion.contains(face1Id)).isFalse()
        assertThat(emptyCellUnion.contains(emptyCellUnion)).isTrue()

        // Intersects(...)
        assertThat(emptyCellUnion.intersects(face1Id)).isFalse()
        assertThat(emptyCellUnion.intersects(emptyCellUnion)).isFalse()

        // Union(...)
        val cellUnion = emptyCellUnion.union(emptyCellUnion)
        assertThat(cellUnion.isEmpty()).isTrue()

        // Intersection(...)
        var intersection = emptyCellUnion.intersection(face1Id)
        assertThat(intersection.isEmpty()).isTrue()
        intersection = emptyCellUnion.intersection(emptyCellUnion)
        assertThat(intersection.isEmpty()).isTrue()

        // Difference(...)
        val difference = emptyCellUnion.difference(emptyCellUnion)
        assertThat(difference.numCells()).isEqualTo(0)

        // Expand(...)
        emptyCellUnion.expand(S1Angle.radians(1), 20)
        assertThat(emptyCellUnion.isEmpty()).isTrue()
        emptyCellUnion.expand(10)
        assertThat(emptyCellUnion.isEmpty()).isTrue()
    }

    @Test
    fun clear() {
        val face1Id = S2CellId.fromFace(1)
        val face1Union = S2CellUnion(face1Id)

        assertThat(face1Union.numCells()).isEqualTo(1)
        assertThat(face1Union.cellIds().size).isEqualTo(1)

        face1Union.clear()
        assertThat(face1Union.numCells()).isEqualTo(0)
        assertThat(face1Union.cellIds().size).isEqualTo(0)
    }

    @Test
    fun leafCellsCovered() {
        var cellUnion = S2CellUnion()
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(0UL)

        val ids = mutableListOf<S2CellId>()
        // One leaf cell on face 0.
        ids.add(S2CellId.fromFace(0).childBegin(S2CellId.kMaxLevel))
        cellUnion = S2CellUnion(ids)
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(1UL)

        // Face 0 itself (which includes the previous leaf cell).
        ids.add(S2CellId.fromFace(0))
        cellUnion = S2CellUnion(ids)
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(1UL shl 60)
        // Five faces.
        cellUnion.expand(0)
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(5UL shl 60)
        // Whole world.
        cellUnion.expand(0)
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(6UL shl 60)

        // Add some disjoint cells.
        ids.add(S2CellId.fromFace(1).childBegin(1))
        ids.add(S2CellId.fromFace(2).childBegin(2))
        ids.add(S2CellId.fromFace(2).childEnd(2).previous())
        ids.add(S2CellId.fromFace(3).childBegin(14))
        ids.add(S2CellId.fromFace(4).childBegin(27))
        ids.add(S2CellId.fromFace(4).childEnd(15).previous())
        ids.add(S2CellId.fromFace(5).childBegin(30))
        cellUnion = S2CellUnion(ids)
        val expected = 1UL + (1UL shl 6) + (1UL shl 30) + (1UL shl 32) + (2UL shl 56) + (1UL shl 58) + (1UL shl 60)
        assertThat(cellUnion.leafCellsCovered()).isEqualTo(expected)
    }

    @Test
    fun cellUnion() {
        val x = listOf("5/1322210002222", "5/132221001", "5/132221002", "5/132221003")
                .map { S2CellId.fromDebugString(it) }
        val y = listOf("5/13222100232000", "5/13222100232003", "5/13222100233110", "5/13222100233111")
                .map { S2CellId.fromDebugString(it) }
        val out = mutableListOf<S2CellId>()
        S2CellUnion.getIntersection(x, y, out)
        println(out)
        assertThat(out).isEqualTo(y)
    }

    companion object {

        fun checkCovering(region: S2Region, covering: S2CellUnion, checkTight: Boolean, id: S2CellId = S2CellId()) {
            if (!id.isValid) {
                for (face in 0..5) {
                    checkCovering(region, covering, checkTight, S2CellId.fromFace(face))
                }
                return
            }

            if (!region.mayIntersect(S2Cell(id))) {
                // If region does not intersect id, then neither should the covering.
                if (checkTight) assertThat(!covering.intersects(id)).isTrue()

            } else if (!covering.contains(id)) {
                // The region may intersect id, but we can't assert that the covering
                // intersects id because we may discover that the region does not actually
                // intersect upon further subdivision.  (MayIntersect is not exact.)
                assertThat(!region.contains(S2Cell(id))).isTrue()
                assertThat(!id.isLeaf).isTrue()
                val end = id.childEnd()
                var child = id.childBegin()
                while (child != end) {
                    checkCovering(region, covering, checkTight, child)
                    child = child.next()
                }
            }
        }

    }

}
