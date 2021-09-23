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
package dilivia.s2.index

import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser.makeIndex
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.S2TextParser.parsePoints
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.shape.S2Shape
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import java.util.*


class S2MinDistanceTargetsUnitTest {

    @Test
    fun pointTargetUpdateMinDistanceToEdgeWhenEqual() {
        // Verifies that UpdateMinDistance only returns true when the new distance
        // is less than the old distance (not less than or equal to).
        val target = S2MinDistancePointTarget(makePoint("1:0"))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val edge = parsePoints("0:-1, 0:1")
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isTrue()
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isFalse()
    }

    @Test
    fun pointTargetUpdateMinDistanceToCellWhenEqual() {
        // Verifies that UpdateMinDistance only returns true when the new distance
        // is less than the old distance (not less than or equal to).
        val target = S2MinDistancePointTarget(makePoint("1:0"));
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val cell = S2Cell(S2CellId.fromPoint(makePoint("0:0")))
        assertThat(target.updateMinDistance(cell, dist)).isTrue()
        assertThat(target.updateMinDistance(cell, dist)).isFalse()
    }

    @Test
    fun edgeTargetUpdateMinDistanceToEdgeWhenEqual() {
        val target = S2MinDistanceEdgeTarget(makePoint("1:0"), makePoint("1:1"))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val edge = parsePoints("0:-1, 0:1")
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isTrue()
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isFalse()
    }

    @Test
    fun edgeTargetUpdateMinDistanceToCellWhenEqual() {
        val target = S2MinDistanceEdgeTarget(makePoint("1:0"), makePoint("1:1"));
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val cell = S2Cell(S2CellId.fromPoint(makePoint("0:0")))
        assertThat(target.updateMinDistance(cell, dist)).isTrue()
        assertThat(target.updateMinDistance(cell, dist)).isFalse()
    }

    @Test
    fun cellTargetUpdateMinDistanceToEdgeWhenEqual() {
        val target = S2MinDistanceCellTarget(S2Cell(S2CellId.fromPoint(makePoint("0:1"))))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val edge = parsePoints("0:-1, 0:1")
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isTrue()
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isFalse()
    }

    @Test
    fun cellTargetUpdateMinDistanceToCellWhenEqual() {
        val target = S2MinDistanceCellTarget(S2Cell(S2CellId.fromPoint(makePoint("0:1"))))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val cell = S2Cell(S2CellId.fromPoint(makePoint("0:0")))
        assertThat(target.updateMinDistance(cell, dist)).isTrue()
        assertThat(target.updateMinDistance(cell, dist)).isFalse()
    }

    @Test
    fun cellUnionTargetUpdateMinDistanceToEdgeWhenEqual() {
        val target = S2MinDistanceCellUnionTarget(S2CellUnion(S2CellId.fromPoint(makePoint("0:1"))))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val edge = parsePoints("0:-1, 0:1")
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isTrue()
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isFalse()
    }

    @Test
    fun cellUnionTargetUpdateMinDistanceToCellWhenEqual() {
        val target = S2MinDistanceCellUnionTarget(S2CellUnion(S2CellId.fromPoint(makePoint("0:1"))))
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val cell = S2Cell(S2CellId.fromPoint(makePoint("0:0")))
        assertThat(target.updateMinDistance(cell, dist)).isTrue()
        assertThat(target.updateMinDistance(cell, dist)).isFalse()
    }

    @Test
    fun shapeIndexTargetUpdateMinDistanceToEdgeWhenEqual() {
        val target_index = makeIndex("1:0 # #")
        val target = S2MinDistanceShapeIndexTarget(target_index)
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val edge = parsePoints("0:-1, 0:1");
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isTrue()
        assertThat(target.updateMinDistance(edge[0], edge[1], dist)).isFalse()
    }

    @Test
    fun shapeIndexTargetUpdateMinDistanceToCellWhenEqual() {
        val target_index = makeIndex("1:0 # #");
        val target = S2MinDistanceShapeIndexTarget(target_index)
        val dist = S2MinDistance(S1ChordAngle.infinity())
        val cell = S2Cell(S2CellId.fromPoint(makePoint("0:0")))
        assertThat(target.updateMinDistance(cell, dist)).isTrue()
        assertThat(target.updateMinDistance(cell, dist)).isFalse()
    }

    private fun getContainingShapes(target: S2MinDistanceTarget, index: S2ShapeIndex, max_shapes: Int): List<Int> {
        val shape_ids = TreeSet<Int>()
        target.visitContainingShapes(index, object : S2DistanceTarget.ShapeVisitor {
            override fun visit(containing_shape: S2Shape, target_point: S2Point): Boolean {
                shape_ids.add(containing_shape.id)
                return shape_ids.size < max_shapes
            }
        })
        return shape_ids.toList()
    }

    // Given two sorted vectors "x" and "y", returns true if x is a subset of y
// and x.size() == x_size.
    private fun isSubsetOfSize(x: List<Int>, y: List<Int>, x_size: Int): Boolean {
        if (x.size != x_size) return false;
        return y.containsAll(x)
    }

    @Test
    fun pointTargetVisitContainingShapes() {
        // Only shapes 2 and 4 should contain the target point.
        val index = makeIndex(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0"
        );
        val target = S2MinDistancePointTarget(makePoint("1:1"));
        assertThat(isSubsetOfSize(getContainingShapes(target, index, 1), listOf(2, 4), 1)).isTrue()
        assertThat(getContainingShapes(target, index, 5)).isEqualTo(listOf(2, 4))
    }

    @Test
    fun edgeTargetVisitContainingShapes() {
        // Only shapes 2 and 4 should contain the target point.
        val index = makeIndex(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0"
        );
        val target = S2MinDistanceEdgeTarget(makePoint("1:2"), makePoint("2:1"));
        assertThat(isSubsetOfSize(getContainingShapes(target, index, 1), listOf(2, 4), 1)).isTrue()
        assertThat(getContainingShapes(target, index, 5)).isEqualTo((listOf(2, 4)))
    }

    @Test
    fun CellTargetVisitContainingShapes() {
        val index = makeIndex("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
        // Only shapes 2 and 4 should contain a very small cell near 1:1.
        val cellid1 = S2CellId.fromPoint(makePoint("1:1"))
        val target1 = S2MinDistanceCellTarget(S2Cell(cellid1))
        assertThat(isSubsetOfSize(getContainingShapes(target1, index, 1), listOf(2, 4), 1)).isTrue()
        assertThat(getContainingShapes(target1, index, 5)).isEqualTo((listOf(2, 4)))

        // For a larger cell that properly contains one or more index cells, all
        // shapes that intersect the first such cell in S2CellId order are returned.
        // In the test below, this happens to again be the 1st and 3rd polygons
        // (whose shape_ids are 2 and 4).
        val cellid2 = cellid1.parent(5);
        val target2 = S2MinDistanceCellTarget(S2Cell(cellid2))
        assertThat(getContainingShapes(target2, index, 5)).isEqualTo((listOf(2, 4)))
    }

    @Test
    fun cellUnionTargetVisitContainingShapes() {
        val index = makeIndex("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
        // Shapes 2 and 4 contain the leaf cell near 1:1, while shape 3 contains the
        // leaf cell near 7:7.
        val cellid1 = S2CellId.fromPoint(makePoint("1:1"))
        val cellid2 = S2CellId.fromPoint(makePoint("7:7"))
        val target1 = S2MinDistanceCellUnionTarget(S2CellUnion(cellid1, cellid2))
        assertThat(isSubsetOfSize(getContainingShapes(target1, index, 1), listOf(2, 3, 4), 1)).isTrue()
        assertThat(getContainingShapes(target1, index, 5)).isEqualTo((listOf(2, 3, 4)))
    }

    @Test
    fun shapeIndexTargetVisitContainingShapes() {
        // Create an index containing a repeated grouping of one point, one
        // polyline, and one polygon.
        val index = makeIndex(
            "1:1 | 4:4 | 7:7 | 10:10 # 1:1, 1:2 | 4:4, 4:5 | 7:7, 7:8 " +
                    "| 10:10, 10:11 # 0:0, 0:3, 3:0 | 3:3, 3:6, 6:3 | 6:6, 6:9, 9:6 " +
                    "| 9:9, 9:12, 12:9"
        );

        // Construct a target consisting of one point, one polyline, and one polygon
        // with two loops where only the second loop is contained by a polygon in
        // the index above.
        val target_index = makeIndex("1:1 # 4:5, 5:4 # 20:20, 20:21, 21:20; 10:10, 10:11, 11:10");

        val target = S2MinDistanceShapeIndexTarget(target_index)
        // These are the shape_ids of the 1st, 2nd, and 4th polygons of "index"
        // (noting that the 4 points are represented by one S2PointVectorShape).
        assertThat(getContainingShapes(target, index, 5)).isEqualTo((listOf(5, 6, 8)))
    }

    @Test
    fun shapeIndexTargetVisitContainingShapesEmptyAndFull() {
        // Verify that VisitContainingShapes never returns empty polygons and always
        // returns full polygons (i.e., those containing the entire sphere).

        // Creating an index containing one empty and one full polygon.
        val index = makeIndex("# # empty | full");

        // Check only the full polygon is returned for a point target.
        val point_index = makeIndex("1:1 # #");
        val point_target = S2MinDistanceShapeIndexTarget(point_index)
        assertThat(getContainingShapes(point_target, index, 5)).isEqualTo((listOf(1)))

        // Check only the full polygon is returned for a full polygon target.
        val full_polygon_index = makeIndex("# # full");
        val full_target = S2MinDistanceShapeIndexTarget(full_polygon_index);
        assertThat(getContainingShapes(full_target, index, 5)).containsExactly(1)

        // Check that nothing is returned for an empty polygon target.  (An empty
        // polygon has no connected components and does not intersect anything, so
        // according to the API of getContainingShapes nothing should be returned.)
        val empty_polygon_index = makeIndex("# # empty");
        val empty_target = S2MinDistanceShapeIndexTarget(empty_polygon_index);
        assertThat(getContainingShapes(empty_target, index, 5)).isEmpty()
    }

}
