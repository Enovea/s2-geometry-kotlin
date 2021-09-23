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
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeClipping.kFaceClipErrorUVCoord
import dilivia.s2.edge.S2EdgeClipping.kIntersectsRectErrorUVDist
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2ShapeIndexRegion.Companion.makeS2ShapeIndexRegion
import dilivia.s2.shape.S2LaxLoopShape
import dilivia.s2.shape.S2Shape
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test


class S2ShapeIndexRegionUnitTest {

    companion object {

        fun makeCellId(str: String): S2CellId = S2CellId.fromDebugString(str)

        // Pad by at least twice the maximum error for reliable results.
        val kPadding = 2 * (kFaceClipErrorUVCoord + kIntersectsRectErrorUVDist)

        fun newPaddedCell(id: S2CellId, padding_uv: Double): S2Shape {
          val (face, i, j, _) = id.toFaceIJOrientation()
            val uv = S2CellId.ijLevelToBoundUV(i, j, id.level()).expanded(padding_uv)
            val vertices = mutableListOf<S2Point>()
            repeat(4) { k ->
                vertices.add(S2Coords.faceUvToXyz(face, uv.getVertex(k)).normalize())
            }
            return S2LaxLoopShape(vertices)
        }
    }


    @Test
    fun getCapBound() {
        val id = S2CellId.fromDebugString("3/0123012301230123012301230123")

        // Add a polygon that is slightly smaller than the cell being tested.
        val index = MutableS2ShapeIndex()
        index.add(newPaddedCell(id, -kPadding))
        val cellBound = S2Cell(id).capBound
        val indexBound = makeS2ShapeIndexRegion(index).capBound
        assertThat(indexBound.contains(cellBound)).isTrue()

        // Note that S2CellUnion::GetCapBound returns a slightly larger bound than
        // S2Cell::GetBound even when the cell union consists of a single S2CellId.
        assertThat(indexBound.radius() <= cellBound.radius() * 1.00001).isTrue()
    }

    @Test
    fun getRectBound() {
        val id = S2CellId.fromDebugString("3/0123012301230123012301230123")

        // Add a polygon that is slightly smaller than the cell being tested.
        val index = MutableS2ShapeIndex()
        index.add(newPaddedCell(id, -kPadding))
        val cellBound = S2Cell (id).rectBound
        val indexBound = makeS2ShapeIndexRegion (index).rectBound
        assertThat(cellBound).isEqualTo(indexBound)
    }

    @Test
    fun getCellUnionBoundMultipleFaces() {
        val ids = mutableListOf(makeCellId("3/00123"), makeCellId("2/11200013"))
        val index = MutableS2ShapeIndex()
        for (id in ids) index.add(newPaddedCell(id, -kPadding))
        val covering = mutableListOf<S2CellId>()
        makeS2ShapeIndexRegion(index).getCellUnionBound(covering)
        ids.sort()
        assertThat(covering).isEqualTo(ids)
    }

    @Test
    fun getCellUnionBoundOneFace() {
        // This tests consists of 3 pairs of S2CellIds.  Each pair is located within
        // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
        // We expect GetCellUnionBound to compute the smallest cell that bounds the
        // pair on each face.
        val input = mutableListOf(
            makeCellId("5/010"), makeCellId("5/0211030"),
            makeCellId("5/110230123"), makeCellId("5/11023021133"),
            makeCellId("5/311020003003030303"), makeCellId("5/311020023"),
        )
        val expected = mutableListOf(
            makeCellId("5/0"), makeCellId("5/110230"), makeCellId("5/3110200")
        )
        val index = MutableS2ShapeIndex()
        for (id in input) {
            // Add each shape 3 times to ensure that the S2ShapeIndex subdivides.
            repeat (3) {
                index.add(newPaddedCell(id, -kPadding))
            }
        }
        val actual = mutableListOf<S2CellId>()
        makeS2ShapeIndexRegion(index).getCellUnionBound(actual)
        assertThat(actual).isEqualTo(expected)
    }

    @Test
    fun containsCellMultipleShapes() {
        val id = S2CellId.fromDebugString("3/0123012301230123012301230123")

        // Add a polygon that is slightly smaller than the cell being tested.
        val index = MutableS2ShapeIndex()
        index.add(newPaddedCell(id, -kPadding))
        assertThat(makeS2ShapeIndexRegion(index).contains(S2Cell(id))).isFalse()

        // Add a second polygon that is slightly larger than the cell being tested.
        // Note that Contains() should return true if *any* shape contains the cell.
        index.add(newPaddedCell(id, kPadding))
        assertThat(makeS2ShapeIndexRegion(index).contains(S2Cell(id))).isTrue()

        // Verify that all children of the cell are also contained.
        var child = id.childBegin()
        while (child != id.childEnd()) {
            assertThat(makeS2ShapeIndexRegion(index).contains(S2Cell(child))).isTrue()
          child = child.next()
        }
    }

    @Test
    fun intersectsShrunkenCell() {
        val target = S2CellId . fromDebugString ("3/0123012301230123012301230123")

        // Add a polygon that is slightly smaller than the cell being tested.
        val index = MutableS2ShapeIndex()
        index.add(newPaddedCell(target, -kPadding))
        val region = makeS2ShapeIndexRegion (index)

        // Check that the index intersects the cell itself, but not any of the
        // neighboring cells.
        assertThat(region.mayIntersect(S2Cell(target))).isTrue()
        val nbrs = mutableListOf<S2CellId>()
        target.appendAllNeighbors(target.level(), nbrs)
        for (id in nbrs) {
            assertThat(region.mayIntersect(S2Cell(id))).isFalse()
        }
    }

    @Test
    fun intersectsExactCell() {
        val target = S2CellId . fromDebugString ("3/0123012301230123012301230123")

        // Adds a polygon that exactly follows a cell boundary.
        val index = MutableS2ShapeIndex()
        index.add(newPaddedCell(target, 0.0))
        val region = makeS2ShapeIndexRegion (index)

        // Check that the index intersects the cell and all of its neighbors.
        val ids = mutableListOf(target)
        target.appendAllNeighbors(target.level(), ids)
        for (id in ids) {
            assertThat(region.mayIntersect(S2Cell(id))).isTrue()
        }
    }

}
