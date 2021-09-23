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
package dilivia.s2

import dilivia.PreConditions.requireGE
import dilivia.PreConditions.requireLT
import dilivia.math.R2Rect
import dilivia.math.vectors.R2Point
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2Coords.kPosToOrientation
import dilivia.s2.region.S2Cap
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import java.util.logging.Logger
import kotlin.math.IEEErem
import kotlin.math.abs
import kotlin.math.min

/**
 */
class S2CellIdUnitTest {

    @Test
    fun defaultConstructor() {
        val id = S2CellId()
        assertThat(id.id).isEqualTo(0UL)
        assertThat(id.isValid).isFalse()
    }

    @Test
    fun cellIdHash() {
        repeat(15) {
            assertThat(getCellId(0, 90).hashCode()).isEqualTo(805306369)
        }
    }

    @Test
    fun faceDefinitions() {
        assertThat(getCellId(0, 0).face()).isEqualTo(0)
        assertThat(getCellId(0, 90).face()).isEqualTo(1)
        assertThat(getCellId(90, 0).face()).isEqualTo(2)
        assertThat(getCellId(0, 180).face()).isEqualTo(3)
        assertThat(getCellId(0, -90).face()).isEqualTo(4)
        assertThat(getCellId(-90, 0).face()).isEqualTo(5)
    }


    @Test
    fun fromFace() {
        for (face in 0..5) {
            assertThat(S2CellId.fromFace(face)).isEqualTo(S2CellId.fromFacePosLevel(face, 0UL, 0))
        }
    }

    @Test
    fun parentChildRelationships() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678UL, S2CellId.kMaxLevel - 4)
        assertThat(id.isValid).isTrue()
        assertThat(id.face()).isEqualTo(3)
        assertThat(id.pos()).isEqualTo(0x12345700UL)
        assertThat(id.level()).isEqualTo(S2CellId.kMaxLevel - 4)
        assertThat(id.isLeaf).isFalse()

        assertThat(id.childBegin(id.level() + 2).pos()).isEqualTo(0x12345610UL)
        assertThat(id.childBegin().pos()).isEqualTo(0x12345640UL)
        assertThat(id.parent().pos()).isEqualTo(0x12345400UL)
        assertThat(id.parent(id.level() - 2).pos()).isEqualTo(0x12345000UL)

        // Check ordering of children relative to parents.
        assertThat(id.childBegin()).isLessThan(id)
        assertThat(id.childEnd()).isGreaterThan(id)
        assertThat(id.childBegin().next().next().next().next()).isEqualTo(id.childEnd())
        assertThat(id.childBegin(S2CellId.kMaxLevel)).isEqualTo(id.rangeMin())
        assertThat(id.childEnd(S2CellId.kMaxLevel)).isEqualTo(id.rangeMax().next())

        // Check that cells are represented by the position of their center
        // along the Hilbert curve.
        assertThat(id.rangeMin().id + id.rangeMax().id).isEqualTo(2UL * id.id)
    }

    @Test
    fun sentinelRangeMinMax() {
        assertThat(S2CellId.sentinel.rangeMin()).isEqualTo(S2CellId.sentinel)
        assertThat(S2CellId.sentinel.rangeMax()).isEqualTo(S2CellId.sentinel)
    }

    @Test
    fun centerSiTi() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678UL, S2CellId.kMaxLevel)
        // Check that the (si, ti) coordinates of the center end in a 1 followed by (30 - level) 0s.

        // Leaf level, 30.
        var faceSiTi = id.centerSiTi
        assertThat(faceSiTi.si and 1U).isEqualTo(1U shl 0)
        assertThat(faceSiTi.ti and 1U).isEqualTo(1U shl 0)

        // Level 29.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 1).centerSiTi
        assertThat(faceSiTi.si and 3U).isEqualTo(1U shl 1)
        assertThat(faceSiTi.ti and 3U).isEqualTo(1U shl 1)

        // Level 28.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 2).centerSiTi
        assertThat(faceSiTi.si and 7U).isEqualTo(1U shl 2)
        assertThat(faceSiTi.ti and 7U).isEqualTo(1U shl 2)

        // Level 20.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 10).centerSiTi
        assertThat(faceSiTi.si and ((1U shl 11) - 1U)).isEqualTo(1U shl 10)
        assertThat(faceSiTi.ti and ((1U shl 11) - 1U)).isEqualTo(1U shl 10)

        // Level 10.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 20).centerSiTi
        assertThat(faceSiTi.si and ((1U shl 21) - 1U)).isEqualTo(1U shl 20)
        assertThat(faceSiTi.ti and ((1U shl 21) - 1U)).isEqualTo(1U shl 20)

        // Level 0.
        faceSiTi = id.parent(0).centerSiTi
        assertThat(faceSiTi.si and ((1U shl 31) - 1U)).isEqualTo(1U shl 30)
        assertThat(faceSiTi.ti and ((1U shl 31) - 1U)).isEqualTo(1U shl 30)
    }

    @Test
    fun wrapping() {
        // Check wrapping from beginning of Hilbert curve to end and vice versa.
        assertThat(S2CellId.begin(0).prevWrap()).isEqualTo(S2CellId.end(0).previous())

        assertThat(S2CellId.begin(S2CellId.kMaxLevel).prevWrap()).isEqualTo(
            S2CellId.fromFacePosLevel(5, 0UL.inv() shr S2CellId.kFaceBits, S2CellId.kMaxLevel)
        )
        assertThat(S2CellId.begin(S2CellId.kMaxLevel).advanceWrap(-1)).isEqualTo(
            S2CellId.fromFacePosLevel(5, 0UL.inv() shr S2CellId.kFaceBits, S2CellId.kMaxLevel)
        )

        assertThat(S2CellId.end(4).previous().nextWrap()).isEqualTo(S2CellId.begin(4))
        assertThat(S2CellId.end(4).advance(-1).advanceWrap(1)).isEqualTo(S2CellId.begin(4))

        assertThat(S2CellId.end(S2CellId.kMaxLevel).previous().nextWrap()).isEqualTo(
            S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel)

        )
        assertThat(S2CellId.end(S2CellId.kMaxLevel).advance(-1).advanceWrap(1)).isEqualTo(
            S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel)
        )
    }

    @Test
    fun advance() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678.toULong(), S2CellId.kMaxLevel - 4)
        // Check basic properties of advance().
        assertThat(S2CellId.begin(0).advance(7)).isEqualTo(S2CellId.end(0))
        assertThat(S2CellId.begin(0).advance(12)).isEqualTo(S2CellId.end(0))
        assertThat(S2CellId.end(0).advance(-7)).isEqualTo(S2CellId.begin(0))
        assertThat(S2CellId.end(0).advance(-12000000)).isEqualTo(S2CellId.begin(0))
        val numLevel5Cells = 6 shl (2 * 5)
        assertThat(S2CellId.end(5).advance(500L - numLevel5Cells)).isEqualTo(S2CellId.begin(5).advance(500))
        assertThat(id.childBegin(S2CellId.kMaxLevel).advance(256)).isEqualTo(id.next().childBegin(S2CellId.kMaxLevel))
        assertThat(
            S2CellId.fromFacePosLevel(1, 0UL, S2CellId.kMaxLevel).advance(4L shl (2 * S2CellId.kMaxLevel))
        ).isEqualTo(
            S2CellId.fromFacePosLevel(5, 0UL, S2CellId.kMaxLevel)
        )

        // Check basic properties of advance_wrap().
        assertThat(S2CellId.begin(0).advanceWrap(7)).isEqualTo(S2CellId.fromFace(1))
        assertThat(S2CellId.begin(0).advanceWrap(12)).isEqualTo(S2CellId.begin(0))
        assertThat(S2CellId.fromFace(5).advanceWrap(-7)).isEqualTo(S2CellId.fromFace(4))
        assertThat(S2CellId.begin(0).advanceWrap(-12000000)).isEqualTo(S2CellId.begin(0))
        assertThat(S2CellId.begin(5).advanceWrap(-11788)).isEqualTo(S2CellId.begin(5).advanceWrap(6644))
        assertThat(id.childBegin(S2CellId.kMaxLevel).advanceWrap(256)).isEqualTo(
            id.next().childBegin(S2CellId.kMaxLevel)
        )
        assertThat(
            S2CellId.fromFacePosLevel(5, 0UL, S2CellId.kMaxLevel).advanceWrap(2L shl (2 * S2CellId.kMaxLevel))
        ).isEqualTo(
            S2CellId.fromFacePosLevel(1, 0UL, S2CellId.kMaxLevel)
        )
    }

    @Test
    fun distanceFromBegin() {
        assertThat(S2CellId.end(0).distanceFromBegin()).isEqualTo(6UL)
        assertThat(S2CellId.end(S2CellId.kMaxLevel).distanceFromBegin()).isEqualTo(
            6UL * (1L shl (2 * S2CellId.kMaxLevel)).toULong()
        )

        assertThat(S2CellId.begin(0).distanceFromBegin()).isEqualTo(0UL)
        assertThat(S2CellId.begin(S2CellId.kMaxLevel).distanceFromBegin()).isEqualTo(0UL)

        val id = S2CellId.fromFacePosLevel(3, 0x12345678UL, S2CellId.kMaxLevel - 4)
        assertThat(S2CellId.begin(id.level()).advance(id.distanceFromBegin().toLong())).isEqualTo(id)
    }

    @Test
    fun maximumTile() {
        // This method is tested more thoroughly in s2cell_union_test.cc.
        for (iter in 0 until 1000) {
            val id = S2Random.randomCellId(10)

            // Check that "limit" is returned for tiles at or beyond "limit".
            assertThat(id.maximumTile(id)).isEqualTo(id)
            assertThat(id.child(0).maximumTile(id)).isEqualTo(id)
            assertThat(id.child(1).maximumTile(id)).isEqualTo(id)
            assertThat(id.next().maximumTile(id)).isEqualTo(id)
            assertThat(id.maximumTile(id.child(0))).isEqualTo(id.child(0))

            // Check that the tile size is increased when possible.
            assertThat(id.child(0).maximumTile(id.next())).isEqualTo(id)
            assertThat(id.child(0).maximumTile(id.next().child(0))).isEqualTo(id)
            assertThat(id.child(0).maximumTile(id.next().child(1).child(0))).isEqualTo(id)
            assertThat(id.child(0).child(0).maximumTile(id.next())).isEqualTo(id)
            assertThat(id.child(0).child(0).child(0).maximumTile(id.next())).isEqualTo(id)

            // Check that the tile size is decreased when necessary.
            assertThat(id.maximumTile(id.child(0).next())).isEqualTo(id.child(0))
            assertThat(id.maximumTile(id.child(0).next().child(0))).isEqualTo(id.child(0))
            assertThat(id.maximumTile(id.child(0).next().child(1))).isEqualTo(id.child(0))
            assertThat(id.maximumTile(id.child(0).child(0).next())).isEqualTo(id.child(0).child(0))
            assertThat(id.maximumTile(id.child(0).child(0).child(0).next())).isEqualTo(id.child(0).child(0).child(0))

            // Check that the tile size is otherwise unchanged.
            assertThat(id.maximumTile(id.next())).isEqualTo(id)
            assertThat(id.maximumTile(id.next().child(0))).isEqualTo(id)
            assertThat(id.maximumTile(id.next().child(1).child(0))).isEqualTo(id)
        }
    }

    @Test
    fun getCommonAncestorLevel() {
        // Two identical cell ids.
        assertThat(S2CellId.fromFace(0).getCommonAncestorLevel(S2CellId.fromFace(0))).isEqualTo(0)
        assertThat(
            S2CellId.fromFace(0).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(0).childBegin(30))
        ).isEqualTo(30)

        // One cell id is a descendant of the other.
        assertThat(S2CellId.fromFace(0).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(0))).isEqualTo(0)
        assertThat(S2CellId.fromFace(5).getCommonAncestorLevel(S2CellId.fromFace(5).childEnd(30).previous())).isEqualTo(
            0
        )

        // Two cells that have no common ancestor.
        assertThat(S2CellId.fromFace(0).getCommonAncestorLevel(S2CellId.fromFace(5))).isEqualTo(-1)
        assertThat(
            S2CellId.fromFace(2).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(3).childEnd(20))
        ).isEqualTo(-1)

        // Two cells that have a common ancestor distinct from both of them.
        assertThat(
            S2CellId.fromFace(5).childBegin(9).next().childBegin(15)
                .getCommonAncestorLevel(S2CellId.fromFace(5).childBegin(9).childBegin(20))
        ).isEqualTo(8)
        assertThat(
            S2CellId.fromFace(0).childBegin(2).childBegin(30)
                .getCommonAncestorLevel(S2CellId.fromFace(0).childBegin(2).next().childBegin(5))
        ).isEqualTo(1)
    }

    @Test
    fun inverses() {
        // Check the conversion of random leaf cells to S2LatLngs and back.
        for (i in 0 until 200000) {
            val id = S2Random.randomCellId(S2CellId.kMaxLevel)
            assertThat(id.isLeaf).isTrue()
            assertThat(id.level()).isEqualTo(S2CellId.kMaxLevel)
            val center = id.toLatLng()
            assertThat(S2CellId.fromLatLng(center).id).isEqualTo(id.id)
        }
    }

    @Test
    fun tokens() {
        // Test random cell ids at all levels.
        for (i in 0 until 10000) {
            val id = S2Random.randomCellId()
            val token = id.toToken();
            assertThat(token.length).isLessThanOrEqualTo(16)
            assertThat(S2CellId.fromToken(token)).isEqualTo(id)
        }
        // Check that invalid cell ids can be encoded, and round-trip is
        // the identity operation.
        var token = S2CellId.none.toToken()
        assertThat(S2CellId.fromToken(token)).isEqualTo(S2CellId.none)

        // Sentinel is invalid.
        token = S2CellId.sentinel.toToken();
        assertThat(S2CellId.fromToken(token)).isEqualTo(S2CellId.sentinel)

        // Check an invalid face.
        token = S2CellId.fromFace(7).toToken();
        assertThat(S2CellId.fromToken(token)).isEqualTo(S2CellId.fromFace(7))

        // Check that supplying tokens with non-alphanumeric characters
        // returns S2CellId.none().
        assertThat(S2CellId.fromToken("876b e99")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromToken("876bee99\n")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromToken("876[ee99")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromToken(" 876bee99")).isEqualTo(S2CellId.none)
    }

    @Test
    fun containment() {
        // Test contains() and intersects().
        val parentMap = mutableMapOf<S2CellId, S2CellId>()
        val cells = mutableListOf<S2CellId>()
        for (face in 0..5) {
            expandCell(S2CellId.fromFace(face), cells, parentMap)
        }
        for (end_id in cells) {
            for (begin_id in cells) {
                var contained = true
                var id = begin_id
                while (id != end_id) {
                    if (parentMap[id] == null) {
                        contained = false
                        break
                    }
                    id = parentMap.getValue(id)
                }
                assertThat(end_id.contains(begin_id)).isEqualTo(contained)
                assertThat(begin_id in end_id.rangeMin()..end_id.rangeMax()).isEqualTo(contained)
                assertThat(end_id.contains(begin_id) || begin_id.contains(end_id)).isEqualTo(end_id.intersects(begin_id))
            }
        }
    }

    private fun expandCell(parent: S2CellId, cells: MutableList<S2CellId>, parent_map: MutableMap<S2CellId, S2CellId>) {
        cells.add(parent)
        if (parent.level() == kMaxExpandLevel) return
        val (face, _, _, orientation) = parent.toFaceIJOrientation(true)
        assertThat(face).isEqualTo(parent.face())

        var child = parent.childBegin()
        var pos = 0
        while (child != parent.childEnd()) {
            parent_map[child] = parent
            // Do some basic checks on the children.
            assertThat(parent.child(pos)).isEqualTo(child)
            assertThat(child.childPosition()).isEqualTo(pos)
            // Test child_position(level) on all the child's ancestors.
            var ancestor = child
            while (ancestor.level() >= 1) {
                assertThat(ancestor.childPosition()).isEqualTo(child.childPosition(ancestor.level()))
                ancestor = parent_map[ancestor]!!
            }
            assertThat(child.childPosition(child.level())).isEqualTo(pos)
            assertThat(child.level()).isEqualTo(parent.level() + 1)
            assertThat(child.isLeaf).isFalse()

            val (childFace, _, _, child_orientation) = child.toFaceIJOrientation(true)
            assertThat(childFace).isEqualTo(face)
            assertThat(child_orientation).isEqualTo(orientation!! xor kPosToOrientation[pos])
            expandCell(child, cells, parent_map)

            child = child.next()
            ++pos
        }
    }

    @Test
    fun continuity() {
        // Make sure that sequentially increasing cell ids form a continuous
        // path over the surface of the sphere, i.e. there are no
        // discontinuous jumps from one region to another.

        val maxDist = S2Coords.projection.kMaxEdge.getValue(kMaxWalkLevel)
        val end = S2CellId.end(kMaxWalkLevel)
        var id = S2CellId.begin(kMaxWalkLevel)
        while (id != end) {
            assertThat(id.toPointRaw().angle(id.nextWrap().toPointRaw())).isLessThanOrEqualTo(maxDist)
            assertThat(id.advanceWrap(1)).isEqualTo(id.nextWrap())
            assertThat(id.nextWrap().advanceWrap(-1)).isEqualTo(id)

            // Check that the ToPointRaw() returns the center of each cell
            // in (s,t) coordinates.
            val (_, u, v) = S2Coords.xyzToFaceUv(id.toPointRaw())
            val kCellSize = 1.0 / (1 shl kMaxWalkLevel)
            assertThat(S2Coords.uvToSt(u).IEEErem(0.5 * kCellSize)).isCloseTo(0.0, Offset.offset(1e-15))
            assertThat(S2Coords.uvToSt(v).IEEErem(0.5 * kCellSize)).isCloseTo(0.0, Offset.offset(1e-15))

            id = id.next()
        }
    }

    @Test
    fun coverage() {
        // Make sure that random points on the sphere can be represented to the
        // expected level of accuracy, which in the worst case is sqrt(2/3) times
        // the maximum arc length between the points on the sphere associated with
        // adjacent values of "i" or "j".  (It is sqrt(2/3) rather than 1/2 because
        // the cells at the corners of each face are stretched -- they have 60 and
        // 120 degree angles.)
        val maxDist = 0.5 * S2Coords.projection.kMaxDiag.getValue(S2CellId.kMaxLevel)
        for (i in 0 until 1000000) {
            val p = S2Random.randomPoint()
            val q = S2CellId.fromPoint(p).toPointRaw()
            assertThat(p.angle(q)).isLessThanOrEqualTo(maxDist)
        }
    }

    @Test
    fun neighbors() {
        // Check the edge neighbors of face 1.
        val outFaces = intArrayOf(5, 3, 2, 0)
        val faceNbrs = S2CellId.fromFace(1).getEdgeNeighbors()
        for (i in 0..3) {
            assertThat(faceNbrs[i].isFace).isTrue()
            assertThat(faceNbrs[i].face()).isEqualTo(outFaces[i])
        }

        // Check the edge neighbors of the corner cells at all levels.  This case is
        // trickier because it requires projecting onto adjacent faces.
        val kMaxIJ = S2CellId.kMaxSize - 1
        for (level in 1..S2CellId.kMaxLevel) {
            val id = S2CellId.fromFaceIJ(1, 0, 0).parent(level)
            val nbrs = id.getEdgeNeighbors()
            // These neighbors were determined manually using the face and axis
            // relationships defined in s2coords.cc.
            val sizeIj = S2CellId.getSizeIJ(level)
            assertThat(nbrs[0]).isEqualTo(S2CellId.fromFaceIJ(5, kMaxIJ, kMaxIJ).parent(level))
            assertThat(nbrs[1]).isEqualTo(S2CellId.fromFaceIJ(1, sizeIj, 0).parent(level))
            assertThat(nbrs[2]).isEqualTo(S2CellId.fromFaceIJ(1, 0, sizeIj).parent(level))
            assertThat(nbrs[3]).isEqualTo(S2CellId.fromFaceIJ(0, kMaxIJ, 0).parent(level))
        }

        // Check the vertex neighbors of the center of face 2 at level 5.
        val nbrs = mutableListOf<S2CellId>()
        S2CellId.fromPoint(S2Point(0, 0, 1)).appendVertexNeighbors(5, nbrs)
        nbrs.sort()
        for (i in 0..3) {
            assertThat(nbrs[i]).isEqualTo(
                S2CellId.fromFaceIJ(
                    2,
                    (1 shl 29) - if (i < 2) 1 else 0,
                    (1 shl 29) - if (i == 0 || i == 3) 1 else 0
                ).parent(5),
            )
        }
        nbrs.clear()

        // Check the vertex neighbors of the corner of faces 0, 4, and 5.
        var id = S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel)
        id.appendVertexNeighbors(0, nbrs)
        nbrs.sort()
        assertThat(nbrs.size).isEqualTo(3)
        assertThat(nbrs[0]).isEqualTo(S2CellId.fromFace(0))
        assertThat(nbrs[1]).isEqualTo(S2CellId.fromFace(4))
        assertThat(nbrs[2]).isEqualTo(S2CellId.fromFace(5))

        // Check that AppendAllNeighbors produces results that are consistent
        // with AppendVertexNeighbors for a bunch of random cells.
        for (i in 0 until 1000) {
            id = S2Random.randomCellId()
            if (id.isLeaf) id = id.parent()

            // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell ids,
            // so it's not reasonable to use large values of "diff".
            val maxDiff = min(5, S2CellId.kMaxLevel - id.level() - 1)
            val level = id.level() + S2Random.randomInt(maxDiff + 1)
            testAllNeighbors(id, level)
        }
    }

    private fun testAllNeighbors(id: S2CellId, level: Int) {
        requireGE(level, id.level())
        requireLT(level, S2CellId.kMaxLevel)

        // We compute AppendAllNeighbors, and then add in all the children of "id"
        // at the given level.  We then compare this against the result of finding
        // all the vertex neighbors of all the vertices of children of "id" at the
        // given level.  These should give the same result.
        val all = mutableListOf<S2CellId>()
        val expected = mutableListOf<S2CellId>()
        id.appendAllNeighbors(level, all)
        val end = id.childEnd(level + 1)
        var c = id.childBegin(level + 1)
        while (c != end) {
            all.add(c.parent())
            c.appendVertexNeighbors(level, expected)
            c = c.next()
        }
        // Sort the results and eliminate duplicates.
        assertThat(all.toSortedSet()).isEqualTo(expected.toSortedSet())
    }

    @Test
    fun expandedByDistanceUV() {
        val maxDistDegrees = 10.0
        S2Random.reset(364)
        for (iter in 0 until 100) {
            val id = S2Random.randomCellId()
            val distDegrees = S2Random.randomDouble(-maxDistDegrees, maxDistDegrees)
            testExpandedByDistanceUV(id, S1Angle.degrees(distDegrees))
        }
    }

    private fun testExpandedByDistanceUV(id: S2CellId, distance: S1Angle) {
        val bound = id.boundUV
        val expanded = S2CellId.expandedByDistanceUV(bound, distance)
        for (iter in 0 until 100) {
            // Choose a point on the boundary of the rectangle.
            val face = S2Random.randomInt(6)
            val centerUv = sampleBoundary(bound)
            val center = S2Coords.faceUvToXyz(face, centerUv.x, centerUv.y).normalize()

            // Now sample a point from a disc of radius (2 * distance).
            val p = S2Random.samplePoint(S2Cap.fromCenterAngle(center, 2.0 * distance.abs()))

            // Find the closest point on the boundary to the sampled point.
            val uv = S2Coords.faceXyztoUv(face, p) ?: continue

            val closestUv = projectToBoundary(uv, bound)
            val closest = S2Coords.faceUvToXyz(face, closestUv[0], closestUv[1]).normalize()
            val actualDist = S1Angle(p, closest)

            if (distance >= S1Angle.zero()) {
                // "expanded" should contain all points in the original bound, and also
                // all points within "distance" of the boundary.
                if (bound.contains(uv) || actualDist < distance) {
                    assertThat(expanded.contains(uv)).isTrue()
                }
            } else {
                // "expanded" should not contain any points within "distance" of the
                // original boundary.
                if (actualDist < -distance) {
                    assertThat(expanded.contains(uv)).isFalse()
                }
            }
        }
    }

    @Test
    fun toStringMethod() {
        assertThat(S2CellId.fromFace(3).toString()).isEqualTo("3/")
        assertThat(S2CellId.fromFace(4).rangeMin().toString()).isEqualTo("4/000000000000000000000000000000")
        assertThat(S2CellId.none.toString()).isEqualTo("Invalid: 0000000000000000")
    }

    @Test
    fun fromDebugString() {
        assertThat(S2CellId.fromDebugString("3/")).isEqualTo(S2CellId.fromFace(3))
        assertThat(S2CellId.fromDebugString("0/21")).isEqualTo(S2CellId.fromFace(0).child(2).child(1))
        assertThat(S2CellId.fromDebugString("4/000000000000000000000000000000")).isEqualTo(
            S2CellId.fromFace(4).rangeMin()
        )
        assertThat(S2CellId.fromDebugString("4/0000000000000000000000000000000")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString("")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString("7/")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString(" /")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString("3:0")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString("3/ 12")).isEqualTo(S2CellId.none)
        assertThat(S2CellId.fromDebugString("3/1241")).isEqualTo(S2CellId.none)
    }

    @Test
    fun outputOperator() {
        val cell = S2CellId(0xbb04000000000000UL)
        assertThat(cell.toString()).isEqualTo("5/31200")
    }

    companion object {

        private val logger = Logger.getLogger(S2CellIdUnitTest::class.qualifiedName!!)

        const val kMaxExpandLevel = 3

        const val kMaxWalkLevel = 8

        fun getCellId(lat_degrees: Int, lng_degrees: Int): S2CellId =
            getCellId(lat_degrees.toDouble(), lng_degrees.toDouble())

        fun getCellId(lat_degrees: Double, lng_degrees: Double): S2CellId {
            val id = S2CellId.fromLatLng(S2LatLng.fromDegrees(lat_degrees, lng_degrees))
            logger.info { "CellId from $lat_degrees - $lng_degrees : $id = ${id.id} " }
            return id
        }

        // Returns a random point on the boundary of the given rectangle.
        fun sampleBoundary(rect: R2Rect): R2Point {
            val uv = doubleArrayOf(0.0, 0.0)
            val d = S2Random.randomInt(2)
            uv[d] = S2Random.randomDouble(rect[d][0], rect[d][1])
            uv[1 - d] = if (S2Random.oneIn(2)) rect[1 - d][0] else rect[1 - d][1]
            return R2Point(uv)
        }

        // Returns the closest point to "uv" on the boundary of "rect".
        fun projectToBoundary(uv: R2Point, rect: R2Rect): R2Point {
            val du0 = abs(uv[0] - rect[0][0])
            val du1 = abs(uv[0] - rect[0][1])
            val dv0 = abs(uv[1] - rect[1][0])
            val dv1 = abs(uv[1] - rect[1][1])
            val dmin = min(min(du0, du1), min(dv0, dv1))
            if (du0 == dmin) return R2Point(rect[0][0], rect[1].project(uv[1]))
            if (du1 == dmin) return R2Point(rect[0][1], rect[1].project(uv[1]))
            if (dv0 == dmin) return R2Point(rect[0].project(uv[0]), rect[1][0])
            assertThat(dmin).withFailMessage("Bug in ProjectToBoundary").isEqualTo(dv1)
            return R2Point(rect[0].project(uv[0]), rect[1][1])
        }

    }

}
