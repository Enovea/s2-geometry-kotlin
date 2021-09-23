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
package dilivia.s2.coords

import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.S2Random.randomCellId
import dilivia.s2.S2Random.randomInt
import dilivia.s2.coords.S2Coords.kIJtoPos
import dilivia.s2.coords.S2Coords.kInvertMask
import dilivia.s2.coords.S2Coords.kPosToIJ
import dilivia.s2.coords.S2Coords.kSwapMask
import org.apache.commons.math3.util.FastMath.abs
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test

class S2CoordsUnitTest {

    private fun swapAxes(ij: Int): Int = ((ij shr 1) and 1) + ((ij and 1) shl 1)

    private fun invertBits(ij: Int): Int = ij xor 3

    @Test
    fun traversalOrder() {
        for (r in 0..3) {
            for (i in 0..3) {
                // Check consistency with respect to swapping axes.
                assertThat(kIJtoPos[r xor kSwapMask][swapAxes(i)]).isEqualTo(kIJtoPos[r][i])
                assertThat(swapAxes(kPosToIJ[r xor kSwapMask][i])).isEqualTo(kPosToIJ[r][i])

                // Check consistency with respect to reversing axis directions.
                assertThat(kIJtoPos[r xor kInvertMask][invertBits(i)]).isEqualTo(kIJtoPos[r][i])
                assertThat(invertBits(kPosToIJ[r xor kInvertMask][i])).isEqualTo(kPosToIJ[r][i])

                // Check that the two tables are inverses of each other.
                assertThat(kIJtoPos[r][kPosToIJ[r][i]]).isEqualTo(i)
                assertThat(kPosToIJ[r][kIJtoPos[r][i]]).isEqualTo(i)
            }
        }
    }

    @Test
    fun stUvConversions() {
        // Check boundary conditions.
        var s = 0.0
        while (s <= 1) {
            val u = S2Coords.stToUv(s)
            assertThat(u).isEqualTo(2 * s - 1)
            s += 0.5
        }
        for (u in -1..1) {
            s = S2Coords.uvToSt(u.toDouble())
            assertThat(s).isEqualTo(0.5 * (u + 1))
        }
        // Check that UVtoST and STtoUV are inverses.
        var x = 0.0
        while (x <= 1.0) {
            assertThat(S2Coords.uvToSt(S2Coords.stToUv(x))).isCloseTo(x, Offset.offset(1e-15))
            assertThat(S2Coords.stToUv(S2Coords.uvToSt(2 * x - 1))).isCloseTo(2 * x - 1, Offset.offset(1e-15))
            x += 0.0001
        }
    }

    @Test
    fun faceUVtoXYZ() {
        // Check that each face appears exactly once.
        var sum = S2Point()
        for (face in 0..5) {
            val center = S2Coords.faceUvToXyz(face, 0.0, 0.0)
            assertThat(center).isEqualTo(S2Coords.norm(face))
            assertThat(abs(center[center.largestAbsComponent()])).isEqualTo(1.0)
            sum = sum + center.abs()
        }
        assertThat(sum).isEqualTo(S2Point(2, 2, 2))

        // Check that each face has a right-handed coordinate system.
        for (face in 0..5) {
            assertThat(
                S2Coords.uAxis(face)
                    .crossProd(S2Coords.vAxis(face))
                    .dotProd(S2Coords.faceUvToXyz(face, 0.0, 0.0))
            ).isEqualTo(1.0)
        }

        // Check that the Hilbert curves on each face combine to form a
        // continuous curve over the entire cube.
        for (face in 0..5) {
            // The Hilbert curve on each face starts at (-1,-1) and terminates
            // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
            val sign = if (face and kSwapMask != 0) -1.0 else 1.0
            assertThat(S2Coords.faceUvToXyz(face, sign, -sign)).isEqualTo(
                S2Coords.faceUvToXyz(
                    (face + 1) % 6,
                    -1.0,
                    -1.0
                )
            )
        }
    }

    @Test
    fun faceXYZtoUVW() {
        for (face in 0..5) {
            assertThat(S2Coords.faceXyzToUvw(face, S2Point(0, 0, 0))).isEqualTo(S2Point(0, 0, 0))
            assertThat(S2Coords.faceXyzToUvw(face, S2Coords.uAxis(face))).isEqualTo(S2Point(1, 0, 0))
            assertThat(S2Coords.faceXyzToUvw(face, -S2Coords.uAxis(face))).isEqualTo(S2Point(-1, 0, 0))
            assertThat(S2Coords.faceXyzToUvw(face, S2Coords.vAxis(face))).isEqualTo(S2Point(0, 1, 0))
            assertThat(S2Coords.faceXyzToUvw(face, -S2Coords.vAxis(face))).isEqualTo(S2Point(0, -1, 0))
            assertThat(S2Coords.faceXyzToUvw(face, S2Coords.norm(face))).isEqualTo(S2Point(0, 0, 1))
            assertThat(S2Coords.faceXyzToUvw(face, -S2Coords.norm(face))).isEqualTo(S2Point(0, 0, -1))
        }
    }

    @Test
    fun xyzToFaceSiTi() {
        // Check the conversion of random cells to center points and back.
        for (level in 0..S2CellId.kMaxLevel) {
            for (i in 0 until 1000) {
                val id = randomCellId(level)

                val (actual_level, faceSiTi) = S2Coords.xyzToFaceSiTi(id.toPoint())
                assertThat(actual_level).isEqualTo(level)
                val actual_id =
                    S2CellId.fromFaceIJ(faceSiTi.face, (faceSiTi.si / 2U).toInt(), (faceSiTi.ti / 2U).toInt())
                        .parent(level)
                assertThat(actual_id).isEqualTo(id)

                // Now test a point near the cell center but not equal to it.
                val p_moved = id.toPoint() + S2Point(1e-13, 1e-13, 1e-13)
                val (actual_level_moved, faceSiTi_moved) = S2Coords.xyzToFaceSiTi(p_moved)
                assertThat(actual_level_moved).isEqualTo(-1)
                assertThat(faceSiTi_moved).isEqualTo(faceSiTi)

                // Finally, test some random (si,ti) values that may be at different
                // levels, or not at a valid level at all (for example, si == 0).
                val face_random = randomInt(S2CellId.kNumFaces)
                var si_random = 0U
                var ti_random = 0U
                val mask = -1 shl (S2CellId.kMaxLevel - level)
                do {
                    si_random = (randomInt() and mask).toUInt()
                    ti_random = (randomInt() and mask).toUInt()
                } while (si_random > S2Coords.kMaxSiTi || ti_random > S2Coords.kMaxSiTi)
                val p_random = S2Coords.faceSiTiToXyz(face_random, si_random, ti_random)
                val (actual_level_random, faceSiTi_random) = S2Coords.xyzToFaceSiTi(p_random)
                if (faceSiTi_random.face != face_random) {
                    // The chosen point is on the edge of a top-level face cell.
                    assertThat(actual_level_random).isEqualTo(-1)
                    assertThat(
                        faceSiTi_random.si == 0U || faceSiTi_random.si == S2Coords.kMaxSiTi
                                || faceSiTi_random.ti == 0U || faceSiTi_random.ti == S2Coords.kMaxSiTi
                    ).isTrue()
                } else {
                    assertThat(faceSiTi_random.si).isEqualTo(si_random)
                    assertThat(faceSiTi_random.ti).isEqualTo(ti_random)
                    if (actual_level_random >= 0) {
                        assertThat(
                            S2CellId.fromFaceIJ(
                                faceSiTi_random.face,
                                (faceSiTi_random.si / 2U).toInt(),
                                (faceSiTi_random.ti / 2U).toInt()
                            ).parent(actual_level_random).toPoint()
                        ).isEqualTo(p_random)
                    }
                }
            }
        }
    }

    @Test
    fun uvNorms() {
        // Check that GetUNorm and GetVNorm compute right-handed normals for
        // an edge in the increasing U or V direction.
        for (face in 0..5) {
            var x = -1.0
            while (x <= 1.0) {
                assertThat(
                    S2Coords.faceUvToXyz(face, x, -1.0).crossProd(S2Coords.faceUvToXyz(face, x, 1.0))
                        .angle(S2Coords.getUNorm(face, x))
                ).isZero()
                assertThat(
                    S2Coords.faceUvToXyz(face, -1.0, x).crossProd(S2Coords.faceUvToXyz(face, 1.0, x))
                        .angle(S2Coords.getVNorm(face, x))
                ).isZero()
                x += 1.0 / 1024.0
            }
        }
    }

    @Test
    fun uvwAxis() {
        for (face in 0..5) {
            // Check that axes are consistent with FaceUVtoXYZ.
            assertThat(S2Coords.uAxis(face)).isEqualTo(
                S2Coords.faceUvToXyz(face, 1.0, 0.0) - S2Coords.faceUvToXyz(face, 0.0, 0.0)
            )
            assertThat(S2Coords.vAxis(face)).isEqualTo(
                S2Coords.faceUvToXyz(face, 0.0, 1.0) - S2Coords.faceUvToXyz(face, 0.0, 0.0)
            )
            assertThat(S2Coords.norm(face)).isEqualTo(S2Coords.faceUvToXyz(face, 0.0, 0.0))

            // Check that every face coordinate frame is right-handed.
            assertThat(
                S2Coords.uAxis(face).crossProd(S2Coords.vAxis(face)).dotProd(S2Coords.norm(face))
            ).isEqualTo(1.0)

            // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
            assertThat(S2Coords.uvwAxis(face, 0)).isEqualTo(S2Coords.uAxis(face))
            assertThat(S2Coords.uvwAxis(face, 1)).isEqualTo(S2Coords.vAxis(face))
            assertThat(S2Coords.uvwAxis(face, 2)).isEqualTo(S2Coords.norm(face))
        }
    }

    @Test
    fun uvwFace() {
        // Check that GetUVWFace is consistent with GetUVWAxis.
        for (face in 0..5) {
            for (axis in 0..2) {
                assertThat(S2Coords.uvwFace(face, axis, 0)).isEqualTo(S2Coords.face(-S2Coords.uvwAxis(face, axis)))
                assertThat(S2Coords.uvwFace(face, axis, 1)).isEqualTo(S2Coords.face(S2Coords.uvwAxis(face, axis)) )
            }
        }
    }

}
