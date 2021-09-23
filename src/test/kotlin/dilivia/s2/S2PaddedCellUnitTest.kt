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

import dilivia.math.R1Interval
import dilivia.math.R2Rect
import dilivia.math.vectors.R2Point
import dilivia.s2.region.S2Cell
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import kotlin.math.pow

class S2PaddedCellUnitTest {

    private val logger = KotlinLogging.logger {  }

    private fun compareS2CellToPadded(cell: S2Cell, pcell: S2PaddedCell, padding: Double) {
        logger.debug { "Compare S2Cell $cell with padded (padding=$padding) $pcell" }
        assertThat(pcell.id).isEqualTo(cell.id())
        assertThat(pcell.level).isEqualTo(cell.level())
        assertThat(pcell.padding).isEqualTo(padding)
        logger.trace { "Cell bound          = ${cell.boundUV()}, ${S2PaddedCell(cell.id(), 0.0).bound}" }
        logger.trace { "Expanded cell bound = ${cell.boundUV().expanded(padding)}" }
        logger.trace { "Padded cell bound   = ${pcell.bound}" }
        assertThat(pcell.bound).isEqualTo(cell.boundUV().expanded(padding))
        val centerUv = cell.id().centerUV
        assertThat(pcell.middle()).isEqualTo(R2Rect.fromPoint(centerUv).expanded(padding))
        assertThat(pcell.getCenter()).isEqualTo(cell.getCenter())
    }

    @Test
    fun s2CellMethods() {
        // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
        val kIters = 1000
        repeat(kIters) { iter ->
            logger.trace { "iteration ${iter + 1}" }
            val id = S2Random.randomCellId()
            val padding = 1e-15.pow(S2Random.randomDouble())
            val cell = S2Cell(id)
            val pcell = S2PaddedCell(id, padding)
            //compareS2CellToPadded(cell, pcell, padding)
            if (!id.isLeaf) {
                val children = Array(4) { S2Cell() }
                assertThat(cell.subdivide(children)).isTrue()
                for (pos in 0..3) {
                  val (i, j) = pcell.getChildIJ(pos)
                    compareS2CellToPadded(children[pos], S2PaddedCell(pcell, i, j), padding)
                }
            }
        }
    }

    @Test
    fun getEntryExitVertices() {
        val kIters = 1000
        repeat(kIters) {
            val id = S2Random.randomCellId()
            // Check that entry/exit vertices do not depend on padding.
            assertThat(S2PaddedCell(id, 0.0).getEntryVertex()).isEqualTo(S2PaddedCell(id, 0.5).getEntryVertex())
            assertThat(S2PaddedCell(id, 0.0).getExitVertex()).isEqualTo(S2PaddedCell(id, 0.5).getExitVertex())

            // Check that the exit vertex of one cell is the same as the entry vertex
            // of the immediately following cell.  (This also tests wrapping from the
            // end to the start of the S2CellId curve with high probability.)
            assertThat(S2PaddedCell(id, 0.0).getExitVertex()).isEqualTo(S2PaddedCell(id.nextWrap(), 0.0).getEntryVertex())

            // Check that the entry vertex of a cell is the same as the entry vertex
            // of its first child, and similarly for the exit vertex.
            if (!id.isLeaf) {
                assertThat(S2PaddedCell(id, 0.0).getEntryVertex()).isEqualTo(S2PaddedCell(id.child(0), 0.0).getEntryVertex())
                assertThat(S2PaddedCell(id, 0.0).getExitVertex()).isEqualTo(S2PaddedCell(id.child(3), 0.0).getExitVertex())
            }
        }
    }

    private fun sampleInterval(x: R1Interval): Double {
        return S2Random.randomDouble(x.lo, x.hi)
    }

    @Test
    fun shrinkToFit() {
        val kIters = 1000
        repeat(kIters) {
            // Start with the desired result and work backwards.
            val result = S2Random.randomCellId()
            val resultUv = result.boundUV
            val sizeUv = resultUv.size

            // Find the biggest rectangle that fits in "result" after padding.
            // (These calculations ignore numerical errors.)
            val maxPadding = 0.5 * kotlin.math.min(sizeUv[0], sizeUv[1])
            val padding = maxPadding * S2Random.randomDouble()
            val maxRect = resultUv.expanded(-padding)

            // Start with a random subset of the maximum rectangle.
            val a = R2Point(sampleInterval(maxRect[0]), sampleInterval(maxRect[1]))
            val b = R2Point(sampleInterval(maxRect[0]), sampleInterval(maxRect[1]))
            if (!result.isLeaf) {
                // If the result is not a leaf cell, we must ensure that no child of
                // "result" also satisfies the conditions of ShrinkToFit().  We do this
                // by ensuring that "rect" intersects at least two children of "result"
                // (after padding).
                val axis = S2Random.randomInt(2)
                val center = result.centerUV[axis]

                // Find the range of coordinates that are shared between child cells
                // along that axis.
                val shared = R1Interval(center - padding, center + padding)
                val mid = sampleInterval(shared.intersection(maxRect[axis]))
                a[axis] = sampleInterval(R1Interval(maxRect[axis].lo, mid))
                b[axis] = sampleInterval(R1Interval(mid, maxRect[axis].hi))
            }
            val rect = R2Rect.fromPointPair(a, b)

            // Choose an arbitrary ancestor as the S2PaddedCell.
            val initialId = result.parent(S2Random.randomInt(result.level() + 1))
            assertThat(S2PaddedCell(initialId, padding).shrinkToFit(rect)).isEqualTo(result)
        }
    }

}
