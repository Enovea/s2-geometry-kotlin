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

import dilivia.s2.S1Angle
import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.S2CellId
import dilivia.s2.S2Debug
import dilivia.s2.S2Earth
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.coords.S2Coords
import org.assertj.core.api.Assertions
import org.junit.jupiter.api.Test

class S2PolygonInitToSimplifiedInCellUnitTest {

    // Creates a polygon from loops specified as a comma separated list of u:v
    // coordinates relative to a cell. The loop "0:0, 1:0, 1:1, 0:1" is
    // counter-clockwise.
    fun makeCellPolygon(cell: S2Cell, strs: List<String>): S2Polygon {
        val loops = mutableListOf<S2Loop>()
        for (str in strs) {
            val points = S2TextParser.parseLatLngs(str)
            val loop_vertices = mutableListOf<S2Point>()
            val uv = cell.boundUV()
            for (p in points) {
            val u = p.lat().degrees()
                val v = p.lng().degrees()
            loop_vertices.add(
                S2Coords.faceUvToXyz(
                    face = cell.face(),
                    u = uv[0][0] * (1 - u) + uv[0][1] * u,
                    v = uv[1][0] * (1 - v) + uv[1][1] * v
                ).normalize()
            )
        }
            loops.add(S2Loop(loop_vertices))
        }
        return S2Polygon(loops)
    }

    @Test
    fun pointsOnCellBoundaryKept() {
        val cell = S2Cell(S2CellId.fromToken("89c25c"))
        val polygon = makeCellPolygon(cell, listOf("0.1:0, 0.2:0, 0.2:0.5"))
        val tolerance = S1Angle(polygon.loop(0).vertex(0), polygon.loop(0).vertex(1)) * 1.1
        val simplified = S2Polygon()
                simplified.initToSimplified(polygon, IdentitySnapFunction(tolerance))
        Assertions.assertThat(simplified.isEmpty()).isTrue()
        val simplified_in_cell = S2Polygon()
                simplified_in_cell.initToSimplifiedInCell(polygon, cell, tolerance)
        Assertions.assertThat(simplified_in_cell.boundaryEquals(polygon)).isTrue()
        Assertions.assertThat(simplified_in_cell.numVertices()).isEqualTo(3)
        Assertions.assertThat(simplified.getSnapLevel()).isEqualTo(-1)
    }

    @Test
    fun pointsInsideCellSimplified() {
        val cell_id = S2CellId.fromToken("89c25c")
        val cell = S2Cell(cell_id)
        val polygon = makeCellPolygon(cell, listOf("0.3:0, 0.4:0, 0.4:0.5, 0.4:0.8, 0.2:0.8"))
        val tolerance = S1Angle(polygon.loop(0).vertex(0), polygon.loop(0).vertex(1)) * 1.1
        val simplified = S2Polygon()
                simplified.initToSimplifiedInCell(polygon, cell, tolerance)
        Assertions.assertThat(simplified.boundaryNear(polygon, S1Angle.radians(1e-15))).isTrue()
        Assertions.assertThat(simplified.numVertices()).isEqualTo(4)
        Assertions.assertThat(simplified.getSnapLevel()).isEqualTo(-1)
    }

    @Test
    fun cellCornerKept() {
        val cell = S2Cell(S2CellId.fromToken("00001"))
        val input = makeCellPolygon(cell, listOf("1:0, 1:0.05, 0.99:0"))
        val tolerance = 0.02 * S1Angle(cell.getVertex(0), cell.getVertex(1))
        val simplified = S2Polygon()
                simplified.initToSimplifiedInCell(input, cell, tolerance)
        Assertions.assertThat(simplified.boundaryNear(input, S1Angle.radians(1e-15))).isTrue()
    }

    @Test
    fun narrowStripRemoved() {
        val cell = S2Cell(S2CellId.fromToken("00001"))
        val input = makeCellPolygon(cell, listOf("0.9:0, 0.91:0, 0.91:1, 0.9:1"))
        val tolerance = 0.02 * S1Angle(cell.getVertex(0), cell.getVertex(1))
        val simplified = S2Polygon()
                simplified.initToSimplifiedInCell(input, cell, tolerance)
        Assertions.assertThat(simplified.isEmpty()).isTrue()
    }

    @Test
    fun narrowGapRemoved() {
        val cell = S2Cell(S2CellId.fromToken("00001"))
        val input = makeCellPolygon(cell, listOf("0.7:0, 0.75:0, 0.75:1, 0.7:1", "0.76:0, 0.8:0, 0.8:1, 0.76:1"))
        val expected = makeCellPolygon(cell, listOf("0.7:0, 0.8:0, 0.8:1, 0.7:1"))
        val tolerance = 0.02 * S1Angle(cell.getVertex(0), cell.getVertex(1))
        val simplified = S2Polygon()
                simplified.initToSimplifiedInCell(input, cell, tolerance)
        Assertions.assertThat(simplified.boundaryNear(expected, S1Angle.radians(1e-15))).isTrue()
    }

    @Test
    fun closelySpacedEdgeVerticesKept() {
        val cell = S2Cell(S2CellId.fromToken("00001"))
        val input = makeCellPolygon(cell, listOf("0:0.303, 0:0.302, 0:0.301, 0:0.3, 0.1:0.3, 0.1:0.4"))
        val tolerance = 0.02 * S1Angle(cell.getVertex(0), cell.getVertex(1))
        val simplified = S2Polygon()
                simplified.initToSimplifiedInCell(input, cell, tolerance)
        Assertions.assertThat(simplified.boundaryApproxEquals(input, S1Angle.radians(1e-15))).isTrue()
    }

    @Test
    fun polylineAssemblyBug() {
        val cell = S2Cell(S2CellId.fromToken("5701"))
        val polygon = makePolygon(
            "55.8699252:-163.9412145, " + // South-west corner of 5701
                    "54.7672352:-166.7579678, " + // North-east corner of 5701
                    /* Offending part: a tiny triangle near south-east corner */
                    "54.7109214:-164.6376338, " + // forced vertex, on edge 4
                    "54.7140193:-164.6398404, " +
                    "54.7113202:-164.6374015"
        )  // forced vertex, on edge 4
        val tolerance = S1Angle.radians(2.138358e-05)  // 136.235m
        val max_dist = S1Angle.radians(2.821947e-09)  // 18mm
        val simplified_in_cell = S2Polygon()
                simplified_in_cell.initToSimplifiedInCell(polygon, cell, tolerance, max_dist)
        Assertions.assertThat(simplified_in_cell.isEmpty()).isFalse()
    }

    @Test
    fun interiorEdgesSnappedToBoundary() {
        val polygon = S2TextParser.makePolygon(
            "37.8011672:-122.3247322, 37.8011648:-122.3247399, " +
                    "37.8011647:-122.3247403, 37.8011646:-122.3247408, " +
                    "37.8011645:-122.3247411, 37.8011633:-122.3247449, " +
                    "37.8011621:-122.3247334"
        )
        val cell = S2Cell(S2CellId.fromDebugString("4/001013300"))
        val snap_radius = S1Angle.radians(S2Earth.metersToRadians(1.0))
        val boundary_tolerance = S1Angle.radians(0.5 * S2Coords.projection.kMaxWidth.getValue(S2CellId.kMaxLevel - 1)) + IntLatLngSnapFunction.minSnapRadiusForExponent(
            7
        )
        val simplified_polygon = S2Polygon()
                simplified_polygon.debugOverride = S2Debug.DISABLE
        simplified_polygon.initToSimplifiedInCell(polygon, cell, snap_radius, boundary_tolerance)
        val error = S2Error()
        Assertions.assertThat(simplified_polygon.findValidationError(error)).isFalse()
    }


}
