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
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2TextParser
import dilivia.s2.S2TextParser.makeIndex
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2BooleanOperation
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.region.S2CellUnionUnitTest.Companion.checkCovering
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2ShapeIndexBufferedRegionUnitTest {

  private val logger = KotlinLogging.logger {  }

    @Test
    fun emptyIndex() {
        // Test buffering an empty S2ShapeIndex.
        val index = MutableS2ShapeIndex()
        val radius = S1ChordAngle(S1Angle.degrees(2))
        val region = S2ShapeIndexBufferedRegion(index, radius)
        val coverer = S2RegionCoverer()
        val covering = coverer.getCovering(region)
        assertThat(covering.isEmpty()).isTrue()
    }

    @Test
    fun fullPolygon() {
        // Test buffering an S2ShapeIndex that contains a full polygon.
        val index = makeIndex("# # full");
        val radius = S1ChordAngle(S1Angle.degrees(2))
        val region = S2ShapeIndexBufferedRegion(index, radius)
        val coverer = S2RegionCoverer()
        val covering = coverer.getCovering(region)
        assertThat(covering.numCells()).isEqualTo(6)
        for (id in covering) {
            assertThat(id.isFace).isTrue()
        }
    }

    @Test
    fun fullAfterBuffering() {
        // Test a region that becomes the full polygon after buffering.
        val index = makeIndex("0:0 | 0:90 | 0:180 | 0:-90 | 90:0 | -90:0 # #");
        val radius = S1ChordAngle(S1Angle.degrees(60))
        val region = S2ShapeIndexBufferedRegion(index, radius)
        val coverer = S2RegionCoverer(maxCells = 1000)
        val covering = coverer.getCovering(region)
        assertThat(covering.numCells()).isEqualTo(6)
        for (id in covering) {
            assertThat(id.isFace).isTrue();
        }
    }

    @Test
    fun pointZeroRadius() {
        // Test that buffering a point using a zero radius produces a non-empty
        // covering.  (This requires using "less than or equal to" distance tests.)
        val index = makeIndex("34:25 # #");
        val region = S2ShapeIndexBufferedRegion(index, S1ChordAngle.zero())
        val coverer = S2RegionCoverer()
        val covering = coverer.getCovering(region)
        assertThat(covering.numCells()).isEqualTo(1)
        for (id in covering) {
            assertThat(id.isLeaf).isTrue()
        }
    }

    @Test
    fun bufferedPointVsCap() {
        // Compute an S2Cell covering of a buffered S2Point, then make sure that the
        // covering is equivalent to the corresponding S2Cap.
        val index = makeIndex("3:5 # #")
        val point = makePoint ("3:5")
        val radius = S1ChordAngle(S1Angle.degrees(2))
        val region = S2ShapeIndexBufferedRegion(index, radius)
        val coverer = S2RegionCoverer(maxCells = 50)
        val covering = coverer.getCovering(region)
        val equivalentCap = S2Cap(point, radius)
        checkCovering(equivalentCap, covering, true)
    }

    // Verifies that an arbitrary S2ShapeIndex is buffered correctly, by first
    // converting the covering to an S2Polygon and then checking that (a) the
    // S2Polygon contains the original geometry and (b) the distance between the
    // original geometry and the boundary of the S2Polygon is at least "radius".
    //
    // The "radius" parameter is an S1Angle for convenience.
    // TODO(ericv): Add Degrees, Radians, etc, methods to S1ChordAngle?
    fun testBufferIndex(index_str: String, radius_angle: S1Angle, coverer: S2RegionCoverer)
    {
        val index = makeIndex(index_str);
        val radius = S1ChordAngle(radius_angle);
        val region = S2ShapeIndexBufferedRegion(index, radius)
        val covering = coverer.getCovering(region)
      logger.info { "Covering uses ${covering.numCells()} cells vs. max of ${coverer.getMaxCells()}" }
      logger.debug { """
        |
        |--------------------------
        |S2Polygon: 
        |${ covering.joinToString("\n") { S2TextParser.toString(S2Polygon(S2Cell()))}}
        |--------------------------
      """.trimMargin() }
        
        // Compute an S2Polygon representing the union of the cells in the covering.
        val covering_polygon = S2Polygon()
        covering_polygon.initToCellUnionBorder(covering)
        val covering_index = MutableS2ShapeIndex()
        covering_index.add(S2Polygon.Shape(polygon = covering_polygon))

        // (a) Check that the covering contains the original index.
        assertThat(S2BooleanOperation.contains(covering_index, index)).isTrue()

        // (b) Check that the distance between the boundary of the covering and the
        // the original indexed geometry is at least "radius".
        val query = S2ClosestEdgeQuery(covering_index)
        query.options.includeInteriors = false
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(index)
        assertThat(query.isDistanceLess(target, radius)).isFalse()
    }

    @Test
    fun pointSet() {
        // Test buffering a set of points.
        val coverer = S2RegionCoverer(maxCells = 100)
        testBufferIndex("10:20 | 10:23 | 10:26 # #", S1Angle.degrees(5), coverer)
    }

    @Test
    fun polyline() {
        // Test buffering a polyline.
      val coverer = S2RegionCoverer(maxCells = 100)
        testBufferIndex("# 10:5, 20:30, -10:60, -60:100 #", S1Angle.degrees(2), coverer)
    }

    @Test
    fun polygonWithHole() {
        // Test buffering a polygon with a hole.
      val coverer = S2RegionCoverer(maxCells = 100)
        testBufferIndex("# # 10:10, 10:100, 70:0; 11:11, 69:0, 11:99", S1Angle.degrees(2), coverer)
    }

}
