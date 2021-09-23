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
package dilivia.s2.index.shape

import dilivia.math.M_PI
import dilivia.math.M_PI_4
import dilivia.s2.S1Angle
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2TextParser.makeIndex
import dilivia.s2.shape.S2PointVectorShape
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test

//
// Note that the "real" testing of these methods is in s2loop_measures_test
// and s2polyline_measures_test.  This file only checks the handling of
// multiple shapes and shapes of different dimensions.

class S2ShapeIndexMeasuresUnitTest {

    @Test
    fun getDimensionEmpty() {
        assertEquals(-1, makeIndex("# #").getDimension())
    }

    @Test
    fun testGetDimensionPoints() {
        assertEquals(0, makeIndex("0:0 # #").getDimension())

        // Create an index with an empty point set.
        val index = MutableS2ShapeIndex()
        index.add(S2PointVectorShape())
        assertEquals(0, index.getDimension())
    }

    @Test
    fun testGetDimensionPointsAndLines() {
        assertEquals(1, makeIndex("0:0 # 1:1, 1:2 #").getDimension());

        // Note that a polyline with one vertex has no edges, so it is effectively
        // empty for the purpose of testing GetDimension().
        assertEquals(1, makeIndex("0:0 # 1:1 #").getDimension());
    }

    @Test
    fun getDimensionPointsLinesAndPolygons() {
        assertEquals(2, makeIndex("0:0 # 1:1, 2:2 # 3:3, 3:4, 4:3").getDimension());

        assertEquals(2, makeIndex("# # empty").getDimension());
    }

    @Test
    fun getNumPointsEmpty() {
        assertEquals(0, makeIndex("# #").getNumPoints())
    }

    @Test
    fun getNumPointsTwoPoints() {
        assertEquals(2, makeIndex("0:0 | 1:0 # #").getNumPoints());
    }

    @Test
    fun getNumPointsLineAndPolygon() {
        assertEquals(0, makeIndex("# 1:1, 1:2 # 0:3, 0:5, 2:5").getNumPoints())
    }

    @Test
    fun getLengthEmpty() {
        assertEquals(S1Angle.zero(), makeIndex("# #").getLength())
    }

    @Test
    fun getLengthTwoLines() {
        assertEquals(S1Angle.degrees(2), makeIndex("4:4 # 0:0, 1:0 | 1:0, 2:0 # 5:5, 5:6, 6:5").getLength())
    }

    @Test
    fun getPerimeterEmpty() {
        assertEquals(S1Angle.zero(), makeIndex("# #").getPerimeter());
    }

    @Test
    fun getPerimeterDegeneratePolygon() {
        assertEquals(4.0, makeIndex("4:4 # 0:0, 1:0 | 2:0, 3:0 # 0:1, 0:2, 0:3").getPerimeter().degrees(), 1e-15);
    }

    @Test
    fun getAreaEmpty() {
        assertEquals(0.0, makeIndex("# #").getArea())
    }

    @Test
    fun getAreaTwoFullPolygons() {
        assertEquals(8 * M_PI, makeIndex("# # full | full").getArea())
    }

    @Test
    fun getApproxAreaEmpty() {
        assertEquals(0.0, makeIndex("# #").getApproxArea())
    }

    @Test
    fun getApproxAreaTwoFullPolygons() {
        assertEquals(8 * M_PI, makeIndex("# # full | full").getApproxArea())
    }

    @Test
    fun GetCentroidEmpty() {
        assertEquals(S2Point(0, 0, 0), makeIndex("# #").getCentroid())
    }

    @Test
    fun getCentroidPoints() {
        assertEquals(S2Point(1, 1, 0), makeIndex("0:0 | 0:90 # #").getCentroid());
    }

    @Test
    fun getCentroidPolyline() {
        // Checks that points are ignored when computing the centroid.
        assertTrue(S2PointUtil.approxEquals(
                S2Point(1, 1, 0),
            makeIndex("5:5 | 6:6 # 0:0, 0:90 #").getCentroid()))
    }

    @Test
    fun getCentroidPolygon() {
        // Checks that points and polylines are ignored when computing the centroid.
        assertTrue(
            S2PointUtil.approxEquals(
                S2Point(M_PI_4, M_PI_4, M_PI_4),
                makeIndex("5:5 # 6:6, 7:7 # 0:0, 0:90, 90:0").getCentroid()))
    }

}  // namespace
