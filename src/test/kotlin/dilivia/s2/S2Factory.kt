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

import dilivia.math.M_PI
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import kotlin.math.cos
import kotlin.math.sin

object S2Factory {

    fun parseVertices(str: String): List<S2Point> {
        val vertices: MutableList<S2Point> = mutableListOf()
        if (str != "") {
            for (token in str.split(',')) {
                val colon = token.indexOf(':')
                require(colon != -1) { "Illegal string:$token. Should look like '35:20'" }
                val lat: Double = token.substring(0, colon).toDouble()
                val lng: Double = token.substring(colon + 1).toDouble()
                vertices.add(S2LatLng.fromDegrees(lat, lng).toPoint())
            }
        }
        return vertices
    }

    @JvmStatic
    fun makePoint(str: String): S2Point {
        val vertices = parseVertices(str)
        return vertices.first()
    }

    @JvmStatic
    fun makeLoop(str: String): S2Loop {
        val vertices = parseVertices(str)
        return S2Loop(vertices)
    }

    @JvmStatic
    @JvmOverloads
    fun makePolyline(str: String, debugOverride: S2Debug = S2Debug.ALLOW): S2Polyline {
        val vertices = parseVertices(str)
        return S2Polyline(vertices.toMutableList(), debugOverride)
    }

    @JvmStatic
    fun makeRegularPoints(center: S2Point, radius: S1Angle, num_vertices: Int): List<S2Point> {
        val loop = S2Loop.makeRegularLoop(center, radius, num_vertices)
        val points = ArrayList<S2Point>(loop.numVertices)
        for (i in 0 until loop.numVertices) {
            points.add(loop.vertex(i))
        }
        return points;
    }


    fun concentricLoopsPolygon(center: S2Point, num_loops: Int, num_vertices_per_loop: Int, polygon: S2Polygon) {
        val m = S2PointUtil.getFrame(center)
        val loops = mutableListOf<S2Loop>()
        for (li in 0 until num_loops) {
            val vertices = mutableListOf<S2Point>()
            val radius = 0.005 * (li + 1) / num_loops
            val radian_step = 2 * M_PI / num_vertices_per_loop
            for (vi in 0 until num_vertices_per_loop) {
                val angle = vi * radian_step
                val p = S2Point(radius * cos(angle), radius * sin(angle), 1.0)
                vertices.add(S2PointUtil.fromFrame(m, p.normalize()))
            }
            loops.add(S2Loop(vertices))
        }
        polygon.initNested(loops)
    }

}
