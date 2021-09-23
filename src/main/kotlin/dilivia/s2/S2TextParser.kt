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

import dilivia.PreConditions.requireEQ
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2LatLngRect
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.S2LaxPolygonShape
import dilivia.s2.shape.S2LaxPolylineShape
import dilivia.s2.shape.S2PointVectorShape
import dilivia.s2.shape.S2Shape
import mu.KotlinLogging
import java.math.RoundingMode
import java.text.DecimalFormat
import java.text.DecimalFormatSymbols
import java.util.*

object S2TextParser {

    private val logger = KotlinLogging.logger { }

    // Parses a string of one or more latitude-longitude coordinates in degrees,
    // and return the corresponding vector of S2LatLng points.
    // Examples of the input format:
    //     ""                            // no points
    //     "-20:150"                     // one point
    //     "-20:150, -20:151, -19:150"   // three points
    fun parseLatLngs(str: String): List<S2LatLng> {
        val latlngs = mutableListOf<S2LatLng>()
        check(parseLatLngs(str, latlngs)) { ": str == \"$str\"" }
        return latlngs
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun parseLatLngs(str: String, latlngs: MutableList<S2LatLng>): Boolean {
        val ps = str.split(',')
                .filter { pointStr -> pointStr.isNotBlank() }
                .map { pointStr -> pointStr.trim() }
                .map { pointStr ->
                    val coordsStr = pointStr.split(':')
                    if (coordsStr.size != 2) return false
                    Pair(coordsStr[0].trim(), coordsStr[1].trim())
                }

        for (p in ps) {
            val lat = p.first.toDoubleOrNull() ?: return false
            val lng = p.second.toDoubleOrNull() ?: return false
            latlngs.add(S2LatLng.fromDegrees(lat, lng))
        }
        return true
    }

    // Parses a string in the same format as ParseLatLngs, and return the
    // corresponding vector of S2Point values.
    fun parsePoints(str: String): List<S2Point> {
        val vertices = mutableListOf<S2Point>()
        check(parsePoints(str, vertices)) { ": str == \"$str\"" }
        return vertices
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun parsePoints(str: String, vertices: MutableList<S2Point>): Boolean {
        val latlngs = mutableListOf<S2LatLng>()
        if (!parseLatLngs(str, latlngs)) return false
        for (latlng in latlngs) {
            vertices.add(latlng.toPoint())
        }
        return true
    }

    fun makePoint(str: String): S2Point {
        val vertices = mutableListOf<S2Point>()
        check(parsePoints(str, vertices) && vertices.size == 1) { ": str == \"$str\"" }
        return vertices[0];
    }

    fun makeLatLngRect(str: String): S2LatLngRect {
        val latlngs = mutableListOf<S2LatLng>()
        check(parseLatLngs(str, latlngs) && !latlngs.isEmpty())
        var rect = S2LatLngRect.fromPoint(latlngs[0])
        for (i in 1 until latlngs.size) {
            rect = rect.addPoint(latlngs[i])
        }
        return rect
    }

    fun makeLoop(str: String, debug: S2Debug = S2Debug.ALLOW): S2Loop {
        return when (str) {
            "empty" -> S2Loop(S2Loop.kEmpty)
            "full" -> S2Loop(S2Loop.kFull)
            else -> {
                val vertices = mutableListOf<S2Point>()
                check(parsePoints(str, vertices)) { ": str == \"$str\"" }
                S2Loop(vertices, debugOverride = debug)
            }
        }
    }

    // Like MakePolyline, but returns an S2LaxPolylineShape instead.
    fun makeLaxPolyline(str: String): S2LaxPolylineShape {
        val vertices = mutableListOf<S2Point>()
        check(parsePoints(str, vertices)) { ": str == \"$str\"" }
        return S2LaxPolylineShape(vertices)
    }

    // Parses an S2CellId in the format "f/dd..d" where "f" is a digit in the
    // range [0-5] representing the S2CellId face, and "dd..d" is a string of
    // digits in the range [0-3] representing each child's position with respect
    // to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents S2CellId::FromFace(4), and "3/02" represents
    // S2CellId::FromFace(3).child(0).child(2).
    //
    // This function is a wrapper for S2CellId::FromDebugString().
    fun makeCellId(str: String): S2CellId {
        val cellId = S2CellId.fromDebugString(str)
        check(cellId != S2CellId.none) { "Invalid cell id: str == \"$str\"" }
        return cellId;
    }

    // Parses a comma-separated list of S2CellIds in the format above, and returns
    // the corresponding S2CellUnion.  (Note that S2CellUnions are automatically
    // normalized by sorting, removing duplicates, and replacing groups of 4 child
    // cells by their parent cell.)
    fun makeCellUnion(str: String): S2CellUnion {
        val cellIds = mutableListOf<S2CellId>()
        str.split(",").map { it.trim() }.forEach { cellStr -> cellIds.add(makeCellId(cellStr)) }
        return S2CellUnion(cellIds)
    }


    // Returns a MutableS2ShapeIndex containing the points, polylines, and loops
    // (in the form of a single polygon) described by the following format:
    //
    //   point1|point2|... # line1|line2|... # polygon1|polygon2|...
    //
    // Examples:
    //   1:2 | 2:3 # #                     // Two points
    //   # 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
    //   # # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
    //   5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each
    //   # # empty                         // One empty polygon
    //   # # empty | full                  // One empty polygon, one full polygon
    //
    // Loops should be directed so that the region's interior is on the left.
    // Loops can be degenerate (they do not need to meet S2Loop requirements).
    //
    // CAVEAT: Because whitespace is ignored, empty polygons must be specified
    //         as the string "empty" rather than as the empty string ("").
    fun makeIndex(str: String): MutableS2ShapeIndex {
        val index = MutableS2ShapeIndex()
        check(makeIndex(str, index)) { ": str == \"$str\"" }
        return index;
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun makeIndex(str: String, index: MutableS2ShapeIndex): Boolean {
        val strs = str.split('#').map { it.trim() }
        check(3 == strs.size) { "Must contain two # characters: $str" }

        val points = mutableListOf<S2Point>()
        for (point_str in strs[0].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val point = makePoint(point_str)
            points.add(point)
        }
        if (points.isNotEmpty()) {
            index.add(S2PointVectorShape(points = points))
        }

        for (line_str in strs[1].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val laxPolyline = makeLaxPolyline(line_str)
            index.add(laxPolyline)
        }

        for (polygon_str in strs[2].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val laxPolygon = makeLaxPolygon(polygon_str)
            index.add(laxPolygon)
        }
        return true
    }

   fun makePolygon(str: String, debug: S2Debug = S2Debug.ALLOW): S2Polygon = internalMakePolygon(str, debug, true)

    fun makeLaxPolygon(str: String): S2LaxPolygonShape {
        val loopStrs = str.trim().let {
            var s = it
            while (s.endsWith(";")) s = s.removeSuffix(";").trim()
            s
        }.split(';')
        val loops = mutableListOf<MutableList<S2Point>>()
        for (loop_str in loopStrs) {
            if (loop_str == "full") {
                loops.add(mutableListOf<S2Point>())
            } else if (loop_str != "empty") {
                val points = mutableListOf<S2Point>()
                check(parsePoints(loop_str, points)) { "Fail to parse loop: $loop_str" }
                if (points.isNotEmpty()) {
                    loops.add(points)
                }
            }
        }
        return S2LaxPolygonShape(loops)
    }


    fun makeVerbatimPolygon(str: String): S2Polygon {
        logger.debug { "Create verbatim Polygon: $str" }
        return internalMakePolygon(str, S2Debug.ALLOW, false)
    }

    private fun internalMakePolygon(str: String, debug: S2Debug, normalizeLoops: Boolean): S2Polygon {
        val loopStrs = str.trim().let {
            var s = it
            while (s.endsWith(";")) s = s.removeSuffix(";").trim()
            s
        }.split(';').map { it.trim() }.filter { it.isNotBlank() }
        val loops = mutableListOf<S2Loop>()
        for (loop_str in loopStrs) {
            val s = loop_str.trim()
            val loop = makeLoop(s, debug)
            // Don't normalize loops that were explicitly specified as "full".
            if (normalizeLoops && !loop.isFull()) loop.normalize()
            loops.add(loop)
            logger.trace { "Internal Make Polygon: create loop $loop_str => ${loop.toDebugString()}" }
        }
        return S2Polygon(loops, debug)
    }

    fun toString(polygon: S2Polygon, loopSeparator: String = " ; "): String {
        if (polygon.isEmpty()) {
            return "empty"
        } else if (polygon.isFull()) {
            return "full"
        }
        var out = ""
        for (i in 0 until polygon.numLoops()) {
            if (i > 0) out += loopSeparator
            val loop = polygon.loop(i)
            out += loop.vertices().joinToString(", ") { p -> toString(p) }
        }
        return out
    }

    fun toString(loop: S2Loop): String {
        if (loop.isEmpty()) {
            return "empty"
        } else if (loop.isFull()) {
            return "full"
        }
        var out = ""
        if (loop.numVertices > 0) {
            out += appendVertices(loop.vertices(), loop.numVertices, out)
            //out += loop.vertices().joinToString(", ", prefix = "[", postfix = "]")
        }
        return out
    }

    fun toString(polyline: S2Polyline): String {
        var out = ""
        if (polyline.numVertices > 0) {
            out += appendVertices(polyline.vertices(), polyline.numVertices, out)
        }
        return out
    }

    fun toString(ll: S2LatLng): String {
        val df = DecimalFormat("#.##############", DecimalFormatSymbols.getInstance(Locale.ENGLISH))
        df.roundingMode = RoundingMode.HALF_UP
        val lat = df.format(ll.latDegrees())
        val lng = df.format(ll.lngDegrees())
        return "$lat:$lng"
    }

    fun toString(p: S2Point): String {
        val ll = S2LatLng.fromPoint(p)
        return toString(ll)
    }

    fun toString(points: List<S2Point>): String {
        var out = ""
        points.forEachIndexed { index, point ->
            if (index > 0) out += ", "
            out += toString(point)
        }
        return out
    }


    fun toString(polygon: S2LaxPolygonShape, loopSeparator: String = "|"): String {
        var out = ""
        for (i in 0 until polygon.numLoops()) {
            if (i > 0) out += loopSeparator
            val n = polygon.numLoopVertices(i)
            if (n == 0) {
                out += "full"
            } else {
                out += toString((0 until n).map { j -> polygon.loopVertex(i, j) })
            }
        }
        return out
    }


    fun toString(index: S2ShapeIndex): String {
        var out = ""
        for (dim in 0..2) {
            if (dim > 0) out += "#"
            var count = 0
            for (shape: S2Shape? in index) {
                if (shape == null || shape.dimension != dim) continue
                out += if (count > 0) " | " else if (dim > 0) " " else ""
                for (i in 0 until shape.numChains) {
                    if (i > 0) out += if (dim == 2) "; " else " | "
                    val chain = shape.chain(i)
                    if (chain.length == 0) {
                        requireEQ(dim, 2)
                        out += "full"
                    } else {
                        out = appendVertex(shape.edge(chain.start).v0, out)
                    }
                    var limit = chain.start + chain.length
                    if (dim != 1) --limit
                    for (e in chain.start until limit) {
                        out += ", ";
                        out = appendVertex(shape.edge(e).v1, out)
                    }
                    ++count
                }
            }
            // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
            if (dim == 1 || (dim == 0 && count > 0)) out += " ";
        }
        return out
    }

    fun appendVertices(v: List<S2Point>, n: Int, out: String): String {
        var result = out
        for (i in 0 until n) {
            if (i > 0) result += ", "
            result += appendVertex(v[i], out)
        }
        return result
    }

    fun appendVertex(p: S2Point, out: String): String {
        val ll = S2LatLng.fromPoint(p)
        return appendVertex(ll, out)
    }


    fun appendVertex(ll: S2LatLng, out: String): String = out + toString(ll)

}
