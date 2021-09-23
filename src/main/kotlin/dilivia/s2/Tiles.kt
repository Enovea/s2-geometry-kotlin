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

import dilivia.s2.S1Angle.Companion.times
import dilivia.s2.region.S2LatLngRect
import kotlin.math.ceil

/**
 * A class that provides a uniform (square) tiling system for a specified bounding box and tile size.
 * A unique tile Id is assigned for each tile based on the following rules:
 *    Tile numbers start at 0 at the min y, x (lower left)
 *    Tile numbers increase by column (x,longitude) then by row (y,latitude)
 *    Tile numbers increase along each row by increasing x,longitude.
 * Contains methods for converting x,y or lat,lng into tile Id and
 * vice-versa.  Methods for relative tiles (using row and column offsets).
 * are also provided. Also includes a method to get a list of tiles covering
 * a bounding box.
 *
 * @constructor A bottom left coord, with tile_size and number of rows and columns
 * @param   bottomLeftCorner       Bottom left coord of the tileset
 * @param   tileSize    The size of a tile in both dimensions
 * @param   columns      Number of tiles in the x axis
 * @param   rows         Number of tiles in the y axis
 *
 * @author fmeurisse
 */
class Tiles(bottomLeftCorner: S2LatLng, tileSize: S1Angle, columns: Int, rows: Int) {
    /** Bounding box of the tiling system. */
    val bounds: S2LatLngRect = S2LatLngRect.fromPointPair(
        p1 = bottomLeftCorner,
        p2 = bottomLeftCorner.plus(S2LatLng.fromDegrees(rows * tileSize.degrees(), columns * tileSize.degrees()))
    )

    /** Tile size.  Tiles are square (equal y and x size). */
    val tileSize: S1Angle = tileSize

    /** Number of rows (latitude) */
    val nrows: Int = rows

    /** Number of longitude (x or longitude). */
    val ncolumns: Int = columns

    /**
     * Get the "row" based on latitude.
     * @param   latitude   latitude coordinate
     * @return  Returns the tile row. Returns -1 if outside the tile system bounds.
     */
    fun row(latitude: S1Angle): Int {
        latitude.normalize()
        // Return -1 if outside the tile system bounds
        if (!this.bounds.lat.contains(latitude.radians)) {
            return -1
        }

        // If equal to the max y return the largest row
        val latHi = this.bounds.latHi()
        val latLo = this.bounds.latLo()
        return if (latitude == latHi) nrows - 1 else ((latitude - latLo) / tileSize).radians.toInt()
    }

    /**
     * Get the "column" based on longitude.
     * @param   longitude   longitude coordinate
     * @return  Returns the tile column. Returns -1 if outside the tile system bounds.
     */
    fun column(longitude: S1Angle): Int {
        longitude.normalize()
        // Return -1 if outside the tile system bounds
        if (!this.bounds.lng.contains(longitude.radians)) {
            return -1;
        }

        // If equal to the max x return the largest column
        val lngHi = this.bounds.lngHi()
        val lngLo = this.bounds.lngLo()
        return if (longitude == lngHi) ncolumns - 1 else {
            val col = ((longitude - lngLo) / tileSize).radians
            if (col >= 0.0) col.toInt() else (col.toInt() - 1)
        }
    }

    /**
     * Convert latitude,longitude to a tile Index.
     * @param   latitude   Latitude of the point.
     * @param   longitude  Longitude of the point
     * @return  Returns the tile Index. -1 (error is returned if the x,y is
     *          outside the bounding box of the tiling sytem).
     */
    fun tileIndex(latitude: S1Angle, longitude: S1Angle): Int {
        latitude.normalize()
        longitude.normalize()
        // Return -1 if totally outside the extent.
        if (!this.bounds.lat.contains(latitude.radians) || !this.bounds.lng.contains(longitude.radians)) {
            return -1
        }

        // Find the tileid by finding the latitude row and longitude column
        return tileIndex(column(longitude), row(latitude))
    }

    fun tileIndex(point: S2LatLng): Int = tileIndex(point.lat(), point.lng())

    fun tileIndex(point: S2Point): Int = tileIndex(S2LatLng.fromPoint(point))

    fun tileIndex(col: Int, row: Int): Int = (row * ncolumns) + col

    /**
     * Get the tile row, col based on tile Id.
     * @param  tileid  Tile Id.
     * @return  Returns a pair indicating {row, col}
     */
    fun tileRowColumn(index: Int): Pair<Int, Int> = (index / ncolumns) to (index % ncolumns)

    fun tileBottomLeftCorner(index: Int): S2LatLng {
        val row = index / ncolumns
        val col = index - (row * ncolumns)
        return S2LatLng.fromLatLng(this.bounds.latLo() + row * tileSize, this.bounds.lngLo() + col * tileSize)
    }

    fun tileBounds(index: Int): S2LatLngRect {
        val bottomLeftCorner = tileBottomLeftCorner(index)
        return S2LatLngRect(bottomLeftCorner, bottomLeftCorner + S2LatLng.fromLatLng(tileSize, tileSize))
    }

    fun tileCenter(index: Int): S2LatLng {
        return tileBounds(index).center
    }

    fun tileCount(): Int = nrows * ncolumns

    fun tileList(box: S2LatLngRect): List<Int> {
        val tiles = mutableListOf<Int>()
        tileList(box, tiles)
        return tiles
    }

    fun tileList(box: S2LatLngRect, tiles: MutableList<Int>) {
        tiles.clear()

    }

    companion object {
        /**
         * Get a maximum tileid given a bounds and a tile size.
         * @param bounds      the region for which to compute the maximum tile id
         * @param tile_size   the size of a tile within the region
         * @return the highest tile number within the region
         */
        fun maxTileIndex(bounds: S2LatLngRect, tileSize: S1Angle): Int {
            val cols = ceil(bounds.lng.length / tileSize.radians).toInt()
            val rows = ceil(bounds.lat.length / tileSize.radians).toInt()
            return (cols * rows) - 1
        }


    }
}
