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

import dilivia.PreConditions.checkEQ
import dilivia.math.DoubleType
import dilivia.math.M_1_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.R1Interval
import dilivia.math.R2Rect
import dilivia.s2.*
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2Coords.kPosToIJ
import dilivia.s2.coords.S2Coords.kPosToOrientation
import dilivia.s2.edge.S2EdgeCrosser
import dilivia.s2.edge.S2EdgeDistances
import org.apache.commons.math3.util.FastMath.asin
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it supports efficient containment and
 * intersection tests.  However, it is also a more expensive representation (currently 48 bytes rather than 8).
 *
 * @constructor An S2Cell always corresponds to a particular S2CellId. The other constructors are just convenience
 * methods.
 *
 * This class is a port of the S2Cell class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2Cell : S2Region, Comparable<S2Cell> {

    private var id: S2CellId

    /** The face of the cell. */
    private var face: Byte

    /** The level of the cell. */
    private var level: Byte

    /** The orientation of the cell. */
    private var orientation: Byte

    /** the bounds of this cell in (u,v)-space. */
    private var boundUV: R2Rect

    private constructor(id: S2CellId, face: Int, level: Int, orientation: Int, boundUV: R2Rect) {
        this.id = id
        this.face = face.toByte()
        this.level = level.toByte()
        this.orientation = orientation.toByte()
        this.boundUV = boundUV
    }

    constructor() {
        this.id = S2CellId.none
        this.face = 0
        this.level = 0
        this.orientation = 0
        this.boundUV = R2Rect.empty()
    }

    /**
     *
     */
    constructor(id: S2CellId) {
        this.id = id
        val faceIJ = id.toFaceIJOrientation(true)
        this.face = faceIJ.face.toByte()
        this.orientation = faceIJ.orientation!!.toByte()  // Compress int to a byte.
        this.level = id.level().toByte()
        this.boundUV = S2CellId.ijLevelToBoundUV(faceIJ.i, faceIJ.j, this.level.toInt())
    }

    // Convenience constructors.
    constructor(p: S2Point) : this(S2CellId.fromPoint(p))

    // The S2LatLng must be normalized.
    constructor(ll: S2LatLng) : this(S2CellId.fromLatLng(ll))

    fun id(): S2CellId = id

    fun face(): Int = face.toInt()

    fun level(): Int = level.toInt()

    fun orientation(): Int = orientation.toInt()

    fun boundUV(): R2Rect = boundUV

    fun isLeaf(): Boolean = level == S2CellId.kMaxLevel.toByte()

    // These are equivalent to the S2CellId methods, but have a more efficient
    // implementation since the level has been precomputed.
    fun getSizeIJ(): Int = TODO()
    fun getSizeST(): Double = TODO()

    // Returns the k-th vertex of the cell (k = 0,1,2,3).  Vertices are returned
    // in CCW order (lower left, lower right, upper right, upper left in the UV
    // plane).  The points returned by GetVertexRaw are not normalized.
    // For convenience, the argument is reduced modulo 4 to the range [0..3].
    fun getVertex(k: Int): S2Point = getVertexRaw(k).normalize()
    fun getVertexRaw(k: Int): S2Point {
        // Vertices are returned in the order SW, SE, NE, NW.
        val vertex = boundUV.getVertex(k)
        return S2Coords.faceUvToXyz(face(), vertex[0], vertex[1])
    }

    // Returns the inward-facing normal of the great circle passing through the
    // edge from vertex k to vertex k+1 (mod 4).  The normals returned by
    // GetEdgeRaw are not necessarily unit length.  For convenience, the
    // argument is reduced modulo 4 to the range [0..3].
    fun getEdge(k: Int): S2Point = getEdgeRaw(k).normalize()
    fun getEdgeRaw(k: Int): S2Point = when (k and 3) {
        0 -> S2Coords.getVNorm(face(), boundUV[1][0]) // South
        1 -> S2Coords.getUNorm(face(), boundUV[0][1]) // East
        2 -> -S2Coords.getVNorm(face(), boundUV[1][1]) // North
        else -> -S2Coords.getUNorm(face(), boundUV[0][0]) // West
    }

    // If this is not a leaf cell, sets children[0..3] to the four children of
    // this cell (in traversal order) and return true.  Otherwise returns false.
    // This method is equivalent to the following:
    //
    // for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
    //   children[pos] = S2Cell(id);
    //
    // except that it is more than two times faster.
    /**
     * Return the inward-facing normal of the great circle passing through the
     * edge from vertex k to vertex k+1 (mod 4). The normals returned by
     * GetEdgeRaw are not necessarily unit length.
     *
     *
     * If this is not a leaf cell, set children[0..3] to the four children of
     * this cell (in traversal order) and return true. Otherwise returns false.
     * This method is equivalent to the following:
     *
     *
     * for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
     * children[i] = S2Cell(id);
     *
     *
     * except that it is more than two times faster.
     */
    fun subdivide(children: Array<S2Cell>): Boolean {
        // This function is equivalent to just iterating over the child cell ids
        // and calling the S2Cell constructor, but it is about 2.5 times faster.
        if (id.isLeaf) return false

        // Compute the cell midpoint in uv-space.
        val uvMid = id.centerUV

        // Create four children with the appropriate bounds.
        var childId = id.childBegin()
        var pos = 0
        while (pos < 4) {
            val child = children[pos]
            child.face = face;
            child.level = (level + 1).toByte()
            child.orientation = (orientation.toInt() xor kPosToOrientation[pos]).toByte()
            child.id = childId;
            // We want to split the cell in half in "u" and "v".  To decide which
            // side to set equal to the midpoint value, we look at cell's (i,j)
            // position within its parent.  The index for "i" is in bit 1 of ij.
            val ij = kPosToIJ[orientation.toInt()][pos]
            val i = ij shr 1
            val j = ij and 1
            val bound = arrayOf(
                    doubleArrayOf(0.0, 0.0),
                    doubleArrayOf(0.0, 0.0),
            )
            bound[0][i] = boundUV[0][i];
            bound[0][1 - i] = uvMid[0];
            bound[1][j] = boundUV[1][j];
            bound[1][1 - j] = uvMid[1];
            child.boundUV = R2Rect(R1Interval(bound[0][0], bound[0][1]), R1Interval(bound[1][0], bound[1][1]))
            ++pos
            childId = childId.next()
        }
        return true
    }

    // Returns the direction vector corresponding to the center in (s,t)-space of
    // the given cell.  This is the point at which the cell is divided into four
    // subcells; it is not necessarily the centroid of the cell in (u,v)-space
    // or (x,y,z)-space.  The point returned by GetCenterRaw is not necessarily
    // unit length.
    fun getCenter(): S2Point = getCenterRaw().normalize()
    fun getCenterRaw(): S2Point = id.toPointRaw()

    // Returns the average area of cells at this level in steradians.  This is
    // accurate to within a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is
    // extremely cheap to compute.
    fun averageArea(): Double = averageArea(level.toInt())

    // Returns the approximate area of this cell in steradians.  This method is
    // accurate to within 3% percent for all cell sizes and accurate to within
    // 0.1% for cells at level 5 or higher (i.e. squares 350km to a side or
    // smaller on the Earth's surface).  It is moderately cheap to compute.
    fun approxArea(): Double {
        // All cells at the first two levels have the same area.
        if (level < 2) return averageArea(level())

        // First, compute the approximate area of the cell when projected
        // perpendicular to its normal.  The cross product of its diagonals gives
        // the normal, and the length of the normal is twice the projected area.
        val flatArea = 0.5 * (getVertex(2) - getVertex(0))
                .crossProd(getVertex(3) - getVertex(1))
                .norm()

        // Now, compensate for the curvature of the cell surface by pretending
        // that the cell is shaped like a spherical cap.  The ratio of the
        // area of a spherical cap to the area of its projected disc turns out
        // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
        // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
        // Here we set Pi*r*r == flat_area to find the equivalent disc.
        return flatArea * 2 / (1 + sqrt(1 - min(M_1_PI * flatArea, 1.0)));
    }

    // Returns the area of this cell as accurately as possible.  This method is
    // more expensive but it is accurate to 6 digits of precision even for leaf
    // cells (whose area is approximately 1e-18).
    fun exactArea(): Double {
        // There is a straightforward mathematical formula for the exact surface
        // area (based on 4 calls to asin), but as the cell size gets small this
        // formula has too much cancellation error.  So instead we compute the area
        // as the sum of two triangles (which is very accurate at all cell levels).
        val v0 = getVertex(0)
        val v1 = getVertex(1)
        val v2 = getVertex(2)
        val v3 = getVertex(3)
        return S2Measures.area(v0, v1, v2) + S2Measures.area(v0, v2, v3)
    }

    // Returns the distance from the cell to the given point.  Returns zero if
    // the point is inside the cell.
    fun getDistance(target: S2Point): S1ChordAngle = getDistanceInternal(target, true /*to_interior*/)

    // Return the distance from the cell boundary to the given point.
    fun getBoundaryDistance(target: S2Point): S1ChordAngle = getDistanceInternal(target, false /*to_interior*/);

    // Returns the maximum distance from the cell (including its interior) to the
    // given point.
    fun getMaxDistance(target: S2Point): S1ChordAngle {
        // First check the 4 cell vertices.  If all are within the hemisphere
        // centered around target, the max distance will be to one of these vertices.
        val targetUvw = S2Coords.faceXyzToUvw(face(), target)
        val maxDist = maxOf(
                maxOf(vertexChordDist(targetUvw, 0, 0), vertexChordDist(targetUvw, 1, 0)),
                maxOf(vertexChordDist(targetUvw, 0, 1), vertexChordDist(targetUvw, 1, 1))
        )

        if (maxDist <= S1ChordAngle.right()) {
            return maxDist;
        }

        // Otherwise, find the minimum distance d_min to the antipodal point and the
        // maximum distance will be Pi - d_min.
        return S1ChordAngle.straight() - getDistance(-target)
    }

    // Returns the minimum distance from the cell to the given edge AB.  Returns
    // zero if the edge intersects the cell interior.
    fun getDistance(a: S2Point, b: S2Point): S1ChordAngle {
        // Possible optimizations:
        //  - Currently the (cell vertex, edge endpoint) distances are computed
        //    twice each, and the length of AB is computed 4 times.
        //  - To fix this, refactor GetDistance(target) so that it skips calculating
        //    the distance to each cell vertex.  Instead, compute the cell vertices
        //    and distances in this function, and add a low-level UpdateMinDistance
        //    that allows the XA, XB, and AB distances to be passed in.
        //  - It might also be more efficient to do all calculations in UVW-space,
        //    since this would involve transforming 2 points rather than 4.

        // First, check the minimum distance to the edge endpoints A and B.
        // (This also detects whether either endpoint is inside the cell.)
        val minDist = minOf(getDistance(a), getDistance(b))
        if (minDist == S1ChordAngle.zero()) return minDist

        // Otherwise, check whether the edge crosses the cell boundary.
        // Note that S2EdgeCrosser needs pointers to vertices.
        val v = (0..3).map { i -> getVertex(i) }
        val crosser = S2EdgeCrosser(a, b, v[3])
        if (v.any { p -> crosser.crossingSign(p) >= 0 }) return S1ChordAngle.zero()

        // Finally, check whether the minimum distance occurs between a cell vertex
        // and the interior of the edge AB.  (Some of this work is redundant, since
        // it also checks the distance to the endpoints A and B again.)
        //
        // Note that we don't need to check the distance from the interior of AB to
        // the interior of a cell edge, because the only way that this distance can
        // be minimal is if the two edges cross (already checked above).
        for (i in 0..3) {
            S2EdgeDistances.updateMinDistance(v[i], a, b, minDist)
        }
        return minDist;
    }

    // Returns the maximum distance from the cell (including its interior) to the
    // given edge AB.
    fun getMaxDistance(a: S2Point, b: S2Point): S1ChordAngle {
        // If the maximum distance from both endpoints to the cell is less than Pi/2
        // then the maximum distance from the edge to the cell is the maximum of the
        // two endpoint distances.
        val maxDist = maxOf(getMaxDistance(a), getMaxDistance(b))
        if (maxDist <= S1ChordAngle.right()) {
            return maxDist
        }

        return S1ChordAngle.straight() - getDistance(-a, -b)
    }

    // Returns the distance from the cell to the given cell.  Returns zero if
    // one cell contains the other.
    fun getDistance(target: S2Cell): S1ChordAngle {
        // If the cells intersect, the distance is zero.  We use the (u,v) ranges
        // rather S2CellId::intersects() so that cells that share a partial edge or
        // corner are considered to intersect.
        if (face == target.face && boundUV.intersects(target.boundUV)) {
            return S1ChordAngle.zero()
        }

        // Otherwise, the minimum distance always occurs between a vertex of one
        // cell and an edge of the other cell (including the edge endpoints).  This
        // represents a total of 32 possible (vertex, edge) pairs.
        //
        // TODO(ericv): This could be optimized to be at least 5x faster by pruning
        // the set of possible closest vertex/edge pairs using the faces and (u,v)
        // ranges of both cells.
        val va = (0..3).map { i -> getVertex(i) }
        val vb = (0..3).map { i -> target.getVertex(i) }

        val minDist = S1ChordAngle.infinity()
        for (i in 0..3) {
            for (j in 0..3) {
                S2EdgeDistances.updateMinDistance(va[i], vb[j], vb[(j + 1) and 3], minDist)
                S2EdgeDistances.updateMinDistance(vb[i], va[j], va[(j + 1) and 3], minDist)
            }
        }
        return minDist
    }

    // Returns the maximum distance from the cell (including its interior) to the
    // given target cell.
    fun getMaxDistance(target: S2Cell): S1ChordAngle {
        // Need to check the antipodal target for intersection with the cell. If it
        // intersects, the distance is S1ChordAngle::Straight().
        if (face() == oppositeFace(target.face()) && boundUV.intersects(oppositeUV(target.boundUV))) {
            return S1ChordAngle.straight()
        }

        // Otherwise, the maximum distance always occurs between a vertex of one
        // cell and an edge of the other cell (including the edge endpoints).  This
        // represents a total of 32 possible (vertex, edge) pairs.
        //
        // TODO(user): When the maximum distance is at most Pi/2, the maximum is
        // always attained between a pair of vertices, and this could be made much
        // faster by testing each vertex pair once rather than the current 4 times.
        val va = (0..3).map { i -> getVertex(i) }
        val vb = (0..3).map { i -> target.getVertex(i) }
        val maxDist = S1ChordAngle.negative()
        for (i in 0..3) {
            for (j in 0..3) {
                S2EdgeDistances.updateMaxDistance(va[i], vb[j], vb[(j + 1) and 3], maxDist)
                S2EdgeDistances.updateMaxDistance(vb[i], va[j], va[(j + 1) and 3], maxDist)
            }
        }
        return maxDist;
    }

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public override fun clone(): S2Region {
        val clone = S2Cell()
        clone.face = face
        clone.level = level
        clone.orientation = orientation
        clone.boundUV = boundUV.clone()
        return clone
    }


    // Use the cell center in (u,v)-space as the cap axis. This vector is
    // very close to GetCenter() and faster to compute. Neither one of these
    // vectors yields the bounding cap with minimal surface area, but they
    // are both pretty close.
    //
    // It's possible to show that the two vertices that are furthest from
    // the (u,v)-origin never determine the maximum cap size (this is a
    // possible future optimization).
    override val capBound: S2Cap
        get() {
            // Use the cell center in (u,v)-space as the cap axis.  This vector is
            // very close to GetCenter() and faster to compute.  Neither one of these
            // vectors yields the bounding cap with minimal surface area, but they
            // are both pretty close.
            //
            // It's possible to show that the two vertices that are furthest from
            // the (u,v)-origin never determine the maximum cap size (this is a
            // possible future optimization).
            val center = S2Coords.faceUvToXyz(face(), boundUV.center).normalize()
            var cap = S2Cap.fromPoint(center)
            for (k in 0..3) {
                cap = cap.addPoint(getVertex(k))
            }
            return cap;
        }

    // attained at the vertices. Furthermore, the latitude range is
    // determined by one pair of diagonally opposite vertices and the
    // longitude range is determined by the other pair.
    //
    // We first determine which corner (i,j) of the cell has the largest
    // absolute latitude. To maximize latitude, we want to find the point in
    // the cell that has the largest absolute z-coordinate and the smallest
    // absolute x- and y-coordinates. To do this we look at each coordinate
    // (u and v), and determine whether we want to minimize or maximize that
    // coordinate based on the axis direction and the cell's (u,v) quadrant.
    // 35.26 degrees
    override val rectBound: S2LatLngRect
        get() {
            if (level > 0) {
                // Except for cells at level 0, the latitude and longitude extremes are
                // attained at the vertices.  Furthermore, the latitude range is
                // determined by one pair of diagonally opposite vertices and the
                // longitude range is determined by the other pair.
                //
                // We first determine which corner (i,j) of the cell has the largest
                // absolute latitude.  To maximize latitude, we want to find the point in
                // the cell that has the largest absolute z-coordinate and the smallest
                // absolute x- and y-coordinates.  To do this we look at each coordinate
                // (u and v), and determine whether we want to minimize or maximize that
                // coordinate based on the axis direction and the cell's (u,v) quadrant.
                val u = boundUV[0][0] + boundUV[0][1]
                val v = boundUV[1][0] + boundUV[1][1];
                val i = if (S2Coords.uAxis(face())[2] == 0.0) {
                    if (u < 0) 1 else 0
                } else {
                    if (u > 0) 1 else 0
                }
                val j = if (S2Coords.vAxis(face())[2] == 0.0) {
                    if (v < 0) 1 else 0
                } else {
                    if (v > 0) 1 else 0
                }
                val lat = R1Interval.fromPointPair(getLatitude(i, j), getLatitude(1 - i, 1 - j))
                val lng = S1Interval.fromPointPair(getLongitude(i, 1 - j), getLongitude(1 - i, j))

                // We grow the bounds slightly to make sure that the bounding rectangle
                // contains S2LatLng(P) for any point P inside the loop L defined by the
                // four *normalized* vertices.  Note that normalization of a vector can
                // change its direction by up to 0.5 * DBL_EPSILON radians, and it is not
                // enough just to add Normalize() calls to the code above because the
                // latitude/longitude ranges are not necessarily determined by diagonally
                // opposite vertex pairs after normalization.
                //
                // We would like to bound the amount by which the latitude/longitude of a
                // contained point P can exceed the bounds computed above.  In the case of
                // longitude, the normalization error can change the direction of rounding
                // leading to a maximum difference in longitude of 2 * DBL_EPSILON.  In
                // the case of latitude, the normalization error can shift the latitude by
                // up to 0.5 * DBL_EPSILON and the other sources of error can cause the
                // two latitudes to differ by up to another 1.5 * DBL_EPSILON, which also
                // leads to a maximum difference of 2 * DBL_EPSILON.
                return S2LatLngRect(lat, lng)
                        .expanded(S2LatLng.fromRadians(2 * DoubleType.epsilon, 2 * DoubleType.epsilon))
                        .polarClosure()
            }

            // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
            checkEQ((if (face() < 3) 1.0 else -1.0), S2Coords.norm(face())[face() % 3])

            val bound: S2LatLngRect = when (face()) {
                0 -> S2LatLngRect(R1Interval(-M_PI_4, M_PI_4), S1Interval(-M_PI_4, M_PI_4))
                1 -> S2LatLngRect(R1Interval(-M_PI_4, M_PI_4), S1Interval(M_PI_4, 3 * M_PI_4))
                2 -> S2LatLngRect(R1Interval(kPoleMinLat, M_PI_2), S1Interval.full())
                3 -> S2LatLngRect(R1Interval(-M_PI_4, M_PI_4), S1Interval(3 * M_PI_4, -3 * M_PI_4))
                4 -> S2LatLngRect(R1Interval(-M_PI_4, M_PI_4), S1Interval(-3 * M_PI_4, -M_PI_4))
                else -> S2LatLngRect(R1Interval(-M_PI_2, -kPoleMinLat), S1Interval.full())
            }
            // Finally, we expand the bound to account for the error when a point P is
            // converted to an S2LatLng to test for containment.  (The bound should be
            // large enough so that it contains the computed S2LatLng of any contained
            // point, not just the infinite-precision version.)  We don't need to expand
            // longitude because longitude is calculated via a single call to atan2(),
            // which is guaranteed to be semi-monotonic.  (In fact the Gnu implementation
            // is also correctly rounded, but we don't even need that here.)
            return bound.expanded(S2LatLng.fromRadians(DoubleType.epsilon, 0.0))
        }

    // The point 'p' does not need to be normalized.
    override fun contains(cell: S2Cell): Boolean {
        return id.contains(cell.id)
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        return id.intersects(cell.id)
    }

    // Returns true if the cell contains the given point "p".  Note that unlike
    // S2Loop/S2Polygon, S2Cells are considered to be closed sets.  This means
    // that points along an S2Cell edge (or at a vertex) belong to the adjacent
    // cell(s) as well.
    //
    // If instead you want every point to be contained by exactly one S2Cell,
    // you will need to convert the S2Cells to S2Loops (which implement point
    // containment this way).
    //
    // The point "p" does not need to be normalized.
    override operator fun contains(p: S2Point): Boolean {
        // We can't just call XYZtoFaceUV, because for points that lie on the
        // boundary between two faces (i.e. u or v is +1/-1) we need to return
        // true for both adjacent cells.
        val uv = S2Coords.faceXyztoUv(face(), p) ?: return false

        // Expand the (u,v) bound to ensure that
        //
        //   S2Cell(S2CellId(p)).Contains(p)
        //
        // is always true.  To do this, we need to account for the error when
        // converting from (u,v) coordinates to (s,t) coordinates.  At least in the
        // case of S2_QUADRATIC_PROJECTION, the total error is at most DBL_EPSILON.
        return boundUV.expanded(DoubleType.epsilon).contains(uv)
    }

    // Returns the latitude or longitude of the cell vertex given by (i,j),
    // where "i" and "j" are either 0 or 1.
    private fun getLatitude(i: Int, j: Int): Double {
        val p = S2Coords.faceUvToXyz(face(), boundUV[0][i], boundUV[1][j])
        return S2LatLng.latitude(p).radians
    }

    private fun getLongitude(i: Int, j: Int): Double {
        val p = S2Coords.faceUvToXyz(face(), boundUV[0][i], boundUV[1][j])
        return S2LatLng.longitude(p).radians
    }

    // Return the squared chord distance from point P to corner vertex (i,j).
    private fun vertexChordDist(p: S2Point, i: Int, j: Int): S1ChordAngle {
        val vertex = S2Point(boundUV[0][i], boundUV[1][j], 1.0).normalize()
        return S1ChordAngle.between(p, vertex)
    }

    // Given a point P and either the lower or upper edge of the S2Cell (specified
    // by setting "v_end" to 0 or 1 respectively), return true if P is closer to
    // the interior of that edge than it is to either endpoint.
    private fun uedgeIsClosest(target: S2Point, vEnd: Int): Boolean {
        val u0 = boundUV[0][0]
        val u1 = boundUV[0][1]
        val v = boundUV[1][vEnd]
        // These are the normals to the planes that are perpendicular to the edge
        // and pass through one of its two endpoints.
        val dir0 = S2Point(v * v + 1, -u0 * v, -u0)
        val dir1 = S2Point(v * v + 1, -u1 * v, -u1)
        return target.dotProd(dir0) > 0 && target.dotProd(dir1) < 0;
    }

    // Given a point P and either the left or right edge of the S2Cell (specified
    // by setting "u_end" to 0 or 1 respectively), return true if P is closer to
    // the interior of that edge than it is to either endpoint.
    private fun vedgeIsClosest(target: S2Point, uEnd: Int): Boolean {
        val v0 = boundUV[1][0]
        val v1 = boundUV[1][1]
        val u = boundUV[0][uEnd]
        // See comments above.
        val dir0 = S2Point(-u * v0, u * u + 1, -v0)
        val dir1 = S2Point(-u * v1, u * u + 1, -v1)
        return target.dotProd(dir0) > 0 && target.dotProd(dir1) < 0
    }

    // Returns the distance from the given point to the interior of the cell if
    // "to_interior" is true, and to the boundary of the cell otherwise.
    private fun getDistanceInternal(targetXYZ: S2Point, toInterior: Boolean): S1ChordAngle {
        // All calculations are done in the (u,v,w) coordinates of this cell's face.
        val target = S2Coords.faceXyzToUvw(face(), targetXYZ)

        // Compute dot products with all four upward or rightward-facing edge
        // normals.  "dirIJ" is the dot product for the edge corresponding to axis
        // I, endpoint J.  For example, dir01 is the right edge of the S2Cell
        // (corresponding to the upper endpoint of the u-axis).
        val dir00 = target[0] - target[2] * boundUV[0][0]
        val dir01 = target[0] - target[2] * boundUV[0][1]
        val dir10 = target[1] - target[2] * boundUV[1][0]
        val dir11 = target[1] - target[2] * boundUV[1][1]
        var inside = true
        if (dir00 < 0) {
            inside = false  // Target is to the left of the cell
            if (vedgeIsClosest(target, 0)) return edgeDistance(-dir00, boundUV[0][0])
        }
        if (dir01 > 0) {
            inside = false  // Target is to the right of the cell
            if (vedgeIsClosest(target, 1)) return edgeDistance(dir01, boundUV[0][1])
        }
        if (dir10 < 0) {
            inside = false  // Target is below the cell
            if (uedgeIsClosest(target, 0)) return edgeDistance(-dir10, boundUV[1][0])
        }
        if (dir11 > 0) {
            inside = false  // Target is above the cell
            if (uedgeIsClosest(target, 1)) return edgeDistance(dir11, boundUV[1][1])
        }
        if (inside) {
            if (toInterior) return S1ChordAngle.zero()
            // Although you might think of S2Cells as rectangles, they are actually
            // arbitrary quadrilaterals after they are projected onto the sphere.
            // Therefore the simplest approach is just to find the minimum distance to
            // any of the four edges.
            return minOf(
                    minOf(edgeDistance(-dir00, boundUV[0][0]), edgeDistance(dir01, boundUV[0][1])),
                    minOf(edgeDistance(-dir10, boundUV[1][0]), edgeDistance(dir11, boundUV[1][1]))
            )
        }
        // Otherwise, the closest point is one of the four cell vertices.  Note that
        // it is *not* trivial to narrow down the candidates based on the edge sign
        // tests above, because (1) the edges don't meet at right angles and (2)
        // there are points on the far side of the sphere that are both above *and*
        // below the cell, etc.
        return minOf(
                minOf(vertexChordDist(target, 0, 0), vertexChordDist(target, 1, 0)),
                minOf(vertexChordDist(target, 0, 1), vertexChordDist(target, 1, 1))
        )
    }

    // Return the latitude or longitude of the cell vertex given by (i,j),
    // where "i" and "j" are either 0 or 1.
    override fun toString(): String {
        return "[$face, $level, $orientation, $id]"
    }

    override fun compareTo(other: S2Cell): Int = id.compareTo(other.id)

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S2Cell) return false

        if (id != other.id) return false

        return true
    }

    override fun hashCode(): Int = id.hashCode()

    companion object {
        private const val kMaxCellSize = 1 shl S2CellId.kMaxLevel

        // We grow the bounds slightly to make sure that the bounding rectangle
        // also contains the normalized versions of the vertices. Note that the
        // maximum result magnitude is Pi, with a floating-point exponent of 1.
        // Therefore adding or subtracting 2**-51 will always change the result.
        private const val kMaxError = 1.0 / (1L shl 51)

        // The 4 cells around the equator extend to +/-45 degrees latitude at the
        // midpoints of their top and bottom edges. The two cells covering the
        // poles extend down to +/-35.26 degrees at their vertices.
        // adding kMaxError (as opposed to the C version) because of asin and atan2
        // roundoff errors
        private val kPoleMinLat = asin(sqrt(1.0 / 3.0)) - kMaxError

        // Returns the cell corresponding to the given S2 cube face.
        fun fromFace(face: Int): S2Cell {
            return S2Cell(S2CellId.fromFace(face))
        }

        // Returns a cell given its face (range 0..5), Hilbert curve position within
        // that face (an unsigned integer with S2CellId::kPosBits bits), and level
        // (range 0..kMaxLevel).  The given position will be modified to correspond
        // to the Hilbert curve position at the center of the returned cell.  This
        // is a static function rather than a constructor in order to indicate what
        // the arguments represent.
        @JvmStatic
        fun fromFacePosLevel(face: Int, pos: Byte, level: Int): S2Cell {
            return S2Cell(S2CellId.fromFacePosLevel(face, pos.toULong(), level))
        }

        /**
         * Return the average area of cells at this level. This is accurate to within
         * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
         * compute.
         */
        @JvmStatic
        fun averageArea(level: Int): Double {
            return S2Coords.projection.kAvgArea.getValue(level)
        }

        // Given the dot product of a point P with the normal of a u- or v-edge at the
        // given coordinate value, return the distance from P to that edge.
        private fun edgeDistance(dirIJ: Double, uv: Double): S1ChordAngle {
            // Let P by the target point and let R be the closest point on the given
            // edge AB.  The desired distance PR can be expressed as PR^2 = PQ^2 + QR^2
            // where Q is the point P projected onto the plane through the great circle
            // through AB.  We can compute the distance PQ^2 perpendicular to the plane
            // from "dirIJ" (the dot product of the target point P with the edge
            // normal) and the squared length the edge normal (1 + uv**2).
            val pq2 = (dirIJ * dirIJ) / (1 + uv * uv)

            // We can compute the distance QR as (1 - OQ) where O is the sphere origin,
            // and we can compute OQ^2 = 1 - PQ^2 using the Pythagorean theorem.
            // (This calculation loses accuracy as angle POQ approaches Pi/2.)
            val qr = 1 - sqrt(1 - pq2)
            return S1ChordAngle.fromLength2(pq2 + qr * qr)
        }

        private fun oppositeFace(face: Int): Int {
            return if (face >= 3) face - 3 else face + 3
        }

        // The antipodal UV is the transpose of the original UV, interpreted within
        // the opposite face.
        private fun oppositeUV(uv: R2Rect): R2Rect {
            return R2Rect(uv[1], uv[0])
        }

    }
}
