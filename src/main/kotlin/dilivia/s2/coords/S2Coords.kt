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

import dilivia.PreConditions.checkState
import dilivia.PreConditions.requireArgument
import dilivia.math.vectors.R2Point
import dilivia.s2.S2Point
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.round

enum class Projections {
    LINEAR, TAN, QUADRATIC
}

/**
 * S2 coordinate systems.
 *
 * This object contains the various coordinate systems used throughout the library. Most importantly, S2 defines a
 * framework for decomposing the unit sphere into a hierarchy of "cells".  Each cell is a quadrilateral bounded by four
 * geodesics. The top level of the hierarchy is obtained by projecting the six faces of a cube onto the unit sphere,
 * and lower levels are obtained by subdividing each cell into four children recursively. Cells are numbered such that
 * sequentially increasing cells follow a continuous space-filling curve over the entire sphere. The transformation is
 * designed to make the cells at each level fairly uniform in size.
 *
 *
 * S2Cell Decomposition
 * ---------------------
 *
 * The following methods define the cube-to-sphere projection used by the S2Cell decomposition.
 *
 * In the process of converting a latitude-longitude pair to a 64-bit cell id, the following coordinate systems
 * are used:
 *
 *  (id)
 *    An S2CellId is a 64-bit encoding of a face and a Hilbert curve position on that face.
 *    The Hilbert curve position implicitly encodes both the position of a cell and its subdivision level
 *    (see S2CellId).
 *
 *  (face, i, j)
 *    Leaf-cell coordinates.  "i" and "j" are integers in the range [0,(2**30)-1] that identify a particular
 *    leaf cell on the given face.
 *    The (i, j) coordinate system is right-handed on each face, and the faces are oriented such that Hilbert
 *    curves connect continuously from one face to the next.
 *
 *  (face, s, t)
 *    Cell-space coordinates.  "s" and "t" are real numbers in the range [0,1] that identify a point on the
 *    given face.  For example, the point (s, t) = (0.5, 0.5) corresponds to the center of the top-level face
 *    cell.  This point is also a vertex of exactly four cells at each subdivision level greater than zero.
 *
 *  (face, si, ti)
 *    Discrete cell-space coordinates.  These are obtained by multiplying "s" and "t" by 2**31 and rounding
 *    to the nearest unsigned integer. Discrete coordinates lie in the range [0,2**31].  This coordinate system
 *    can represent the edge and center positions of all cells with no loss of precision (including non-leaf cells).
 *    In binary, each coordinate of a level-k cell center ends with a 1 followed by (30 - k) 0s.
 *    The coordinates of its edges end with (at least) (31 - k) 0s.
 *
 *  (face, u, v)
 *    Cube-space coordinates in the range [-1,1].  To make the cells at each level more uniform in size after
 *    they are projected onto the sphere, we apply a nonlinear transformation of the form u=f(s), v=f(t).
 *    The (u, v) coordinates after this transformation give the actual coordinates on the cube face
 *    (modulo some 90 degree rotations) before it is projected onto the unit sphere.
 *
 *  (face, u, v, w)
 *    Per-face coordinate frame.  This is an extension of the (face, u, v) cube-space coordinates that adds a
 *    third axis "w" in the direction of the face normal.  It is always a right-handed 3D coordinate system.
 *    Cube-space coordinates can be converted to this frame by setting w=1, while (u,v,w) coordinates can be
 *    projected onto the cube face by dividing by w, i.e. (face, u/w, v/w).
 *
 *  (x, y, z)
 *    Direction vector (S2Point).  Direction vectors are not necessarily unit length, and are often chosen
 *    to be points on the biunit cube [-1,+1]x[-1,+1]x[-1,+1].  They can be be normalized to obtain the
 *    corresponding point on the unit sphere.
 *
 *  (lat, lng)
 *    Latitude and longitude (S2LatLng).  Latitudes must be between -90 and 90 degrees inclusive, and
 *    longitudes must be between -180 and 180 degrees inclusive.
 *
 * Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are right-handed on all six faces.
 *
 * This class is a port of s2coords and s2metrics of the Google S2 Geometry project
 * (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2Coords {

    // The canonical Hilbert traversal order looks like an inverted 'U':
    // the subcells are visited in the order (0,0), (0,1), (1,1), (1,0).
    // The following tables encode the traversal order for various
    // orientations of the Hilbert curve (axes swapped and/or directions
    // of the axes reversed).

    // Together these flags define a cell orientation.  If 'kSwapMask'
    // is true, then canonical traversal order is flipped around the
    // diagonal (i.e. i and j are swapped with each other).  If
    // 'kInvertMask' is true, then the traversal order is rotated by 180
    // degrees (i.e. the bits of i and j are inverted, or equivalently,
    // the axis directions are reversed).
    internal val kSwapMask = 0x01
    internal val kInvertMask = 0x02

    // kIJtoPos[orientation][ij] -> pos
    //
    // Given a cell orientation and the (i,j)-index of a subcell (0=(0,0),
    // 1=(0,1), 2=(1,0), 3=(1,1)), return the order in which this subcell is
    // visited by the Hilbert curve (a position in the range [0..3]).
    internal val kIJtoPos = arrayOf(
        //        (0,0) (0,1) (1,0) (1,1)
        intArrayOf(0, 1, 3, 2),  // canonical order
        intArrayOf(0, 3, 1, 2),  // axes swapped
        intArrayOf(2, 3, 1, 0),  // bits inverted
        intArrayOf(2, 1, 3, 0),  // swapped & inverted
    ) //[4][4];

    // kPosToIJ[orientation][pos] -> ij
    //
    // Return the (i,j) index of the subcell at the given position 'pos' in the
    // Hilbert curve traversal order with the given orientation.  This is the
    // inverse of the previous table:
    //
    //   kPosToIJ[r][kIJtoPos[r][ij]] == ij
    internal val kPosToIJ = arrayOf(
        //         0  1  2  3
        intArrayOf(0, 1, 3, 2),    // canonical order:    (0,0), (0,1), (1,1), (1,0)
        intArrayOf(0, 2, 3, 1),    // axes swapped:       (0,0), (1,0), (1,1), (0,1)
        intArrayOf(3, 2, 0, 1),    // bits inverted:      (1,1), (1,0), (0,0), (0,1)
        intArrayOf(3, 1, 0, 2),    // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    ) //[4][4];

    // kPosToOrientation[pos] -> orientation_modifier
    //
    // Return a modifier indicating how the orientation of the child subcell
    // with the given traversal position [0..3] is related to the orientation
    // of the parent cell.  The modifier should be XOR-ed with the parent
    // orientation to obtain the curve orientation in the child.
    val kPosToOrientation = intArrayOf(
        kSwapMask,
        0,
        0,
        kInvertMask + kSwapMask,
    )

    // The U,V,W axes for each face.
    val kFaceUVWAxes //[6][3][3]
            = arrayOf(
        arrayOf(
            intArrayOf(0, 1, 0),
            intArrayOf(0, 0, 1),
            intArrayOf(1, 0, 0)
        ),
        arrayOf(
            intArrayOf(-1, 0, 0),
            intArrayOf(0, 0, 1),
            intArrayOf(0, 1, 0)
        ),
        arrayOf(
            intArrayOf(-1, 0, 0),
            intArrayOf(0, -1, 0),
            intArrayOf(0, 0, 1)
        ),
        arrayOf(
            intArrayOf(0, 0, -1),
            intArrayOf(0, -1, 0),
            intArrayOf(-1, 0, 0)
        ),
        arrayOf(
            intArrayOf(0, 0, -1),
            intArrayOf(1, 0, 0),
            intArrayOf(0, -1, 0)
        ),
        arrayOf(
            intArrayOf(0, 1, 0),
            intArrayOf(1, 0, 0),
            intArrayOf(0, 0, -1)
        )
    )

    // The precomputed neighbors of each face (see GetUVWFace).
    val kFaceUVWFaces //[6][3][2]
            = arrayOf(
        arrayOf(intArrayOf(4, 1), intArrayOf(5, 2), intArrayOf(3, 0)),
        arrayOf(intArrayOf(0, 3), intArrayOf(5, 2), intArrayOf(4, 1)),
        arrayOf(intArrayOf(0, 3), intArrayOf(1, 4), intArrayOf(5, 2)),
        arrayOf(intArrayOf(2, 5), intArrayOf(1, 4), intArrayOf(0, 3)),
        arrayOf(intArrayOf(2, 5), intArrayOf(3, 0), intArrayOf(1, 4)),
        arrayOf(intArrayOf(4, 1), intArrayOf(3, 0), intArrayOf(2, 5))
    )

    /** The sphere to unit cube projection. */
    var projection: S2Projection = S2QuadraticProjection

    /** This is the number of levels needed to specify a leaf cell. */
    const val kMaxCellLevel = 30

    /** The maximum index of a valid leaf cell plus one. The range of valid leaf cell indices is [0..kLimitIJ-1]. */
    val kLimitIJ = 1 shl kMaxCellLevel  // == S2CellId::kMaxSize

    /** The maximum value of an si- or ti-coordinate.  The range of valid (si,ti) values is [0..kMaxSiTi]. */
    val kMaxSiTi = 1.toUInt() shl (kMaxCellLevel + 1)

    /**
     * Convert an s- or t-value to the corresponding u- or v-value.  This is a non-linear transformation
     * from [-1,1] to [-1,1] that attempts to make the cell sizes more uniform.
     *
     * @param s a s- or t-value in [-1,1]
     * @return The corresponding u- or v-value (in [-1,1]).
     */
    fun stToUv(s: Double): Double = projection.stToUv(s)

    /**
     * The inverse of the stToUv transformation.  Note that it is not always true that uvToSt(stToUv(x)) == x
     * due to numerical errors.
     *
     * @param u u- or v-value (in [-1,1])
     * @return The corresponding s- or t-value (in [-1,1]).
     */
    fun uvToSt(u: Double): Double = projection.uvToSt(u)

    /**
     * Convert the i- or j-index of a leaf cell to the minimum corresponding s-or t-value contained by that cell.
     * The argument must be in the range [0..2**30], i.e. up to one position beyond the normal range of valid leaf
     * cell indices.
     *
     * @param i a i- or j-index of a leaf cell
     * @return The corresponding minimun s- or t-value contained by the cell.
     */
    fun ijToStMin(i: Int): Double {
        requireArgument { i in 0..kLimitIJ }
        return (1.0 / kLimitIJ) * i
    }

    /**
     * Get the i- or j-index of the leaf cell containing the given s- or t-value.
     * If the argument is outside the range spanned by valid leaf cell indices, return the index of the closest valid
     * leaf cell (i.e., return values are clamped to the range of valid leaf cell indices).
     *
     * @param s a s- or t-value
     * @return The corresponding cell index.
     */
    fun stToIj(s: Double): Int = max(0, min(kLimitIJ - 1, round(kLimitIJ * s - 0.5).toInt()))

    /**
     * Convert an si- or ti-value to the corresponding s- or t-value.
     *
     * @param si a si- or ti-value
     * @return The corresponding s- or t- value.
     */
    fun siTiToSt(si: UInt): Double {
        requireArgument { si <= kMaxSiTi }
        return (1.0 / kMaxSiTi.toDouble()) * si.toDouble()
    }

    /**
     * Get the si- or ti-coordinate that is nearest to the given s- or t-value.
     * The result may be outside the range of valid (si,ti)-values.
     *
     * @param s a s- or t-value
     * @return The corresponding si- or ti- value.
     */
    fun stToSiTi(s: Double): UInt {
        // kMaxSiTi == 2^31, so the result doesn't fit in an int32 when s == 1.
        return round(s * kMaxSiTi.toDouble()).toUInt()
    }

    /**
     * Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
     *
     * @param faceUV a (face, u, v) coordinates
     * @return The corresponding sphere point
     */
    fun faceUvToXyz(faceUV: FaceUV): S2Point = faceUvToXyz(faceUV.face, faceUV.u, faceUV.v)

    /**
     * Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
     *
     * @param face The face number
     * @param uv The (u,v) coordinates
     * @return The corresponding sphere point
     */
    fun faceUvToXyz(face: Int, uv: R2Point): S2Point = faceUvToXyz(face, uv[0], uv[1])

    /**
     * Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
     *
     * @param face The face number
     * @param u the u coordinate
     * @param v the v coordinate
     * @return The corresponding sphere point
     */
    fun faceUvToXyz(face: Int, u: Double, v: Double): S2Point {
        return when (face) {
            0 -> S2Point(1.0, u, v)
            1 -> S2Point(-u, 1.0, v)
            2 -> S2Point(-u, -v, 1.0)
            3 -> S2Point(-1.0, -v, -u)
            4 -> S2Point(v, -1.0, -u)
            else -> S2Point(v, u, -1.0)
        }
    }

    /**
     * If the dot product of p with the given face normal is positive, return the corresponding u and v values
     * (which may lie outside the range [-1,1]).  Otherwise return null.
     *
     * @param face The face number.
     * @param p The direction vector to convert.
     * @return The (u,v) coordinate of the vector if the dot product of p with the given face normal is positive else
     * null.
     */
    fun faceXyztoUv(face: Int, p: S2Point): R2Point? {
        if (face < 3) {
            if (p[face] <= 0) return null
        } else {
            if (p[face - 3] >= 0) return null
        }
        return validFaceXyzToUv(face, p)
    }

    /**
     * Given a *valid* face for the given point p (meaning that dot product of p with the face normal is positive),
     * set the corresponding u and v values (which may lie outside the range [-1,1]).
     *
     * @param face The face number
     * @param p A point
     * @param uv The (u,v) coordinate point to update.
     */
    fun validFaceXyzToUv(face: Int, p: S2Point, uv: R2Point) {
        requireArgument { p.dotProd(norm(face)) > 0 }
        when (face) {
            0 -> {
                uv[0] = p[1] / p[0]; uv[1] = p[2] / p[0]; }
            1 -> {
                uv[0] = -p[0] / p[1]; uv[1] = p[2] / p[1]; }
            2 -> {
                uv[0] = -p[0] / p[2]; uv[1] = -p[1] / p[2]; }
            3 -> {
                uv[0] = p[2] / p[0]; uv[1] = p[1] / p[0]; }
            4 -> {
                uv[0] = p[2] / p[1]; uv[1] = -p[0] / p[1]; }
            else -> {
                uv[0] = -p[1] / p[2]; uv[1] = -p[0] / p[2]; }
        }
    }

    /**
     * Given a *valid* face for the given point p (meaning that dot product of p with the face normal is positive),
     * return the corresponding u and v values (which may lie outside the range [-1,1]).
     *
     * @param face The face number
     * @param p A point
     * @return The (u,v) coordinate point to update.
     */
    fun validFaceXyzToUv(face: Int, p: S2Point): R2Point {
        val uv = R2Point()
        validFaceXyzToUv(face, p, uv)
        return uv
    }

    /**
     * Transform the given point P to the (u,v,w) coordinate frame of the given face
     * (where the w-axis represents the face normal).
     *
     * @param face The face number
     * @param p The point to transform.
     * @return the (u,v,w) coordinate frame
     */
    fun faceXyzToUvw(face: Int, p: S2Point): S2Point {
        // The result coordinates are simply the dot products of P with the (u,v,w)
        // axes for the given face (see kFaceUVWAxes).
        return when (face) {
            0 -> S2Point(p.y, p.z, p.x)
            1 -> S2Point(-p.x, p.z, p.y)
            2 -> S2Point(-p.x, -p.y, p.z)
            3 -> S2Point(-p.z, -p.y, -p.x)
            4 -> S2Point(-p.z, p.x, -p.y)
            else -> S2Point(p.y, p.x, -p.z)
        }
    }

    /**
     * Get the face containing the given direction vector.
     * (For points on the boundary between faces, the result is arbitrary but repeatable.)
     *
     * @param p A point.
     * @return The face containing the given direction vector.
     */
    fun face(p: S2Point): Int {
        var face = p.largestAbsComponent()
        if (p[face] < 0) face += 3
        return face
    }

    /**
     * Convert a direction vector (not necessarily unit length) to (face, u, v) coordinates.
     *
     * @param p A point
     * @return The corresponding (face, u, v) coordinates.
     */
    fun xyzToFaceUv(p: S2Point): FaceUV {
        val face = face(p)
        val uv = validFaceXyzToUv(face, p)
        return FaceUV(face = face, u = uv.x, v = uv.y)
    }

    /**
     * Convert a direction vector (not necessarily unit length) to (face, si, ti) coordinates and, if p is exactly
     * equal to the center of a cell, return the level of this cell (-1 otherwise).
     *
     * @param p A point
     * @return cell level and (face, si, ti) coordinates pair.
     */
    fun xyzToFaceSiTi(p: S2Point): Pair<Int, FaceSiTi> {
        val (face, u, v) = xyzToFaceUv(p)
        val si = stToSiTi(projection.uvToSt(u))
        val ti = stToSiTi(projection.uvToSt(v))
        val faceSiTi = FaceSiTi(face, si, ti)
        // If the levels corresponding to si,ti are not equal, then p is not a cell
        // center.  The si,ti values 0 and kMaxSiTi need to be handled specially
        // because they do not correspond to cell centers at any valid level; they
        // are mapped to level -1 by the code below.
        val level = kMaxCellLevel - (si or kMaxSiTi).countTrailingZeroBits()
        if (level < 0 || level != kMaxCellLevel - (ti or kMaxSiTi).countTrailingZeroBits()) {
            return -1 to faceSiTi
        }
        checkState { level <= kMaxCellLevel }
        // In infinite precision, this test could be changed to ST == SiTi. However,
        // due to rounding errors, UVtoST(XYZtoFaceUV(FaceUVtoXYZ(STtoUV(...)))) is
        // not idempotent. On the other hand, center_raw is computed exactly the same
        // way p was originally computed (if it is indeed the center of an S2Cell):
        // the comparison can be exact.
        val center = faceSiTiToXyz(faceSiTi).normalize()
        return (if (p == center) level else -1) to faceSiTi
    }

    /**
     * Convert (face, si, ti) coordinates to a direction vector (not necessarily unit length).
     *
     * @param faceSiTi (face, si, ti) coordinates
     * @return The corresponding direction vector.
     */
    fun faceSiTiToXyz(faceSiTi: FaceSiTi): S2Point = faceSiTiToXyz(faceSiTi.face, faceSiTi.si, faceSiTi.ti)

    /**
     * Convert (face, si, ti) coordinates to a direction vector (not necessarily unit length).
     *
     * @param face The face number
     * @param si The si coordinate.
     * @param ti The ti coordinate.
     * @return The corresponding direction vector.
     */
    fun faceSiTiToXyz(face: Int, si: UInt, ti: UInt): S2Point {
        val u = projection.stToUv(siTiToSt(si))
        val v = projection.stToUv(siTiToSt(ti))
        return faceUvToXyz(face, u, v)
    }

    /**
     * Get the right-handed normal (not necessarily unit length) for an edge in the direction of the positive v-axis
     * at the given u-value on the given face.  (This vector is perpendicular to the plane through the sphere origin
     * that contains the given edge.)
     *
     * @param face The face number.
     * @param u u coordinate.
     * @return the right-handed normal of the edge.
     */
    fun getUNorm(face: Int, u: Double): S2Point = when (face) {
        0 -> S2Point(u, -1.0, 0.0)
        1 -> S2Point(1.0, u, 0.0)
        2 -> S2Point(1.0, 0.0, u)
        3 -> S2Point(-u, 0.0, 1.0)
        4 -> S2Point(0.0, -u, 1.0)
        else -> S2Point(0.0, -1.0, -u)
    }

    // Return the right-handed normal (not necessarily unit length) for an
    // edge in the direction of the positive u-axis at the given v-value on
    // the given face.

    /**
     * Get the right-handed normal (not necessarily unit length) for an edge in the direction of the positive u-axis
     * at the given v-value on the given face.
     *
     * @param face The face number.
     * @param v v coordinate.
     * @return the right-handed normal of the edge.
     */
    fun getVNorm(face: Int, v: Double): S2Point {
        return when (face) {
            0 -> S2Point(-v, 0.0, 1.0)
            1 -> S2Point(0.0, -v, 1.0)
            2 -> S2Point(0.0, -1.0, -v)
            3 -> S2Point(v, -1.0, 0.0)
            4 -> S2Point(1.0, v, 0.0)
            else -> S2Point(1.0, 0.0, v)
        }
    }

    // Return the unit-length normal, u-axis, or v-axis for the given face.
    fun norm(face: Int): S2Point = uvwAxis(face, 2)

    fun uAxis(face: Int): S2Point = uvwAxis(face, 0)

    fun vAxis(face: Int): S2Point = uvwAxis(face, 1)

    /**
     * Get the given axis of the given face (u=0, v=1, w=2).
     *
     * @param face The face number.
     * @param axis The axis (u=0, v=1, w=2)
     * @return The uvw axis.
     */
    fun uvwAxis(face: Int, axis: Int): S2Point {
        val p = kFaceUVWAxes[face][axis]
        return S2Point(p[0], p[1], p[2])
    }

    /**
     * With respect to the (u,v,w) coordinate system of a given face, return the face that lies in the given
     * direction (negative=0, positive=1) of the given axis (u=0, v=1, w=2).  For example, getUVWFace(4, 0, 1)
     * returns the face that is adjacent to face 4 in the positive u-axis direction.
     *
     * @param face The face number.
     * @param axis The axis (u=0, v=1, w=2)
     * @param direction The direction (negative=0, positive=1)
     * @return the face.
     */
    fun uvwFace(face: Int, axis: Int, direction: Int): Int {
        requireArgument { face in 0..5 }
        requireArgument { axis in 0..2 }
        requireArgument { direction in 0..1 }
        return kFaceUVWFaces[face][axis][direction]
    }

}

data class FaceIJ(val face: Int, val i: Int, val j: Int, val orientation: Int?)
data class FaceSiTi(val face: Int, val si: UInt, val ti: UInt)
data class FaceUV(val face: Int, val u: Double, val v: Double) {
    fun getUV(): R2Point = R2Point(u, v)
}

