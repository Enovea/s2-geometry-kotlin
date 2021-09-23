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

import dilivia.PreConditions
import dilivia.math.M_PI
import dilivia.math.matrix.Matrix3x3Double
import dilivia.s2.S2PointUtil.fromFrame
import dilivia.s2.S2PointUtil.getFrame
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2LatLngRect
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.asin
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.pow
import org.apache.commons.math3.util.FastMath.sin
import org.apache.commons.math3.util.FastMath.sqrt
import kotlin.random.Random
import kotlin.random.nextUInt
import kotlin.random.nextULong

/**
 *
 */
object S2Random {

    var random = Random(1)

    // Reset the generator state using the given seed.
    @JvmStatic
    fun reset(seed: Int) {
        random = Random(seed)
    }

    // Return a uniformly distributed 64-bit unsigned integer.
    @JvmStatic
    fun randomULong(): ULong = random.nextULong()

    // Return a uniformly distributed 32-bit unsigned integer.
    @JvmStatic
    fun randomUInt(): UInt = random.nextUInt()

    @JvmStatic
    fun randomUInt(until: UInt) = random.nextUInt(until)

    @JvmStatic
    fun randomUShort(): UShort = random.nextUInt(0U, UShort.MAX_VALUE.toUInt()).toUShort()

    @JvmStatic
    fun randomDouble(): Double = random.nextDouble()

    // Return a uniformly distributed integer in the range [0,n).
    @JvmStatic
    fun randomInt(n: Int): Int = random.nextInt(n)

    @JvmStatic
    fun randomInt(min: Int, limit: Int): Int = random.nextInt(min, limit)

    @JvmStatic
    fun randomInt(): Int = random.nextInt()

    @JvmStatic
    fun randomLong(): Long = random.nextLong()

    @JvmStatic
    fun randomLong(min: Long, limit: Long): Long = random.nextLong(min, limit)

    @JvmStatic
    fun randomLong(n: Long): Long = random.nextLong(n)

    // Return a uniformly distributed "double" in the range [min, limit).
    @JvmStatic
    fun randomDouble(min: Double, limit: Double): Double = random.nextDouble(min, limit)

    // Return true with probability 1 in n.
    @JvmStatic
    fun oneIn(n: Int): Boolean = randomInt(n) == 0

    // Skewed: pick "base" uniformly from range [0,max_log] and then
    // return "base" random bits.  The effect is to pick a number in the
    // range [0,2^max_log-1] with bias towards smaller numbers.

    // Pick "base" uniformly from range [0,maxLog] and then return
    // "base" random bits. The effect is to pick a number in the range
    // [0,2^maxLog-1] with bias towards smaller numbers.
    @JvmStatic
    fun skewed(maxLog: Int): Int {
        val base = abs(random.nextInt()) % (maxLog + 1)
        // if (!base) return 0; // if 0==base, we & with 0 below.
        //
        // this distribution differs slightly from ACMRandom's Skewed,
        // since 0 occurs approximately 3 times more than 1 here, and
        // ACMRandom's Skewed never outputs 0.
        return random.nextInt() and (1 shl base) - 1
    }

    // Return a random unit-length vector.
    @JvmStatic
    fun randomPoint(): S2Point {
        // The order of evaluation of function arguments is unspecified,
        // so we may not just call S2Point with three RandDouble-based args.
        // Use temporaries to induce sequence points between calls.
        // The order of evaluation of function arguments is unspecified,
        // so we may not just call S2Point with three RandDouble-based args.
        // Use temporaries to induce sequence points between calls.
        val x: Double = randomDouble(-1.0, 1.0)
        val y: Double = randomDouble(-1.0, 1.0)
        val z: Double = randomDouble(-1.0, 1.0)
        val point = S2Point(x, y, z)
        point.normalize()
        return point
    }

    // Return a right-handed coordinate frame (three orthonormal vectors).
    @JvmStatic
    fun randomFrame(): Triple<S2Point, S2Point, S2Point> {
        return randomFrameAt(randomPoint())
    }

    // Given a unit-length z-axis, compute x- and y-axes such that (x,y,z) is a
    // right-handed coordinate frame (three orthonormal vectors).
    @JvmStatic
    fun randomFrameAt(z: S2Point): Triple<S2Point, S2Point, S2Point> {
        val x = z.crossProd(randomPoint()).normalized()
        val y = z.crossProd(x).normalized()
        return Triple(x, y, z)
    }

    // Return a cap with a random axis such that the log of its area is
    // uniformly distributed between the logs of the two given values.
    // (The log of the cap angle is also approximately uniformly distributed.)
    @JvmStatic
    fun randomCap(min_area: Double, max_area: Double): S2Cap {
        val capArea = max_area * pow(min_area / max_area, random.nextDouble())
        PreConditions.requireGE(capArea, min_area)
        PreConditions.requireLE(capArea, max_area)

        // The surface area of a cap is 2*Pi times its height.
        return S2Cap.fromCenterArea(randomPoint(), capArea)
    }

    @JvmStatic
    fun fromCols(cols: Triple<S2Point, S2Point, S2Point>): Matrix3x3Double =
            Matrix3x3Double.fromCols(cols.first, cols.second, cols.third)

    // Return a point chosen uniformly at random (with respect to area)
    // from the given cap.
    @JvmStatic
    fun samplePoint(cap: S2Cap): S2Point {
        // We consider the cap axis to be the "z" axis.  We choose two other axes to
        // complete the coordinate frame.
        val m = getFrame(cap.center)

        // The surface area of a spherical cap is directly proportional to its
        // height.  First we choose a random height, and then we choose a random
        // point along the circle at that height.
        val h = random.nextDouble() * cap.height
        val theta = 2 * M_PI * random.nextDouble()
        val r = sqrt(h * (2 - h))  // Radius of circle.

        // The result should already be very close to unit-length, but we might as
        // well make it accurate as possible.
        return fromFrame(m, S2Point(cos(theta) * r, sin(theta) * r, 1 - h)).normalize()
    }

    // Return a point chosen uniformly at random (with respect to area on the
    // sphere) from the given latitude-longitude rectangle.
    @JvmStatic
    fun samplePoint(rect: S2LatLngRect): S2Point {
        // First choose a latitude uniformly with respect to area on the sphere.
        val sinLo = sin(rect.lat.lo)
        val sinHi = sin(rect.lat.hi)
        val lat = asin(random.nextDouble(sinLo, sinHi))

        // Now choose longitude uniformly within the given range.
        val lng = rect.lng.lo + random.nextDouble() * rect.lng.length
        return S2LatLng.fromRadians(lat, lng).normalized().toPoint()
    }

    // Return a random cell id at the given level or at a randomly chosen
    // level.  The distribution is uniform over the space of cell ids,
    // but only approximately uniform over the surface of the sphere.
    @JvmStatic
    fun randomCellId(level: Int): S2CellId {
        val face = random.nextInt(S2CellId.kNumFaces)
        val pos = random.nextLong().toULong() and ((1UL shl S2CellId.kPosBits) - 1UL)
        return S2CellId.fromFacePosLevel(face, pos, level)
    }
    @JvmStatic
    fun randomCellId(): S2CellId {
        return randomCellId(random.nextInt(S2CellId.kMaxLevel + 1))
    }

}
