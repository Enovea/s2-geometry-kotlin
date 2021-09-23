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
import dilivia.math.M_PI_2
import dilivia.math.vectors.R2Point
import dilivia.math.vectors.R2VectorDouble
import dilivia.s2.S1Angle.Companion.degrees
import dilivia.s2.S1Angle.Companion.e5
import dilivia.s2.S1Angle.Companion.e6
import dilivia.s2.S1Angle.Companion.e7
import dilivia.s2.S1Angle.Companion.radians
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.FastMath.IEEEremainder
import org.apache.commons.math3.util.FastMath.abs
import org.apache.commons.math3.util.FastMath.atan2
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.max
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sin
import org.apache.commons.math3.util.FastMath.sqrt

/**
 * This class represents a point on the unit sphere as a pair of latitude-longitude coordinates. Like the rest of the
 * "geometry" package, the intent is to represent spherical geometry as a mathematical abstraction, so functions that
 * are specifically related to the Earth's geometry (e.g. easting/northing conversions) should be put elsewhere.
 *
 * This class is a port of the S2LatLng class of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class S2LatLng : Comparable<S2LatLng>,Cloneable {

    /**  Latitude of the point in radians. */
    var latRadians: Double
    /** Longitude of the point in radians. */
    var lngRadians: Double

    /**
     * This is internal to avoid ambiguity about which units are expected.
     */
    private constructor(latRadians: Double, lngRadians: Double) {
        this.latRadians = latRadians
        this.lngRadians = lngRadians
    }

    constructor(): this(0,0)

    constructor(latRadians: Int, lngRadians: Int): this(latRadians.toDouble(), lngRadians.toDouble())

    private constructor(p: R2Point): this(p[0], p[1])

    public override fun clone(): S2LatLng {
        return S2LatLng(latRadians, lngRadians)
    }
    /**
     * Get the latitude of this point as a new S1Angle.
     * @return A S1Angle instance that represents the latitude of the point.
     */
    fun lat(): S1Angle {
        return radians(latRadians)
    }

    /**
     * Get the latitude of this point as degrees.
     * @return the latitude in degrees.
     */
    @Strictfp
    fun latDegrees(): Double {
        return 180.0 / M_PI * latRadians
    }

    @Strictfp
    fun latDegrees(degrees: Double) {
        latRadians = degrees * M_PI / 180.0
    }

    /**
     * Get the longitude of this point as a new S1Angle.
     * @return A S1Angle instance that represents the longitude of the point.
     */
    fun lng(): S1Angle {
        return radians(lngRadians)
    }

    /**
     * Get the longitude of this point as degrees.
     * @return The longitude in degrees.
     */
    fun lngDegrees(): Double {
        return 180.0 / M_PI * lngRadians
    }

    fun lngDegrees(degrees: Double) {
        lngRadians = degrees * M_PI / 180.0
    }

    var coords: R2VectorDouble
        get() = R2VectorDouble(x = latRadians, y = lngRadians)
        set(value) {
            latRadians = value.x
            lngRadians = value.y
        }

    /**
     * Indicates if the point is valid. Is true if the latitude is between -90 and 90 degrees inclusive and the
     * longitude is between -180 and 180 degrees inclusive.
     */
    val isValid: Boolean
        get() = abs(lat().radians) <= M_PI_2 && abs(lng().radians) <= M_PI

    /**
     * Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts a multiple of 360 degrees to the
     * longitude if necessary to reduce it to the range [-180, 180].
     * If the current point is valid then the returned point will have the same coordinates.
     *
     * @return a new S2LatLng based on this instance for which isValid will be `true`.
     */
    fun normalized(): S2LatLng {
        // [-S2.M_PI, S2.M_PI] inclusive, which is what we want here.
        // drem(x, 2 * S2.M_PI) reduces its argument to the range
        return S2LatLng(
                latRadians = max(-M_PI_2, min(M_PI_2, lat().radians)),
                lngRadians = IEEEremainder(lng().radians,2 * M_PI)
        )
    }

    fun normalize(): S2LatLng {
        latRadians = max(-M_PI_2, min(M_PI_2, lat().radians))
        lngRadians = IEEEremainder(lng().radians, 2 * M_PI)
        return this
    }

    /**
     * Converts a normalized S2LatLng to the equivalent unit-length vector. The maximum error in the result is
     * 1.5 * DBL_EPSILON.  (This does not include the error of converting degrees, E5, E6, or E7 to radians.)
     *
     * @return The unit-length vector that represents this point.
     */
    fun toPoint(): S2Point {
        val phi = lat().radians
        val theta = lng().radians
        val cosphi = cos(phi)
        return dilivia.s2.S2Point(cos(theta) * cosphi, sin(theta) * cosphi, sin(phi))
    }

    /**
     * Compute the distance between this point and another S2LatLng instance (measured along the surface of the sphere).
     * This implements the Haversine formula, which is numerically stable for small distances but only gets about 8
     * digits of precision for very large distances (e.g. antipodal points). Note that 8 digits is still accurate to
     * within about 10cm for a sphere the size of the Earth.
     *
     * This is equivalent to
     *
     *      S1Angle(toPoint(), o.toPoint())
     *
     * except that this function is slightly faster, and is also somewhat less accurate for distances approaching 180
     * degrees.  Both S2LatLngs must be normalized.
     *
     * sinlat = sin((lat2 - lat1) / 2)
     * sinlng = sin((lng2 - lng1) / 2)
     * x = sinlat^2 + sinlng^2 * cos(lat1) * cos(lat2)
     * d = 2 * atan(sqrt(x) / sqrt(1 - x))
     *
     * @param o A S2LatLng point.
     * @return The distance between this point and o measured along the surface of the unit sphere.
     */
    fun getDistance(o: S2LatLng): S1Angle {
        val lat1 = lat().radians
        val lat2 = o.lat().radians
        val lng1 = lng().radians
        val lng2 = o.lng().radians
        val dlat = sin(0.5 * (lat2 - lat1))
        val dlng = sin(0.5 * (lng2 - lng1))
        val x = dlat * dlat + dlng * dlng * cos(lat1) * cos(lat2)
        return radians(2 * atan2(sqrt(x), sqrt(max(0.0, 1.0 - x))))
    }

    /**
     * Sum this point with an other. The result is not normalized.
     *
     * @param other Another point.
     * @return The sum of the 2 points.
     */
    operator fun plus(other: S2LatLng): S2LatLng = fromRadians(
            latRadians = latRadians + other.latRadians,
            lngRadians = lngRadians + other.lngRadians
    )

    /**
     * Subtract this point with an other. The result is not normalized.
     *
     * @param other Another point.
     * @return The subtraction of the 2 points.
     */
    operator fun minus(other: S2LatLng): S2LatLng = fromRadians(
            latRadians = latRadians - other.latRadians,
            lngRadians = lngRadians - other.lngRadians
    )

    /**
     * Multiply the lat-lng pair by a constant.
     * @param other A constant.
     * @return A new S2LatLng instance that represents this point multiply by other.
     */
    operator fun times(other: Double): S2LatLng = fromRadians(
            latRadians = other * latRadians,
            lngRadians = other * lngRadians
    )

    override fun equals(other: Any?): Boolean {
        if (other is S2LatLng) {
            val o = other
            return latRadians == o.latRadians && lngRadians == o.lngRadians
        }
        return false
    }

    override fun hashCode(): Int {
        var value: Long = 17
        value += 37 * value + java.lang.Double.doubleToLongBits(latRadians)
        value += 37 * value + java.lang.Double.doubleToLongBits(lngRadians)
        return (value xor (value ushr 32)).toInt()
    }

    override fun compareTo(other: S2LatLng): Int {
        val latComparison = latRadians.compareTo(other.latRadians)
        return if (latComparison != 0) latComparison else lngRadians.compareTo(other.lngRadians)
    }

    /**
     * Returns true if both the latitude and longitude of the given point are
     * within `maxError` radians of this point.
     */
    /**
     * Returns true if the given point is within `1e-9` radians of this
     * point. This corresponds to a distance of less than `1cm` at the
     * surface of the Earth.
     */
    @JvmOverloads
    fun approxEquals(o: S2LatLng, maxError: Double = 1e-9): Boolean {
        return (abs(latRadians - o.latRadians) < maxError && abs(lngRadians - o.lngRadians) < maxError)
    }

    override fun toString(): String {
        return "(${latDegrees()}, ${lngDegrees()})"
    }

    fun toStringDegrees(): String {
        return "(" + latDegrees() + ", " + lngDegrees() + ")"
    }

    companion object {

        /** The center point the lat/lng coordinate system.  */
        @JvmStatic
        val center = S2LatLng(0.0, 0.0)

        /**
         * Creates a S2LatLng instance from a lat S1Angle and lng S1Angle.
         *
         * @param lat The latitude angle.
         * @param lng The longitude angle.
         * @return The S2LatLng instance.
         */
        @JvmStatic
        fun fromLatLng(lat: S1Angle, lng: S1Angle): S2LatLng = S2LatLng(lat.radians, lng.radians)

        /**
         * Convert a point (not necessarily normalized) to an S2LatLng.
         * We use atan2 to compute the latitude because the input vector is not necessarily unit length, and atan2 is
         * much more accurate than asin near the poles.
         * Note that atan2(0, 0) is defined to be zero.
         *
         * @param p The point to convert.
         * @return A S2LatLng instance that represents the given xyz point.
         */
        @JvmStatic
        @Strictfp
        fun fromPoint(p: S2Point) : S2LatLng = S2LatLng(
                latRadians = FastMath.atan2(p.z, FastMath.sqrt(p.x * p.x + p.y * p.y)),
                lngRadians = FastMath.atan2(p.y, p.x)
        )

        /** A S2LatLng for which isValid will return false. */
        @JvmStatic
        val invalid = S2LatLng(latRadians = M_PI, lngRadians =  2 * M_PI)

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in radians.
         * @param latRadians The latitude of the point in radians.
         * @param lngRadians The longitude of the point in radians.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        fun fromRadians(latRadians: Double, lngRadians: Double): S2LatLng {
            return S2LatLng(latRadians, lngRadians)
        }

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in degrees.
         * @param latDegrees The latitude of the point in degrees.
         * @param lngDegrees The longitude of the point in degrees.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        @Strictfp
        fun fromDegrees(latDegrees: Double, lngDegrees: Double): S2LatLng {
            return fromLatLng(degrees(latDegrees), degrees(lngDegrees))
        }

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in degrees.
         * @param latDegrees The latitude of the point in degrees.
         * @param lngDegrees The longitude of the point in degrees.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        fun fromDegrees(latDegrees: Int, lngDegrees: Int): S2LatLng {
            return fromLatLng(degrees(latDegrees.toDouble()), degrees(lngDegrees.toDouble()))
        }

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in 1e5 degrees.
         * @param latE5 The latitude of the point in 1e5 degrees.
         * @param lngE5 The longitude of the point in 1e5 degrees.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        fun fromE5(latE5: Int, lngE5: Int): S2LatLng {
            return fromLatLng(e5(latE5), e5(lngE5))
        }

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in 1e6 degrees.
         * @param latE6 The latitude of the point in 1e5 degrees.
         * @param lngE6 The longitude of the point in 1e5 degrees.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        fun fromE6(latE6: Int, lngE6: Int): S2LatLng {
            return fromLatLng(e6(latE6), e6(lngE6))
        }

        /**
         * Factory that builds a S2LatLng instance from latitude and longitude in 1e7 degrees.
         * @param latE7 The latitude of the point in 1e5 degrees.
         * @param lngE7 The longitude of the point in 1e5 degrees.
         * @return A S2LatLng instance.
         */
        @JvmStatic
        fun fromE7(latE7: Int, lngE7: Int): S2LatLng {
            return fromLatLng(e7(latE7), e7(lngE7))
        }

        /**
         * Compute the latitude of a given point.
         *
         *  lat = atan(z / sqrt(x^2 + y^2))
         *
         * @param p A point (not necessarily on the unit sphere).
         * @return The latitude of the point.
         */
        @JvmStatic
        fun latitude(p: S2Point): S1Angle {
            // We use atan2 rather than asin because the input vector is not necessarily unit length, and atan2 is much
            // more accurate than asin near the poles.
            return radians(atan2(p[2], sqrt(p[0] * p[0] + p[1] * p[1])))
        }

        /**
         * Compute the longitude of a given point.
         *
         *  lng = atan(y / x)
         *
         * @param p A point (not necessarily on the unit sphere).
         * @return The latitude of the point.
         */
        @JvmStatic
        fun longitude(p: S2Point): S1Angle {
            // Note that atan2(0, 0) is defined to be zero.
            return radians(atan2(p[1], p[0]))
        }

        /**
         *
         */
        operator fun Double.times(p: S2LatLng): S2LatLng = p * this

    }
}

