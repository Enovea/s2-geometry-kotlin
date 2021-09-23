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
import org.apache.commons.math3.util.FastMath.atan2
import org.apache.commons.math3.util.FastMath.cos
import org.apache.commons.math3.util.FastMath.min
import org.apache.commons.math3.util.FastMath.sin
import tech.units.indriya.quantity.Quantities
import tech.units.indriya.unit.Units.METRE
import javax.measure.MetricPrefix.KILO
import javax.measure.Quantity
import javax.measure.quantity.Length

/**
 * The earth modeled as a sphere.
 *
 * This class is a port of the s2earth methods of the Google S2 Geometry project (https://github.com/google/s2geometry).
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
object S2Earth {

    // These functions convert between distances on the unit sphere
    // (expressed as angles subtended from the sphere's center) and
    // distances on the Earth's surface.  This is possible only because
    // the Earth is modeled as a sphere; otherwise a given angle would
    // correspond to a range of distances depending on where the
    // corresponding line segment was located.
    //
    // Note that you will lose precision if you use the ToDistance() method,
    // since Meters is a single-precision type.  If you need more precision,
    // use one of the direct conversion methods below.
    fun toAngle(distance: Quantity<Length>): S1Angle = S1Angle.radians(toRadians(distance))
    fun toChordAngle(distance: Quantity<Length>): S1ChordAngle = S1ChordAngle(toAngle(distance))
    fun toDistance(angle: S1Angle): Quantity<Length> = Quantities.getQuantity(toMeters(angle), METRE)
    fun toDistance(cangle: S1ChordAngle): Quantity<Length> = Quantities.getQuantity(toMeters(cangle), METRE)

    // Convenience functions.  These methods also return a double-precision
    // result, unlike the generic ToDistance() method.
    fun toRadians(distance: Quantity<Length>): Double = distance.to(METRE).value.toDouble() / radiusMeters
    fun toMeters(angle: S1Angle): Double = angle.radians * radiusMeters
    fun toMeters(cangle: S1ChordAngle): Double = toMeters(cangle.toAngle())
    fun toKm(angle: S1Angle): Double = angle.radians * radiusKm
    fun toKm(cangle: S1ChordAngle): Double = toKm(cangle.toAngle())
    fun kmToRadians(km: Double) = km / radiusKm
    fun radiansToKm(radians: Double) = radians * radiusKm
    fun metersToRadians(meters: Double) = meters / radiusMeters
    fun radiansToMeters(radians: Double) = radians * radiusMeters

    // These functions convert between areas on the unit sphere and areas on the Earth's surface.
    // Note that the area of a region on the unit sphere is equal to the
    // solid angle it subtends from the sphere's center (measured in steradians).
    fun squareKmToSteradians(km2: Double): Double = km2 / (radiusKm * radiusKm)
    fun squareMetersToSteradians(m2: Double): Double = m2 / (radiusMeters * radiusMeters)
    fun steradiansToSquareKm(steradians: Double): Double = steradians * radiusKm * radiusKm
    fun steradiansToSquareMeters(steradians: Double): Double = steradians * radiusMeters * radiusMeters

    // Convenience function for the frequent case where you need to call
    // ToRadians in order to convert an east-west distance on the globe to
    // radians. The output is a function of how close to the poles you are
    // (i.e. at the bulge at the equator, one unit of longitude represents a
    // much farther distance). The function will never return more than 2*PI
    // radians, even if you're trying to go 100 million miles west at the north
    // pole.
    fun toLongitudeRadians(distance: Quantity<Length>, latitudeRadians: Double): Double {
        val scalar = cos(latitudeRadians)
        return if (scalar == 0.0) M_PI * 2
        else min(toRadians(distance) / scalar, M_PI * 2)
    }

    // Computes the initial bearing from a to b. This is the bearing an observer
    // at point a has when facing point b. A bearing of 0 degrees is north, and it
    // increases clockwise (90 degrees is east, etc).
    // If a == b, a == -b, or a is one of the Earths' poles, the return value is
    // undefined.
    fun getInitialBearing(a: S2LatLng, b: S2LatLng): S1Angle {
        val lat1 = a.lat().radians
        val cosLat2 = cos(b.lat().radians)
        val latDiff = b.lat().radians - a.lat().radians
        val lngDiff = b.lng().radians - a.lng().radians
        val x = sin(latDiff) + sin(lat1) * cosLat2 * 2 * haversine(lngDiff)
        val y = sin(lngDiff) * cosLat2
        return S1Angle.radians(atan2(y, x))
    }

    // Returns the distance between two points.  Example:
    // double miles = Miles(geostore::S2Earth::GetDistance(a, b)).value();
    //
    // Note that these methods only have single-precision accuracy, since
    // Meters is a single-precision type.  If you ned more precision, use one
    // of the methods below.
    fun getDistance(a: S2Point, b: S2Point): Quantity<Length> = toDistance(S1Angle(a, b))
    fun getDistance(a: S2LatLng, b: S2LatLng): Quantity<Length> = toDistance(a.getDistance(b))

    // Convenience functions.  These methods also return a double-precision
    // result, unlike the generic GetDistance() method.
    fun getDistanceKm(a: S2Point, b: S2Point): Double = getDistance(a, b).to(KILO(METRE)).value.toDouble()
    fun getDistanceKm(a: S2LatLng, b: S2LatLng): Double = getDistance(a, b).to(KILO(METRE)).value.toDouble()
    fun getDistanceMeters(a: S2Point, b: S2Point) = getDistance(a, b).to(METRE).value.toDouble()
    fun getDistanceMeters(a: S2LatLng, b: S2LatLng): Double = getDistance(a, b).to(METRE).value.toDouble()

    // Returns the Earth's mean radius, which is the radius of the equivalent
    // sphere with the same surface area.  According to NASA, this value is
    // 6371.01 +/- 0.02 km.  The equatorial radius is 6378.136 km, and the polar
    // radius is 6356.752 km.  They differ by one part in 298.257.
    //
    // Reference: http://ssd.jpl.nasa.gov/phys_props_earth.html, which quotes
    // Yoder, C.F. 1995. "Astrometric and Geodetic Properties of Earth and the
    // Solar System" in Global Earth Physics, A Handbook of Physical Constants,
    // AGU Reference Shelf 1, American Geophysical Union, Table 2.
    val radius = Quantities.getQuantity(6371010.0, METRE)

    // Convenience functions.
    val radiusKm: Double = radius.to(KILO(METRE)).value.toDouble()
    val radiusMeters: Double = radius.to(METRE).value.toDouble()

    // Returns the altitude of the lowest known point on Earth. The lowest known
    // point on Earth is the Challenger Deep with an altitude of -10898 meters
    // above the surface of the spherical earth.
    val lowestAltitude = Quantities.getQuantity(-10898, METRE)

    // Convenience functions.
    val lowestAltitudeKm = lowestAltitude.to(KILO(METRE)).value.toDouble()
    val lowestAltitudeMeters = lowestAltitude.to(METRE).value.toDouble()

    // Returns the altitude of the highest known point on Earth. The highest
    // known point on Earth is Mount Everest with an altitude of 8846 meters
    // above the surface of the spherical earth.
    val highestAltitude = Quantities.getQuantity(8846, METRE)

    // Convenience functions.
    val highestAltitudeKm = highestAltitude.to(KILO(METRE)).value.toDouble()
    val highestAltitudeMeters = highestAltitude.to(METRE).value.toDouble()

}


operator fun <Q : Quantity<Q>> Quantity<Q>.div(other: Quantity<Q>): Quantity<*> = this.divide(other)

// http://en.wikipedia.org/wiki/Haversine_formula
// Haversine(x) has very good numerical stability around zero.
// Haversine(x) == (1-cos(x))/2 == sin(x/2)^2; must be implemented with the
// second form to reap the numerical benefits.
fun haversine(radians: Double): Double {
    val sinHalf = sin(radians / 2)
    return sinHalf * sinHalf
}
