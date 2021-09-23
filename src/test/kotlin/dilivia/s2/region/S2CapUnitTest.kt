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

import dilivia.math.DoubleType
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng
import dilivia.s2.S2LatLng.Companion.fromDegrees
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.coords.S2Coords
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import kotlin.math.atan
import kotlin.math.sqrt
import kotlin.random.Random

class S2CapUnitTest {

    private fun getLatLngPoint(latDegrees: Double, lngDegrees: Double): S2Point {
        return fromDegrees(latDegrees, lngDegrees).toPoint()
    }

    private fun getLatLngPoint(latDegrees: Int, lngDegrees: Int): S2Point = getLatLngPoint(latDegrees.toDouble(), lngDegrees.toDouble())

    @Test
    fun basic() {
        // Test basic properties of empty and full caps.
        val empty = S2Cap.empty
        val full = S2Cap.full
        assertThat(empty.isValid).isTrue()
        assertThat(empty.isEmpty).isTrue()
        assertThat(empty.complement.isFull).isTrue()
        assertThat(full.isValid).isTrue()
        assertThat(full.isFull).isTrue()
        assertThat(full.complement.isEmpty).isTrue()
        assertThat(full.height).isEqualTo(2.0)
        assertThat(full.radius().degrees()).isEqualTo(180.0)

        // Test the S1Angle constructor using out-of-range arguments.
        assertThat(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(-20)).isEmpty).isTrue()
        assertThat(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(5)).isFull).isTrue()
        assertThat(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.infinity()).isFull).isTrue()

        // Check that the default S2Cap is identical to Empty().
        val defaultEmpty = S2Cap()
        assertThat(defaultEmpty.isValid).isTrue()
        assertThat(defaultEmpty.isEmpty).isTrue()
        assertThat(defaultEmpty.center).isEqualTo(empty.center)
        assertThat(defaultEmpty.height).isEqualTo(empty.height)

        // Containment and intersection of empty and full caps.
        assertThat(empty.contains(empty)).isTrue()
        assertThat(full.contains(empty)).isTrue()
        assertThat(full.contains(full)).isTrue()
        assertThat(empty.interiorIntersects(empty)).isFalse()
        assertThat(full.interiorIntersects(full)).isTrue()
        assertThat(full.interiorIntersects(empty)).isFalse()

        // Singleton cap containing the x-axis.
        val xaxis = S2Cap.fromPoint(S2Point(1, 0, 0))
        assertThat(xaxis.contains(S2Point(1, 0, 0))).isTrue()
        assertThat(xaxis.contains(S2Point(1.0, 1e-20, 0.0))).isFalse()
        assertThat(xaxis.radius().radians).isEqualTo(0.0)

        // Singleton cap containing the y-axis.
        val yaxis = S2Cap.fromPoint(S2Point(0, 1, 0))
        assertThat(yaxis.contains(xaxis.center)).isFalse()
        assertThat(xaxis.height).isEqualTo(0.0)

        // Check that the complement of a singleton cap is the full cap.
        val xcomp = xaxis.complement
        assertThat(xcomp.isValid).isTrue()
        assertThat(xcomp.isFull).isTrue()
        assertThat(xcomp.contains(xaxis.center)).isTrue()

        // Check that the complement of the complement is *not* the original.
        assertThat(xcomp.complement.isValid).isTrue()
        assertThat(xcomp.complement.isEmpty).isTrue()
        assertThat(xcomp.complement.contains(xaxis.center)).isFalse()

        // Check that very small caps can be represented accurately.
        // Here "kTinyRad" is small enough that unit vectors perturbed by this
        // amount along a tangent do not need to be renormalized.
        val kTinyRad = 1e-10
        val tiny = S2Cap.fromCenterAngle(S2Point(1, 2, 3).normalize(), S1Angle.radians(kTinyRad))
        val tangent = tiny.center.crossProd(S2Point(3, 2, 1)).normalize()
        assertThat(tiny.contains(tiny.center + (0.99 * kTinyRad) * tangent)).isTrue()
        assertThat(tiny.contains(tiny.center + 1.01 * kTinyRad * tangent)).isFalse()

        // Basic tests on a hemispherical cap.
        val hemi = S2Cap.fromCenterHeight(S2Point(1, 0, 1).normalize(), 1.0)
        assertThat(hemi.complement.center).isEqualTo(-hemi.center)
        assertThat(hemi.complement.height).isEqualTo(1.0)
        assertThat(hemi.contains(S2Point(1, 0, 0))).isTrue()
        assertThat(hemi.complement.contains(S2Point(1, 0, 0))).isFalse()
        assertThat(hemi.contains(S2Point(1.0, 0.0, -(1 - kEps)).normalize())).isTrue()
        assertThat(hemi.interiorContains(S2Point(1.0, 0.0, -(1 + kEps)).normalize())).isFalse()

        // A concave cap.  Note that the error bounds for point containment tests
        // increase with the cap angle, so we need to use a larger error bound
        // here.  (It would be painful to do this everywhere, but this at least
        // gives an example of how to compute the maximum error.)
        val center = getLatLngPoint(80.0, 10.0)
        val radius = S1ChordAngle(S1Angle.degrees(150))
        val maxError = (radius.getS2PointConstructorMaxError() + radius.getS1AngleConstructorMaxError() + 3 * DoubleType.epsilon)  // getLatLngPoint() error
        val concave = S2Cap(center, radius)
        val concaveMin = S2Cap(center, radius.plusError(-maxError))
        val concaveMax = S2Cap(center, radius.plusError(maxError))
        assertThat(concaveMax.contains(getLatLngPoint(-70, 10))).isTrue()
        assertThat(concaveMin.contains(getLatLngPoint(-70, 10))).isFalse()
        assertThat(concaveMax.contains(getLatLngPoint(-50, -170))).isTrue()
        assertThat(concaveMin.contains(getLatLngPoint(-50, -170))).isFalse()

        // Cap containment tests.
        assertThat(empty.contains(xaxis)).isFalse()
        assertThat(empty.interiorIntersects(xaxis)).isFalse()
        assertThat(full.contains(xaxis)).isTrue()
        assertThat(full.interiorIntersects(xaxis)).isTrue()
        assertThat(xaxis.contains(full)).isFalse()
        assertThat(xaxis.interiorIntersects(full)).isFalse()
        assertThat(xaxis.contains(xaxis)).isTrue()
        assertThat(xaxis.interiorIntersects(xaxis)).isFalse()
        assertThat(xaxis.contains(empty)).isTrue()
        assertThat(xaxis.interiorIntersects(empty)).isFalse()
        assertThat(hemi.contains(tiny)).isTrue()
        assertThat(hemi.contains(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(M_PI_4 - kEps)))).isTrue()
        assertThat(hemi.contains(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(M_PI_4 + kEps)))).isFalse()
        assertThat(concave.contains(hemi)).isTrue()
        assertThat(concave.interiorIntersects(hemi.complement)).isTrue()
        assertThat(concave.contains(S2Cap.fromCenterHeight(-concave.center, 0.1))).isFalse()
    }

    @Test
    fun addEmptyCapToNonEmptyCap() {
        val nonEmptyCap = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(10))
        val initialArea = nonEmptyCap.area
        assertThat(nonEmptyCap.addCap(S2Cap.empty).area).isEqualTo(initialArea)
    }

    @Test
    fun addNonEmptyCapToEmptyCap() {
        val empty = S2Cap.empty
        val nonEmptyCap = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(10))
        assertThat(empty.addCap(nonEmptyCap).area).isEqualTo(nonEmptyCap.area)
    }

    @Test
    fun getRectBound() {
        // Empty and full caps.
        assertThat(S2Cap.empty.rectBound.isEmpty).isTrue()
        assertThat(S2Cap.full.rectBound.isFull).isTrue()

        val kDegreeEps = 1e-13
        // Maximum allowable error for latitudes and longitudes measured in
        // degrees.  (assertEquals isn't sufficient.)

        // Cap that includes the south pole.
        var rect = S2Cap.fromCenterAngle(getLatLngPoint(-45, 57), S1Angle.degrees(50)).rectBound
        assertThat(rect.latLo().degrees()).isCloseTo(-90.0, Offset.offset(kDegreeEps))
        assertThat(rect.latHi().degrees()).isCloseTo(5.0, Offset.offset(kDegreeEps))
        assertThat(rect.lng.isFull).isTrue()

        // Cap that is tangent to the north pole.
        rect = S2Cap.fromCenterAngle(S2Point(1, 0, 1).normalize(), S1Angle.radians(M_PI_4 + 1e-16)).rectBound
        assertThat(rect.lat.lo).isCloseTo(0.0, Offset.offset(kEps))
        assertThat(rect.lat.hi).isCloseTo(M_PI_2, Offset.offset(kEps))
        assertThat(rect.lng.isFull).isTrue()

        rect = S2Cap.fromCenterAngle(S2Point(1, 0, 1).normalize(), S1Angle.degrees(45 + 5e-15)).rectBound
        assertThat(rect.latLo().degrees()).isCloseTo(0.0, Offset.offset(kDegreeEps))
        assertThat(rect.latHi().degrees()).isCloseTo(90.0, Offset.offset(kDegreeEps))
        assertThat(rect.lng.isFull).isTrue()

        // The eastern hemisphere.
        rect = S2Cap.fromCenterAngle(S2Point(0, 1, 0), S1Angle.radians(M_PI_2 + 2e-16)).rectBound
        assertThat(rect.latLo().degrees()).isCloseTo(-90.0, Offset.offset(kDegreeEps))
        assertThat(rect.latHi().degrees()).isCloseTo(90.0, Offset.offset(kDegreeEps))
        assertThat(rect.lng.isFull).isTrue()

        // A cap centered on the equator.
        rect = S2Cap.fromCenterAngle(getLatLngPoint(0, 50), S1Angle.degrees(20)).rectBound
        assertThat(rect.latLo().degrees()).isCloseTo(-20.0, Offset.offset(kDegreeEps))
        assertThat(rect.latHi().degrees()).isCloseTo(20.0, Offset.offset(kDegreeEps))
        assertThat(rect.lngLo().degrees()).isCloseTo(30.0, Offset.offset(kDegreeEps))
        assertThat(rect.lngHi().degrees()).isCloseTo(70.0, Offset.offset(kDegreeEps))

        // A cap centered on the north pole.
        rect = S2Cap.fromCenterAngle(getLatLngPoint(90, 123), S1Angle.degrees(10)).rectBound
        assertThat(rect.latLo().degrees()).isCloseTo(80.0, Offset.offset(kDegreeEps))
        assertThat(rect.latHi().degrees()).isCloseTo(90.0, Offset.offset(kDegreeEps))
        assertThat(rect.lng.isFull).isTrue()
    }

    @Test
    fun cellMethods() {
        // For each cube face, we construct some cells on
        // that face and some caps whose positions are relative to that face,
        // and then check for the expected intersection/containment results.

        // The distance from the center of a face to one of its vertices.
        val kFaceRadius = atan(sqrt(2.0))

        for (face in 0..5) {
            // The cell consisting of the entire face.
            val rootCell = S2Cell.fromFace(face)

            // A leaf cell at the midpoint of the v=1 edge.
            val edgeCell = S2Cell(S2Coords.faceUvToXyz(face, 0.0, 1 - kEps))

            // A leaf cell at the u=1, v=1 corner.
            val cornerCell = S2Cell(S2Coords.faceUvToXyz(face, 1 - kEps, 1 - kEps))

            // Quick check for full and empty caps.
            assertThat(S2Cap.full.contains(rootCell)).isTrue()
            assertThat(S2Cap.empty.mayIntersect(rootCell)).isFalse()

            // Check intersections with the bounding caps of the leaf cells that are
            // adjacent to 'corner_cell' along the Hilbert curve.  Because this corner
            // is at (u=1,v=1), the curve stays locally within the same cube face.
            var first = cornerCell.id()
            repeat(3) { first = first.previous() }
            var last = cornerCell.id()
            repeat(3) { last = last.next()}
            var id = first
            while (id <= last) {
                val cell = S2Cell(id)
                assertThat(cell.capBound.contains(cornerCell)).isEqualTo(id == cornerCell.id())
                assertThat(cell.capBound.mayIntersect(cornerCell)).isEqualTo(id.parent().contains(cornerCell.id()))
                id = id.next()
            }

            val antiFace = (face + 3) % 6  // Opposite face.
            for (capFace in 0..5) {
                // A cap that barely contains all of 'capFace'.
                val center = S2Coords.norm(capFace)
                val covering = S2Cap.fromCenterAngle(center, S1Angle.radians(kFaceRadius + kEps))
                assertThat(covering.contains(rootCell)).isEqualTo(capFace == face)
                assertThat(covering.mayIntersect(rootCell)).isEqualTo(capFace != antiFace)
                assertThat(covering.contains(edgeCell)).isEqualTo(center.dotProd(edgeCell.getCenter()) > 0.1)
                assertThat(covering.contains(edgeCell)).isEqualTo(covering.mayIntersect(edgeCell))
                assertThat(covering.contains(cornerCell)).isEqualTo(capFace == face)
                assertThat(covering.mayIntersect(cornerCell)).isEqualTo(center.dotProd(cornerCell.getCenter()) > 0)

                // A cap that barely intersects the edges of 'capFace'.
                val bulging = S2Cap.fromCenterAngle(center, S1Angle.radians(M_PI_4 + kEps))
                assertThat(bulging.contains(rootCell)).isFalse()
                assertThat(bulging.mayIntersect(rootCell)).isEqualTo(capFace != antiFace)
                assertThat(bulging.contains(edgeCell))
                    .withFailMessage("capFace = $capFace == $face != $bulging.contains($edgeCell) == ${bulging.contains(edgeCell)}")
                    .isEqualTo(capFace == face)
                assertThat(bulging.mayIntersect(edgeCell)).isEqualTo(center.dotProd(edgeCell.getCenter()) > 0.1)
                assertThat(bulging.contains(cornerCell)).isFalse()
                assertThat(bulging.mayIntersect(cornerCell)).isFalse()

                // A singleton cap.
                val singleton = S2Cap.fromCenterAngle (center, S1Angle.zero())
                assertThat(singleton.mayIntersect(rootCell)).isEqualTo(capFace == face)
                assertThat(singleton.mayIntersect(edgeCell)).isFalse()
                assertThat(singleton.mayIntersect(cornerCell)).isFalse()
            }
        }
    }

    @Test
    fun getCellUnionBoundLevel1Radius() {
        // Check that a cap whose radius is approximately the width of a level 1
        // S2Cell can be covered by only 3 faces.
        val cap = S2Cap.fromCenterAngle(S2Point(1, 1, 1).normalize(), S1Angle.radians(S2Coords.projection.kMinWidth.getValue(1)))
        val covering = mutableListOf<S2CellId>()
        cap.getCellUnionBound(covering)
        assertThat(covering.size).isEqualTo(3)
    }

    @Test
    fun expanded() {
        assertThat(S2Cap.empty.expanded(S1Angle.radians(2)).isEmpty).isTrue()
        assertThat(S2Cap.full.expanded(S1Angle.radians(2)).isFull).isTrue()
        val cap50 = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(50))
        val cap51 = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(51))
        assertThat(cap50.expanded(S1Angle.radians(0)).approxEquals(cap50)).isTrue()
        assertThat(cap50.expanded(S1Angle.degrees(1)).approxEquals(cap51)).isTrue()
        assertThat(cap50.expanded(S1Angle.degrees(129.99)).isFull).isFalse()
        assertThat(cap50.expanded(S1Angle.degrees(130.01)).isFull).isTrue()
    }

    @Test
    fun getCentroid() {
        // Empty and full caps.
        assertThat(S2Cap.empty.centroid).isEqualTo(S2Point())
        assertThat(S2Cap.full.centroid.norm() <= 1e-15).isTrue()

        // Random caps.
        for (i in 0 until 100) {
            val center = S2Random.randomPoint()
            val height = Random.nextDouble(0.0, 2.0)
            val cap = S2Cap.fromCenterHeight(center, height)
            val centroid = cap.centroid
            val expected = center * (1.0 - height / 2.0) * cap.area
            assertThat((expected - centroid).norm() <= 1e-15).isTrue()
        }
    }

    @Test
    fun union() {
        // Two caps which have the same center but one has a larger radius.
        val a = S2Cap.fromCenterAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.2))
        val b = S2Cap.fromCenterAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.3))
        assertThat(b.contains(a)).isTrue()
        assertThat(a.union(b)).isEqualTo(b)

        // Two caps where one is the full cap.
        assertThat(a.union(S2Cap.full).isFull).isTrue()

        // Two caps where one is the empty cap.
        assertThat(a.union(S2Cap.empty)).isEqualTo(a)

        // Two caps which have different centers, one entirely encompasses the other.
        val c = S2Cap.fromCenterAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(1.5))
        assertThat(c.contains(a)).isTrue()
        assertThat(c.center).isEqualTo(a.union(c).center)
        assertThat(c.radius()).isEqualTo(a.union(c).radius())

        // Two entirely disjoint caps.
        val d = S2Cap.fromCenterAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(0.1))
        assertThat(d.contains(a)).isFalse()
        assertThat(d.intersects(a)).isFalse()
        assertThat(a.union(d).approxEquals(d.union(a))).isTrue()
        assertThat(50.4588).isCloseTo(S2LatLng.fromPoint(a.union(d).center).lat().degrees(), Offset.offset(0.001))
        assertThat(S2LatLng.fromPoint(a.union(d).center).lng().degrees()).isCloseTo(10.4525, Offset.offset(0.001))
        assertThat(a.union(d).radius().degrees()).isCloseTo(0.7425, Offset.offset(0.001))

        // Two partially overlapping caps.
        val e = S2Cap.fromCenterAngle(getLatLngPoint(50.3, 10.3), S1Angle.degrees(0.2))
        assertThat(e.contains(a)).isFalse()
        assertThat(e.intersects(a)).isTrue()
        assertThat(a.union(e).approxEquals(e.union(a))).isTrue()
        assertThat(S2LatLng.fromPoint(a.union(e).center).lat().degrees()).isCloseTo(50.1500, Offset.offset(0.001))
        assertThat(S2LatLng.fromPoint(a.union(e).center).lng().degrees()).isCloseTo(10.1495, Offset.offset(0.001))
        assertThat(a.union(e).radius().degrees()).isCloseTo(0.3781, Offset.offset(0.001))

        // Two very large caps, whose radius sums to in excess of 180 degrees, and
        // whose centers are not antipodal.
        val f = S2Cap.fromCenterAngle(S2Point(0, 0, 1).normalize(), S1Angle.degrees(150))
        val g = S2Cap.fromCenterAngle(S2Point(0, 1, 0).normalize(), S1Angle.degrees(150))
        assertThat(f.union(g).isFull).isTrue()

        // Two non-overlapping hemisphere caps with antipodal centers.
        val hemi = S2Cap.fromCenterHeight(S2Point(0, 0, 1).normalize(), 1.0)
        assertThat(hemi.union(hemi.complement).isFull).isTrue()
    }

    /*
    fun testEncodeDecode() {
        S2Cap cap = S2Cap . fromCenterHeight (S2Point(3, 2, 1).normalize(), 1)
        Encoder encoder
                cap.Encode(& encoder)
        Decoder decoder (encoder.base(), encoder.length())
        S2Cap decoded_cap
                assertThat(decoded_cap.Decode(& decoder)).isTrue()
        assertEquals(cap, decoded_cap)
    }
    */

    companion object {
        // About 9 times the double-precision roundoff relative error.
        const val kEps = 1e-15
    }
}
