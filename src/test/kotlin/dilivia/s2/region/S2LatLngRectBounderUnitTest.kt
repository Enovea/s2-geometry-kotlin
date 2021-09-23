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
import dilivia.math.M_PI
import dilivia.math.M_PI_2
import dilivia.math.M_PI_4
import dilivia.math.R1Interval
import dilivia.math.vectors.times
import dilivia.s2.S1Angle
import dilivia.s2.S1Interval
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import dilivia.s2.edge.S2EdgeDistances
import org.assertj.core.api.Assertions
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import java.lang.Math.pow
import kotlin.math.asin
import kotlin.math.pow
import kotlin.math.sqrt


class S2LatLngRectBounderUnitTest {
  
  fun getEdgeBound(a: S2Point, b: S2Point): S2LatLngRect {
    val bounder = S2LatLngRectBounder()
    bounder.addPoint(a);
    bounder.addPoint(b);
    return bounder.getBound()
  }

  fun getEdgeBound(x1: Double, y1: Double, z1: Double, x2: Double, y2: Double, z2: Double): S2LatLngRect {
    return getEdgeBound(S2Point(x1, y1, z1).normalize(), S2Point(x2, y2, z2).normalize());
  }

  fun getEdgeBound(x1: Int, y1: Int, z1: Int, x2: Int, y2: Int, z2: Int): S2LatLngRect {
    return getEdgeBound(S2Point(x1, y1, z1).normalize(), S2Point(x2, y2, z2).normalize());
  }

  val kRectError = S2LatLngRectBounder.maxErrorForTests()

  @Test
  fun rectBounderMaxLatitudeSimple() {
    // Check cases where the min/max latitude is attained at a vertex.
    val kCubeLat = asin (1 / sqrt(3.0));  // 35.26 degrees
    assertThat(getEdgeBound(1, 1, 1, 1, -1, -1).approxEquals(S2LatLngRect(R1Interval(-kCubeLat, kCubeLat), S1Interval(-M_PI_4, M_PI_4)), kRectError)).isTrue()    // NOLINT        
    assertThat(getEdgeBound(1, -1, 1, 1, 1, -1).approxEquals(S2LatLngRect(R1Interval(-kCubeLat, kCubeLat), S1Interval(-M_PI_4, M_PI_4)), kRectError)).isTrue()   // NOLINT        

    // Check cases where the min/max latitude occurs in the edge interior.
    // These tests expect the result to be pretty close to the middle of the
    // allowable error range (i.e., by adding 0.5 * kRectError).

    // Max latitude, CW edge
    assertThat(getEdgeBound(1, 1, 1, 1, -1, 1).lat.hi     ).isEqualTo(M_PI_4 + 0.5 * kRectError.lat().radians );  // NOLINT
    // Max latitude, CCW edge
    assertThat(getEdgeBound(1, -1, 1, 1, 1, 1).lat.hi     ).isEqualTo(M_PI_4 + 0.5 * kRectError.lat().radians );  // NOLINT
    // Min latitude, CW edge
    assertThat(getEdgeBound(1, -1, -1, -1, -1, -1).lat.lo ).isEqualTo(-M_PI_4 - 0.5 * kRectError.lat().radians);  // NOLINT
    // Min latitude, CCW edge
    assertThat(getEdgeBound(-1, 1, -1, -1, -1, -1).lat.lo ).isEqualTo(-M_PI_4 - 0.5 * kRectError.lat().radians);  // NOLINT

    // Check cases where the edge passes through one of the poles.
    assertThat(getEdgeBound(0.3, 0.4, 1.0, -0.3, -0.4, 1.0).lat.hi  ).isEqualTo(M_PI_2 );  // NOLINT
    assertThat(getEdgeBound(0.3, 0.4, -1.0, -0.3, -0.4, -1.0).lat.lo).isEqualTo(-M_PI_2);  // NOLINT
  }

  @Test
  fun rectBounderMaxLatitudeRandom() {
    // Check that the maximum latitude of edges is computed accurately to within
    // 3 * DoubleType.epsilon (the expected maximum error).  We concentrate on maximum
    // latitudes near the equator and north pole since these are the extremes.

    val kIters = 100;
    repeat(kIters) {
      // Construct a right-handed coordinate frame (U,V,W) such that U points
      // slightly above the equator, V points at the equator, and W is slightly
      // offset from the north pole.
      val randomPoint = S2Random.randomPoint()
      randomPoint[2] = DoubleType.epsilon * 1e-6 * 1e12.pow(S2Random.randomDouble());  // log is uniform
      val u = randomPoint.normalize()
      val v = S2PointUtil.robustCrossProd(S2Point(0, 0, 1), u).normalize()
      val w = S2PointUtil.robustCrossProd(u, v).normalize()

      // Construct a line segment AB that passes through U, and check that the
      // maximum latitude of this segment matches the latitude of U.
      val a = (u - S2Random.randomDouble() * v).normalize()
      val b = (u + S2Random.randomDouble() * v).normalize()
      val ab_bound = getEdgeBound(a, b);
      assertThat(ab_bound.lat.hi).isCloseTo(S2LatLng.latitude(u).radians, Offset.offset(kRectError.lat().radians))

      // Construct a line segment CD that passes through W, and check that the
      // maximum latitude of this segment matches the latitude of W.
      val c = (w - S2Random.randomDouble() * v).normalize();
      val d = (w + S2Random.randomDouble() * v).normalize();
      val cd_bound = getEdgeBound(c, d);
      assertThat(cd_bound.lat.hi).isCloseTo(S2LatLng.latitude(w).radians, Offset.offset(kRectError.lat().radians))
    }
  }

  fun perturbATowardsB(a: S2Point, b: S2Point): S2Point {
    val choice = S2Random.randomDouble();
    if (choice < 0.1) {
      return a;
    }
    if (choice < 0.3) {
      // Return a point that is exactly proportional to A and that still
      // satisfies S2Point.isUnitLength().
      while (true) {
        val b =(2 - a.norm() + 5 * (S2Random.randomDouble() -0.5) * DoubleType.epsilon) * a;
        if (b != a && S2PointUtil.isUnitLength(b))
          return b;
      }
    }
    if (choice < 0.5) {
      // Return a point such that the distance squared to A will underflow.
      return S2EdgeDistances.interpolateAtDistance(S1Angle.radians(1e-300), a, b);
    }
    // Otherwise return a point whose distance from A is near DoubleType.epsilon such
    // that the log of the pdf is uniformly distributed.
    val distance = DoubleType.epsilon * 1e-5 * pow(1e6, S2Random.randomDouble());
    return S2EdgeDistances.interpolateAtDistance(S1Angle.radians(distance), a, b);
  }

  fun randomPole(): S2Point {
    return S2Point(0, 0, if(S2Random.oneIn(2)) 1 else -1);
  }

  fun pointNearPole(): S2Point {
    return perturbATowardsB(randomPole(), S2Random.randomPoint())
  }

  fun pointNearEquator(): S2Point {
    return perturbATowardsB(S2Point(S2Random.randomDouble(), S2Random.randomDouble(), 0.0).normalize(), randomPole());
  }

  @Test
  fun rectBounderNearlyIdenticalOrAntipodalPoints() {
    // Test pairs of points that are either:
    //  - identical
    //  - nearly or exactly proportional, e.g. (1,0,0) vs. (1+2e-16, 0, 0)
    //  - very close to each other
    // Furthermore we want to test cases where the two points are:
    //  - on a nearly-polar great circle
    //  - on a nearly-equatorial great circle
    //  - near the poles, but on any great circle
    //  - near the equator, but on any great circle
    //  - positioned arbitrarily
    // Also test the corresponding situations for antipodal points, i.e. by
    // negating one of the points so that they are almost 180 degrees apart.
    val kIters = 10000;
    repeat(kIters) {
    lateinit var a: S2Point
      lateinit var b: S2Point
    when(S2Random.randomInt(5)) {
    0 -> {
      // Two nearby points on a nearly-polar great circle.
      a = S2Random.randomPoint()
      b = perturbATowardsB(a, pointNearPole())
    }
    1 -> {
      // Two nearby points on a nearly-equatorial great circle.
      a = pointNearEquator();
      b = perturbATowardsB(a, pointNearEquator());
    }
    2 -> {
      // Two nearby points near a pole, but on any great circle.
      a = pointNearPole();
      b = perturbATowardsB(a, S2Random.randomPoint());
    }
    3 -> {
      // Two nearby points near the equator, but on any great circle.
      a = pointNearEquator();
      b = perturbATowardsB(a, S2Random.randomPoint());
    }
    4 -> {
      // Two nearby points anywhere on the sphere.
      a = S2Random.randomPoint();
      b = perturbATowardsB(a, S2Random.randomPoint());
    }
      else -> Assertions.fail("Error")
  }
    // The two points are chosen to be so close to each other that the min/max
    // latitudes are nearly always achieved at the edge endpoints.  The only
    // thing we need to watch out for is that the latitude error bound is
    // slightly larger if the min/max latitude occurs in the edge interior.
    val expected_bound = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(a), S2LatLng.fromPoint(b));
    val bound = getEdgeBound (a, b)
    assertThat(bound.contains(expected_bound)).isTrue();
    assertThat(expected_bound.expanded(kRectError).polarClosure().contains(bound)).isTrue();

    // If the two points are close enough and one point is negated (antipodal
    // points), the bound should be the entire sphere.
    if ((a - b).crossProd(a + b).norm() <= 6.110 * DoubleType.epsilon) {
      assertThat(getEdgeBound(a, -b)).isEqualTo(S2LatLngRect.full())
    }
  }
  }

  fun getSubregionBound(x_lat: Double, x_lng: Double, y_lat: Double, y_lng: Double): S2LatLngRect {
    val input = S2LatLngRect.fromPointPair(S2LatLng.fromRadians(x_lat, x_lng), S2LatLng.fromRadians(y_lat, y_lng));
    val output = S2LatLngRectBounder.expandForSubregions(input)

    // Test that the bound is actually expanded.
    assertThat(output.contains(input)).isTrue();
    if ( input . lat  == S2LatLngRect.fullLat()) {
    assertThat(input.lat.contains(output.lat)).isFalse()
  }
    return output;
  }


  @Test
  fun rectBounderExpandForSubregions() {
    // First we check the various situations where the bound contains
    // nearly-antipodal points.  The tests are organized into pairs where the
    // two bounds are similar except that the first bound meets the
    // nearly-antipodal criteria while the second does not.

    // Cases where the bound does not straddle the equator (but almost does),
    // and spans nearly 180 degrees in longitude.
    assertThat(getSubregionBound(3e-16, 0.0, 1e-14, M_PI).isFull).isTrue();
    assertThat(getSubregionBound(9e-16, 0.0, 1e-14, M_PI).isFull).isFalse()
    assertThat(getSubregionBound(1e-16, 7e-16, 1e-14, M_PI).isFull).isTrue();
    assertThat(getSubregionBound(3e-16, 14e-16, 1e-14, M_PI).isFull).isFalse()
    assertThat(getSubregionBound(1e-100, 14e-16, 1e-14, M_PI).isFull).isTrue();
    assertThat(getSubregionBound(1e-100, 22e-16, 1e-14, M_PI).isFull).isFalse()

    // Cases where the bound spans at most 90 degrees in longitude, and almost
    // 180 degrees in latitude.  Note that DoubleType.epsilon is about 2.22e-16, which
    // implies that the double-precision value just below Pi/2 can be written as
    // (M_PI_2 - 2e-16).
    assertThat(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 7e-16, 0.0).isFull).isTrue();
    assertThat(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 30e-16, 0.0).isFull).isFalse()
    assertThat(getSubregionBound(-M_PI_2 + 4e-16, 0.0, M_PI_2 - 2e-16, 1e-7).isFull).isTrue();
    assertThat(getSubregionBound(-M_PI_2 + 30e-16, 0.0, M_PI_2, 1e-7).isFull).isFalse()
    assertThat(getSubregionBound(-M_PI_2 + 4e-16, 0.0, M_PI_2 - 4e-16, M_PI_2).isFull).isTrue();
    assertThat(getSubregionBound(-M_PI_2, 0.0, M_PI_2 - 30e-16, M_PI_2).isFull).isFalse()

    // Cases where the bound straddles the equator and spans more than 90
    // degrees in longitude.  These are the cases where the critical distance is
    // between a corner of the bound and the opposite longitudinal edge.  Unlike
    // the cases above, here the bound may contain nearly-antipodal points (to
    // within 3.055 * DoubleType.epsilon) even though the latitude and longitude ranges
    // are both significantly less than (Pi - 3.055 * DoubleType.epsilon).
    assertThat(getSubregionBound(-M_PI_2, 0.0, M_PI_2 - 1e-8, M_PI - 1e-7).isFull).isTrue();
    assertThat(getSubregionBound(-M_PI_2, 0.0, M_PI_2 - 1e-7, M_PI - 1e-7).isFull).isFalse()
    assertThat(getSubregionBound(-M_PI_2 + 1e-12, -M_PI + 1e-4, M_PI_2, 0.0).isFull).isTrue();
    assertThat(getSubregionBound(-M_PI_2 + 1e-11, -M_PI + 1e-4, M_PI_2, 0.0).isFull).isTrue();

    // Now we test cases where the bound does not contain nearly-antipodal
    // points, but it does contain points that are approximately 180 degrees
    // apart in latitude.
    assertThat(getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 2e-16).approxEquals(S2LatLngRect(R1Interval(1.5, 1.5), S1Interval.full()), kRectError)).isTrue();
    assertThat(getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 7e-16).approxEquals(S2LatLngRect(R1Interval(1.5, 1.5), S1Interval(-M_PI_2, M_PI_2 - 7e-16)), kRectError)).isTrue();

    // Test the full and empty bounds.
    assertThat(S2LatLngRectBounder.expandForSubregions(S2LatLngRect.full()).isFull).isTrue();
    assertThat(S2LatLngRectBounder.expandForSubregions(S2LatLngRect.empty()).isEmpty).isTrue();

    // Check for cases where the bound is expanded to include one of the poles.
    assertThat(getSubregionBound(-M_PI_2 + 1e-15, 0.0, -M_PI_2 + 1e-15, 0.0).approxEquals(S2LatLngRect(R1Interval(-M_PI_2, -M_PI_2 + 1e-15), S1Interval.full()), kRectError)).isTrue();
    assertThat(getSubregionBound(M_PI_2 - 1e-15, 0.0, M_PI_2 - 1e-15, 0.0).approxEquals(S2LatLngRect(R1Interval(M_PI_2 - 1e-15, M_PI_2), S1Interval.full()), kRectError)).isTrue();
  }

}
