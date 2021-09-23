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
@file:Suppress("UsePropertyAccessSyntax")

package dilivia.s2.region

import dilivia.math.M_PI
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2Debug
import dilivia.s2.S2Error
import dilivia.s2.S2Factory
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2PointUtil
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.layers.S2PolygonLayer
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.coords.S2Coords
import dilivia.s2.edge.S2EdgeCrossings.kIntersectionMergeRadius
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.shape.S2ClosestEdgeQuery
import dilivia.s2.region.S2CellUnionUnitTest.Companion.checkCovering
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.api.Assertions.fail
import org.assertj.core.data.Offset
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.MethodSource
import java.util.stream.Stream
import kotlin.math.min

// A set of nested loops around the point 0:0 (lat:lng).
// Every vertex of kNear0 is a vertex of kNear1.
private const val kNearPoint: String = "0:0"
private const val kNear0: String = "-1:0, 0:1, 1:0, 0:-1;"
private const val kNear1: String = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;"
private const val kNear2: String = "-1:-2, -2:5, 5:-2;"
private const val kNear3: String = "-2:-2, -3:6, 6:-3;"
private const val kNearHemi: String = "0:-90, -90:0, 0:90, 90:0;"

// A set of nested loops around the point 0:180 (lat:lng).
// Every vertex of kFar0 and kFar2 belongs to kFar1, and all
// the loops except kFar2 are non-convex.
private const val kFar0: String = "0:179, 1:180, 0:-179, 2:-180;"
private const val kFar1: String = "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;"
private const val kFar2: String = "3:-178, 3:178, -1:179, -1:-179;"
private const val kFar3: String = "-3:-178, 4:-177, 4:177, -3:178, -2:179;"
private const val kFarHemi: String = "0:-90, 60:90, -60:90;"

// A set of nested loops around the point -90:0 (lat:lng).
private const val kSouthPoint: String = "-89.9999:0.001"
private const val kSouth0a: String = "-90:0, -89.99:0.01, -89.99:0;"
private const val kSouth0b: String = "-90:0, -89.99:0.03, -89.99:0.02;"
private const val kSouth0c: String = "-90:0, -89.99:0.05, -89.99:0.04;"
private const val kSouth1: String = "-90:0, -89.9:0.1, -89.9:-0.1;"
private const val kSouth2: String = "-90:0, -89.8:0.2, -89.8:-0.2;"
private const val kSouthHemi: String = "0:-180, 0:60, 0:-60;"

// Two different loops that surround all the Near and Far loops except
// for the hemispheres.
private const val kNearFar1: String =
    "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, 1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;"
private const val kNearFar2: String = "-2:15, -2:170, -8:-175, 8:-175, 2:170, 2:15, 8:-4, -8:-4;"

// Loops that result from intersection of other loops.
private const val kFarHSouthH: String = "0:-180, 0:90, -60:90, 0:-90;"

// Rectangles that form a cross, with only shared vertices, no crossing edges.
// Optional holes outside the intersecting region.
private const val kCross1: String = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1;"
private const val kCross1SideHole: String = "-1.5:0.5, -1.2:0.5, -1.2:-0.5, -1.5:-0.5;"
private const val kCross2: String = "1:-2, 1:-1, 1:1, 1:2, -1:2, -1:1, -1:-1, -1:-2;"
private const val kCross2SideHole: String = "0.5:-1.5, 0.5:-1.2, -0.5:-1.2, -0.5:-1.5;"
private const val kCrossCenterHole: String = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;"

// Two rectangles that intersect, but no edges cross and there's always
// local containment (rather than crossing) at each shared vertex.
// In this ugly ASCII art, 1 is A+B, 2 is B+C:
//      +---+---+---+
//      | A | B | C |
//      +---+---+---+
private const val kOverlap1: String = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;"
private const val kOverlap1SideHole: String = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;"
private const val kOverlap2: String = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;"
private const val kOverlap2SideHole: String = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;"
private const val kOverlapCenterHole: String = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;"

// An empty polygon.
private const val kEmpty: String = ""

// By symmetry, the intersection of the two polygons has almost half the area
// of either polygon.
private const val kOverlap3: String = "-10:10, 0:10, 0:-10, -10:-10, -10:0"
private const val kOverlap4: String = "-10:0, 10:0, 10:-10, -10:-10"

// Some standard polygons to use in the tests.
private val empty: S2Polygon = S2Polygon()
private val full: S2Polygon = makePolygon("full")
private val near0: S2Polygon = makePolygon(kNear0)
private val near10: S2Polygon = makePolygon(kNear0 + kNear1)
private val near30: S2Polygon = makePolygon(kNear3 + kNear0)
private val near_32: S2Polygon = makePolygon(kNear2 + kNear3)
private val near_3210: S2Polygon = makePolygon(kNear0 + kNear2 + kNear3 + kNear1)
private val near_H3210: S2Polygon = makePolygon(kNear0 + kNear2 + kNear3 + kNearHemi + kNear1)

private val far_10: S2Polygon = makePolygon(kFar0 + kFar1)
private val far_21: S2Polygon = makePolygon(kFar2 + kFar1)
private val far_321: S2Polygon = makePolygon(kFar2 + kFar3 + kFar1)
private val far_H20: S2Polygon = makePolygon(kFar2 + kFarHemi + kFar0)
private val far_H3210: S2Polygon = makePolygon(kFar2 + kFarHemi + kFar0 + kFar1 + kFar3)

private val south_0ab: S2Polygon = makePolygon(kSouth0a + kSouth0b)
private val south_2: S2Polygon = makePolygon(kSouth2)
private val south_210b: S2Polygon = makePolygon(kSouth2 + kSouth0b + kSouth1)
private val south_H21: S2Polygon = makePolygon(kSouth2 + kSouthHemi + kSouth1)
private val south_H20abc: S2Polygon = makePolygon(kSouth2 + kSouth0b + kSouthHemi + kSouth0a + kSouth0c)

private val nf1N10F2S10abc: S2Polygon =
    makePolygon(kSouth0c + kFar2 + kNear1 + kNearFar1 + kNear0 + kSouth1 + kSouth0b + kSouth0a)

private val nf2N2F210S210ab: S2Polygon =
    makePolygon(kFar2 + kSouth0a + kFar1 + kSouth1 + kFar0 + kSouth0b + kNearFar2 + kSouth2 + kNear2)

private val f32_n0: S2Polygon = makePolygon(kFar2 + kNear0 + kFar3)
private val n32_s0b: S2Polygon = makePolygon(kNear3 + kSouth0b + kNear2)

private val cross1: S2Polygon = makePolygon(kCross1)
private val cross1_side_hole: S2Polygon = makePolygon(kCross1 + kCross1SideHole)
private val cross1_center_hole: S2Polygon = makePolygon(kCross1 + kCrossCenterHole)
private val cross2: S2Polygon = makePolygon(kCross2)
private val cross2_side_hole: S2Polygon = makePolygon(kCross2 + kCross2SideHole)
private val cross2_center_hole: S2Polygon = makePolygon(kCross2 + kCrossCenterHole)

private val overlap1: S2Polygon = makePolygon(kOverlap1)
private val overlap1_side_hole: S2Polygon = makePolygon(kOverlap1 + kOverlap1SideHole)
private val overlap1_center_hole: S2Polygon = makePolygon(kOverlap1 + kOverlapCenterHole)
private val overlap2: S2Polygon = makePolygon(kOverlap2)
private val overlap2_side_hole: S2Polygon = makePolygon(kOverlap2 + kOverlap2SideHole)
private val overlap2_center_hole: S2Polygon = makePolygon(kOverlap2 + kOverlapCenterHole)

private val far_H: S2Polygon = makePolygon(kFarHemi)
private val south_H: S2Polygon = makePolygon(kSouthHemi)
private val far_H_south_H: S2Polygon = makePolygon(kFarHSouthH)

private val far_south = run {
    val p = S2Polygon()
    p.initToIntersection(far_H, south_H)
    checkEqual(p, far_H_south_H, S1Angle.radians(1e-15))
    p
}

private val logger = KotlinLogging.logger(S2PolygonUnitTest::class.java.name)

fun makePolygon(str: String): S2Polygon {
    val polygon = S2TextParser.makeVerbatimPolygon(str)

    // Check that InitToSnapped() is idempotent.
    val snapped1 = S2Polygon()
    val snapped2 = S2Polygon()
    snapped1.initToSnapped(polygon)
    snapped2.initToSnapped(snapped1)
    assertThat(snapped1).isEqualTo(snapped2)

    return polygon
}

fun checkContains(aStr: String, bStr: String) {
    val a = makePolygon(aStr)
    val b = makePolygon(bStr)
    assertThat(a.contains(b)).isTrue()
    assertThat(a.approxContains(b, S1Angle.radians(1e-15))).isTrue()
    assertThat(a.approxDisjoint(b, S1Angle.radians(1e-15))).isFalse()
}

fun checkContainsPoint(aStr: String, bStr: String) {
    val a = S2TextParser.makePolygon(aStr)
    assertThat(a.contains(S2TextParser.makePoint(bStr)))
        .withFailMessage("$aStr did not contain $bStr")
        .isTrue()
}


fun checkEqual(a: S2Polygon, b: S2Polygon, max_error: S1Angle = S1Angle.zero()) {
    if (a.boundaryApproxEquals(b, max_error)) return
    val builder = S2Builder(S2Builder.Options())
    val a2 = S2Polygon()
    val b2 = S2Polygon()
    builder.startLayer(S2PolygonLayer(a2))
    builder.addPolygon(a)
    val error = S2Error()
    assertThat(builder.build(error)).isTrue()
    builder.startLayer(S2PolygonLayer(b2))
    builder.addPolygon(b)
    assertThat(builder.build(error)).isTrue()
    assertThat(a2.boundaryApproxEquals(b2, max_error))
        .withFailMessage(
            "\na: ${S2TextParser.toString(a)}" +
                    "\nb: ${S2TextParser.toString(b)}" +
                    "\na2: ${S2TextParser.toString(a2)}" +
                    "\nb2: ${S2TextParser.toString(b2)}"
        )
        .isTrue()
}

fun checkComplementary(a: S2Polygon, b: S2Polygon) {
    val b1 = S2Polygon()
    b1.initToComplement(b)
    checkEqual(a, b1)
}


// Given a pair of polygons where A contains B, check that various identities
// involving union, intersection, and difference operations hold true.
fun testOneNestedPair(a: S2Polygon, b: S2Polygon) {
    assertThat(a.contains(b)).isTrue()
    assertThat(a.intersects(b)).isEqualTo(!b.isEmpty())
    assertThat(b.intersects(a)).isEqualTo(!b.isEmpty())

    val c = S2Polygon()
    c.initToUnion(a, b)
    checkEqual(c, a)

    val d = S2Polygon()
    d.initToIntersection(a, b)
    checkEqual(d, b)

    val e = S2Polygon()
    e.initToDifference(b, a)
    assertThat(e.isEmpty()).isTrue()

    val f = S2Polygon()
    val g = S2Polygon()
    f.initToDifference(a, b)
    g.initToSymmetricDifference(a, b)
    checkEqual(f, g)
}


// Given a pair of disjoint polygons A and B, check that various identities
// involving union, intersection, and difference operations hold true.
fun testOneDisjointPair(a: S2Polygon, b: S2Polygon) {
    assertThat(a.intersects(b)).isFalse()
    assertThat(b.intersects(a)).isFalse()
    assertThat(a.contains(b)).isEqualTo(b.isEmpty())
    assertThat(b.contains(a)).isEqualTo(a.isEmpty())

    val ab = S2Polygon()
    val builder = S2Builder(S2Builder.Options())
    builder.startLayer(S2PolygonLayer(ab))
    builder.addPolygon(a)
    builder.addPolygon(b)
    val error = S2Error()
    assertThat(builder.build(error)).isTrue()

    val c = S2Polygon()
    c.initToUnion(a, b)
    checkEqual(c, ab)

    val d = S2Polygon()
    d.initToIntersection(a, b)
    assertThat(d.isEmpty()).isTrue()

    val e = S2Polygon()
    e.initToDifference(a, b)
    checkEqual(e, a)

    val f = S2Polygon()
    f.initToDifference(b, a)
    checkEqual(f, b)

    val g = S2Polygon()
    g.initToSymmetricDifference(a, b)
    checkEqual(g, ab)
}

// Given polygons A and B whose union covers the sphere, check that various
// identities involving union, intersection, and difference hold true.
fun testOneCoveringPair(a: S2Polygon, b: S2Polygon) {
    assertThat(a.contains(b)).isEqualTo(a.isFull())
    assertThat(b.contains(a)).isEqualTo(b.isFull())

    val c = S2Polygon()
    c.initToUnion(a, b)
    assertThat(c.isFull()).isTrue()
}

// Given polygons A and B such that both A and its complement intersect both B
// and its complement, check that various identities involving union,
// intersection, and difference hold true.
fun testOneOverlappingPair(a: S2Polygon, b: S2Polygon) {
    assertThat(a.contains(b)).isFalse()
    assertThat(b.contains(a)).isFalse()
    assertThat(a.intersects(b)).isTrue()

    val c = S2Polygon()
    c.initToUnion(a, b)
    assertThat(c.isFull()).isFalse()

    val d = S2Polygon()
    d.initToIntersection(a, b)
    assertThat(d.isEmpty()).isFalse()

    val e = S2Polygon()
    e.initToDifference(b, a)
    assertThat(e.isEmpty()).isFalse()

    val f = S2Polygon()
    val g = S2Polygon()
    val h = S2Polygon()
    f.initToDifference(a, b)
    g.initToUnion(e, f)
    h.initToSymmetricDifference(a, b)
    checkEqual(g, h)
}

// Given a pair of polygons where A contains B, test various identities
// involving A, B, and their complements.
fun testNestedPair(a: S2Polygon, b: S2Polygon) {
    logger.trace { "Test nested pair:\na = ${S2TextParser.toString(a)}\nb = ${S2TextParser.toString(b)}" }

    val a1 = S2Polygon()
    val b1 = S2Polygon()
    a1.initToComplement(a)
    b1.initToComplement(b)

    logger.trace {
        "Test nested pair:\na complement = ${S2TextParser.toString(a1)}\nb complement = ${
            S2TextParser.toString(
                b1
            )
        }"
    }

    testOneNestedPair(a, b)
    testOneNestedPair(b1, a1)
    testOneDisjointPair(a1, b)
    testOneCoveringPair(a, b1)
}

// Given a pair of disjoint polygons A and B, test various identities
// involving A, B, and their complements.
fun testDisjointPair(a: S2Polygon, b: S2Polygon) {
    val a1 = S2Polygon()
    val b1 = S2Polygon()
    a1.initToComplement(a)
    b1.initToComplement(b)

    testOneDisjointPair(a, b)
    testOneCoveringPair(a1, b1)
    testOneNestedPair(a1, b)
    testOneNestedPair(b1, a)
}

// Given polygons A and B such that both A and its complement intersect both B
// and its complement, test various identities involving these four polygons.
fun testOverlappingPair(a: S2Polygon, b: S2Polygon) {
    val a1 = S2Polygon()
    val b1 = S2Polygon()
    a1.initToComplement(a)
    b1.initToComplement(b)

    testOneOverlappingPair(a, b)
    testOneOverlappingPair(a1, b1)
    testOneOverlappingPair(a1, b)
    testOneOverlappingPair(a, b1)
}

// "a1" is the complement of "a", and "b1" is the complement of "b".
fun testOneComplementPair(a: S2Polygon, a1: S2Polygon, b: S2Polygon, b1: S2Polygon) {
    // Check DeMorgan's Law and that subtraction is the same as intersection
    // with the complement.  This function is called multiple times in order to
    // test the various combinations of complements.

    val aAndB1 = S2Polygon()
    aAndB1.initToIntersection(a, b1)
    val a1OrB = S2Polygon()
    a1OrB.initToUnion(a1, b)
    val aMinusB = S2Polygon()
    aMinusB.initToDifference(a, b)

    checkComplementary(a1OrB, aAndB1)
    checkEqual(aMinusB, aAndB1)
}

// Test identities that should hold for any pair of polygons A, B and their
// complements.
fun testComplements(a: S2Polygon, b: S2Polygon) {
    val a1 = S2Polygon()
    val b1 = S2Polygon()
    a1.initToComplement(a)
    b1.initToComplement(b)

    testOneComplementPair(a, a1, b, b1)
    testOneComplementPair(a1, a, b, b1)
    testOneComplementPair(a, a1, b1, b)
    testOneComplementPair(a1, a, b1, b)

    // There is a lot of redundancy if we do this test for each complementary
    // pair, so we just do it once instead.
    val aXorB1 = S2Polygon()
    val a1XorB = S2Polygon()
    aXorB1.initToSymmetricDifference(a, b1)
    a1XorB.initToSymmetricDifference(a1, b)
    checkEqual(aXorB1, a1XorB)
}

fun testDestructiveUnion(a: S2Polygon, b: S2Polygon) {
    val c = S2Polygon()
    c.initToUnion(a, b)
    val polygons = mutableListOf<S2Polygon>()
    polygons.add(a.clone())
    polygons.add(b.clone())
    val cDestructive = S2Polygon.destructiveUnion(polygons)
    checkEqual(c, cDestructive)
}

fun testRelation(
    a: S2Polygon,
    b: S2Polygon,
    contains: Boolean,
    contained: Boolean,
    intersects: Boolean,
    description: String = "Test relation (contains = $contains, contained = $contained, intersects = $intersects): \na = ${
        S2TextParser.toString(
            a
        )
    }\nb = ${S2TextParser.toString(b)}"
) {
    S2PolygonUnitTest.logger.debug { description }
    assertThat(a.contains(b)).isEqualTo(contains)
    assertThat(b.contains(a)).isEqualTo(contained)
    assertThat(a.intersects(b)).isEqualTo(intersects)
    if (contains) testNestedPair(a, b)
    if (contained) testNestedPair(b, a)
    if (!intersects) testDisjointPair(a, b)
    if (intersects && !(contains or contained)) {
        testOverlappingPair(a, b)  // See testOverlappingPair for definition
    }
    testDestructiveUnion(a, b)
    testComplements(a, b)
}

fun makeLoops(loop_vertices: List<List<S2Point>>): List<S2Loop> {
    val result = mutableListOf<S2Loop>()
    loop_vertices.forEach { vertices ->
        result.add(S2Loop(vertices))
        val error = S2Error()
        assertThat(result.last().findValidationError(error))
            .withFailMessage("Loop ${result.size - 1}")
            .isFalse()
    }
    return result
}


fun polylineIntersectionSharedEdgeTest(p: S2Polygon, startVertex: Int, direction: Int) {
    S2PolygonUnitTest.logger.debug { "Polyline intersection shared edge test start= $startVertex, direction=$direction" }
    val points = listOf(p.loop(0).vertex(startVertex), p.loop(0).vertex(startVertex + direction))
    val polyline = S2Polyline(points)
    var polylines: List<S2Polyline>
    if (direction < 0) {
        polylines = p.intersectWithPolyline(polyline)
        assertThat(polylines.size).isEqualTo(0)
        polylines = p.subtractFromPolyline(polyline)
        assertThat(polylines.size).isEqualTo(1)
        assertThat(polylines[0].numVertices).isEqualTo(2)
        assertThat(polylines[0].vertex(0)).isEqualTo(points[0])
        assertThat(polylines[0].vertex(1)).isEqualTo(points[1])
        assertThat(p.intersects(polyline)).isFalse()
        assertThat(p.contains(polyline)).isFalse()
    } else {
        polylines = p.intersectWithPolyline(polyline)
        assertThat(polylines.size).isEqualTo(1)
        assertThat(polylines[0].numVertices).isEqualTo(2)
        assertThat(polylines[0].vertex(0)).isEqualTo(points[0])
        assertThat(polylines[0].vertex(1)).isEqualTo(points[1])
        polylines = p.subtractFromPolyline(polyline)
        assertThat(polylines.size).isEqualTo(0)
        assertThat(p.intersects(polyline)).isTrue()
        assertThat(p.contains(polyline)).isTrue()
    }
}


fun checkCoveringIsConservative(polygon: S2Polygon, cells: List<S2CellId>) {
    // Check that Contains(S2Cell) and mayIntersect(S2Cell) are implemented
    // conservatively, by comparing against the Contains/Intersect result with
    // the "cell polygon" defined by the four cell vertices.  Please note that
    // the cell polygon is *not* an exact representation of the S2Cell: cell
    // vertices are rounded from their true mathematical positions, which leads
    // to tiny cracks and overlaps between the cell polygons at different cell
    // levels.  That is why Contains(S2Cell) and mayIntersect(S2Cell) cannot be
    // implemented by simply converting the cell to an S2Polygon.  But it is
    // still useful to do this as a sanity check.  In particular:
    //
    //  - If Contains(cell) is true, the polygon must contain the cell polygon.
    //  - If the polygon intersects the cell polygon, then mayIntersect(cell)
    //    must return true.
    //
    for (cell_id in cells) {
        val cell = S2Cell(cell_id)
        val cellPoly = S2Polygon(cell)
        if (polygon.contains(cell)) {
            assertThat(polygon.contains(cellPoly)).isTrue()
        }
        if (polygon.intersects(cellPoly)) {
            assertThat(polygon.mayIntersect(cell)).isTrue()
        }
    }
}

// Remove a random polygon from "pieces" and return it.
fun choosePiece(pieces: MutableList<S2Polygon>): S2Polygon {
    val i = S2Random.randomInt(pieces.size)
    return pieces.removeAt(i)
}

fun splitAndAssemble(polygon: S2Polygon) {
    // Normalize the polygon's loop structure by rebuilding it with S2Builder.
    val builder = S2Builder(S2Builder.Options())
    val expected = S2Polygon()
    builder.startLayer(S2PolygonLayer(expected))
    builder.addPolygon(polygon)
    val error = S2Error()
    assertThat(builder.build(error)).isTrue()

    val kIters = 3 // google.DEBUG_MODE ? 3 : 10
    repeat(kIters) {
        val coverer = S2RegionCoverer()
        // Compute the minimum level such that the polygon's bounding
        // cap is guaranteed to be cut.
        val diameter = 2 * polygon.capBound.radius().radians
        val minLevel = S2Coords.projection.kMaxWidth.getLevelForMaxValue(diameter)

        // Now choose a level that has up to 500 cells in the covering.
        val max = 4 // google.DEBUG_MODE ? 4 : 6
        val level = minLevel + S2Random.randomInt(max)
        coverer.setMinLevel(minLevel)
        coverer.setMaxLevel(level)
        coverer.setMaxCells(500)

        val cells = mutableListOf<S2CellId>()
        coverer.getCovering(polygon, cells)
        val covering = S2CellUnion(cells)
        checkCovering(polygon, covering, false)
        checkCoveringIsConservative(polygon, cells)
        logger.debug { "${cells.size} cells in covering" }
        val pieces = mutableListOf<S2Polygon>()
        for ((i, cell_id) in cells.withIndex()) {
            val cell = S2Cell(cell_id)
            val window = S2Polygon(cell)
            val piece = S2Polygon()
            piece.initToIntersection(polygon, window)
            logger.debug {
                "\nPiece $i:\n  Window: ${S2TextParser.toString(window)}" +
                        "\n  Piece: ${S2TextParser.toString(piece)}"
            }
            pieces.add(piece)
        }

        // Now we repeatedly remove two random pieces, compute their union, and
        // insert the result as a new piece until only one piece is left.
        //
        // We don't use S2Polygon.DestructiveUnion() because it joins the pieces
        // in a mostly deterministic order.  We don't just call random_shuffle()
        // on the pieces and repeatedly join the last two pieces in the vector
        // because this always joins a single original piece to the current union
        // rather than doing the unions according to a random tree structure.
        while (pieces.size > 1) {
            val a = choosePiece(pieces)
            val b = choosePiece(pieces)
            val c = S2Polygon()
            c.initToUnion(a, b)
            S2PolygonUnitTest.logger.debug {
                "\nJoining piece a: ${S2TextParser.toString(a)}" +
                        "\n  With piece b: ${S2TextParser.toString(b)}" +
                        "\n  To get piece c: ${S2TextParser.toString(c)}"
            }
            pieces.add(c)
        }
        val result = pieces[0]
        pieces.removeLast()

        // The moment of truth!
        assertThat(expected.boundaryNear(result, S1Angle.radians(2e-15)))
            .withFailMessage(
                "\nActual:\n${S2TextParser.toString(result)}" +
                        "\nExpected:\n${S2TextParser.toString(expected)}"
            )
            .isTrue()

        // Check that ApproxEquals produces the same result.
        if (!expected.approxEquals(result, kIntersectionMergeRadius)) {
            val symmetricDifference = S2Polygon()
            symmetricDifference.initToApproxSymmetricDifference(expected, result, kIntersectionMergeRadius)
            fail<Void>(S2TextParser.toString(symmetricDifference))
        }
    }
}

// Helper function for testing the distance methods.  "boundary_x" is the
// expected result of projecting "x" onto the polygon boundary.  For
// convenience it can be set to S2Point() to indicate that (boundary_x == x).
fun testDistanceMethods(polygon: S2Polygon, x: S2Point, boundary_x: S2Point) {
    // This error is not guaranteed by the implementation but is okay for tests.
    val kMaxError = S1Angle.radians(1e-15)

    if (boundary_x == S2Point()) boundary_x.plusAssign(x)
    assertThat(S1Angle(boundary_x, polygon.projectToBoundary(x))).isLessThanOrEqualTo(kMaxError)

    if (polygon.isEmpty() || polygon.isFull()) {
        assertThat(polygon.getDistanceToBoundary(x)).isEqualTo(S1Angle.infinity())
    } else {
        // EXPECT_NEAR only works with doubles.
        assertThat(polygon.getDistanceToBoundary(x).degrees()).isCloseTo(
            S1Angle(x, boundary_x).degrees(),
            Offset.offset(kMaxError.degrees())
        )
    }
    if (polygon.contains(x)) {
        assertThat(polygon.getDistance(x)).isEqualTo(S1Angle.zero())
        assertThat(polygon.project(x)).isEqualTo(x)
    } else {
        assertThat(polygon.getDistance(x)).isEqualTo(polygon.getDistanceToBoundary(x))
        assertThat(polygon.project(x)).isEqualTo(polygon.projectToBoundary(x))
    }
}

class S2PolygonUnitTest {

    @Test
    fun init() {
        checkContains(kNear1, kNear0)
        checkContains(kNear2, kNear1)
        checkContains(kNear3, kNear2)
        checkContains(kNearHemi, kNear3)
        checkContains(kFar1, kFar0)
        checkContains(kFar2, kFar1)
        checkContains(kFar3, kFar2)
        checkContains(kFarHemi, kFar3)
        checkContains(kSouth1, kSouth0a)
        checkContains(kSouth1, kSouth0b)
        checkContains(kSouth1, kSouth0c)
        checkContains(kSouthHemi, kSouth2)
        checkContains(kNearFar1, kNear3)
        checkContains(kNearFar1, kFar3)
        checkContains(kNearFar2, kNear3)
        checkContains(kNearFar2, kFar3)

        checkContainsPoint(kNear0, kNearPoint)
        checkContainsPoint(kNear1, kNearPoint)
        checkContainsPoint(kNear2, kNearPoint)
        checkContainsPoint(kNear3, kNearPoint)
        checkContainsPoint(kNearHemi, kNearPoint)
        checkContainsPoint(kSouth0a, kSouthPoint)
        checkContainsPoint(kSouth1, kSouthPoint)
        checkContainsPoint(kSouth2, kSouthPoint)
        checkContainsPoint(kSouthHemi, kSouthPoint)
    }

    @Test
    fun overlapFractions() {
        var a = makePolygon(kEmpty)
        var b = makePolygon(kEmpty)
        var result = S2Polygon.getOverlapFractions(a, b)
        assertThat(result.first).isEqualTo(1.0)
        assertThat(result.second).isEqualTo(1.0)

        b = makePolygon(kOverlap3)
        result = S2Polygon.getOverlapFractions(a, b)
        assertThat(result.first).isEqualTo(1.0)
        assertThat(result.second).isEqualTo(0.0)

        a = makePolygon(kOverlap4)
        result = S2Polygon.getOverlapFractions(a, b)
        assertThat(result.first).isCloseTo(0.5, Offset.offset(1e-14))
        assertThat(result.second).isCloseTo(0.5, Offset.offset(1e-14))
    }

    @Test
    fun originNearPole() {
        // S2Polygon operations are more efficient if S2PointUtil.origin() is near a pole.
        // (Loops that contain a pole tend to have very loose bounding boxes because
        // they span the full longitude range.  S2Polygon canonicalizes all loops so
        // that they don't contain S2PointUtil.origin(), thus by placing S2PointUtil.origin() near a
        // pole we minimize the number of canonical loops which contain that pole.)
        assertThat(S2LatLng.latitude(S2PointUtil.origin()).degrees()).isGreaterThanOrEqualTo(80.0)
    }

    @Test
    fun approxContainsAndDisjoint() {
        // We repeatedly choose a random cell id and intersect its bounding polygon
        // "A" with the bounding polygon "B" of one its child cells.  The result may
        // not be contained by either A or B, because the vertices of B near the
        // edge midpoints of A may be slightly outside A, and even when the crossing
        // edges are intersected, the intersection point may also be slightly
        // outside A and/or B.
        //
        // We repeat the test many times and expect that some fraction of the exact
        // tests should fail, while all of the approximate test should succeed.
        val kIters = 1000
        var exactContains = 0
        var exactDisjoint = 0
        S2Random.reset(1)
        repeat(kIters) { iter ->
            S2Random.reset(2 * iter)
            val id = S2Random.randomCellId(10)
            val parentPolygon = S2Polygon((S2Cell(id)))
            val childPolygon = S2Polygon(S2Cell(id.child(0)))

            // Get the intersection.  There is no guarantee that the intersection will
            // be contained by A or B.  Similarly, the intersection may slightly
            // overlap an adjacent disjoint polygon C.
            val intersection = S2Polygon()
            intersection.initToIntersection(parentPolygon, childPolygon)
            if (parentPolygon.contains(intersection)) {
                ++exactContains
            }
            assertThat(parentPolygon.approxContains(intersection, kIntersectionMergeRadius)).isTrue()

            val adjacentPolygon = S2Polygon(S2Cell(id.child(1)))
            if (!adjacentPolygon.intersects(intersection)) {
                ++exactDisjoint
            }
            assertThat(adjacentPolygon.approxDisjoint(intersection, kIntersectionMergeRadius)).isTrue()
        }
        // All of the approximate results are true, so we check that at least some
        // of the exact results are false in order to make sure that this test
        // actually tests something.
        //
        // There are two vertices in each child cell that have a 50% chance of being
        // outside the parent cell.  When a vertex is outside, an intersection point
        // is computed near that vertex that also has a 50% chance of being
        // outside.  Snapping used to choose one of these vertices at random, but
        // currently the vertex whose S2CellId is smaller is always chosen.  For the
        // exact containment test, it turns out that one vertex is adjacent to a
        // lower-numbered S2CellId and the other is adjacent to a higher-numbered
        // S2CellId, which means that one vertex will always be chosen outside the
        // parent if possible, and the other will always be chosen inside if
        // possible.  This works out to an expectation that 0.5 * 0.75 = 37.5% of
        // the exact containment tests will succeed.
        //
        // For the exact disjoint test, there is one shared vertex that might be
        // replaced by a computed intersection point.  The shared vertex is inside
        // the parent 50% of the time.  Otherwise there is a 50% chance that the
        // intersection point will not be chosen for snapping because it has a
        // higher S2CellId that the shared vertex, and otherwise there is still a
        // 50% chance that intersection point is on the side of the shared edge that
        // results in no intersection.  This works out to an expectation that
        // (1 - 0.5 * 0.5 * 0.5) = 87.5% of the exact disjoint tests will succeed.
        logger.debug { "Exact contains: $exactContains" }
        logger.debug { "Exact disjoint: $exactDisjoint" }
        assertThat(exactContains.toDouble()).isLessThan(0.40 * kIters)  // about 37.5% succeed
        assertThat(exactDisjoint.toDouble()).isLessThan(0.98 * kIters)  // TODO(fmeurisse) check why it is not about 87.5% succeed
    }

    @Test
    fun relations() {
        testRelation(near10, empty, contains = true, contained = false, intersects = false)
        testRelation(near10, near10, contains = true, contained = true, intersects = true)
        testRelation(full, near10, contains = true, contained = false, intersects = true)
        testRelation(near10, near30, contains = false, contained = true, intersects = true)
        testRelation(near10, near_32, contains = false, contained = false, intersects = false)
        testRelation(near10, near_3210, contains = false, contained = true, intersects = true)
        testRelation(near10, near_H3210, contains = false, contained = false, intersects = false)
        testRelation(near30, near_32, contains = true, contained = false, intersects = true)
        testRelation(near30, near_3210, contains = true, contained = false, intersects = true)
        testRelation(near30, near_H3210, contains = false, contained = false, intersects = true)
        testRelation(near_32, near_3210, contains = false, contained = true, intersects = true)
        testRelation(near_32, near_H3210, contains = false, contained = false, intersects = false)
        testRelation(near_3210, near_H3210, contains = false, contained = false, intersects = false)

        testRelation(far_10, far_21, contains = false, contained = false, intersects = false)
        testRelation(far_10, far_321, contains = false, contained = true, intersects = true)
        testRelation(far_10, far_H20, contains = false, contained = false, intersects = false)
        testRelation(far_10, far_H3210, contains = false, contained = false, intersects = false)
        testRelation(far_21, far_321, contains = false, contained = false, intersects = false)
        testRelation(far_21, far_H20, contains = false, contained = false, intersects = false)
        testRelation(far_21, far_H3210, contains = false, contained = true, intersects = true)
        testRelation(far_321, far_H20, contains = false, contained = false, intersects = true)
        testRelation(far_321, far_H3210, contains = false, contained = false, intersects = true)
        testRelation(far_H20, far_H3210, contains = false, contained = false, intersects = true)

        testRelation(south_0ab, south_2, contains = false, contained = true, intersects = true)
        testRelation(south_0ab, south_210b, contains = false, contained = false, intersects = true)
        testRelation(south_0ab, south_H21, contains = false, contained = true, intersects = true)
        testRelation(south_0ab, south_H20abc, contains = false, contained = true, intersects = true)
        testRelation(south_2, south_210b, contains = true, contained = false, intersects = true)
        testRelation(south_2, south_H21, contains = false, contained = false, intersects = true)
        testRelation(south_2, south_H20abc, contains = false, contained = false, intersects = true)
        testRelation(south_210b, south_H21, contains = false, contained = false, intersects = true)
        testRelation(south_210b, south_H20abc, contains = false, contained = false, intersects = true)
        testRelation(south_H21, south_H20abc, contains = true, contained = false, intersects = true)

        testRelation(nf1N10F2S10abc, nf2N2F210S210ab, contains = false, contained = false, intersects = true)
        testRelation(nf1N10F2S10abc, near_32, contains = true, contained = false, intersects = true)
        testRelation(nf1N10F2S10abc, far_21, contains = false, contained = false, intersects = false)
        testRelation(nf1N10F2S10abc, south_0ab, contains = false, contained = false, intersects = false)
        testRelation(nf1N10F2S10abc, f32_n0, contains = true, contained = false, intersects = true)

        testRelation(nf2N2F210S210ab, near10, contains = false, contained = false, intersects = false)
        testRelation(nf2N2F210S210ab, far_10, contains = true, contained = false, intersects = true)
        testRelation(nf2N2F210S210ab, south_210b, contains = true, contained = false, intersects = true)
        testRelation(nf2N2F210S210ab, south_0ab, contains = true, contained = false, intersects = true)
        testRelation(nf2N2F210S210ab, n32_s0b, contains = true, contained = false, intersects = true)

        testRelation(cross1, cross2, contains = false, contained = false, intersects = true)
        testRelation(cross1_side_hole, cross2, contains = false, contained = false, intersects = true)
        testRelation(cross1_center_hole, cross2, contains = false, contained = false, intersects = true)
        testRelation(cross1, cross2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(cross1, cross2_center_hole, contains = false, contained = false, intersects = true)
        testRelation(cross1_side_hole, cross2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(cross1_center_hole, cross2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(cross1_side_hole, cross2_center_hole, contains = false, contained = false, intersects = true)
        testRelation(cross1_center_hole, cross2_center_hole, contains = false, contained = false, intersects = true)

        // These cases, when either polygon has a hole, test a different code path
        // from the other cases.
        testRelation(overlap1, overlap2, contains = false, contained = false, intersects = true)
        testRelation(overlap1_side_hole, overlap2, contains = false, contained = false, intersects = true)
        testRelation(overlap1_center_hole, overlap2, contains = false, contained = false, intersects = true)
        testRelation(overlap1, overlap2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(overlap1, overlap2_center_hole, contains = false, contained = false, intersects = true)
        testRelation(overlap1_side_hole, overlap2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(overlap1_center_hole, overlap2_side_hole, contains = false, contained = false, intersects = true)
        testRelation(overlap1_side_hole, overlap2_center_hole, contains = false, contained = false, intersects = true)
        testRelation(overlap1_center_hole, overlap2_center_hole, contains = false, contained = false, intersects = true)
    }

    @Test
    fun emptyAndFull() {
        assertThat(empty.isEmpty()).isTrue()
        assertThat(full.isEmpty()).isFalse()
        assertThat(empty.isFull()).isFalse()
        assertThat(full.isFull()).isTrue()

        testNestedPair(empty, empty)
        testNestedPair(full, empty)
        testNestedPair(full, full)
    }


    @ParameterizedTest
    @MethodSource("testCases")
    fun operations(test: TestCase) {
        logger.debug { "Polygon operation test case $test" }
        val a = makePolygon(test.a)
        val b = makePolygon(test.b)
        val expectedAAndB = makePolygon(test.aAndB)
        val expectedAOrB = makePolygon(test.aOrB)
        val expectedAMinusB = makePolygon(test.aMinusB)
        val expectedAXorB = makePolygon(test.aXorB)

        // The intersections in the "expected" data were computed in lat-lng
        // space, while the actual intersections are computed using geodesics.
        // The error due to this depends on the length and direction of the line
        // segment being intersected, and how close the intersection is to the
        // endpoints of the segment.  The worst case is for a line segment between
        // two points at the same latitude, where the intersection point is in the
        // middle of the segment.  In this case the error is approximately
        // (p * t^2) / 8, where "p" is the absolute latitude in radians, "t" is
        // the longitude difference in radians, and both "p" and "t" are small.
        // The test cases all have small latitude and longitude differences.
        // If "p" and "t" are converted to degrees, the following error bound is
        // valid as long as (p * t^2 < 150).

        val kMaxError = S1Angle.radians(1e-4)

        val aAndB = S2Polygon()
        aAndB.initToIntersection(a, b)
        checkEqual(aAndB, expectedAAndB, kMaxError)

        val aOrB = S2Polygon()
        aOrB.initToUnion(a, b)
        checkEqual(aOrB, expectedAOrB, kMaxError)
        testDestructiveUnion(a, b)

        val aMinusB = S2Polygon()
        aMinusB.initToDifference(a, b)
        checkEqual(aMinusB, expectedAMinusB, kMaxError)

        val aXorB = S2Polygon()
        aXorB.initToSymmetricDifference(a, b)
        checkEqual(aXorB, expectedAXorB, kMaxError)
    }

    @Test
    fun intersectionSnapFunction() {
        // This tests that an intersection point is rounded to the nearest allowable
        // vertex position (using E0 coordinates, i.e. integer lat/lng values).
        val a = makePolygon("0:0, 0:10, 1:10, 1:0")
        val b = makePolygon("0:0, 0:10, 3:0")
        val expected = makePolygon("0:0, 0:10, 1:7, 1:0")
        val actual = S2Polygon()
        actual.initToIntersection(a, b, IntLatLngSnapFunction(0))  // E0 coords
        checkEqual(expected, actual)
    }

    @Test
    fun intersectionPreservesLoopOrder() {
        val a = makePolygon("0:0, 0:10, 10:10, 10:0")
        val b = makePolygon("1:1, 1:9, 9:5; 2:2, 2:8, 8:5")
        val actual = S2Polygon()
        actual.initToIntersection(a, b)
        assertThat(S2TextParser.toString(actual)).isEqualTo(S2TextParser.toString(b))
    }


    // The "Bug" tests are regression tests from previous versions of the algorithm.
    @Test
    fun bug1() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10531193335759943, -0.80522214810955617, 0.58354664670985534),
                S2Point(-0.10531194840431297, -0.80522215192439039, 0.58354663873039425),
                S2Point(-0.10531192794033867, -0.80522217497559767, 0.58354661061568747),
                S2Point(-0.10531191284235047, -0.80522217121852058, 0.58354661852470402)
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10531174240075937, -0.80522236320875284, 0.58354638436119843),
                S2Point(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
                S2Point(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),  // A
                S2Point(-0.10531191288915481, -0.80522217116640804, 0.5835466185881667),   // B
                S2Point(-0.10531191288915592, -0.8052221711664066, 0.58354661858816803),   // B
                S2Point(-0.10531192039151964, -0.80522217309710431, 0.58354661457010204),  // A
                S2Point(-0.10531192794033779, -0.80522217497559878, 0.58354661061568636),
                S2Point(-0.1053117575499668, -0.80522236690813498, 0.58354637652254981),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    @Test
    fun bug2() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10618951389689163, -0.80546461394606728, 0.58305277875939732),
                S2Point(-0.10618904764039243, -0.8054645437464607, 0.58305296065497536),
                S2Point(-0.10618862643748632, -0.80546451917975415, 0.58305307130470341),
                S2Point(-0.10617606798507535, -0.80544758470051458, 0.58307875187433833),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10618668131028208, -0.80544613076731553, 0.58307882755616247),
                S2Point(-0.10618910658843225, -0.80546454998744921, 0.58305294129732887),
                S2Point(-0.10618904764039225, -0.80546454374646081, 0.58305296065497536),
                S2Point(-0.10618898834264634, -0.80546453817003949, 0.58305297915823251),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            ).isFalse()
    }

    @Test
    fun bug3() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10703494861068318, -0.80542232562508131, 0.58295659972299307),
                S2Point(-0.10703494998722708, -0.80542232255642865, 0.58295660370995028),
                S2Point(-0.10703495367938694, -0.80542232008675829, 0.58295660644418046),
                S2Point(-0.10703495869785147, -0.80542231887781635, 0.58295660719304865),
                S2Point(-0.10703496369792719, -0.80542231925353791, 0.58295660575589636),
                S2Point(-0.10703496733984781, -0.80542232111324863, 0.58295660251780734),
                S2Point(-0.10703496864776367, -0.80542232395864055, 0.58295659834642488),
                S2Point(-0.10703496727121976, -0.80542232702729322, 0.58295659435946767),
                S2Point(-0.10703496357905991, -0.80542232949696357, 0.5829565916252375),
                S2Point(-0.10703495856059538, -0.80542233070590552, 0.58295659087636931),
                S2Point(-0.10703495356051966, -0.80542233033018396, 0.58295659231352159),
                S2Point(-0.10703494991859903, -0.80542232847047324, 0.58295659555161061),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10703494861068762, -0.80542232562508098, 0.58295659972299274),
                S2Point(-0.10703494998723152, -0.80542232255642832, 0.58295660370994995),
                S2Point(-0.10703495367939138, -0.80542232008675796, 0.58295660644418013),
                S2Point(-0.10703495869785591, -0.80542231887781601, 0.58295660719304832),
                S2Point(-0.10703496369793163, -0.80542231925353758, 0.58295660575589603),
                S2Point(-0.10703496733985225, -0.8054223211132483, 0.58295660251780701),
                S2Point(-0.10703496864776811, -0.80542232395864022, 0.58295659834642455),
                S2Point(-0.1070349672712242, -0.80542232702729288, 0.58295659435946734),
                S2Point(-0.10703496357906438, -0.80542232949696346, 0.58295659162523727),
                S2Point(-0.10703495856059982, -0.80542233070590519, 0.58295659087636897),
                S2Point(-0.1070349535605241, -0.80542233033018362, 0.58295659231352126),
                S2Point(-0.10703494991860348, -0.8054223284704729, 0.58295659555161028),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty()).withFailMessage(
            "\nS2Polygon: ${S2TextParser.toString(a)}" +
                    "\nS2Polygon: ${S2TextParser.toString(b)}"
        ).isFalse()
    }

    @Test
    fun bug4() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10667065556339718, -0.80657502337947207, 0.58142764201754193),
                S2Point(-0.10667064691895933, -0.80657502457251051, 0.58142764194845853),
                S2Point(-0.10667064691930939, -0.80657502457246333, 0.58142764194845975),
                S2Point(-0.10667065556339746, -0.80657502337947395, 0.5814276420175396),
                S2Point(-0.10667077559567185, -0.80657589269604968, 0.58142641405029793),
                S2Point(-0.10667077059539463, -0.80657589232162286, 0.58142641548708696),
                S2Point(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
                S2Point(-0.10667063169531328, -0.80657498170361974, 0.58142770421053058),
                S2Point(-0.10667064898418178, -0.8065749793175444, 0.58142770434869739),
            ),
            listOf(
                S2Point(-0.10667064691897719, -0.80657502457250896, 0.58142764194845697),
                S2Point(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
                S2Point(-0.10667064691861985, -0.80657502457255736, 0.58142764194845586),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10667064691896312, -0.80657502457251107, 0.58142764194845697),
                S2Point(-0.10667064691896297, -0.80657502457251007, 0.58142764194845853),
                S2Point(-0.10667064033974753, -0.80657498051058207, 0.58142770427961399),
                S2Point(-0.10667064076268165, -0.80657498045444342, 0.58142770427989865),
                S2Point(-0.10667051785242875, -0.80657409963649807, 0.58142894872603923),
                S2Point(-0.1066707756642685, -0.80657588679775971, 0.58142642222003538),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Loop 1: Edge 1 crosses edge 3
        assertThat(c.isEmpty()).withFailMessage(
            "\nS2Polygon: ${S2TextParser.toString(a)}" +
                    "\nS2Polygon: ${S2TextParser.toString(b)}"
        ).isFalse()
    }

    @Test
    fun bug5() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10574444273627338, -0.80816264611829447, 0.57938868667714882),
                S2Point(-0.10574444845633162, -0.80816268110163325, 0.57938863683652475),
                S2Point(-0.10574444825833453, -0.80816268112970524, 0.57938863683350494),
                S2Point(-0.10574444253827629, -0.80816264614636646, 0.57938868667412902),
                S2Point(-0.10574408792844124, -0.80816047738475361, 0.57939177648757634),
                S2Point(-0.10574408812643833, -0.80816047735668162, 0.57939177649059592),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.1057440881264381, -0.80816047735668017, 0.57939177649059825),
                S2Point(-0.10574408802743954, -0.80816047737071606, 0.57939177648908835),
                S2Point(-0.10574408812649677, -0.8081604773570521, 0.57939177649006868),
                S2Point(-0.10574408812649701, -0.80816047735705354, 0.57939177649006646),
                S2Point(-0.10574408802703171, -0.80816047737077379, 0.57939177648908202),
                S2Point(-0.10574408792844098, -0.80816047738475194, 0.57939177648757834),
                S2Point(-0.10574408792838257, -0.80816047738438168, 0.5793917764881058),
                S2Point(-0.1057440879283823, -0.80816047738438002, 0.57939177648810791),
                S2Point(-0.10574407993470979, -0.80816042849578984, 0.57939184613891748),
                S2Point(-0.10574408013270691, -0.80816042846771807, 0.57939184614193739),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Loop 0 edge 8 crosses loop 1 edge 0
        assertThat(c.isEmpty()).withFailMessage(
            "\nS2Polygon: ${S2TextParser.toString(a)}" +
                    "\nS2Polygon: ${S2TextParser.toString(b)}"
        ).isFalse()
    }

    @Test
    fun bug6() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10618849949725141, -0.80552159562437586, 0.58297423747304822),
                S2Point(-0.10618849959636036, -0.80552159561106063, 0.58297423747339361),
                S2Point(-0.10618849949722192, -0.80552159562415893, 0.5829742374733532),
                S2Point(-0.10618834540082922, -0.80552043435619214, 0.58297587011440333),
                S2Point(-0.10618834559910612, -0.80552043432999554, 0.58297587011448437),
                S2Point(-0.10618849969546933, -0.80552159559774539, 0.58297423747373922),
                S2Point(-0.10618849969546955, -0.80552159559774716, 0.582974237473737),
                S2Point(-0.10618849969549882, -0.80552159559796233, 0.58297423747343424),
                S2Point(-0.10618849959710704, -0.80552159561096182, 0.58297423747339394),
                S2Point(-0.10618849949725161, -0.80552159562437742, 0.58297423747304589),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10618856154870562, -0.80552206324314812, 0.58297358004005528),
                S2Point(-0.10618849949722212, -0.80552159562416048, 0.58297423747335086),
                S2Point(-0.10618849969549901, -0.80552159559796388, 0.58297423747343191),
                S2Point(-0.10618856174698249, -0.8055220632169513, 0.58297358004013622),
                S2Point(-0.10618857104277038, -0.80552213326985989, 0.58297348155149287),
                S2Point(-0.10618857084449349, -0.80552213329605649, 0.58297348155141182),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Loop 0 edge 0 crosses loop 1 edge 4
        assertThat(c.isEmpty()).withFailMessage(
            "\nS2Polygon: ${S2TextParser.toString(a)}" +
                    "\nS2Polygon: ${S2TextParser.toString(b)}"
        ).isFalse()
    }

    @Test
    fun bug7() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10651728339354898, -0.80806023027835039, 0.57938996589599123),
                S2Point(-0.10651728368541774, -0.80806023024121265, 0.57938996589412783),
                S2Point(-0.10651743884289547, -0.80806147782022508, 0.5793881973990701),
                S2Point(-0.1065172793067945, -0.80806153133252501, 0.5793881520963412),
                S2Point(-0.10651707335497011, -0.80806158532388361, 0.57938811465868356),
                S2Point(-0.10651593657771009, -0.80806167503227055, 0.57938819853274059),
                S2Point(-0.10651567693742285, -0.80806182530835402, 0.57938803667826444),
                S2Point(-0.10651496089498214, -0.80806213485510237, 0.57938773659696563),
                S2Point(-0.10651453461919227, -0.80806229235522298, 0.57938759530083062),
                S2Point(-0.10651448583749658, -0.80806230280784852, 0.57938758969074455),
                S2Point(-0.10651428153471061, -0.80806061225022852, 0.57938998503506256),
                S2Point(-0.10651428161845182, -0.8080606122395747, 0.57938998503452654),
                S2Point(-0.10651427761078044, -0.80806057978063328, 0.57939003104095654),
                S2Point(-0.10651427761077951, -0.80806057978062562, 0.57939003104096709),
                S2Point(-0.10651387099203104, -0.8080572864940091, 0.5793946988282096),
                S2Point(-0.10651387099202798, -0.80805728649398445, 0.57939469882824468),
                S2Point(-0.10651386444607201, -0.80805723347699177, 0.57939477397218053),
                S2Point(-0.10651386444607169, -0.8080572334769891, 0.57939477397218409),
                S2Point(-0.106513765993723, -0.80805643609199118, 0.57939590414857456),
                S2Point(-0.10651376671438624, -0.8080564359989727, 0.57939590414581921),
                S2Point(-0.10651368187839319, -0.80805575808078389, 0.57939686520139033),
                S2Point(-0.10651465698432123, -0.80805552598235797, 0.57939700963750851),
                S2Point(-0.1065149024434091, -0.80805548225095913, 0.57939702550292815),
                S2Point(-0.10651504788182964, -0.80805555533715756, 0.5793968968362615),
                S2Point(-0.10651511658091152, -0.80805559604710031, 0.57939682743066534),
                S2Point(-0.10651517919248171, -0.80805562751022852, 0.57939677204023521),
                S2Point(-0.10651528575974038, -0.80805561374213786, 0.57939677165077275),
                S2Point(-0.10651648823358072, -0.80805539171529139, 0.57939686023850034),
                S2Point(-0.10651666406737116, -0.80805537863686483, 0.57939684615295572),
                S2Point(-0.10651674780673852, -0.80805605121551227, 0.57939589274577097),
                S2Point(-0.10651674667750256, -0.80805605136137271, 0.57939589274994641),
                S2Point(-0.10651678418140036, -0.80805634336988752, 0.57939547860450136),
                S2Point(-0.10651680240261223, -0.80805648524178364, 0.57939527739240138),
                S2Point(-0.10651680240261237, -0.80805648524178486, 0.57939527739239993),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10651727337444802, -0.80806023111043901, 0.57938996657744879),
                S2Point(-0.10651727440799089, -0.80806022882029649, 0.57938996958144073),
                S2Point(-0.10651679374955145, -0.80805648637258243, 0.57939527740611751),
                S2Point(-0.10651677552833975, -0.80805634450068775, 0.57939547861821594),
                S2Point(-0.10651673802444192, -0.80805605249217261, 0.57939589276366099),
                S2Point(-0.10651674651102909, -0.80805605138312775, 0.5793958927502102),
                S2Point(-0.10651673915225639, -0.80805605233507238, 0.57939589277542292),
                S2Point(-0.10651665541288889, -0.80805537975642383, 0.57939684618260878),
                S2Point(-0.10651667272185343, -0.80805537751730583, 0.57939684612330267),
                S2Point(-0.1065167564612207, -0.8080560500959526, 0.57939589271611924),
                S2Point(-0.1065167553320342, -0.80805605024202609, 0.57939589271998793),
                S2Point(-0.10651679283446101, -0.80805634223908773, 0.57939547859078699),
                S2Point(-0.10651681105567287, -0.80805648411098374, 0.57939527737868723),
                S2Point(-0.10651680240318392, -0.80805648524170914, 0.5793952773924006),
                S2Point(-0.10651680240261234, -0.80805648524178475, 0.57939527739239982),
                S2Point(-0.1065168110556733, -0.80805648411098718, 0.57939527737868224),
                S2Point(-0.10651729169518892, -0.80806022641135866, 0.57938996976297907),
                S2Point(-0.10651729210462238, -0.80806022661896348, 0.579389969398166),
                S2Point(-0.1065172934126499, -0.80806022944626155, 0.57938996521453356),
                S2Point(-0.10651729203606744, -0.80806023249651726, 0.57938996121349717),
                S2Point(-0.1065172883437291, -0.80806023495241674, 0.57938995846713126),
                S2Point(-0.10651728332499401, -0.80806023615590394, 0.5793899577113224),
                S2Point(-0.10651727832462815, -0.80806023578450537, 0.57938995914858893),
                S2Point(-0.10651727468247554, -0.80806023393773707, 0.57938996239381635),
            ),
            listOf(
                S2Point(-0.10651680240204828, -0.80805648524185858, 0.57939527739240082),
                S2Point(-0.10651679861449742, -0.80805648573682254, 0.57939527739840524),
                S2Point(-0.10651680240261419, -0.80805648524178353, 0.57939527739240138),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Loop 0: Edge 33 crosses edge 35
        assertThat(c.isEmpty()).withFailMessage(
            "\nS2Polygon: ${S2TextParser.toString(a)}" +
                    "\nS2Polygon: ${S2TextParser.toString(b)}"
        ).isFalse()
    }

    @Test
    fun bug8() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10703872198218529, -0.80846112144645677, 0.57873424566545062),
                S2Point(-0.10703872122182066, -0.80846111957630917, 0.57873424841857957),
                S2Point(-0.10703873813385757, -0.80846111582010538, 0.57873425053786276),
                S2Point(-0.1070387388942222, -0.80846111769025297, 0.57873424778473381),
                S2Point(-0.10703873050793056, -0.80846111955286837, 0.57873424673382978),
                S2Point(-0.1070387388942227, -0.80846111769025419, 0.57873424778473193),
                S2Point(-0.10703919382477994, -0.80846223660916783, 0.57873260056976505),
                S2Point(-0.10703917691274406, -0.80846224036537406, 0.57873259845047831),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10703917691274355, -0.80846224036537273, 0.57873259845047997),
                S2Point(-0.1070391853685064, -0.8084622384873289, 0.57873259951008804),
                S2Point(-0.10703919381027188, -0.80846223657409677, 0.57873260062144094),
                S2Point(-0.10703919381027233, -0.80846223657409788, 0.57873260062143939),
                S2Point(-0.10703918536876245, -0.80846223848727206, 0.57873259951012024),
                S2Point(-0.10703919382478132, -0.80846223660917116, 0.57873260056976017),
                S2Point(-0.10703957146434441, -0.80846316542623331, 0.57873123320737097),
                S2Point(-0.10703955455230836, -0.8084631691824391, 0.57873123108808489),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(a)}" }
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(b)}" }
        val c = S2Polygon()
        c.initToUnion(a, b)
        //  Loop 1: Edge 1 crosses edge 3
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(c)}" }
    }

    @Test
    fun bug9() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10639937100501309, -0.80810205676564995, 0.57935329437301375),
                S2Point(-0.10639937101137514, -0.80810205688156922, 0.57935329421015713),
                S2Point(-0.10639937101137305, -0.80810205688156944, 0.57935329421015713),
                S2Point(-0.106399371005011, -0.80810205676565017, 0.57935329437301375),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10639937099530022, -0.8081020567669569, 0.57935329437297489),
                S2Point(-0.10639937102108385, -0.80810205688026293, 0.5793532942101961),
                S2Point(-0.10639937102108181, -0.80810205688026326, 0.5793532942101961),
                S2Point(-0.10639937099529816, -0.80810205676695701, 0.57935329437297478),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    @Test
    fun bug10() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10592889932808099, -0.80701394501854917, 0.58095400922339757),
                S2Point(-0.10592787800899696, -0.8070140771413753, 0.58095401191158469),
                S2Point(-0.1059270044681431, -0.80701419014619669, 0.58095401421031945),
                S2Point(-0.10592685562894633, -0.80701420940058122, 0.58095401460194696),
                S2Point(-0.10592685502239066, -0.80701420947920588, 0.58095401460332308),
                S2Point(-0.10592681668594067, -0.80701421444855337, 0.5809540146902914),
                S2Point(-0.10592586497682262, -0.8070143378130904, 0.58095401684902004),
                S2Point(-0.10592586434121586, -0.80701433789547994, 0.58095401685046155),
                S2Point(-0.10592585898876766, -0.80701428569270217, 0.58095409034224832),
                S2Point(-0.10592585898876755, -0.80701428569270128, 0.58095409034224987),
                S2Point(-0.10592571912106936, -0.8070129215545373, 0.58095601078971082),
                S2Point(-0.10592571912106795, -0.80701292155452331, 0.58095601078973025),
                S2Point(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
                S2Point(-0.10592546630689463, -0.80701045544795602, 0.58095948256771723),
                S2Point(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
                S2Point(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
                S2Point(-0.1059267844512099, -0.80700966174311928, 0.58096034476466896),
                S2Point(-0.10592686088387009, -0.80700965393230761, 0.58096034167862642),
                S2Point(-0.10592691331665709, -0.80700961093727019, 0.58096039184274961),
                S2Point(-0.10592705773734933, -0.80700947507458121, 0.58096055423665138),
                S2Point(-0.10592721940752658, -0.80700934249808198, 0.58096070892049412),
                S2Point(-0.10592756003095027, -0.80700933299293154, 0.58096066001769275),
                S2Point(-0.10592832507751106, -0.80700935762745474, 0.58096048630521868),
                S2Point(-0.1059284165295875, -0.80701007424011018, 0.58095947418602778),
                S2Point(-0.10592841614913188, -0.80701007428931704, 0.58095947418704452),
                S2Point(-0.10592864947042728, -0.8070119434176124, 0.58095683523192998),
                S2Point(-0.1059286884898481, -0.80701225600079662, 0.58095639390519271),
                S2Point(-0.10592868927069989, -0.80701225581371527, 0.58095639402269295),
                S2Point(-0.10592869427137827, -0.80701225619024619, 0.58095639258785126),
                S2Point(-0.10592869791375134, -0.80701225804491505, 0.58095638934738025),
                S2Point(-0.10592869922184817, -0.80701226088076483, 0.5809563851695615),
                S2Point(-0.10592869922184843, -0.80701226088076705, 0.58095638516955805),
                S2Point(-0.10592869784516552, -0.80701226393793402, 0.58095638117383475),
                S2Point(-0.10592869415258396, -0.80701226639725276, 0.58095637843085768),
                S2Point(-0.10592868991437976, -0.80701226741266929, 0.58095637779310561),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10592564460843924, -0.80700972122716552, 0.58096046996257766),
                S2Point(-0.10592539435053176, -0.80700975987840939, 0.58096046190138972),
                S2Point(-0.10592547496472972, -0.80701045435596641, 0.58095948250602925),
                S2Point(-0.10592546630689462, -0.80701045544795591, 0.58095948256771723),
                S2Point(-0.10592546630693271, -0.80701045544826022, 0.58095948256728758),
                S2Point(-0.1059254749287661, -0.80701045440038255, 0.5809594824508878),
                S2Point(-0.10592572778318898, -0.80701292050174633, 0.58095601067279068),
                S2Point(-0.1059257191207934, -0.80701292155455673, 0.58095601078973391),
                S2Point(-0.1059257194541381, -0.80701292151405679, 0.58095601078521419),
                S2Point(-0.10592572778319062, -0.80701292050176254, 0.58095601067276803),
                S2Point(-0.10592586765088864, -0.80701428463992497, 0.58095409022530931),
                S2Point(-0.10592585899855227, -0.80701428569151201, 0.58095409034211776),
                S2Point(-0.10592585898857355, -0.80701428569272593, 0.58095409034225098),
                S2Point(-0.10592586765088888, -0.80701428463992686, 0.58095409022530675),
                S2Point(-0.10592587247896063, -0.80701433172842685, 0.58095402393347073),
                S2Point(-0.10592681605007616, -0.80701420941876889, 0.58095402179319922),
                S2Point(-0.10592685438651758, -0.80701420444942229, 0.58095402170623067),
                S2Point(-0.10592685499307326, -0.80701420437079774, 0.58095402170485466),
                S2Point(-0.10592685562894634, -0.80701420940058122, 0.58095401460194696),
                S2Point(-0.10592685499689927, -0.80701420437030225, 0.58095402170484534),
                S2Point(-0.10592700383609792, -0.80701418511591771, 0.58095402131321794),
                S2Point(-0.10592787737695626, -0.80701407211109533, 0.58095401901448296),
                S2Point(-0.10592889869604118, -0.80701393998826909, 0.58095401632629584),
                S2Point(-0.10592889996012077, -0.80701395004882903, 0.58095400212049919),
                S2Point(-0.10592787864104941, -0.80701408217165349, 0.58095400480868631),
                S2Point(-0.10592787800903029, -0.80701407714164064, 0.58095401191120999),
                S2Point(-0.10592787864103763, -0.80701408217165482, 0.5809540048086862),
                S2Point(-0.10592700510019466, -0.80701419517647521, 0.58095400710742118),
                S2Point(-0.1059270044681431, -0.80701419014619669, 0.58095401421031934),
                S2Point(-0.10592700510018833, -0.8070141951764761, 0.58095400710742118),
                S2Point(-0.10592685626275877, -0.80701421443063182, 0.58095400749904391),
                S2Point(-0.10592685565826369, -0.80701421450898914, 0.58095400750041526),
                S2Point(-0.10592685502239063, -0.80701420947920566, 0.58095401460332308),
                S2Point(-0.10592685565826078, -0.80701421450898947, 0.58095400750041526),
                S2Point(-0.10592681732181129, -0.80701421947833718, 0.58095400758738369),
                S2Point(-0.10592681668594069, -0.80701421444855348, 0.58095401469029151),
                S2Point(-0.10592681732180521, -0.80701421947833796, 0.58095400758738369),
                S2Point(-0.10592586561269894, -0.80701434284287321, 0.58095400974611222),
                S2Point(-0.10592586497746249, -0.80701433781815202, 0.58095401684187198),
                S2Point(-0.10592586561268771, -0.80701434284287465, 0.58095400974611222),
                S2Point(-0.10592586497708102, -0.80701434292526464, 0.58095400974755396),
                S2Point(-0.10592586434121586, -0.80701433789548005, 0.58095401685046166),
                S2Point(-0.10592585567909471, -0.80701433894825569, 0.58095401696740323),
                S2Point(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011),
                S2Point(-0.10592571045894811, -0.80701292260731206, 0.58095601090665361),
                S2Point(-0.10592571912060067, -0.80701292155459425, 0.58095601078971715),
                S2Point(-0.10592571878923682, -0.80701292159485349, 0.58095601079421),
                S2Point(-0.10592571045894694, -0.80701292260730051, 0.58095601090666993),
                S2Point(-0.10592545760452345, -0.80701045650593073, 0.58095948268477515),
                S2Point(-0.10592545764454649, -0.80701045650106651, 0.58095948268423492),
                S2Point(-0.10592537647753246, -0.80700975726109381, 0.58096046879584118),
                S2Point(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
                S2Point(-0.10592538413784101, -0.80700975119062324, 0.58096047583161736),
                S2Point(-0.10592564339592514, -0.80700971114934217, 0.58096048418271495),
                S2Point(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
                S2Point(-0.10592564496449927, -0.80700971099098684, 0.58096048411668999),
                S2Point(-0.10592678502227458, -0.80700965660628099, 0.58096035179610783),
                S2Point(-0.10592678388014524, -0.80700966687995779, 0.58096033773323019),
            ),
            listOf(
                S2Point(-0.10592585898876757, -0.80701428569270128, 0.58095409034224987),
                S2Point(-0.10592585897888845, -0.80701428569390288, 0.58095409034238166),
                S2Point(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011),
            ),
            listOf(
                S2Point(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
                S2Point(-0.10592546623958927, -0.8070104554564449, 0.58095948256819674),
                S2Point(-0.10592546626662946, -0.80701045545303429, 0.580959482568004),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(a)}" }
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(b)}" }
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Inconsistent loop orientations detected
        logger.debug { "\nS2Polygon: ${S2TextParser.toString(c)}" }
    }

    @Test
    fun bug11() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
                S2Point(-0.10727349807040805, -0.80875763112192245, 0.57827631000568813),
                S2Point(-0.10727349807040625, -0.80875763112192278, 0.57827631000568813),
            ),
            listOf(
                S2Point(-0.1072729603486537, -0.80875606054879057, 0.57827860629945249),
                S2Point(-0.10727299870478688, -0.80875633377729705, 0.57827821705818028),
                S2Point(-0.10727299875560981, -0.80875633413933223, 0.57827821654242495),
                S2Point(-0.10727309272230967, -0.80875700360375646, 0.57827726282438607),
                S2Point(-0.10727318660000487, -0.80875767243400742, 0.57827631000742785),
                S2Point(-0.10727349802669105, -0.80875763101356435, 0.57827631016534387),
                S2Point(-0.10727349803435525, -0.80875763107087817, 0.57827631008376468),
                S2Point(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
                S2Point(-0.1072734980420204, -0.80875763112819909, 0.57827631000217561),
                S2Point(-0.10727318657570066, -0.80875767255391384, 0.57827630984423972),
                S2Point(-0.10727318651657966, -0.80875767256177711, 0.57827630984420975),
                S2Point(-0.10727318650891528, -0.80875767250445951, 0.57827630992579371),
                S2Point(-0.10727318640981781, -0.80875767251785957, 0.57827630992543622),
                S2Point(-0.10727309252411468, -0.80875700363055636, 0.57827726282367087),
                S2Point(-0.10727299855741491, -0.8087563341661328, 0.57827821654170874),
                S2Point(-0.10727299850659211, -0.8087563338040985, 0.57827821705746318),
                S2Point(-0.10727296014242577, -0.80875606051836801, 0.57827860638025652),
                S2Point(-0.10727296024152315, -0.80875606050496729, 0.57827860638061501),
                S2Point(-0.10727296023340849, -0.8087560604477102, 0.57827860646219797),
                S2Point(-0.10727348576547496, -0.80875598914629976, 0.57827860869282954),
                S2Point(-0.1072734857817042, -0.80875598926081438, 0.57827860852966395),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.1072734857735896, -0.80875598920355718, 0.5782786086112468),
                S2Point(-0.10727348576547457, -0.80875598914629976, 0.57827860869282954),
                S2Point(-0.10727839137361543, -0.80875532356817348, 0.57827862950694298),
                S2Point(-0.10727839137881608, -0.80875532356471602, 0.57827862951081388),
                S2Point(-0.10727839143632178, -0.80875532355090063, 0.5782786295194674),
                S2Point(-0.10727839149361706, -0.80875532355509905, 0.57827862950296649),
                S2Point(-0.1072783915353497, -0.80875532357618651, 0.57827862946573261),
                S2Point(-0.10727839154773799, -0.80875532360290581, 0.57827862942606567),
                S2Point(-0.10727848921795155, -0.80875531035110082, 0.57827862984032907),
                S2Point(-0.1072784892332832, -0.80875531046514559, 0.57827862967798682),
                S2Point(-0.10727971608197531, -0.8087551454635169, 0.57827863284376713),
                S2Point(-0.10727986275126807, -0.80875539440654376, 0.57827825747332484),
                S2Point(-0.10727959167812619, -0.80875599171505064, 0.57827747239052929),
                S2Point(-0.10727974196569352, -0.80875625444235633, 0.57827707706958686),
                S2Point(-0.10727993501555312, -0.80875677560355186, 0.57827631237878363),
                S2Point(-0.10727870858143702, -0.80875693828645479, 0.57827631237896882),
                S2Point(-0.1072787085493927, -0.80875693804871851, 0.5782763127174031),
                S2Point(-0.10727615977928232, -0.80875727704955946, 0.57827631143112901),
                S2Point(-0.10727615977915911, -0.80875727704957578, 0.57827631143112901),
                S2Point(-0.10727349803435751, -0.80875763107088128, 0.57827631008375968),
                S2Point(-0.10727349803435574, -0.80875763107088183, 0.57827631008375979),
                S2Point(-0.10727318656803594, -0.80875767249659658, 0.57827630992582391),
                S2Point(-0.10727318650891531, -0.80875767250445962, 0.57827630992579382),
                S2Point(-0.10727309262321218, -0.80875700361715641, 0.57827726282402847),
                S2Point(-0.10727299865651231, -0.80875633415273218, 0.57827821654206735),
                S2Point(-0.10727299860568951, -0.80875633379069789, 0.57827821705782179),
                S2Point(-0.10727296024152314, -0.80875606050496718, 0.57827860638061501),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    @Test
    fun bug12() {
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10772916872905106, -0.80699542608967267, 0.58064861015531188),
                S2Point(-0.10772916892726483, -0.80699542606300401, 0.58064861015560143),
                S2Point(-0.10772916892726613, -0.80699542606301333, 0.58064861015558844),
                S2Point(-0.10772916872905235, -0.806995426089682, 0.58064861015529889),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.10772916872905348, -0.80699542608969022, 0.58064861015528724),
                S2Point(-0.10772916892726496, -0.80699542606300489, 0.58064861015559999),
                S2Point(-0.10772930108168739, -0.80699639165138115, 0.58064724364290399),
                S2Point(-0.10772930088347589, -0.80699639167806647, 0.58064724364259113),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    @Test
    fun bug13() {
        // This test exercises a rare special case in GetCrossedVertexIndex where
        // two crossing edge chains snap to a different permutation of the same
        // vertices.  In this example one input edge crosses another edge from right
        // to left, the first edge snaps to BCD and the second snaps to ABDC, and
        // triangle BCD is CCW.  Since BCD is to the right of BD, this means that
        // the first edge has not yet crossed the second at vertex B, leaving C or D
        // as the possible crossing vertices.
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.38306437985388492, -0.74921955334206214, 0.54030708099846292),
                S2Point(-0.3830643798552798, -0.74921955334134249, 0.5403070809984718),
                S2Point(-0.38306437985529124, -0.74921955334136414, 0.54030708099843361),
                S2Point(-0.38306437985389635, -0.74921955334208379, 0.54030708099842473),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.38306437985390962, -0.74921955334210588, 0.54030708099838465),
                S2Point(-0.38306437985527797, -0.74921955334134205, 0.54030708099847369),
                S2Point(-0.38306437985527941, -0.74921955334134405, 0.54030708099847014),
                S2Point(-0.38306437985391095, -0.74921955334210777, 0.54030708099838098),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    @Test
    fun bug14() {
        // This test exercises another rare case where the crossing vertices chosen
        // by GetCrossedVertexIndex() are not ordered correctly along the edge being
        // crossed.  This is handled by adding extra edges to the output in order to
        // link up the crossings in the correct order.
        val aVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.3837392878495085, -0.7477800800281974, 0.5418201831546835),
                S2Point(-0.38373928785696076, -0.7477800800212292, 0.54182018315902258),
                S2Point(-0.38373928785701278, -0.74778008002124685, 0.5418201831589613),
                S2Point(-0.38373928785703426, -0.7477800800212544, 0.54182018315893576),
                S2Point(-0.38373947205489456, -0.74778014227795497, 0.5418199667802881),
                S2Point(-0.38373947204434411, -0.74778014228781997, 0.54181996677414512),
                S2Point(-0.38373947205872994, -0.74778014228185352, 0.54181996677219124),
                S2Point(-0.38373947218468357, -0.74778014288930306, 0.54181996584462788),
                S2Point(-0.3837396702525171, -0.74778021044361542, 0.54181973233114322),
                S2Point(-0.38373967023137123, -0.74778021046333043, 0.54181973231891067),
                S2Point(-0.38373947216030285, -0.74778014290791484, 0.54181996583620895),
                S2Point(-0.38373947217087578, -0.74778014289805739, 0.54181996584232528),
                S2Point(-0.38373947215649007, -0.74778014290402395, 0.54181996584427927),
                S2Point(-0.3837394720305386, -0.74778014229658485, 0.5418199667718262),
                S2Point(-0.38373928783585998, -0.74778008004095942, 0.54182018314673686),
                S2Point(-0.38373928784641037, -0.7477800800310942, 0.54182018315287972),
                S2Point(-0.38373928783578648, -0.74778008004093421, 0.54182018314682368),
                S2Point(-0.383739287835765, -0.74778008004092666, 0.54182018314684921),
            ),
        )
        val bVertices: List<List<S2Point>> = listOf(
            listOf(
                S2Point(-0.38373923813692823, -0.7477800632164362, 0.54182024156551456),
                S2Point(-0.3837392878569364, -0.74778008002122087, 0.54182018315905123),
                S2Point(-0.38373928784640354, -0.74778008003106944, 0.54182018315291858),
                S2Point(-0.38373928784638789, -0.74778008003108642, 0.54182018315290648),
                S2Point(-0.38373928784638023, -0.74778008003109453, 0.54182018315290048),
                S2Point(-0.38373928783692102, -0.74778008004124585, 0.54182018314559),
                S2Point(-0.38373928783691913, -0.74778008004124541, 0.54182018314559188),
                S2Point(-0.38373928784636568, -0.74778008003110774, 0.54182018315289271),
                S2Point(-0.38373928784637329, -0.74778008003109953, 0.54182018315289848),
                S2Point(-0.38373928783583561, -0.74778008004095109, 0.5418201831467655),
                S2Point(-0.38373923811582744, -0.74778006323616641, 0.54182024155322883),
                S2Point(-0.38373857650312843, -0.74777983961840766, 0.54182101875399913),
                S2Point(-0.38373857652422921, -0.74777983959867744, 0.54182101876628486),
            ),
        )
        val a = S2Polygon(makeLoops(aVertices))
        val b = S2Polygon(makeLoops(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        // Given edges do not form loops (indegree != outdegree)
        assertThat(c.isEmpty())
            .withFailMessage(
                "\nS2Polygon: ${S2TextParser.toString(a)}" +
                        "\nS2Polygon: ${S2TextParser.toString(b)}"
            )
            .isFalse()
    }

    // This tests polygon-polyline intersections.
    // It covers the same edge cases as TestOperations and also adds some
    // extra tests for shared edges.
    @Test
    fun polylineIntersectionCross1() {
        for (v in 0..2) {
            polylineIntersectionSharedEdgeTest(cross1, v, 1)
            polylineIntersectionSharedEdgeTest(cross1, v + 1, -1)
            polylineIntersectionSharedEdgeTest(cross1_side_hole, v, 1)
            polylineIntersectionSharedEdgeTest(cross1_side_hole, v + 1, -1)
        }
    }

    @ParameterizedTest
    @MethodSource("testCases")
    fun polylineIntersection(test: TestCase) {

        // See comments in TestOperations about the vlue of this constant.
        val kMaxError = S1Angle.radians(1e-4)

        // This duplicates some of the tests in TestOperations by
        // converting the outline of polygon A to a polyline then intersecting
        // it with the polygon B. It then converts B to a polyline and intersects
        // it with A. It then feeds all of the results into a polygon builder and
        // tests that the output is equal to doing an intersection between A and B.
        logger.debug { "Polyline intersection test case: $test " }
        val a = makePolygon(test.a)
        val b = makePolygon(test.b)
        val expectedAAndB = makePolygon(test.aAndB)

        val points = mutableListOf<S2Point>()
        val polylines = mutableListOf<S2Polyline>()
        for (ab in 0..1) {
            val tmp = if (ab == 0) b else a
            val tmp2 = if (ab == 0) a else b

            for (l in 0 until tmp.numLoops()) {
                points.clear()
                if (tmp.loop(l).isHole()) {
                    for (v in (tmp.loop(l).numVertices downTo 0)) {
                        points.add(tmp.loop(l).vertex(v))
                    }
                } else {
                    for (v in 0..tmp.loop(l).numVertices) {
                        points.add(tmp.loop(l).vertex(v))
                    }
                }
                val polyline = S2Polyline(points)
                val intersection = tmp2.intersectWithPolyline(polyline)
                polylines.addAll(intersection)
            }
        }

        val builder = S2Builder(S2Builder.Options())
        val aAndB = S2Polygon()
        builder.startLayer(S2PolygonLayer(aAndB))
        for (polyline in polylines) {
            builder.addPolyline(polyline)
        }

        val error = S2Error()

        val isBuilt = builder.build(error)
        assertThat(isBuilt)
            .withFailMessage(error.text)
            .isTrue()
        checkEqual(aAndB, expectedAAndB, kMaxError)
    }

    @Test
    fun splitting() {
        // It takes too long to test all the polygons in debug mode, so we just pick
        // out some of the more interesting ones.
        splitAndAssemble(near10)
        splitAndAssemble(near_H3210)
        splitAndAssemble(far_H3210)
        splitAndAssemble(south_0ab)
        splitAndAssemble(south_210b)
        splitAndAssemble(south_H20abc)
        splitAndAssemble(nf1N10F2S10abc)
        splitAndAssemble(nf2N2F210S210ab)
        splitAndAssemble(far_H)
        splitAndAssemble(south_H)
        splitAndAssemble(far_H_south_H)
    }

    @Test
    fun initToCellUnionBorder() {
        // Test S2Polygon.InitToCellUnionBorder().
        // The main thing to check is that adjacent cells of different sizes get
        // merged correctly.  To do this we generate two random adjacent cells,
        // convert to polygon, and make sure the polygon only has a single loop.
        repeat(200) { iter ->
            logger.trace { "Iteration $iter" }

            // Choose a random non-leaf cell.
            val bigCell = S2Random.randomCellId(S2Random.randomInt(S2CellId.kMaxLevel))
            // Get all neighbors at some smaller level.
            val smallLevel = bigCell.level() + S2Random.randomInt(min(16, S2CellId.kMaxLevel - bigCell.level()))
            val neighbors = mutableListOf<S2CellId>()
            bigCell.appendAllNeighbors(smallLevel, neighbors)
            // Pick one at random.
            val smallCell = neighbors[S2Random.randomInt(neighbors.size)]
            // If it's diagonally adjacent, bail out.
            val edgeNeighbors = bigCell.getEdgeNeighbors()
            var diagonal = true
            for (i in 0..3) {
                if (edgeNeighbors[i].contains(smallCell)) {
                    diagonal = false
                }
            }
            logger.debug { "$iter: big_cell $bigCell small_cell $smallCell" + (if (diagonal) "  diagonal - bailing out!" else "") }
            if (!diagonal) {

                val cells = mutableListOf<S2CellId>()
                cells.add(bigCell)
                cells.add(smallCell)
                val cellUnion = S2CellUnion(cells)
                assertThat(cellUnion.numCells()).isEqualTo(2)
                val poly = S2Polygon()
                poly.initToCellUnionBorder(cellUnion)
                assertThat(poly.numLoops()).isEqualTo(1)
                // If the conversion were perfect we could test containment, but due to
                // rounding the polygon won't always exactly contain both cells.  We can
                // at least test intersection.
                assertThat(poly.mayIntersect(S2Cell(bigCell))).isTrue()
                assertThat(poly.mayIntersect(S2Cell(smallCell))).isTrue()
            }
        }
    }

    @Test
    fun unionWithAmbgiuousCrossings() {
        val aVertices = listOf(
            S2Point(0.044856812877680216, -0.80679210859571904, 0.5891301722422051),
            S2Point(0.044851868273159699, -0.80679240802900054, 0.5891301386444033),
            S2Point(0.044854246527738666, -0.80679240292188514, 0.58912996457145106)
        )
        val bVertices = listOf(
            S2Point(0.044849715793028468, -0.80679253837178111, 0.58913012401412856),
            S2Point(0.044855344598821352, -0.80679219751320641, 0.589130162266992),
            S2Point(0.044854017712818696, -0.80679210327223405, 0.58913039235179754)
        )
        val a = S2Polygon(S2Loop(aVertices))
        val b = S2Polygon(S2Loop(bVertices))
        val c = S2Polygon()
        c.initToUnion(a, b)
        assertThat(c.isEmpty()).isFalse()
    }

    @Test
    fun initToSloppySupportsEmptyPolygons() {
        val emptyPolygon = S2Polygon()
        val polygon = S2Polygon()
        polygon.initToSnapped(emptyPolygon)
        // InitToSloppy is further tested by SnapSplitsPolygon.
    }

    @Test
    fun initToSnappedDoesNotRotateVertices() {
        // This particular example came from MapFacts, but in fact InitToSnapped
        // used to cyclically rotate the vertices of all "hole" loops.
        val polygon = S2TextParser.makePolygon(
            "49.9305505:-124.8345463, 49.9307448:-124.8299657, " +
                    "49.9332101:-124.8301996, 49.9331224:-124.8341368; " +
                    "49.9311087:-124.8327042, 49.9318176:-124.8312621, " +
                    "49.9318866:-124.8334451"
        )
        val polygon2 = S2Polygon()
        polygon2.initToSnapped(polygon)

        // Check that the first vertex is the same when converted to E7.
        assertThat(S2LatLng.latitude(polygon2.loop(0).vertex(0)).e7()).isEqualTo(
            S2LatLng.latitude(
                polygon.loop(0).vertex(0)
            ).e7()
        )
        assertThat(S2LatLng.longitude(polygon2.loop(0).vertex(0)).e7()).isEqualTo(
            S2LatLng.longitude(
                polygon.loop(0).vertex(0)
            ).e7()
        )

        // Check that snapping twice doesn't rotate the vertices.
        val polygon3 = S2Polygon()
        polygon3.initToSnapped(polygon2)
        assertThat(polygon2 == polygon3).isTrue()
    }

    @Test
    fun initToSnappedWithSnapLevel() {
        val polygon = S2TextParser.makePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0")
        for (level in 0..S2CellId.kMaxLevel) {
            val snappedPolygon = S2Polygon()
            snappedPolygon.initToSnapped(polygon, level)
            assertThat(snappedPolygon.isValid()).isTrue()
            val mergeRadius =
                minOf(S1Angle.radians(S2Coords.projection.kMaxDiag.getValue(level)), SnapFunction.kMaxSnapRadius())
            assertThat(snappedPolygon.approxContains(polygon, mergeRadius)).isTrue()
        }
    }

    @Test
    fun initToSnappedIsValid_A() {
        val poly = S2TextParser.makePolygon(
            "53.1328020478452:6.39444903453293, 53.1328019:6.394449, " +
                    "53.1327091:6.3961766, 53.1313753:6.3958652, 53.1312825:6.3975924, " +
                    "53.132616:6.3979042, 53.1326161348736:6.39790423150577"
        )
        logger.debug { "\nInput: ${S2TextParser.toString(poly)}" }
        assertThat(poly.isValid()).isTrue()
        val polySnapped = S2Polygon()
        polySnapped.debugOverride = S2Debug.DISABLE
        polySnapped.initToSnapped(poly)
        logger.debug { "\nSnapped: ${S2TextParser.toString(polySnapped)}" }
        val error = S2Error()
        assertThat(polySnapped.findValidationError(error)).isFalse()
    }

    @Test
    fun initToSnappedIsValid_B() {
        val poly = S2TextParser.makePolygon(
            "51.6621651:4.9858102, 51.6620965:4.9874227, 51.662028:4.9890355, " +
                    "51.6619796006122:4.99017864445347, 51.6622335420397:4.98419752545216, " +
                    "51.6622334:4.9841975; 51.66189957578:4.99206198576131, " +
                    "51.6618911:4.9922612, 51.6618224:4.9938741, 51.6605122:4.993639, " +
                    "51.6604437:4.9952519, 51.6603751:4.9968648, 51.6603064:4.9984777, " +
                    "51.6602379:5.0000907, 51.660169:5.0017037, 51.6601003:5.0033165, " +
                    "51.6600318:5.0049298, 51.659963:5.0065427, 51.6598943:5.0081561, " +
                    "51.6612044207178:5.00839208571886, 51.6612732068132:5.00677860122814, " +
                    "51.6612732:5.0067786, 51.6613418:5.0051654, 51.6614106:5.0035525, " +
                    "51.6614793:5.0019393, 51.6615479:5.0003263, " +
                    "51.6615946694783:4.99923124520759, 51.6616389353165:4.99819106536521, " +
                    "51.6616852:4.9971, 51.6617538:4.995487, " +
                    "51.661753964726:4.99548702962593"
        )
        logger.debug { "\nInput: ${S2TextParser.toString(poly)}" }
        assertThat(poly.isValid()).isTrue()
        val polySnapped = S2Polygon()
        polySnapped.debugOverride = S2Debug.DISABLE
        polySnapped.initToSnapped(poly)
        logger.debug { "\nSnapped: ${S2TextParser.toString(polySnapped)}" }
        val error = S2Error()
        assertThat(polySnapped.findValidationError(error)).isFalse()
    }

    @Test
    fun initToSnappedIsValid_C() {
        val poly = S2TextParser.makePolygon(
            "53.5316236236404:19.5841192796855, 53.5416584:19.5915903, " +
                    "53.5416584189104:19.5915901888287; 53.5416584:19.5915903, " +
                    "53.5363122:19.62299, 53.5562817:19.6378935, 53.5616342:19.606474; " +
                    "53.5616342:19.606474, 53.5916039:19.6288326, 53.5912689:19.6307982, " +
                    "53.5925176:19.6317308, 53.5928526:19.6297652, 53.6015949:19.6362943, " +
                    "53.6015950436033:19.6362944072725, 53.6015950814439:19.6362941852262, " +
                    "53.5616342380536:19.6064737764314"
        )
        logger.debug { "\nInput: ${S2TextParser.toString(poly)}" }
        assertThat(poly.isValid()).isTrue()
        val polySnapped = S2Polygon()
        polySnapped.debugOverride = S2Debug.DISABLE
        polySnapped.initToSnapped(poly)
        logger.debug { "\nSnapped: ${S2TextParser.toString(polySnapped)}" }
        val error = S2Error()
        assertThat(polySnapped.findValidationError(error)).isFalse()
    }

    @Test
    fun initToSnappedIsValid_D() {
        val poly = S2TextParser.makePolygon(
            "52.0909316:4.8673826, 52.0909317627574:4.86738262858533, " +
                    "52.0911338452911:4.86248482549567, 52.0911337:4.8624848, " +
                    "52.0910665:4.8641176, 52.090999:4.8657502"
        )
        logger.debug { "\nInput: ${S2TextParser.toString(poly)}" }
        assertThat(poly.isValid()).isTrue()
        val polySnapped = S2Polygon()
        polySnapped.debugOverride = S2Debug.DISABLE
        polySnapped.initToSnapped(poly)
        logger.debug { "\nSnapped: ${S2TextParser.toString(polySnapped)}" }
        val error = S2Error()
        assertThat(polySnapped.findValidationError(error)).isFalse()
    }

    @Test
    fun multipleInit() {
        val polygon = S2TextParser.makePolygon("0:0, 0:2, 2:0")
        assertThat(polygon.numLoops()).isEqualTo(1)
        assertThat(polygon.numVertices()).isEqualTo(3)
        val bound1 = polygon.rectBound

        val loops = mutableListOf<S2Loop>()
        loops.add(S2TextParser.makeLoop("10:0, -10:-20, -10:20"))
        loops.add(S2TextParser.makeLoop("40:30, 20:10, 20:50"))
        polygon.initNested(loops)
        assertThat(polygon.isValid()).isTrue()
        assertThat(polygon.numLoops()).isEqualTo(2)
        assertThat(polygon.numVertices()).isEqualTo(6)
        assertThat(bound1 != polygon.rectBound).isTrue()
    }

    @Test
    fun initSingleLoop() {
        val polygon = S2Polygon(S2Loop(S2Loop.kEmpty))
        assertThat(polygon.isEmpty()).isTrue()
        polygon.init(S2Loop(S2Loop.kFull))
        assertThat(polygon.isFull()).isTrue()
        polygon.init(S2TextParser.makeLoop("0:0, 0:10, 10:0"))
        assertThat(polygon.numVertices()).isEqualTo(3)
    }

    // This test checks that S2Polygons created directly from S2Cells behave
    // identically to S2Polygons created from the vertices of those cells; this
    // previously was not the case, because S2Cells calculate their bounding
    // rectangles slightly differently, and S2Polygons created from them just
    // copied the S2Cell bounds.
    @Test
    fun testS2CellConstructorAndContains() {
        val latlng = S2LatLng.fromLatLng(S1Angle.e6(40565459), S1Angle.e6(-74645276))
        val cell = S2Cell(latlng)
        val cellAsPolygon = S2Polygon(cell)
        val empty = S2Polygon()
        val polygonCopy = S2Polygon()
        polygonCopy.initToUnion(cellAsPolygon, empty)
        assertThat(polygonCopy.contains(cellAsPolygon)).isTrue()
        assertThat(cellAsPolygon.contains(polygonCopy)).isTrue()
    }

    @Test
    fun project() {
        val polygon = makePolygon(kNear0 + kNear2)
        var projected: S2Point

        // The point inside the polygon should be projected into itself.
        var point: S2Point = S2TextParser.makePoint("1.1:0")
        projected = polygon.project(point)
        assertThat(S2PointUtil.approxEquals(point, projected)).isTrue()

        // The point is on the outside of the polygon.
        point = S2TextParser.makePoint("5.1:-2")
        projected = polygon.project(point)
        assertThat(S2PointUtil.approxEquals(S2TextParser.makePoint("5:-2"), projected)).isTrue()

        // The point is inside the hole in the polygon.
        point = S2TextParser.makePoint("-0.49:-0.49")
        projected = polygon.project(point)
        assertThat(
            S2PointUtil.approxEquals(
                S2TextParser.makePoint("-0.5:-0.5"),
                projected,
                S1Angle.radians(1e-6)
            )
        ).isTrue()

        point = S2TextParser.makePoint("0:-3")
        projected = polygon.project(point)
        assertThat(S2PointUtil.approxEquals(S2TextParser.makePoint("0:-2"), projected)).isTrue()
    }

    @Test
    fun getDistance() {
        // The empty and full loops don't have boundaries.
        testDistanceMethods(empty, S2Point(0, 1, 0), S2Point())
        testDistanceMethods(full, S2Point(0, 1, 0), S2Point())

        // A polygon consisting of two nested rectangles centered around
        // S2LatLng(0,0).  Note that because lines of latitude are curved on the
        // sphere, it is not straightforward to project points onto any edge except
        // along the equator.  (The equator is the only line of latitude that is
        // also a geodesic.)
        val nested = S2TextParser.makePolygon("3:1, 3:-1, -3:-1, -3:1; 4:2, 4:-2, -4:-2, -4:2;")

        // All points on the boundary of the polygon should be at distance zero.
        for (i in 0 until nested.numLoops()) {
            val loop = nested.loop(i)
            for (j in 0 until loop.numVertices) {
                // A vertex.
                testDistanceMethods(nested, loop.vertex(j), S2Point())
                // A point along an edge.
                testDistanceMethods(
                    nested,
                    S2EdgeDistances.interpolate(S2Random.randomDouble(), loop.vertex(j), loop.vertex(j + 1)),
                    S2Point()
                )
            }
        }
        // A point outside the outer shell that projects to an edge.
        testDistanceMethods(nested, S2LatLng.fromDegrees(0.0, -4.7).toPoint(), S2LatLng.fromDegrees(0, -2).toPoint())
        // A point outside the outer shell that projects to a vertex.
        testDistanceMethods(nested, S2LatLng.fromDegrees(6, -3).toPoint(), S2LatLng.fromDegrees(4, -2).toPoint())
        // A point inside the polygon that projects to an outer edge.
        testDistanceMethods(nested, S2LatLng.fromDegrees(0.0, 1.7).toPoint(), S2LatLng.fromDegrees(0, 2).toPoint())
        // A point inside the polygon that projects to an inner vertex.
        testDistanceMethods(nested, S2LatLng.fromDegrees(-3.3, -1.3).toPoint(), S2LatLng.fromDegrees(-3, -1).toPoint())
        // A point inside the inner hole.
        testDistanceMethods(nested, S2LatLng.fromDegrees(0.0, 0.1).toPoint(), S2LatLng.fromDegrees(0, 1).toPoint())
    }

    @Test
    fun area() {
        assertThat(empty.getArea()).isEqualTo(0.0)
        assertThat(full.getArea()).isEqualTo(4 * M_PI)
        assertThat(south_H.getArea()).isEqualTo(2 * M_PI)
        assertThat(far_H_south_H.getArea()).isCloseTo(M_PI, Offset.offset(1e-15))

        val twoShells = makePolygon(kCross1SideHole + kCrossCenterHole)
        assertThat(twoShells.getArea()).isEqualTo(twoShells.loop(0).area + twoShells.loop(1).area)

        val holeyShell = makePolygon(kCross1 + kCrossCenterHole)
        assertThat(holeyShell.getArea()).isEqualTo(holeyShell.loop(0).area - holeyShell.loop(1).area)
    }

    @Test
    fun uninitializedIsValid() {
        val polygon = S2Polygon()
        assertThat(polygon.isValid()).isTrue()
    }

    @Test
    fun fullPolygonShape() {
        val shape = S2Polygon.Shape(polygon = full)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isTrue()
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.chain(0).start).isEqualTo(0)
        assertThat(shape.chain(0).length).isEqualTo(0)
        assertThat(shape.getReferencePoint().contained).isTrue()
    }

    @Test
    fun emptyPolygonShape() {
        val shape = S2Polygon.Shape(polygon = empty)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    fun testPolygonShape(polygon: S2Polygon) {
        require(!polygon.isFull())
        val shape = S2Polygon.Shape(polygon = polygon)
        assertThat(shape.polygon).isEqualTo(polygon)
        assertThat(shape.numEdges).isEqualTo(polygon.numVertices())
        assertThat(shape.numChains).isEqualTo(polygon.numLoops())
        var e = 0
        for (i in 0 until polygon.numLoops()) {
            val loop_i = polygon.loop(i)
            assertThat(shape.chain(i).start).isEqualTo(e)
            assertThat(shape.chain(i).length).isEqualTo(loop_i.numVertices)
            for (j in 0 until loop_i.numVertices) {
                val edge = shape.edge(e)
                assertThat(edge.v0).isEqualTo(loop_i.orientedVertex(j))
                assertThat(edge.v1).isEqualTo(loop_i.orientedVertex(j + 1))
                ++e
            }
        }
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isEqualTo(polygon.contains(S2PointUtil.origin()))
    }

    @Test
    fun oneLoopPolygonShape() {
        testPolygonShape(near0)
    }

    @Test
    fun severalLoopPolygonShape() {
        testPolygonShape(near_3210)
    }

    @Test
    fun manyLoopPolygonShape() {
        val kNumLoops = 100
        val kNumVerticesPerLoop = 6
        val polygon = S2Polygon()
        S2Factory.concentricLoopsPolygon(S2Point(1, 0, 0), kNumLoops, kNumVerticesPerLoop, polygon)
        testPolygonShape(polygon)
    }

    @Test
    fun pointInBigLoop() {
        // This code used to demonstrate a bug in S2ShapeIndex.
        val center = S2LatLng.fromRadians(0.3, 2.0)
        val radius = S1Angle.degrees(80)
        val poly = S2Polygon(S2Loop.makeRegularLoop(center.toPoint(), radius, 10))
        assertThat(poly.mayIntersect(S2Cell(S2CellId.fromLatLng(center)))).isTrue()
    }


    @Test
    fun indexContainsOnePolygonShape() {
        val index = near0.index
        assertThat(index.nextNewShapeId()).isEqualTo(1)
        val shape = (index.shape(0) as S2Polygon.Shape)
        assertThat(shape.polygon).isEqualTo(near0)
    }

    @Test
    fun polygonPolygonDistance() {
        // Verify that the example code for S2Polygon.index() actually works.
        val polygon1 = near0.clone()
        val polygon2 = far_10.clone()
        val query = S2ClosestEdgeQuery(polygon1.index)
        val target = S2ClosestEdgeQuery.ShapeIndexTarget(polygon2.index)
        val distance = query.getDistance(target);
        assertThat(distance).isGreaterThan(S1ChordAngle(S1Angle.degrees(175)))
    }

    @Test
    fun test() {
        testRelation(nf1N10F2S10abc, nf2N2F210S210ab, contains = false, contained = false, intersects = true)
    }

    companion object {
        val logger = KotlinLogging.logger(S2PolygonUnitTest::class.java.name)

        @JvmStatic
        fun testCases(): Stream<TestCase> = listOf(
            // Two triangles that share an edge.
            TestCase(
                a = "4:2, 3:1, 3:3;",
                b = "3:1, 2:2, 3:3;",
                aAndB = "",  // and
                aOrB = "4:2, 3:1, 2:2, 3:3;",  // or
                aMinusB = "4:2, 3:1, 3:3;",  // minus
                aXorB = "4:2, 3:1, 2:2, 3:3;"  // xor
            ),
            // Two vertical bars and a horizontal bar connecting them.
            TestCase(
                a = "0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;",
                b = "1:1, 1:4, 2:4, 2:1;",
                aAndB = "1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;",
                aOrB = "0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;",
                aMinusB = "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;",
                aXorB = "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;1:2, 1:3, 2:3, 2:2"
            ),
            // Two vertical bars and two horizontal bars.
            TestCase(
                a = "1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;",
                b = "-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;",
                aAndB = "1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91; -1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",
                aOrB = "-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, -1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, 2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88;  0:90, 0:91, 1:91, 1:90;",
                aMinusB = "1:88, 1:89, 2:89, 2:88;   1:90, 1:91, 2:91, 2:90; 1:92, 1:93, 2:93, 2:92;   -1:88, -1:89, 0:89, 0:88;  -1:90, -1:91, 0:91, 0:90;   -1:92, -1:93, 0:93, 0:92;",
                aXorB = "1:88, 1:89, 2:89, 2:88;   -1:88, -1:89, 0:89, 0:88;   1:90, 1:91, 2:91, 2:90;   -1:90, -1:91, 0:91, 0:90;  1:92, 1:93, 2:93, 2:92;   -1:92, -1:93, 0:93, 0:92;  -2:89, -2:90, -1:90, -1:89;   -2:91, -2:92, -1:92, -1:91;  0:89, 0:90, 1:90, 1:89;   0:91, 0:92, 1:92, 1:91;  2:89, 2:90, 3:90, 3:89;   2:91, 2:92, 3:92, 3:91;"
            ),
            // Two interlocking square doughnuts.
            TestCase(
                a = "-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;",
                b = "-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;",
                aAndB = "-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;",  // and
                aOrB = "-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93;  0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92;   -2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",
                aMinusB = "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93;   -1:-90, -1:-89, 0:-89, 0:-90;",
                aXorB = "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93;   -3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88, -2:-88, -2:-90, -1:-90, -1:-91;   -1:-90, -1:-89, 0:-89, 0:-90;   1:-91, 0:-91, 0:-90, 1:-90;"
            ),
            // An incredibly thin triangle intersecting a square, such that the two
            // intersection points of the triangle with the square are identical.
            // This results in a degenerate loop that needs to be handled correctly.
            TestCase(
                "10:44, 10:46, 12:46, 12:44;",
                "11:45, 89:45.00000000000001, 90:45;",
                "",  // Empty intersection!
                // Original square with extra vertex, and triangle disappears (due to
                // default vertex_merge_radius of S2.kIntersectionMergeRadius).
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;",  // or
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;",  // minus
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;",  // xor
            ),
        ).stream()
    }

    data class TestCase(
        val a: String,
        val b: String,
        val aAndB: String,
        val aOrB: String,
        val aMinusB: String,
        val aXorB: String,
    )

}


/*


// Verifies that S2Polygon does not destroy or replace pointers to S2Loop, so
// caller can rely on using raw pointers.
fun LoopPointers() {
  vector<unique_ptr<S2Loop>> loops
  loops.emplace_back(S2TextParser.MakeLoop("4:4, 4:6, 6:6, 6:4"))
  loops.emplace_back(S2TextParser.MakeLoop("3:3, 3:7, 7:7, 7:3"))
  loops.emplace_back(S2TextParser.MakeLoop("2:2, 2:8, 8:8, 8:2"))
  loops.emplace_back(S2TextParser.MakeLoop("1:1, 1:9, 9:9, 9:1"))
  loops.emplace_back(S2TextParser.MakeLoop("10:10, 15:15, 20:10"))
  loops.emplace_back(S2TextParser.MakeLoop("-1:-1, -9:-1, -9:-9, -1:-9"))
  loops.emplace_back(S2TextParser.MakeLoop("-5:-5, -6:-5, -6:-6, -5:-6"))

  std.set<const S2Loop*> loops_raw_ptrs
  for (auto& loop : loops) {
    loops_raw_ptrs.insert(loop)
  }
  S2Polygon polygon(std.move(loops))

  // Check that loop pointers didn't change (but could've gotten reordered).
  assertThat(polygon.num_loops()).isEqualTo(loops_raw_ptrs.size())
  for (int i = 0; i < polygon.num_loops(); i++) {
     EXPECT_EQ(1, loops_raw_ptrs.count(polygon.loop(i))) << "loop " << i
  }
}



*/

