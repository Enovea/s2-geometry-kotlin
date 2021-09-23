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

import dilivia.s2.S2WedgeRelations.WedgeRelation.WEDGE_EQUALS
import dilivia.s2.S2WedgeRelations.WedgeRelation.WEDGE_IS_DISJOINT
import dilivia.s2.S2WedgeRelations.WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED
import dilivia.s2.S2WedgeRelations.WedgeRelation.WEDGE_PROPERLY_CONTAINS
import dilivia.s2.S2WedgeRelations.WedgeRelation.WEDGE_PROPERLY_OVERLAPS
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test


class S2WedgeRelationsUnitTest {

    fun testWedge(
        a0: S2Point,
        ab1: S2Point,
        a2: S2Point,
        b0: S2Point,
        b2: S2Point,
        contains: Boolean,
        intersects: Boolean,
        wedge_relation: S2WedgeRelations.WedgeRelation
    ) {
        a0.normalize()
        ab1.normalize()
        a2.normalize()
        b0.normalize()
        b2.normalize()
        assertThat(S2WedgeRelations.wedgeContains(a0, ab1, a2, b0, b2)).isEqualTo(contains)
        assertThat(S2WedgeRelations.wedgeIntersects(a0, ab1, a2, b0, b2)).isEqualTo(intersects)
        assertThat(S2WedgeRelations.getWedgeRelation(a0, ab1, a2, b0, b2)).isEqualTo(wedge_relation)
    }

    @Test
    fun wedges() {
        // For simplicity, all of these tests use an origin of (0, 0, 1).
        // This shouldn't matter as long as the lower-level primitives are
        // implemented correctly.

        // Intersection in one wedge.
        testWedge(
            S2Point(-1, 0, 10),
            S2Point(0, 0, 1),
            S2Point(1, 2, 10),
            S2Point(0, 1, 10),
            S2Point(1, -2, 10),
            false,
            true,
            WEDGE_PROPERLY_OVERLAPS
        )
        // Intersection in two wedges.
        testWedge(
            S2Point(-1, -1, 10),
            S2Point(0, 0, 1),
            S2Point(1, -1, 10),
            S2Point(1, 0, 10),
            S2Point(-1, 1, 10),
            false,
            true,
            WEDGE_PROPERLY_OVERLAPS
        )

        // Normal containment.
        testWedge(
            S2Point(-1, -1, 10),
            S2Point(0, 0, 1),
            S2Point(1, -1, 10),
            S2Point(-1, 0, 10),
            S2Point(1, 0, 10),
            true,
            true,
            WEDGE_PROPERLY_CONTAINS
        )
        // Containment with equality on one side.
        testWedge(
            S2Point(2, 1, 10),
            S2Point(0, 0, 1),
            S2Point(-1, -1, 10),
            S2Point(2, 1, 10),
            S2Point(1, -5, 10),
            true,
            true,
            WEDGE_PROPERLY_CONTAINS
        )
        // Containment with equality on the other side.
        testWedge(
            S2Point(2, 1, 10),
            S2Point(0, 0, 1),
            S2Point(-1, -1, 10),
            S2Point(1, -2, 10),
            S2Point(-1, -1, 10),
            true,
            true,
            WEDGE_PROPERLY_CONTAINS
        )

        // Containment with equality on both sides.
        testWedge(
            S2Point(-2, 3, 10),
            S2Point(0, 0, 1),
            S2Point(4, -5, 10),
            S2Point(-2, 3, 10),
            S2Point(4, -5, 10),
            true,
            true,
            WEDGE_EQUALS
        )

        // Disjoint with equality on one side.
        testWedge(
            S2Point(-2, 3, 10),
            S2Point(0, 0, 1),
            S2Point(4, -5, 10),
            S2Point(4, -5, 10),
            S2Point(-2, -3, 10),
            false,
            false,
            WEDGE_IS_DISJOINT
        )
        // Disjoint with equality on the other side.
        testWedge(
            S2Point(-2, 3, 10),
            S2Point(0, 0, 1),
            S2Point(0, 5, 10),
            S2Point(4, -5, 10),
            S2Point(-2, 3, 10),
            false,
            false,
            WEDGE_IS_DISJOINT
        )
        // Disjoint with equality on both sides.
        testWedge(
            S2Point(-2, 3, 10),
            S2Point(0, 0, 1),
            S2Point(4, -5, 10),
            S2Point(4, -5, 10),
            S2Point(-2, 3, 10),
            false,
            false,
            WEDGE_IS_DISJOINT
        )

        // B contains A with equality on one side.
        testWedge(
            S2Point(2, 1, 10),
            S2Point(0, 0, 1),
            S2Point(1, -5, 10),
            S2Point(2, 1, 10),
            S2Point(-1, -1, 10),
            false,
            true,
            WEDGE_IS_PROPERLY_CONTAINED
        )
        // B contains A with equality on the other side.
        testWedge(
            S2Point(2, 1, 10),
            S2Point(0, 0, 1),
            S2Point(1, -5, 10),
            S2Point(-2, 1, 10),
            S2Point(1, -5, 10),
            false,
            true,
            WEDGE_IS_PROPERLY_CONTAINED
        )
    }

}
