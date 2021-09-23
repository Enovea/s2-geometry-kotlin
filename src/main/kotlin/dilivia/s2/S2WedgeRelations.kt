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

// Defines functions for determining the relationship between two angles
// ("wedges") that share a common vertex.
object S2WedgeRelations {
    // Detailed relation from one wedge A to another wedge B.
    enum class WedgeRelation {
        WEDGE_EQUALS,                 // A and B are equal.
        WEDGE_PROPERLY_CONTAINS,      // A is a strict superset of B.
        WEDGE_IS_PROPERLY_CONTAINED,  // A is a strict subset of B.
        WEDGE_PROPERLY_OVERLAPS,      // A-B, B-A, and A intersect B are non-empty.
        WEDGE_IS_DISJOINT,            // A and B are disjoint.
    }
    
    // Given an edge chain (x0, x1, x2), the wedge at x1 is the region to the
    // left of the edges.  More precisely, it is the set of all rays from x1x0
    // (inclusive) to x1x2 (exclusive) in the *clockwise* direction.
    //
    // The following functions compare two *non-empty* wedges that share the
    // same middle vertex: A=(a0, ab1, a2) and B=(b0, ab1, b2).

    // Returns the relation from wedge A to B.
    // REQUIRES: A and B are non-empty.
    fun getWedgeRelation(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): WedgeRelation {
        // There are 6 possible edge orderings at a shared vertex (all
        // of these orderings are circular, i.e. abcd == bcda):
        //
        //  (1) a2 b2 b0 a0: A contains B
        //  (2) a2 a0 b0 b2: B contains A
        //  (3) a2 a0 b2 b0: A and B are disjoint
        //  (4) a2 b0 a0 b2: A and B intersect in one wedge
        //  (5) a2 b2 a0 b0: A and B intersect in one wedge
        //  (6) a2 b0 b2 a0: A and B intersect in two wedges
        //
        // We do not distinguish between 4, 5, and 6.
        // We pay extra attention when some of the edges overlap.  When edges
        // overlap, several of these orderings can be satisfied, and we take
        // the most specific.
        if (a0 == b0 && a2 == b2) return WEDGE_EQUALS

        if (S2Predicates.orderedCCW(a0, a2, b2, ab1)) {
            // The cases with this vertex ordering are 1, 5, and 6,
            // although case 2 is also possible if a2 == b2.
            if (S2Predicates.orderedCCW(b2, b0, a0, ab1)) return WEDGE_PROPERLY_CONTAINS;

            // We are in case 5 or 6, or case 2 if a2 == b2.
            return if(a2 == b2) WEDGE_IS_PROPERLY_CONTAINED else WEDGE_PROPERLY_OVERLAPS;
        }

        // We are in case 2, 3, or 4.
        if (S2Predicates.orderedCCW(a0, b0, b2, ab1)) return WEDGE_IS_PROPERLY_CONTAINED;
        return if(S2Predicates.orderedCCW(a0, b0, a2, ab1)) WEDGE_IS_DISJOINT else WEDGE_PROPERLY_OVERLAPS
    }

    // Returns true if wedge A contains wedge B.  Equivalent to but faster than
    // GetWedgeRelation() == WEDGE_PROPERLY_CONTAINS || WEDGE_EQUALS.
    // REQUIRES: A and B are non-empty.
    fun wedgeContains(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        // For A to contain B (where each loop interior is defined to be its left
        // side), the CCW edge order around ab1 must be a2 b2 b0 a0.  We split
        // this test into two parts that test three vertices each.
        return (S2Predicates.orderedCCW(a2, b2, b0, ab1) && S2Predicates.orderedCCW(b0, a0, a2, ab1))
    }
    
    // Returns true if wedge A intersects wedge B.  Equivalent to but faster
    // than GetWedgeRelation() != WEDGE_IS_DISJOINT.
    // REQUIRES: A and B are non-empty.
    fun wedgeIntersects(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        // For A not to intersect B (where each loop interior is defined to be
        // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
        // Note that it's important to write these conditions as negatives
        // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
        // results when two vertices are the same.
        return !(S2Predicates.orderedCCW(a0, b2, b0, ab1) && S2Predicates.orderedCCW(b0, a2, a0, ab1))
    }

}
