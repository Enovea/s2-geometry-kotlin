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
package dilivia.s2.builder

import dilivia.PreConditions.checkEQ
import dilivia.PreConditions.checkGT
import dilivia.PreConditions.checkState
import dilivia.collections.assign
import dilivia.collections.assignWith
import dilivia.s2.S2PointUtil
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.shape.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.S2Shape.Companion.containsBruteForce
import mu.KotlinLogging

object S2PolygonBoundariesBuilder {

    private val logger = KotlinLogging.logger {  }

    // The purpose of this function is to construct polygons consisting of
    // multiple loops.  It takes as input a collection of loops whose boundaries
    // do not cross, and groups them into polygons whose interiors do not
    // intersect (where the boundary of each polygon may consist of multiple
    // loops).
    //
    // some of those islands have lakes, then the input to this function would
    // islands, and their lakes.  Each loop would actually be present twice, once
    // in each direction (see below).  The output would consist of one polygon
    // representing each lake, one polygon representing each island not including
    // islands or their lakes, and one polygon representing the rest of the world
    //
    // This method is intended for internal use; external clients should use
    // S2Builder, which has more convenient interface.
    //
    // The input consists of a set of connected components, where each component
    // consists of one or more loops.  The components must satisfy the following
    // properties:
    //
    //  - The loops in each component must form a subdivision of the sphere (i.e.,
    //    they must cover the entire sphere without overlap), except that a
    //    component may consist of a single loop if and only if that loop is
    //    degenerate (i.e., its interior is empty).
    //
    //  - The boundaries of different components must be disjoint (i.e. no
    //    crossing edges or shared vertices).
    //
    //  - No component should be empty, and no loop should have zero edges.
    //
    // The output consists of a set of polygons, where each polygon is defined by
    // the collection of loops that form its boundary.  This function does not
    // actually construct any S2Shapes; it simply identifies the loops that belong
    // to each polygon.
    fun buildPolygonBoundaries(components: List<List<S2Shape>>,
                               polygons: MutableList<MutableList<S2Shape>>) {
        polygons.clear()
        if (components.isEmpty()) return

        // Since the loop boundaries do not cross, a loop nesting hierarchy can be
        // defined by choosing any point on the sphere as the "point at infinity".
        // Loop A then contains loop B if (1) A contains the boundary of B and (2)
        // loop A does not contain the point at infinity.
        //
        // We choose S2::Origin() for this purpose.  The loop nesting hierarchy then
        // determines the face structure.  Here are the details:
        //
        // 1. Build an S2ShapeIndex of all loops that do not contain S2::Origin().
        //    This leaves at most one unindexed loop per connected component
        //    (the "outer loop").
        //
        // 2. For each component, choose a representative vertex and determine
        //    which indexed loops contain it.  The "depth" of this component is
        //    defined as the number of such loops.
        //
        // 3. Assign the outer loop of each component to the containing loop whose
        //    depth is one less.  This generates a set of multi-loop polygons.
        //
        // 4. The outer loops of all components at depth 0 become a single face.

        val index = MutableS2ShapeIndex()
        // A map from shape.id() to the corresponding component number.
        val component_ids = mutableListOf<Int>()
        val outer_loops = mutableListOf<S2Shape>()
        var outerLoopId = 0
        for (i in components.indices) {
            val component = components[i]
            for (loop in component) {
                if (component.size > 1 && !containsBruteForce(loop, S2PointUtil.origin())) {
                    index.add(loop);
                    component_ids.add(i)
                } else {
                    loop.id = --outerLoopId
                    outer_loops.add(loop)
                }
            }
            // Check that there is exactly one outer loop in each component.
            checkEQ(i + 1, outer_loops.size) { "Component is not a subdivision" }
        }

        logger.trace { "Outer loops: $outer_loops\nIndex\n${index.toDebugString()}" }

        // Find the loops containing each component.
        val ancestors = ArrayList<List<S2Shape>>(components.size)
        ancestors.assign(components.size, emptyList())
        val contains_query = makeS2ContainsPointQuery(index)
        for (i in 0 until outer_loops.size) {
            val loop = outer_loops[i]
            checkGT(loop.numEdges, 0)
            ancestors[i] = contains_query.getContainingShapes(loop.edge(0).v0)
        }

        logger.trace { """
            |
            |${ancestors.mapIndexed { index, shapes -> "Outer Loop $index -> ${shapes.map { s -> s.id }}" }.joinToString("\n")}
        """.trimMargin() }

        // Assign each outer loop to the component whose depth is one less.
        // Components at depth 0 become a single face.
        val children = mutableMapOf<S2Shape?, MutableList<S2Shape>>()
        for (i in 0 until outer_loops.size) {
            var ancestor: S2Shape? = null
            val depth = ancestors[i].size
            if (depth > 0) {
                for (candidate in ancestors[i]) {
                    if (ancestors[component_ids[candidate.id]].size == depth - 1) {
                        checkState { ancestor == null }
                        ancestor = candidate
                    }
                }
                check(ancestor != null)
            }
            if (!children.containsKey(ancestor)) children[ancestor] = mutableListOf()
            children.getValue(ancestor).add(outer_loops[i])
        }

        logger.trace { """
            |Children:
            |------------------------------
            |${children.entries
                .sortedBy { entry -> entry.key?.id }
                .joinToString("\n") { entry -> "${entry.key?.id} => ${entry.value.map { s -> "${s.id};${s.edge(0).v0.x }" }}"}}
            |------------------------------
        """.trimMargin() }

        // There is one face per loop that is not an outer loop, plus one for the
        // outer loops of components at depth 0.
        polygons.assignWith(index.nextNewShapeId() + 1) { mutableListOf() }
        for (i in 0 until index.nextNewShapeId()) {
            var polygon = polygons[i]
            val loop = index.shape(i)!!
            val loopChildren = children[loop]
            if (loopChildren != null) {
                polygons[i] = loopChildren
                polygon = polygons[i]
            }
            polygon.add(loop)
        }
        if (children.containsKey(null)) polygons[polygons.lastIndex] = children.getValue(null)

        // Explicitly release the shapes from the index so they are not deleted.
        index.removeAll()
    }
}
