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

import dilivia.ComparisonChain
import dilivia.s2.S2Error
import dilivia.s2.S2TextParser
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.builder.layers.Layer
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2FindPolygonDegeneraciesUnitTest {

    data class TestDegeneracy(val edge_str: String, val is_hole: Boolean) : Comparable<TestDegeneracy> {

        override fun compareTo(other: TestDegeneracy): Int = ComparisonChain.start()
            .compare(edge_str, other.edge_str)
            .compare(is_hole, other.is_hole)
            .result()

    }

    class DegeneracyCheckingLayer(val expected: List<TestDegeneracy>) : Layer() {

        override fun graphOptions(): GraphOptions = GraphOptions(
            EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS, DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS
        )

        override fun build(g: Graph, error: S2Error) {
            val degeneracies = S2FindPolygonDegeneracies.findPolygonDegeneracies(g, error)
            // Convert the output into a human-readable format.
            val actual = mutableListOf<TestDegeneracy>()
            for (degeneracy in degeneracies) {
                val edge = g.edge(degeneracy.edgeId)
                val points = listOf(g.vertex(edge.first), g.vertex(edge.second))
                actual.add(TestDegeneracy(S2TextParser.toString(points), degeneracy.isHole));
            }
            actual.sort()
            assertThat(actual).containsExactly(*expected.sorted().toTypedArray())
            assertThat(S2FindPolygonDegeneracies.isFullyDegenerate(g)).isEqualTo(degeneracies.size == g.numEdges)
        }

        override fun toString(): String {
            return expected.joinToString(" ")
        }

    }

    fun expectDegeneracies(polygon_str: String, expected: List<TestDegeneracy>) {
        val builder = S2Builder(S2Builder.Options())
        builder.startLayer(DegeneracyCheckingLayer(expected))
        val polygon = S2TextParser.makeLaxPolygon(polygon_str)
        builder.addIsFullPolygonPredicate { g, error -> polygon.getReferencePoint().contained }

        for (i in 0 until polygon.numEdges) {
            val edge = polygon.edge(i)
            builder.addEdge(edge.v0, edge.v1)
        }
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
    }

    @Test
    fun EmptyPolygon() {
        expectDegeneracies("", emptyList());
    }

    @Test
    fun NoDegeneracies() {
        expectDegeneracies("0:0, 0:1, 1:0", emptyList());
    }

    @Test
    fun PointShell() {
        expectDegeneracies("0:0", listOf(TestDegeneracy("0:0, 0:0", false)))
    }

    @Test
    fun SiblingPairShells() {
        expectDegeneracies(
            "0:0, 0:1, 1:0; 1:0, 0:1, 0:0",
            listOf(
                TestDegeneracy("0:0, 0:1", false),
                TestDegeneracy("0:1, 0:0", false),
                TestDegeneracy("0:1, 1:0", false),
                TestDegeneracy("1:0, 0:1", false),
                TestDegeneracy("0:0, 1:0", false),
                TestDegeneracy("1:0, 0:0", false)
            )
        )
    }

    @Test
    fun AttachedSiblingPairShells() {
        expectDegeneracies(
            "0:0, 0:1, 1:0; 1:0, 2:0",
            listOf(
                TestDegeneracy("1:0, 2:0", false),
                TestDegeneracy("2:0, 1:0", false)
            )
        )
    }

    @Test
    fun AttachedSiblingPairHoles() {
        expectDegeneracies(
            "0:0, 0:3, 3:0; 0:0, 1:1",
            listOf(
                TestDegeneracy("0:0, 1:1", true),
                TestDegeneracy("1:1, 0:0", true)
            )
        )
    }

    @Test
    fun AttachedSiblingPairShellsAndHoles() {
        expectDegeneracies(
            "0:0, 0:3, 3:0; 3:0, 1:1; 3:0, 5:5",
            listOf(
                TestDegeneracy("3:0, 1:1", true),
                TestDegeneracy("1:1, 3:0", true),
                TestDegeneracy("3:0, 5:5", false),
                TestDegeneracy("5:5, 3:0", false),
            )
        )
    }

    @Test
    fun DegenerateShellsOutsideLoop() {
        expectDegeneracies(
            "0:0, 0:3, 3:3, 3:0; 4:4, 5:5; 6:6",
            listOf(
                TestDegeneracy("4:4, 5:5", false),
                TestDegeneracy("5:5, 4:4", false),
                TestDegeneracy("6:6, 6:6", false)
            )
        )
    }

    @Test
    fun DegenerateHolesWithinLoop() {
        expectDegeneracies(
            "0:0, 0:5, 5:5, 5:0; 1:1, 2:2; 3:3",
            listOf(
                TestDegeneracy("1:1, 2:2", true),
                TestDegeneracy("2:2, 1:1", true),
                TestDegeneracy("3:3, 3:3", true)
            )
        )
    }

    @Test
    fun PointHoleWithinFull() {
        expectDegeneracies("full; 0:0", listOf(TestDegeneracy("0:0, 0:0", true)))
    }

    @Test
    fun SiblingPairHolesWithinFull() {
        expectDegeneracies(
            "full; 0:0, 0:1, 1:0; 1:0, 0:1, 0:0",
            listOf(
                TestDegeneracy("0:0, 0:1", true), TestDegeneracy("0:1, 0:0", true),
                TestDegeneracy("0:1, 1:0", true),
                TestDegeneracy("1:0, 0:1", true), TestDegeneracy("0:0, 1:0", true), TestDegeneracy("1:0, 0:0", true)
            )
        )
    }

}
