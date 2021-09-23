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
package dilivia.s2.shape

import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2LaxPolylineShapeUnitTest {

    @Test
    fun noVertices() {
        val vertices: List<S2Point> = emptyList()
        val shape = S2LaxPolylineShape(vertices)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun oneVertex() {
        val vertices = listOf(S2Point(1, 0, 0))
        val shape = S2LaxPolylineShape(vertices)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
    }

    @Test
    fun edgeAccess() {
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1")
        val shape = S2LaxPolylineShape(vertices)
        assertThat(shape.numEdges).isEqualTo(2)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.chain(0).start).isEqualTo(0)
        assertThat(shape.chain(0).length).isEqualTo(2)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        val edge0 = shape.edge(0)
        assertThat(edge0.v0).isEqualTo(vertices[0])
        assertThat(edge0.v1).isEqualTo(vertices[1])
        val edge1 = shape.edge(1)
        assertThat(edge1.v0).isEqualTo(vertices[1])
        assertThat(edge1.v1).isEqualTo(vertices[2])
    }

}
