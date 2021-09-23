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

import dilivia.s2.S2PointUtil
import dilivia.s2.S2TextParser
import dilivia.s2.region.S2Loop
import dilivia.s2.shape.S2Shape.Companion.containsBruteForce
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2LaxLoopShapeUnitTest {

    @Test
    fun laxLoopShapeEmptyLoop() {
        // Test S2Loop constructor.
        val shape = S2LaxLoopShape()
        shape.init(S2Loop(S2Loop.kEmpty))
        assertThat(shape.numVertices()).isEqualTo(0)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun laxLoopShapeNonEmptyLoop() {
        // Test vector<S2Point> constructor.
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val shape = S2LaxLoopShape(vertices)
        assertThat(shape.numVertices()).isEqualTo(vertices.size)
        assertThat(shape.numEdges).isEqualTo(vertices.size)
        assertThat(shape.numChains).isEqualTo(1);
        assertThat(shape.chain(0).start).isEqualTo(0);
        assertThat(shape.chain(0).length).isEqualTo(vertices.size)
        for (i in 0..vertices.lastIndex) {
            assertThat(shape.vertex(i)).isEqualTo(vertices[i])
            val edge = shape.edge(i)
            assertThat(edge.v0).isEqualTo(vertices[i])
            assertThat(edge.v1).isEqualTo(vertices[(i + 1) % vertices.size])
        }
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun laxClosedPolylineShapeNoInterior() {
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val shape = S2LaxClosedPolylineShape(vertices)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun vertexIdLaxLoopShapeEmptyLoop() {
        val shape = S2VertexIdLaxLoopShape(intArrayOf(), emptyArray())
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numVertices()).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun vertexIdLaxLoopShapeInvertedLoop() {
        val vertex_array = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val vertex_ids = intArrayOf(0, 3, 2, 1)  // Inverted.
        val shape = S2VertexIdLaxLoopShape(vertex_ids, vertex_array.toTypedArray())
        assertThat(shape.numEdges).isEqualTo(4)
        assertThat(shape.numVertices()).isEqualTo(4)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.chain(0).start).isEqualTo(0)
        assertThat(shape.chain(0).length).isEqualTo(4)
        assertThat(shape . vertex (0)).isEqualTo(vertex_array [0])
        assertThat(shape . vertex (1)).isEqualTo(vertex_array [3])
        assertThat(shape . vertex (2)).isEqualTo(vertex_array [2])
        assertThat(shape . vertex (3)).isEqualTo(vertex_array [1])
        assertThat(shape.dimension).isEqualTo(2)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        assertThat(containsBruteForce(shape, S2PointUtil.origin())).isTrue()
    }

}
