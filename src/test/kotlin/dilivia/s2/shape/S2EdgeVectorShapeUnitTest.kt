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
import dilivia.s2.S2Random
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2EdgeVectorShapeUnitTest {

    @Test
    fun empty() {
        val shape = S2EdgeVectorShape()
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun edgeAccess() {
        val shape = S2EdgeVectorShape()
        S2Random.reset(0)
        val kNumEdges = 100
        repeat(kNumEdges) {
            val a = S2Random.randomPoint()  // Control the evaluation order
            shape.add(a, S2Random.randomPoint())
        }
        assertThat(shape.numEdges).isEqualTo(kNumEdges)
        assertThat(shape.numChains).isEqualTo(kNumEdges)
        assertThat(shape.dimension).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        S2Random.reset(0)
        repeat(kNumEdges) { i ->
            assertThat(shape.chain(i).start).isEqualTo(i)
            assertThat(shape.chain(i).length).isEqualTo(1)
            val edge = shape.edge(i)
            assertThat(edge.v0).isEqualTo(S2Random.randomPoint())
            assertThat(edge.v1).isEqualTo(S2Random.randomPoint())
        }
    }

    @Test
    fun singletonConstructor() {
        val a = S2Point(1, 0, 0)
        val b = S2Point(0, 1, 0)
        val shape = S2EdgeVectorShape(a, b)
        assertThat(shape.numEdges).isEqualTo(1)
        assertThat(shape.numChains).isEqualTo(1)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        val edge = shape.edge(0)
        assertThat(edge.v0).isEqualTo(a)
        assertThat(edge.v1).isEqualTo(b)
    }
}
