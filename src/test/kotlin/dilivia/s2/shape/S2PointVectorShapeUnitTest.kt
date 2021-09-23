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

class S2PointVectorShapeUnitTest {

    @Test
    fun empty() {
        val points = emptyList<S2Point>()
        val shape = S2PointVectorShape(points = points)
        assertThat(shape.numEdges).isEqualTo(0)
        assertThat(shape.numChains).isEqualTo(0)
        assertThat(shape.dimension).isEqualTo(0)
        assertThat(shape.isEmpty()).isTrue()
        assertThat(shape.isFull()).isFalse()
        assertThat(shape.getReferencePoint().contained).isFalse()
    }

    @Test
    fun constructionAndAccess() {
        val points = mutableListOf<S2Point>()
        S2Random.reset(1)
        val kNumPoints = 100
        repeat(kNumPoints) {
            points.add(S2Random.randomPoint())
        }
        val shape = S2PointVectorShape(points = points)

        assertThat(shape.numEdges).isEqualTo(kNumPoints)
        assertThat(shape.numChains).isEqualTo(kNumPoints)
        assertThat(shape.dimension).isEqualTo(0)
        assertThat(shape.isEmpty()).isFalse()
        assertThat(shape.isFull()).isFalse()
        S2Random.reset(1)
        repeat(kNumPoints) { i ->
            assertThat(shape.chain(i).start).isEqualTo(i)
            assertThat(shape.chain(i).length).isEqualTo(1)
            val edge = shape.edge(i)
            val pt = S2Random.randomPoint()
            assertThat(edge.v0).isEqualTo(pt)
            assertThat(edge.v1).isEqualTo(pt)
        }
    }

}
