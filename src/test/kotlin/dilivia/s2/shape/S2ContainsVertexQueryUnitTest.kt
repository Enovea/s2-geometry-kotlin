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

import dilivia.s2.S1Angle
import dilivia.s2.S2TextParser
import dilivia.s2.region.S2Loop
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2ContainsVertexQueryUnitTest {

    @Test
    fun undetermined() {
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("1:2"))
        q.addEdge(S2TextParser.makePoint("3:4"), 1)
        q.addEdge(S2TextParser.makePoint("3:4"), -1)
        assertThat(q.containsSign()).isEqualTo(0)
    }

    fun testContainedWithDuplicates() {
        // The S2::Ortho reference direction points approximately due west.
        // Containment is determined by the unmatched edge immediately clockwise.
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("0:0"))
        q.addEdge(S2TextParser.makePoint("3:-3"), -1)
        q.addEdge(S2TextParser.makePoint("1:-5"), 1)
        q.addEdge(S2TextParser.makePoint("2:-4"), 1)
        q.addEdge(S2TextParser.makePoint("1:-5"), -1)
        assertThat(q.containsSign()).isEqualTo(1)
    }

    fun testNotContainedWithDuplicates() {
        // The S2::Ortho reference direction points approximately due west.
        // Containment is determined by the unmatched edge immediately clockwise.
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("1:1"))
        q.addEdge(S2TextParser.makePoint("1:-5"), 1)
        q.addEdge(S2TextParser.makePoint("2:-4"), -1)
        q.addEdge(S2TextParser.makePoint("3:-3"), 1)
        q.addEdge(S2TextParser.makePoint("1:-5"), -1)
        assertThat(q.containsSign()).isEqualTo(-1)
    }

    fun testMatchesLoopContainment() {
        // Check that the containment function defined is compatible with S2Loop
        // (which at least currently does not use this class).
        val loop = S2Loop.makeRegularLoop(S2TextParser.makePoint("89:-179"), S1Angle.degrees(10), 1000)
        for (i in 1..loop.numVertices) {
            val q = S2ContainsVertexQuery(loop.vertex(i))
            q.addEdge(loop.vertex(i - 1), -1)
            q.addEdge(loop.vertex(i + 1), 1)
            assertThat(loop.contains(loop.vertex(i))).isEqualTo(q.containsSign() > 0)
        }
    }
}
  
