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
package dilivia.s2.index.shape

import dilivia.s2.S2TextParser
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2CountEdgesUnitTest {

    @Test
    fun countEdgesUpToStopsEarly() {
        val index = S2TextParser.makeIndex(
            "0:0 | 0:1 | 0:2 | 0:3 | 0:4 # 1:0, 1:1 | 1:2, 1:3 | 1:4, 1:5, 1:6 #"
        )
        // Verify the test parameters.
        assertThat(index.nextNewShapeId()).isEqualTo(4)
        assertThat(index.shape(0)?.numEdges).isEqualTo(5)
        assertThat(index.shape(1)?.numEdges).isEqualTo(1)
        assertThat(index.shape(2)?.numEdges).isEqualTo(1)
        assertThat(index.shape(3)?.numEdges).isEqualTo(2)

        assertThat(S2CountEdges.countEdges(index)).isEqualTo(9)
        assertThat(S2CountEdges.countEdgesUpTo(index, 1)).isEqualTo(5)
        assertThat(S2CountEdges.countEdgesUpTo(index, 5)).isEqualTo(5)
        assertThat(S2CountEdges.countEdgesUpTo(index, 6)).isEqualTo(6)
        assertThat(S2CountEdges.countEdgesUpTo(index, 8)).isEqualTo(9)
    }

}
