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
package dilivia.s2.index

import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class SequenceLexiconUnitTest {

    @Test
    fun int64() {
        val lex = SequenceLexicon<Long>()
        assertThat( lex.add(listOf())).isEqualTo(0)
        assertThat( lex.add(listOf( 5L ))).isEqualTo(1)
        assertThat( lex.add(listOf())).isEqualTo(0)
        assertThat(lex.add(listOf( 5L, 5L ))).isEqualTo(2)
        assertThat(lex.add(listOf( 5L, 0L, -3L ))).isEqualTo(3)
        assertThat( lex.add(listOf( 5L ))).isEqualTo(1)
        assertThat( lex.add(listOf( 0x7fffffffffffffffL ))).isEqualTo(4)
        assertThat(lex.add(listOf( 5, 0, -3 ))).isEqualTo(3)
        assertThat( lex.add(listOf())).isEqualTo(0)
        assertThat( lex.size()).isEqualTo(5)
        expectSequence(listOf(), lex.sequence(0))
        expectSequence(listOf( 5L ), lex.sequence(1))
        expectSequence(listOf( 5L, 5L ), lex.sequence(2))
        expectSequence(listOf( 5L, 0L, -3L ), lex.sequence(3))
        expectSequence(listOf( 0x7fffffffffffffffL ), lex.sequence(4))
    }

    fun testClear() {
        val lex = SequenceLexicon<Long>()
        assertThat( lex.add(listOf( 1L ))).isEqualTo(0)
        assertThat( lex.add(listOf( 2L ))).isEqualTo(1)
        lex.clear()
        assertThat( lex.add(listOf( 2L ))).isEqualTo(0)
        assertThat( lex.add(listOf( 1L ))).isEqualTo(1)
    }

    companion object {
        fun <T> expectSequence(expected: List<T>, actual: SequenceLexicon<T>.Sequence) {
            assertThat( actual.size()).isEqualTo(expected.size)
            assertThat( actual.values()).isEqualTo(expected)
        }
    }

}

