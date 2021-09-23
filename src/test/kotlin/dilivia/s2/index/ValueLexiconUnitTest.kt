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

import dilivia.s2.S2Point
import dilivia.s2.index.lexicon.InMemoryValueLexicon
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class ValueLexiconUnitTest {

  @Test
  fun duplicateValues() {
    val lex = InMemoryValueLexicon<Long>()
    assertThat( lex.add(5L)).isEqualTo(0)
    assertThat( lex.add(0L)).isEqualTo(1)
    assertThat( lex.add(0L)).isEqualTo(1)
    assertThat( lex.add(-3L)).isEqualTo(2)
    assertThat( lex.add(5L)).isEqualTo(0)
    assertThat( lex.add(0L)).isEqualTo(1)
    assertThat( lex.add(0x7fffffffffffffffL)).isEqualTo(3)
    assertThat( lex.add(-0x80L shl (14*4))).isEqualTo(4)
    assertThat( lex.add(0x7fffffffffffffffL)).isEqualTo(3)
    assertThat( lex.add(-0x80L shl (14*4))).isEqualTo(4)
    assertThat( lex.size()).isEqualTo(5)
    assertThat( lex.value(0)).isEqualTo(5L)
    assertThat( lex.value(1)).isEqualTo(0L)
    assertThat( lex.value(2)).isEqualTo(-3L)
    assertThat( lex.value(3)).isEqualTo(0x7fffffffffffffff)
    assertThat( lex.value(4)).isEqualTo(-0x80L shl (14*4))
  }

  @Test
  fun clear() {
    val lex = InMemoryValueLexicon<Long>()
    assertThat( lex.add(1)).isEqualTo(0)
    assertThat( lex.add(2)).isEqualTo(1)
    assertThat( lex.add(1)).isEqualTo(0)
    lex.clear()
    assertThat( lex.add(2)).isEqualTo(0)
    assertThat( lex.add(1)).isEqualTo(1)
    assertThat( lex.add(2)).isEqualTo(0)
  }

  @Test
  fun floatEquality() {
    val lex = InMemoryValueLexicon<S2Point>()
    val a = S2Point (1.0, 0.0, 0.0)
    val b = S2Point (1.0, -0.0, 0.0)
    val c = S2Point(1.0, 0.0, -0.0)
    assertThat( lex.add(a)).isEqualTo(0)
    assertThat( lex.add(b)).isEqualTo(0)
    assertThat( lex.add(c)).isEqualTo(0)
    assertThat( lex.size()).isEqualTo(1)
  }

}
