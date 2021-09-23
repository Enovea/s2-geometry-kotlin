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
package dilivia.collections

fun <A, B> Pair<A, B>.reverse(): Pair<B, A> = Pair(this.second, this.first)

operator fun <A> Pair<A, A>.get(i: Int): A = when(i) {
    0 -> this.first
    1 -> this.second
    else -> throw IllegalArgumentException("i out of range 0..1")
}

