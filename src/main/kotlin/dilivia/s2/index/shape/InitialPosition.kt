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

/**
 * When passed to an Iterator constructor, specifies whether the iterator should be positioned :
 * - at the beginning of the index (BEGIN),
 * - the end of the index (END),
 * - or arbitrarily (UNPOSITIONED).
 *
 * By default dilivia.iterators are unpositioned, since this avoids an extra seek in this situation where one of the seek
 * methods (such as Locate) is immediately called.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
enum class InitialPosition { BEGIN, END, UNPOSITIONED }
