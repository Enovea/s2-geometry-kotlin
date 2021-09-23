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
package dilivia.s2.coords

import dilivia.s2.S2Point


// The XYZ and face,si,ti coordinates of an S2Point and, if this point is equal
// to the center of an S2Cell, the level of this cell (-1 otherwise).
data class S2XYZFaceSiTi(val xyz: S2Point, val face: Int, val si: UInt, val ti: UInt, val cellLevel: Int)
