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
import dilivia.s2.region.S2Loop

// S2LaxClosedPolylineShape is like S2LaxPolylineShape except that the last
// vertex is implicitly joined to the first.  It is also like S2LaxLoopShape
// except that it does not have an interior (which makes it more efficient to
// index).
class S2LaxClosedPolylineShape : S2LaxLoopShape {

    // Constructs an empty loop.
    constructor(): super()

    // Constructs an S2LaxLoopShape with the given vertices.
    constructor(vertices: List<S2Point>): super(vertices)

    // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
    constructor(loop: S2Loop): super(loop)

    override val dimension: Int = 1

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

}
