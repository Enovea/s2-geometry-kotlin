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

import dilivia.s2.S1Angle
import dilivia.s2.S2Earth
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.index.shape.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Loop
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

typealias EdgeIdVector = MutableList<ShapeEdgeId>

class S2ContainsPointQueryTest {

    private val logger = KotlinLogging.logger { }

    @Test
    fun containsPointQueryVertexModelOpen() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.OPEN)
        val q = makeS2ContainsPointQuery(index, options)
        assertThat(q.contains(S2TextParser.makePoint("0:0"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("-1:1"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("1:1"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:2"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:3"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:5"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:7"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("2:6"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("10:10"))).isFalse()

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertThat(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6"))).isFalse()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5"))).isFalse()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7"))).isFalse()
    }

    @Test
    fun containsPointQueryVertexModelSemiOpen() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.SEMI_OPEN)
        val q = makeS2ContainsPointQuery(index, options)
        assertThat(q.contains(S2TextParser.makePoint("0:0"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("-1:1"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("1:1"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:2"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:5"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:7"))).isTrue()  // Contained vertex.
        assertThat(q.contains(S2TextParser.makePoint("2:6"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("10:10"))).isFalse()

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertThat(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6"))).isFalse()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5"))).isFalse()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7"))).isTrue()
    }

    @Test
    fun containsPointQueryVertexModelClosed() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.CLOSED)
        val q = makeS2ContainsPointQuery(index, options)
        assertThat(q.contains(S2TextParser.makePoint("0:0"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("-1:1"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("1:1"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("0:2"))).isFalse()
        assertThat(q.contains(S2TextParser.makePoint("0:5"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("0:7"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("2:6"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q.contains(S2TextParser.makePoint("10:10"))).isFalse()

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertThat(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6"))).isFalse()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6"))).isTrue()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5"))).isTrue()
        assertThat(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7"))).isTrue()
    }

    @Test
    fun containsPointQueryGetContainingShapes() {
        // Also tests shapeContains().
        val kNumVerticesPerLoop = 10
        val kMaxLoopRadius = S1Angle.radians(S2Earth.kmToRadians(10.0))
        val centerCap = S2Cap.fromCenterAngle(S2Random.randomPoint(), kMaxLoopRadius)

        logger.trace { "Test params: num vertices per loop = $kNumVerticesPerLoop, max loop radius = $kMaxLoopRadius, center cap = $centerCap" }

        val index = MutableS2ShapeIndex()
        repeat(100) { i ->
            val center = S2Random.samplePoint(centerCap)
            val radius = kMaxLoopRadius * S2Random.randomDouble()
            val loop = S2Loop.makeRegularLoop(center, radius, kNumVerticesPerLoop)
            logger.trace { "Create loop $i: center = $center, radius = ${radius.radians} => $loop" }
            index.add(S2Loop.Shape(loop = loop))
        }

        logger.trace { "Index:\n-------------------\n${index.toDebugString()}\n--------------------" }

        val query = makeS2ContainsPointQuery(index)
        repeat(100) {
            val p = S2Random.samplePoint(centerCap)
            val expected = mutableListOf<S2Shape>()
            for (shape in index) {
                val loop = (shape as S2Loop.Shape).loop
                val loopContainsPoint = loop.contains(p)
                val shapeContainsPoint = query.shapeContains(shape, p)
                if (loopContainsPoint) {
                    assertThat(shapeContainsPoint)
                        .withFailMessage("loop contains point $p but query result is false: id = ${shape.id}, loop = $loop").isTrue()
                    expected.add(shape)
                } else {
                    assertThat(shapeContainsPoint).isFalse()
                }
            }
            val actual = query.getContainingShapes(p)
            assertThat(actual).isEqualTo(expected)
        }
    }

    fun expectIncidentEdgeIds(expected: EdgeIdVector, index: MutableS2ShapeIndex, p: S2Point) {
        val actual: EdgeIdVector = mutableListOf()
        val q = makeS2ContainsPointQuery(index)
        assertThat(q.visitIncidentEdges(p) { e ->
            actual.add(e.id)
            true
        }).isTrue()
        assertThat(actual).isEqualTo(expected)
    }

    @Test
    fun containsPointQueryVisitIncidentEdges() {
        val index = S2TextParser.makeIndex("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2")
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(0, 0)), index, S2TextParser.makePoint("0:0"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(0, 1), ShapeEdgeId(1, 0)), index, S2TextParser.makePoint("1:1"))
        expectIncidentEdgeIds(
            mutableListOf(ShapeEdgeId(1, 0), ShapeEdgeId(2, 0), ShapeEdgeId(2, 2)),
            index,
            S2TextParser.makePoint("1:2")
        )
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(2, 0), ShapeEdgeId(2, 1)), index, S2TextParser.makePoint("1:3"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(2, 1), ShapeEdgeId(2, 2)), index, S2TextParser.makePoint("2:2"))
    }

}
