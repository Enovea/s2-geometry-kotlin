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
package dilivia.s2.builder.layers

import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.S2TextParser.makePoint
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.S2Builder
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.shape.S2PointVectorShape
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test

class S2PointVectorLayerUnitTest {

    private fun verifyS2PointVectorLayerResults(
        label_set_ids: S2PointVectorLabelSetIds,
        label_set_lexicon: IdSetLexicon,
        output: MutableList<S2Point>,
        str_expected_points: String,
        expected_labels: List<List<Int>>
    ) {
        val expectedPoints = S2TextParser.parsePoints(str_expected_points)
        assertThat(label_set_ids.size).isEqualTo(expected_labels.size)
        for (i in 0 until output.size) {
            assertThat(output[i]).isEqualTo(expectedPoints[i])
            assertThat(label_set_lexicon.idSet(label_set_ids[i]).size()).isEqualTo(expected_labels[i].size)
            var k = 0
            for (label in label_set_lexicon.idSet(label_set_ids[i])) {
                assertThat(label).isEqualTo(expected_labels[i][k++])
            }
        }
    }

    fun addPoint(p: S2Point, builder: S2Builder) {
        builder.addEdge(p, p)
    }

    @Test
    fun mergeDuplicates() {
        val builder = S2Builder(S2Builder.Options())
        val output = mutableListOf<S2Point>()
        val labelSetLexicon = IdSetLexicon()
        val labelSetIds: S2PointVectorLabelSetIds = mutableListOf()
        builder.startLayer(
            S2PointVectorLayer(
                output,
                labelSetIds,
                labelSetLexicon,
                S2PointVectorLayer.Options(DuplicateEdges.MERGE)
            )
        )
        builder.setLabel(1)
        addPoint(makePoint("0:1"), builder)
        addPoint(makePoint("0:2"), builder)
        builder.setLabel(2);
        addPoint(makePoint("0:1"), builder)
        addPoint(makePoint("0:4"), builder)
        addPoint(makePoint("0:5"), builder)
        builder.clearLabels()
        addPoint(makePoint("0:5"), builder)
        addPoint(makePoint("0:6"), builder)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        val expected_labels = listOf(listOf(1, 2), listOf(1), listOf(2), listOf(2), emptyList())
        val expected_points = "0:1, 0:2, 0:4, 0:5, 0:6"
        verifyS2PointVectorLayerResults(labelSetIds, labelSetLexicon, output, expected_points, expected_labels)
    }

    @Test
    fun keepDuplicates() {
        val builder = S2Builder(S2Builder.Options())
        val output = mutableListOf<S2Point>()
        val labelSetLexicon = IdSetLexicon()
        val labelSetIds: S2PointVectorLabelSetIds = mutableListOf()
        builder.startLayer(
            S2PointVectorLayer(
                output,
                labelSetIds,
                labelSetLexicon,
                S2PointVectorLayer.Options(DuplicateEdges.KEEP)
            )
        )

        builder.setLabel(1);
        addPoint(makePoint("0:1"), builder)
        addPoint(makePoint("0:2"), builder)
        builder.setLabel(2);
        addPoint(makePoint("0:1"), builder)
        addPoint(makePoint("0:4"), builder)
        addPoint(makePoint("0:5"), builder)
        builder.clearLabels()
        addPoint(makePoint("0:5"), builder)
        addPoint(makePoint("0:6"), builder)
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()

        val expected_labels = listOf(listOf(1), listOf(2), listOf(1), listOf(2), listOf(2), emptyList(), emptyList())
        val expected_points = "0:1, 0:1, 0:2, 0:4, 0:5, 0:5, 0:6"

        verifyS2PointVectorLayerResults(labelSetIds, labelSetLexicon, output, expected_points, expected_labels)
    }

    @Test
    fun error() {
        val builder = S2Builder(S2Builder.Options())
        val output = mutableListOf<S2Point>()
        builder.startLayer(S2PointVectorLayer(output, options = S2PointVectorLayer.Options(DuplicateEdges.KEEP)))

        addPoint(makePoint("0:1"), builder)
        builder.addEdge(makePoint("0:3"), makePoint("0:4"))
        addPoint(makePoint("0:5"), builder)
        val error = S2Error()
        assertThat(builder.build(error)).isFalse()
        assertThat(error.code).isEqualTo(S2Error.INVALID_ARGUMENT)
        assertThat(error.text).isEqualTo("Found non-degenerate edges")

        assertThat(2).isEqualTo(output.size)
        assertThat(output[0]).isEqualTo(makePoint("0:1"))
        assertThat(output[1]).isEqualTo(makePoint("0:5"))
    }

    @Test
    fun indexedS2PointVectorLayerAddsShapes() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PointVectorLayer(index))
        val point0_str = "0:0"
        val point1_str = "2:2"
        builder.addPoint(makePoint(point0_str))
        builder.addPoint(makePoint(point1_str))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(1)
        val shape = index.shape(0) as S2PointVectorShape
        assertThat(shape.numPoints()).isEqualTo(2);
        assertThat(S2TextParser.toString(shape.point(0))).isEqualTo(point0_str)
        assertThat(S2TextParser.toString(shape.point(1))).isEqualTo(point1_str)
    }

    @Test
    fun indexedS2PointVectorLayerAddsEmptyShape() {
        val builder = S2Builder(S2Builder.Options())
        val index = MutableS2ShapeIndex()
        builder.startLayer(IndexedS2PointVectorLayer(index))
        val error = S2Error()
        assertThat(builder.build(error)).isTrue()
        assertThat(index.nextNewShapeId()).isEqualTo(0)
    }

}
