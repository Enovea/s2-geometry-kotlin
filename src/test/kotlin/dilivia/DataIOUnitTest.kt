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
package dilivia

import org.assertj.core.api.BDDAssertions
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.api.Test
import java.io.ByteArrayInputStream
import java.io.ByteArrayOutputStream
import java.io.DataInputStream
import java.io.DataOutputStream

internal class DataIOUnitTest {

    @DisplayName("Pack / Unpack int")
    @Test
    fun packAndUnpackInt() {

        repeat(10000) {
            val value = Random.randomInt()
            val output = ByteArrayOutputStream().use { bytesOutput ->
                DataOutputStream(bytesOutput).use { dataOutput ->
                    DataIO.packInt(dataOutput, value)
                }
                bytesOutput.toByteArray()
            }

            val packedLongSize = DataIO.packedIntSize(value)
            BDDAssertions.then(output)
                .withFailMessage("Packed size ${output.size} for value $value != expected $packedLongSize")
                .hasSize(packedLongSize)

            DataInputStream(ByteArrayInputStream(output)).use { dataInput ->
                val unpacked = DataIO.unpackInt(dataInput)
                BDDAssertions.then(unpacked)
                    .withFailMessage("Unpacked value $unpacked != initial value $value")
                    .isEqualTo(value)
            }
        }

    }

    @DisplayName("Pack / Unpack long")
    @Test
    fun packAndUnpackLong() {

        repeat(10000) {
            val value = Random.randomLong()
            val output = ByteArrayOutputStream().use { bytesOutput ->
                DataOutputStream(bytesOutput).use { dataOutput ->
                    DataIO.packLong(dataOutput, value)
                }
                bytesOutput.toByteArray()
            }

            val packedLongSize = DataIO.packedLongSize(value)
            BDDAssertions.then(output)
                .withFailMessage("Packed size ${output.size} for value $value != expected $packedLongSize")
                .hasSize(packedLongSize)

            DataInputStream(ByteArrayInputStream(output)).use { dataInput ->
                val unpacked = DataIO.unpackLong(dataInput)
                BDDAssertions.then(unpacked)
                    .withFailMessage("Unpacked value $unpacked != initial value $value")
                    .isEqualTo(value)
            }
        }

    }


    @DisplayName("Pack / Unpack delta long list")
    @Test
    fun packAndUnpackDeltaLongList() {

        repeat(10000) {

            val size = Random.randomInt(1, 50)
            val values = (0 until size).map { Random.randomLong(0, 4000000) }

            val output = ByteArrayOutputStream().use { bytesOutput ->
                DataOutputStream(bytesOutput).use { dataOutput ->
                    DataIO.packDeltaLongList(dataOutput, values)
                }
                bytesOutput.toByteArray()
            }

            DataInputStream(ByteArrayInputStream(output)).use { dataInput ->
                val unpacked = DataIO.unpackDeltaLongList(dataInput)
                BDDAssertions.then(unpacked)
                    .withFailMessage("Unpacked value $unpacked != initial value $values")
                    .isEqualTo(values)
            }
        }

    }


    @DisplayName("Pack / Unpack delta int list")
    @Test
    fun packAndUnpackDeltaIntList() {

        repeat(10000) {

            val size = Random.randomInt(1, 50)
            val values = (0 until size).map { Random.randomInt(-40000, 40000) }

            val output = ByteArrayOutputStream().use { bytesOutput ->
                DataOutputStream(bytesOutput).use { dataOutput ->
                    DataIO.packDeltaIntList(dataOutput, values)
                }
                bytesOutput.toByteArray()
            }

            DataInputStream(ByteArrayInputStream(output)).use { dataInput ->
                val unpacked = DataIO.unpackDeltaIntList(dataInput)
                BDDAssertions.then(unpacked)
                    .withFailMessage("Unpacked value $unpacked != initial value $values")
                    .isEqualTo(values)
            }
        }

    }
}
