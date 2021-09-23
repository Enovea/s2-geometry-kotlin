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

import org.assertj.core.api.Assertions
import org.junit.jupiter.api.Test
import java.io.ByteArrayInputStream
import java.io.ByteArrayOutputStream
import java.io.DataInputStream
import java.io.DataOutputStream

internal class DataIOTest {

    @Test
    fun readUnsignedLong() {

        repeat(1000) { i ->
            Random.reset(i)

            val value = Random.randomULong()
            val bytes = ByteArrayOutputStream().use { output ->
                DataOutputStream(output).use { dataOutput ->
                    DataIO.writeUInt64(dataOutput, value)
                    output.toByteArray()
                }
            }

            ByteArrayInputStream(bytes).use { input -> DataInputStream(input).use { dataInput ->
                Assertions.assertThat(DataIO.readUInt64(dataInput)).isEqualTo(value)
            } }
        }

    }
}
