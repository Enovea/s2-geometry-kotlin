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

import java.io.ByteArrayInputStream
import java.io.ByteArrayOutputStream
import java.io.DataInput
import java.io.DataInputStream
import java.io.DataOutput
import java.io.DataOutputStream
import java.io.IOException
import java.lang.Byte.toUnsignedInt

object DataIO {

    fun readUInt64(bytes: ByteArray, offset: Int = 0): ULong {
        require(offset + 8 <= bytes.size)
        return (toUnsignedInt(bytes[offset + 0]).toULong() shl 56) +
                (toUnsignedInt(bytes[offset + 1]).toULong() shl 48) +
                (toUnsignedInt(bytes[offset + 2]).toULong() shl 40) +
                (toUnsignedInt(bytes[offset + 3]).toULong() shl 32) +
                (toUnsignedInt(bytes[offset + 4]).toULong() shl 24) +
                (toUnsignedInt(bytes[offset + 5]).toULong() shl 16) +
                (toUnsignedInt(bytes[offset + 6]).toULong() shl 8) +
                (toUnsignedInt(bytes[offset + 7]).toULong() shl 0)
    }

    fun writeUInt64(bytes: ByteArray, value: ULong, offset: Int = 0) {
        bytes[offset + 0] = (value shr 56).toByte()
        bytes[offset + 1] = (value shr 48).toByte()
        bytes[offset + 2] = (value shr 40).toByte()
        bytes[offset + 3] = (value shr 32).toByte()
        bytes[offset + 4] = (value shr 24).toByte()
        bytes[offset + 5] = (value shr 16).toByte()
        bytes[offset + 6] = (value shr 8).toByte()
        bytes[offset + 7] = (value shr 0).toByte()
    }

    fun readUInt64(input: DataInput): ULong {
        val bytes = ByteArray(8)
        input.readFully(bytes)
        return readUInt64(bytes)
    }

    fun writeUInt64(output: DataOutput, value: ULong) {
        output.writeByte((value shr 56).toInt())
        output.writeByte((value shr 48).toInt())
        output.writeByte((value shr 40).toInt())
        output.writeByte((value shr 32).toInt())
        output.writeByte((value shr 24).toInt())
        output.writeByte((value shr 16).toInt())
        output.writeByte((value shr 8).toInt())
        output.writeByte((value shr 0).toInt())
    }

    fun readInt32(bytes: ByteArray, offset: Int = 0): Int {
        require(offset + 4 <= bytes.size)
        return (bytes[0].toInt() shl 24) +
                (bytes[1].toInt() shl 16) +
                (bytes[2].toInt() shl 8) +
                (bytes[3].toInt() shl 0)
    }

    fun readUInt32(bytes: ByteArray, offset: Int = 0): UInt {
        require(offset + 4 <= bytes.size)
        return (toUnsignedInt(bytes[offset + 0]).toUInt() shl 24) +
                (toUnsignedInt(bytes[offset + 1]).toUInt() shl 16) +
                (toUnsignedInt(bytes[offset + 2]).toUInt() shl 8) +
                (toUnsignedInt(bytes[offset + 3]).toUInt() shl 0)
    }

    fun writeInt32(bytes: ByteArray, value: Int, offset: Int = 0) {
        require(offset + 4 <= bytes.size)
        bytes[offset + 0] = (value shr 24).toByte()
        bytes[offset + 1] = (value shr 16).toByte()
        bytes[offset + 2] = (value shr 8).toByte()
        bytes[offset + 3] = (value shr 0).toByte()
    }

    fun writeUInt32(bytes: ByteArray, value: UInt, offset: Int = 0) {
        require(offset + 4 <= bytes.size)
        bytes[offset + 0] = (value shr 24).toByte()
        bytes[offset + 1] = (value shr 16).toByte()
        bytes[offset + 2] = (value shr 8).toByte()
        bytes[offset + 3] = (value shr 0).toByte()
    }

    fun readUInt32(input: DataInput): UInt {
        val bytes = ByteArray(4)
        input.readFully(bytes)
        return  readUInt32(bytes)
    }

    fun writeUInt32(output: DataOutput, value: UInt) {
        output.writeByte((value shr 24).toInt())
        output.writeByte((value shr 16).toInt())
        output.writeByte((value shr 8).toInt())
        output.writeByte((value shr 0).toInt())
    }

    fun readInt24(bytes: ByteArray, offset: Int = 0): Int {
        require(offset + 3 <= bytes.size)
        return (bytes[offset + 0].toInt() shl 16) +
                (bytes[offset + 1].toInt() shl 8) +
                (bytes[offset + 2].toInt() shl 0)
    }

    fun writeInt24(bytes: ByteArray, value: Int, offset: Int = 0) {
        require(offset + 3 <= bytes.size)
        bytes[offset + 0] = (value shr 16).toByte()
        bytes[offset + 1] = (value shr 8).toByte()
        bytes[offset + 2] = (value shr 0).toByte()
    }

    fun readUInt24(bytes: ByteArray, offset: Int = 0): UInt {
        require(offset + 3 <= bytes.size)
        return (toUnsignedInt(bytes[offset + 0]).toUInt() shl 16) +
                (toUnsignedInt(bytes[offset + 1]).toUInt() shl 8) +
                (toUnsignedInt(bytes[offset + 2]).toUInt() shl 0)
    }

    fun writeUInt24(bytes: ByteArray, value: UInt, offset: Int = 0) {
        require(offset + 3 <= bytes.size)
        bytes[offset + 0] = (value shr 16).toByte()
        bytes[offset + 1] = (value shr 8).toByte()
        bytes[offset + 2] = (value shr 0).toByte()
    }

    fun readUInt24(input: DataInput): UInt {
        val bytes = ByteArray(3)
        input.readFully(bytes)
        return readUInt24(bytes)
    }

    fun writeUInt24(output: DataOutput, value: UInt) {
        output.writeByte((value shr 16).toInt())
        output.writeByte((value shr 8).toInt())
        output.writeByte((value shr 0).toInt())
    }

    fun readUInt16(bytes: ByteArray, offset: Int = 0): UShort {
        require(offset + 2 <= bytes.size)
        return ((toUnsignedInt(bytes[offset]) shl 8) +
                (toUnsignedInt(bytes[offset + 1]) shl 0)).toUShort()
    }

    fun readUInt16(input: DataInput): UShort {
        val bytes = ByteArray(2)
        input.readFully(bytes)
        return  ((toUnsignedInt(bytes[0]).toUInt() shl 8) +
                (toUnsignedInt(bytes[1]).toUInt() shl 0)).toUShort()
    }

    fun writeUInt16(output: DataOutput, value: UShort) {
        output.writeByte((value.toUInt() shr 8).toInt())
        output.writeByte((value.toUInt() shr 0).toInt())
    }

    fun writeUInt16(bytes: ByteArray, value: UShort, offset: Int = 0) {
        bytes[offset + 0] = (value.toInt() shr 8).toByte()
        bytes[offset + 1] = (value.toInt() shr 0).toByte()
    }

    /**
     * Unpack int value from the input stream.
     *
     * @param `in` The input stream.
     * @return The long value.
     * a
     * @throws java.io.IOException in case of IO error
     */
    @Throws(IOException::class)
    fun unpackInt(input: DataInput): Int {
        var ret = 0
        var v: Int
        do {
            v = input.readByte().toInt()
            ret = (ret shl 7) or (v and 0x7F)
        } while (v and 0x80 == 0)
        return ret
    }

    /**
     * Pack int into an output stream.
     * It will occupy 1-5 bytes depending on value (lower values occupy smaller space)
     *
     * @param out DataOutput to put value into
     * @param value to be serialized, must be non-negative
     * @throws java.io.IOException in case of IO error
     */
    @Throws(IOException::class)
    fun packInt(out: DataOutput, value: Int) {
        // Optimize for the common case where value is small. This is particular important where our caller
        // is SerializerBase.SER_STRING.serialize because most chars will be ASCII characters and hence in this range.
        // credit Max Bolingbroke https://github.com/jankotek/MapDB/pull/489
        var shift = value and 0x7F.inv() //reuse variable
        if (shift != 0) {
            shift = 31 - Integer.numberOfLeadingZeros(value)
            shift -= shift % 7 // round down to nearest multiple of 7
            while (shift != 0) {
                out.writeByte((value ushr shift and 0x7F).toInt())
                //$DELAY$
                shift -= 7
            }
        }
        out.writeByte((value and 0x7F or 0x80).toInt())
    }

    @Throws(IOException::class)
    fun packDeltaIntList(out: DataOutput, value: List<Int>) {
        packInt(out, value.size)
        if (value.isNotEmpty()) {
            var prev = value.first()
            packInt(out, prev)
            for (i in 1..value.lastIndex) {
                val curr = value[i]
                packInt(out, curr - prev)
                prev = curr
            }
        }
    }

    @Throws(IOException::class)
    fun packDeltaLongList(out: DataOutput, value: List<Long>) {
        packInt(out, value.size)
        if (value.isNotEmpty()) {
            var prev = value.first()
            packLong(out, prev)
            for (i in 1..value.lastIndex) {
                val curr = value[i]
                packLong(out, curr - prev)
                prev = curr
            }
        }
    }

    @Throws(IOException::class)
    fun packDeltaLongList(value: List<Long>): ByteArray {
        return ByteArrayOutputStream().use { bytes ->
            DataOutputStream(bytes).use { output ->
                packDeltaLongList(output, value)
            }
            bytes.toByteArray()
        }
    }

    fun unpackDeltaIntList(input: DataInput): List<Int> {
        val size = unpackInt(input)
        val ret = ArrayList<Int>(size)
        var prev = 0
        for (i in 0 until size) {
            prev += unpackInt(input)
            ret.add(prev)
        }
        return ret
    }

    fun unpackDeltaLongList(input: DataInput): List<Long> {
        val size = unpackInt(input)
        val ret = ArrayList<Long>(size)
        var prev = 0L
        for (i in 0 until size) {
            prev += unpackLong(input)
            ret.add(prev)
        }
        return ret
    }

    fun unpackDeltaLongList(bytes: ByteArray): List<Long> {
        return ByteArrayInputStream(bytes).use { DataInputStream(it).use { input ->
            unpackDeltaLongList(input)
        } }
    }

    /**
     * Pack long into output.
     * It will occupy 1-10 bytes depending on value (lower values occupy smaller space)
     *
     * @param out DataOutput to put value into
     * @param value to be serialized, must be non-negative
     *
     * @throws java.io.IOException in case of IO error
     */
    @Throws(IOException::class)
    fun packLong(out: DataOutput, value: Long) {
        var shift = 63 - java.lang.Long.numberOfLeadingZeros(value)
        shift -= shift % 7 // round down to nearest multiple of 7
        while (shift != 0) {
            out.writeByte((value ushr shift and 0x7F).toInt())
            shift -= 7
        }
        out.writeByte((value and 0x7F or 0x80).toInt())
    }

    fun packLong(value: Long): ByteArray {
        val result = ByteArray(packedLongSize(value))
        var idx = 0
        var shift = 63 - java.lang.Long.numberOfLeadingZeros(value)
        shift -= shift % 7 // round down to nearest multiple of 7
        while (shift != 0) {
            result[idx++] = (value ushr shift and 0x7F).toByte()
            shift -= 7
        }
        result[idx] = (value and 0x7F or 0x80).toByte()
        return result
    }

    @Throws(IOException::class)
    fun unpackLong(input: DataInput): Long {
        var ret: Long = 0
        var v: Int
        do {
            v = input.readByte().toInt()
            ret = (ret shl 7) or (v and 0x7F).toLong()
        } while (v and 0x80 == 0)
        return ret
    }

    /**
     * Calculate how much bytes packed long consumes.
     *
     * @param value to calculate
     * @return number of bytes used in packed form
     */
    fun packedLongSize(value: Long): Int {
        var shift = 63 - java.lang.Long.numberOfLeadingZeros(value)
        shift -= shift % 7 // round down to nearest multiple of 7
        var ret = 1
        while (shift != 0) {
            shift -= 7
            ret++
        }
        return ret
    }

    fun packedIntSize(value: Int): Int {
        var shift = 31 - java.lang.Integer.numberOfLeadingZeros(value)
        shift -= shift % 7 // round down to nearest multiple of 7
        var ret = 1
        while (shift != 0) {
            shift -= 7
            ret++
        }
        return ret
    }
}
