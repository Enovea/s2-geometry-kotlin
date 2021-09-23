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
package dilivia.math.vectors

/**
 * Interface that defines vector of dimension 3.
 *
 * @param T The type of numbers
 * @param V The implementation class of the vector.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
interface R3Vector<T, V> : RVector<T, V> where T : Number, T : Comparable<T>, V : R3Vector<T, V> {

    /** The x component of the vector (indice 0). */
    var x: T
        set(value) {
            this[0] = value
        }
        get() = this[0]

    /** The y component of the vector (indice 1). */
    var y: T
        set(value) {
            this[1] = value
        }
        get() = this[1]

    /** The z component of the vector (indice 2). */
    var z: T
        set(value) {
            this[2] = value
        }
        get() = this[2]

    /**
     * Computes the cross product of this vector with another one.
     *
     * @param other A R2Vector.
     * @return The value of this x other.
     */
    fun crossProd(other: R3Vector<T, V>): V = newInstance { i ->
        when (i) {
            /*x =*/ 0 -> numberType.minus(numberType.times(y, other.z), numberType.times(z, other.y))
            /*y =*/ 1 -> numberType.minus(numberType.times(z, other.x), numberType.times(x, other.z))
            /*z =*/ 2 -> numberType.minus(numberType.times(x, other.y), numberType.times(y, other.x))
            else -> throw IllegalStateException("Should not happend.")
        }
    }

    /**
     * Gets a vector orthogonal to the current one with the same norm and counterclockwise to it.
     *
     * @return An orthogonal vector instance.
     */
    fun ortho(): V {
        val k: Int = largestAbsComponent()
        val temp = when (k) {
            1 -> newInstance { i -> if (i == 0) numberType.one() else numberType.zero() }
            2 -> newInstance { i -> if (i == 1) numberType.one() else numberType.zero() }
            else -> newInstance { i -> if (i == 2) numberType.one() else numberType.zero() }
        }
        val ortho = crossProd(temp)
        ortho.normalize()
        return ortho
    }

    /**
     * return the index of the largest component (abs)
     */
    fun largestAbsComponent(): Int {
        val temp = abs()
        return if (temp[0] > temp[1]) {
            if (temp[0] > temp[2]) 0 else 2
        } else {
            if (temp[1] > temp[2]) 1 else 2
        }
    }

    /**
     * Gets the angle between "this" and v in radians. If either vector is zero-length, or nearly zero-length,
     * the result will be zero, regardless of the other value.
     *
     * @param v Another vector.
     * @return The angle between this vector and v.
     */
    fun angle(v: V): T = numberType.atan2(crossProd(v).norm(), dotProd(v))

}


