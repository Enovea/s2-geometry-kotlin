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

import dilivia.PreConditions.requireArgument
import dilivia.math.NumberType

// Exception messages
const val check_message_dimensions_not_equals = "Vector dimensions are not the same %d != %d."
const val check_message_component_idx_out_of_range = "Component index %d is not in range 0..%d"

/**
 * Represents a vector of numbers. This interface defines common operations on vectors.
 *
 * @param T The type of numbers
 * @param V The implementation class of the vector.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
interface RVector<T, V> : Cloneable, Comparable<V> where T : Number, T : Comparable<T>, V : RVector<T, V> {

    val numberType: NumberType<T>

    /** The size of the vector. It is defined as the number of components of the vector. */
    val size: Int

    /**
     * Gets the i-th component of the vector.
     *
     * @param index The index of the desired component.
     * @return The value of the i-th component.
     * @throws IllegalArgumentException If index if outside the range 0..(size - 1).
     */
    operator fun get(index: Int): T

    /**
     * Sets the value of the i-th component of the vector.
     *
     * @param index The index of the component.
     * @param value The new value of the component.
     */
    operator fun set(index: Int, value: T)

    /**
     * Gets the sum of two vectors.
     *
     * @param other The vector of the same type to sum.
     * @return A new vector instance that represents the sum of this vector and the other one.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun plus(other: V): V {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return newInstance { i -> numberType.plus(this[i], other[i]) }
    }

    /**
     * Gets the sum of this vector with a scalar value. The scalar value is added to each component of the vector.
     *
     * @param other The scalar value to add.
     * @return A new vector instance that represents the sum of this vector with the scalar.
     */
    operator fun plus(other: T): V = newInstance { i -> numberType.plus(this[i], other) }

    /**
     * Assigns this vector with the result of the sum of this vector with another one. This method does not create a
     * new vector instance bu modify this one.
     *
     * @param other The vector to add.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun plusAssign(other: V) {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        (0 until size).forEach { i -> this[i] = numberType.plus(this[i], other[i]) }
    }

    /**
     * Assigns this vector with the result of the sum of this vector with a scalar value. This method does not create a
     * new vector instance bu modify this one.
     *
     * @param other The scalar to add.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun plusAssign(other: T) = (0 until size).forEach { i -> this[i] = numberType.plus(this[i], other) }

    /**
     * Defines the unary minus operation on a vector. A new instance of this vector is returned by multiplying each
     * component by -1.
     *
     * @return The unary minus instance of this vector.
     */
    operator fun unaryMinus(): V = newInstance { i -> numberType.unaryMinus(this[i]) }

    /**
     * Gets the subtraction of two vectors.
     *
     * @param other The vector to subtract.
     * @return A new vector instance that represents the subtraction of this vector and the other one.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun minus(other: V): V {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return newInstance { i -> numberType.minus(this[i], other[i]) }
    }

    /**
     * Assigns this vector with the result of the subtraction of this vector with another one. This method does not
     * create a new vector instance bu modify this one.
     *
     * @param other The vector to subtract.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun minusAssign(other: V) {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        (0 until size).forEach { i -> this[i] = numberType.minus(this[i], other[i]) }
    }

    /**
     * Gets the subtraction of this vector by a scalar value. The scalar value is subtracted to each component of
     * the vector.
     *
     * @param other The scalar value to subtract.
     * @return A new vector instance that represents th
     */
    operator fun minus(other: T): V = newInstance { i -> numberType.minus(this[i], other) }

    /**
     * Assigns this vector with the result of the subtraction of this vector with a scalar value. This method does not
     * create a new vector instance but modify this one.
     *
     * @param other The scalar to subtract.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun minusAssign(other: T) = (0 until size).forEach { i -> this[i] = numberType.minus(this[i], other) }

    /**
     * Gets the multiplication of this vector by another one. This method returns a new instance of this vector
     * multiplied by the other vector. The two vector must have the same dimensions.
     * Each component of the result is the multiplication of the components of the two vector components at the same
     * indice.
     *
     * @param other The vector to multiply.
     * @return A new vector instance that represents the multiplication of this vector by the other one.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun times(other: V): V {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return newInstance { i -> numberType.times(this[i], other[i]) }
    }

    /**
     * Gets the multiplication of this vector by a scalar value.
     *
     * @param other The scalar value.
     * @return A new vector instance that is the result of the multiplication of this vector by the specified scalar
     * value.
     */
    operator fun times(other: T): V = newInstance { i -> numberType.times(this[i], other) }

    /**
     * Assigns this vector with the result of the subtraction of this vector with another one. This method does not
     * create a new vector instance but modify this one.
     *
     * @param other The vector to subtract.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun timesAssign(other: V) {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        (0 until size).forEach { i -> this[i] = numberType.times(this[i], other[i]) }
    }

    /**
     * Assigns this vector with the result of the multiplication of this by a scalar value. This method does not
     * create a new vector instance bu modify this one.
     *
     * @param other A scalar value.
     */
    operator fun timesAssign(other: T) = (0 until size).forEach { i -> this[i] = numberType.times(this[i], other) }

    /**
     * Gets the result of the division of this vector by another one. The division is defined as the vector with
     * each component equals to the division of the component i of this vector divided by the component i of the
     * second vector.
     *
     * @param other A vector.
     * @return The result of the division.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun div(other: V): V {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return newInstance { i -> numberType.div(this[i], other[i]) }
    }

    /**
     * Gets the division of this vector by a scalar value. Each component of this vector is divided by the scalar value.
     *
     * @param other A scalar value.
     * @return The result of the division.
     */
    operator fun div(other: T): V = newInstance { i -> numberType.div(this[i], other) }

    /**
     * Assigns to this vector the result of the division of this vector by another one. This method does not
     * create a new vector instance but modify this one.
     *
     * @param other A vector.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    operator fun divAssign(other: V) {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        (0 until size).forEach { i -> this[i] = numberType.div(this[i], other[i]) }
    }

    /**
     * Assigns this vector with the result of the division of this by a scalar value. This method does not
     * create a new vector instance but modify this one.
     *
     * @param other A scalar value.
     */
    operator fun divAssign(other: T) = (0 until size).forEach { i -> this[i] = numberType.div(this[i], other) }

    /**
     * Gets the dot product of this vector with the other one.
     *
     * @param other A vector.
     * @return The dot product.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    fun dotProd(other: RVector<T, V>): T {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return numberType.sum((0 until size).map { i: Int -> numberType.times(this[i], other[i]) })
    }

    /**
     * Computes the squared Euclidean norm (the dot product with itself).
     *
     * @return the squared Euclidean norm.
     */
    fun norm2(): T = this.dotProd(this)

    /**
     * Computes the Euclidean norm.
     *
     * @return the Euclidean norm.
     */
    fun norm(): T = numberType.sqrt(norm2())

    /**
     * Compose a vector from the sqrt of each component.
     *
     * @return A new vector instance with each component equals to the square root of the component of this vector.
     */
    fun sqrt(): V = newInstance { i -> numberType.sqrt(this[i]) }

    /**
     * Normalizes this vector.
     *
     * @return this
     * @throws IllegalStateException If the norm of this vector is 0.
     */
    fun normalize(): V

    /**
     * Gets the normalized version of this vector vector if the norm is nonzero.
     *
     * @return A new vector instance that represents the normalized version of this vector.
     * @throws IllegalStateException If the norm of this vector is 0.
     */
    fun normalized(): V {
        val normalized = clone()
        normalized.normalize()
        return normalized
    }

    /**
     * Compose a vector from the absolute value of each component.
     *
     * @return A new vector instance with x_i = |x_i|
     */
    fun abs(): V = newInstance { i -> numberType.abs(this[i]) }

    /**
     * Creates and returns a copy of this object.
     *
     * @return A cloned instance of this vector.
     */
    public override fun clone(): V = newInstance { i -> this[i] }

    /**
     * Compares vectors according to their lexicographical order.
     *
     * @return A negative value if this vector is before the other one, 0 if they are equals and a positive value if
     * this vector is after the other one.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same.
     */
    override fun compareTo(other: V): Int {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        for (i in 0 until size) {
            val compareTo = numberType.compare(this[i], other[i])
            if (i == size - 1 || compareTo != 0) return compareTo
        }
        throw IllegalStateException("Return should have been called in the for loop.")
    }

    /**
     * Checks if this vector is approximately equals to another with a margin tolerance.
     *
     * @param other A vector
     * @return true if the vectors are approximately equals.
     * @throws IllegalArgumentException If dimensions (size) of the two vectors are not the same
     */
    fun approxEquals(other: V, margin: Double = 1e-15): Boolean {
        requireArgument({ this.size == other.size }, { check_message_dimensions_not_equals.format(this.size, other.size) })
        return (0 until size).all { i -> numberType.abs(numberType.minus(this[i], other[i])) < numberType.cast(margin) }
    }

    /**
     * Creates a new vector instance from an array of numbers.
     *
     * @param coords The components of the new vector instance.
     * @return The vector instance.
     */
    fun newInstance(coords: Array<T>): V

    /**
     * Creates a new vector instance from a lambda function that associated a value to each component index.
     *
     * @param init The components initializer.
     * @return The vector instance.
     */
    fun newInstance(init: (Int) -> T): V

}


