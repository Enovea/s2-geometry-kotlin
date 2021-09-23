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
package dilivia.math.matrix

import dilivia.PreConditions.checkNE
import dilivia.PreConditions.requireBetween
import dilivia.math.DoubleType
import dilivia.math.NumberType
import dilivia.math.vectors.R2VectorDouble
import dilivia.math.vectors.R3VectorDouble
import org.apache.commons.math3.util.FastMath

//
// A simple class to handle 3x3 matrices
// The aim of this class is to be able to manipulate 3x3 matrices
// and 3D vectors as naturally as possible and make calculations
// readable.
// For that reason, the operators +, -, * are overloaded.
// (Reading a = a + b*2 - c is much easier to read than
// a = Sub(Add(a, Mul(b,2)),c)   )
//
// Please be careful about overflows when using those matrices wth integer types
// The calculations are carried with VType. eg : if you are using uint8 as the
// base type, all values will be modulo 256.

/**
 * Interface that represents a 3x3 matrix of numbers.
 *
 * @author Fabien Meurisse
 */
interface Matrix3x3<T, V> where T : Number, T : Comparable<T>, V : Matrix3x3<T, V> {

    /** Matrix elements. */
    val m: Array<Array<T>>

    /** Scalar type operations. */
    val numberType: NumberType<T>

    /**
     * Change the value of all the coefficients of the matrix.
     *
     * @param m00 value of element (0,0)
     * @param m01 value of element (0,1)
     * @param m02 value of element (0,2)
     * @param m10 value of element (1,0)
     * @param m11 value of element (1,1)
     * @param m12 value of element (1,2)
     * @param m20 value of element (2,0)
     * @param m21 value of element (2,1)
     * @param m22 value of element (2,2)
     * @return This matrix after the modification (not a copy).
     */
    fun set(
        m00: T, m01: T, m02: T,
        m10: T, m11: T, m12: T,
        m20: T, m21: T, m22: T
    ): V {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22
        @Suppress("UNCHECKED_CAST")
        return this as V
    }

    // Return matrix element (i,j) with 0<=i<=2 0<=j<=2
    fun get(i: Int, j: Int): T {
        requireBetween(i, 0, 3)
        requireBetween(j, 0, 3)
        return m[i][j]
    }

    // Return matrix element (i/3,i%3) with 0<=i<=8 (access concatenated rows).
    operator fun get(i: Int): T {
        requireBetween(i, 0, 9)
        return get(i / 3, i % 3)
    }

    // Matrix addition
    operator fun plus(mb: V): V {
        @Suppress("UNCHECKED_CAST") val result = newInstance(this as V)
        result += mb
        return result
    }

    operator fun plusAssign(mb: V) {
        for (i in 0..2) {
            for (j in 0..2) {
                m[i][j] = numberType.plus(m[i][j], mb.get(i, j))
            }
        }
    }

    // Matrix subtraction
    operator fun minus(mb: V): V {
        @Suppress("UNCHECKED_CAST") val result = newInstance(this as V)
        result -= mb
        return result
    }

    operator fun minusAssign(mb: V) {
        for (i in 0..2) {
            for (j in 0..2) {
                m[i][j] = numberType.minus(m[i][j], mb.get(i, j))
            }
        }
    }

    // Change the sign of all the coefficients in the matrix
    operator fun unaryMinus(): V {
        for (i in 0..2) {
            for (j in 0..2) {
                m[i][j] = numberType.unaryMinus(m[i][j])
            }
        }
        @Suppress("UNCHECKED_CAST")
        return this as V
    }

    // Matrix multiplication by a scalar
    operator fun times(k: T): V {
        @Suppress("UNCHECKED_CAST") val result = newInstance(this as V)
        result *= k
        return result
    }

    operator fun timesAssign(k: T) {
        for (i in 0..2) {
            for (j in 0..2) {
                m[i][j] = numberType.times(k, m[i][j])
            }
        }
    }

    // Matrix multiplication
    operator fun times(mb: V): V =
        newInstance { i, j ->
            var cij = numberType.zero()
            for (k in 0..2) {
                cij = numberType.plus(cij, numberType.times(m[i][k], mb.get(k, j)))
            }
            cij
        }

    // Return the trace of the matrix
    fun trace(): T = numberType.sum((0..2).map { m[it][it] })

    // Return the transposed matrix
    fun transpose(): V = newInstance { i, j -> m[j][i] }

    // Return the Frobenius norm of the matrix: sqrt(sum(aij^2))
    fun frobeniusNorm(): T {
        var sum = numberType.zero()
        for (i in 0..2) {
            for (j in 0..2) {
                sum = numberType.plus(sum, numberType.times(m[i][j], m[i][j]))
            }
        }
        return numberType.sqrt(sum)
    }

    fun newInstance(matrix3x3: V): V

    /**
     */
    fun newInstance(coords: Array<Array<T>>): V

    /**
     */
    fun newInstance(initializer: (Int, Int) -> T): V
}

class Matrix3x3Double(override val m: Array<Array<Double>>) : Matrix3x3<Double, Matrix3x3Double> {

    override val numberType: DoubleType = DoubleType

    // Initialize the matrix to 0
    constructor() : this(m = Array(3) { arrayOf<Double>(0.0, 0.0, 0.0) })

    // Constructor explicitly setting the values of all the coefficient of
    // the matrix
    constructor(
        m00: Double, m01: Double, m02: Double,
        m10: Double, m11: Double, m12: Double,
        m20: Double, m21: Double, m22: Double
    ) : this(
        arrayOf(
            arrayOf<Double>(m00, m01, m02),
            arrayOf<Double>(m10, m11, m12),
            arrayOf<Double>(m20, m21, m22)
        )
    )

    constructor(other: Matrix3x3Double) : this((0..2).map { other.m[it].copyOf() }.toTypedArray())

    // Multiplication of a matrix by a vector
    operator fun times(v: R3VectorDouble): R3VectorDouble {
        return v.newInstance(
            arrayOf(
                m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
                m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
                m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
            )
        )
    }

    // Return the determinant of the matrix
    fun det(): Double =
        m[0][0] * m[1][1] * m[2][2] +
                m[0][1] * m[1][2] * m[2][0] +
                m[0][2] * m[1][0] * m[2][1] -
                m[2][0] * m[1][1] * m[0][2] -
                m[2][1] * m[1][2] * m[0][0] -
                m[2][2] * m[1][0] * m[0][1]

    // Return the transposed of the matrix of the cofactors
    // (Useful for inversion for example)
    fun comatrixTransposed(): Matrix3x3Double = Matrix3x3Double(
        m[1][1] * m[2][2] - m[2][1] * m[1][2],
        m[2][1] * m[0][2] - m[0][1] * m[2][2],
        m[0][1] * m[1][2] - m[1][1] * m[0][2],

        m[1][2] * m[2][0] - m[2][2] * m[1][0],
        m[2][2] * m[0][0] - m[0][2] * m[2][0],
        m[0][2] * m[1][0] - m[1][2] * m[0][0],

        m[1][0] * m[2][1] - m[2][0] * m[1][1],
        m[2][0] * m[0][1] - m[0][0] * m[2][1],
        m[0][0] * m[1][1] - m[1][0] * m[0][1]
    )

    // Matrix inversion
    fun inverse(): Matrix3x3Double {
        val det = det()
        checkNE(det, 0.0) { "Can't inverse. Determinant = 0." }
        return comatrixTransposed() * (1.0 / det)
    }

    // Return the vector 3D at row i
    fun row(i: Int): R3VectorDouble {
        requireBetween(i, 0, 3)
        return R3VectorDouble(m[i][0], m[i][1], m[i][2])
    }

    // Return the vector 3D at col i
    fun col(i: Int): R3VectorDouble {
        requireBetween(i, 0, 3)
        return R3VectorDouble(m[0][i], m[1][i], m[2][i])
    }

    // Set the vector in row i to be v1
    fun setRow(i: Int, v: R3VectorDouble) {
        requireBetween(i, 0, 3)
        m[i][0] = v[0];
        m[i][1] = v[1];
        m[i][2] = v[2];
    }

    // Set the vector in column i to be v1
    fun setCol(i: Int, v: R3VectorDouble) {
        requireBetween(i, 0, 3)
        m[0][i] = v[0];
        m[1][i] = v[1];
        m[2][i] = v[2];
    }

    // Return a matrix M close to the original but verifying MtM = I
    // (useful to compensate for errors in a rotation matrix)
    fun orthogonalize(): Matrix3x3Double {
        val r1 = row(0).normalize()
        val r2 = (row(2).crossProd(r1)).normalize()
        val r3 = (r1.crossProd(r2)).normalize()
        return fromRows(r1, r2, r3)
    }

    // Returns v.Transpose() * (*this) * u
    fun mulBothSides(v: R3VectorDouble, u: R3VectorDouble): Double {
        return (this * u).dotProd(v)
    }

    // Use the 3x3 matrix as a projective transform for 2d points
    fun project(v: R2VectorDouble): R2VectorDouble {
        val temp = this * R3VectorDouble(v[0], v[1], 1.0)
        return R2VectorDouble(temp[0] / temp[2], temp[1] / temp[2])
    }


    // Return true is one of the elements of the matrix is NaN
    fun isNaN(): Boolean = (0..8).any { get(it).isNaN() }

    override fun newInstance(coords: Array<Array<Double>>): Matrix3x3Double = Matrix3x3Double(coords)

    override fun newInstance(initializer: (Int, Int) -> Double): Matrix3x3Double = Matrix3x3Double(
        Array(3) { i -> Array(3) { j -> initializer(i, j) } }
    )

    override fun newInstance(matrix3x3: Matrix3x3Double): Matrix3x3Double = Matrix3x3Double(
        Array(3) { i -> Array(3) { j -> matrix3x3.get(i, j) } }
    )

    override fun equals(other: Any?): Boolean = when {
        this === other -> true
        other !is Matrix3x3Double -> false
        else -> (0..8).all { i -> this[i] == other[i] }
    }

    override fun hashCode(): Int {
        return m.contentDeepHashCode()
    }

    override fun toString(): String {
        var out = ""
        for (i in 0..2) {
            out += if (i == 0) {
                "["
            } else {
                " "
            }
            for (j in 0..2) {
                out += "" + m[i][j] + " "
            }
            if (i == 2) {
                out += "]"
            } else {
                out += "\n"
            }
        }
        return out
    }

    companion object {


        // Create a matrix from 3 row vectors
        fun fromRows(v1: R3VectorDouble, v2: R3VectorDouble, v3: R3VectorDouble): Matrix3x3Double = Matrix3x3Double(
            v1[0], v1[1], v1[2],
            v2[0], v2[1], v2[2],
            v3[0], v3[1], v3[2]
        )

        // Create a matrix from 3 column vectors
        fun fromCols(v1: R3VectorDouble, v2: R3VectorDouble, v3: R3VectorDouble): Matrix3x3Double = Matrix3x3Double(
            v1[0], v2[0], v3[0],
            v1[1], v2[1], v3[1],
            v1[2], v2[2], v3[2]
        )


        @JvmStatic
        fun fromCols(cols: Triple<R3VectorDouble, R3VectorDouble, R3VectorDouble>): Matrix3x3Double =
            fromCols(cols.first, cols.second, cols.third)

        // Return the identity matrix
        fun identity(): Matrix3x3Double = Matrix3x3Double(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        )

        // Return a matrix full of zeros
        fun zero(): Matrix3x3Double = Matrix3x3Double()

        // Return a diagonal matrix with the coefficients in v
        fun diagonal(v: R3VectorDouble): Matrix3x3Double {
            return Matrix3x3Double(
                v[0], 0.0, 0.0,
                0.0, v[1], 0.0,
                0.0, 0.0, v[2]
            )
        }

        // Return the matrix vvT
        fun sym3(v: R3VectorDouble): Matrix3x3Double = Matrix3x3Double(
            v[0] * v[0], v[0] * v[1], v[0] * v[2],
            v[1] * v[0], v[1] * v[1], v[1] * v[2],
            v[2] * v[0], v[2] * v[1], v[2] * v[2]
        )

        // Return a matrix M such that:
        // for each u,  M * u = v.CrossProd(u)
        fun antiSym3(v: R3VectorDouble): Matrix3x3Double = Matrix3x3Double(
            0.0, -v[2], v[1],
            v[2], 0.0, -v[0],
            -v[1], v[0], 0.0
        )

        // Returns matrix that rotates |rot| radians around axis rot.
        fun rodrigues(rot: R3VectorDouble): Matrix3x3Double {
            val theta = rot.norm()
            val w = rot.normalize()
            val wv = antiSym3(w)
            val ident = identity()
            val a = sym3(w)
            return (1 - FastMath.cos(theta)) * a + FastMath.sin(theta) * wv + FastMath.cos(theta) * ident;

        }

    }
}

operator fun Double.times(m: Matrix3x3Double): Matrix3x3Double = m * this

