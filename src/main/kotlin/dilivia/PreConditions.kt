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

/**
 *
 * @author Fabien Meurisse
 * @since 1.0
 */
object PreConditions {

    var enabled: Boolean = javaClass.desiredAssertionStatus()

    inline fun checkState(assertion: () -> Boolean, lazyMessage: () -> Any) {
        if (enabled) {
            check(assertion(), lazyMessage)
        }
    }

    inline fun checkState(assertion: () -> Boolean) {
        if (enabled) {
            check(assertion())
        }
    }

    inline fun <T : Collection<V>, V> checkIsEmpty(value: T, lazyMessage: () -> Any = { "Collection is not empty." }) {
        checkState({ value.isEmpty() }, lazyMessage)
    }

    inline fun <T> checkEQ(value: T, expected: T, lazyMessage: () -> Any) {
        checkState({ value == expected }, lazyMessage)
    }

    fun <T> checkEQ(value: T, expected: T) {
        checkState { value == expected }
    }

    inline fun <T> checkNE(value: T, notExpected: T, lazyMessage: () -> Any) {
        checkState({ value != notExpected }, lazyMessage)
    }

    fun <T> checkNE(value: T, notExpected: T) {
        checkState { value != notExpected }
    }

    inline fun <T: Comparable<T>>  checkGT(value: T, limit: T, lazyMessage: () -> Any) {
        checkState({ value > limit }, lazyMessage)
    }

    fun <T: Comparable<T>> checkGT(value: T, limit: T) {
        checkState { value > limit }
    }

    inline fun <T: Comparable<T>>  checkGE(value: T, limit: T, lazyMessage: () -> Any) {
        checkState({ value >= limit }, lazyMessage)
    }

    fun <T: Comparable<T>> checkGE(value: T, limit: T) {
        checkState { value >= limit }
    }

    inline fun <T: Comparable<T>>  checkLT(value: T, limit: T, lazyMessage: () -> Any) {
        checkState({ value < limit }, lazyMessage)
    }

    fun <T: Comparable<T>> checkLT(value: T, limit: T) {
        checkState { value < limit }
    }

    inline fun <T: Comparable<T>>  checkLE(value: T, limit: T, lazyMessage: () -> Any) {
        checkState({ value <= limit }, lazyMessage)
    }

    fun <T: Comparable<T>> checkLE(value: T, limit: T) {
        checkState { value <= limit }
    }

    inline fun requireArgument(assertion: () -> Boolean, lazyMessage: () -> Any) {
        if (enabled) {
            require(assertion(), lazyMessage)
        }
    }

    inline fun requireArgument(assertion: () -> Boolean) {
        if (enabled) {
            require(assertion())
        }
    }

    inline fun <T: Comparable<T>> requireBetween(value: T, minInclusive: T, maxExclusive: T, lazyMessage: () -> Any) {
        requireArgument({ minInclusive <= value && value < maxExclusive }, lazyMessage)
    }

    fun <T: Comparable<T>>  requireBetween(value: T, minInclusive: T, maxExclusive: T) {
        requireArgument { minInclusive <= value && value < maxExclusive }
    }

    inline fun <T: Comparable<T>> requireEQ(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value == limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireEQ(value: T, limit: T) {
        requireArgument { value == limit }
    }

    inline fun <T: Comparable<T>> requireNE(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value != limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireNE(value: T, limit: T) {
        requireArgument { value != limit }
    }
    inline fun <T: Comparable<T>> requireGE(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value >= limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireGE(value: T, limit: T) {
        requireArgument { value >= limit }
    }

    inline fun <T: Comparable<T>> requireLE(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value <= limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireLE(value: T, limit: T) {
        requireArgument { value <= limit }
    }

    inline fun <T: Comparable<T>>  requireLT(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value < limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireLT(value: T, limit: T) {
        requireArgument { value < limit }
    }

    inline fun <T: Comparable<T>>  requireGT(value: T, limit: T, lazyMessage: () -> Any) {
        requireArgument({ value > limit }, lazyMessage)
    }

    fun <T: Comparable<T>> requireGT(value: T, limit: T) {
        requireArgument { value > limit }
    }

}
