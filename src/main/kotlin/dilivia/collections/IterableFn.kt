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
package dilivia.collections

fun <T> containsAll(iterable: Iterable<T>, elements: Collection<T>): Boolean {
    var allContains = false
    val containedElements = mutableSetOf<T>()
    val iterator = iterable.iterator()
    while (!allContains && iterator.hasNext()) {
        val next = iterator.next()
        if (next in elements) {
            containedElements.add(next)
            allContains = containedElements.containsAll(elements)
        }
    }
    return allContains
}

fun <T> retainAll(iterable: MutableIterable<T>, elements: Collection<T>): Boolean {
    var modified = false
    val iterator = iterable.iterator()
    while (iterator.hasNext()) {
        val next = iterator.next()
        if (next !in elements) {
            iterator.remove()
            modified = true
        }
    }
    return modified
}

fun <T> removeAll(iterable: MutableIterable<T>, elements: Collection<T>): Boolean {
    var modified = false
    val iterator = iterable.iterator()
    while (iterator.hasNext()) {
        val next = iterator.next()
        if (next in elements) {
            iterator.remove()
            modified = true
        }
    }
    return modified
}
