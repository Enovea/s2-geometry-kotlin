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

import java.math.BigDecimal

/**
 * Abstract base class for vector of big decimals implementations.
 *
 * @param V The implementation class of the vector.
 * @since 1.0
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 */
abstract class AbstractRBigDecimalVector<V : AbstractRBigDecimalVector<V>>(coords: Array<BigDecimal>) : AbstractRVector<BigDecimal, V>(coords) {

    operator fun BigDecimal.plus(other: V): V = other + this

}

operator fun <V : AbstractRBigDecimalVector<V>> BigDecimal.times(other: V): V = other * this
