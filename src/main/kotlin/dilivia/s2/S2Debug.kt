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
package dilivia.s2

var debug: Boolean = System.getProperty("debug", "false").toBoolean()

// Class that allows the -ea validity checks to be enabled or disabled
// for specific objects (e.g., see S2Polygon).
enum class S2Debug {
    ALLOW,    // Validity checks are controlled by --s2debug
    DISABLE   // No validity checks even when --s2debug is true
}
