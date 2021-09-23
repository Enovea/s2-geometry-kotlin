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

//
// S2Error is a simple class consisting of an error code and a human-readable
// error message.


// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
data class S2Error(var code: Int = OK, var text: String = ""): Exception("S2Error $code: $text") {

    fun isOk(): Boolean = code == OK

    fun isNotOk(): Boolean = !isOk()

    fun init(code: Int, text: String): S2Error {
        this.code = code
        this.text = text
        return this
    }

    fun init(error: S2Error): S2Error = init(error.code, error.text)

    fun clear() {
        init(OK, "")
    }

    companion object Code {
        const val OK = 0                  // No error.

        ////////////////////////////////////////////////////////////////////
        // Generic errors not specific to geometric objects:

        const val UNKNOWN = 1000              // Unknown error.
        const val UNIMPLEMENTED = 1001        // Operation is not implemented.
        const val OUT_OF_RANGE = 1002         // Argument is out of range.
        const val INVALID_ARGUMENT = 1003     // Invalid argument (other than a range error).
        const val FAILED_PRECONDITION = 1004  // Object is not in the required state.
        const val INTERNAL = 1005             // An internal invariant has failed.
        const val DATA_LOSS = 1006            // Data loss or corruption.
        const val RESOURCE_EXHAUSTED = 1007   // A resource has been exhausted.

        ////////////////////////////////////////////////////////////////////
        // Error codes in the following range can be defined by clients:

        const val USER_DEFINED_START = 1000000
        const val USER_DEFINED_END   = 9999999

        ////////////////////////////////////////////////////////////////////
        // Errors that apply to more than one type of geometry:

        const val NOT_UNIT_LENGTH = 1     // Vertex is not unit length.
        const val DUPLICATE_VERTICES = 2  // There are two identical vertices.
        const val ANTIPODAL_VERTICES = 3  // There are two antipodal vertices.

        ////////////////////////////////////////////////////////////////////
        // S2Loop errors:

        const val LOOP_NOT_ENOUGH_VERTICES = 100  // Loop with fewer than 3 vertices.
        const val LOOP_SELF_INTERSECTION = 101    // Loop has a self-intersection.

        ////////////////////////////////////////////////////////////////////
        // S2Polygon errors:

        const val POLYGON_LOOPS_SHARE_EDGE = 200  // Two polygon loops share an edge.
        const val POLYGON_LOOPS_CROSS = 201       // Two polygon loops cross.
        const val POLYGON_EMPTY_LOOP = 202        // Polygon has an empty loop.
        const val POLYGON_EXCESS_FULL_LOOP = 203  // Non-full polygon has a full loop.

        // InitOriented() was called and detected inconsistent loop orientations.
        const val POLYGON_INCONSISTENT_LOOP_ORIENTATIONS = 204

        // Loop depths don't correspond to any valid nesting hierarchy.
        const val POLYGON_INVALID_LOOP_DEPTH = 205

        // Actual polygon nesting does not correspond to the nesting hierarchy
        // encoded by the loop depths.
        const val POLYGON_INVALID_LOOP_NESTING = 206

        ////////////////////////////////////////////////////////////////////
        // S2Builder errors:

        // The S2Builder snap function moved a vertex by more than the specified
        // snap radius.
        const val BUILDER_SNAP_RADIUS_TOO_SMALL = 300

        // S2Builder expected all edges to have siblings (as specified by
        // S2Builder::GraphOptions::SiblingPairs::REQUIRE), but some were missing.
        const val BUILDER_MISSING_EXPECTED_SIBLING_EDGES = 301

        // S2Builder found an unexpected degenerate edge.  For example,
        // Graph::GetLeftTurnMap() does not support degenerate edges.
        const val BUILDER_UNEXPECTED_DEGENERATE_EDGE = 302

        // S2Builder found a vertex with (indegree != outdegree), which means
        // that the given edges cannot be assembled into loops.
        const val BUILDER_EDGES_DO_NOT_FORM_LOOPS = 303

        // The edges provided to S2Builder cannot be assembled into a polyline.
        const val BUILDER_EDGES_DO_NOT_FORM_POLYLINE = 304

        // There was an attempt to assemble a polygon from degenerate geometry
        // without having specified a predicate to decide whether the output is
        // the empty polygon (containing no points) or the full polygon
        // (containing all points).
        const val BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED = 305
    }

}
