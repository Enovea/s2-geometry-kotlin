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
package dilivia.s2.builder.graph

import dilivia.s2.builder.EdgeType


enum class DegenerateEdges { DISCARD, DISCARD_EXCESS, KEEP }
enum class DuplicateEdges { MERGE, KEEP }
enum class SiblingPairs { DISCARD, DISCARD_EXCESS, KEEP, REQUIRE, CREATE }


// This class is only needed by S2Builder::Layer implementations.  A layer is
// responsible for assembling an S2Builder::Graph of snapped edges into the
// desired output format (e.g., an S2Polygon).  The GraphOptions class allows
// each Layer type to specify requirements on its input graph: for example, if
// DegenerateEdges::DISCARD is specified, then S2Builder will ensure that all
// degenerate edges are removed before passing the graph to the S2Layer::Build
// method.
data class GraphOptions(

        // Specifies whether the S2Builder input edges should be treated as
        // undirected.  If true, then all input edges are duplicated into pairs
        // consisting of an edge and a sibling (reverse) edge.  The layer
        // implementation is responsible for ensuring that exactly one edge from
        // each pair is used in the output, i.e. *only half* of the graph edges will
        // be used.  (Note that some values of the sibling_pairs() option
        // automatically take care of this issue by removing half of the edges and
        // changing edge_type() to DIRECTED.)
        //
        // DEFAULT: EdgeType::DIRECTED
        var edgeType: EdgeType = EdgeType.DIRECTED,

        // Controls how degenerate edges (i.e., an edge from a vertex to itself) are
        // handled.  Such edges may be present in the input, or they may be created
        // when both endpoints of an edge are snapped to the same output vertex.
        // The options available are:
        //
        // DISCARD: Discards all degenerate edges.  This is useful for layers that
        //          do not support degeneracies, such as S2PolygonLayer.
        //
        // DISCARD_EXCESS: Discards all degenerate edges that are connected to
        //                 non-degenerate edges.  (Any remaining duplicate edges can
        //                 be merged using DuplicateEdges::MERGE.)  This is useful
        //                 for simplifying polygons while ensuring that loops that
        //                 collapse to a single point do not disappear.
        //
        // KEEP: Keeps all degenerate edges.  Be aware that this may create many
        //       redundant edges when simplifying geometry (e.g., a polyline of the
        //       form AABBBBBCCCCCCDDDD).  DegenerateEdges::KEEP is mainly useful
        //       for algorithms that require an output edge for every input edge.
        //
        // DEFAULT: DegenerateEdges::KEEP
        var degenerateEdges: DegenerateEdges = DegenerateEdges.KEEP,

        // Controls how duplicate edges (i.e., edges that are present multiple
        // times) are handled.  Such edges may be present in the input, or they can
        // be created when vertices are snapped together.  When several edges are
        // merged, the result is a single edge labelled with all of the original
        // input edge ids.
        //
        // DEFAULT: DuplicateEdges::KEEP
        var duplicateEdges: DuplicateEdges = DuplicateEdges.KEEP,

        // Controls how sibling edge pairs (i.e., pairs consisting of an edge and
        // its reverse edge) are handled.  Layer types that define an interior
        // (e.g., polygons) normally discard such edge pairs since they do not
        // affect the result (i.e., they define a "loop" with no interior).  The
        // various options include:
        //
        // DISCARD: Discards all sibling edge pairs.
        //
        // DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept
        //                 if the result would otherwise be empty.  This is useful
        //                 for polygons with degeneracies (S2LaxPolygonShape), and
        //                 for simplifying polylines while ensuring that they are
        //                 not split into multiple disconnected pieces.
        //
        // KEEP: Keeps sibling pairs.  This can be used to create polylines that
        //       double back on themselves, or degenerate loops (with a layer type
        //       such as S2LaxPolygonShape).
        //
        // REQUIRE: Requires that all edges have a sibling (and returns an error
        //          otherwise).  This is useful with layer types that create a
        //          collection of adjacent polygons (a polygon mesh).
        //
        // CREATE: Ensures that all edges have a sibling edge by creating them if
        //         necessary.  This is useful with polygon meshes where the input
        //         polygons do not cover the entire sphere.  Such edges always
        //         have an empty set of labels.
        //
        // If edge_type() is EdgeType::UNDIRECTED, a sibling edge pair is considered
        // to consist of four edges (two duplicate edges and their siblings), since
        // only two of these four edges will be used in the final output.
        //
        // Furthermore, since the options REQUIRE and CREATE guarantee that all
        // edges will have siblings, S2Builder implements these options for
        // undirected edges by discarding half of the edges in each direction and
        // changing the edge_type() to EdgeType::DIRECTED.  For example, two
        // undirected input edges between vertices A and B would first be converted
        // into two directed edges in each direction, and then one edge of each pair
        // would be discarded leaving only one edge in each direction.
        //
        // Degenerate edges are considered not to have siblings.  If such edges are
        // present, they are passed through unchanged by SiblingPairs::DISCARD.  For
        // SiblingPairs::REQUIRE or SiblingPairs::CREATE with undirected edges, the
        // number of copies of each degenerate edge is reduced by a factor of two.
        //
        // Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and
        // REQUIRE/CREATE in the case of undirected edges) have the side effect that
        // when duplicate edges are present, all of the corresponding edge labels
        // are merged together and assigned to the remaining edges.  (This avoids
        // the problem of having to decide which edges are discarded.)  Note that
        // this merging takes place even when all copies of an edge are kept, and
        // that even labels attached to duplicate degenerate edges are merged.  For
        // example, consider the graph {AB1, AB2, BA3, CD4, CD5} (where XYn denotes
        // an edge from X to Y with label "n").  With SiblingPairs::DISCARD, we need
        // to discard one of the copies of AB.  But which one?  Rather than choosing
        // arbitrarily, instead we merge the labels of all duplicate edges (even
        // ones where no sibling pairs were discarded), yielding {AB12, CD45, CD45}
        // (assuming that duplicate edges are being kept).
        //
        // DEFAULT: SiblingPairs::KEEP
        var siblingPairs: SiblingPairs = SiblingPairs.KEEP,

        // This is a specialized option that is only needed by clients want to work
        // with the graphs for multiple layers at the same time (e.g., in order to
        // check whether the same edge is present in two different graphs).  [Note
        // that if you need to do this, usually it is easier just to build a single
        // graph with suitable edge labels.]
        //
        // When there are a large number of layers, then by default S2Builder builds
        // a minimal subgraph for each layer containing only the vertices needed by
        // the edges in that layer.  This ensures that layer types that iterate over
        // the vertices run in time proportional to the size of that layer rather
        // than the size of all layers combined.  (For example, if there are a
        // million layers with one edge each, then each layer would be passed a
        // graph with 2 vertices rather than 2 million vertices.)
        //
        // If this option is set to false, this optimization is disabled.  Instead
        // the graph passed to this layer will contain the full set of vertices.
        // (This is not recommended when the number of layers could be large.)
        //
        // DEFAULT: true
        var allowVertexFiltering: Boolean = true
)

