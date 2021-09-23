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
package dilivia.s2.index.shape

import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2Factory
import dilivia.s2.S2Point
import dilivia.s2.edge.S2EdgeClipping
import dilivia.s2.edge.S2EdgeClipping.intersectsRect
import dilivia.s2.edge.S2EdgeClipping.kIntersectsRectErrorUVDist
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.shape.S2EdgeVectorShape
import dilivia.s2.shape.S2Shape
import org.assertj.core.api.Assertions.assertThat
import org.assertj.core.api.Assertions.fail
import org.junit.jupiter.api.BeforeEach
import org.junit.jupiter.api.Disabled
import org.junit.jupiter.api.Test

class MutableS2ShapeIndexUnitTest {
    // This test harness owns a MutableS2ShapeIndex for convenience.
    private lateinit var index: MutableS2ShapeIndex

    @BeforeEach
    fun setUp() {
        index = MutableS2ShapeIndex()
    }

    // Verifies that that every cell of the index contains the correct edges, and
    // that no cells are missing from the index.  The running time of this
    // function is quadratic in the number of edges.
    private fun quadraticValidate() {
        // Iterate through a sequence of nonoverlapping cell ids that cover the
        // sphere and include as a subset all the cell ids used in the index.  For
        // each cell id, verify that the expected set of edges is present.

        // "min_cellid" is the first S2CellId that has not been validated yet.
        var min_cellid = S2CellId.begin(S2CellId.kMaxLevel)
        val it = index.cellIterator(InitialPosition.BEGIN)
        while (true) {
            // Generate a list of S2CellIds ("skipped cells") that cover the gap
            // between the last cell we validated and the next cell in the index.
            var skipped: S2CellUnion = S2CellUnion()
            if (!it.done()) {
                val cellid = it.id()
                assertThat(cellid >= min_cellid).isTrue()
                skipped = S2CellUnion.fromBeginEnd(min_cellid, cellid.rangeMin())
                min_cellid = cellid.rangeMax().next()
            } else {
                // Validate the empty cells beyond the last cell in the index.
                skipped = S2CellUnion.fromBeginEnd(min_cellid, S2CellId.end(S2CellId.kMaxLevel))
            }
            // Iterate through all the shapes, simultaneously validating the current
            // index cell and all the skipped cells.
            var short_edges = 0  // number of edges counted toward subdivision
            for (id in 0 until index.nextNewShapeId()) {
                val shape = (index.shape(id) ?: fail("")) as S2Shape?
                var clipped: S2ClippedShape? = null
                if (!it.done()) clipped = it.cell().findClipped(id)

                // First check that contains_center() is set correctly.
                for (skipped_id in skipped) {
                    validateInterior(shape, skipped_id, false)
                }
                if (!it.done()) {
                    val contains_center = clipped != null && clipped.containsCenter
                    validateInterior(shape, it.id(), contains_center)
                }
                // If this shape has been released, it should not be present at all.
                if (shape == null) {
                    assertThat(clipped).isNull()
                    continue
                }
                // Otherwise check that the appropriate edges are present.
                for (e in 0 until shape.numEdges) {
                    val edge = shape.edge(e)
                    for (j in 0 until skipped.numCells()) {
                        validateEdge(edge.v0, edge.v1, skipped.cellId(j), false)
                    }
                    if (!it.done()) {
                        val has_edge = clipped != null && clipped.containsEdge(e)
                        validateEdge(edge.v0, edge.v1, it.id(), has_edge)
                        val max_level = index.getEdgeMaxLevel(edge)
                        if (has_edge && it.id().level() < max_level) {
                            ++short_edges
                        }
                    }
                }
            }
            assertThat(short_edges <= index.options().maxEdgesPerCell).isTrue()
            if (it.done()) break
            it.next()
        }
    }

    // Given an edge and a cell id, determines whether or not the edge should be
    // present in that cell and verify that this matches "indexhas_edge".
    // Verify that "indexhas_edge" is true if and only if the edge AB intersects
    // the given cell id.
    private fun validateEdge(a: S2Point, b: S2Point, id: S2CellId, indexhas_edge: Boolean) {
        // Expand or shrink the padding slightly to account for errors in the
        // function we use to test for intersection (IntersectsRect).
        var padding = AbstractMutableS2ShapeIndex.kCellPadding
        padding += (if (indexhas_edge) 1 else -1) * kIntersectsRectErrorUVDist
        val bound = id.boundUV.expanded(padding)
        val clippedEdge = S2EdgeClipping.clipToPaddedFace(a, b, id.face(), padding)
        assertThat(indexhas_edge).isEqualTo(clippedEdge != null && intersectsRect(clippedEdge.first, clippedEdge.second, bound))
    }

    // Given a shape and a cell id, determines whether or not the shape contains
    // the cell center and verify that this matches "indexcontains_center".
    private fun validateInterior(shape: S2Shape?, id: S2CellId, indexcontains_center: Boolean) {
        if (shape == null) {
            assertThat(indexcontains_center).isFalse()
        } else {
            assertThat(indexcontains_center).isEqualTo(S2Shape.containsBruteForce(shape, id.toPoint()))
        }
    }


    private fun testIteratorMethods(index: MutableS2ShapeIndex) {
        val it = index.cellIterator(InitialPosition.BEGIN)
        assertThat(it.prev()).isFalse()
        it.finish()
        assertThat(it.done()).isTrue()
        val ids = mutableListOf<S2CellId>()
        var it2 = index.cellIterator()
        var min_cellid = S2CellId.begin(S2CellId.kMaxLevel)
        it.begin()
        while (!it.done()) {
            val cellid = it.id()
            val skipped = S2CellUnion.fromBeginEnd(min_cellid, cellid.rangeMin())
            for (skipped_id in skipped) {
                assertThat(it2.locate(skipped_id.toPoint())).isFalse()
                assertThat(it2.locate(skipped_id)).isEqualTo(CellRelation.DISJOINT)
                it2.begin()
                it2.seek(skipped_id)
                assertThat(it2.id()).isEqualTo(cellid)
            }
            if (ids.isNotEmpty()) {
                it2 = it.clone()
                assertThat(it2.prev()).isTrue()
                assertThat(it2.id()).isEqualTo(ids.last())
                it2.next()
                assertThat(it2.id()).isEqualTo(cellid)
                it2.seek(ids.last())
                assertThat(it2.id()).isEqualTo(ids.last())
            }
            it2.begin()
            assertThat(it.center()).isEqualTo(cellid.toPoint())
            assertThat(it2.locate(it.center())).isTrue()
            assertThat(it2.id()).isEqualTo(cellid)
            it2.begin()
            assertThat(it2.locate(cellid)).isEqualTo(CellRelation.INDEXED)
            assertThat(it2.id()).isEqualTo(cellid)
            if (!cellid.isFace) {
                it2.begin()
                assertThat(it2.locate(cellid.parent())).isEqualTo(CellRelation.SUBDIVIDED)
                assertThat(it2.id() <= cellid).isTrue()
                assertThat(it2.id() >= cellid.parent().rangeMin()).isTrue()
            }
            if (!cellid.isLeaf) {
                for (i in 0..3) {
                    it2.begin()
                    assertThat(it2.locate(cellid.child(i))).isEqualTo(CellRelation.INDEXED)
                    assertThat(it2.id()).isEqualTo(cellid)
                }
            }
            ids.add(cellid)
            min_cellid = cellid.rangeMax().next()
            it.next()
        }
    }

    @Test
    fun noEdges() {
        val it = index.cellIterator(InitialPosition.BEGIN)
        assertThat(it.done()).isTrue()
        testIteratorMethods(index)
    }

    @Test
    fun oneEdge() {
        assertThat(index.add(S2EdgeVectorShape(S2Point(1, 0, 0), S2Point(0, 1, 0)))).isEqualTo(0)
        quadraticValidate()
        testIteratorMethods(index)
    }

    @Test
    fun shrinkToFitOptimization() {
        // This used to trigger a bug in the ShrinkToFit optimization.  The loop
        // below contains almost all of face 0 except for a small region in the
        // 0/00000 subcell.  That subcell is the only one that contains any edges.
        // This caused the index to be built only in that subcell.  However, all the
        // other cells on that face should also have index entries, in order to
        // indicate that they are contained by the loop.
        val loop = S2Loop.makeRegularLoop(S2Point(1.0, 0.5, 0.5).normalize(), S1Angle.degrees(89), 100)
        index.add(S2Loop.Shape(loop = loop))
        quadraticValidate()
    }

    @Test
    fun loopsSpanningThreeFaces() {
        val polygon = S2Polygon()
        val kNumEdges = 100  // Validation is quadratic
        // Construct two loops consisting of kNumEdges vertices each, centered
        // around the cube vertex at the start of the Hilbert curve.
        S2Factory.concentricLoopsPolygon(S2Point(1, -1, -1).normalize(), 2, kNumEdges, polygon)
        val loops = polygon.loops()
        for (loop in loops) {
        index.add(S2Loop.Shape(loop = loop))
    }
        quadraticValidate()
        testIteratorMethods(index)
        //TestEncodeDecode()
    }

    @Test
    fun manyIdenticalEdges() {
        val kNumEdges = 100;  // Validation is quadratic
        val a = S2Point(0.99, 0.99, 1.0).normalize()
        val b = S2Point(-0.99, -0.99, 1.0).normalize()
        for (i in 0 until kNumEdges) {
            assertThat(index.add(S2EdgeVectorShape(a, b))).isEqualTo(i)
        }
        quadraticValidate()
        testIteratorMethods(index)
        // Since all edges span the diagonal of a face, no subdivision should
        // have occurred (with the default index options).
        val it = index.cellIterator(InitialPosition.BEGIN)
        while (!it.done()) {
            assertThat(it.id().level()).isEqualTo(0)
            it.next()
        }
    }

    @Test
    fun degenerateEdge() {
        // This test verifies that degenerate edges are supported.  The following
        // point is a cube face vertex, and so it should be indexed in 3 cells.
        val a = S2Point (1, 1, 1).normalize()
        val shape = S2EdgeVectorShape()
        shape.add(a, a)
        index.add(shape)
        quadraticValidate()
        // Check that exactly 3 index cells contain the degenerate edge.
        var count = 0
        val it = index.cellIterator(InitialPosition.BEGIN)
        while (!it.done()) {
            assertThat(it.id().isLeaf).isTrue()
            assertThat(it.cell().numClipped).isEqualTo(1)
            assertThat(it.cell().clipped(0).numEdges).isEqualTo(1)
            it.next(); ++count
        }
        assertThat(count).isEqualTo(3)
    }

    @Test
    fun manyTinyEdges() {
        // This test adds many edges to a single leaf cell, to check that
        // subdivision stops when no further subdivision is possible.
        val kNumEdges = 100;  // Validation is quadratic
        // Construct two points in the same leaf cell.
        val a = S2CellId.fromPoint(S2Point(1, 0, 0)).toPoint()
        val b =(a + S2Point(0.0, 1e-12, 0.0)).normalize()
        val shape = S2EdgeVectorShape()
        repeat(kNumEdges) { shape.add(a, b) }
        index.add(shape)
        quadraticValidate()
        // Check that there is exactly one index cell and that it is a leaf cell.
        val it = index.cellIterator(InitialPosition.BEGIN)
        assertThat(!it.done()).isTrue()
        assertThat(it.id().isLeaf).isTrue()
        it.next()
        assertThat(it.done()).isTrue()
    }

    @Test
    @Disabled("failed - TODO")
    fun testSimpleUpdates() {
//        // Add 5 loops one at a time, then release them one at a time,
//        // validating the index at each step.
//        S2Polygon polygon
//                S2Testing::ConcentricLoopsPolygon(S2Point(1, 0, 0), 5, 20, & polygon)
//        for (int i = 0; i < polygon.num_loops(); ++i) {
//            index.Add(make_unique < S2Loop::Shape > (polygon.loop(i)))
//            quadraticValidate()
//        }
//        for (int id = 0; id < polygon.num_loops(); ++id) {
//            index.Release(id)
//            quadraticValidate()
//            TestEncodeDecode()
//        }
    }

    @Test
    @Disabled("failed - TODO")
    fun randomUpdates() {
//        // Allow the seed to be varied from the command line.
//        S2Testing::rnd.Reset(FLAGS_s2_random_seed)
//
//        // A few polylines.
//        index.Add(make_unique < S2Polyline::OwningShape > (
//                MakePolyline("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")))
//        index.Add(make_unique < S2Polyline::OwningShape > (
//                MakePolyline("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")))
//        index.Add(make_unique < S2Polyline::OwningShape > (
//                MakePolyline("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")))
//
//        // A loop that used to trigger an indexing bug.
//        index.Add(make_unique < S2Loop::OwningShape > (S2Loop::MakeRegularLoop(
//                S2Point(1, 0.5, 0.5).normalize(), S1Angle::Degrees(89), 20)))
//
//        // Five concentric loops.
//        S2Polygon polygon5
//                S2Testing::ConcentricLoopsPolygon(S2Point(1, -1, -1).normalize(),
//                        5, 20, & polygon5)
//        for (int i = 0; i < polygon5.num_loops(); ++i) {
//            index.Add(make_unique < S2Loop::Shape > (polygon5.loop(i)))
//        }
//
//        // Two clockwise loops around S2Cell cube vertices.
//        index.Add(make_unique < S2Loop::OwningShape > (S2Loop::MakeRegularLoop(
//                S2Point(-1, 1, 1).normalize(), S1Angle::Radians(M_PI - 0.001), 10)))
//        index.Add(make_unique < S2Loop::OwningShape > (S2Loop::MakeRegularLoop(
//                S2Point(-1, -1, -1).normalize(), S1Angle::Radians(M_PI - 0.001), 10)))
//
//        // A shape with no edges and no interior.
//        index.Add(make_unique < S2Loop::OwningShape > (
//                make_unique<S2Loop>(S2Loop::kEmpty())))
//
//        // A shape with no edges that covers the entire sphere.
//        index.Add(make_unique < S2Loop::OwningShape > (
//                make_unique<S2Loop>(S2Loop::kFull())))
//
//        vector<unique_ptr<S2Shape>> released
//                vector<int> added (index.num_shape_ids())
//        std::iota(added.begin(), added.end(), 0)
//        quadraticValidate()
//        TestEncodeDecode()
//        for (int iter = 0; iter < 100; ++iter) {
//            S2_VLOG(1) < < "Iteration: " << iter
//            // Choose some shapes to add and release.
//            int num_updates = 1+S2Testing::rnd.Skewed(5)
//            for (int n = 0; n < num_updates; ++n) {
//            if (S2Testing::rnd.OneIn(2) && !added.empty()) {
//                int i = S2Testing ::rnd.Uniform(added.size())
//                S2_VLOG(1) < < "  Released shape " << added[i]
//                << " (" << index.shape(added[i]) << ")"
//                released.push_back(index.Release(added[i]))
//                added.erase(added.begin() + i)
//            } else if (!released.empty()) {
//                int i = S2Testing ::rnd.Uniform(released.size())
//                S2Shape * shape = released[i].get()
//                index.Add(std::move(released[i]));  // Changes shape->id().
//                released.erase(released.begin() + i)
//                added.push_back(shape->id())
//                S2_VLOG(1) < < "  Added shape " << shape->id()
//                << " (" << shape << ")"
//            }
//        }
//            quadraticValidate()
//            TestEncodeDecode()
//        }
    }
}

// TODO(fmeurisse) Test multi threading
/*
// A test that repeatedly updates "index" in one thread and attempts to
// concurrently read the index from several other threads.  When all threads
// have finished reading, the first thread makes another update.
//
// Note that we only test concurrent read access, since MutableS2ShapeIndex
// requires all updates to be single-threaded and not concurrent with any
// reads.
class LazyUpdatesTest : public ::testing::Test {
    public:
    LazyUpdatesTest() : num_updates_(0), num_readers_left_(0) {
    }

    // The function executed by each reader thread.
    void ReaderThread ()

    protected:
    class ReaderThreadPool {
        public :
        ReaderThreadPool(LazyUpdatesTest* test, int num_threads)
        : threads_(make_unique<std::thread[]>(num_threads)),
        num_threads_(num_threads)
        {
            for (int i = 0; i < num_threads_; ++i) {
            threads_[i] = std::thread(& LazyUpdatesTest ::ReaderThread, test)
        }
        }
        ~ReaderThreadPool()
        {
            for (int i = 0; i < num_threads_; ++i) threads_[i].join()
        }

        private :
        unique_ptr<std::thread[]> threads_
        int num_threads_
    }

    MutableS2ShapeIndex index
            // The following fields are guarded by lock_.
            absl::Mutex lock_
            int num_updates_
            int num_readers_left_

            // Signalled when a new update is ready to be processed.
            absl::CondVar update_ready_
            // Signalled when all readers have processed the latest update.
            absl::CondVar all_readers_done_
}

void LazyUpdatesTest::ReaderThread() {
    lock_.Lock()
    for (int last_update = 0;; last_update = num_updates_) {
        while (num_updates_ == last_update) {
            update_ready_.Wait(& lock_)
        }
        if (num_updates_ < 0) break

        // The index is built on demand the first time we attempt to use it.
        // We intentionally release the lock so that many threads have a chance
        // to access the MutableS2ShapeIndex in parallel.
        lock_.Unlock()
        for (MutableS2ShapeIndex:: Iterator it (&index, S2ShapeIndex::BEGIN)
        !it.done(); it.Next()) {
        continue;  // NOLINT
    }
        lock_.Lock()
        if (--num_readers_left_ == 0) {
            all_readers_done_.Signal()
        }
    }
    lock_.Unlock()
}

TEST_F(LazyUpdatesTest, ConstMethodsThreadSafe) {
    // Ensure that lazy updates are thread-safe.  In other words, make sure that
    // nothing bad happens when multiple threads call "const" methods that
    // cause pending updates to be applied.

    // The number of readers should be large enough so that it is likely that
    // several readers will be running at once (with a multiple-core CPU).
    const int kNumReaders = 8
    ReaderThreadPool pool (this, kNumReaders)
    lock_.Lock()
    const int kIters = 100
    for (int iter = 0; iter < kIters; ++iter) {
        // Loop invariant: lock_ is held and num_readers_left_ == 0.
        S2_DCHECK_EQ(0, num_readers_left_)
        // Since there are no readers, it is safe to modify the index.
        index.Clear()
        int num_vertices = 4 * S2Testing::rnd.Skewed(10);  // Up to 4K vertices
        unique_ptr<S2Loop> loop (S2Loop::MakeRegularLoop(
                S2Testing::RandomPoint(), S2Testing::KmToAngle(5), num_vertices))
        index.Add(make_unique < S2Loop::Shape > (loop.get()))
        num_readers_left_ = kNumReaders
        ++num_updates_
        update_ready_.SignalAll()
        while (num_readers_left_ > 0) {
            all_readers_done_.Wait(& lock_)
        }
    }
    // Signal the readers to exit.
    num_updates_ = -1
    update_ready_.SignalAll()
    lock_.Unlock()
    // ReaderThreadPool destructor waits for all threads to complete.
}

TEST(MutableS2ShapeIndex, MixedGeometry) {
    // This test used to trigger a bug where the presence of a shape with an
    // interior could cause shapes that don't have an interior to suddenly
    // acquire one.  This would cause extra S2ShapeIndex cells to be created
    // that are outside the bounds of the given geometry.
    vector<unique_ptr<S2Polyline>> polylines
            polylines.push_back(MakePolyline("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6"))
    polylines.push_back(MakePolyline("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6"))
    polylines.push_back(MakePolyline("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6"))
    MutableS2ShapeIndex index
            for (auto& polyline : polylines) {
        index.Add(make_unique < S2Polyline::OwningShape > (std::move(polyline)))
    }
    S2Loop loop (S2Cell(S2CellId::Begin(S2CellId::kMaxLevel)))
    index.Add(make_unique < S2Loop::Shape > (& loop))
    MutableS2ShapeIndex::Iterator it (&index)
    // No geometry intersects face 1, so there should be no index cells there.
    assertThat(it.Locate(S2CellId::FromFace(1))).isEqualTo(S2ShapeIndex::DISJOINT)
}

TEST(S2Shape, user_data) {
    struct MyData {
        int x, y
        MyData(int _x, int _y) : x(_x), y(_y) {}
    }
    class MyEdgeVectorShape : public S2EdgeVectorShape {
        public:
        explicit MyEdgeVectorShape (const MyData & data)
        : S2EdgeVectorShape(), data_(data) {
    }
        const void * user_data () const override { return &data_; }
        void * mutable_user_data() override { return &data_; }

        private:
        MyData data_
    }
    MyEdgeVectorShape shape (MyData(3, 5))
    MyData * data = static_cast < MyData * >(shape.mutable_user_data())
    S2_DCHECK_EQ(3, data->x)
    data->y = 10
    S2_DCHECK_EQ(10, static_cast < const MyData * > (shape.user_data())->y)
}
*/
