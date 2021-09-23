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
package dilivia.s2.index.point

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2Point
import dilivia.s2.edge.S2EdgeDistances
import dilivia.s2.index.Delta
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceCellTarget
import dilivia.s2.index.S2MinDistanceEdgeTarget
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistancePointTarget
import dilivia.s2.index.S2MinDistanceShapeIndexTarget
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.index.point.S2ClosestPointQueryBase.Result
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2Region

//
// See S2ClosestPointQuery (defined below) for an overview.

// Given a set of points stored in an S2PointIndex, S2ClosestPointQuery
// provides methods that find the closest point(s) to a given query point
// or query edge.  Example usage:
//
// void Test(const vector<S2Point>& index_points,
//           const vector<S2Point>& target_points) {
//   // The template argument allows auxiliary data to be attached to each
//   // point (in this case, the array index).
//   S2PointIndex<int> index;
//   for (int i = 0; i < index_points.size(); ++i) {
//     index.Add(index_points[i], i);
//   }
//   S2ClosestPointQuery<int> query(&index);
//   query.mutable_options()->set_max_results(5);
//   for (const S2Point& target_point : target_points) {
//     S2ClosestPointQueryPointTarget target(target_point);
//     for (const auto& result : query.FindClosestPoints(&target)) {
//       // The Result class contains the following methods:
//       //   distance() is the distance to the target.
//       //   point() is the indexed point.
//       //   data() is the auxiliary data.
//       DoSomething(target_point, result);
//     }
//   }
// }
//
// You can find either the k closest points, or all points within a given
// radius, or both (i.e., the k closest points up to a given maximum radius).
// E.g. to find all the points within 5 kilometers, call
//
//   query.mutable_options()->set_max_distance(
//       S2Earth::ToAngle(util::units::Kilometers(5)));
//
// By default *all* points are returned, so you should always specify either
// max_results() or max_distance() or both.  There is also a FindClosestPoint()
// convenience method that returns only the closest point.
//
// You can restrict the results to an arbitrary S2Region, for example:
//
//   S2LatLngRect rect(...);
//   query.mutable_options()->set_region(&rect);  // Does *not* take ownership.
//
// To find the closest points to a query edge rather than a point, use:
//
//   S2ClosestPointQueryEdgeTarget target(v0, v1);
//   query.FindClosestPoints(&target);
//
// Similarly you can find the closest points to an S2Cell by using an
// S2ClosestPointQuery::CellTarget, and you can find the closest points to an
// arbitrary collection of points, polylines, and polygons by using an
// S2ClosestPointQuery::ShapeIndexTarget.
//
// The implementation is designed to be fast for both small and large
// point sets.
class S2ClosestPointQuery<T : Comparable<T>> {
    // See S2ClosestPointQueryBase for full documentation.
    private lateinit var options: Options
    private val base: S2ClosestPointQueryBase<S2MinDistance, T> = S2ClosestPointQueryBase(S2MinDistanceFactory)

    // Convenience constructor that calls Init().  Options may be specified here
    // or changed at any time using the mutable_options() accessor method.
    constructor(index: S2PointIndex<T>, options: Options = Options()) {
        init(index, options)
    }

    // Default constructor; requires Init() to be called.
    constructor() {
        TODO()
    }

    // Initializes the query.  Options may be specified here or changed at any
    // time using the mutable_options() accessor method.
    //
    // REQUIRES: "index" must persist for the lifetime of this object.
    // REQUIRES: ReInit() must be called if "index" is modified.
    fun init(index: S2PointIndex<T>, options: Options = Options()) {
        this.options = options
        base.init(index)
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    fun reInit() {
        base.reInit()
    }

    // Returns a reference to the underlying S2PointIndex.
    fun index(): S2PointIndex<T> = base.index()

    fun options(): Options = options

    // Returns the closest points to the given target that satisfy the current
    // options.  This method may be called multiple times.
    fun findClosestPoints(target: S2MinDistanceTarget): List<Result<S2MinDistance, T>> = base.findClosestPoints(target, options)

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    fun findClosestPoints(target: S2MinDistanceTarget, results: MutableList<Result<S2MinDistance, T>>) {
        base.findClosestPoints(target, options, results)
    }

    //////////////////////// Convenience Methods ////////////////////////

    // Returns the closest point to the target.  If no point satisfies the search
    // criteria, then a Result object with distance() == Infinity() and
    // is_empty() == true is returned.
    fun findClosestPoint(target: S2MinDistanceTarget): Result<S2MinDistance, T> {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResult(1)
        return base.findClosestPoint(target, tmpOptions)
    }

    // Returns the minimum distance to the target.  If the index or target is
    // empty, returns S1ChordAngle::Infinity().
    //
    // Use IsDistanceLess() if you only want to compare the distance against a
    // threshold value, since it is often much faster.
    fun getDistance(target: S2MinDistanceTarget): S1ChordAngle = findClosestPoint(target).distance.value

    // Returns true if the distance to "target" is less than "limit".
    //
    // This method is usually much faster than GetDistance(), since it is much
    // less work to determine whether the minimum distance is above or below a
    // threshold than it is to calculate the actual minimum distance.
    fun isDistanceLess(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResult(1)
        tmpOptions.maxDistance = S2MinDistance(limit)
        tmpOptions.maxError = S1ChordAngle.straight()
        return !base.findClosestPoint(target, tmpOptions).isEmpty()
    }

    // Like IsDistanceLess(), but also returns true if the distance to "target"
    // is exactly equal to "limit".
    fun isDistanceLessOrEqual(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmp_options = options.clone()
        tmp_options.setMaxResult(1)
        tmp_options.setInclusiveMaxDistance(limit)
        tmp_options.maxError  = S1ChordAngle.straight()
        return !base.findClosestPoint(target, tmp_options).isEmpty()
    }

    // Like IsDistanceLessOrEqual(), except that "limit" is increased by the
    // maximum error in the distance calculation.  This ensures that this
    // function returns true whenever the true, exact distance is less than
    // or equal to "limit".
    //
    // For example, suppose that we want to test whether two geometries might
    // intersect each other after they are snapped together using S2Builder
    // (using the IdentitySnapFunction with a given "snap_radius").  Since
    // S2Builder uses exact distance predicates (s2predicates.h), we need to
    // measure the distance between the two geometries conservatively.  If the
    // distance is definitely greater than "snap_radius", then the geometries
    // are guaranteed to not intersect after snapping.
    fun isConservativeDistanceLessOrEqual(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResult(1)
        tmpOptions.setConservativeMaxDistance(limit)
        tmpOptions.maxError = S1ChordAngle.straight()
        return !base.findClosestPoint(target, tmpOptions).isEmpty()
    }

    // Options that control the set of points returned.  Note that by default
    // *all* points are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    //
    // This class is also available as S2ClosestPointQuery<Data>::Options.
    // (It is defined here to avoid depending on the "Data" template argument.)
    class Options() : S2ClosestPointQueryBase.Options<S2MinDistance>(S2MinDistanceFactory), Cloneable {

        constructor(maxResult: Int = kMaxMaxResults,
                    maxDistance: S2MinDistance = S2MinDistanceFactory.infinity(),
                    maxError: Delta = Delta.zero(),
                    region: S2Region? = null,
                    useBruteForce: Boolean = false): this() {
            this.setMaxResult(maxResult)
            this.maxDistance = maxDistance
            this.maxError = maxError
            this.region = region
            this.useBruteForce = useBruteForce
        }

        // See S2ClosestPointQueryBaseOptions for the full set of options.

        // Specifies that only points whose distance to the target is less than
        // "max_distance" should be returned.
        //
        // Note that points whose distance is exactly equal to "max_distance" are
        // not returned.  Normally this doesn't matter, because distances are not
        // computed exactly in the first place, but if such points are needed then
        // see set_inclusive_max_distance() below.
        //
        // DEFAULT: Distance::Infinity()
        fun setMaxDistance(maxDistance: S1ChordAngle) {
            this.maxDistance = S2MinDistance(maxDistance)
        }

        // Like set_max_distance(), except that points whose distance is exactly
        // equal to "max_distance" are also returned.  Equivalent to calling
        // set_max_distance(max_distance.Successor()).
        fun setInclusiveMaxDistance(maxDistance: S1ChordAngle) {
            setMaxDistance(maxDistance.successor())
        }

        // Like set_inclusive_max_distance(), except that "max_distance" is also
        // increased by the maximum error in the distance calculation.  This ensures
        // that all points whose true distance is less than or equal to
        // "max_distance" will be returned (along with some points whose true
        // distance is slightly greater).
        //
        // Algorithms that need to do exact distance comparisons can use this
        // option to find a set of candidate points that can then be filtered
        // further (e.g., using s2pred::CompareDistance).
        fun setConservativeMaxDistance(maxDistance: S1ChordAngle): Unit {
            this.maxDistance = S2MinDistance(maxDistance.plusError(S2EdgeDistances.getUpdateMinDistanceMaxError(maxDistance)).successor())
        }

        // Versions of set_max_distance that take an S1Angle argument.  (Note that
        // these functions require a conversion, and that the S1ChordAngle versions
        // are preferred.)
        fun setMaxDistance(maxDistance: S1Angle) {

        }

        fun setInclusiveMaxDistance(maxDistance: S1Angle) {
            setInclusiveMaxDistance(S1ChordAngle(maxDistance))
        }
        fun setConservativeMaxDistance(maxDistance: S1Angle): Unit {
            setConservativeMaxDistance(S1ChordAngle(maxDistance))
        }

        // See S2ClosestPointQueryBaseOptions for documentation.
        fun setMaxError(maxError: S1Angle) {
            this.maxError = S1ChordAngle(maxError)
        }

        public override fun clone(): Options {
            return Options(getMaxResult(), maxDistance.clone(), maxError, region, useBruteForce)
        }

    }

    // Target subtype that computes the closest distance to a point.
    //
    // This class is also available as S2ClosestPointQuery<Data>::PointTarget.
    // (It is defined here to avoid depending on the "Data" template argument.)
    class S2ClosestPointQueryPointTarget(point: S2Point) : S2MinDistancePointTarget(point) {

        override fun maxBruteForceIndexSize(): Int {
            // Using BM_FindClosest (which finds the single closest point), the
            // break-even points are approximately X, Y, and Z points for grid,
            // fractal, and regular loop geometry respectively.
            //
            // TODO(ericv): Adjust using benchmarks.
            return 150;
        }


    }

    // Target subtype that computes the closest distance to an edge.
    //
    // This class is also available as S2ClosestPointQuery<Data>::EdgeTarget.
    // (It is defined here to avoid depending on the "Data" template argument.)
    class S2ClosestPointQueryEdgeTarget(a: S2Point, b: S2Point) : S2MinDistanceEdgeTarget(a, b) {

        override fun maxBruteForceIndexSize(): Int {
            // Using BM_FindClosestToEdge (which finds the single closest point), the
            // break-even points are approximately X, Y, and Z points for grid,
            // fractal, and regular loop geometry respectively.
            //
            // TODO(ericv): Adjust using benchmarks.
            return 100;
        }

    }

    // Target subtype that computes the closest distance to an S2Cell
    // (including the interior of the cell).
    //
    // This class is also available as S2ClosestPointQuery<Data>::CellTarget.
    // (It is defined here to avoid depending on the "Data" template argument.)
    class S2ClosestPointQueryCellTarget(cell: S2Cell) : S2MinDistanceCellTarget(cell) {

        override fun maxBruteForceIndexSize(): Int {
            // Using BM_FindClosestToCell (which finds the single closest point), the
            // break-even points are approximately X, Y, and Z points for grid,
            // fractal, and regular loop geometry respectively.
            //
            // TODO(ericv): Adjust using benchmarks.
            return 5 //50;
        }

    }

    // Target subtype that computes the closest distance to an S2ShapeIndex
    // (an arbitrary collection of points, polylines, and/or polygons).
    //
    // By default, distances are measured to the boundary and interior of
    // polygons in the S2ShapeIndex rather than to polygon boundaries only.
    // If you wish to change this behavior, you may call
    //
    //   target.set_include_interiors(false);
    //
    // (see S2MinDistanceShapeIndexTarget for details).
    //
    // This class is also available as S2ClosestPointQuery<Data>::ShapeIndexTarget.
    // (It is defined here to avoid depending on the "Data" template argument.)
    class S2ClosestPointQueryShapeIndexTarget(index: S2ShapeIndex) : S2MinDistanceShapeIndexTarget(index) {

        override fun maxBruteForceIndexSize(): Int {
            // For BM_FindClosestToSameSizeAbuttingIndex (which uses a nearby
            // S2ShapeIndex target of similar complexity), the break-even points are
            // approximately X, Y, and Z points for grid, fractal, and regular loop
            // geometry respectively.
            //
            // TODO(ericv): Adjust using benchmarks.
            return 30;
        }

    }

}
