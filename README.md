![Build Status](https://github.com/Enovea/s2-geometry-kotlin/actions/workflows/maven.yml/badge.svg)


# Overview

S2 is a library for spherical geometry that aims to have the same robustness,
flexibility, and performance as the best planar geometry libraries.

This is a library for manipulating geometric shapes. Unlike many geometry
libraries, S2 is primarily designed to work with _spherical geometry_, i.e.,
shapes drawn on a sphere rather than on a planar 2D map. (In fact, the name S2
is derived from the mathematical notation for the unit sphere *S²*.) This makes
it especially suitable for working with geographic data.

More details about S2 in general are available on the S2 Geometry Website
[s2geometry.io](https://s2geometry.io/).

## Scope

The library provides the following:

*   Representations of angles, intervals, latitude-longitude points, unit
    vectors, and so on, and various operations on these types.

*   Geometric shapes over the unit sphere, such as spherical caps ("discs"),
    latitude-longitude rectangles, polylines, and polygons. These are
    collectively known as "regions".

*   A hierarchical decomposition of the sphere into regions called "cells". The
    hierarchy starts with the six faces of a projected cube and recursively
    subdivides them in a quadtree-like fashion.

*   Robust constructive operations (e.g., union) and boolean predicates (e.g.,
    containment) for arbitrary collections of points, polylines, and polygons.

*   Fast in-memory indexing of collections of points, polylines, and polygons.

*   Algorithms for measuring distances and finding nearby objects.

*   Robust algorithms for snapping and simplifying geometry (with accuracy and
    topology guarantees).

*   A collection of efficient yet exact mathematical predicates for testing
    relationships among geometric objects.

*   Support for spatial indexing, including the ability to approximate regions
    as collections of discrete "S2 cells". This feature makes it easy to build
    large distributed spatial indexes.

On the other hand, the following are outside the scope of S2:

*   Planar geometry.

*   Conversions to/from common GIS formats.

### Robustness

What do we mean by "robust"?

In the S2 library, the core operations are designed to be 100% robust. This
means that each operation makes strict mathematical guarantees about its output,
and is implemented in such a way that it meets those guarantees for all possible
valid inputs. For example, if you compute the intersection of two polygons, not
only is the output guaranteed to be topologically correct (up to the creation of
degeneracies), but it is also guaranteed that the boundary of the output stays
within a user-specified tolerance of true, mathematically exact result.

Robustness is very important when building higher-level algorithms, since
unexpected results from low-level operations can be very difficult to handle. S2
achieves this goal using a combination of techniques from computational
geometry, including *conservative error bounds*, *exact geometric predicates*,
and *snap rounding*.

The implementation attempts to be precise both in terms of mathematical
definitions (e.g. whether regions include their boundaries, and how degeneracies
are handled) and numerical accuracy (e.g. minimizing cancellation error).

Note that the intent of this library is to represent geometry as a mathematical
abstraction. For example, although the unit sphere is obviously a useful
approximation for the Earth's surface, functions that are specifically related
to geography are not part of the core library (e.g. easting/northing
conversions, ellipsoid approximations, geodetic vs. geocentric coordinates,
etc).

For an analogous library in C++, see https://github.com/google/s2geometry, in
Go, see https://github.com/golang/geo, and Python, see
https://github.com/google/s2geometry/tree/master/src/python

# Status of the Go Library

This library is principally a port of the
[C++ S2 library](https://github.com/google/s2geometry), adapting to Kotlin idioms
where it makes sense. We detail the progress of this port below relative to that
C++ library.

## ℝ¹ - One-dimensional Cartesian coordinates

Full parity with C++.

## ℝ² - Two-dimensional Cartesian coordinates

Full parity with C++.

## ℝ³ - Three-dimensional Cartesian coordinates

Full parity with C++.

## S¹ - Circular Geometry

Full parity with C++.

## S² - Spherical Geometry

Full parity with C++.

## Encode/Decode

**Not Started Yet.**
Encoding and decoding of S2 types is fully implemented and interoperable with
C++ and Go.
