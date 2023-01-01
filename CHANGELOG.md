# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [2.1.0] - 2023-01-01

### Added
 - Added `crate::FloatTriangulations` for additional methods of triangulations over `f32` or `f64`
 - Added `FloatTriangulation::get_edges_in_rectangle`
 - Added `FloatTriangulation::get_edges_in_circle`
 - Added `FloatTriangulation::get_vertices_in_rectangle`
 - Added `FloatTriangulation::get_vertices_in_circle`
 - Added `index` function for any fixed and dynamic handle, returning the internal index of the element.
 - Added missing `DelaunayTriangulation::Clone` implementation
 - Added `DirectedEdgeHandle::positions`
 - Added `Triangulation::clear`

### Bugfixes
 - (breaking fix) `ConstrainedDelaunayTriangulation::can_add_constraint_edge` accidentally returned the wrong result (`true` when it should have returned `false` and vice versa) #75
 - Fixes a crash that could occur when inserting an element into a CDT (#78)

## [2.0.0] - 2022-01-29
 This release is focussed on API refinement, performance improvements and usability improvements.
  
 ## General API refinement
 Many features of Spade 1 showed too little benefit to legitimize keeping them. Spade 2.0 attempts to slim down the API surface to a few important functions and types and leaves out the more complicated bits to make it more simple to use.

### Removed features (planned to be re-introduced in an upcoming release)
 - Removed interpolation functions

### Removed features (not planned to re-introduced)
 - Removed support for integer coordinates. Only `f64` and `f32` can be used as scalar types.
 - Removed support for `cgmath` and `nalgebra` crate to slim down the crate dependencies.
 - Removed `PointN`. Using Spade now requires to convert any position into `Spade::Point2`.
 - Removed r-tree data structure (it has been moved into the `rstar` crate)
 - Removed geometric primitives (`SimpleCircle`, `SimpleEdge`, `SimpleTriangle`)
 - Removed support for custom calculation kernels. Spade now uses the precise calculation kernel by default (previously called `FloatKernel`).
 - Removed support for Delaunay triangulations backed by r-trees. Use the newly introduced `HierarchyHintGenerator` instead.

### Added features
 - Added basic support for Voronoi diagrams (refer to the documentation of struct `DelaunayTriangulation`)
 - Added struct `HierarchyHintGenerator`, a hint generator with small performance and memory overhead which is optimized for uniformly distributed data.
 - Improved conciseness of handle types (refer to the `handles` module for more information):
   - Introduced `UndirectedEdgeHandle` to refer to undirected edges of the triangulation
   - `FaceHandle` is now parameterized to denote if it refers to an inner face or to the single outer face
   - Added some methods on the individual handle types
   - Added iterators for all handle types (see `handles` module)
   - Fixed directed edge handles can now be `rev`ersed
 - Both `DelaunayTriangulation` and `CDT` now allow to also annotate undirected edges with arbitrary data
 - Added `Triangulation::convex_hull` which allows to directly iterate over the triangulation's convex hull
 - Added `Triangulation::bulk_load` for efficient bulk creation of triangulations.

### Changed
 - Performance improvements
  - The newly introduced `HierarchyHintGenerator` should outperform the r-tree based hint generator consistently for uniformly distributed data sets (2x - 3x times faster)
  - `LastUsedVertexHintGenerator` and its predecessor `DelaunayWalkLocate` are still comparable.
  - The new `Triangulation::bulk_load` will be more efficient for creating triangulations from a given set of vertices compared to inserting them incrementally.
 - Cleaned up crate dependencies: Spade 2.0 depends only on 4 other crates without any transitive dependencies. This should make it more viable to include Spade into other projects.
 - Renamed any `o_next` / `o_prev` function to `next` / `prev`.
 - Renamed `sym` to `rev`
 - Many other smaller changes - they are omitted from this changelog.

 ### Migration notes
 Keep in mind that not all features have been preserved before attempting to upgrade. Other than that, upgrading should be mostly consist of renaming:

  - `FloatDelaunayTriangulation` and `FloatCdt` have been replaced by `DelaunayTriangulation` and `ConstrainedDelaunayTriangulation`
  - `IntDelaunayTriangulation` is not supported anymore. Use `f64` or `f32` as scalar type instead if possible.
  - `DelaunayWalkLocate` has been renamed to `LastUsedVertexHintGenerator`
  - `DelaunayTreeLocate` has been replaced by `HierarchyHintGenerator`

## [1.8.2] - 2020-04-01
### Bugfixes
 - Removing elements from an rtree could leave the tree in an inconsistent state (#55). This made some nearest neighbor queries return incorrect results.

## [1.8.0] - 2019-06-20
### Changed
 - Bumped compatible nalgebra version to 0.18
 - Updated edition to 2018
 - Fixed all clippy findings
 - Cargo fmt'ed source code
### Bugfixes
 - SimpleTriangle now overwrites both Hash and PartialEq

## [1.7.0] - 2019-02-08
### Changed
 - Updated README.md
 - Merge #41: Use smallvec in barycentric interpolation
 - Bumped compatible cgmath version to 0.17.*
 - Bumped compatible nalgebra version to 0.17.*
  
## [1.6.0] - 2018-11-1
###
 - `SimpleEdge`, `SimpleCircle`, `BoundingRect` now derive `Clone`, `Copy`, `PartialEq`, `Eq`, `PartialOrd`, `Ord` and `Hash`
 - `CdtEdge` is now public - the type cannot be used, but may be necessary for method signatures containing CDTs.
 - Added some standard `derive`s for various other types that are applicable.
### Changed
 - Bumped compatible nalgebra version to 0.16.*
 - Bumped compatible num version to `>=0.1, <=0.2.*`. This may be a breaking change due to the way cargo resolves dependencies. An upgrade to num 0.2 is recommended.

## [1.5.1] - 2018-06-12
### Bugfixes
 - `nearest_neighbor` sometimes incorrectly returned `None`

## [1.5.0] - 2018-06-07
### Added
 - Added `RTree::nearest_neighbor_iterator`
### Changed
 - Improved performance of single nearest neighbor queries
 - Bumped compatible nalgebra version to 0.15.*
 - Bumped compatible cgmath version to 0.16.*

## [1.4.0] - 2018-02-25
### Added
 - Added cargo feature 'serde_serialize', rtrees, triangulations and primitives now support serialization with serde!
 - `RTree`s implement `Debug`
 - New constructor: `BoundingRect::from_points`
### Changed
 - Bumped compatible nalgebra version to 0.14.*
 ### Bugfixes
  - `BoundingRect` now implements `SpatialObject`

## [1.3.0] - 2018-01-06
### Deprecated
 - `spade::delaunay::DelaunayTriangulation::lookup` is deprecated, use `locate_vertex` instead
 - `spade::delaunay::DelaunayTriangulation::lookup_and_remove` is deprecated, use `locate_and_remove` instead
### Changed
 - Bumped compatible `cgmath` and `nalgebra` versions.
 - Spade's various kernels are not instantiable anymore and implement `Clone`.
### Added
 - Spade now implements constrained delaunay triangulations!
 - `locate` and `nearest_neighbor` now also work with degenerate triangulations
 - Added constrained triangulation example
 - Added bulk loading for r-trees
 - `DelaunayTriangulation` now implements `Clone`
 - New method for triangulations: `get_edge_from_neighbors` to get an existing edge from its two adjacent points.
  
### Bugfixes
 - Fixed `SimpleCircle` distance calculation for more than 2 dimensions.
 

## [1.2.0] - 2017-05-13
### Changed
  - Bumped compatible `cgmath` and `nalgebra` versions. Unfortunately, due to the way Cargo handles "external dependencies" (thus, dependencies whose types are part of spade's public API like `cgmath` and `nalgebra` Points), this must be considered a breaking change.
  - `FloatKernel` now works happily with `f32` coordinates. They will be casted to double precision before the kernel evaluates any query, thus, no performance gain results from this. Only the space requirements will differ.
  
## [1.1.0] - 2017-04-12
### Deprecated
  - `spade::delaunay::RTreeDelaunayLocate<T>` is deprecated, use `spade::delaunay::DelaunayTreeLocate<T>` instead
  - `spade::delaunay::TriangulationWalkLocate<T>` is deprecated, use `spade::delaunay::DelaunayWalkLocate` instead (without any type argument)
### Changed
  - Insertion into a delaunay triangulation now uses `SmallVec` from the `smallvec` crate for better performance.
  - Improved interpolation performance - natural neighbor interpolation methods will be significantly faster now.
### Added
  - Added struct `spade::delaunay::DelaunayWalkLocate`. The struct will now keep track of the last query made and use it for the next query. A query can be an interpolation, lookup, locate or nearest neighbor query, insertion will also update the hint to the inserted vertex. This means: Subsequent queries that are close to each other will run in O(1) without the need of giving an explicit hint. This behaviour was only implemented for insertion until now.
  - Added method `spade::delaunay::DelaunayTriangulation::nearest_neighbor`
  - Added method `spade::delaunay::DelaunayTriangulation::locate_vertex`
  - Added method `spade::primitives::SimpleEdge::length2`
  - Added an interpolation benchmark
  
## [1.0.0] - 2017-03-02
A lot has changed for the 1.0. release, only larger changes are shown.

### Changed
  - Changed project license from Apache 2.0 to dual MIT / Apache 2.0
  - `VectorN` renamed to `PointN`
  - Bumped supported nalgebra version to `0.11.*`
  - `Point2`, `Point3` ... from cgmath and nalgebra crates *do* implement `VectorN`
  - `Vector2`, `Vector3` ... from cgmath and nalgebra crates *do not* implement `VectorN`
  - Moved all kernels from the root to `spade::kernels`
  - Moved delaunay triangulation and handle types to `spade::delaunay`
  - Moved the rtree to `spade::rtree::RTree`
  - Renamed `DelaunayTriangulation::lookup_in_triangulation` to `locate`
  - Renamed `DelaunayTriangulation::handle(..)` to `vertex(..)`
  
### Added
 - Added support for vertex removal in delaunay triangulations!
 - Added `DelaunayTriangulation::barycentric_interpolation(..)`
 - Added support for different lookup methods for delaunay triangulation (see `DelaunayLookupStructure`)
 - Added a user guide! Check it out [here](https://stoeoef.gitbooks.io/spade-user-manual/content/)
 - `locate` and `lookup` can now be used with a hint (see `{insert|locate}_with_hint`)
 - Added type shorthands `spade::delaunay::{Float|Int}Triangulation` and quick creation methods `with_tree_lookup` and `with_walk_lookup`
 - Added struct `EdgeHandle` and type `FixedEdgeHandle`.
 - Added struct `FaceHandle` and type `FixedFaceHandle`.
 - Added methods `DelaunayTriangulation::edge(..)` and `::face(..)`
 - Added `DelaunayTriangulation::infinite_face()`
 - Added `DelaunayTriangulation::is_degenerate()`
 - Added `DelaunayTriangulation::num_{edges|triangles|faces}()`
 - `CCWIterator` and `ONextIterator` now both implement `DoubleEndedIterator`
 
### Removed
 - Removed support for pointer like types (e.g. inserting `Box<Point2<_>>`) to simplify type signatures.
 - Removed `DelaunayTriangulation::lookup_in_circle(..)` and `lookup_in_rect(..)`. These methods will likely be added again in a later release.
 
## [0.3.0] - 2016-12-10
### Changed
  - `VectorN` trait trimmed down, now contains only a small set of functions.
  
### Added
  - `[S; 2]`, `[S; 3]` and `[S; 4]` now implement `VectorN` and `TwoDimensional` or `ThreeDimensional` when appropriate. This change allows to insert fixed size arrays that encode positions directly into `RTree` or `DelaunayInterpolation`.
  
### Removed
  - Removed dependencies on crate `rand` and crate `noise`
  
## [0.2.1] - 2016-10-11
### Changed
  - Function signatures of `nn_interpolation_c1_*` slightly modified.
  - Sibson's c1 interpolant now comes with a flatness factor
### Fixes
  - Wrong documentation link in crate description
  - Fixed signatures of `estimate_normal` and `estimate_gradient`
  
## [0.2.0] - 2016-10-08
### Added
  - `DelaunayTriangulation`: `estimate_normal`, `estimate_normals`, `estimate_gradient`, `estimate_gradients`
  - Added Sibson's c1 interpolant, `DelaunayTriangulation::nn_interpolation_c1_sibson`
  - Added Farin's c1 interpolant, `DelaunayTriangulation::nn_interpolation_c1_farin`
  - trait `ThreeDimensional`
  
### Changed
 - Type signatures of `RTree` and `DelaunayTriangulation` have now an additional parameter, `B`.
 This allows to insert pointer like objects (that is, an object `B: Borrow<T>`) into the tree.

## [0.1.1] - 2016-09-28
### Added
 - Documentaion to all functions and types intended for public use
 - `RTree::lookup_mut(..)`
 - `RTree::contains(..)`
 - `DelaunayTriangulation::handle_mut(..)`
 - `DelaunayTriangulation::lookup_mut(..)`
 - `DelaunayKernel::point_on_edge(..)`
 - `SimpleTriangle::nearest_point_on_edge(..)`
 - types `TwoDimensional`, `HasPosition2D`
 
### Removed
 - `SimpleEdge::point_on_edge(..)`
 - `SimpleTriangle::is_ordered_ccw(..)`

### Fixed
 - Potential crashes when inserting points into a `DelaunayTriangulation`,
 even though `FloatKernel` was used.

### Changed
 - cgmath dependency bumped to 0.12.*
 - DelaunayTriangulations and some primitives now will only work with two
 dimensional coordinates. Using higher dimensions actually yielded unspecified
 results.

## 0.1.0 - 2016-09-23
Initial commit

[2.1.0]: https://github.com/Stoeoef/spade/compare/v2.0.0...v2.1.0

[2.0.0]: https://github.com/Stoeoef/spade/compare/v1.8.2...v2.0.0

[1.8.2]: https://github.com/Stoeoef/spade/compare/v1.8.0...v1.8.2

[1.8.0]: https://github.com/Stoeoef/spade/compare/v1.7.0...v1.8.0

[1.7.0]: https://github.com/Stoeoef/spade/compare/v1.6.0...v1.7.0

[1.6.0]: https://github.com/Stoeoef/spade/compare/v1.5.1...v1.6.0

[1.5.1]: https://github.com/Stoeoef/spade/compare/v1.5.0...v1.5.1

[1.5.0]: https://github.com/Stoeoef/spade/compare/v1.4.0...v1.5.0

[1.4.0]: https://github.com/Stoeoef/spade/compare/v1.3.0...v1.4.0

[1.3.0]: https://github.com/Stoeoef/spade/compare/v1.2.0...v1.3.0

[1.2.0]: https://github.com/Stoeoef/spade/compare/v1.1.0...v1.2.0

[1.1.0]: https://github.com/Stoeoef/spade/compare/v1.0.0...v1.1.0

[1.0.0]: https://github.com/Stoeoef/spade/compare/v0.3.0...v1.0.0

[0.3.0]: https://github.com/Stoeoef/spade/compare/v0.2.1...v0.3.0

[0.2.1]: https://github.com/Stoeoef/spade/compare/v0.2.0...v0.2.1

[0.2.0]: https://github.com/Stoeoef/spade/compare/v0.1.1...v0.2.0

[0.1.1]: https://github.com/Stoeoef/spade/compare/v0.1.0...v0.1.1
