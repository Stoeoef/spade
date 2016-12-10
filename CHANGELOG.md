# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

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

[0.3.0] https://github.com/Stoeoef/spade/compare/v0.2.1...v0.3.0
[0.2.1] https://github.com/Stoeoef/spade/compare/v0.2.0...v0.2.1
[0.2.0] https://github.com/Stoeoef/spade/compare/v0.1.1...v0.2.0
[0.1.1] https://github.com/Stoeoef/spade/compare/v0.1.0...v0.1.1
