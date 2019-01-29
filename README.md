[![Build Status](https://travis-ci.org/Stoeoef/spade.svg?branch=master)](https://travis-ci.org/Stoeoef/spade)
[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)
# spade

 * [Documentation](https://docs.rs/spade/)
 * [Using spade](#using-spade)
 * [Examples](#examples)
 * [Project state](#project-state)
 * [Performance](#performance)
 * [License](#license)

Spade (SPAtial DatastructurEs, obviously!) implements a few nifty data structures for spatial access operations:

- An n-dimensional [r*-tree](https://en.wikipedia.org/wiki/R*_tree) for efficient nearest-neighbor and point lookup queries. Note that a faster successor is available with the [rstar crate](https://crates.io/crates/rstar).
- 2D [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation), optionally backed by an r-tree for faster insertion and nearest neighbor lookup
- 2D constrained Delaunay triangulation (CDT)

Some other noteworthy features:
- [natural neighbor interpolation](https://en.wikipedia.org/wiki/Natural_neighbor) on this triangulation
- Precise and adaptive calculation methods to avoid rounding issues
- supports [serde](https://crates.io/crates/serde) with the `serde_serialize` feature

All structures are purely written in rust, the package currently supports vectors from the [nalgebra](http://nalgebra.org/) and [cgmath](https://github.com/brendanzab/cgmath) packages. However, using these
packages is not required.

# Compatibility note
Spade complies with semantic versioning, and since it is past its 1.0 version, current minor version changes will be backward compatible. However, due to the way cargo resolves dependencies, there might be issues when using spade combined with cgmath or nalgebra: every time spade updates these libraries, the using code must be update too, even if spade would still work happily with an older version. To avoid this, consider switching to fixed size arrays as points, implementing your own point type or do some creative hacking into your cargo.lock to force cargo into using the correct cgmath / nalgebra version.

# Documentation
The documentation can be found on [docs.rs](https://docs.rs/spade/).
There is also a [user guide](https://stoeoef.gitbooks.io/spade-user-manual/content/) available.

# Project state
Spade is being passively maintained, please file all bugs that you can find! However, I don't plan any major release at the moment.
Spade's r-tree has been split off into the smaller [rstar crate](https://crates.io/crates/rstar) which is the recommended replacement. rstar compiles faster, runs faster and is more actively developed. 

# Examples
## R-Tree
This image shows the structure of an r*-tree with some points inserted in a circular pattern.
Points are shown as blue dots, the tree's directory nodes are displayed as boxes of different colors (depending on their depth in the tree).
Note that the implementation tries prevent any boxes from overlapping, resulting in faster query performance. You can find this example in `/examples/interactivedemo`, run it with `cargo run rtree`.

![An example R-Tree with a few inserted points](/images/rtree_demo.png?raw=true)

## CDT
CDT's are usual Delaunay triangulations with a few "fixed" edges:
![An example of a CDT](/images/cdt_demo.png?raw=true)

## Interpolation
The [user guide](https://stoeoef.gitbooks.io/spade-user-manual/) has a an [own chapter](https://stoeoef.gitbooks.io/spade-user-manual/content/interpolation.html) about interpolation, along with some nice images.
An example showcasing spade's interpolation features can be found in `/examples/nninterpolation`, run it with `cargo run`.

# Performance
Neither spade's triangulation nor r-tree were optimized to be exceptionally fast. The library focussed on a rich and high quality feature set and an early 1.0 release at the cost of performance. Compared to any GCed language, spade's performance is likely better (feel free to contribute a benchmark!) but it is _by far not_ in the same ballpark as its C/C++ contenders like CGAL. However, for many use cases, spade will hopefully be fast enough and a viable rust only alternative.

The [user guide](https://stoeoef.gitbooks.io/spade-user-manual/content/triangulation-performance.html) contains detailed graphs, benchmarks and more information about the delaunay triangulation's performance.

# License
Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
