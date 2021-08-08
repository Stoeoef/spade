[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)

# spade

 * [Documentation](https://docs.rs/spade/)
 * [Using spade](#using-spade)
 * [Examples](#examples)
 * [Project state](#project-state)
 * [Performance](#performance)
 * [License](#license)

Delaunay triangulations for the rust ecosystem.

- 2D [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation), optionally backed by a hierarchy
 structure for improved nearest neighbor and insertion performance.
- 2D constrained Delaunay triangulation (CDT)
- Precise calculation kernel out of the box to prevent incorrect geometries due to rounding issues
- Supports extracting the Voronoi diagram

# Project goals

Project goals, in the order of their importance:

 1. Robustness - all data structures should behave correctly. An incorrect result, even if triggered only under rare circumstances, is not acceptable. This is why Spade uses a precise calculation kernel by default.
 2. Easy to use - favor an easy to use API over an API that exposes all bells and whistles.
 3. Performance - Delaunay triangulations are often a low level component of an application. Optimization in this area pays off greatly.
 4. Small footprint - Spade should be a sensible library to include in your project that doesn't require too many dependencies. Bigger dependencies will be feature gated.

# Project state and contributing
Looking for co-maintainers! I don't have the best track record when it comes to being an open-source maintainer and would love to see this burden being distributed on more shoulders! If you are interested in seeing a project like this being around for a long while, please get in touch.

If you want to contribute, please consider opening an issue first. I'm happy to help out with any questions!


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
