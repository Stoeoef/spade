[![Build Status](https://travis-ci.org/Stoeoef/spade.svg?branch=master)](https://travis-ci.org/Stoeoef/spade)
[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)
# spade

 * [Documentation](https://docs.rs/spade/)
 * [Using spade](#using-spade)
 * [Examples](#examples)
 * [Feedback](#feedback)
 * [Performance](#performance)
 * [License](#license)

Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.

The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R*_tree)) for efficient nearest-neighbor and point lookup queries.

The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for faster insertion times and nearest neighbor lookup.
The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) which allows for a smooth interpolation on the resulting grid.

All structures are purely written in (safe) rust, the package currently supports vectors from the [nalgebra](http://nalgebra.org/) and [cgmath](https://github.com/brendanzab/cgmath) crates.

# Project state
The next goal is to release a 1.0 version, there's already an appropriate milestone for it. Check its progress to see if we're there yet.
Also, there is a small user guide being produced right now which should be a helpful ressource when starting to use the library.

The algorithms and features for the 1.0 version are already there, the public interface will change though (mostly due to renaming / repackaging).

# Documentation
The documentation can be found under [docs.rs](https://docs.rs/spade/).
There is also a (yet unfinished) [user guide](https://stoeoef.gitbooks.io/spade-user-manual/content/) available.

# Using spade
Add this to your `Cargo.toml`:
```
spade = "0.3.*"
```
# Feedback and contributing
Do you miss a feature? Many features may be easy to implement, I might just have missed that use case. Please do post an issue on GitHub. If you think there's a feature missing and you are interested to implement it yourself, don't be shy to mention it - I'd be happy to help you getting started. Just post an appropriate issue on GitHub.

# Examples
_Note: If you have opened this on docs.rs, you won't see any images. Use the README.md on the GitHub page._
## R-Tree
This image shows the structure of an r*-tree with some points inserted in a circular pattern.
Points are shown as blue dots, the tree's directory nodes are displayed as boxes of different colors (depending on their depth in the tree).
Note that the implementation tries prevent any boxes from overlapping, resulting in faster query performance. You can find this example in `/examples/interactivedemo`, run it with `cargo run`.

![An example R-Tree with a few inserted points](/images/rtree_demo.png?raw=true)

## Interpolation
The [user guide](https://stoeoef.gitbooks.io/spade-user-manual/) has a an [own chapter](https://stoeoef.gitbooks.io/spade-user-manual/content/interpolation.html) about interpolation, along with some nice images.
An example showcasing spade's interpolation features can be found in `/examples/nninterpolation`, run it with `cargo run`.

# Performance
The [user guide](https://stoeoef.gitbooks.io/spade-user-manual/content/triangulation-performance.html) contains detailed graphs and information about the delaunay triangulation's performance.

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
