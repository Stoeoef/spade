[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)
![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/Stoeoef/spade/Continuous%20integration/master)

# spade

Delaunay triangulations for the rust ecosystem.

- 2D [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation), optionally backed by a hierarchy
 structure for improved nearest neighbor and insertion performance.
- Allows both incremental and bulk loading creation of triangulations
- Support for vertex removal
- 2D constrained Delaunay triangulation (CDT)
- [Delaunay refinement](https://en.wikipedia.org/wiki/Delaunay_refinement)
- Uses precise geometric predicates to prevent incorrect geometries due to rounding issues
- Supports extracting the Voronoi diagram

---------------------------
<img src="images/basic_voronoi.svg" width=60% style="margin-left: auto;margin-right:auto;display:block">

---------------------------

# Project goals

Project goals, in the order of their importance:

 1. Robustness - all data structures should behave correctly. An incorrect result, even if triggered only under rare circumstances, is not acceptable. This is why Spade uses a precise calculation kernel by default.
 2. Easy to use - favor an easy to use API over an API that exposes all bells and whistles.
 3. Performance - Delaunay triangulations are often a low level component of an application. Optimization in this area pays off greatly.
 4. Small footprint - Spade should be a sensible library to include in your project that doesn't require too many dependencies. Bigger dependencies will be feature gated when possible.

# Roadmap

For Spade 2.x:
 - Add back the removed interpolation methods (natural neighbor interpolation, #67)

For Spade 3:
 - Possibly base `spade` on `nalgebra` as underlying vector and matrix library. Not much else planned yet!

# Project state and contributing

Looking for co-maintainers! Projects with just a single maintainer can be a little unreliable due to the single point of failure. I would love to see this burden being distributed on more shoulders! This is less about *implementing* things but rather about general maintenance tasks, e.g. package updates, minor bug fixes, reviewing PRs, etc...

If you want to contribute, please consider opening an issue first. I'm happy to help out with any questions!

# Performance and comparison to other Delaunay crates

Refer to [the delaunay_compare readme](./delaunay_compare/README.md) for some benchmarks and a comparison with other crates.

# License
Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.