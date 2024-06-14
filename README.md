[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)
[![License](https://img.shields.io/crates/l/spade)](https://crates.io/crates/spade)
![GitHub Workflow Status (branch)](https://img.shields.io/github/actions/workflow/status/Stoeoef/spade/spade_actions.yaml?branch=master)
[![dependency status](https://deps.rs/repo/github/Stoeoef/spade/status.svg)](https://deps.rs/repo/github/Stoeoef/spade)

# spade

Delaunay triangulations for the rust ecosystem.

- 2D [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation), optionally backed
  by a hierarchy structure for improved nearest neighbor and insertion performance.
- Allows both incremental and bulk loading creation of triangulations
- Support for vertex removal
- 2D constrained Delaunay triangulation (CDT)
- [Delaunay refinement](https://en.wikipedia.org/wiki/Delaunay_refinement)
- Uses precise geometric predicates to prevent incorrect geometries due to rounding issues
- Supports extracting the Voronoi diagram
- Natural neighbor interpolation

---------------------------

<img src="images/basic_voronoi.svg" width=60% style="margin-left: auto;margin-right:auto;display:block;object-fit:cover" alt="A Delaunay triangulation and the dual Voronoi Graph">

---------------------------

# Project goals

Project goals, in the order of their importance:

1. Robustness - all data structures should behave correctly. An incorrect result, even if triggered only under rare
   circumstances, is not acceptable. This is why Spade uses a precise calculation kernel by default.
2. Easy to use - favor an easy-to-use API over an API that exposes all bells and whistles.
3. Performance - Delaunay triangulations are often a low level component of an application. Optimization in this area
   pays off greatly.
4. Small footprint - Spade should be a sensible library to include in your project that doesn't require too many
   dependencies. Bigger dependencies will be feature gated when possible.

# Roadmap

For Spade 3:

- Possibly API simplification by un-supporting non-f64 outputs.

# Performance and comparison to other crates

Refer to [the delaunay_compare readme](./delaunay_compare/README.md) for some benchmarks and a comparison with other
crates.

# License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.