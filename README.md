[![Build Status](https://travis-ci.org/Stoeoef/spade.svg?branch=master)](https://travis-ci.org/Stoeoef/spade)
[![Docs](https://docs.rs/spade/badge.svg)](https://docs.rs/spade/)
[![Crates.io](https://img.shields.io/crates/v/spade.svg)](https://crates.io/crates/spade)
# spade

 * [Documentation](https://docs.rs/spade/)
 * [Using spade](#using-spade)
 * [Examples](#examples)
 * [Feedback](#feedback)
 * [Performance](#performance)

Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.

The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R*_tree)) for efficient nearest-neighbor and point lookup queries.

The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for faster insertion times and nearest neighbor lookup.
The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) which allows for a smooth interpolation on the resulting grid.

All structures are purely written in (safe) rust, the package currently supports vectors from the [nalgebra](http://nalgebra.org/) and [cgmath](https://github.com/brendanzab/cgmath) crates.

# Documentation
The documentation can be found under [docs.rs](https://docs.rs/spade/)

# Using spade
Add this to your `Cargo.toml`:
```
spade = "0.3.*"
```
# Feedback
Do you miss a feature? Many features may be easy to implement, the crate's main author might just have missed that use case. Please do post an issue on GitHub. Please do report bugs as well. If you're in for an adventure, pull requests of any kind are very welcome.

# Examples
_Note: If you have opened this on docs.rs, you won't see any images. Use the README.md on the GitHub page._
## R-Tree
This image shows the structure of an r*-tree with some points inserted in a circular pattern.
Points are shown as blue dots, the tree's directory nodes are displayed as boxes of different colors (depending on their depth in the tree).
Note that the implementation tries prevent any boxes from overlapping, resulting in faster query performance. You can find this example in `/examples/rtreevis`, run it with `cargo run`.

![An example R-Tree with a few inserted points](/images/rtree_demo.png?raw=true)

## Natural neighbor interpolation
The delaunay triangulation of a height field is shown with orange lines, the green grid shows the natural neighbor interpolated heights. In contrast to barycentric interpolation (the gray polygons), natural neighbor interpolation is smooth (C¹ differentiable, except for the data points themselves). You can find this example in `/examples/nninterpolation`, run it with `cargo run`.

![Delaunay triangulation with a grid showing interpolated values](/images/nninterpolation.png?raw=true)

## C¹ interpolants
Naive natural neighbor interpolation is not smooth at the data points. Other interpolants have been proposed which, given a gradient for each data point, will approximate those gradients and yield a completely smooth surface. Spade can approximate those gradients for you.
The following picture presents a side by side comparison of all implemented interpolants.
High resolutions: [left](https://raw.githubusercontent.com/Stoeoef/spade/master/images/InterpolationC0.png), [middle](https://raw.githubusercontent.com/Stoeoef/spade/master/images/InterpolationSibson.png), [right](https://raw.githubusercontent.com/Stoeoef/spade/master/images/InterpolationFarin.png)
![Side by side comparison of interpolation methods](/images/InterpolationMethods.png?raw=true)

# Performance
The following measurements were taken on an Intel Core i7-3517u.
![Performance of opererations on the r-tree implementation](/images/rtree_analysis.png?raw_true)

Insertion performance for various delaunay kernels:
![Performance of opererations on the r-tree implementation](/images/delaunay_analysis.png?raw_true)
