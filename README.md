# spade

Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.

The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R*_tree)) for efficient nearest-neighbor and point lookup queries.

The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for faster insertion times and nearest neighbor lookup.
The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) which allows for a smooth interpolation on the resulting grid.

All classes are purely written in rust, the package currently supports vectors from the [nalgebra](http://nalgebra.org/) and [cgmath](https://github.com/brendanzab/cgmath) crates.

# Documentation
[Link to documentation](https://stoeoef.github.io/spade/)

# Examples
## R-Tree
This image shows the structure of an r*-tree with some points inserted in a circular pattern.
Points are shown as blue dots, the tree's directory nodes are displayed as boxes of different colors (depending on their depth in the tree).
Note that the implementation tries prevent any boxes from overlapping, resulting in faster query performance.
![An example R-Tree with a few inserted points](/images/rtree_demo.png?raw=true)

## Natural neighbor interpolation
The delaunay triangulation of a height field is shown with orange lines, the green grid shows the natural neighbor interpolated heights. In contrast to barycentric interpolation (the gray polygons), natural neighbor interpolation is smooth (C¹ differentiable, except for the data points themselves).
![Delaunay triangulation with a grid showing interpolated values](/images/nninterpolation.png?raw=true)

# Performance
The following measurements were taken on an Intel Core i7-3517u.
![Performance of opererations on the r-tree implementation](/images/rtree_analysis.png?raw_true)

Insertion performance for various delaunay kernels:
![Performance of opererations on the r-tree implementation](/images/delaunay_analysis.png?raw_true)
