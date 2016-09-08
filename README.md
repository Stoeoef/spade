# spade

Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.

The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R*_tree)) for efficient nearest-neighbor and point lookup queries.

The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for fast nearest neighbor lookup.
The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) which allows for a smooth interpolation on the resulting grid.

All classes are purely written in rust, the package currently supports vectors from the [nalgebra](http://nalgebra.org/) and [cgmath](https://github.com/brendanzab/cgmath).

# Documentation
[Link to documentation](https://stoeoef.github.io/spade/)
