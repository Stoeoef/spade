// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! # Spade
//! Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.
//! The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R-tree)) for efficient nearest-neighbor and point lookup queries.
//! The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for fast 
//! nearest neighbor lookup. The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) 
//! which offers smooth interpolation in a delaunay triangulation.
//! All classes are purely written in (safe) rust.
//! The datastructures take either fixed size arrays or vectors of the nalgebra or cgmath crate as input.
//!
//! # Features
//! * An n-dimensional R-Tree: `RTree`
//! * A 2D delaunay triangulation: `DelaunayTriangulation`
//!   * Supports integral and floating point vectors as input
//!   * Uses exact predicates to avoid floating point rounding issues, see: `FloatKernel`
//!   * Natural neighbor interpolation

extern crate num;
extern crate cgmath;
extern crate nalgebra;
extern crate clamp;

#[cfg(test)]
extern crate approx;
#[cfg(test)]
extern crate rand;

#[cfg(test)]
mod testutils;

mod traits;
mod vector_traits;
mod misc;
mod boundingvolume;
mod bigvec;
mod exactpred;

pub mod delaunay;
pub mod kernels;
pub mod primitives;
pub mod rtree;

pub use traits::*;
pub use boundingvolume::*;
pub use vector_traits::{VectorN, TwoDimensional, ThreeDimensional};
