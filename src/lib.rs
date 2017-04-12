// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! # Spade
//! Spade (SPAtial Data structurEs, obviously!) implements a few nifty data structures optimized for spatial access operations.
//! The first major data structure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R-tree)) for efficient nearest-neighbor and point lookup queries.
//! The second data structures implements a two dimensional delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)). 
//! The triangulation also implements various interpolation methods like nearest neighbor
//! interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)).
//! All classes are purely written in rust.
//! The data structures take either fixed size arrays or points of the nalgebra or cgmath crate as input.
//! There is also a [user guide](https://stoeoef.gitbooks.io/spade-user-manual/content/) available
//! to complement this reference.
//!
//! # Features
//! * An n-dimensional r*-tree: `spade::rtree::RTree`
//! * A 2D delaunay triangulation: `spade::delaunay::DelaunayTriangulation`
//!   * Supports integral and floating point coordinates as input
//!   * Uses exact predicates to avoid floating point rounding issues, see: `spade::kernels::FloatKernel`
//!   * Natural neighbor interpolation
//!   * Can be backed up by an r*-tree to improve performance when inserting randomly distributed points

#![warn(missing_docs)]

extern crate num;
extern crate cgmath;
extern crate nalgebra;
extern crate clamp;
extern crate smallvec;

#[cfg(test)]
extern crate approx;
#[cfg(test)]
extern crate rand;

#[cfg(test)]
mod testutils;

mod traits;
mod point_traits;
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
pub use point_traits::{PointN, TwoDimensional, ThreeDimensional};
