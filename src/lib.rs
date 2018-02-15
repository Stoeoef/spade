// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! # Spade
//! Spade (SPAtial Data structurEs, obviously!) implements a few nifty data structures optimized for spatial access operations.
//!
//! # Features
//! * An n-dimensional r*-tree: `spade::rtree::RTree`
//! * A 2D Delaunay triangulation: `spade::delaunay::DelaunayTriangulation`
//!   * Supports integral and floating point coordinates as input
//!   * Uses exact predicates to avoid floating point rounding issues, see `spade::kernels::FloatKernel`
//!   * Natural neighbor interpolation
//!   * Can be backed up by an r*-tree to improve performance when inserting randomly distributed points
//!   * Supports vertex removal
//! * A 2D constrained Delaunay triangulation: `spade::delaunay::ConstrainedDelaunayTriangulation`
//!
//! # Supported point types
//! Spade works well with points from the `nalgebra` and `cgmath` packages. Also, fixed size arrays of size 2, 3 and 4 are
//! supported. Also own vector types can be defined.
//! Please note that, due to the way cargo resolves dependencies, there might be issues when using spade combined with cgmath 
//! or nalgebra: every time spade updates these libraries, the using code must be update too, even if spade would still work
//!  with the older version. To avoid this, consider switching to fixed size arrays as points until 
//! [public / private dependencies make their way into cargo](https://github.com/rust-lang/rust/issues/44663).



#![warn(missing_docs)]

extern crate num;
extern crate cgmath;
extern crate nalgebra;
extern crate clamp;
extern crate smallvec;
extern crate pdqselect;

#[cfg(test)]
extern crate approx;
#[cfg(test)]
extern crate rand;

#[cfg(test)]
mod testutils;

mod traits;
mod point_traits;
mod misc;
mod boundingrect;
mod bigvec;
mod exactpred;

pub mod delaunay;
pub mod kernels;
pub mod primitives;
pub mod rtree;

pub use traits::*;
pub use boundingrect::*;
pub use point_traits::{PointN, TwoDimensional, ThreeDimensional};
