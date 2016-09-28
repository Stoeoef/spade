// Copyright 2016 The Spade Developers. For a full listing of the authors,
// refer to the Cargo.toml file at the top-level directory of this distribution.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! # Spade
//! Spade (SPAtial DatastructurEs, obviously!) implements a few nifty datastructures optimized for spatial access operations.
//! The first major datastructure is an n-dimensional r*-tree ([wikipedia](https://en.wikipedia.org/wiki/R-tree)) for efficient nearest-neighbor and point lookup queries.
//! The second datastructures implements a 2D delaunay triangulation ([wikipedia](https://en.wikipedia.org/wiki/Delaunay_triangulation)) backed by an r-tree for fast 
//! nearest neighbor lookup. The triangulation also implements natural neighbor interpolation ([wikipedia](https://en.wikipedia.org/wiki/Natural_neighbor)) 
//! which offers a smooth interpolation in a delaunay triangulation.
//! All classes are purely written in (safe) rust, the package currently supports vectors from the nalgebra and cgmath crates.
//!
//! # Features
//! * An n-dimensional R-Tree: `RTree`
//! * A 2D delaunay triangulation: `DelaunayTriangulation`
//!   * Supports integral and floating point vectors as input
//!   * Uses exact predicates to avoid floating point rounding issues, see: `FloatKernel`
//!   * Natural neighbor interpolation

#![warn(missing_docs)]
extern crate num;
extern crate cgmath;
extern crate nalgebra;
extern crate clamp;
#[cfg(test)]
extern crate approx;

extern crate rand;

pub mod testutils;

mod traits;
mod rtree;
mod misc;
mod boundingvolume;
mod bigvec;
mod exactpred;
mod kernels;

pub mod primitives;

mod delaunay;
mod planarsubdivision;

pub use planarsubdivision::{VertexHandle, AllEdgesIterator, AllVerticesIterator, FixedVertexHandle};
pub use rtree::*;
pub use delaunay::*;
pub use traits::*;
pub use boundingvolume::*;
pub use kernels::*;
