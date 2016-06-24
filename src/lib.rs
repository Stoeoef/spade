// Copyright 2016 The RTree Developers. For a full listing of the authors,
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

//! A rust implementation of two dimensional r*-trees
//!
//! [R-trees](https://en.wikipedia.org/wiki/R-tree) provide efficient nearest-neighbor searches for
//! a set of various objects. [R*-trees](https://en.wikipedia.org/wiki/R*_tree) (&quot;R-Star-Trees&quot;) 
//! are a common variant of r-trees and use more advanced heuristics to improve query performance. This
//! package implements only r*-trees at the moment, the original r-tree might be implemented in the future.
//! Instead of linear time complexity, r-trees yield logarithmic complexity
//! for look-up operations and nearest neighbor queries. Inserting into an r-tree runs in O(log(n)) time on average.
//! The main data structure is `RTree`. If you want to insert some simple geometric primitives,
//! consider to take a look at the `primitives` module. If your object is not among those, consider
//! implementing the `SpatialObject` trait.
//!
//! # Basic Example
//!
//! ```
//! use rtree::RTree;
//! use rtree::primitives::SimplePoint;
//!
//! let mut rtree = RTree::new();
//! // Insert two points
//! rtree.insert(SimplePoint::new([0.5, 0.5f32]));
//! rtree.insert(SimplePoint::new([1.0, 1.0f32]));
//!
//! if rtree.lookup([0.5, 0.5]).is_some() {
//!   println!("We'fe found a point at [0.5, 0.5]!");
//! }
//! 
//! let nearest = rtree.nearest_neighbor([1.5, 1.5]).unwrap();
//! println!("nearest neighbor at [1.5, 1.5]: {:?}", nearest);
//!
//! // Iterate over all elements
//! for point in rtree.iter() {
//!   println!("Found point: {:?}", point);
//! }
//! ```


extern crate cgmath;
extern crate num;
#[cfg(test)]
extern crate rand;
#[cfg(test)]
#[macro_use]
extern crate glium;

pub use rtree::*;
pub use traits::*;
pub use boundingvolume::*;

mod rtree;
mod traits;
mod boundingvolume;
pub mod primitives;

#[cfg(test)]
mod testutils;

mod misc;

