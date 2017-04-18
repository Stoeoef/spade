// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! A two dimensional delaunay triangulation.

mod dcel;
mod delaunay2d;
mod delaunay_locate;
mod cdt;
mod delaunay_basic;

pub use self::delaunay2d::*;
pub use self::cdt::ConstrainedDelaunayTriangulation;
pub use self::dcel::{FixedVertexHandle, FixedEdgeHandle, FixedFaceHandle,
                     VertexHandle, EdgeHandle, FaceHandle,
                     CCWIterator, ONextIterator};
pub use self::delaunay_locate::{DelaunayTreeLocate, DelaunayWalkLocate,
                                DelaunayLocateStructure};
#[allow(deprecated)]
pub use self::delaunay_locate::{TriangulationWalkLocate, RTreeDelaunayLocate};
