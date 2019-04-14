// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! A two dimensional Delaunay triangulation.

mod cdt;
mod dcel;
mod delaunay2d;
mod delaunay_basic;
mod delaunay_locate;
mod line_intersection_iterator;

pub use self::cdt::{CdtEdge, ConstrainedDelaunayTriangulation, FloatCDT};
pub use self::dcel::{
    CCWIterator, EdgeHandle, FaceHandle, FixedEdgeHandle, FixedFaceHandle, FixedVertexHandle,
    ONextIterator, VertexHandle,
};
pub use self::delaunay2d::*;
pub use self::delaunay_locate::{DelaunayLocateStructure, DelaunayTreeLocate, DelaunayWalkLocate};
#[allow(deprecated)]
pub use self::delaunay_locate::{RTreeDelaunayLocate, TriangulationWalkLocate};
