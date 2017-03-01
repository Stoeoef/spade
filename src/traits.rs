// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use cgmath as cg;

use num::{Signed, BigInt, BigRational, zero};
use num::rational::Ratio;
use boundingvolume::BoundingRect;
use bigvec::{AdaptiveInt};
use point_traits::{PointN, PointNExtensions, TwoDimensional};
use std::fmt::Debug;

/// Number types that can be used with spade.
/// 
/// Number types that can be used for spade's datastructures will need to
/// implement this trait. Can be an integer or floating point type.
/// Note that `copy` is not necessary to support big integers from the `num`
/// crate.
pub trait SpadeNum: Signed + Clone + Debug + PartialOrd { }

/// Floating point types that can be used with spade.
///
/// Used by all operations that require precise division.
pub trait SpadeFloat: SpadeNum + cg::BaseFloat { }

impl SpadeNum for i32 { }
impl SpadeNum for i64 { }
impl SpadeNum for f32 { }
impl SpadeNum for f64 { }

impl SpadeFloat for f32 { }
impl SpadeFloat for f64 { }

impl SpadeNum for BigInt { }
impl SpadeNum for BigRational { }
impl SpadeNum for AdaptiveInt { }
impl SpadeNum for Ratio<AdaptiveInt> { }

/// Describes objects that can be located by r-trees.
///
/// See the `primitives` module for some basic implementations which can also serve 
/// as useful examples for own implementations.
pub trait SpatialObject {
    /// The object's point type.
    type Point: PointN;

    /// Returns the object's minimal bounding rectangle.
    ///
    /// The minimal bounding rectangle is the smallest axis aligned rectangle that completely
    /// contains the object.
    /// <b>Note:</b> The rectangle must be as small as possible, otherwise some queries
    /// might fail.
    fn mbr(&self) -> BoundingRect<Self::Point>;

    /// Returns the squared euclidean distance from the object's contour.
    /// Returns a value samller than zero if the point is contained within the object.
    fn distance2(&self, point: &Self::Point) -> <Self::Point as PointN>::Scalar;

    /// Returns true if a given point is contained in this object.
    fn contains(&self, point: &Self::Point) -> bool {
        self.distance2(point) <= zero()
    }
}

/// An object that has a position.
/// 
/// Describes a point like object that has a well defined position.
/// Since this trait also implements `SpatialObject`, `HasPosition` can serve as a quick and
/// easy implementation of your own pointlike objects that can be inserted into r-trees
/// or delaunay triangulations:
///
/// ```
/// extern crate cgmath;
/// extern crate spade;
///
/// use cgmath::Point3;
/// use spade::HasPosition;
/// use spade::rtree::RTree;
///
/// struct MyPoint {
///   position: Point3<f32>,
///   my_data: i32,
/// }
/// impl HasPosition for MyPoint {
///    type Point = Point3<f32>;
///    fn position(&self) -> Point3<f32> {
///      self.position
///    }
/// }
/// fn main() {
///   let mut tree = RTree::new();
///   tree.insert(MyPoint { position: Point3::new(0.0, 1.0, 0.0), my_data: 42 });
///   assert_eq!(tree.lookup(&Point3::new(0.0, 1.0, 0.0)).unwrap().my_data, 42);
/// }
/// ```
pub trait HasPosition {
    /// The object's point type
    type Point: PointN;
    /// Return's the object's position.
    fn position(&self) -> Self::Point;
}

/// An object with a two dimensional position.
///
/// This trait is similar to `HasPosition`, but is only implemented for two dimensional points.
pub trait HasPosition2D: HasPosition where Self::Point: TwoDimensional { }

impl <V: HasPosition> HasPosition2D for V where V::Point: TwoDimensional { }

impl <V> HasPosition for V where V: PointN {
    type Point = V;
    fn position(&self) -> V {
        self.clone()
    }
}

impl <S>  SpatialObject for S where S: HasPosition {
    type Point = S::Point;

    fn mbr(&self) -> BoundingRect<S::Point> {
        BoundingRect::from_point(self.position())
    }

    fn distance2(&self, point: &S::Point) -> <S::Point as PointN>::Scalar {
        self.position().sub(point).length2()
    }

    fn contains(&self, point: &S::Point) -> bool {
        self.position() == *point
    }
}



