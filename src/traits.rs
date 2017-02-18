// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use nalgebra as na;
use cgmath as cg;

use num::{Signed, BigInt, BigRational, zero};
use num::rational::Ratio;
use boundingvolume::BoundingRect;
use bigvec::{AdaptiveInt};
use vector_traits::{VectorN, VectorNExtensions, TwoDimensional};
use std::fmt::Debug;

/// Scalar that can be used for spade's datastructures will need to implement
/// this trait. Can be an integer or floating point type.
/// Note that `copy` is not necessary.
pub trait SpadeNum: Signed + Clone + Debug + PartialOrd { }
/// Trait for `SpadeNum`s that are also a floating point number.
/// Used by all operations that require precise division.
pub trait SpadeFloat: SpadeNum + cg::BaseFloat + na::BaseFloat { }

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
    /// The object's vector type.
    type Vector: VectorN;

    /// Returns the object's minimal bounding rectangle.
    ///
    /// The minimal bounding rectangle is the smallest axis aligned rectangle that completely
    /// contains the object.
    /// <b>Note:</b> The rectangle must be as small as possible, otherwise some queries
    /// might fail.
    fn mbr(&self) -> BoundingRect<Self::Vector>;

    /// Returns the squared euclidean distance from the object's contour.
    /// Returns a value samller than zero if the point is contained within the object.
    fn distance2(&self, point: &Self::Vector) -> <Self::Vector as VectorN>::Scalar;

    /// Returns true if a given point is contained in this object.
    fn contains(&self, point: &Self::Vector) -> bool {
        self.distance2(point) <= zero()
    }
}

/// An object with a well defined location.
/// Since this trait also implements `SpatialObject`, `HasPosition` can serve as a quick and
/// easy implementation of your own pointlike objects that can be inserted into r-trees:
///
/// ```
/// extern crate cgmath;
/// extern crate spade;
///
/// use cgmath::Vector3;
/// use spade::HasPosition;
/// use spade::rtree::RTree;
///
/// struct MyPoint {
///   position: Vector3<f32>,
///   my_data: i32,
/// }
/// impl HasPosition for MyPoint {
///    type Vector = Vector3<f32>;
///    fn position(&self) -> Vector3<f32> {
///      self.position
///    }
/// }
/// fn main() {
///   let mut tree = RTree::new();
///   tree.insert(MyPoint { position: Vector3::new(0.0, 1.0, 0.0), my_data: 42 });
///   assert_eq!(tree.lookup(&Vector3::new(0.0, 1.0, 0.0)).unwrap().my_data, 42);
/// }
/// ```
pub trait HasPosition {
    /// The object's vector type
    type Vector: VectorN;
    /// Return's the object's position.
    fn position(&self) -> Self::Vector;
}

/// Like `HasPosition`, but is only implemented for two dimensional positions.
pub trait HasPosition2D: HasPosition where Self::Vector: TwoDimensional { }

impl <V: HasPosition> HasPosition2D for V where V::Vector: TwoDimensional { }

impl <V> HasPosition for V where V: VectorN {
    type Vector = V;
    fn position(&self) -> V {
        self.clone()
    }
}

impl <S>  SpatialObject for S where S: HasPosition {
    type Vector = S::Vector;

    fn mbr(&self) -> BoundingRect<S::Vector> {
        BoundingRect::from_point(self.position())
    }

    fn distance2(&self, point: &S::Vector) -> <S::Vector as VectorN>::Scalar {
        self.position().sub(point).length2()
    }

    fn contains(&self, point: &S::Vector) -> bool {
        self.position() == *point
    }
}



