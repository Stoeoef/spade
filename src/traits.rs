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

use cgmath as cg;
use num::{Signed, BigInt, BigRational, zero};
use num::rational::Ratio;
use misc::{min_inline, max_inline};
use boundingvolume::BoundingRect;
use bigvec::{AdaptiveInt};
use std::ops::{Add, Sub, Index, IndexMut, Div, Mul};
use std::fmt::Debug;
use nalgebra as na;
use nalgebra::{Repeat};

pub trait SpadeNum: Signed + Clone + Debug + PartialOrd { }
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

/// Abstraction over vectors in different dimensions.
pub trait VectorN where Self: Clone,
Self: Add<Self, Output=Self>,
Self: Sub<Self, Output=Self>,
Self: Div<<Self as VectorN>::Scalar, Output=Self>,
Self: Mul<<Self as VectorN>::Scalar, Output=Self>,
Self: Index<usize, Output=<Self as VectorN>::Scalar>,
Self: IndexMut<usize, Output=<Self as VectorN>::Scalar>,
Self: Debug,
Self: PartialEq {
    type Scalar: SpadeNum;

    /// Creates a new vector with all compoenents set to a certain value.
    fn from_value(value: Self::Scalar) -> Self;

    /// Creates a new vector with all components initialized to zero.
    fn new() -> Self {
        Self::from_value(zero())
    }

    /// The (fixed) number of dimensions of this vector trait.
    fn dimensions() -> usize;

    /// Applies a binary operation component wise.
    fn component_wise<F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar>(&self, rhs: &Self, f: F) -> Self {
        let mut result = self.clone();
        for i in 0 .. Self::dimensions() {
            result[i] = f(self[i].clone(), rhs[i].clone());
        }
        result
    }

    /// Maps an unary operation to all compoenents.
    fn map<F: Fn(Self::Scalar) -> O::Scalar, O: VectorN>(&self, f: F) -> O {
        let mut result = O::new();
        for i in 0 .. Self::dimensions() {
            result[i]  = f(self[i].clone());
        }
        result
    }

    /// Returns a new vector containing the minimum values of this and another vector (componentwise)
    fn min_vec(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| min_inline(l, r))
    }

    /// Returns a new vector containing the maximum values of this and another vector (componentwise)
    fn max_vec(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| max_inline(l, r))
    }

    /// Fold operation over all vector components.
    fn fold<T, F: Fn(T, Self::Scalar) -> T>(&self, mut acc: T, f: F) -> T {
        for i in 0 .. Self::dimensions() {
            acc = f(acc, self[i].clone());
        }
        acc
    }

    /// Checks if a property holds for all components.
    fn all_comp_wise<F: Fn(Self::Scalar, Self::Scalar) -> bool>(&self, rhs: &Self, f: F) -> bool {
        for i in 0 .. Self::dimensions() {
            if !f(self[i].clone(), rhs[i].clone()) {
                return false;
            }
        }
        true
    }

    /// Returns the vector's dot product.
    fn dot(&self, rhs: &Self) -> Self::Scalar {
        self.component_wise(rhs, |l, r| l * r).fold(zero(), |acc, val| acc + val)
    }

    /// Returns the vector's squared length.
    fn length2(&self) -> Self::Scalar {
        self.dot(&self)
    }
}

/// Only implemented by two dimensional vectors.
/// Some datastructures won't work for if 3 or
/// more dimensional vectors are used, this trait
/// will ensure that they all have the correct
/// dimension.
pub trait TwoDimensional : VectorN { }

impl <S: SpadeNum + cg::BaseNum> TwoDimensional for cg::Vector2<S> { }
impl <S: SpadeNum + na::BaseNum> TwoDimensional for na::Vector2<S> { }

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector2<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 2 }
    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector3<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 3 }
    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector4<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }
    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector2<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 2 }
    fn from_value(value: Self::Scalar) -> Self {
        na::Vector2::repeat(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector3<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 3 }
    fn from_value(value: Self::Scalar) -> Self {
        na::Vector3::repeat(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector4<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }
    fn from_value(value: Self::Scalar) -> Self {
        na::Vector4::repeat(value)
    }
}

/// Describes objects that can be located by r-trees.
///
/// See the `primitives` module for some basic implementations which can also serve 
/// as useful examples for own implementations.
pub trait SpatialObject {
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
/// use spade::{RTree, HasPosition};
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
    type Vector: VectorN;
    fn position(&self) -> Self::Vector;
}

/// Like `HasPosition`, but will only work for two dimensional vectors.
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
        (self.position() - point.clone()).length2()
    }

    fn contains(&self, point: &S::Vector) -> bool {
        self.position() == *point
    }
}

