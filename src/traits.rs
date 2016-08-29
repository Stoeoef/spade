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

use cgmath::{BaseNum, BaseFloat, Vector2, Vector3, Vector4, Array, Zero};
use num::{Bounded, Signed};
use misc::{min_inline, max_inline};
use boundingvolume::BoundingRect;
use std::ops::{Add, Sub, Index, IndexMut, Div, Mul};

pub trait RTreeNum: BaseNum + Bounded + Signed { }
pub trait RTreeFloat: RTreeNum + BaseFloat { }

impl <S> RTreeNum for S where S: BaseNum + Bounded + Signed { }
impl <S> RTreeFloat for S where S: RTreeNum + BaseFloat { }

/// Abstraction over Vectors in different dimensions.
/// Implement this if you need to use your own vector classes.
pub trait VectorN where Self: Copy, Self: Clone,
Self: Add<Self, Output=Self>,
Self: Sub<Self, Output=Self>,
Self: Div<<Self as VectorN>::Scalar, Output=Self>,
Self: Mul<<Self as VectorN>::Scalar, Output=Self>,
Self: Index<usize, Output=<Self as VectorN>::Scalar>,
Self: IndexMut<usize, Output=<Self as VectorN>::Scalar>,
Self: ::std::fmt::Debug,
Self: PartialEq {
    type Scalar: RTreeNum;

    /// Creates a new vector with all compoenents set to a certain value.
    fn from_value(value: Self::Scalar) -> Self;

    /// Creates a new vector with all components initialized to zero.
    fn new() -> Self {
        Self::from_value(Self::Scalar::zero())
    }

    /// The (fixed) number of dimensions of this vector trait.
    fn dimensions() -> usize;

    /// Applies a binary operation component wise.
    fn component_wise<F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar>(&self, rhs: &Self, f: F) -> Self {
        let mut result = self.clone();
        for i in 0 .. Self::dimensions() {
            result[i] = f(self[i], rhs[i]);
        }
        result
    }

    /// Maps an unary operation to all compoenents.
    fn map<F: Fn(Self::Scalar) -> Self::Scalar>(&self, f: F) -> Self {
        let mut result = self.clone();
        for i in 0 .. Self::dimensions() {
            result[i]  = f(self[i]);
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
            acc = f(acc, self[i]);
        }
        acc
    }

    /// Checks if a property holds for all components.
    fn all_comp_wise<F: Fn(Self::Scalar, Self::Scalar) -> bool>(&self, rhs: &Self, f: F) -> bool {
        for i in 0 .. Self::dimensions() {
            if !f(self[i], rhs[i]) {
                return false;
            }
        }
        true
    }

    /// Returns the vector's dot product.
    fn dot(&self, rhs: &Self) -> Self::Scalar {
        self.component_wise(rhs, |l, r| l * r).fold(Self::Scalar::zero(), 
                                                    |acc, val| acc + val)
    }

    /// Returns the vector's squared length.
    fn length2(&self) -> Self::Scalar {
        self.dot(&self)
    }

}

impl<S: RTreeNum> VectorN for Vector2<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 2 }
    fn from_value(value: Self::Scalar) -> Self {
        Array::from_value(value)
    }
}

impl<S: RTreeNum> VectorN for Vector3<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 3 }
    fn from_value(value: Self::Scalar) -> Self {
        Array::from_value(value)
    }
}

impl<S: RTreeNum> VectorN for Vector4<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }
    fn from_value(value: Self::Scalar) -> Self {
        Array::from_value(value)
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
    // contains the object.
    /// <b>Note:</b> The rectangle must be as small as possible, otherwise some queries
    /// might fail.
    fn mbr(&self) -> BoundingRect<Self::Vector>;

    /// Returns the distance from the object's contour.
    /// Note that this is not necessarily the euclidean distance,
    /// this functions result's will only be used for comparison.
    /// Returns zero if the point is contained within the object.
    fn distance(&self, point: Self::Vector) -> <Self::Vector as VectorN>::Scalar;

    /// Returns true if a given point is contained in this object.
    fn contains(&self, point: Self::Vector) -> bool {
        self.distance(point) <= <Self::Vector as VectorN>::Scalar::zero()
    }
}

pub trait HasPosition {
    type Vector: VectorN;
    fn position(&self) -> Self::Vector;
}

impl <V> HasPosition for V where V: VectorN {
    type Vector = V;
    fn position(&self) -> V {
        *self
    }
}

impl <S>  SpatialObject for S where S: HasPosition {
    type Vector = S::Vector;

    fn mbr(&self) -> BoundingRect<S::Vector> {
        BoundingRect::from_point(self.position())
    }

    fn distance(&self, point: S::Vector) -> <S::Vector as VectorN>::Scalar {
        (self.position() - point).length2()
    }
}
