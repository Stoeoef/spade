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


use nalgebra as na;
use cgmath as cg;
use nalgebra::{Repeat};

use std::fmt::Debug;
use traits::SpadeNum;
use num::{zero};
use misc::{min_inline, max_inline};

/// Abstraction over vectors in different dimensions.
pub trait VectorN where Self: Clone,
Self: Debug,
Self: PartialEq {
    /// The vector's internal scalar type.
    type Scalar: SpadeNum;

    /// Creates a new vector with all compoenents set to a certain value.
    fn from_value(value: Self::Scalar) -> Self;

    /// Returns the nth element of the vector.
    fn nth(&self, index: usize) -> &Self::Scalar;
    /// Returns a mutable reference to the nth element of the vector.
    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar;

    /// Creates a new vector with all components initialized to zero.
    fn new() -> Self {
        Self::from_value(zero())
    }

    /// The (fixed) number of dimensions of this vector trait.
    fn dimensions() -> usize;

    fn add(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| l + r)
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| l - r)
    }

    fn div(&self, scalar: Self::Scalar) -> Self {
        self.map(|x| x / scalar.clone())
    }

    fn mul(&self, scalar: Self::Scalar) -> Self {
        self.map(|x| x * scalar.clone())
    }

    /// Applies a binary operation component wise.
    fn component_wise<F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar>(&self, rhs: &Self, f: F) -> Self {
        let mut result = self.clone();
        for i in 0 .. Self::dimensions() {
            *result.nth_mut(i) = f(self.nth(i).clone(), rhs.nth(i).clone());
        }
        result
    }

    /// Maps an unary operation to all compoenents.
    fn map<F: Fn(Self::Scalar) -> O::Scalar, O: VectorN>(&self, f: F) -> O {
        let mut result = O::new();
        for i in 0 .. Self::dimensions() {
            *result.nth_mut(i)  = f(self.nth(i).clone());
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
            acc = f(acc, self.nth(i).clone());
        }
        acc
    }

    /// Checks if a property holds for all components.
    fn all_comp_wise<F: Fn(Self::Scalar, Self::Scalar) -> bool>(&self, rhs: &Self, f: F) -> bool {
        for i in 0 .. Self::dimensions() {
            if !f(self.nth(i).clone(), rhs.nth(i).clone()) {
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

/// A two dimensional Vector.
/// Some datastructures will only work if two dimensional vectors are given,
/// this trait makes sure that only such vectors can be passed.
pub trait TwoDimensional : VectorN { }

impl <S: SpadeNum + cg::BaseNum> TwoDimensional for cg::Vector2<S> { }
impl <S: SpadeNum + na::BaseNum> TwoDimensional for na::Vector2<S> { }
impl <S: SpadeNum + Copy> TwoDimensional for [S; 2] { }

/// A three dimensional Vector.
/// Some algorithms will only work with three dimensional vectors, this trait makes
/// sure that only such vectors can be used.
pub trait ThreeDimensional : VectorN {
    /// The cross product of this vector and another.
    fn cross(&self, other: &Self) -> Self;
}

impl <S: SpadeNum + cg::BaseNum> ThreeDimensional for cg::Vector3<S> {
    fn cross(&self, other: &Self) -> Self {
        cg::Vector3::cross(*self, *other)
    }
}

impl <S: SpadeNum + na::BaseNum> ThreeDimensional for na::Vector3<S> {
    fn cross(&self, other: &Self) -> Self {
        na::cross(self, other)
    }
}

impl <S: SpadeNum + Copy> ThreeDimensional for [S; 3] {
    fn cross(&self, other: &Self) -> Self {
        [
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0]
        ]
    }
}

impl <S: SpadeNum + Copy> VectorN for [S; 2] {
    type Scalar = S;
    fn dimensions() -> usize { 2 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }
    
    fn from_value(value: Self::Scalar) -> Self {
        [value; 2]
    }
}

impl <S: SpadeNum + Copy> VectorN for [S; 3] {
    type Scalar = S;
    fn dimensions() -> usize { 3 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }
    
    fn from_value(value: Self::Scalar) -> Self {
        [value; 3]
    }
}

impl <S: SpadeNum + Copy> VectorN for [S; 4] {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }
    
    fn from_value(value: Self::Scalar) -> Self {
        [value; 4]
    }
}

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector2<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 2 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector3<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 3 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + cg::BaseNum> VectorN for cg::Vector4<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector2<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 2 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        na::Vector2::repeat(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector3<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 3 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        na::Vector3::repeat(value)
    }
}

impl<S: SpadeNum + na::BaseNum> VectorN for na::Vector4<S> {
    type Scalar = S;
    
    fn dimensions() -> usize { 4 }

    fn nth(&self, index: usize) -> &S { &self[index] }
    fn nth_mut(&mut self, index: usize) -> &mut S { &mut self[index] }

    fn from_value(value: Self::Scalar) -> Self {
        na::Vector4::repeat(value)
    }
}
