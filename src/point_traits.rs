// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use cgmath as cg;
use nalgebra as na;

use crate::misc::{max_inline, min_inline};
use crate::traits::SpadeNum;
use num::zero;
use std::fmt::Debug;

/// Abstraction over a point with a fixed number of dimensions.
///
/// Spade will work with any point type implementing this trait, at the
/// moment points of the `cgmath` and `nalgebra` crates are supported.
/// Also, the trait is implemented for fixed arrays of length 2, 3 and 4, allowing
/// to use spade's datastructures with fixed size arrays as point coordinates.
/// That means that the trait's methods are also implemented for
/// these array types, thus be careful when importing `PointN`.
///
/// Implement this if you want spade to support your own point types.
/// Also, an empty implementation of `TwoDimensional` or `ThreeDimensional`
/// should be given if appropriate.
pub trait PointN
where
    Self: Clone + Debug + PartialEq,
{
    /// The points's internal scalar type.
    type Scalar: SpadeNum;

    /// The (fixed) number of dimensions of this point type.
    fn dimensions() -> usize;

    /// Creates a new point with all components set to a certain value.
    fn from_value(value: Self::Scalar) -> Self;

    /// Returns the nth element of this point.
    fn nth(&self, index: usize) -> &Self::Scalar;
    /// Returns a mutable reference to the nth element of this point.
    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar;
}

/// Adds some private methods to the `PointN` trait.
pub trait PointNExtensions: PointN {
    /// Creates a new point with all components initialized to zero.
    fn new() -> Self {
        Self::from_value(zero())
    }

    /// Adds two points.
    fn add(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| l + r)
    }

    /// Substracts two points.
    fn sub(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, |l, r| l - r)
    }

    /// Divides this point with a scalar value.
    fn div(&self, scalar: Self::Scalar) -> Self {
        self.map(|x| x / scalar.clone())
    }

    /// Multiplies this point with a scalar value.
    fn mul(&self, scalar: Self::Scalar) -> Self {
        self.map(|x| x * scalar.clone())
    }

    /// Applies a binary operation component wise.
    fn component_wise<F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar>(
        &self,
        rhs: &Self,
        f: F,
    ) -> Self {
        let mut result = self.clone();
        for i in 0..Self::dimensions() {
            *result.nth_mut(i) = f(self.nth(i).clone(), rhs.nth(i).clone());
        }
        result
    }

    /// Maps an unary operation to all components.
    fn map<F: Fn(Self::Scalar) -> O::Scalar, O: PointN>(&self, f: F) -> O {
        let mut result = O::new();
        for i in 0..Self::dimensions() {
            *result.nth_mut(i) = f(self.nth(i).clone());
        }
        result
    }

    /// Returns a new point containing the minimum values of this and another point (componentwise)
    fn min_point(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, min_inline)
    }

    /// Returns a new point containing the maximum values of this and another point (componentwise)
    fn max_point(&self, rhs: &Self) -> Self {
        self.component_wise(rhs, max_inline)
    }

    /// Fold operation over all point components.
    fn fold<T, F: Fn(T, Self::Scalar) -> T>(&self, mut acc: T, f: F) -> T {
        for i in 0..Self::dimensions() {
            acc = f(acc, self.nth(i).clone());
        }
        acc
    }

    /// Checks if a property holds for all components of this and another point.
    fn all_comp_wise<F: Fn(Self::Scalar, Self::Scalar) -> bool>(&self, rhs: &Self, f: F) -> bool {
        for i in 0..Self::dimensions() {
            if !f(self.nth(i).clone(), rhs.nth(i).clone()) {
                return false;
            }
        }
        true
    }

    /// Returns the point's dot product.
    fn dot(&self, rhs: &Self) -> Self::Scalar {
        self.component_wise(rhs, |l, r| l * r)
            .fold(zero(), |acc, val| acc + val)
    }

    /// Returns the point's squared length.
    fn length2(&self) -> Self::Scalar {
        self.dot(&self)
    }

    fn lex_compare(&self, other: &Self) -> std::cmp::Ordering {
        for i in 0..Self::dimensions() {
            let left = self.nth(i);
            let right = other.nth(i);
            let compare = left.partial_cmp(right).unwrap();
            if compare != std::cmp::Ordering::Equal {
                return compare;
            }
        }
        std::cmp::Ordering::Equal
    }
}

impl<T> PointNExtensions for T where T: PointN {}

/// A two dimensional Point.
///
/// Some datastructures will only work if two dimensional points are given,
/// this trait makes sure that only such points can be passed.
pub trait TwoDimensional: PointN {}

impl<S: SpadeNum + cg::BaseNum> TwoDimensional for cg::Point2<S> {}
impl<S: SpadeNum + na::Scalar> TwoDimensional for na::Point2<S> {}
impl<S: SpadeNum + Copy> TwoDimensional for [S; 2] {}

/// A three dimensional Point.
///
/// Some algorithms will only work with three dimensional points, this trait makes
/// sure that only such points can be used.
pub trait ThreeDimensional: PointN {
    /// The cross product of this point and another.
    fn cross(&self, other: &Self) -> Self {
        let mut result = Self::new();
        *result.nth_mut(0) =
            self.nth(1).clone() * other.nth(2).clone() - self.nth(2).clone() * other.nth(1).clone();
        *result.nth_mut(1) =
            self.nth(2).clone() * other.nth(0).clone() - self.nth(0).clone() * other.nth(2).clone();
        *result.nth_mut(2) =
            self.nth(0).clone() * other.nth(1).clone() - self.nth(1).clone() * other.nth(0).clone();
        result
    }
}

impl<S: SpadeNum + cg::BaseNum> ThreeDimensional for cg::Point3<S> {}

impl<S: SpadeNum + na::Scalar> ThreeDimensional for na::Point3<S> {}

impl<S: SpadeNum + Copy> ThreeDimensional for [S; 3] {}

impl<S: SpadeNum + Copy> PointN for [S; 2] {
    type Scalar = S;
    fn dimensions() -> usize {
        2
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        [value; 2]
    }
}

impl<S: SpadeNum + Copy> PointN for [S; 3] {
    type Scalar = S;
    fn dimensions() -> usize {
        3
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        [value; 3]
    }
}

impl<S: SpadeNum + Copy> PointN for [S; 4] {
    type Scalar = S;

    fn dimensions() -> usize {
        4
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        [value; 4]
    }
}

impl<S: SpadeNum + cg::BaseNum> PointN for cg::Point2<S> {
    type Scalar = S;

    fn dimensions() -> usize {
        2
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + cg::BaseNum> PointN for cg::Point3<S> {
    type Scalar = S;

    fn dimensions() -> usize {
        3
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        cg::Array::from_value(value)
    }
}

impl<S: SpadeNum + na::Scalar> PointN for na::Point2<S> {
    type Scalar = S;

    fn dimensions() -> usize {
        2
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        na::Point2::new(value, value)
    }
}

impl<S: SpadeNum + na::Scalar> PointN for na::Point3<S> {
    type Scalar = S;

    fn dimensions() -> usize {
        3
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        na::Point3::new(value, value, value)
    }
}

impl<S: SpadeNum + na::Scalar + na::Scalar> PointN for na::Point4<S> {
    type Scalar = S;

    fn dimensions() -> usize {
        4
    }

    fn nth(&self, index: usize) -> &S {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut S {
        &mut self[index]
    }

    fn from_value(value: Self::Scalar) -> Self {
        na::Point4::new(value, value, value, value)
    }
}
