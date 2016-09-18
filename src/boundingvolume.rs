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

use num::{Signed, zero, one};
use traits::{VectorN};
use misc::max_inline;

/// An axis aligned minimal bounding rectangle.
///
/// A axis aligned minimal bounding rectangle is the smallest rectangle that completly surrounds an object and is aligned along the two axis in the plane.
#[derive(Clone, PartialEq, Debug)]
pub struct BoundingRect<V: VectorN> {
    lower: V,
    upper: V,
}

impl <V> BoundingRect<V> where V: VectorN {

    /// Creates a bounding rectangle that contains exactly one point.
    ///
    /// This will create a bounding rectangle with `lower == upper == point`.
    pub fn from_point(point: V) -> BoundingRect<V> {
        BoundingRect {
            lower: point.clone(),
            upper: point,
        }
    }

    /// Creates a bounding rectangle that contains two points.
    pub fn from_corners(corner1: &V, corner2: &V) -> BoundingRect<V> {
        BoundingRect {
            lower: corner1.min_vec(&corner2),
            upper: corner1.max_vec(&corner2),
        }
    }

    /// Returns the lower corner of the bounding rectangle.
    ///
    /// The lower corner has the smaller coordinates.
    pub fn lower(&self) -> V {
        self.lower.clone()
    }

    /// Returns the upper corner of the bounding rectangle.
    ///
    /// The upper corner has the larger coordinates.
    pub fn upper(&self) -> V {
        self.upper.clone()
    }

    /// Checks if a point is contained within the bounding rectangle.
    ///
    /// A point lying exactly on the bounding rectangle's border is also contained.
    #[inline]
    pub fn contains_point(&self, point: &V) -> bool {
        self.lower.all_comp_wise(&point, |l, r| l <= r) &&
            self.upper.all_comp_wise(&point, |l, r| l >= r)
    }
    
    /// Checks if another bounding rectangle is completley contained withing this rectangle.
    ///
    /// A rectangle is contained if and only if all four corner are contained (see `contains_point`).
    #[inline]
    pub fn contains_rect(&self, rect: &BoundingRect<V>) -> bool {
        self.lower.all_comp_wise(&rect.lower(), |l, r| l <= r) &&
            self.upper.all_comp_wise(&rect.upper(), |l, r| l >= r)
    }

    /// Enlarges this bounding rectangle to contain a point.
    ///
    /// If the point is already contained, nothing will be changed.
    /// Otherwise, this will enlarge `self` to be just large enough
    /// to contain the new point.
    #[inline]
    pub fn add_point(&mut self, point: V) {
        self.lower = self.lower.min_vec(&point);
        self.upper = self.upper.max_vec(&point);
    }

    /// Enlarges this bounding rectangle to contain a rectangle.
    ///
    /// If the rectangle is already contained, nothing will be changed.
    /// Otherwise, this will enlarge `self` to be just large enough
    // to contain the new rectangle.
    #[inline]
    pub fn add_rect(&mut self, rect: &BoundingRect<V>) {
        self.lower = self.lower.min_vec(&rect.lower);
        self.upper = self.upper.max_vec(&rect.upper);
    }

    /// Returns the rectangle's area.
    pub fn area(&self) -> V::Scalar {
        let diag = self.upper() - self.lower();
        diag.fold(one(), |acc, value| max_inline(acc * value, zero()))
    }

    /// Returns half of the rectangle's margin, thus `width + height`.
    pub fn half_margin(&self) -> V::Scalar {
        let diag = self.upper() - self.lower();
        diag.fold(zero(), |acc, value| max_inline(acc + value, zero()))
    }

    /// Returns the rectangle's center.
    pub fn center(&self) -> V {
        let two = one::<V::Scalar>() + one::<V::Scalar>();
        self.lower() + (self.upper() - self.lower()) / two
    }

    /// Returns the intersection of this and another bounding rectangle.
    ///
    /// If the rectangles do not intersect, a bounding rectangle with an area and
    /// margin of zero is returned.
    pub fn intersect(&self, other: &BoundingRect<V>) -> BoundingRect<V> {
        BoundingRect {
            lower: self.lower.max_vec(&other.lower),
            upper: self.upper.min_vec(&other.upper),
        }
    }

    /// Returns true if this and another bounding rectangle intersect each other.
    /// If the rectangles just "touch" each other at one side, true is returned.
    pub fn intersects(&self, other: &BoundingRect<V>) -> bool {
        self.lower.all_comp_wise(&other.upper(), |l, r| l <= r) &&
            self.upper.all_comp_wise(&other.lower(), |l, r| l >= r)
    }

    #[doc(hidden)]
    pub fn min_point(&self, point: &V) -> V {
        self.upper.min_vec(&self.lower.max_vec(&point))
    }

    #[doc(hidden)]
    pub fn min_dist2(&self, point: &V) -> V::Scalar {
        (self.min_point(point) - point.clone()).length2()
    }

    #[doc(hidden)]
    pub fn max_dist2(&self, point: &V) -> V::Scalar {
        let l = self.lower();
        let u = self.upper();
        let d1 = (l - point.clone()).map(|v| v.abs());
        let d2 = (u - point.clone()).map(|v| v.abs());
        let max_delta = d1.max_vec(&d2);
        max_delta.length2()
    }

    #[doc(hidden)]
    pub fn min_max_dist2(&self, point: &V) -> V::Scalar {
        let l = self.lower();
        let u = self.upper();
        let (mut min, mut max) = (V::new(), V::new());
        for i in 0 .. V::dimensions() {
            if (l[i].clone() - point[i].clone()).abs() < (u[i].clone() - point[i].clone()).abs() { 
                min[i] = l[i].clone();
                max[i] = u[i].clone();
            } else {
                min[i] = u[i].clone();
                max[i] = l[i].clone();
            }
        }
        let mut result = zero();
        for i in 0 .. V::dimensions() {
            let mut p = min.clone();
            // Only set one component to the maximum distance
            p[i] = max[i].clone();
            let new_dist = (p - point.clone()).length2();
            if new_dist < result || i == 0 {
                result = new_dist
            }
        }
        result
    }
}
