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

use cgmath::{BaseNum, Vector2, zero};
use cgmath::conv::array2;
use num::{Bounded, Signed};
use boundingvolume::BoundingRect;
use misc::{length2};

pub trait RTreeNum: BaseNum + Bounded + Signed { }

impl <S> RTreeNum for S where S: BaseNum + Bounded + Signed {}

/// Describes objects that can be located by r-trees.
///
/// See the `primitives` module for some basic implementations.
/// # Example
/// ```
/// use rtree::{SpatialObject, BoundingRect};
/// 
/// struct Circle {
///  radius: f32,
///  origin: [f32; 2],
/// }
/// 
/// impl SpatialObject for Circle {
///   type Scalar = f32;
///
///   fn mbr(&self) -> BoundingRect<f32> {
///     let (x, y) = (self.origin[0], self.origin[1]);
///     let r = self.radius;
///     BoundingRect::from_corners(&[x - r, y - r], &[x + r, y + r])
///   }
///
///   fn distance(&self, point: [f32; 2]) -> f32 {
///     let dx = self.origin[0] - point[0];
///     let dy = self.origin[1] - point[1];
///     ((dx * dx + dy * dy).sqrt() - self.radius).max(0.)
///   }
/// }
/// ```
pub trait SpatialObject {
    type Scalar: RTreeNum;

    /// Returns the object's minimal bounding rectangle.
    ///
    /// The minimal bounding rectangle is the smallest axis aligned rectangle that completely
    // contains the object.
    /// <b>Note:</b> The rectangle must be as small as possible, otherwise some queries
    /// might fail.
    fn mbr(&self) -> BoundingRect<Self::Scalar>;

    /// Returns the distance from the object's contour.
    /// Note that this is not necessarily the euclidean distance,
    /// this functions result's will only be used for comparison.
    /// Returns zero if the point is contained within the object.
    fn distance(&self, point: [Self::Scalar; 2]) -> Self::Scalar;

    /// Returns true if a given point is contained in this object.
    fn contains(&self, point: [Self::Scalar; 2]) -> bool {
        self.distance(point) == zero()
    }
}

impl <S: RTreeNum> SpatialObject for Vector2<S> {
    type Scalar = S;
    fn mbr(&self) -> BoundingRect<S> {
        BoundingRect::from_point(array2(*self))
    }

    fn distance(&self, point: [Self::Scalar; 2]) -> S {
        let point: Vector2<_> = point.into();
        length2(&(point - *self))
    }
}
