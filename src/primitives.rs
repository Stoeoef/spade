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

//! Contains some useful primitives that can be inserted into r-trees.
//!
//! Use these objects if only the geometrical properties (position and size)
//! are important. If additional data needs to be stored per object, consider
//! implementing `SpatialObject`.


use cgmath::{Vector3, Zero, One};
use traits::{RTreeFloat, RTreeNum, SpatialObject, VectorN};
use num::{Float, one, zero, Signed};
use boundingvolume::BoundingRect;

#[derive(Clone, Debug)]
/// An edge defined by it's two end points.
pub struct SimpleEdge<V: VectorN> {
    /// The edge's origin.
    pub from: V,
    /// The edge's destination.
    pub to: V,
}

pub struct EdgeSideInfo<S> {
    signed_side: S,
}

impl <S> EdgeSideInfo<S> where S: RTreeNum  {
    /// Returns `true` if the query point lies on the left side of the query edge.
    pub fn is_on_left_side(&self) -> bool {
        self.signed_side > S::zero()
    }

    /// Returns `true` if the query point lies on the right side of the query edge.
    pub fn is_on_right_side(&self) -> bool {
        self.signed_side < S::zero()
    }

    /// Returns `true` if the query point lies on the left side or on the query edge.
    pub fn is_on_left_side_or_on_line(&self) -> bool {
        self.signed_side >= S::zero()
    }

    /// Returns `true` if the query point lies on the right side or on the query edge.
    pub fn is_on_right_side_or_on_line(&self) -> bool {
        self.signed_side <= S::zero()
    }

    /// Returns `true` if the query point lies on the line.
    pub fn is_on_line(&self) -> bool {
        self.signed_side.abs() == zero()
    }
}

impl <V> SimpleEdge<V> where V: VectorN {
    /// Creates a new edge from `from` to `to`.
    pub fn new(from: V, to: V) -> SimpleEdge<V> {
        SimpleEdge {
            from: from,
            to: to,
        }
    }

    /// Determines on which side of this edge a given point lies.
    pub fn side_query(&self, q: &V) -> EdgeSideInfo<V::Scalar> {
        let (a, b) = (&self.from, &self.to);
        let signed_side = (b[0].clone() - a[0].clone()) * (q[1].clone() - a[1].clone()) 
            - (b[1].clone() - a[1].clone()) * (q[0].clone() - a[0].clone());
        EdgeSideInfo { signed_side: signed_side }
    }

    pub fn is_point_on_edge(&self, query_point: &V) -> bool {
        self.is_projection_on_edge(query_point)
            && self.side_query(query_point).is_on_line()
    }

    pub fn is_projection_on_edge(&self, query_point: &V) -> bool {
        let (p1, p2) = (self.from.clone(), self.to.clone());
        let dir = p2 - p1.clone();
        let s = (query_point.clone() - p1).dot(&dir);
        zero::<V::Scalar>() <= s && s <= dir.length2()

    }
}

impl <V> SimpleEdge<V> where V: VectorN, V::Scalar: RTreeFloat {
    /// Yields the nearest point on this edge.
    pub fn nearest_point(&self, query_point: &V) -> V {
        let (p1, p2) = (self.from.clone(), self.to.clone());
        let dir = p2.clone() - p1.clone();
        let s = self.project_point(query_point);
        if V::Scalar::zero() < s && s < one() {
            p1 + dir * s

        } else {
            if s <= V::Scalar::zero() {
                p1
            } else {
                p2
            }
        }
    }

    pub fn projection_distance2(&self, query_point: &V) -> V::Scalar {
        let s = self.project_point(query_point);
        let p = self.from.clone() + (self.to.clone() - self.from.clone()) * s;
        p.distance2(query_point)
    }

    /// Projects a point on this line and returns its relative position.
    /// 
    /// This method will return a value between 0. and 1. (linearly interpolated) if the projected
    /// point lies between `self.from` and `self.to`, a value close to zero (due to rounding errors)
    /// if the projected point is equal to `self.from` and a value smaller than zero if the projected
    /// point lies "before" `self.from`. Analogously, a value close to 1. or greater than 1. is
    /// returned if the projected point is equal to or lies behind `self.to`.
    pub fn project_point(&self, query_point: &V) -> V::Scalar {
        let (p1, p2) = (self.from.clone(), self.to.clone());
        let dir = p2 - p1.clone();
        (query_point.clone() - p1).dot(&dir) / dir.length2()
    }
}

impl <V: VectorN> SpatialObject for SimpleEdge<V> where V::Scalar: RTreeFloat {
    type Vector = V;

    fn mbr(&self) -> BoundingRect<V> {
        BoundingRect::from_corners(&self.from, &self.to)
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        let nn = self.nearest_point(point);
        (point.clone() - nn).length2()
    }
}

/// A triangle, defined by it's three points.
#[derive(Clone)]
pub struct SimpleTriangle<V: VectorN> {
    v0: V,
    v1: V,
    v2: V,
}

impl <V> SimpleTriangle<V> where V: VectorN {
    /// Checks if the given points are ordered counter clock wise.
    pub fn is_ordered_ccw(v0: &V, v1: &V, v2: &V) -> bool {
        let edge = SimpleEdge::new(v0.clone(), v1.clone());
        edge.side_query(v2).is_on_left_side_or_on_line()
    }

    pub fn new(v0: V, v1: V, v2: V) -> SimpleTriangle<V> {
        SimpleTriangle { v0: v0, v1: v1, v2: v2 }
    }

    pub fn vertices(&self) -> [&V; 3] {
        [&self.v0, &self.v1, &self.v2]
    }

    pub fn double_area(&self) -> V::Scalar {
        let b = self.v1.clone() - self.v0.clone();
        let c = self.v2.clone() - self.v0.clone();
        (b[0].clone() * c[1].clone() - b[1].clone() * c[0].clone()).abs()
    }

}

impl <V> PartialEq for SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    fn eq(&self, rhs: &SimpleTriangle<V>) -> bool {
        let vl = self.vertices();
        let vr = rhs.vertices();
        if let Some(index) = vr.iter().position(|v| *v == vl[0]) {
            let r1 = vr[(index + 1) % 3];
            let r2 = vr[(index + 2) % 3];
            vl[1] == r1 && vl[2] == r2
        } else {
            false
        }
    }
}

impl <V> SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    pub fn circumcenter(&self) -> V {
        let one: V::Scalar = One::one();
        let two = one + one;
        let b = self.v1.clone() - self.v0.clone();
        let c = self.v2.clone() - self.v0.clone();
        // Calculate circumcenter position
        let d = two * (b[0] * c[1] - c[0] * b[1]);
        let len_b = b.dot(&b);
        let len_c = c.dot(&c);
        let x = (len_b * c[1] - len_c * b[1]) / d;
        let y = (-len_b * c[0] + len_c * b[0]) / d;
        let mut result = V::new();
        result[0] = x;
        result[1] = y;
        result + self.v0.clone()
    }

    pub fn barycentric_interpolation(&self, coord: &V) -> Vector3<V::Scalar> {
        let (v1, v2, v3) = (self.v0.clone(), self.v1.clone(), self.v2.clone());
        let (x, y) = (coord[0], coord[1]);
        let (x1, x2, x3) = (v1[0], v2[0], v3[0]);
        let (y1, y2, y3) = (v1[1], v2[1], v3[1]);
        let det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        let lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
        let lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
        let lambda3 = one::<V::Scalar>() - lambda1 - lambda2;
        Vector3::new(lambda1, lambda2, lambda3)
    }
}


impl <V> SpatialObject for SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    type Vector = V;
    
    fn mbr(&self) -> BoundingRect<V> {
        let mut result = BoundingRect::from_corners(&self.v0, &self.v1);
        result.add_point(self.v2.clone());
        result
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        let ordered_ccw = SimpleTriangle::is_ordered_ccw(&self.v0, &self.v1, &self.v2);
        for i in 0 .. 3 {
            let edge = SimpleEdge::new(self.vertices()[i].clone(), 
                                       self.vertices()[(i + 1) % 3].clone());
            if edge.side_query(point).is_on_right_side() == ordered_ccw {
                return edge.distance2(point);
            }
        }
        // The point is within the triangle
        zero()
    }
}

pub struct SimpleCircle<V: VectorN> {
    pub origin: V,
    pub radius: V::Scalar,
}

impl <V> SimpleCircle<V> where V: VectorN, V::Scalar: RTreeFloat {
    pub fn new(origin: V, radius: V::Scalar) -> SimpleCircle<V> {
        SimpleCircle {
            origin: origin,
            radius: radius,
        }
    }
}

impl <V> SpatialObject for SimpleCircle<V> where V: VectorN, V::Scalar: RTreeFloat {
    type Vector = V;

    fn mbr(&self) -> BoundingRect<V> {
        let r = V::from_value(self.radius);
        BoundingRect::from_corners(&(self.origin.clone() - r.clone()), &(self.origin.clone() + r.clone()))
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        let dx = self.origin[0] - point[0];
        let dy = self.origin[1] - point[1];
        let dist = ((dx * dx + dy * dy).sqrt() - self.radius).max(zero());
        dist * dist
    }

    // Since containment checks do not require the calculation of the square root
    // we can redefine this method.
    fn contains(&self, point: &V) -> bool {
        let dx = self.origin[0] - point[0];
        let dy = self.origin[1] - point[1];
        let r2 = self.radius * self.radius;
        dx * dx + dy * dy <= r2
    }
}

#[cfg(test)]
mod test {
    use super::{SimpleEdge, SimpleTriangle};
    use traits::SpatialObject;
    use cgmath::{Vector2, ApproxEq};
    
    #[test]
    fn test_edge_distance() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.distance2(&Vector2::new(1.0, 0.0)).approx_eq(&0.5));
        assert!(e.distance2(&Vector2::new(0.0, 1.)).approx_eq(&0.5));
        assert!(e.distance2(&Vector2::new(-1.0, -1.)).approx_eq(&2.));
        assert!(e.distance2(&Vector2::new(2.0, 2.0)).approx_eq(&2.));
    }

    #[test]
    fn test_edge_side() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.side_query(&Vector2::new(1.0, 0.0)).is_on_right_side());
        assert!(e.side_query(&Vector2::new(0.0, 1.0)).is_on_left_side());
        assert!(e.side_query(&Vector2::new(0.5, 0.5)).is_on_line());
    }

    #[test]
    fn test_triangle_distance() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1., 0.);
        let v3 = Vector2::new(0., 1.);
        let t = SimpleTriangle::new(v1, v2, v3);
        assert_eq!(t.distance2(&Vector2::new(0.25, 0.25)), 0.);
        assert!(t.distance2(&Vector2::new(-1., -1.)).approx_eq(&2.));
        assert!(t.distance2(&Vector2::new(0., -1.)).approx_eq(&1.));
        assert!(t.distance2(&Vector2::new(-1., 0.)).approx_eq(&1.));
        assert!(t.distance2(&Vector2::new(1., 1.)).approx_eq(&0.5));
        assert!(t.distance2(&Vector2::new(0.5, 0.5)).approx_eq(&0.0));
        assert!(t.distance2(&Vector2::new(0.6, 0.6)) > 0.001);
    }

    #[test]
    fn test_triangle_is_ordered_ccw() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1f32, 0.);
        let v3 = Vector2::new(0f32, 1.);
        assert!(SimpleTriangle::is_ordered_ccw(&v1, &v2, &v3));
        assert!(!SimpleTriangle::is_ordered_ccw(&v2, &v1, &v3));
        assert!(SimpleTriangle::is_ordered_ccw(&v3, &v1, &v2));
    }        

}
