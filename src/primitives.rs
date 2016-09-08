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
use traits::{RTreeFloat, SpatialObject, VectorN};
use num::{Float, one, zero};
use boundingvolume::BoundingRect;

#[derive(Clone, Debug)]
/// An edge defined by it's two end points.
pub struct SimpleEdge<V: VectorN> {
    /// The edge's origin.
    pub from: V,
    /// The edge's destination.
    pub to: V,
}

impl <V: VectorN> SimpleEdge<V> where V::Scalar: RTreeFloat {
    /// Creates a new edge from `from` to `to`.
    pub fn new(from: V, to: V) -> SimpleEdge<V> {
        SimpleEdge {
            from: from,
            to: to,
        }
    }

    /// Yields the nearest point on this edge.
    pub fn nearest_point(&self, query_point: V) -> V {
        let (p1, p2) = (self.from, self.to);
        let dir = p2 - p1;
        let s: V::Scalar = (query_point - p1).dot(&dir) / dir.length2();
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
    /// Returns a value indicating on which side of this edge a given point lies.
    ///
    /// Returns a value smaller than zero if `query_point` is on the right side,
    /// or greater than zero if `query_point` is on the left side.
    pub fn signed_side(&self, query_point: V) -> V::Scalar {
        let (ref a, ref b) = (self.from, self.to);
        (b[0] - a[0]) * (query_point[1] - a[1]) - (b[1] - a[1])  * (query_point[0] - a[0])
    }

    /// Returns `true` if a given point is on the right side of this edge.
    pub fn is_on_right_side(&self, query_point: V) -> bool {
        self.signed_side(query_point) < V::Scalar::zero()
    }

    /// Returns `false` if a given point is on the left side of this edge.
    pub fn is_on_left_side(&self, query_point: V) -> bool {
        self.signed_side(query_point) > V::Scalar::zero()
    }

    pub fn approx_is_on_right_side(&self, query_point: V, eps: V::Scalar) -> bool {
        self.signed_side(query_point) < eps
    }

    pub fn approx_is_on_left_side(&self, query_point: V, eps: V::Scalar) -> bool {
        self.signed_side(query_point) > -eps
    }

}

impl <V: VectorN> SpatialObject for SimpleEdge<V> where V::Scalar: RTreeFloat {
    type Vector = V;

    fn mbr(&self) -> BoundingRect<V> {
        BoundingRect::from_corners(&self.from, &self.to)
    }

    fn distance(&self, point: V) -> V::Scalar {
        let nn = self.nearest_point(point);
        (point - nn).length2()
    }
}

/// A triangle, defined by it's three points.
#[derive(Clone)]
pub struct SimpleTriangle<V: VectorN> where V::Scalar: RTreeFloat {
    vertices: [V; 3],
}

impl <V> SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    /// Checks if the given points are ordered counter clock wise.
    pub fn is_ordered_ccw(vertices: &[V; 3]) -> bool {
        let edge = SimpleEdge::new(vertices[0].clone(), vertices[1].clone());
        return edge.signed_side(vertices[2]) >= zero()
    }

    pub fn new(mut vertices: [V; 3]) -> SimpleTriangle<V> {
        if !SimpleTriangle::is_ordered_ccw(&vertices) {
            vertices.swap(1, 2);
        }
        SimpleTriangle { vertices: vertices }
    }

    pub fn vertices(&self) -> &[V; 3] {
        &self.vertices
    }

    pub fn area(&self) -> V::Scalar {
        let one: V::Scalar = One::one();
        let two = one + one;
        let a = self.vertices[0];
        let b = self.vertices[1] - a;
        let c = self.vertices[2] - a;
        (b[0] * c[1] - b[1] * c[0]).abs() / two
    }

    pub fn circumcenter(&self) -> V {
        let one: V::Scalar = One::one();
        let two = one + one;
        let a = self.vertices[0];
        let b = self.vertices[1] - a;
        let c = self.vertices[2] - a;
        // Calculate circumcenter position
        let d = two * (b[0] * c[1] - c[0] * b[1]);
        let len_b = b.dot(&b);
        let len_c = c.dot(&c);
        let x = (len_b * c[1] - len_c * b[1]) / d;
        let y = (-len_b * c[0] + len_c * b[0]) / d;
        let mut result = V::new();
        result[0] = x;
        result[1] = y;
        result + a
    }

    pub fn barycentric_interpolation(&self, coord: &V) -> Vector3<V::Scalar> {
        let vertices = self.vertices();
        let (v1, v2, v3) = (vertices[0], vertices[1], vertices[2]);
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

impl <V> PartialEq for SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    fn eq(&self, rhs: &SimpleTriangle<V>) -> bool {
        let vl = self.vertices;
        let vr = rhs.vertices;
        if let Some(index) = vr.iter().position(|v| *v == vl[0]) {
            let r1 = vr[(index + 1) % 3];
            let r2 = vr[(index + 2) % 3];
            vl[1] == r1 && vl[2] == r2
        } else {
            false
        }
    }
}


impl <V> SpatialObject for SimpleTriangle<V> where V: VectorN, V::Scalar: RTreeFloat {
    type Vector = V;
    
    fn mbr(&self) -> BoundingRect<V> {
        let mut result = BoundingRect::new();
        for vertex in self.vertices.iter() {
            result.add_point(vertex);
        }
        result
    }

    fn distance(&self, point: V) -> V::Scalar {
        for i in 0 .. 3 {
            let edge = SimpleEdge::new(self.vertices[i], self.vertices[(i + 1) % 3]);
            if edge.signed_side(point) <= zero() {
                return edge.distance(point);
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
        BoundingRect::from_corners(&(self.origin - r), &(self.origin + r))
    }

    fn distance(&self, point: V) -> V::Scalar {
        let dx = self.origin[0] - point[0];
        let dy = self.origin[1] - point[1];
        ((dx * dx + dy * dy).sqrt() - self.radius).max(zero())
    }

    // Since containment checks do not require the calculation of the square root
    // we can redefine this method.
    fn contains(&self, point: V) -> bool {
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
        assert!(e.distance(Vector2::new(1.0, 0.0)).approx_eq(&0.5));
        assert!(e.distance(Vector2::new(0.0, 1.)).approx_eq(&0.5));
        assert!(e.distance(Vector2::new(-1.0, -1.)).approx_eq(&2.));
        assert!(e.distance(Vector2::new(2.0, 2.0)).approx_eq(&2.));
    }

    #[test]
    fn test_edge_side() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.signed_side(Vector2::new(1.0, 0.0)) < 0.0);
        assert!(e.signed_side(Vector2::new(0.0, 1.0)) > 0.0);
        assert!(e.signed_side(Vector2::new(0.5, 0.5)).approx_eq(&0.0));
    }

    #[test]
    fn test_triangle_distance() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1., 0.);
        let v3 = Vector2::new(0., 1.);
        let t = SimpleTriangle::new([v1, v2, v3]);
        assert_eq!(t.distance(Vector2::new(0.25, 0.25)), 0.);
        assert!(t.distance(Vector2::new(-1., -1.)).approx_eq(&2.));
        assert!(t.distance(Vector2::new(0., -1.)).approx_eq(&1.));
        assert!(t.distance(Vector2::new(-1., 0.)).approx_eq(&1.));
        assert!(t.distance(Vector2::new(1., 1.)).approx_eq(&0.5));
        assert!(t.distance(Vector2::new(0.5, 0.5)).approx_eq(&0.0));
        assert!(t.distance(Vector2::new(0.6, 0.6)) > 0.001);
    }

    #[test]
    fn test_triangle_is_ordered_ccw() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1f32, 0.);
        let v3 = Vector2::new(0f32, 1.);
        assert!(SimpleTriangle::is_ordered_ccw(&[v1, v2, v3]));
        assert!(!SimpleTriangle::is_ordered_ccw(&[v2, v1, v3]));
        assert!(SimpleTriangle::is_ordered_ccw(&[v3, v1, v2]));
    }        

}
