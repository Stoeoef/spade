// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Contains some useful primitives that can be inserted into r-trees.
//!
//! Use these objects if only the geometrical properties (position and size)
//! are important. If additional data needs to be stored per object, consider
//! implementing `SpatialObject`.


use cgmath::{Point3, Zero, One};
use traits::{SpadeFloat, SpadeNum, SpatialObject};
use point_traits::{PointN, PointNExtensions, TwoDimensional};
use num::{Float, one, zero, Signed};
use boundingvolume::BoundingRect;
use kernels::{TrivialKernel, DelaunayKernel};

/// An edge defined by it's two end points.
#[derive(Clone, Debug)]
pub struct SimpleEdge<V: PointN> {
    /// The edge's origin.
    pub from: V,
    /// The edge's destination.
    pub to: V,
}


/// Yields information about on which side of a line a point lies.
#[derive(Debug, Clone, Copy)]
pub struct EdgeSideInfo<S> {
    signed_side: S,
}

impl <S> PartialEq for EdgeSideInfo<S>
    where S: SpadeNum
{
    fn eq(&self, other: &EdgeSideInfo<S>) -> bool {
        if self.is_on_line() || other.is_on_line() {
            self.is_on_line() && other.is_on_line()
        } else {
            self.is_on_right_side() == other.is_on_right_side()
        }
    }
}

impl <S> EdgeSideInfo<S> where S: SpadeNum  {
    #[doc(hidden)]
    pub fn from_determinant(s: S) -> EdgeSideInfo<S> {
        EdgeSideInfo { signed_side: s }
    }

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

    /// Returns `true` if the query point lies on an infinite line going
    /// through the edge's start and end point.
    pub fn is_on_line(&self) -> bool {
        self.signed_side.abs() == zero()
    }

    /// Inverts this edge side information.
    /// If this information indicates the position of a point P
    /// relative to an edge A -> B, the inverted information will
    /// indicate the relative position of P relative to the edge
    /// B -> A.
    pub fn reversed(&self) -> EdgeSideInfo<S> {
        EdgeSideInfo { signed_side: -self.signed_side.clone() }
    }
}

impl <V> SimpleEdge<V> where V: PointN {
    /// Creates a new edge from `from` to `to`.
    pub fn new(from: V, to: V) -> SimpleEdge<V> {
        SimpleEdge {
            from: from,
            to: to,
        }
    }

    /// Projects a point onto the infinite line going through the
    /// edge's start and end point and returns `true` if the projected
    /// points lies between `from` and `to`.
    pub fn is_projection_on_edge(&self, query_point: &V) -> bool {
        let (p1, p2) = (&self.from, &self.to);
        let dir = p2.sub(p1);
        let s = query_point.sub(p1).dot(&dir);
        zero::<V::Scalar>() <= s && s <= dir.length2()
    }

    /// Returns the edge's squared length.
    pub fn length2(&self) -> V::Scalar {
        let diff = self.from.sub(&self.to);
        diff.dot(&diff)
    }
}

impl <V> SimpleEdge<V> where V: TwoDimensional {
    /// Determines on which side of this edge a given point lies.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    ///
    /// use nalgebra::Point2;
    /// use spade::kernels::TrivialKernel;
    /// use spade::primitives::SimpleEdge;
    ///
    /// # fn main() {
    /// let e = SimpleEdge::new(Point2::new(0f32, 0.), Point2::new(1., 1.));
    /// assert!(e.side_query::<TrivialKernel>(&Point2::new(1.0, 0.0)).is_on_right_side());
    /// assert!(e.side_query::<TrivialKernel>(&Point2::new(0.0, 1.0)).is_on_left_side());
    /// assert!(e.side_query::<TrivialKernel>(&Point2::new(0.5, 0.5)).is_on_line());
    /// # }
    /// ```
    pub fn side_query<K: DelaunayKernel<V::Scalar>>(&self, q: &V) -> EdgeSideInfo<V::Scalar> {
        K::side_query(&self, q)
    }

    pub fn intersects_edge_non_colinear<K>(&self, other: &SimpleEdge<V>) -> bool 
        where K: DelaunayKernel<V::Scalar>
    {
        let other_from = self.side_query::<K>(&other.from);
        let other_to = self.side_query::<K>(&other.to);
        let self_from = other.side_query::<K>(&self.from);
        let self_to = other.side_query::<K>(&self.to);

        debug_assert!(![&other_from, &other_to, &self_from, &self_to].iter()
                      .all(|q| q.is_on_line()), 
                      "Given edge is colinear.");

        other_from != other_to && self_from != self_to
    }
}

impl <V> SimpleEdge<V> where V: PointN, V::Scalar: SpadeFloat {
    /// Yields the nearest point on this edge.
    pub fn nearest_point(&self, query_point: &V) -> V {
        let (p1, p2) = (&self.from, &self.to);
        let dir = p2.sub(p1);
        let s = self.project_point(query_point);
        if V::Scalar::zero() < s && s < one() {
            p1.add(&dir.mul(s))

        } else {
            if s <= V::Scalar::zero() {
                p1.clone()
            } else {
                p2.clone()
            }
        }
    }

    /// Returns the squared distance of a given point to its
    /// projection onto the infinite line going through this edge's start 
    /// and end point.
    pub fn projection_distance2(&self, query_point: &V) -> V::Scalar {
        let s = self.project_point(query_point);
        let p = self.from.add(&self.to.sub(&self.from).mul(s));
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
        let (ref p1, ref p2) = (self.from.clone(), self.to.clone());
        let dir = p2.sub(p1);
        query_point.sub(p1).dot(&dir) / dir.length2()
    }
}

impl <V: PointN> SpatialObject for SimpleEdge<V> where V::Scalar: SpadeFloat {
    type Point = V;

    fn mbr(&self) -> BoundingRect<V> {
        BoundingRect::from_corners(&self.from, &self.to)
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        let nn = self.nearest_point(point);
        point.sub(&nn).length2()
    }
}

/// A triangle, defined by it's three points.
#[derive(Clone)]
pub struct SimpleTriangle<V: PointN> {
    v0: V,
    v1: V,
    v2: V,
}

impl <V> SimpleTriangle<V> where V: PointN {
    /// Checks if the given points are ordered counter clock wise.
    pub fn new(v0: V, v1: V, v2: V) -> SimpleTriangle<V> {
        SimpleTriangle { v0: v0, v1: v1, v2: v2 }
    }

    /// Returns the triangle's vertices.
    pub fn vertices(&self) -> [&V; 3] {
        [&self.v0, &self.v1, &self.v2]
    }
}

impl <V: TwoDimensional> SimpleTriangle<V> where V: TwoDimensional {

    /// Returns the triangle's doubled area.
    pub fn double_area(&self) -> V::Scalar {
        let b = self.v1.sub(&self.v0);
        let c = self.v2.sub(&self.v0);
        (b.nth(0).clone() * c.nth(1).clone() - b.nth(1).clone() * c.nth(0).clone()).abs()
    }

}

impl <V> PartialEq for SimpleTriangle<V> where V: PointN {
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

impl <V> SimpleTriangle<V> where V: PointN, V::Scalar: SpadeFloat {

    /// Returns the nearest point lying on any of the triangle's edges.
    pub fn nearest_point_on_edge(&self, pos: &V) -> V {
        let e0 = SimpleEdge::new(self.v0.clone(), self.v1.clone());
        let e1 = SimpleEdge::new(self.v1.clone(), self.v2.clone());
        let e2 = SimpleEdge::new(self.v2.clone(), self.v0.clone());
        let p0 = e0.nearest_point(pos);
        let p1 = e1.nearest_point(pos);
        let p2 = e2.nearest_point(pos);
        let d0 = p0.distance2(pos);
        let d1 = p1.distance2(pos);
        let d2 = p2.distance2(pos);
        if d0 <= d1 && d0 <= d2 {
            return p0;
        }
        if d1 <= d0 && d1 <= d2 {
            return p1;
        }
        p2
    }
}

impl <V> SimpleTriangle<V> where V: TwoDimensional, V::Scalar: SpadeFloat {

    /// Returns the position of the triangle's circumcenter.
    pub fn circumcenter(&self) -> V {
        let one: V::Scalar = One::one();
        let two = one + one;
        let b = self.v1.sub(&self.v0);
        let c = self.v2.sub(&self.v0);
        // Calculate circumcenter position
        let d = two * (*b.nth(0) * *c.nth(1) - *c.nth(0) * *b.nth(1));
        let len_b = b.dot(&b);
        let len_c = c.dot(&c);
        let x = (len_b * *c.nth(1) - len_c * *b.nth(1)) / d;
        let y = (-len_b * *c.nth(0) + len_c * *b.nth(0)) / d;
        let mut result = V::new();
        *result.nth_mut(0) = x;
        *result.nth_mut(1) = y;
        result.add(&self.v0)
    }

    /// Returns the barycentric coordinates of a point.
    pub fn barycentric_interpolation(&self, coord: &V) -> Point3<V::Scalar> {
        let (v1, v2, v3) = (self.v0.clone(), self.v1.clone(), self.v2.clone());
        let (x, y) = (*coord.nth(0), *coord.nth(1));
        let (x1, x2, x3) = (*v1.nth(0), *v2.nth(0), *v3.nth(0));
        let (y1, y2, y3) = (*v1.nth(1), *v2.nth(1), *v3.nth(1));
        let det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        let lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
        let lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
        let lambda3 = one::<V::Scalar>() - lambda1 - lambda2;
        Point3::new(lambda1, lambda2, lambda3)
    }
}


impl <V> SpatialObject for SimpleTriangle<V> where V: TwoDimensional, V::Scalar: SpadeFloat {
    type Point = V;
    
    fn mbr(&self) -> BoundingRect<V> {
        let mut result = BoundingRect::from_corners(&self.v0, &self.v1);
        result.add_point(self.v2.clone());
        result
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        
        let ordered_ccw = TrivialKernel::is_ordered_ccw(
            &self.v0, &self.v1, &self.v2);
        for i in 0 .. 3 {
            let edge = SimpleEdge::new(self.vertices()[i].clone(), 
                                       self.vertices()[(i + 1) % 3].clone());
            if edge.side_query::<TrivialKernel>(point).is_on_right_side() == ordered_ccw {
                return edge.distance2(point);
            }
        }
        // The point lies within the triangle
        zero()
    }
}

/// An n-dimensional circle, defined by its origin and radius.
pub struct SimpleCircle<V: PointN> {
    /// The circle's center.
    pub center: V,
    /// The circle's radius.
    pub radius: V::Scalar,
}

impl <V> SimpleCircle<V> where V: PointN, V::Scalar: SpadeFloat {

    /// Create a new circle.
    pub fn new(center: V, radius: V::Scalar) -> SimpleCircle<V> {
        SimpleCircle {
            center: center,
            radius: radius,
        }
    }
}

impl <V> SpatialObject for SimpleCircle<V> where V: PointN, V::Scalar: SpadeFloat {
    type Point = V;

    fn mbr(&self) -> BoundingRect<V> {
        let r = V::from_value(self.radius);
        BoundingRect::from_corners(&self.center.sub(&r), &self.center.add(&r))
    }

    fn distance2(&self, point: &V) -> V::Scalar {
        let dx = *self.center.nth(0) - *point.nth(0);
        let dy = *self.center.nth(1) - *point.nth(1);
        let dist = ((dx * dx + dy * dy).sqrt() - self.radius).max(zero());
        dist * dist
    }

    // Since containment checks do not require the calculation of the square root
    // we can redefine this method.
    fn contains(&self, point: &V) -> bool {
        let dx = *self.center.nth(0) - *point.nth(0);
        let dy = *self.center.nth(1) - *point.nth(1);
        let r2 = self.radius * self.radius;
        dx * dx + dy * dy <= r2
    }
}

#[cfg(test)]
mod test {
    use super::{SimpleEdge, SimpleTriangle};
    use traits::SpatialObject;
    use kernels::{TrivialKernel, FloatKernel};
    use cgmath::{Point2};
    use approx::ApproxEq;

    #[test]
    fn test_edge_distance() {
        let e = SimpleEdge::new(Point2::new(0f32, 0.), Point2::new(1., 1.));
        assert!(e.distance2(&Point2::new(1.0, 0.0)).relative_eq(&0.5, 1e-10, 1e-10));
        assert!(e.distance2(&Point2::new(0.0, 1.)).relative_eq(&0.5, 1e-10, 1e-10));
        assert!(e.distance2(&Point2::new(-1.0, -1.)).relative_eq(&2., 1e-10, 1e-10));
        assert!(e.distance2(&Point2::new(2.0, 2.0)).relative_eq(&2., 1e-10, 1e-10));
    }

    #[test]
    fn test_edge_side() {
        let e = SimpleEdge::new(Point2::new(0f32, 0.), Point2::new(1., 1.));
        assert!(e.side_query::<TrivialKernel>(&Point2::new(1.0, 0.0)).is_on_right_side());
        assert!(e.side_query::<TrivialKernel>(&Point2::new(0.0, 1.0)).is_on_left_side());
        assert!(e.side_query::<TrivialKernel>(&Point2::new(0.5, 0.5)).is_on_line());
    }

    #[test]
    fn test_intersects_middle() {
        let e1 = SimpleEdge::new(Point2::new(0f32, 0f32), Point2::new(5f32, 5f32));
        let e2 = SimpleEdge::new(Point2::new(-1.5, 1.), Point2::new(1.0, -1.5));
        let e3 = SimpleEdge::new(Point2::new(0.5, 4.), Point2::new(0.5, -4.));
        assert!(!e1.intersects_edge_non_colinear::<TrivialKernel>(&e2));
        assert!(!e2.intersects_edge_non_colinear::<TrivialKernel>(&e1));
        assert!(e1.intersects_edge_non_colinear::<TrivialKernel>(&e3));
        assert!(e3.intersects_edge_non_colinear::<TrivialKernel>(&e1));
        assert!(e2.intersects_edge_non_colinear::<TrivialKernel>(&e3));
        assert!(e3.intersects_edge_non_colinear::<TrivialKernel>(&e2));
    }
    
    #[test]
    fn test_intersects_end_points() {
        // Check for intersection of one endpoint touches another edge
        let e1 = SimpleEdge::new(Point2::new(0.33f64, 0.33f64), Point2::new(1.0, 0.0));
        let e2 = SimpleEdge::new(Point2::new(0.33, -1.0), Point2::new(0.33, 1.0));
        assert!(e1.intersects_edge_non_colinear::<FloatKernel>(&e2));
        assert!(e2.intersects_edge_non_colinear::<FloatKernel>(&e1));
        let e3 = SimpleEdge::new(Point2::new(0.0, -1.0), Point2::new(2.0, 1.0));
        assert!(e1.intersects_edge_non_colinear::<FloatKernel>(&e3));
        assert!(e3.intersects_edge_non_colinear::<FloatKernel>(&e1));
        // Check for intersection if only end points overlap
        let e4 = SimpleEdge::new(Point2::new(0.33, 0.33), Point2::new(0.0, 2.0));
        assert!(e1.intersects_edge_non_colinear::<FloatKernel>(&e4));
        assert!(e4.intersects_edge_non_colinear::<FloatKernel>(&e1));
        // // Check for intersection if edges are colinear
        // let e5 = SimpleEdge::new(Point2::new(1.0, 0.0), Point2::new(3.0, 2.0));
        // let e6 = SimpleEdge::new(Point2::new(-1., -2.), Point2::new(-2., -3.));
        // assert!(e3.intersects_edge_non_colinear::<FloatKernel>(&e5));
        // assert!(e5.intersects_edge_non_colinear::<FloatKernel>(&e3));
        // // These edges are colinear but not intersecting
        // assert!(!e3.intersects_edge_non_colinear::<FloatKernel>(&e6));
        // assert!(!e6.intersects_edge_non_colinear::<FloatKernel>(&e3));
    }

    #[test]
    fn test_triangle_distance() {
        let v1 = Point2::new(0f32, 0.);
        let v2 = Point2::new(1., 0.);
        let v3 = Point2::new(0., 1.);
        let t = SimpleTriangle::new(v1, v2, v3);
        assert_eq!(t.distance2(&Point2::new(0.25, 0.25)), 0.);
        assert!(t.distance2(&Point2::new(-1., -1.)).relative_eq(&2., 1e-10, 1e-10));
        assert!(t.distance2(&Point2::new(0., -1.)).relative_eq(&1., 1e-10, 1e-10));
        assert!(t.distance2(&Point2::new(-1., 0.)).relative_eq(&1., 1e-10, 1e-10));
        assert!(t.distance2(&Point2::new(1., 1.)).relative_eq(&0.5, 1e-10, 1e-10));
        assert!(t.distance2(&Point2::new(0.5, 0.5)).relative_eq(&0.0, 1e-10, 1e-10));
        assert!(t.distance2(&Point2::new(0.6, 0.6)) > 0.001);
    }
}
