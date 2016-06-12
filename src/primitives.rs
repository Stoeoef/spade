//! Contains some useful primitives that can be inserted into r-trees.
//!
//! Use these objects if only the geometrical properties (position and size)
//! are important. If additional data needs to be stored per object, consider
//! implementing `SpatialObject`.


use cgmath::{Vector2, Vector3, BaseFloat, InnerSpace};
use cgmath::conv::array2;
use traits::{RTreeNum, SpatialObject};
use num::{Bounded, zero, one};
use boundingvolume::BoundingRect;
use misc::{length2};


/// A two dimensional point.
#[derive(Clone, PartialEq, Debug)]
pub struct SimplePoint<S: RTreeNum> {
    coords: Vector2<S>
}

impl <S: RTreeNum + Bounded> SimplePoint<S> {
    pub fn new(coords: [S; 2]) -> SimplePoint<S> {
        SimplePoint { coords: coords.into() }
    }

    pub fn coords(&self) -> [S; 2] {
        array2(self.coords)
    }
}

impl <S: RTreeNum + Bounded> SpatialObject for SimplePoint<S> {
    type Scalar = S;

    fn mbr(&self) -> BoundingRect<S> {
        BoundingRect::from_point(array2(self.coords))
    }

    fn distance(&self, point: [S; 2]) -> S {
        let point: Vector2<_> = point.clone().into();
        length2(&(point - self.coords))
    }
}
        
#[derive(Clone, Debug)]
/// An edge defined by it's two end points.
pub struct SimpleEdge<S: BaseFloat + Bounded> {
    /// The edge's origin.
    pub from: Vector2<S>,
    /// The edge's destination.
    pub to: Vector2<S>,
}

impl <S: BaseFloat + Bounded> SimpleEdge<S> {
    /// Creates a new edge from `from` to `to`.
    pub fn new(from: Vector2<S>, to: Vector2<S>) -> SimpleEdge<S> {
        SimpleEdge {
            from: from,
            to: to,
        }
    }

    /// Yields the nearest point on this edge.
    pub fn nearest_point(&self, query_point: [S; 2]) -> [S; 2] {
        let (p1, p2) = (self.from, self.to);
        let dir = p2 - p1;
        let point: Vector2<_> = query_point.into();
        let s: S = (point - p1).dot(dir) / dir.magnitude2();
        if S::zero() < s && s < one() {
            array2(p1 + dir * s)
        } else {
            if s <= zero() {
                array2(p1)
            } else {
                array2(p2)
            }
        }
    }
    /// Returns a value indicating on which side of this edge a given point lies.
    ///
    /// Returns a value smaller than zero if `query_point` is on the right side,
    /// or greater than zero if `query_point` is on the left side.
    pub fn signed_side(&self, query_point: [S; 2]) -> S {
        let (ref a, ref b) = (self.from, self.to);
        let point: Vector2<_> = query_point.into();
        (b.x - a.x) * (point.y - a.y) - (b.y - a.y)  * (point.x - a.x)
    }

    /// Returns `true` if a given point is on the right side of this edge.
    fn is_on_right_side(&self, query_point: [S; 2]) -> bool {
        self.signed_side(query_point) < zero()
    }

    /// Returns `false` if a given point is on the left side of this edge.
    fn is_on_left_side(&self, query_point: [S; 2]) -> bool {
        self.signed_side(query_point) > zero()
    }
}

impl <S: BaseFloat + RTreeNum> SpatialObject for SimpleEdge<S> {
    type Scalar = S;

    fn mbr(&self) -> BoundingRect<S> {
        BoundingRect::from_corners(&array2(self.from), &array2(self.to))
    }

    fn distance(&self, point: [S; 2]) -> S {
        let point_vec: Vector2<_> = point.into();
        let nn: Vector2<_> = self.nearest_point(point).into();
        (point_vec - nn).magnitude2()
    }
}

/// A triangle, defined by it's three points.
#[derive(Clone)]
pub struct SimpleTriangle<S: BaseFloat + RTreeNum> {
    vertices: [Vector2<S>; 3],
}

impl <S> SimpleTriangle<S> where S: BaseFloat + RTreeNum {
    /// Checks if the given points are ordered counter clock wise.
    pub fn is_ordered_ccw(vertices: &[Vector2<S>; 3]) -> bool {
        let edge = SimpleEdge::new(vertices[0].clone(), vertices[1].clone());
        return edge.signed_side(array2(vertices[2])) >= zero()
    }

    pub fn new(mut vertices: [Vector2<S>; 3]) -> SimpleTriangle<S> {
        if !SimpleTriangle::is_ordered_ccw(&vertices) {
            vertices.swap(1, 2);
        }
        SimpleTriangle { vertices: vertices }
    }

    pub fn vertices(&self) -> &[Vector2<S>; 3] {
        &self.vertices
    }

    pub fn barycentric_interpolation(&self, coord: &Vector2<S>) -> Vector3<S> {
        let vertices = self.vertices();
        let (v1, v2, v3) = (vertices[0], vertices[1], vertices[2]);
        let (x, y) = (coord.x, coord.y);
        let (x1, x2, x3) = (v1.x, v2.x, v3.x);
        let (y1, y2, y3) = (v1.y, v2.y, v3.y);
        let det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        let lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
        let lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
        let lambda3 = one::<S>() - lambda1 - lambda2;
        Vector3::new(lambda1, lambda2, lambda3)
    }
}

impl <S> PartialEq for SimpleTriangle<S> where S: BaseFloat + RTreeNum {
    fn eq(&self, rhs: &SimpleTriangle<S>) -> bool {
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


impl <S> SpatialObject for SimpleTriangle<S> where S: BaseFloat + RTreeNum {
    type Scalar = S;
    
    fn mbr(&self) -> BoundingRect<S> {
        let mut result = BoundingRect::new();
        for vertex in self.vertices.iter() {
            result.add_point(&array2(*vertex));
        }
        result
    }

    fn distance(&self, point: [S; 2]) -> S {
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

// struct SimpleCircle<S: RTreeNum> {
//     pub origin: [S; 2],
//     pub radius: S,
// }

// impl <S: RTreeNum> SpatialObject for Circle<S> {
//     type Scalar = S;

//     fn mbr(&self) -> BoundingRectangle<S> {
//         let (x, y) = (self.origin[0], self.origin[1]);
//         let r = self.radius;
//         BoundingRect::from_corners(&[x - r, y - r], &[x + r, y + r])
//     }

//     fn distance(&self, point: [S; 2]) -> S {
//         let dx = self.origin[0] - point[0];
//         let dy = self.origin[1] - point[1];
//         let r2 = self.radius * self.radius;
//         (dx * dx + dy * dy - r2).max(zero())
//     }
// }

#[cfg(test)]
mod test {
    use super::{SimpleEdge, SimpleTriangle};
    use traits::SpatialObject;
    use cgmath::{Vector2, ApproxEq};
    
    #[test]
    fn test_edge_distance() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.distance([1.0, 0.0]).approx_eq(&0.5));
        assert!(e.distance([0.0, 1.]).approx_eq(&0.5));
        assert!(e.distance([-1.0, -1.]).approx_eq(&2.));
        assert!(e.distance([2.0, 2.0]).approx_eq(&2.));
    }

    #[test]
    fn test_edge_side() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.signed_side([1.0, 0.0]) < 0.0);
        assert!(e.signed_side([0.0, 1.0]) > 0.0);
        assert!(e.signed_side([0.5, 0.5]).approx_eq(&0.0));
    }

    #[test]
    fn test_triangle_distance() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1., 0.);
        let v3 = Vector2::new(0., 1.);
        let t = SimpleTriangle::new([v1, v2, v3]);
        assert_eq!(t.distance([0.25, 0.25]), 0.);
        assert!(t.distance([-1., -1.]).approx_eq(&2.));
        assert!(t.distance([0., -1.]).approx_eq(&1.));
        assert!(t.distance([-1., 0.]).approx_eq(&1.));
        assert!(t.distance([1., 1.]).approx_eq(&0.5));
        assert!(t.distance([0.5, 0.5]).approx_eq(&0.0));
        assert!(t.distance([0.6, 0.6]) > 0.001);
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
