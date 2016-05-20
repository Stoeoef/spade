use cgmath::{Vector2, Vector3, BaseFloat, EuclideanVector};
use traits::{SpatialObject, SpatialEdge};
use num::{zero, one};
use boundingvolume::BoundingRect;


#[derive(Clone)]
pub struct SimpleEdge<S: BaseFloat> {
    pub from: Vector2<S>,
    pub to: Vector2<S>,
}

impl <S: BaseFloat> SimpleEdge<S> {
    pub fn new(from: Vector2<S>, to: Vector2<S>) -> SimpleEdge<S> {
        SimpleEdge {
            from: from,
            to: to,
        }
    }
}

impl <S: BaseFloat> SpatialEdge for SimpleEdge<S> {

    fn from(&self) -> Vector2<S> {
        self.from
    }

    fn to(&self) -> Vector2<S> {
        self.to
    }
}

impl <S: BaseFloat> SpatialObject for SimpleEdge<S> {
    
    type Scalar = S;

    fn mbr(&self) -> BoundingRect<S> {
        BoundingRect::from_corners(&self.from, &self.to)
    }

    fn distance2(&self, point: &Vector2<S>) -> S {
        (point - self.nearest_point(point)).length2()
    }
}



// impl <T> SpatialObject for T where T: SpatialEdge {
//     type Scalar = T::Scalar;
//     fn mbr(&self) -> BoundingRect<T::Scalar> {
//         BoundingRect::from_corners(&self.from(), &self.to())
//     }

//     fn distance2(&self, point: &Vector2<T::Scalar>) -> T::Scalar {
//         (point - self.nearest_point(point)).length2()
//     }
// }

pub trait SpatialTriangle<S: BaseFloat> {
    fn vertices(&self) -> [Vector2<S>; 3];
}

#[derive(Clone)]
pub struct SimpleTriangle<S: BaseFloat> {
    pub vertices: [Vector2<S>; 3],
}

impl <S> SimpleTriangle<S> where S: BaseFloat {
    pub fn is_ordered_ccw(vertices: &[Vector2<S>; 3]) -> bool {
        let edge = SimpleEdge::new(vertices[0].clone(), vertices[1].clone());
        return edge.signed_side(&vertices[2]) >= zero()
    }

    pub fn new(mut vertices: [Vector2<S>; 3]) -> SimpleTriangle<S> {
        if !SimpleTriangle::is_ordered_ccw(&vertices) {
            vertices.swap(1, 2);
        }
        SimpleTriangle::new_ordered(vertices)
    }

    pub fn new_ordered(vertices: [Vector2<S>; 3]) -> SimpleTriangle<S> {
        assert!(SimpleTriangle::is_ordered_ccw(&vertices));
        SimpleTriangle {
            vertices: vertices,
        }
    }
}

impl <S> PartialEq for SimpleTriangle<S> where S: BaseFloat {
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


impl <S> SpatialObject for SimpleTriangle<S> where S: BaseFloat {
    type Scalar = S;
    
    fn mbr(&self) -> BoundingRect<S> {
        let mut result = BoundingRect::from_corners(&self.vertices[0], &self.vertices[1]);
        result.add_point(&self.vertices[2]);
        result
    }

    fn distance2(&self, point: &Vector2<S>) -> S {
        for i in 0 .. 3 {
            let edge = SimpleEdge::new(self.vertices[i], self.vertices[(i + 1) % 3]);
            if edge.signed_side(point) <= zero() {
                return edge.distance2(point);
            }
        }
        // The point is within the triangle
        zero()
    }

}

impl <S> SpatialTriangle<S> where S: BaseFloat {
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

#[cfg(test)]
mod test {
    use super::{SimpleEdge, SimpleTriangle};
    use traits::{SpatialObject, SpatialEdge};
    use cgmath::{Vector2, ApproxEq};
    
    #[test]
    fn test_edge_distance2() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.distance2(&Vector2::new(1.0, 0.0)).approx_eq(&0.5));
        assert!(e.distance2(&Vector2::new(0.0, 1.)).approx_eq(&0.5));
        assert!(e.distance2(&Vector2::new(-1.0, -1.)).approx_eq(&2.));
        assert!(e.distance2(&Vector2::new(2.0, 2.0)).approx_eq(&2.));
    }

    #[test]
    fn test_edge_side() {
        let e = SimpleEdge::new(Vector2::new(0f32, 0.), Vector2::new(1., 1.));
        assert!(e.signed_side(&Vector2::new(1.0, 0.0)) < 0.0);
        assert!(e.signed_side(&Vector2::new(0.0, 1.0)) > 0.0);
        assert!(e.signed_side(&Vector2::new(0.5, 0.5)).approx_eq(&0.0));
    }

    #[test]
    fn test_triangle_distance2() {
        let v1 = Vector2::new(0f32, 0.);
        let v2 = Vector2::new(1., 0.);
        let v3 = Vector2::new(0., 1.);
        let t = SimpleTriangle::new([v1, v2, v3]);
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
        assert!(SimpleTriangle::is_ordered_ccw(&[v1, v2, v3]));
        assert!(!SimpleTriangle::is_ordered_ccw(&[v2, v1, v3]));
        assert!(SimpleTriangle::is_ordered_ccw(&[v3, v1, v2]));
    }        

}
