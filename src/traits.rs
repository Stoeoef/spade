use cgmath::{BaseFloat, Vector2, Vector, EuclideanVector, one, zero};
use num::traits::Zero;
use boundingvolume::BoundingRect;

pub trait SpatialObject {
    type Scalar: BaseFloat;
    fn mbr(&self) -> BoundingRect<Self::Scalar>;

    fn distance2(&self, point: &Vector2<Self::Scalar>) -> Self::Scalar {
        self.mbr().min_dist(point)
    }

    fn contains(&self, point: &Vector2<Self::Scalar>) -> bool {
        self.distance2(point) == zero()
    }
}

pub trait SpatialEdge : SpatialObject {

    fn from(&self) -> Vector2<Self::Scalar>;
    fn to(&self) -> Vector2<Self::Scalar>;

    fn nearest_point(&self, point: &Vector2<Self::Scalar>) -> Vector2<Self::Scalar> {
        let (p1, p2) = (self.from(), self.to());
        let dir = p2 - p1;
        let s: Self::Scalar = (point - p1).dot(dir) / dir.length2();
        if Self::Scalar::zero() < s && s < one() {
            p1 + dir * s
        } else {
            if s <= zero() {
                p1
            } else {
                p2
            }
        }
    }

    /**
     * Returns a value smaller than zero if point is on the right side of this edge,
     * zero if c is on the edge and a positive value if point is to the left.
     **/
    fn signed_side(&self, point: &Vector2<Self::Scalar>) -> Self::Scalar {
        let (ref a, ref b) = (self.from(), self.to());
        (b.x - a.x) * (point.y - a.y) - (b.y - a.y)  * (point.x - a.x)
    }

    fn is_on_right_side(&self, point: &Vector2<Self::Scalar>) -> bool {
        self.signed_side(point) < zero()
    }

    fn is_on_left_side(&self, point: &Vector2<Self::Scalar>) -> bool {
        self.signed_side(point) > zero()
    }
}


pub trait PointObject : SpatialObject {
    fn position(&self) -> Vector2<Self::Scalar>;
}


impl <S: BaseFloat> PointObject for Vector2<S> {

    fn position(&self) -> Vector2<S> {
        self.clone()
    }
}

impl <S: BaseFloat> SpatialObject for Vector2<S> {
    type Scalar = S;
    fn mbr(&self) -> BoundingRect<S> {
        BoundingRect::from_point(self.position())
    }

    fn distance2(&self, point: &Vector2<S>) -> S {
        // Overwrite default implementation since this will be faster
        (point - self.position()).length2()
    }
}
