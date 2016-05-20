use cgmath::{Vector2, BaseFloat, EuclideanVector};
use num::{Float, zero, cast};

#[derive(Clone, PartialEq)]
pub struct BoundingRect<S: BaseFloat> {
    lower: Vector2<S>,
    upper: Vector2<S>,
}

// A call to l.min(r) does not seem to be inlined, thus we define it ourselves
// This does improve performance significantly, especially for larger node sizes
#[inline]
fn fmin<'a, S: BaseFloat>(l: &'a S, r: &'a S) -> &'a S {
    if l < r {
        l
    } else {
        r
    }
}

#[inline]
fn fmax<'a, S: BaseFloat>(l: &'a S, r: &'a S) -> &'a S {
    if l > r {
        l
    } else {
        r
    }
}

impl <S: BaseFloat> BoundingRect<S> {
    pub fn new() -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(Float::infinity(), Float::infinity()),
            upper: Vector2::new(Float::neg_infinity(), Float::neg_infinity()),
        }
    }

    pub fn from_point(point: Vector2<S>) -> BoundingRect<S> {
        BoundingRect {
            lower: point,
            upper: point,
        }
    }

    pub fn from_corners(corner1: &Vector2<S>, corner2: &Vector2<S>) -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(corner1.x.min(corner2.x), corner1.y.min(corner2.y)),
            upper: Vector2::new(corner1.x.max(corner2.x), corner1.y.max(corner2.y)),
        }
    }

    #[inline]
    pub fn lower(&self) -> Vector2<S> {
        self.lower
    }

    #[inline]
    pub fn upper(&self) -> Vector2<S> {
        self.upper
    }

    #[inline]
    pub fn contains_point(&self, point: &Vector2<S>) -> bool {
        self.lower.x <= point.x && self.lower.y <= point.y
            && self.upper.x >= point.x && self.upper.y >= point.y
    }
    
    #[inline]
    pub fn contains_rect(&self, rect: &BoundingRect<S>) -> bool {
        self.contains_point(&rect.lower()) && self.contains_point(&rect.upper())
    }

    #[inline]
    pub fn add_point(&mut self, point: &Vector2<S>) {
        self.lower = Vector2::new(self.lower.x.min(point.x), 
                                  self.lower.y.min(point.y));
        self.upper = Vector2::new(self.upper.x.max(point.x),
                                  self.upper.y.max(point.y));
    }

    #[inline]
    pub fn add_rect(&mut self, rect: &BoundingRect<S>) {
        self.lower = Vector2::new(*fmin(&self.lower.x, &rect.lower.x),
                                  *fmin(&self.lower.y, &rect.lower.y));
        self.upper = Vector2::new(*fmax(&self.upper.x, &rect.upper.x),
                                  *fmax(&self.upper.y, &rect.upper.y));
    }

    pub fn area(&self) -> S {
        let diag = self.upper - self.lower;
        (diag.x * diag.y).max(zero())
    }

    pub fn half_margin(&self) -> S {
        let diag = self.upper - self.lower;
        (diag.x + diag.y).max(zero())
    }

    pub fn center(&self) -> Vector2<S> {
        self.lower + (self.upper - self.lower) / cast::<_, S>(2i32).unwrap()
    }

    pub fn intersect(&self, other: &BoundingRect<S>) -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(self.lower.x.max(other.lower.x), 
                                self.lower.y.max(other.lower.y)),
            upper: Vector2::new(self.upper.x.min(other.upper.x),
                                self.upper.y.min(other.upper.y)),
        }
    }

    pub fn min_dist(&self, point: &Vector2<S>) -> S {
        let l = self.lower();
        let u = self.upper();
        let x = u.x.min(point.x.max(l.x));
        let y = u.y.min(point.y.max(l.y));
        (Vector2::new(x, y) - point).length2()
    }

    pub fn min_max_dist(&self, point: &Vector2<S>) -> S {
        let l = self.lower();
        let u = self.upper();
        let (min_x, max_x) = if (l.x - point.x).abs() < (u.x - point.x).abs() { (l.x, u.x) } else { (u.x, l.x) };
        let (min_y, max_y) = if (l.y - point.y).abs() < (u.y - point.y).abs() { (l.y, u.y) } else { (u.y, l.y) };
        let p1 = Vector2::new(min_x, max_y);
        let p2 = Vector2::new(max_x, min_y);
        (p1 - point).length2().max((p2 - point).length2())
    }
}
