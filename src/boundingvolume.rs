use cgmath::{Vector2, BaseFloat, EuclideanVector};
use num::{Float, zero, cast};
use misc::{fmin, fmax, fclamp};

#[derive(Clone, PartialEq)]
pub struct BoundingRect<S: BaseFloat> {
    lower: Vector2<S>,
    upper: Vector2<S>,
}

impl <S: BaseFloat> BoundingRect<S> {
    #[inline]
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
        let rl = rect.lower();
        let ru = rect.upper();
        self.lower.x <= rl.x && self.lower.y <= rl.y
            && self.upper.x >= ru.x && self.upper.y >= ru.y
    }

    #[inline]
    pub fn add_point(&mut self, point: &Vector2<S>) {
        self.lower = Vector2::new(*fmin(&self.lower.x, &point.x),
                                  *fmin(&self.lower.y, &point.y));
        self.upper = Vector2::new(*fmax(&self.upper.x, &point.x),
                                  *fmax(&self.upper.y, &point.y));
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
        *fmax(&(diag.x * diag.y), &zero())
    }

    pub fn half_margin(&self) -> S {
        let diag = self.upper - self.lower;
        *fmax(&(diag.x + diag.y), &zero())
    }

    pub fn center(&self) -> Vector2<S> {
        self.lower + (self.upper - self.lower) / cast::<_, S>(2i32).unwrap()
    }

    pub fn intersect(&self, other: &BoundingRect<S>) -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(*fmax(&self.lower.x, &other.lower.x),
                                *fmax(&self.lower.y, &other.lower.y)),
            upper: Vector2::new(*fmin(&self.upper.x, &other.upper.x),
                                *fmin(&self.upper.y, &other.upper.y)),
        }
    }

    pub fn min_dist(&self, point: &Vector2<S>) -> S {
        let l = self.lower();
        let u = self.upper();
        let x = *fclamp(&l.x, &u.x, &point.x);
        let y = *fclamp(&l.y, &u.y, &point.y);
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
