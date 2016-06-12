use cgmath::{BaseNum, Vector2};
use cgmath::conv::array2;
use num::{Bounded, zero, cast};
use misc::{length2, clamp};
use traits::RTreeNum;

/// An axis aligned minimal bounding rectangle.
///
/// A axis aligned minimal bounding rectangle is the smallest rectangle that completly surrounds an object and is aligned along the two axis in the plane.
#[derive(Clone, PartialEq, Debug)]
pub struct BoundingRect<S: BaseNum + RTreeNum> {
    lower: Vector2<S>,
    upper: Vector2<S>,
}

impl <S: BaseNum + RTreeNum> BoundingRect<S> {

    /// Creates an uninitialized bounding rectangle.
    ///
    /// This will create a new bounding rectangle that contains no objects.
    /// Note that `lower` and `upper` are set to positive / negative infinity,
    /// thus some operations (like `area`, `margin`) will yield unspecified results.
    pub fn new() -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(Bounded::max_value(), Bounded::max_value()),
            upper: Vector2::new(Bounded::min_value(), Bounded::min_value()),
        }
    }

    /// Creates a bounding rectangle that contains exactly one point.
    ///
    /// This will create a bounding rectangle with `lower == upper == point`.
    pub fn from_point(point: [S; 2]) -> BoundingRect<S> {
        BoundingRect {
            lower: point.into(),
            upper: point.into(),
        }
    }

    /// Creates a bounding rectangle that contains two points.
    pub fn from_corners(corner1: &[S; 2], corner2: &[S; 2]) -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(corner1[0].partial_min(corner2[0]), corner1[1].partial_min(corner2[1])),
            upper: Vector2::new(corner1[0].partial_max(corner2[0]), corner1[1].partial_max(corner2[1])),
        }
    }

    /// Returns the lower corner of the bounding rectangle.
    ///
    /// The lower corner has the smaller coordinates.
    pub fn lower(&self) -> [S; 2] {
        array2(self.lower)
    }

    /// Returns the upper corner of the bounding rectangle.
    ///
    /// The upper corner has the larger coordinates.
    pub fn upper(&self) -> [S; 2] {
        array2(self.upper)
    }

    /// Checks if a point is contained within the bounding rectangle.
    ///
    /// A point lying exactly on the bounding rectangle's border is also contained.
    #[inline]
    pub fn contains_point(&self, point: &[S; 2]) -> bool {
        self.lower.x <= point[0] && self.lower.y <= point[1]
            && self.upper.x >= point[0] && self.upper.y >= point[1]
    }
    
    /// Checks if another bounding rectangle is completley contained withing this rectangle.
    ///
    /// A rectangle is contained if and only if all four corner are contained (see `contains_point`).
    #[inline]
    pub fn contains_rect(&self, rect: &BoundingRect<S>) -> bool {
        let rl = rect.lower;
        let ru = rect.upper;
        self.lower.x <= rl.x && self.lower.y <= rl.y
            && self.upper.x >= ru.x && self.upper.y >= ru.y
    }

    /// Enlarges this bounding rectangle to contain a point.
    ///
    /// If the point is already contained, nothing will be changed.
    /// Otherwise, this will enlarge `self` to be just large enough
    /// to contain the new point.
    #[inline]
    pub fn add_point(&mut self, point: &[S; 2]) {
        self.lower = Vector2::new(self.lower.x.partial_min(point[0]),
                                  self.lower.y.partial_min(point[1]));
        self.lower = Vector2::new(self.upper.x.partial_max(point[0]),
                                  self.upper.y.partial_max(point[1]));
    }

    /// Enlarges this bounding rectangle to contain a rectangle.
    ///
    /// If the rectangle is already contained, nothing will be changed.
    /// Otherwise, this will enlarge `self` to be just large enough
    // to contain the new rectangle.
    #[inline]
    pub fn add_rect(&mut self, rect: &BoundingRect<S>) {
        self.lower = Vector2::new(self.lower.x.partial_min(rect.lower.x),
                                  self.lower.y.partial_min(rect.lower.y));
        self.upper = Vector2::new(self.upper.x.partial_max(rect.upper.x),
                                  self.upper.y.partial_max(rect.upper.y));
    }

    /// Returns the rectangle's area.
    pub fn area(&self) -> S {
        let diag = self.upper - self.lower;
        (diag.x * diag.y).partial_max(zero())
    }

    /// Returns half of the rectangle's margin, thus `width + height`.
    pub fn half_margin(&self) -> S {
        let diag = self.upper - self.lower;
        (diag.x + diag.y).partial_max(zero())
    }

    /// Returns the rectangle's center.
    pub fn center(&self) -> Vector2<S> {
        self.lower + (self.upper - self.lower) / cast::<_, S>(2i32).unwrap()
    }

    /// Returns the intersection of this and another bounding rectangle.
    ///
    /// If the rectangles do not intersect, a bounding rectangle with an area and
    /// margin of zero is returned.
    pub fn intersect(&self, other: &BoundingRect<S>) -> BoundingRect<S> {
        BoundingRect {
            lower: Vector2::new(self.lower.x.partial_max(other.lower.x),
                                self.lower.y.partial_max(other.lower.y)),
            upper: Vector2::new(self.upper.x.partial_min(other.upper.x),
                                self.upper.y.partial_min(other.upper.y)),
        }
    }

    #[doc(hidden)]
    pub fn min_dist2(&self, point: &Vector2<S>) -> S {
        let l = self.lower;
        let u = self.upper;
        let x = clamp(l.x, u.x, point.x);
        let y = clamp(l.y, u.y, point.y);
        length2(&(Vector2::new(x, y) - point))
    }

    #[doc(hidden)]
    pub fn min_max_dist2(&self, point: &Vector2<S>) -> S {
        let l = self.lower;
        let u = self.upper;
        let (min_x, max_x) = if (l.x - point.x).abs() < (u.x - point.x).abs() { (l.x, u.x) } else { (u.x, l.x) };
        let (min_y, max_y) = if (l.y - point.y).abs() < (u.y - point.y).abs() { (l.y, u.y) } else { (u.y, l.y) };
        let p1 = Vector2::new(min_x, max_y);
        let p2 = Vector2::new(max_x, min_y);
        length2(&(p1 - point)).partial_max(length2(&(p2 - point)))
    }
}
