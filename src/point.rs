use num_traits::{Num, Signed};

#[cfg(feature = "serde")]
use serde_crate::{Deserialize, Serialize};

/// A coordinate type that can be used with a triangulation.
///
/// Internally, most calculations are performed after converting the type into a `f64`.
/// However, changing this to `f32` will reduce the required storage space.
///
/// This type should usually be either `f32` or `f64`.
pub trait SpadeNum: Num + PartialOrd + Into<f64> + Copy + Signed + std::fmt::Debug {}

impl<T> SpadeNum for T where T: Num + PartialOrd + Into<f64> + Copy + Signed + std::fmt::Debug {}

/// A two dimensional point.
///
/// This is the basic type used for defining positions.
#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Default)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
pub struct Point2<S> {
    /// The point's x coordinate
    pub x: S,
    /// The point's y coordinate
    pub y: S,
}

impl<S> Point2<S> {
    /// Creates a new point.
    #[inline]
    pub const fn new(x: S, y: S) -> Self {
        Point2 { x, y }
    }
}

impl<S: SpadeNum> Point2<S> {
    /// Returns the squared distance of this point and another point.
    #[inline]
    pub fn point_distance_2(&self, other: Point2<S>) -> S {
        self.sub(other).length2()
    }

    pub(crate) fn mul(&self, factor: S) -> Self {
        Point2 {
            x: self.x * factor,
            y: self.y * factor,
        }
    }

    pub(crate) fn add(&self, other: Self) -> Self {
        Point2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    pub(crate) fn distance2(&self, other: Self) -> S {
        let delta_x = self.x - other.x;
        let delta_y = self.y - other.y;
        delta_x * delta_x + delta_y * delta_y
    }

    pub(crate) fn length2(&self) -> S {
        self.x * self.x + self.y * self.y
    }

    pub(crate) fn sub(&self, other: Self) -> Self {
        Point2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    pub(crate) fn dot(&self, other: Self) -> S {
        self.x * other.x + self.y * other.y
    }
}

impl<S: SpadeNum> From<Point2<S>> for [S; 2] {
    #[inline]
    fn from(point: Point2<S>) -> Self {
        [point.x, point.y]
    }
}

impl<S: SpadeNum> From<Point2<S>> for (S, S) {
    #[inline]
    fn from(point: Point2<S>) -> (S, S) {
        (point.x, point.y)
    }
}

impl<S: SpadeNum> From<[S; 2]> for Point2<S> {
    #[inline]
    fn from(source: [S; 2]) -> Self {
        Self::new(source[0], source[1])
    }
}

impl<S: SpadeNum> From<(S, S)> for Point2<S> {
    #[inline]
    fn from(source: (S, S)) -> Self {
        Self::new(source.0, source.1)
    }
}

/// An object with position.
///
/// Vertices need to implement this trait to allow insertion into triangulations.
///
/// Spade defines [Point2] as a basic point primitive.
pub trait HasPosition {
    /// The number type used by this coordinate type.
    type Scalar: SpadeNum;

    /// Returns the position of this object.
    ///
    /// **Note**: It is assumed that the position doesn't change once it has been
    /// inserted into a triangulation. Failing this requirement will lead to logical errors.
    fn position(&self) -> Point2<Self::Scalar>;
}

impl<S: SpadeNum> HasPosition for Point2<S> {
    type Scalar = S;

    fn position(&self) -> Point2<S> {
        *self
    }
}

/*
TODO
#[derive(Debug, Default, PartialEq, PartialOrd, Clone, Copy)]
struct Normal {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Normal {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Normal { x, y, z }
    }

    pub(crate) fn sub(&self, other: Self) -> Self {
        Normal {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    pub(crate) fn add(&self, other: Self) -> Self {
        Normal {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub(crate) fn mul(&self, factor: f64) -> Self {
        Normal {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
        }
    }

    pub(crate) fn length2(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z + self.z
    }

    pub(crate) fn cross(&self, other: Self) -> Self {
        Normal {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Gradient {
    pub x: f64,
    pub y: f64,
}

impl Gradient {
    pub fn new(x: f64, y: f64) -> Self {
        Gradient { x, y }
    }
}

trait HasGradient {
    fn gradient(&self) -> Gradient;
}

trait HasGradientMut {
    fn set_gradient(&mut self, gradient: Gradient);
}
*/
