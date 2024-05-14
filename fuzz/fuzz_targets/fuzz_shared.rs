use spade::{HasPosition, Point2};

#[derive(Clone, Copy, arbitrary::Arbitrary)]
pub struct FuzzPoint {
    pub x: f64,
    pub y: f64,
}

impl HasPosition for FuzzPoint {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        Point2::new(self.x, self.y)
    }
}

impl core::fmt::Debug for FuzzPoint {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_fmt(format_args!("Point2::new({:?}, {:?})", self.x, self.y))
    }
}
