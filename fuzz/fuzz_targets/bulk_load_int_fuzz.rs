#![no_main]
use libfuzzer_sys::fuzz_target;

use spade::{DelaunayTriangulation, Point2, Triangulation, TriangulationExt};

fuzz_target!(|data: Vec<IntFuzzPoint>| {
    const MAX_VALUE: i32 = 1024;
    for p in &data {
        if p.x > MAX_VALUE || p.x < -MAX_VALUE || p.y > MAX_VALUE || p.y < -MAX_VALUE {
            return;
        }
    }
    let triangulation = DelaunayTriangulation::<_>::bulk_load(
        data.iter()
            .map(|point| Point2::new(point.x as f64, point.y as f64))
            .collect(),
    )
    .unwrap();

    triangulation.sanity_check();
});

#[derive(Clone, Copy, arbitrary::Arbitrary)]
pub struct IntFuzzPoint {
    x: i32,
    y: i32,
}

impl core::fmt::Debug for IntFuzzPoint {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_fmt(format_args!("Point2::new({:?}.0, {:?}.0)", self.x, self.y))
    }
}
