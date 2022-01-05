#![no_main]
use libfuzzer_sys::fuzz_target;

use spade::triangulation::TriangulationExt;
use spade::{DelaunayTriangulation, Point2, Triangulation};

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

    let mut triangulation: DelaunayTriangulation<_> = DelaunayTriangulation::new();
    for point in data {
        triangulation
            .insert(Point2::new(point.x as f64, point.y as f64))
            .unwrap();
    }

    triangulation.sanity_check();
});

#[derive(Clone, Copy, arbitrary::Arbitrary)]
pub struct IntFuzzPoint {
    x: i32,
    y: i32,
}

impl std::fmt::Debug for IntFuzzPoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("Point2::new({:?}.0, {:?}.0)", self.x, self.y))
    }
}
