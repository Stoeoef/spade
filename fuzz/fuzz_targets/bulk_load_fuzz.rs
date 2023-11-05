#![no_main]
use libfuzzer_sys::fuzz_target;

use spade::{DelaunayTriangulation, Point2, Triangulation, TriangulationExt};

fuzz_target!(|data: Vec<FuzzPoint>| {
    for p in &data {
        if spade::validate_coordinate(p.x).is_err() || spade::validate_coordinate(p.y).is_err() {
            return;
        }
        if p.x.abs() > 20.0 || p.y.abs() > 20.0 {
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
pub struct FuzzPoint {
    x: f64,
    y: f64,
}

impl core::fmt::Debug for FuzzPoint {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_fmt(format_args!("Point2::new({:?}, {:?})", self.x, self.y))
    }
}
