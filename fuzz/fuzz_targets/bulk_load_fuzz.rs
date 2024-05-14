#![no_main]
mod fuzz_shared;
use fuzz_shared::FuzzPoint;
use libfuzzer_sys::fuzz_target;
use spade::{DelaunayTriangulation, Triangulation, TriangulationExt};

fuzz_target!(|data: Vec<FuzzPoint>| {
    for p in &data {
        if spade::validate_coordinate(p.x).is_err() || spade::validate_coordinate(p.y).is_err() {
            return;
        }
        if p.x.abs() > 20.0 || p.y.abs() > 20.0 {
            return;
        }
    }
    let triangulation = DelaunayTriangulation::<_>::bulk_load(data.clone()).unwrap();
    triangulation.sanity_check();
});
