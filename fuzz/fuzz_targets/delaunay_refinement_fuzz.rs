#![no_main]
use libfuzzer_sys::fuzz_target;

use spade::{
    AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation,
    TriangulationExt,
};

fuzz_target!(|data: Vec<IntFuzzPoint>| {
    const MAX_VALUE: i32 = 140;
    for p in &data {
        if p.x > MAX_VALUE || p.x < -MAX_VALUE || p.y > MAX_VALUE || p.y < -MAX_VALUE {
            return;
        }
    }
    let mut triangulation = ConstrainedDelaunayTriangulation::<_>::bulk_load(
        data.iter()
            .map(|point| Point2::new(point.x as f64, point.y as f64))
            .collect(),
    )
    .unwrap();

    let parameters = RefinementParameters::new()
        .with_angle_limit(AngleLimit::new_from_radius_to_shortest_edge_ratio(0.7));

    triangulation.refine(parameters);
    triangulation.sanity_check();
});

#[derive(Clone, Copy, arbitrary::Arbitrary)]
pub struct IntFuzzPoint {
    x: i32,
    y: i32,
}

impl std::fmt::Debug for IntFuzzPoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("Point2::new({:?}.0, {:?}.0),", self.x, self.y))
    }
}
