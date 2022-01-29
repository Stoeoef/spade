use crate::{
    ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation, TriangulationExt,
};

fn refinement_parameters() -> RefinementParameters<f64> {
    Default::default()
}

fn fuzz_test(vertices: Vec<Point2<f64>>) {
    let mut triangulation = ConstrainedDelaunayTriangulation::<_>::bulk_load(vertices).unwrap();
    triangulation.refine(refinement_parameters());
    triangulation.sanity_check();
}

#[test]
fn refine_fuzz_1() {
    fuzz_test(vec![
        Point2::new(5.0, 0.0),
        Point2::new(1.0, 2.0),
        Point2::new(1.0, 16.0),
    ]);
}
