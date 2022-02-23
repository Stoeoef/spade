use crate::{
    ConstrainedDelaunayTriangulation, InsertionError, Point2, RefinementParameters, Triangulation,
    TriangulationExt,
};

fn refinement_parameters() -> RefinementParameters<f64> {
    RefinementParameters::new().with_max_additional_vertices(1000)
}

fn fuzz_test(
    vertices: Vec<Point2<f64>>,
    edges: Vec<[Point2<f64>; 2]>,
) -> Result<(), InsertionError> {
    let mut triangulation = ConstrainedDelaunayTriangulation::<_>::bulk_load(vertices)?;
    let result = triangulation.refine(refinement_parameters());
    for edge in edges {
        triangulation.add_constraint_edge(edge[0], edge[1])?;
    }

    assert!(result.refinement_complete);

    triangulation.sanity_check();
    Ok(())
}

#[test]
fn refine_fuzz_1() -> Result<(), InsertionError> {
    fuzz_test(
        vec![
            Point2::new(5.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(1.0, 16.0),
        ],
        vec![],
    )
}

#[test]
pub fn refine_fuzz_2() -> Result<(), InsertionError> {
    fuzz_test(
        vec![],
        vec![
            [Point2::new(0.0, 0.0), Point2::new(60.0, 0.0)],
            [Point2::new(0.0, 0.0), Point2::new(60.0, 10.0)],
            [Point2::new(0.0, 0.0), Point2::new(40.0, 0.0)],
            [Point2::new(0.0, 0.0), Point2::new(0.0, -40.0)],
        ],
    )
}
