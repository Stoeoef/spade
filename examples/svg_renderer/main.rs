pub mod quicksketch;
mod scenario;
mod scenario_list;

type Result = std::result::Result<(), Box<dyn std::error::Error>>;

/// Used for rendering SVGs for documentation. These are inlined (via #[doc = include_str!(...)])
/// into the doc comment of a few items. That makes sure they will be visible even for offline users.
fn main() -> Result {
    scenario_list::circumcircle_scenario()
        .save_to_svg("circumcircle", "images/circumcircle.svg")?;
    scenario_list::lhs_rhs_scenario(false).save_to_svg("rhs", "images/rhs.svg")?;
    scenario_list::lhs_rhs_scenario(true).save_to_svg("ls", "images/lhs.svg")?;
    scenario_list::outer_face_scenario().save_to_svg("outer_face", "images/outer_faces.svg")?;
    scenario_list::basic_voronoi_example()
        .save_to_svg("basic_voronoi", "images/basic_voronoi.svg")?;
    scenario_list::voronoi_edge_details_scenario()
        .save_to_svg("voronoi_edge_details", "images/voronoi_edge_details.svg")?;
    scenario_list::delaunay_directed_edge_details_scenario().save_to_svg(
        "delaunay_directed_edge_details",
        "images/delaunay_directed_edge_details.svg",
    )?;
    scenario_list::delaunay_directed_edge_vertex_and_face_scenario().save_to_svg(
        "delaunay_directed_edge_face_and_vertex",
        "images/delaunay_directed_edge_face_and_vertex.svg",
    )?;
    scenario_list::cdt_scenario().save_to_svg("cdt_scenario", "images/cdt.svg")?;
    scenario_list::circular_iterator_example()
        .save_to_svg("circular_iterator", "images/circular_iterator.svg")?;
    scenario_list::face_adjacent_edges_scenario()
        .save_to_svg("face_adjacent_edges", "images/face_adjacent_edges.svg")?;
    scenario_list::convex_hull_scenario()
        .save_to_svg("convex_hull", "images/convex_hull_scenario.svg")?;
    Ok(())
}

pub fn convert_point(point: spade::Point2<f64>) -> quicksketch::Point {
    quicksketch::Point::new(point.x, point.y)
}
