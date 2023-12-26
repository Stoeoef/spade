pub mod quicksketch;
mod scenario;
mod scenario_list;
use anyhow::Result;

/// Used for rendering SVGs for documentation. These are inlined (via #[doc = include_str!(...)])
/// into the doc comment of a few items. That makes sure they will be visible even for offline users.
fn main() -> Result<()> {
    scenario_list::natural_neighbor_area_scenario(false)?.save_to_svg(
        "natural_neighbor_insertion_cell",
        "images/natural_neighbor_insertion_cell.svg",
    )?;

    scenario_list::natural_neighbor_area_scenario(true)?.save_to_svg(
        "natural_neighbor_polygon",
        "images/natural_neighbor_polygon.svg",
    )?;

    scenario_list::natural_neighbors_scenario().save_to_svg(
        "natural_neighbor_scenario",
        "images/natural_neighbor_scenario.svg",
    )?;

    scenario_list::refinement_maximum_area_scenario(None).save_to_svg(
        "refinement_maximum_area_no_limit",
        "images/refinement_maximum_area_no_limit.svg",
    )?;

    scenario_list::refinement_maximum_area_scenario(Some(200.0)).save_to_svg(
        "refinement_maximum_area_200",
        "images/refinement_maximum_area_200.svg",
    )?;

    scenario_list::refinement_maximum_area_scenario(Some(100.0)).save_to_svg(
        "refinement_maximum_area_100",
        "images/refinement_maximum_area_100.svg",
    )?;

    scenario_list::shape_iterator_scenario(true, true).save_to_svg(
        "shape_iterator_circle_vertices",
        "images/shape_iterator_circle_vertices.svg",
    )?;
    scenario_list::shape_iterator_scenario(true, false).save_to_svg(
        "shape_iterator_circle_edges",
        "images/shape_iterator_circle_edges.svg",
    )?;
    scenario_list::shape_iterator_scenario(false, true).save_to_svg(
        "shape_iterator_rectangle_vertices",
        "images/shape_iterator_rectangle_vertices.svg",
    )?;
    scenario_list::shape_iterator_scenario(false, false).save_to_svg(
        "shape_iterator_rectangle_edges",
        "images/shape_iterator_rectangle_edges.svg",
    )?;

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
    scenario_list::inner_voronoi_vertex_example()
        .save_to_svg("inner_voronoi_vertex", "images/inner_voronoi_vertex.svg")?;
    scenario_list::outer_voronoi_vertex_example()
        .save_to_svg("outer_voronoi_vertex", "images/outer_voronoi_vertex.svg")?;
    scenario_list::dual_edge_example().save_to_svg("dual_edges", "images/dual_edges.svg")?;
    scenario_list::project_point_scenario()
        .save_to_svg("project_point", "images/project_point.svg")?;
    scenario_list::refinement_scenario(true).save_to_svg("refined", "images/refined.svg")?;
    scenario_list::refinement_scenario(false).save_to_svg("unrefined", "images/unrefined.svg")?;
    scenario_list::angle_limit_scenario(0.0)
        .save_to_svg("angle_limit_0", "images/angle_limit_00.svg")?;
    scenario_list::angle_limit_scenario(20.0)
        .save_to_svg("angle_limit_20", "images/angle_limit_20.svg")?;
    scenario_list::angle_limit_scenario(30.0)
        .save_to_svg("angle_limit_30", "images/angle_limit_30.svg")?;
    scenario_list::angle_limit_scenario(34.0)
        .save_to_svg("angle_limit_34", "images/angle_limit_34.svg")?;

    scenario_list::exclude_outer_faces_scenario(false).save_to_svg(
        "exclude_unrefined",
        "images/exclude_outer_faces_unrefined.svg",
    )?;
    scenario_list::exclude_outer_faces_scenario(true)
        .save_to_svg("exclude_refined", "images/exclude_outer_faces_refined.svg")?;

    Ok(())
}

pub fn convert_point(point: spade::Point2<f64>) -> quicksketch::Point {
    quicksketch::Point::new(point.x, point.y)
}
