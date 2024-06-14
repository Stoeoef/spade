use anyhow::Result;

use crate::quicksketch::Sketch;
use crate::scenario_list::TryAddConstraintScenario::{BaseCdt, SplitAll};

pub mod quicksketch;
mod scenario;
mod scenario_list;

use scenario_list::*;

/// Used for rendering SVGs for documentation. These are inlined (via #[doc = include_str!(...)])
/// into the doc comment of a few items. That makes sure they will be visible even for offline users.
fn main() -> Result<()> {
    fn save_svg(doc: Sketch, name: &'static str) -> Result<()> {
        doc.save_to_svg(name.to_string(), format!("images/{name}.svg"))?;
        Ok(())
    }

    save_svg(
        try_add_constraint_scenario(BaseCdt)?,
        "add_constraint_and_split_initial",
    )?;

    save_svg(
        try_add_constraint_scenario(SplitAll)?,
        "add_constraint_and_split_added",
    )?;

    save_svg(
        line_intersection_iterator_scenario()?,
        "line_intersection_iterator_scenario",
    )?;

    save_svg(
        natural_neighbor_area_scenario(false)?,
        "natural_neighbor_insertion_cell",
    )?;

    save_svg(
        natural_neighbor_area_scenario(true)?,
        "natural_neighbor_polygon",
    )?;

    save_svg(natural_neighbors_scenario(), "natural_neighbor_scenario")?;

    save_svg(
        refinement_maximum_area_scenario(None),
        "refinement_maximum_area_no_limit",
    )?;

    save_svg(
        refinement_maximum_area_scenario(Some(200.0)),
        "refinement_maximum_area_200",
    )?;

    save_svg(
        refinement_maximum_area_scenario(Some(100.0)),
        "refinement_maximum_area_100",
    )?;

    save_svg(
        shape_iterator_scenario(true, true),
        "shape_iterator_circle_vertices",
    )?;

    save_svg(
        shape_iterator_scenario(true, false),
        "shape_iterator_circle_edges",
    )?;

    save_svg(
        shape_iterator_scenario(false, true),
        "shape_iterator_rectangle_vertices",
    )?;

    save_svg(
        shape_iterator_scenario(false, false),
        "shape_iterator_rectangle_edges",
    )?;

    save_svg(circumcircle_scenario(), "circumcircle")?;
    save_svg(lhs_rhs_scenario(false), "rhs")?;
    save_svg(lhs_rhs_scenario(true), "lhs")?;
    save_svg(outer_face_scenario(), "outer_faces")?;
    save_svg(basic_voronoi_example(), "basic_voronoi")?;
    save_svg(voronoi_edge_details_scenario(), "voronoi_edge_details")?;
    save_svg(
        delaunay_directed_edge_details_scenario(),
        "delaunay_directed_edge_details",
    )?;
    save_svg(
        delaunay_directed_edge_vertex_and_face_scenario(),
        "delaunay_directed_edge_face_and_vertex",
    )?;
    save_svg(cdt_scenario(), "cdt")?;
    save_svg(circular_iterator_example(), "circular_iterator")?;
    save_svg(face_adjacent_edges_scenario(), "face_adjacent_edges")?;
    save_svg(convex_hull_scenario(), "convex_hull")?;
    save_svg(inner_voronoi_vertex_example(), "inner_voronoi_vertex")?;
    save_svg(outer_voronoi_vertex_example(), "outer_voronoi_vertex")?;
    save_svg(dual_edge_example(), "dual_edges")?;
    save_svg(project_point_scenario(), "project_point")?;
    save_svg(refinement_scenario(true), "refined")?;
    save_svg(refinement_scenario(false), "unrefined")?;
    save_svg(angle_limit_scenario(0.0), "angle_limit_00")?;
    save_svg(angle_limit_scenario(20.0), "angle_limit_20")?;
    save_svg(angle_limit_scenario(30.0), "angle_limit_30")?;
    save_svg(angle_limit_scenario(34.0), "angle_limit_34")?;

    save_svg(
        exclude_outer_faces_scenario(false),
        "exclude_outer_faces_unrefined",
    )?;
    save_svg(
        exclude_outer_faces_scenario(true),
        "exclude_outer_faces_refined",
    )?;

    Ok(())
}

pub fn convert_point(point: spade::Point2<f64>) -> quicksketch::Point {
    quicksketch::Point::new(point.x, point.y)
}
