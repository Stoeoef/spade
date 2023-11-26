use super::quicksketch::{
    ArrowType, HorizontalAlignment, Point, Sketch, SketchColor, SketchElement, SketchFill,
    SketchLayer, StrokeStyle, Vector,
};

use anyhow::{Context, Result};

use cgmath::{Angle, Bounded, Deg, EuclideanSpace, InnerSpace, Point2, Vector2};

use spade::{
    handles::{
        FixedDirectedEdgeHandle,
        VoronoiVertex::{self, Inner, Outer},
    },
    AngleLimit, FloatTriangulation as _, HasPosition, InsertionError, RefinementParameters,
    Triangulation as _,
};

use crate::{
    convert_point,
    scenario::{
        convert_triangulation, Cdt, ConversionOptions, DirectedEdgeType, EdgeMode, FaceType,
        Triangulation, UndirectedEdgeType, VertexType,
    },
};

fn big_triangulation() -> Result<Triangulation, InsertionError> {
    let mut result = Triangulation::new();

    result.insert(VertexType::new(15.0, 39.0))?;
    result.insert(VertexType::new(-5.0, -42.0))?;
    result.insert(VertexType::new(30.0, 2.0))?;
    result.insert(VertexType::new(5.0, -12.0))?;
    result.insert(VertexType::new(-35.0, 32.0))?;
    result.insert(VertexType::new(-2.0, -24.0))?;
    result.insert(VertexType::new(-67.0, 52.0))?;
    result.insert(VertexType::new(-14.0, -52.0))?;
    result.insert(VertexType::new(76.0, 10.0))?;
    result.insert(VertexType::new(12.0, 10.0))?;
    result.insert(VertexType::new(-35.0, 5.0))?;
    result.insert(VertexType::new(33.0, -30.0))?;
    result.insert(VertexType::new(-30.0, -25.0))?;
    result.insert(VertexType::new(-55.0, -25.0))?;
    result.insert(VertexType::new(50.0, 50.0))?;
    result.insert(VertexType::new(45.0, -47.0))?;

    Ok(result)
}

fn add_circular_arrows(sketch: &mut Sketch, center: Point, is_lhs: bool, radius: f64) {
    let angle_to_position =
        |angle: Deg<f64>| center + Vector::new(angle.sin(), -angle.cos()) * radius;

    const WIDTH: f64 = 1.0;
    const COLOR: SketchColor = SketchColor::ROYAL_BLUE;
    const ARROW_TYPE: ArrowType = ArrowType::FilledArrow;

    let zero = Deg(0.0);

    let radii = Vector::new(radius, radius);
    let arrow = SketchElement::path()
        .stroke_width(WIDTH)
        .stroke_color(COLOR);

    let arrow1 = arrow.clone().move_to(angle_to_position(zero)).arc_to(
        radii,
        zero,
        false,
        true,
        angle_to_position(Deg(140.0)),
    );
    let arrow1 = if is_lhs {
        arrow1.with_arrow_end(ARROW_TYPE)
    } else {
        arrow1.with_arrow_start(ARROW_TYPE)
    };

    let arrow2 = arrow.move_to(angle_to_position(Deg(180.0))).arc_to(
        radii,
        zero,
        false,
        true,
        angle_to_position(Deg(320.0)),
    );

    let arrow2 = if is_lhs {
        arrow2.with_arrow_end(ARROW_TYPE)
    } else {
        arrow2.with_arrow_start(ARROW_TYPE)
    };

    sketch.add(arrow1);
    sketch.add(arrow2);
}

pub fn lhs_rhs_scenario(is_lhs: bool) -> Sketch {
    let mut triangulation = Triangulation::new();

    triangulation.insert(VertexType::new(-67.0, 52.0)).unwrap();
    triangulation.insert(VertexType::new(-40.0, -24.0)).unwrap();
    triangulation.insert(VertexType::new(-30.0, 20.0)).unwrap();
    triangulation.insert(VertexType::new(-6.0, -42.0)).unwrap();
    triangulation.insert(VertexType::new(-5.0, 0.0)).unwrap();
    triangulation.insert(VertexType::new(15.0, 42.0)).unwrap();
    triangulation.insert(VertexType::new(30.0, 2.0)).unwrap();

    let mut sketch = convert_triangulation(
        &triangulation,
        &ConversionOptions {
            directed_edge_mode: EdgeMode::Directed { reversed: is_lhs },
            ..Default::default()
        },
    );

    for face in triangulation.inner_faces() {
        let circumcenter = face.center();
        add_circular_arrows(&mut sketch, convert_point(circumcenter), is_lhs, 4.5);
    }

    sketch.set_width(420);

    sketch
}

pub fn circumcircle_scenario() -> Sketch {
    let triangulation = big_triangulation().unwrap();

    let mut sketch = convert_triangulation(&triangulation, &Default::default());

    for face in triangulation.inner_faces().skip(5).take(8) {
        let (circumcenter, radius_squared) = face.circumcircle();
        if radius_squared < 880.0 {
            sketch.add(
                SketchElement::circle(convert_point(circumcenter), radius_squared.sqrt())
                    .stroke_width(0.5)
                    .stroke_color(SketchColor::ROYAL_BLUE),
            );
        }
    }

    sketch.set_width(500);
    sketch
}

pub fn outer_face_scenario() -> Sketch {
    let mut triangulation = Triangulation::default();

    triangulation.insert(VertexType::new(50.0, 10.0)).unwrap();
    triangulation.insert(VertexType::new(40.0, -35.0)).unwrap();
    triangulation.insert(VertexType::new(-50.0, 50.0)).unwrap();
    triangulation.insert(VertexType::new(-20.0, 0.0)).unwrap();
    triangulation.insert(VertexType::new(35.0, 45.0)).unwrap();

    let mut sketch = convert_triangulation(&triangulation, &Default::default());

    for face in triangulation.inner_faces() {
        let center = face.center();
        const FONT_SIZE: f64 = 10.0;

        sketch.add(
            SketchElement::text("inner")
                .position(convert_point(center))
                .horizontal_alignment(HorizontalAlignment::Middle)
                .font_size(FONT_SIZE),
        );

        sketch.add(
            SketchElement::text("outer")
                .position(Point::new(-40.0, -25.0))
                .font_size(FONT_SIZE),
        );
    }

    sketch.set_width(400);
    sketch.set_relative_padding(0.05);

    sketch
}

pub fn basic_voronoi_example() -> Sketch {
    let triangulation = big_triangulation().unwrap();

    let mut sketch = convert_triangulation(&triangulation, &Default::default());
    const LINE_COLOR: SketchColor = SketchColor::ROYAL_BLUE;

    for edge in triangulation.undirected_voronoi_edges() {
        match edge.vertices() {
            [Inner(from), Inner(to)] => {
                sketch.add(
                    SketchElement::line(
                        convert_point(from.circumcenter()),
                        convert_point(to.circumcenter()),
                    )
                    .stroke_color(LINE_COLOR),
                );
            }
            [Inner(from), Outer(edge)] | [Outer(edge), Inner(from)] => {
                let from = convert_point(from.circumcenter());
                let to_direction = edge.direction_vector();
                let to_direction = Vector::new(to_direction.x, to_direction.y);

                sketch.add(
                    SketchElement::line(from, from + to_direction * 4.0)
                        .stroke_color(LINE_COLOR)
                        .stroke_style(StrokeStyle::Dashed),
                );
            }
            [Outer(_), Outer(_)] => {}
        }
    }

    let mut min: Point = Bounded::max_value();
    let mut max: Point = Bounded::min_value();

    for vertex in triangulation.vertices() {
        let position = vertex.position();
        min = min.zip(convert_point(position), f64::min);
        max = max.zip(convert_point(position), f64::max);
    }

    min.y -= 5.0;
    max.y -= 15.0;

    sketch
        .set_width(700)
        .set_relative_padding(0.24)
        .set_view_box_min(min)
        .set_view_box_max(max);
    sketch
}

pub fn voronoi_edge_details_scenario() -> Sketch {
    let mut triangulation = Triangulation::new();

    triangulation.insert(VertexType::new(-50.0, -50.0)).unwrap();
    triangulation.insert(VertexType::new(-41.0, -60.0)).unwrap();
    triangulation.insert(VertexType::new(-41.0, 61.0)).unwrap();
    triangulation.insert(VertexType::new(-40.0, -40.0)).unwrap();
    triangulation.insert(VertexType::new(-30.0, 25.0)).unwrap();
    let vertex = triangulation.insert(VertexType::new(30.0, 0.0)).unwrap();
    triangulation.insert(VertexType::new(40.0, 20.0)).unwrap();
    triangulation.insert(VertexType::new(41.0, -70.0)).unwrap();
    triangulation.insert(VertexType::new(41.0, -80.0)).unwrap();
    triangulation.insert(VertexType::new(41.0, -20.0)).unwrap();
    triangulation.insert(VertexType::new(107.0, -20.0)).unwrap();
    triangulation.insert(VertexType::new(90.0, -40.0)).unwrap();
    triangulation.insert(VertexType::new(90.0, 60.0)).unwrap();

    for undirected_edge in triangulation.fixed_undirected_edges() {
        triangulation
            .undirected_edge_data_mut(undirected_edge)
            .color = SketchColor::DARK_GRAY;
    }

    let mut sketch = convert_triangulation(
        &triangulation,
        &ConversionOptions {
            vertex_stroke_color: SketchColor::DARK_GRAY,
            vertex_color: SketchColor::GRAY,
            ..Default::default()
        },
    );

    const LINE_COLOR: SketchColor = SketchColor::ROYAL_BLUE;

    for face in triangulation.inner_faces() {
        let center = convert_point(face.circumcenter());
        sketch.add(
            SketchElement::circle(center, 2.0)
                .fill(SketchFill::solid(SketchColor::ROYAL_BLUE))
                .stroke_color(SketchColor::BLACK)
                .stroke_width(0.5),
        );
    }

    let example_face = triangulation.vertex(vertex).as_voronoi_face();

    const SHIFT: f64 = -2.5;

    let create_line =
        |from: VoronoiVertex<VertexType, DirectedEdgeType, UndirectedEdgeType, FaceType>,
         to: VoronoiVertex<VertexType, DirectedEdgeType, UndirectedEdgeType, FaceType>| {
            let from = from.as_delaunay_face().unwrap().circumcenter();
            let to = to.as_delaunay_face().unwrap().circumcenter();

            SketchElement::line(convert_point(from), convert_point(to))
                .stroke_color(LINE_COLOR)
                .draw_double_line()
                .with_arrow_start(ArrowType::HalfArrow)
                .shift_from_and_to(SHIFT)
        };

    for edge in triangulation.undirected_voronoi_edges() {
        if let [Inner(from), Inner(to)] = edge.vertices() {
            let directed = edge.as_directed();
            if directed.face() == example_face || directed.rev().face() == example_face {
                // Edges of the example face are drawn manually
                continue;
            }

            sketch.add(create_line(Inner(from), Inner(to)));
        }
    }

    let edge = example_face.adjacent_edges().next().unwrap().next().next();

    let main_edge = create_line(edge.from(), edge.to());

    // main edge and rev
    sketch.add(
        main_edge
            .create_adjacent_text("edge")
            .font_size(5.0)
            .dy(-3.3),
    );

    let rev_text = create_line(edge.to(), edge.from())
        .create_adjacent_text("edge.rev()")
        .font_size(5.0)
        .dy(5.7);
    sketch.add(rev_text);
    sketch.add(main_edge);

    // next
    let next = edge.prev(); // Take prev instead of next since SVG renders everything with a LHS
    let next_edge = create_line(next.from(), next.to());

    sketch.add(
        next_edge
            .create_adjacent_text("edge.next()")
            .font_size(5.0)
            .dy(5.3),
    );
    sketch.add(next_edge);

    // auxiliary edge
    let aux = next.prev();
    sketch.add(create_line(aux.from(), aux.to()));

    // prev
    let prev = aux.prev();
    let prev_line = create_line(prev.from(), prev.to());
    sketch.add(
        prev_line
            .create_adjacent_text("edge.prev()")
            .font_size(5.0)
            .dy(-3.2),
    );
    sketch.add(prev_line);

    sketch.set_view_box_min(Point::new(-11.0, -30.0));
    sketch.set_view_box_max(Point::new(56.0, 27.0));
    sketch.set_width(400);

    sketch
}

fn delaunay_edge_details_triangulation(
) -> Result<(Triangulation, FixedDirectedEdgeHandle), InsertionError> {
    let mut triangulation = Triangulation::new();

    triangulation.insert(VertexType::new(-50.0, -50.0))?;
    triangulation.insert(VertexType::new(-41.0, -60.0))?;
    triangulation.insert(VertexType::new(-41.0, 61.0))?;
    triangulation.insert(VertexType::new(-40.0, -40.0))?;
    triangulation.insert(VertexType::new(-30.0, 25.0))?;
    let to = triangulation.insert(VertexType::new(0.0, 0.0))?;
    let from = triangulation.insert(VertexType::new(40.0, 20.0))?;
    triangulation.insert(VertexType::new(41.0, -70.0))?;
    triangulation.insert(VertexType::new(41.0, -80.0))?;
    triangulation.insert(VertexType::new(41.0, -20.0))?;
    triangulation.insert(VertexType::new(107.0, -20.0))?;
    triangulation.insert(VertexType::new(90.0, -40.0))?;
    triangulation.insert(VertexType::new(90.0, 60.0))?;

    let edge = triangulation
        .get_edge_from_neighbors(from, to)
        .unwrap()
        .fix();
    Ok((triangulation, edge))
}

pub fn delaunay_directed_edge_details_scenario() -> Sketch {
    let (triangulation, edge) = delaunay_edge_details_triangulation().unwrap();
    let edge = triangulation.directed_edge(edge);

    let mut sketch = convert_triangulation(
        &triangulation,
        &ConversionOptions {
            directed_edge_mode: EdgeMode::Directed { reversed: false },
            ..Default::default()
        },
    );

    let from_pos = convert_point(edge.from().position());
    let to_pos = convert_point(edge.to().position());

    const FONT_SIZE: f64 = 5.0;
    const DY: f64 = -2.1;
    let edge_label = SketchElement::line(from_pos, to_pos)
        .create_adjacent_text("e")
        .font_size(FONT_SIZE)
        .dy(DY);

    sketch.add(edge_label);
    let rev_label = SketchElement::line(from_pos, to_pos)
        .create_adjacent_text("e.rev()")
        .font_size(FONT_SIZE)
        .dy(5.0);

    sketch.add(rev_label);

    let next = edge.prev(); // Use prev since SVG uses a left handed coordinate system
    let next_label = SketchElement::line(from_pos, convert_point(next.from().position()))
        .create_adjacent_text("e.next()")
        .font_size(FONT_SIZE)
        .dy(DY);

    sketch.add(next_label);

    let prev = edge.next(); // Use prev since SVG uses a left handed coordinate system
    let prev_label = SketchElement::line(convert_point(prev.to().position()), to_pos)
        .create_adjacent_text("e.prev()")
        .font_size(FONT_SIZE)
        .dy(5.2);

    sketch.add(prev_label);

    sketch.set_view_box_min(Point::new(-5.0, -20.0));
    sketch.set_view_box_max(Point::new(46.0, 17.0));
    sketch.set_width(300);

    sketch
}

pub fn delaunay_directed_edge_vertex_and_face_scenario() -> Sketch {
    let (mut triangulation, edge) = delaunay_edge_details_triangulation().unwrap();

    let face_fixed = triangulation.directed_edge(edge).face().fix();
    triangulation.face_data_mut(face_fixed).fill =
        SketchFill::StripePattern(SketchColor::LIGHT_GRAY, SketchColor::DARK_GRAY, 5.0);

    let mut sketch = convert_triangulation(
        &triangulation,
        &ConversionOptions {
            directed_edge_mode: EdgeMode::Directed { reversed: false },
            ..Default::default()
        },
    );

    let edge = triangulation.directed_edge(edge);
    let face = edge.face().as_inner().unwrap();

    let from_pos = convert_point(edge.from().position());
    let to_pos = convert_point(edge.to().position());

    const FONT_SIZE: f64 = 5.0;
    const DY: f64 = -2.1;
    let edge_label = SketchElement::line(from_pos, to_pos)
        .create_adjacent_text("e")
        .font_size(FONT_SIZE)
        .dy(DY);

    sketch.add(edge_label);

    let face_position = convert_point(face.circumcenter());

    sketch.add(
        SketchElement::text("e.face()")
            .position(face_position)
            .horizontal_alignment(HorizontalAlignment::Middle)
            .font_size(FONT_SIZE),
    );

    let edge_from_position = convert_point(edge.to().position());
    sketch.add(
        SketchElement::text("e.from()")
            .position(edge_from_position + Vector::new(-7.0, 16.0))
            .font_size(FONT_SIZE),
    );

    sketch.add(
        SketchElement::path()
            .move_to(edge_from_position + Vector::new(1.8, 6.0))
            .line_to(edge_from_position + Vector::new(4.0, 12.0))
            .stroke_color(SketchColor::BLACK)
            .with_arrow_start(ArrowType::FilledArrow)
            .stroke_style(StrokeStyle::SmallDashed)
            .stroke_width(0.7),
    );

    let edge_to_position = convert_point(edge.from().position());
    sketch.add(
        SketchElement::text("e.to()")
            .position(edge_to_position + Vector::new(-7.0, 16.0))
            .font_size(FONT_SIZE),
    );

    sketch.add(
        SketchElement::path()
            .move_to(edge_to_position + Vector::new(0.5, 6.0))
            .line_to(edge_to_position + Vector::new(0.5, 11.0))
            .stroke_color(SketchColor::BLACK)
            .with_arrow_start(ArrowType::FilledArrow)
            .stroke_style(StrokeStyle::SmallDashed)
            .stroke_width(0.7),
    );

    sketch.set_view_box_min(Point::new(-5.0, -20.0));
    sketch.set_view_box_max(Point::new(46.0, 35.0));
    sketch.set_width(300);

    sketch
}

pub fn cdt_scenario() -> Sketch {
    let create_cdt = |offset| {
        let mut cdt = Cdt::new();
        let v = |x, y| VertexType::new(x + offset, y);

        cdt.insert(v(-50.0, -50.0)).unwrap();
        cdt.insert(v(-41.0, -60.0)).unwrap();
        cdt.insert(v(-41.0, 61.0)).unwrap();
        cdt.insert(v(-20.0, -40.0)).unwrap();
        cdt.insert(v(-10.0, 15.0)).unwrap();
        cdt.insert(v(0.0, 0.0)).unwrap();
        cdt.insert(v(40.0, 20.0)).unwrap();
        cdt.insert(v(41.0, -60.0)).unwrap();
        cdt.insert(v(41.0, -70.0)).unwrap();
        cdt.insert(v(41.0, -20.0)).unwrap();
        cdt.insert(v(65.0, -20.0)).unwrap();
        cdt.insert(v(70.0, -40.0)).unwrap();
        cdt.insert(v(70.0, 60.0)).unwrap();
        cdt.insert(v(-20.0, 40.0)).unwrap();
        cdt.insert(v(35.0, 35.0)).unwrap();
        cdt.insert(v(35.0, -35.0)).unwrap();
        cdt.insert(v(-35.0, -35.0)).unwrap();
        cdt.insert(v(-35.0, 35.0)).unwrap();
        cdt
    };

    let mut cdt = create_cdt(0.0);
    let v0 = cdt.insert(VertexType::new(35.0, 35.0)).unwrap();
    let v1 = cdt.insert(VertexType::new(35.0, -35.0)).unwrap();
    let v2 = cdt.insert(VertexType::new(-35.0, -35.0)).unwrap();
    let v3 = cdt.insert(VertexType::new(-35.0, 35.0)).unwrap();

    cdt.add_constraint(v0, v1);
    cdt.add_constraint(v1, v2);
    cdt.add_constraint(v2, v3);
    cdt.add_constraint(v3, v0);

    for edge in cdt.fixed_undirected_edges() {
        if cdt.is_constraint_edge(edge) {
            cdt.undirected_edge_data_mut(edge).data_mut().color = SketchColor::DARK_RED;
        }
    }

    let mut sketch = convert_triangulation(&cdt, &Default::default());

    let cdt2 = create_cdt(140.0);
    let sketch2 = convert_triangulation(&cdt2, &Default::default());

    for layer in [
        SketchLayer::EDGES,
        SketchLayer::VERTICES,
        SketchLayer::BACKGROUND,
    ] {
        let layer = layer as usize;
        sketch.items[layer].extend_from_slice(&sketch2.items[layer]);
    }

    sketch.add(
        SketchElement::text("CDT")
            .position(Point::new(-30.0, -67.0))
            .horizontal_alignment(HorizontalAlignment::Middle)
            .font_size(13.0),
    );

    sketch.add(
        SketchElement::text("No CDT")
            .position(Point::new(120.0, -67.0))
            .horizontal_alignment(HorizontalAlignment::Middle)
            .font_size(13.0),
    );

    sketch.set_width(935);
    sketch.set_relative_padding(0.03);
    sketch
}

pub fn circular_iterator_example() -> Sketch {
    let mut triangulation = Triangulation::new();

    let mut v0 = VertexType::new(0.0, 0.0);
    v0.radius = 3.5;
    let v0 = triangulation.insert(v0).unwrap();
    triangulation.insert(VertexType::new(66.0, -5.0)).unwrap();
    triangulation.insert(VertexType::new(6.0, 40.0)).unwrap();
    triangulation.insert(VertexType::new(-55.0, 5.0)).unwrap();
    triangulation.insert(VertexType::new(60.0, 40.0)).unwrap();
    triangulation.insert(VertexType::new(-45.0, 25.0)).unwrap();
    triangulation.insert(VertexType::new(45.0, -40.0)).unwrap();
    triangulation.insert(VertexType::new(-49.0, -30.0)).unwrap();

    let mut sketch = convert_triangulation(
        &triangulation,
        &ConversionOptions {
            directed_edge_mode: EdgeMode::Directed { reversed: false },
            ..Default::default()
        },
    );

    for (index, edge) in triangulation.vertex(v0).out_edges().enumerate() {
        let from = convert_point(edge.from().position());
        let to = convert_point(edge.to().position());
        let dy = match index {
            0 | 2 | 1 | 6 => -2.0,
            _ => 7.0,
        };

        sketch.add(
            SketchElement::line(from, to)
                .create_adjacent_text(format!("e{}", 6 - index))
                .dy(dy)
                .horizontal_alignment(HorizontalAlignment::Middle)
                .font_size(8.0),
        );
    }

    sketch.add(
        SketchElement::text("v")
            .horizontal_alignment(HorizontalAlignment::Middle)
            .dy(2.1)
            .font_size(7.0),
    );

    add_circular_arrows(&mut sketch, Point::new(0.0, 0.0), false, 10.0);

    sketch.set_relative_padding(-0.01);
    sketch.set_width(400);

    sketch
}

pub fn face_adjacent_edges_scenario() -> Sketch {
    let mut d = Triangulation::new();

    d.insert(VertexType::new(0.0, -15.0)).unwrap();
    d.insert(VertexType::new(60.0, 15.0)).unwrap();
    d.insert(VertexType::new(30.0, 45.0)).unwrap();

    let conversion_options = ConversionOptions {
        directed_edge_mode: EdgeMode::Disabled,
        ..Default::default()
    };
    let mut sketch = convert_triangulation(&d, &conversion_options);

    let face = d.inner_faces().next().unwrap();

    for (index, edge) in face.adjacent_edges().into_iter().enumerate() {
        let from = convert_point(edge.from().position());
        let to = convert_point(edge.to().position());

        let line = SketchElement::line(from, to)
            .stroke_color(SketchColor::ROYAL_BLUE)
            .with_arrow_start(ArrowType::FilledArrow)
            .shift_from(-7.0)
            .shift_to(-2.4);
        sketch.add(
            line.create_adjacent_text(format!("e{}", 2 - index))
                .font_size(6.0)
                .dy(if index == 1 { 5.0 } else { -1.2 }),
        );
        sketch.add(line);
    }

    let mut center = convert_point(face.center());
    center.x += 4.0;
    sketch.add(
        SketchElement::text("face")
            .horizontal_alignment(HorizontalAlignment::Middle)
            .position(center)
            .font_size(6.0)
            .dy(4.0),
    );

    sketch.set_width(220);

    sketch
}

pub fn convex_hull_scenario() -> Sketch {
    let d = big_triangulation().unwrap();

    let mut sketch = convert_triangulation(&d, &Default::default());

    for (index, edge) in d.convex_hull().enumerate() {
        let from = convert_point(edge.from().position());
        let to = convert_point(edge.to().position());

        let line = SketchElement::line(from, to)
            .stroke_color(SketchColor::ROYAL_BLUE)
            .with_arrow_start(ArrowType::FilledArrow)
            .stroke_width(1.4)
            .shift_from(-9.0)
            .shift_to(-2.4);
        sketch.add(
            line.create_adjacent_text(format!("e{}", d.convex_hull_size() - index - 1))
                .font_size(7.0)
                .dy(if index == 4 || index == 3 { 6.0 } else { -1.7 }),
        );
        sketch.add(line);
    }

    sketch.set_width(470);

    sketch
}

pub fn inner_voronoi_vertex_example() -> Sketch {
    let mut triangulation = big_triangulation().unwrap();

    for face in triangulation.fixed_inner_faces().skip(13).take(6) {
        triangulation.face_data_mut(face).fill = SketchFill::Solid(SketchColor::STEEL_BLUE);
    }

    let mut sketch = convert_triangulation(&triangulation, &ConversionOptions::default());

    for face in triangulation.inner_faces().skip(13).take(6) {
        let circumcenter = convert_point(face.circumcenter());
        let center = convert_point(face.center());

        sketch.add(
            SketchElement::line(circumcenter, center)
                .with_arrow_end(ArrowType::FilledArrow)
                .stroke_color(SketchColor::SALMON)
                .stroke_width(0.8)
                .shift_to(-2.0),
        );

        sketch.add(
            SketchElement::circle(circumcenter, 2.0)
                .stroke_color(SketchColor::BLACK)
                .stroke_width(0.8)
                .fill(SketchFill::Solid(SketchColor::ROYAL_BLUE)),
        );
    }
    sketch
}

pub fn outer_voronoi_vertex_example() -> Sketch {
    let mut triangulation = big_triangulation().unwrap();

    let edge_selection = triangulation
        .convex_hull()
        .step_by(2)
        .map(|e| e.fix())
        .collect::<Vec<_>>();

    for edge in &edge_selection {
        triangulation
            .undirected_edge_data_mut(edge.as_undirected())
            .color = SketchColor::ROYAL_BLUE
    }

    let mut sketch = convert_triangulation(&triangulation, &ConversionOptions::default());

    for edge in edge_selection {
        let edge = triangulation.directed_edge(edge);
        let from = convert_point(edge.from().position());
        let to = convert_point(edge.to().position());

        let center = EuclideanSpace::centroid(&[from, to]);

        let diff = to - from;
        let direction = Vector2::new(-diff.y, diff.x).normalize();

        let arrow_end = center + direction * 30.0;

        sketch.add(
            SketchElement::line(center, arrow_end)
                .stroke_color(SketchColor::SALMON)
                .stroke_width(0.8)
                .stroke_style(StrokeStyle::Dashed)
                .with_arrow_end(ArrowType::FilledArrow),
        );

        sketch.add(
            SketchElement::circle(center, 1.0)
                .stroke_width(0.7)
                .stroke_color(SketchColor::BLACK)
                .fill(SketchFill::Solid(SketchColor::LIGHT_BLUE)),
        );
    }

    sketch
}

pub fn dual_edge_example() -> Sketch {
    let mut triangulation = big_triangulation().unwrap();

    for edge in triangulation.fixed_undirected_edges().skip(2).take(4) {
        triangulation.undirected_edge_data_mut(edge).color = SketchColor::SALMON;
    }

    let mut sketch = convert_triangulation(&triangulation, &ConversionOptions::default());

    for edge in triangulation.undirected_edges().skip(2).take(4) {
        let edge = edge.as_voronoi_edge().as_directed();

        if let (Inner(from), Inner(to)) = (edge.from(), edge.to()) {
            let from = convert_point(from.circumcenter());
            let to = convert_point(to.circumcenter());
            sketch.add(
                SketchElement::line(from, to)
                    .stroke_color(SketchColor::ROYAL_BLUE)
                    .stroke_width(1.3),
            );
            let circle = |pos| {
                SketchElement::circle(pos, 1.3)
                    .fill(SketchFill::Solid(SketchColor::LIGHT_BLUE))
                    .stroke_color(SketchColor::BLACK)
                    .stroke_width(1.0)
            };

            sketch.add(circle(from));
            sketch.add(circle(to));
        }
    }

    sketch.set_width(500);

    sketch
}

pub fn project_point_scenario() -> Sketch {
    let mut sketch = Sketch::new();

    let from = Point2::new(0.0, 0.0);
    let to = Point2::new(20.0, 4.0);
    let offset = Vector2::new(to.y, -to.x) * 2.0;
    let direction = to - from;

    let far_before = from - direction * 2.0;
    let far_behind = to + direction * 2.0;

    let before_area = SketchElement::path()
        .move_to(far_before + offset)
        .line_to(from + offset)
        .line_to(from - offset)
        .line_to(far_before - offset);

    let edge_area = SketchElement::path()
        .move_to(from + offset)
        .line_to(to + offset)
        .line_to(to - offset)
        .line_to(from - offset);

    let behind_area = SketchElement::path()
        .move_to(far_behind + offset)
        .line_to(to + offset)
        .line_to(to - offset)
        .line_to(far_behind - offset);

    sketch.add(before_area.fill(SketchFill::solid(SketchColor::DARK_SEA_GREEN)));
    sketch.add(edge_area.fill(SketchFill::solid(SketchColor::CORNFLOWER_BLUE)));
    sketch.add(behind_area.fill(SketchFill::solid(SketchColor::SALMON)));

    sketch.add(
        SketchElement::line(from, to - direction * 0.12)
            .with_arrow_end(ArrowType::FilledArrow)
            .stroke_width(0.5),
    );

    const FONT_SIZE: f64 = 2.0;

    sketch.add(
        SketchElement::text("before")
            .position(Point2::new(-7.0, -2.0))
            .font_size(FONT_SIZE),
    );
    sketch.add(
        SketchElement::text("on edge")
            .position(Point2::new(4.0, 7.0))
            .font_size(FONT_SIZE),
    );
    sketch.add(
        SketchElement::text("behind")
            .position(Point2::new(21.0, 5.0))
            .font_size(FONT_SIZE),
    );

    sketch.add(
        SketchElement::line(from + offset, from - offset)
            .stroke_width(0.1)
            .stroke_style(StrokeStyle::SmallDashed)
            .stroke_color(SketchColor::DARK_SLATE_GREY),
    );

    sketch.add(
        SketchElement::line(to + offset, to - offset)
            .stroke_width(0.1)
            .stroke_style(StrokeStyle::SmallDashed)
            .stroke_color(SketchColor::DARK_SLATE_GREY),
    );

    sketch
        .set_view_box_min(from)
        .set_view_box_max(to)
        .set_relative_padding(0.42)
        .set_height(400);
    sketch
}

pub fn shape_iterator_scenario(use_circle_metric: bool, iterate_vertices: bool) -> Sketch {
    let t = big_triangulation().unwrap();

    let center = spade::Point2::new(0.0, 2.0);
    let radius = 30.0;
    let radius_2 = radius * radius;

    let lower = spade::Point2::new(-35.0, -30.0);
    let upper = spade::Point2::new(20.0, 55.0);

    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    match (use_circle_metric, iterate_vertices) {
        (true, true) => vertices = t.get_vertices_in_circle(center, radius_2).collect(),
        (true, false) => edges = t.get_edges_in_circle(center, radius_2).collect(),
        (false, true) => vertices = t.get_vertices_in_rectangle(lower, upper).collect(),
        (false, false) => edges = t.get_edges_in_rectangle(lower, upper).collect(),
    }

    let mut result = convert_triangulation(&t, &Default::default());

    if use_circle_metric {
        result.add(
            SketchElement::circle(convert_point(center), radius)
                .stroke_width(0.5)
                .stroke_color(SketchColor::TEAL),
        );
    } else {
        let l = convert_point(lower);
        let u = convert_point(upper);

        result.add_with_layer(
            SketchElement::path()
                .move_to(l)
                .line_to(Point2::new(l.x, u.y))
                .line_to(u)
                .line_to(Point2::new(u.x, l.y))
                .close()
                .stroke_width(0.7)
                .stroke_color(SketchColor::TEAL),
            SketchLayer::BACKGROUND,
        );
    }

    for vertex in vertices {
        result.add_with_layer(
            SketchElement::circle(convert_point(vertex.position()), 2.0).fill(SketchColor::SALMON),
            SketchLayer::VERTICES,
        );
    }

    for edge in edges {
        let [from, to] = edge.positions().map(convert_point);
        result.add_with_layer(
            SketchElement::line(from, to)
                .stroke_color(SketchColor::SALMON)
                .stroke_width(0.71),
            SketchLayer::EDGES,
        );
    }

    result.set_width(500);
    result
}

pub fn natural_neighbors_scenario() -> Sketch {
    let triangulation = big_triangulation().unwrap();

    let nn = triangulation.natural_neighbor();

    let mut nns = Vec::new();
    let query_point = spade::Point2::new(-1.0, 4.0);
    nn.get_weights(query_point, &mut nns);

    let mut result = convert_triangulation(&triangulation, &Default::default());

    let offsets = [
        Vector2::new(2.0, 5.0),
        Vector2::new(-4.0, 6.0),
        Vector2::new(0.0, 6.0),
        Vector2::new(1.0, 6.0),
        Vector2::new(0.0, 8.0),
        Vector2::new(3.0, 4.0),
        Vector2::new(5.0, 2.0),
    ];

    for (index, (neighbor, weight)) in nns.iter().enumerate() {
        let neighbor = triangulation.vertex(*neighbor);
        let position = convert_point(neighbor.position());
        result.add(
            SketchElement::circle(position, 2.0)
                .fill(SketchFill::Solid(SketchColor::SALMON))
                .stroke_width(0.5)
                .stroke_color(SketchColor::BLACK),
        );
        result.add(
            SketchElement::text(index.to_string())
                .position(position)
                .font_size(3.0)
                .horizontal_alignment(HorizontalAlignment::Middle)
                .dy(1.15),
        );

        result.add(
            SketchElement::text(format!("{:.2}", weight))
                .position(position + offsets[index])
                .font_size(5.0),
        );
    }

    result.add(
        SketchElement::circle(convert_point(query_point), 1.5)
            .fill(SketchFill::Solid(SketchColor::ROYAL_BLUE))
            .stroke_width(0.2)
            .stroke_color(SketchColor::BLACK),
    );

    result.set_view_box_min(Point2::new(-40.0, -30.0));
    result.set_view_box_max(Point2::new(30.0, 45.0));
    result.set_width(400);
    result
}

/// Only used for internal documentation of natural neighbor area calculation
pub fn natural_neighbor_area_scenario(include_faces: bool) -> Result<Sketch> {
    let mut triangulation: Triangulation = Default::default();

    triangulation.insert(VertexType::new(45.0, 30.0))?;
    triangulation.insert(VertexType::new(7.5, 40.0))?;
    triangulation.insert(VertexType::new(-45.0, 42.0))?;
    triangulation.insert(VertexType::new(-55.0, 0.0))?;
    triangulation.insert(VertexType::new(-32.0, -42.0))?;
    triangulation.insert(VertexType::new(-2.0, -32.0))?;
    triangulation.insert(VertexType::new(25.0, -32.0))?;

    triangulation.insert(VertexType::new(70.0, 40.0))?;
    triangulation.insert(VertexType::new(20.0, 65.0))?;
    triangulation.insert(VertexType::new(-60.0, 50.0))?;
    triangulation.insert(VertexType::new(-70.0, -15.5))?;
    triangulation.insert(VertexType::new(-50.0, -60.0))?;
    triangulation.insert(VertexType::new(-9.0, -70.0))?;
    triangulation.insert(VertexType::new(42.0, -50.0))?;

    for edge in triangulation.fixed_undirected_edges() {
        triangulation.undirected_edge_data_mut(edge).color = SketchColor::LIGHT_GRAY;
    }

    for face in triangulation.fixed_inner_faces() {
        triangulation.face_data_mut(face).fill = SketchFill::solid(SketchColor::ANTIQUE_WHITE);
    }

    let query_vertex = VertexType::new(-5.0, -5.0);
    let mut nns = Vec::new();
    triangulation
        .natural_neighbor()
        .get_weights(query_vertex.position(), &mut nns);

    for (vertex, _) in &nns {
        triangulation.vertex_data_mut(*vertex).color = Some(SketchColor::CRIMSON);
    }

    let mut result = convert_triangulation(&triangulation, &Default::default());

    for (index, (vertex, _)) in nns.iter().enumerate() {
        let vertex = triangulation.vertex(*vertex);
        result.add(
            SketchElement::text(format!("{index}"))
                .position(convert_point(vertex.position()))
                .font_size(2.5)
                .horizontal_alignment(HorizontalAlignment::Middle)
                .dy(0.9),
        );
    }

    for edge in triangulation.undirected_voronoi_edges() {
        let [v0, v1] = edge.vertices();

        if let (Some(v0), Some(v1)) = (v0.position(), v1.position()) {
            result.add(
                SketchElement::line(convert_point(v0), convert_point(v1))
                    .stroke_color(SketchColor::SALMON)
                    .stroke_width(0.5),
            );
        }
    }

    let mut circumcenters = Vec::new();
    for face in triangulation.inner_faces() {
        let circumcenter = convert_point(face.circumcenter());
        circumcenters.push(circumcenter);

        result.add(
            SketchElement::circle(circumcenter, 1.0)
                .fill(SketchFill::Solid(SketchColor::ROYAL_BLUE)),
        );
    }

    if !include_faces {
        result.add(
            SketchElement::circle(convert_point(query_vertex.position()), 1.0)
                .fill(SketchFill::Solid(SketchColor::RED)),
        );
    }

    let mut insertion_cell_points = Vec::new();
    let inserted = triangulation.insert(query_vertex)?;
    for edge in triangulation
        .vertex(inserted)
        .as_voronoi_face()
        .adjacent_edges()
    {
        let context = "Edge was infinite - is insertion position correct?";
        let from = edge.from().position().context(context)?;
        let to = edge.to().position().context(context)?;
        insertion_cell_points.push(convert_point(from));

        result.add(
            SketchElement::line(convert_point(from), convert_point(to))
                .stroke_color(SketchColor::ORANGE_RED),
        );
    }
    for cell_point in &insertion_cell_points {
        result.add(
            SketchElement::circle(*cell_point, 1.0)
                .fill(SketchFill::solid(SketchColor::DARK_GREEN)),
        );
    }

    if include_faces {
        let nn3 = convert_point(triangulation.vertex(nns[3].0).position());
        let nn4 = convert_point(triangulation.vertex(nns[4].0).position());
        let nn5 = convert_point(triangulation.vertex(nns[5].0).position());

        let last_edge = SketchElement::line(nn3, nn4)
            .with_arrow_end(ArrowType::FilledArrow)
            .stroke_color(SketchColor::ROYAL_BLUE)
            .shift_from(-2.3)
            .shift_to(-7.0);

        result.add(
            last_edge
                .create_adjacent_text("last_edge")
                .font_size(3.0)
                .dy(3.5),
        );
        result.add(last_edge);
        let stop_edge = SketchElement::line(nn4, nn5)
            .with_arrow_end(ArrowType::FilledArrow)
            .stroke_color(SketchColor::ROYAL_BLUE)
            .shift_from(-2.3)
            .shift_to(-7.0);
        result.add(
            stop_edge
                .create_adjacent_text("stop_edge")
                .font_size(3.0)
                .dy(3.5),
        );
        result.add(stop_edge);

        result.add(
            SketchElement::text("first")
                .position(insertion_cell_points[6] + Vector2::new(0.0, 0.0))
                .dy(3.3)
                .font_size(2.5),
        );
        result.add(
            SketchElement::text("last")
                .position(insertion_cell_points[0] + Vector2::new(-5.5, 0.0))
                .dy(3.0)
                .font_size(2.5),
        );

        result.add(
            SketchElement::text("c2")
                .position(circumcenters[4] + Vector2::new(-3.5, 0.0))
                .dy(-1.0)
                .font_size(2.5),
        );

        result.add(
            SketchElement::text("c1")
                .position(circumcenters[2])
                .dy(-2.0)
                .font_size(2.5),
        );

        result.add(
            SketchElement::text("c0")
                .position(circumcenters[3] + Vector2::new(1.5, 0.0))
                .dy(0.5)
                .font_size(2.5),
        );

        let path = SketchElement::path()
            .move_to(insertion_cell_points[6])
            .line_to(insertion_cell_points[0])
            .line_to(circumcenters[4])
            .line_to(circumcenters[2])
            .line_to(circumcenters[3])
            .close()
            .fill(SketchFill::solid(SketchColor::ORANGE))
            .opacity(0.75);

        result.add_with_layer(path, SketchLayer::EDGES);
    }

    result
        .set_view_box_min(Point2::new(-60.0, -45.0))
        .set_view_box_max(Point2::new(45.0, 60.0));

    Ok(result)
}

pub fn refinement_scenario(do_refine: bool) -> Sketch {
    let mut cdt = create_refinement_cdt();

    let parameters = RefinementParameters::default();

    if do_refine {
        cdt.refine(parameters);
    }

    convert_refinement_cdt(&mut cdt)
}

pub fn exclude_outer_faces_scenario(do_refine: bool) -> Sketch {
    let mut cdt = create_refinement_cdt();

    let num_additional_vertices = if do_refine { 500 } else { 0 };

    let parameters = RefinementParameters::<f64>::default()
        .exclude_outer_faces(true)
        .with_max_additional_vertices(num_additional_vertices);

    let result = cdt.refine(parameters);
    for face in &result.excluded_faces {
        cdt.face_data_mut(*face).fill = SketchFill::solid(SketchColor::TAN);
    }

    convert_refinement_cdt(&mut cdt)
}

fn convert_refinement_cdt(cdt: &mut Cdt) -> Sketch {
    for vertex in cdt.fixed_vertices() {
        cdt.vertex_data_mut(vertex).radius = 0.4;
    }
    for edge in cdt.fixed_undirected_edges() {
        if cdt.is_constraint_edge(edge) {
            cdt.undirected_edge_data_mut(edge).data_mut().color = SketchColor::DARK_RED;
        } else {
            cdt.undirected_edge_data_mut(edge).data_mut().color = SketchColor::DARK_GRAY;
        }
    }
    let mut sketch = convert_triangulation(cdt, &ConversionOptions::default());
    sketch.set_width(360);
    sketch
}

fn create_refinement_cdt() -> Cdt {
    let mut cdt = Cdt::new();
    let v0 = VertexType::new(0.0, 0.0);
    let v1 = VertexType::new(0.0, 100.0);
    let v2 = VertexType::new(100.0, 100.0);
    let v3 = VertexType::new(100.0, 0.0);

    cdt.insert(v0).unwrap();
    cdt.insert(v1).unwrap();
    cdt.insert(v2).unwrap();
    cdt.insert(v3).unwrap();

    let inner_vertices = [
        // "A" Shape
        VertexType::new(10.0, 75.0),
        VertexType::new(20.0, 75.0),
        VertexType::new(25.0, 60.0),
        VertexType::new(35.0, 60.0),
        VertexType::new(40.0, 75.0),
        VertexType::new(50.0, 75.0),
        VertexType::new(30.0, 20.0),
    ];

    cdt.add_constraint_edges(inner_vertices, true).unwrap();

    let inner_vertices = [
        // Inner "A"
        VertexType::new(25.0, 55.0),
        VertexType::new(35.0, 55.0),
        VertexType::new(30.0, 40.0),
    ];

    cdt.add_constraint_edges(inner_vertices, true).unwrap();

    let inner_vertices = [
        // "C" Shape bottom half
        VertexType::new(55.0, 55.0),
        VertexType::new(55.0, 65.0),
        VertexType::new(60.0, 75.0),
        VertexType::new(70.0, 85.0),
        VertexType::new(80.0, 90.0),
        VertexType::new(85.0, 90.0),
        VertexType::new(85.0, 80.0),
        VertexType::new(80.0, 80.0),
        VertexType::new(70.0, 75.0),
        VertexType::new(65.0, 65.0),
        VertexType::new(65.0, 55.0),
        // "C" Shape upper half
        VertexType::new(70.0, 45.0),
        VertexType::new(80.0, 40.0),
        VertexType::new(85.0, 40.0),
        VertexType::new(85.0, 30.0),
        VertexType::new(80.0, 30.0),
        VertexType::new(70.0, 35.0),
        VertexType::new(60.0, 45.0),
    ];

    cdt.add_constraint_edges(inner_vertices, true).unwrap();
    cdt
}

pub fn angle_limit_scenario(angle_limit_degrees: f64) -> Sketch {
    let mut cdt = create_angle_limit_cdt();

    cdt.refine(
        RefinementParameters::default().with_angle_limit(AngleLimit::from_deg(angle_limit_degrees)),
    );
    let mut result = convert_refinement_cdt(&mut cdt);
    result.set_width(200);
    result
}

fn create_angle_limit_cdt() -> Cdt {
    let mut cdt = Cdt::new();

    let num_slices = 22;

    let mut vertices = vec![VertexType::new(0.0, 0.0)];

    for index in 0..num_slices {
        if index == 2 || index == 5 || index == 6 || index == 15 || index == 16 || index == 17 {
            // Add some arbitrary irregularities to make the result look more interesting
            continue;
        }

        let angle = std::f64::consts::PI * 0.9 * index as f64 / num_slices as f64;
        let distance = 50.0;
        let (sin, cos) = angle.sin_cos();
        vertices.push(VertexType::new(sin * distance, cos * distance));
    }

    cdt.add_constraint_edges(vertices, true).unwrap();
    cdt
}

pub fn refinement_maximum_area_scenario(max_area: Option<f64>) -> Sketch {
    let triangulation = big_triangulation().unwrap();
    let mut cdt = Cdt::from(triangulation);

    cdt.refine(RefinementParameters::default().with_max_allowed_area(max_area.unwrap_or(2000.)));

    let mut result = convert_refinement_cdt(&mut cdt);

    result.set_width(290);
    let description = if let Some(max_area) = max_area {
        format!("Max area: {max_area}")
    } else {
        "No max area".to_string()
    };

    result.add(
        SketchElement::text(description)
            .font_size(10.0)
            .angle(Deg(-33.2))
            .position(Point2::new(-75.0, -20.0)),
    );

    result
}
