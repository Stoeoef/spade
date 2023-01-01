use super::{
    ArrowType, PathPoint, Point, Sketch, SketchCircle, SketchColor, SketchElement, SketchLine,
    SketchPath, SketchText, Style,
};
use cgmath::num_traits::zero;
use cgmath::{Bounded, Deg, InnerSpace, Vector2};
use svg::node::element::path::Data;
use svg::node::element::{
    Circle, Definitions, Group, Line, Marker, Path, Pattern, Rectangle, Text,
};
use svg::{node::Text as TextNode, Document};

pub struct SketchConverter {
    unique_prefix: String,
    unique_id: u32,
    patterns: Vec<Pattern>,
    markers: Vec<Marker>,
}

impl SketchConverter {
    const FILLED_ARROW: &'static str = "filled_arrow";

    pub fn new<S: Into<String>>(unique_prefix: S) -> Self {
        Self {
            unique_prefix: unique_prefix.into(),
            unique_id: 0,
            patterns: Vec::new(),
            markers: Vec::new(),
        }
    }

    pub fn get_unique_pattern_id(&mut self) -> String {
        self.unique_id += 1;
        format!("{}_pattern{:04}", self.unique_prefix, self.unique_id)
    }

    pub fn get_unique_rotation(&mut self) -> Deg<f64> {
        const PATTERN_ROTATION: f64 = 63.512;
        Deg((self.unique_id as f64) * PATTERN_ROTATION)
    }

    pub fn add_pattern(&mut self, pattern: Pattern) {
        self.patterns.push(pattern)
    }
}

impl SketchConverter {
    pub fn convert<S: Into<String>>(unique_prefix: S, sketch: &Sketch) -> Document {
        let mut svg = Group::new();
        let mut converter = SketchConverter::new(unique_prefix);

        let (mut view_box_min, view_box_max) =
            Self::get_items_bounding_box(sketch.items.iter().flatten());

        view_box_min = sketch.view_box_min.unwrap_or(view_box_min);
        let view_box_max = sketch.view_box_max.unwrap_or(view_box_max);

        let mut dimensions = view_box_max - view_box_min;

        let offset = dimensions * sketch.relative_padding;
        view_box_min -= offset;
        dimensions += offset * 2.0;

        let background_rectangle = Rectangle::new()
            .set("x", view_box_min.x)
            .set("y", view_box_min.y)
            .set("width", "100%")
            .set("height", "100%")
            .set("style", sketch.style.get_attribute_string(&mut converter));

        for element in sketch.items.iter().flatten() {
            svg = converter.convert_item(svg, element);
        }

        let mut definitions = Definitions::new();
        for pattern in &converter.patterns {
            definitions = definitions.add(pattern.clone());
        }

        for marker in &converter.markers {
            definitions = definitions.add(marker.clone());
        }

        let mut document = Document::new().set("width", sketch.width);
        if let Some(height) = sketch.height {
            document = document.set("height", height)
        }

        document
            .set(
                "viewBox",
                format!(
                    "{} {} {} {}",
                    view_box_min.x, view_box_min.y, dimensions.x, dimensions.y
                ),
            )
            .set("xmlns:xlink", "http://www.w3.org/1999/xlink")
            .add(definitions)
            .add(background_rectangle)
            .add(svg)
    }

    fn convert_item(&mut self, svg: Group, element: &SketchElement) -> Group {
        match element {
            SketchElement::Text(SketchText {
                text,
                position,
                font_size,
                horizontal_alignment,
                angle,
                dy,
                style,
            }) => {
                let text_node = TextNode::new(text);
                let mut text = Text::new()
                    .set("x", position.x)
                    .set("y", position.y)
                    .set("style", style.get_attribute_string(self))
                    .set("font-size", *font_size)
                    .set(
                        "transform",
                        format!(
                            "translate ({} {}) rotate ({}) translate({} {})",
                            position.x, position.y, angle.0, -position.x, -position.y
                        ),
                    )
                    .set("dy", *dy)
                    .add(text_node);
                if let Some(alignment) = horizontal_alignment {
                    text = text.set("text-anchor", alignment.get_svg_value())
                }

                svg.add(text)
            }
            SketchElement::Circle(SketchCircle {
                center,
                radius,
                style,
            }) => svg.add(
                Circle::new()
                    .set("cx", center.x)
                    .set("cy", center.y)
                    .set("r", *radius)
                    .set("style", style.get_attribute_string(self)),
            ),
            SketchElement::Path(SketchPath { data_points, style }) => {
                let mut data = Data::new();
                for point in data_points {
                    match point {
                        PathPoint::MoveTo(point) => data = data.move_to((point.x, point.y)),
                        PathPoint::MoveBy(vector) => data = data.move_by((vector.x, vector.y)),
                        PathPoint::LineTo(point) => data = data.line_to((point.x, point.y)),
                        PathPoint::LineBy(vector) => data = data.line_by((vector.x, vector.y)),
                        PathPoint::ArcTo {
                            radii,
                            rotation,
                            large_arc,
                            sweep,
                            to,
                        } => {
                            data = data.elliptical_arc_to((
                                radii.x,
                                radii.y,
                                rotation.0,
                                i32::from(*large_arc),
                                i32::from(*sweep),
                                to.x,
                                to.y,
                            ));
                        }
                        PathPoint::Close => data = data.close(),
                    }
                }
                let path = Path::new()
                    .set("d", data)
                    .set("style", style.get_attribute_string(self));

                svg.add(path)
            }
            SketchElement::Line(SketchLine {
                from,
                to,
                draw_double_line,
                style,
                shift_from,
                shift_to,
            }) => {
                let normalized = (to - from).normalize();

                let from = from - normalized * *shift_from;
                let to = to + normalized * *shift_to;

                if *draw_double_line {
                    const FILLED_ARROW_CLEARANCE: f64 = 0.7;
                    const HALF_ARROW_CLEARANCE: f64 = 2.5;

                    let mut end_marker_clearance = zero();
                    let mut start_marker_clearance = zero();

                    match style.marker_end {
                        Some(ArrowType::FilledArrow) => {
                            end_marker_clearance = normalized * FILLED_ARROW_CLEARANCE;
                        }
                        Some(ArrowType::HalfArrow) => {
                            start_marker_clearance = normalized * HALF_ARROW_CLEARANCE;
                        }
                        None => {}
                    }

                    match style.marker_start {
                        Some(ArrowType::FilledArrow) => {
                            start_marker_clearance += normalized * FILLED_ARROW_CLEARANCE;
                        }
                        Some(ArrowType::HalfArrow) => {
                            end_marker_clearance += normalized * HALF_ARROW_CLEARANCE;
                        }
                        None => {}
                    }

                    let perpendicular = Vector2::new(normalized.y, -normalized.x);
                    let line_width = style.stroke_width.unwrap_or(1.0);

                    let offset = perpendicular * line_width * 0.8;
                    let with_line = self.add_svg_line(
                        svg,
                        from - offset + end_marker_clearance,
                        to - offset - start_marker_clearance,
                        style,
                    );
                    self.add_svg_line(
                        with_line,
                        to + offset - end_marker_clearance,
                        from + offset + start_marker_clearance,
                        style,
                    )
                } else {
                    self.add_svg_line(svg, from, to, style)
                }
            }
        }
    }

    fn add_svg_line(&mut self, svg: Group, from: Point, to: Point, style: &Style) -> Group {
        svg.add(
            Line::new()
                .set("x1", from.x)
                .set("y1", from.y)
                .set("x2", to.x)
                .set("y2", to.y)
                .set("style", style.get_attribute_string(self)),
        )
    }

    fn get_items_bounding_box<'a, I: Iterator<Item = &'a SketchElement>>(
        items: I,
    ) -> (Point, Point) {
        let mut min = Point::max_value();
        let mut max = Point::min_value();

        for item in items {
            let (item_min, item_max) = Self::get_item_bounding_box(item);
            min = min.zip(item_min, f64::min);
            max = max.zip(item_max, f64::max);
        }
        (min, max)
    }

    fn get_item_bounding_box(item: &SketchElement) -> (Point, Point) {
        match item {
            SketchElement::Text(SketchText {
                text,
                position,
                font_size,
                ..
            }) => {
                const HEIGHT_MULTIPLIER: f64 = 80.0;
                const WIDTH_MULTIPLIER: f64 = 10.5;

                let font_size_scale = font_size / 100.0;

                let char_len = text.chars().count() as f64;
                let dimensions =
                    Vector2::new(WIDTH_MULTIPLIER * char_len, HEIGHT_MULTIPLIER) * font_size_scale;
                (
                    Point::new(position.x, position.y - dimensions.y),
                    Point::new(position.x + dimensions.x, position.y),
                )
            }
            SketchElement::Circle(SketchCircle { center, radius, .. }) => {
                let offset = Vector2::new(*radius, *radius);
                (center - offset, center + offset)
            }
            SketchElement::Path(SketchPath { data_points, .. }) => {
                let mut min = Point::max_value();
                let mut max = Point::min_value();
                let mut last_point = None;
                let mut loop_start = None;

                for data_point in data_points {
                    last_point = match data_point {
                        PathPoint::MoveTo(point) => {
                            loop_start = Some(*point);
                            loop_start
                        }
                        PathPoint::MoveBy(vector) => {
                            loop_start = Some(
                                last_point.expect("MoveBy must not be the first path element")
                                    + vector,
                            );
                            loop_start
                        }
                        PathPoint::LineTo(point) => {
                            if last_point.is_none() {
                                panic!("Path must start with a move command");
                            }
                            Some(*point)
                        }
                        PathPoint::LineBy(vector) => {
                            Some(last_point.expect("Path must start with a move command") + vector)
                        }
                        PathPoint::ArcTo { to, .. } => Some(*to),
                        PathPoint::Close => loop_start,
                    };
                    if let Some(last_point) = last_point {
                        min = min.zip(last_point, f64::min);
                        max = max.zip(last_point, f64::max);
                    }
                }
                (min, max)
            }
            SketchElement::Line(SketchLine { from, to, .. }) => {
                (from.zip(*to, f64::min), from.zip(*to, f64::max))
            }
        }
    }

    pub fn add_arrowhead(
        &mut self,
        style: &Style,
        arrow_type: ArrowType,
        flip_arrow: bool,
    ) -> String {
        self.unique_id += 1;
        let arrow_url = format!(
            "{}_{}{:04}",
            self.unique_prefix,
            Self::FILLED_ARROW,
            self.unique_id
        );

        let mut path = match arrow_type {
            ArrowType::FilledArrow => Path::new().set(
                "d",
                Data::new()
                    .move_to((0.0, 0.0))
                    .line_to((5.0, 1.5))
                    .line_to((0.0, 3.0))
                    .close(),
            ),
            ArrowType::HalfArrow => Path::new().set(
                "d",
                Data::new()
                    .move_to((1.0, 0.0))
                    .line_to((4.0, 0.0))
                    .line_to((0.0, 2.0))
                    .close(),
            ),
        };

        if flip_arrow {
            path = path.set("transform", "scale(-1 1) translate (-5 0)");
        }

        let mut marker = Marker::new()
            .set("id", arrow_url.clone())
            .set("markerWidth", 5)
            .set("markerHeight", 3)
            .set("stroke", "none")
            .set("orient", "auto")
            .set(
                "fill",
                style.stroke_color.unwrap_or(SketchColor::BLACK).to_string(),
            )
            .add(path);

        marker = match arrow_type {
            ArrowType::FilledArrow => marker
                .set("refX", if flip_arrow { 4.9 } else { 0.1 })
                .set("refY", 1.5),
            ArrowType::HalfArrow => marker
                .set("refX", if flip_arrow { 4.0 } else { 1.0 })
                .set("refY", 0.5),
        };

        self.markers.push(marker);
        arrow_url
    }
}
