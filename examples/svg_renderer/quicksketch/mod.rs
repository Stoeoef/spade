//! This module contains a small library which contain utilities for svg conversion.
//! Feel free to read on but be warned that it's unrelated on how to use Spade.
use cgmath::{Angle, Deg, EuclideanSpace, Zero};
use convert::SketchConverter;
use svg::node::element::{Line, Pattern};
mod color;
mod convert;

pub use color::SketchColor;

pub type Point = cgmath::Point2<f64>;
pub type Vector = cgmath::Vector2<f64>;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum SketchFill {
    Solid(SketchColor),
    StripePattern(SketchColor, SketchColor, f64),
}

impl From<SketchColor> for SketchFill {
    fn from(color: SketchColor) -> Self {
        Self::Solid(color)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub enum HorizontalAlignment {
    Left,
    Middle,
    Right,
}

impl HorizontalAlignment {
    pub(crate) fn get_svg_value(&self) -> &'static str {
        match self {
            HorizontalAlignment::Left => "left",
            HorizontalAlignment::Middle => "middle",
            HorizontalAlignment::Right => "right",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub enum StrokeStyle {
    Dashed,
    SmallDashed,
    Regular,
}

impl StrokeStyle {
    pub(crate) fn get_style_attribute(&self) -> &'static str {
        match self {
            StrokeStyle::Dashed => "4 2",
            StrokeStyle::Regular => "1.0",
            StrokeStyle::SmallDashed => "1 0.5",
        }
    }
}

impl SketchFill {
    pub fn stripes(color1: SketchColor, color2: SketchColor, stripe_width: f64) -> Self {
        SketchFill::StripePattern(color1, color2, stripe_width)
    }

    pub fn solid(color: SketchColor) -> Self {
        SketchFill::Solid(color)
    }

    pub(crate) fn get_attribute_string(&self, converter: &mut SketchConverter) -> String {
        match self {
            SketchFill::Solid(color) => color.to_string(),
            SketchFill::StripePattern(color1, color2, scale) => {
                let pattern_id = converter.get_unique_pattern_id();
                let pattern_rotation = converter.get_unique_rotation();

                let transform = format!("rotate({} 0 0) scale({})", pattern_rotation.0, scale);

                converter.add_pattern(
                    Pattern::new()
                        .set("id", pattern_id.clone())
                        .set("width", 1)
                        .set("height", 1)
                        .set("patternTransform", transform)
                        .set("patternUnits", "userSpaceOnUse")
                        .add(
                            Line::new()
                                .set("x1", 0)
                                .set("y1", -0.1)
                                .set("x2", 0)
                                .set("y2", 1.1)
                                .set("stroke-width", 1.01)
                                .set("stroke", color1.to_string()),
                        )
                        .add(
                            Line::new()
                                .set("x1", 1)
                                .set("y1", -0.1)
                                .set("x2", 1)
                                .set("y2", 1.1)
                                .set("stroke-width", 1.01)
                                .set("stroke", color2.to_string()),
                        ),
                );
                format!("url(#{})", pattern_id)
            }
        }
    }
}

impl Default for SketchFill {
    fn default() -> Self {
        SketchFill::Solid(SketchColor::GRAY)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ArrowType {
    FilledArrow,
    HalfArrow,
}

#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Style {
    stroke_width: Option<f64>,
    stroke_color: Option<SketchColor>,
    stroke_style: Option<StrokeStyle>,
    fill: Option<SketchFill>,
    marker_start: Option<ArrowType>,
    marker_end: Option<ArrowType>,
}

impl Style {
    pub fn new() -> Self {
        Default::default()
    }

    fn get_attribute_string(&self, converter: &mut SketchConverter) -> String {
        let stroke_width = self
            .stroke_width
            .map(|width| format!("stroke-width: {}", width));

        let stroke_color = self
            .stroke_color
            .as_ref()
            .map(|color| format!("stroke: {}", color));

        let fill = format!(
            "fill: {}",
            self.fill
                .as_ref()
                .map(|fill| fill.get_attribute_string(converter))
                .unwrap_or_else(|| "none".to_string())
        );

        let stroke_dash_array = self.stroke_style.map(|stroke_style| {
            format!("stroke-dasharray: {}", stroke_style.get_style_attribute())
        });

        let marker_start = self.marker_start.as_ref().map(|arrow_type| {
            format!(
                "marker-start: url(#{})",
                converter.add_arrowhead(self, *arrow_type, true)
            )
        });

        let marker_end = self.marker_end.as_ref().map(|arrow_type| {
            format!(
                "marker-end: url(#{})",
                converter.add_arrowhead(self, *arrow_type, false)
            )
        });

        IntoIterator::into_iter([
            stroke_width,
            stroke_color,
            marker_start,
            marker_end,
            stroke_dash_array,
        ])
        .flatten()
        .chain(std::iter::once(fill))
        .collect::<Vec<_>>()
        .join("; ")
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SketchCircle {
    center: Point,
    radius: f64,
    style: Style,
}

impl SketchCircle {
    pub fn fill(mut self, fill: impl Into<SketchFill>) -> Self {
        self.style.fill = Some(fill.into());
        self
    }

    pub fn stroke_width(mut self, width: f64) -> Self {
        self.style.stroke_width = Some(width);
        self
    }

    pub fn stroke_color(mut self, color: SketchColor) -> Self {
        self.style.stroke_color = Some(color);
        self
    }

    pub fn stroke_style(mut self, stroke_style: StrokeStyle) -> Self {
        self.style.stroke_style = Some(stroke_style);
        self
    }
}

impl Default for SketchCircle {
    fn default() -> Self {
        SketchCircle {
            center: Point::new(0.0, 0.0),
            radius: 10.0,
            style: Default::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SketchText {
    text: String,
    position: Point,
    horizontal_alignment: Option<HorizontalAlignment>,
    font_size: f64,
    angle: Deg<f64>,
    dy: f64,
    style: Style,
}

impl SketchText {
    pub fn position(mut self, position: Point) -> Self {
        self.position = position;
        self
    }

    pub fn stroke_width(mut self, stroke_width: f64) -> Self {
        self.style.stroke_width = Some(stroke_width);
        self
    }

    pub fn stroke_color(mut self, stroke_color: SketchColor) -> Self {
        self.style.stroke_color = Some(stroke_color);
        self
    }

    pub fn fill(mut self, fill: SketchFill) -> Self {
        self.style.fill = Some(fill);
        self
    }

    pub fn horizontal_alignment(mut self, horizontal_alignment: HorizontalAlignment) -> Self {
        self.horizontal_alignment = Some(horizontal_alignment);
        self
    }

    pub fn font_size(mut self, font_size: f64) -> Self {
        self.font_size = font_size;
        self
    }

    pub fn angle(mut self, angle: Deg<f64>) -> Self {
        self.angle = angle;
        self
    }

    pub fn dy(mut self, dy: f64) -> Self {
        self.dy = dy;
        self
    }
}

impl Default for SketchText {
    fn default() -> Self {
        Self {
            text: Default::default(),
            position: Point::new(0.0, 0.0),
            font_size: 20.0,
            horizontal_alignment: None,
            angle: Deg::zero(),
            dy: 0.0,
            style: Default::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub enum PathPoint {
    MoveTo(Point),
    MoveBy(Vector),
    LineTo(Point),
    LineBy(Vector),
    ArcTo {
        radii: Vector,
        rotation: Deg<f64>,
        large_arc: bool,
        sweep: bool,
        to: Point,
    },
    Close,
}

#[derive(Clone, Debug)]
pub struct SketchPath {
    data_points: Vec<PathPoint>,
    style: Style,
}

impl SketchPath {
    pub fn move_to(mut self, point: Point) -> Self {
        self.data_points.push(PathPoint::MoveTo(point));
        self
    }

    pub fn move_by(mut self, offset: Vector) -> Self {
        self.data_points.push(PathPoint::MoveBy(offset));
        self
    }

    pub fn line_to(mut self, point: Point) -> Self {
        self.data_points.push(PathPoint::LineTo(point));
        self
    }

    pub fn line_by(mut self, offset: Vector) -> Self {
        self.data_points.push(PathPoint::LineBy(offset));
        self
    }

    pub fn arc_to(
        mut self,
        radii: Vector,
        rotation: Deg<f64>,
        large_arc: bool,
        sweep: bool,
        to: Point,
    ) -> Self {
        self.data_points.push(PathPoint::ArcTo {
            radii,
            rotation,
            large_arc,
            sweep,
            to,
        });
        self
    }

    pub fn close(mut self) -> Self {
        self.data_points.push(PathPoint::Close);
        self
    }

    pub fn fill(mut self, fill: SketchFill) -> Self {
        self.style.fill = Some(fill);
        self
    }

    pub fn stroke_width(mut self, width: f64) -> Self {
        self.style.stroke_width = Some(width);
        self
    }

    pub fn stroke_color(mut self, color: SketchColor) -> Self {
        self.style.stroke_color = Some(color);
        self
    }

    pub fn with_arrow_end(mut self, arrow_type: ArrowType) -> Self {
        self.style.marker_end = Some(arrow_type);
        self
    }

    pub fn with_arrow_start(mut self, arrow_type: ArrowType) -> Self {
        self.style.marker_start = Some(arrow_type);
        self
    }

    pub fn stroke_style(mut self, stroke_style: StrokeStyle) -> Self {
        self.style.stroke_style = Some(stroke_style);
        self
    }
}

#[derive(Clone, Debug)]
pub struct SketchLine {
    from: Point,
    to: Point,
    shift_from: f64,
    shift_to: f64,
    draw_double_line: bool,
    style: Style,
}

impl SketchLine {
    pub fn draw_double_line(mut self) -> Self {
        self.draw_double_line = true;
        self
    }

    pub fn with_arrow_end(mut self, arrow_type: ArrowType) -> Self {
        self.style.marker_end = Some(arrow_type);
        self
    }

    pub fn with_arrow_start(mut self, arrow_type: ArrowType) -> Self {
        self.style.marker_start = Some(arrow_type);
        self
    }

    pub fn shift_from(mut self, shift_from: f64) -> Self {
        self.shift_from = shift_from;
        self
    }

    pub fn shift_to(mut self, shift_to: f64) -> Self {
        self.shift_to = shift_to;
        self
    }

    pub fn shift_from_and_to(self, shift: f64) -> Self {
        self.shift_from(shift).shift_to(shift)
    }

    pub fn stroke_width(mut self, width: f64) -> Self {
        self.style.stroke_width = Some(width);
        self
    }

    pub fn stroke_color(mut self, color: SketchColor) -> Self {
        self.style.stroke_color = Some(color);
        self
    }

    pub fn stroke_style(mut self, stroke_style: StrokeStyle) -> Self {
        self.style.stroke_style = Some(stroke_style);
        self
    }

    pub fn create_adjacent_text<S: Into<String>>(&self, text: S) -> SketchText {
        let diff = self.from - self.to;
        let mut rotation: Deg<f64> = Deg::atan2(diff.x, diff.y);

        if rotation > Deg(0.0) {
            rotation -= Deg(180.0);
        }

        let middle = EuclideanSpace::centroid(&[self.from, self.to]);
        let text = SketchElement::text(text);

        text.angle(-rotation - Deg(90.0))
            .position(middle)
            .horizontal_alignment(HorizontalAlignment::Middle)
    }
}

#[derive(Clone, Debug)]
pub enum SketchElement {
    Text(SketchText),
    Circle(SketchCircle),
    Path(SketchPath),
    Line(SketchLine),
}

impl SketchElement {
    pub fn circle(center: Point, radius: f64) -> SketchCircle {
        SketchCircle {
            center,
            radius,
            style: Default::default(),
        }
    }

    pub fn line(from: Point, to: Point) -> SketchLine {
        SketchLine {
            from,
            to,
            shift_from: 0.0,
            shift_to: 0.0,
            draw_double_line: false,
            style: Style {
                stroke_width: Some(1.0),
                stroke_color: Some(SketchColor::BLACK),
                ..Default::default()
            },
        }
    }

    pub fn path() -> SketchPath {
        SketchPath {
            data_points: Default::default(),
            style: Default::default(),
        }
    }

    pub fn text<S: Into<String>>(s: S) -> SketchText {
        SketchText {
            text: s.into(),
            style: Style {
                fill: Some(SketchFill::solid(SketchColor::BLACK)),
                ..Default::default()
            },
            ..Default::default()
        }
    }
}

impl From<SketchText> for SketchElement {
    fn from(text: SketchText) -> Self {
        SketchElement::Text(text)
    }
}

impl From<SketchCircle> for SketchElement {
    fn from(circle: SketchCircle) -> Self {
        SketchElement::Circle(circle)
    }
}

impl From<SketchPath> for SketchElement {
    fn from(polygon: SketchPath) -> Self {
        SketchElement::Path(polygon)
    }
}

impl From<SketchLine> for SketchElement {
    fn from(line: SketchLine) -> Self {
        SketchElement::Line(line)
    }
}

#[derive(Clone, Debug)]
pub struct Sketch {
    pub items: [Vec<SketchElement>; SketchLayer::TOP as usize + 1],
    width: u32,
    relative_padding: f64,
    view_box_min: Option<Point>,
    view_box_max: Option<Point>,
    height: Option<u32>,
    style: Style,
}

#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SketchLayer {
    BACKGROUND = 0,
    EDGES = 1,
    VERTICES = 2,
    TOP = 3,
}

impl Sketch {
    pub fn new() -> Self {
        const EMPTY: Vec<SketchElement> = Vec::new();

        Self {
            items: [EMPTY; SketchLayer::TOP as usize + 1],
            relative_padding: 0.1,
            width: 800,
            height: None,
            view_box_min: None,
            view_box_max: None,
            style: Style::default(),
        }
    }

    pub fn from_iterator(iterator: impl IntoIterator<Item = SketchElement>) -> Self {
        let mut result = Sketch::new();

        for element in iterator {
            result.add(element);
        }

        result
    }

    pub fn set_width(&mut self, width: u32) -> &mut Self {
        self.width = width;
        self
    }

    pub fn set_height(&mut self, height: u32) -> &mut Self {
        self.height = Some(height);
        self
    }

    pub fn set_background(&mut self, background: SketchFill) -> &mut Self {
        self.style.fill = Some(background);
        self
    }

    pub fn set_relative_padding(&mut self, relative_padding: f64) -> &mut Self {
        self.relative_padding = relative_padding;
        self
    }

    pub fn set_view_box_min(&mut self, view_box_min: Point) -> &mut Self {
        self.view_box_min = Some(view_box_min);
        self
    }

    pub fn set_view_box_max(&mut self, view_box_max: Point) -> &mut Self {
        self.view_box_max = Some(view_box_max);
        self
    }

    pub fn add<T: Into<SketchElement>>(&mut self, item: T) -> &mut Self {
        self.add_with_layer(item, SketchLayer::TOP)
    }

    pub fn add_with_layer<T: Into<SketchElement>>(
        &mut self,
        item: T,
        layer: SketchLayer,
    ) -> &mut Self {
        self.items[layer as usize].push(item.into());
        self
    }

    pub fn save_to_svg<S: Into<String>, P: AsRef<std::path::Path>>(
        &self,
        unique_prefix: S,
        path: P,
    ) -> std::io::Result<()> {
        let document = SketchConverter::convert(unique_prefix, self);
        svg::save(path, &document)
    }
}

impl Default for Sketch {
    fn default() -> Self {
        Self::new()
    }
}
