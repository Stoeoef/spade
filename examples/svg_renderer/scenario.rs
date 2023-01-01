use super::quicksketch::{ArrowType, Point, Sketch, SketchColor, SketchElement, SketchFill};
use spade::{DelaunayTriangulation, HasPosition, Point2};

use crate::convert_point;

pub struct VertexType {
    position: Point,
    pub radius: f64,
}

impl VertexType {
    pub fn new(x: f64, y: f64) -> Self {
        const DEFAULT_CIRCLE_RADIUS: f64 = 2.0;

        Self {
            position: Point::new(x, y),
            radius: DEFAULT_CIRCLE_RADIUS,
        }
    }
}

impl HasPosition for VertexType {
    type Scalar = f64;

    fn position(&self) -> Point2<f64> {
        Point2::new(self.position.x, self.position.y)
    }
}

pub struct UndirectedEdgeType {
    pub color: SketchColor,
}

impl AsRef<UndirectedEdgeType> for UndirectedEdgeType {
    fn as_ref(&self) -> &UndirectedEdgeType {
        self
    }
}

impl Default for UndirectedEdgeType {
    fn default() -> Self {
        Self {
            color: SketchColor::BLACK,
        }
    }
}

#[derive(Default)]
pub struct DirectedEdgeType {}

#[derive(Clone, Copy, Debug)]
pub struct FaceType {
    pub fill: SketchFill,
}

impl Default for FaceType {
    fn default() -> Self {
        Self {
            fill: SketchFill::Solid(SketchColor::LIGHT_GRAY),
        }
    }
}

pub type Triangulation =
    DelaunayTriangulation<VertexType, DirectedEdgeType, UndirectedEdgeType, FaceType>;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum EdgeMode {
    Disabled,
    Undirected,
    Directed { reversed: bool },
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ConversionOptions {
    pub directed_edge_mode: EdgeMode,
    pub vertex_stroke_color: SketchColor,
    pub vertex_color: SketchColor,
}

impl Default for ConversionOptions {
    fn default() -> Self {
        Self {
            directed_edge_mode: EdgeMode::Undirected,
            vertex_stroke_color: SketchColor::BLACK,
            vertex_color: SketchColor::DIM_GRAY,
        }
    }
}

pub fn convert_triangulation<T>(triangulation: &T, options: &ConversionOptions) -> Sketch
where
    T: spade::Triangulation<Vertex = VertexType, Face = FaceType>,
    T::UndirectedEdge: AsRef<UndirectedEdgeType>,
{
    let mut sketch = Sketch::new();

    sketch.set_background(SketchFill::stripes(
        SketchColor::WHITE,
        SketchColor::LIGHT_GRAY,
        2.0,
    ));

    for face in triangulation.inner_faces() {
        let [p0, p1, p2] = face.positions();
        let p0 = convert_point(p0);
        let p1 = convert_point(p1);
        let p2 = convert_point(p2);
        sketch.add_with_layer(
            SketchElement::path()
                .move_to(p0)
                .line_to(p1)
                .line_to(p2)
                .close()
                .fill(face.data().fill),
            crate::quicksketch::SketchLayer::BACKGROUND,
        );
    }

    for undirected_edge in triangulation.undirected_edges() {
        let [mut from, mut to] = undirected_edge.positions();
        let color = undirected_edge.data().as_ref().color;

        const SHIFT: f64 = -3.0;
        match options.directed_edge_mode {
            EdgeMode::Directed { reversed } => {
                if reversed {
                    std::mem::swap(&mut from, &mut to);
                }

                let line = SketchElement::line(convert_point(from), convert_point(to))
                    .draw_double_line()
                    .shift_to(SHIFT)
                    .shift_from(SHIFT)
                    .stroke_width(0.7)
                    .stroke_color(color);

                let line = if reversed {
                    line.with_arrow_end(ArrowType::HalfArrow)
                } else {
                    line.with_arrow_start(ArrowType::HalfArrow)
                };

                sketch.add_with_layer(line, crate::quicksketch::SketchLayer::EDGES);
            }
            EdgeMode::Undirected => {
                sketch.add_with_layer(
                    SketchElement::line(convert_point(from), convert_point(to)).stroke_color(color),
                    crate::quicksketch::SketchLayer::EDGES,
                );
            }
            EdgeMode::Disabled => {}
        }
    }

    for vertex in triangulation.vertices() {
        sketch.add_with_layer(
            SketchElement::circle(convert_point(vertex.position()), vertex.data().radius)
                .fill(SketchFill::solid(options.vertex_color))
                .stroke_width(0.5)
                .stroke_color(options.vertex_stroke_color),
            crate::quicksketch::SketchLayer::VERTICES,
        );
    }

    sketch
}
