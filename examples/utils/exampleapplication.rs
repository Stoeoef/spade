use utils::{Vertex};
use rtree::{RTree};
use cgmath::{Vector2};
use glium::{DisplayBuild, Surface, VertexBuffer, Program, Display, DrawParameters};
use glium::glutin::{Event, ElementState};
use glium;

const VERTEX_SHADER_SRC: &'static str = r#"
    #version 140
    in vec2 pos;
    in vec3 color;

    out vec3 fragment_color;
    void main() {
    gl_Position = vec4(pos, 0.0, 1.0);
    fragment_color = color;
        }
    "#;

const FRAGMENT_SHADER_SRC: &'static str = r#"
    #version 140
    out vec4 out_color;
    in vec3 fragment_color;
    void main() {
    out_color = vec4(fragment_color, 1.0);
        }
    "#;


pub struct ExampleApplication {
    pub tree: RTree<Vector2<f32>>,
    program: Program,
    pub edges_buffer: VertexBuffer<Vertex>,
    pub vertices_buffer: VertexBuffer<Vertex>,
    pub selection_buffer: VertexBuffer<Vertex>,
    pub display: Display,
}

impl ExampleApplication {
    pub fn new() -> ExampleApplication {
        let display = glium::glutin::WindowBuilder::new()
            .with_dimensions(800, 800)
            .with_title("RTree Demo".to_string())
            .build_glium()
            .unwrap();
        let program = Program::from_source(&display, VERTEX_SHADER_SRC,
                                           FRAGMENT_SHADER_SRC, None).unwrap();
        let tree_edges_buffer = VertexBuffer::new(&display, &[]).unwrap();
        let tree_vertices_buffer = VertexBuffer::new(&display, &[]).unwrap();
        let selection_buffer = VertexBuffer::new(&display, &[]).unwrap();
        ExampleApplication {
            tree: RTree::new(),
            display: display,
            program: program,
            edges_buffer: tree_edges_buffer,
            vertices_buffer: tree_vertices_buffer,
            selection_buffer: selection_buffer,
        }
    }
    
    pub fn default_handle_event(&mut self, event: &Event) -> bool {
        use glium::glutin::VirtualKeyCode::Escape;
        match event {
            &Event::Refresh => self.draw(),
            &Event::Closed | &Event::KeyboardInput(ElementState::Pressed, _, Some(Escape))
                => return true,
            _ => ()
        }
        false
    }

    pub fn draw(&self) {
        let mut target = self.display.draw();
        target.clear_color(1.0, 1.0, 1.0, 1.0);
        let indices = glium::index::NoIndices(glium::index::PrimitiveType::LinesList);
        let parameters = DrawParameters {
            line_width: Some(1.0),
            .. Default::default()
        };


        target.draw(&self.edges_buffer, &indices, &self.program, 
                    &glium::uniforms::EmptyUniforms, &parameters).unwrap();

        let parameters = DrawParameters {
            point_size: Some(3.0),
            line_width: Some(2.0),
            .. Default::default()
        };

        target.draw(&self.selection_buffer, &indices, &self.program,
                    &glium::uniforms::EmptyUniforms, &parameters).unwrap();

        let indices = glium::index::NoIndices(glium::index::PrimitiveType::Points);
        target.draw(&self.vertices_buffer, &indices, &self.program,
                    &glium::uniforms::EmptyUniforms, &parameters).unwrap();

        target.finish().unwrap();
    }
}
