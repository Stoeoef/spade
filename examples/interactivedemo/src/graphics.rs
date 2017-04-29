// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use ::{ExampleTriangulation};
use spade::{BoundingRect, HasPosition};
use spade::rtree::{RTree, RTreeNode};
use spade::delaunay::Subdivision;
use cgmath::{EuclideanSpace, Point2, Point3, Vector2, Vector3, Array};
use cgmath::conv::{array2, array3};
use glium::{Surface, VertexBuffer, Program, Display, DrawParameters};
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

pub struct RenderData {
    program: Program,
    pub edges_buffer: VertexBuffer<Vertex>,
    pub vertices_buffer: VertexBuffer<Vertex>,
    pub selection_buffer: VertexBuffer<Vertex>,
}

impl RenderData {
    pub fn new(display: &Display) -> RenderData {
        let program = Program::from_source(display, VERTEX_SHADER_SRC,
                                           FRAGMENT_SHADER_SRC, None).unwrap();
        let edges_buffer = VertexBuffer::new(display, &[]).unwrap();
        let vertices_buffer = VertexBuffer::new(display, &[]).unwrap();
        let selection_buffer = VertexBuffer::new(display, &[]).unwrap();
        RenderData {
            program: program,
            edges_buffer: edges_buffer,
            vertices_buffer: vertices_buffer,
            selection_buffer: selection_buffer,
        }
    }

    pub fn draw(&self, display: &Display) {
        let mut target = display.draw();
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

    pub fn update_buffers(&mut self, display: &Display, tree: &RTree<Point2<f64>>, 
                          delaunay: &ExampleTriangulation, draw_tree_nodes: bool) {
        let mut edges = Vec::new();
        let vertices = get_tree_edges(&tree, &mut edges);
        if !draw_tree_nodes {
            edges.clear();
        }
        get_delaunay_edges(&delaunay, &mut edges);
        self.edges_buffer = VertexBuffer::new(display, &edges).unwrap();
        self.vertices_buffer = VertexBuffer::new(display, &vertices).unwrap();
    }

    pub fn update_selection(&mut self, display: &Display, points: &Vec<Point2<f64>>) {
        let color = [1.0, 0.0, 0.0];
        let mut vertices = Vec::new();
        for point in points {
            push_cross(&mut vertices, point, color);
        }
        self.selection_buffer = VertexBuffer::new(display, &vertices).unwrap();
    }
}

#[derive(Clone, Copy)]
pub struct Vertex {
    pub pos: [f32; 2],
    pub color: [f32; 3],
}

implement_vertex!(Vertex, pos, color);
impl Vertex {
    pub fn new(pos: [f32; 2], color: [f32; 3]) -> Vertex {
        Vertex { pos: pos, color: color }
    }
}

pub fn push_rectangle(vec: &mut Vec<Vertex>, rect: &BoundingRect<Point2<f64>>, 
                      color: [f32; 3]) {
    let v0: Point2<_> = rect.lower();
    let v2: Point2<_> = rect.upper();
    let v1 = Point2::new(v2.x, v0.y);
    let v3 = Point2::new(v0.x, v2.y);
    vec.extend([v0, v1, v1, v2, v2, v3, v3, v0].iter().cloned().map(
        |v| Vertex::new(array2(v.to_vec().cast()), array3(color))));
}

pub fn push_cross(vec: &mut Vec<Vertex>, pos: &Point2<f64>, color: [f32; 3]) {
    let mut delta =  Vector2::from_value(0.015);
    let v0 = pos.to_vec() + delta;
    let v1 = pos.to_vec() - delta;
    delta.x *= -1.0;
    let v2 = pos.to_vec() + delta;
    let v3 = pos.to_vec() - delta;
    vec.extend([v0, v1, v2, v3].iter().cloned().map(
        |v| Vertex::new(array2(v.cast()), color)));
}

pub fn get_color_for_depth(depth: usize) -> [f32; 3] {
    array3(match depth {
        0 => Vector3::new(1., 1., 1.),
        1 => Vector3::new(1., 0., 1.) * 0.85,
        2 => Vector3::new(0., 1., 1.) * 0.7,
        3 => Vector3::new(0., 0., 1.) * 0.55,
        4 => Vector3::new(1., 0., 0.) * 0.4,
        _ => Vector3::new(0., 1., 1.) * 0.25,
    })
}

fn get_tree_edges(tree: &RTree<Point2<f64>>, buffer: &mut Vec<Vertex>) -> Vec<Vertex> {
    let mut vertices = Vec::new();
    let vertex_color = Point3::new(0.0, 0.0, 1.0);
    let mut to_visit = vec![tree.root()];
    while let Some(cur) = to_visit.pop() {
        for child in cur.children().iter() {
            match child {
                &RTreeNode::Leaf(point) => vertices.push(Vertex::new(
                    array2(point.to_vec().cast()), array3(vertex_color.clone()))),
                &RTreeNode::DirectoryNode(ref data) => {
                    to_visit.push(data);
                    push_rectangle(buffer, &data.mbr(), get_color_for_depth(data.depth()));
                }                
            }
        }
    }
    vertices
}

fn get_delaunay_edges(delaunay: &ExampleTriangulation,
                      edges_buffer: &mut Vec<Vertex>) {
    let base_color = [0.1, 0.1, 0.2];
    let constraint_color = [0.96, 0.1, 0.1];
    for edge in delaunay.edges() {
        let color = if delaunay.is_constraint_edge(edge.fix()) {
            constraint_color
        } else {
            base_color
        };
        let from = edge.from().position().to_vec();
        let to = edge.to().position().to_vec();
        edges_buffer.push(Vertex::new(array2(from.cast()), color));
        edges_buffer.push(Vertex::new(array2(to.cast()), color));
    }

    // let color = [0.8, 0.2, 0.2];
    // // for &(from, to) in &delaunay.constraints {
    // //     let from = delaunay.vertex(from).position().to_vec();
    // //     let to = delaunay.vertex(to).position().to_vec();
    // //     edges_buffer.push(Vertex::new(array2(from.cast()), color));
    // //     edges_buffer.push(Vertex::new(array2(to.cast()), color));
    // // }

    // for &(from, to) in &delaunay.conflicts {
    //     let from = delaunay.vertex(from).position().to_vec();
    //     let to = delaunay.vertex(to).position().to_vec();
    //     let color = [0.8, 0.5, 0.5];
    //     edges_buffer.push(Vertex::new(array2(from.cast()), color));
    //     let color = [0.8, 0.0, 0.0];
    //     edges_buffer.push(Vertex::new(array2(to.cast()), color));
    // }

}
