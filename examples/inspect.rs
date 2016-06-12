extern crate rtree;
extern crate rand;
extern crate cgmath;
#[macro_use]
extern crate glium;

mod utils;
use utils::exampleapplication::ExampleApplication;
use utils::{Vertex, get_color_for_depth, push_rectangle, push_cross};
use rtree::{RTree, RTreeNode};
use cgmath::{Vector2, Vector3};
use cgmath::conv::*;
use glium::{VertexBuffer};
use glium::glutin::{Event};

fn get_tree_edges(tree: &RTree<Vector2<f32>>) -> (Vec<Vertex>, Vec<Vertex>) {
    let mut edges = Vec::new();
    let mut vertices = Vec::new();
    let vertex_color = Vector3::new(0.0, 0.0, 1.0);
    let mut to_visit = vec![tree.root()];
    while let Some(cur) = to_visit.pop() {
        for child in cur.children().iter() {
            match child {
                &RTreeNode::Leaf(ref point) => vertices.push(Vertex::new(
                    array2(point.clone()), array3(vertex_color.clone()))),
                &RTreeNode::DirectoryNode(ref data) => {
                    to_visit.push(data);
                    push_rectangle(&mut edges, &data.mbr(), &get_color_for_depth(data.depth()));
                }                
            }
        }
    }
    (vertices, edges)
}

fn main() {

    let mut app = ExampleApplication::new();
    populate_tree(&mut app.tree);
    let (vertices, edges) = get_tree_edges(&app.tree);

    app.edges_buffer = VertexBuffer::new(&app.display, &edges).unwrap();
    app.vertices_buffer = VertexBuffer::new(&app.display, &vertices).unwrap();

    loop {
        let events: Vec<_> = app.display.poll_events().collect();

        let mut dirty = false;    
        for event in events.into_iter() {
            if app.default_handle_event(&event) {
                return;
            }
            match event {
                Event::MouseMoved(x, y) => {
                    let (w, h) = app.display.get_framebuffer_dimensions();
                    // Transform x, y into the range [-1 , 1]
                    let y = h as i32 - y;
                    let x = (x as f32 / w as f32) * 2. - 1.;
                    let y = (y as f32 / h as f32) * 2. - 1.;
                    let nearest = app.tree.nearest_neighbor(array2(Vector2::new(x, y)));
                    let mut vertices = Vec::new();
                    push_cross(&mut vertices, nearest.unwrap(), &Vector3::new(1.0, 0.0, 0.0));
                    app.selection_buffer = VertexBuffer::new(&app.display, &vertices).unwrap();
                    dirty = true;
                },
                _ => (),
            }
        }
        if dirty {
            app.draw();
        }
    }
}

fn populate_tree(tree: &mut RTree<Vector2<f32>>) {
    let points = utils::random_points_with_seed(1000, [0, 1, 55, 0]);
    // let points = utils::random_walk(1000, 0.1, [0, 1, 55, 0]);
    for point in points.into_iter() {
        tree.insert(point);
    }
}
