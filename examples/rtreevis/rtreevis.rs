// Copyright 2016 The Spade Developers. For a full listing of the authors,
// refer to the Cargo.toml file at the top-level directory of this distribution.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

extern crate spade;
extern crate rand;
extern crate cgmath;
#[macro_use]
extern crate glium;

mod utils;
use utils::exampleapplication::ExampleApplication;
use utils::{Vertex, get_color_for_depth, push_rectangle, push_cross};
use spade::{RTree, RTreeNode, SpadeNum};
use cgmath::{Vector2, Vector3, BaseFloat, BaseNum};
use cgmath::conv::*;
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use glium::{VertexBuffer};
use glium::glutin::{Event, ElementState, MouseButton};
use glium::glutin::VirtualKeyCode;

enum LookupMode {
    Nearest,
    NearestN,
    InCircle,
    CloseN,
}

impl std::fmt::Display for LookupMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", match self {
            &LookupMode::Nearest => "Nearest neighbor",
            &LookupMode::NearestN => "Nearest N neighbors",
            &LookupMode::InCircle => "Contained in Circle",
            &LookupMode::CloseN => "Close neighbor",
        })
    }
}

fn get_tree_edges(tree: &RTree<Vector2<f32>, Vector2<f32>>) -> (Vec<Vertex>, Vec<Vertex>) {
    let mut edges = Vec::new();
    let mut vertices = Vec::new();
    let vertex_color = Vector3::new(0.0, 0.0, 1.0);
    let mut to_visit = vec![tree.root()];
    while let Some(cur) = to_visit.pop() {
        for child in cur.children().iter() {
            match child {
                &RTreeNode::Leaf(ref point) => vertices.push(Vertex::new(
                    array2(*point), array3(vertex_color.clone()))),
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

    let mut last_point = Vector2::new(0., 0.);
    let mut lookup_mode = LookupMode::Nearest;
    let mut draw_tree_nodes = true;

    println!("RTree Demo");
    print_help();
    loop {
        let events: Vec<_> = app.display.poll_events().collect();

        let mut dirty = false;
        for event in events.into_iter() {
            if app.default_handle_event(&event) {
                return;
            }
            match event {
                Event::KeyboardInput(ElementState::Pressed, _, Some(key)) => {
                    match key {
                        VirtualKeyCode::H => {
                            print_help();
                        },
                        VirtualKeyCode::F => {
                            draw_tree_nodes = !draw_tree_nodes;
                            update_buffers(&mut app, draw_tree_nodes);
                            dirty = true;
                        },
                        VirtualKeyCode::M => {
                            match lookup_mode {
                                LookupMode::Nearest => lookup_mode = LookupMode::NearestN,
                                LookupMode::NearestN => lookup_mode = LookupMode::InCircle,
                                LookupMode::InCircle => lookup_mode = LookupMode::CloseN,
                                LookupMode::CloseN => lookup_mode = LookupMode::Nearest,
                            }
                            println!("Changed lookup mode to {}", lookup_mode);
                        },
                        VirtualKeyCode::A | VirtualKeyCode::B => {
                            // Insert some random points
                            let num = if key == VirtualKeyCode::A { 10usize } else { 100 };
                            let s = app.tree.size() as u32;
                            let seed = [s, s * 35, s + 124, s * 91];
                            let new_points = random_points_with_seed(num, seed);
                            for point in new_points.into_iter() {
                                app.tree.insert(point);
                            }
                            update_buffers(&mut app, draw_tree_nodes);
                            dirty = true;
                        },
                        _ => (),
                    }
                },
                Event::MouseInput(ElementState::Pressed, MouseButton::Left) => {
                    app.tree.insert(last_point.clone());
                    update_buffers(&mut app, draw_tree_nodes);
                    dirty = true;
                },                    
                Event::MouseMoved(x, y) => {
                    let (w, h) = app.display.get_framebuffer_dimensions();
                    // Transform x, y into the range [-1 , 1]
                    let y = h as i32 - y;
                    let x = (x as f32 / w as f32) * 2. - 1.;
                    let y = (y as f32 / h as f32) * 2. - 1.;
                    last_point = Vector2::new(x, y);
                    let color = Vector3::new(1.0, 0.0, 0.0);
                    const LOOKUP_RADIUS2: f32 = 0.2 * 0.2;
                    const N: usize = 10;
                    let mut vertices = Vec::new();
                    let mut points = Vec::new();
                    match lookup_mode {
                        LookupMode::Nearest => {
                            points.extend(
                                app.tree.nearest_neighbor(&last_point).iter().cloned());
                        },
                        LookupMode::InCircle => {
                            points.extend(app.tree.lookup_in_circle(
                                &last_point, &LOOKUP_RADIUS2).iter().cloned());
                        },
                        LookupMode::NearestN => {
                            points.extend(app.tree.nearest_n_neighbors(
                                &last_point, N));
                        },
                        LookupMode::CloseN => {
                            points.extend(app.tree.close_neighbor(&last_point).iter().cloned());
                        },
                    }
                    for point in points.iter().cloned() {
                        push_cross(&mut vertices, point, &color);
                    }

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

fn update_buffers(app: &mut ExampleApplication, draw_tree_nodes: bool) {
    let (vertices, edges) = get_tree_edges(&app.tree);
    if draw_tree_nodes {
        app.edges_buffer = VertexBuffer::new(&app.display, &edges).unwrap();
    } else {
        app.edges_buffer = VertexBuffer::new(&app.display, &[]).unwrap();
    }
    app.vertices_buffer = VertexBuffer::new(&app.display, &vertices).unwrap();
}


fn print_help() {
    println!("H - print this help dialog");
    println!("F - toggle draw tree nodes");
    println!("M - change lookup mode");
    println!("A - add 10 random points.");
    println!("B - add 100 random points.");
}

pub fn random_points_in_range<S: SpadeNum + Rand + SampleRange + BaseNum>(range: S, size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range.clone(), range.clone());
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points
}

pub fn random_points_with_seed<S: SpadeNum + BaseFloat + Rand + SampleRange>(size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    random_points_in_range(S::one(), size, seed)
}
