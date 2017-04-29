// Copyright 2017 The Spade Developers.
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

/*
 * This example is an interactive demo showing the features of spade's delaunay
 * triangulation and R-Tree. Press h for help.
 */

extern crate spade;
extern crate rand;
extern crate cgmath;
#[macro_use]
extern crate glium;

mod graphics;
use graphics::{RenderData};
use spade::rtree::{RTree};
use spade::{SpadeNum};
use spade::delaunay::{ConstrainedDelaunayTriangulation, Locateable};
use spade::kernels::{FloatKernel};
use cgmath::{Point2, BaseFloat, BaseNum};
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use glium::{DisplayBuild};
use glium::glutin::{Event, ElementState, MouseButton};
use glium::glutin::VirtualKeyCode;

// type ExampleTriangulation = DelaunayTriangulation<Point2<f64>, FloatKernel>;
type ExampleTriangulation = ConstrainedDelaunayTriangulation<Point2<f64>, FloatKernel>;

#[derive(Clone, Copy)]
pub enum LookupMode {
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

fn main() {

    let display = glium::glutin::WindowBuilder::new()
        .with_dimensions(800, 800)
        .with_title("Interactive Demo".to_string())
        .build_glium()
        .unwrap();

    let mut delaunay = ExampleTriangulation::new();
    let mut segment_start_maybe = None;
    
    let mut rtree = RTree::new();

    let mut render_data = RenderData::new(&display);

    let mut last_point = Point2::new(0., 0.);
    let mut lookup_mode = LookupMode::Nearest;
    let mut draw_tree_nodes = false;

    println!("Interactive Demo");
    print_help();
    loop {
        let events: Vec<_> = display.poll_events().collect();

        let mut dirty = false;
        for event in events.into_iter() {
            match event {
                Event::Refresh => render_data.draw(&display),
                Event::Closed => return,
                Event::KeyboardInput(ElementState::Pressed, _, Some(key)) => {
                    match key {
                        VirtualKeyCode::Escape => return,
                        VirtualKeyCode::H => {
                            print_help();
                        },
                        VirtualKeyCode::F => {
                            draw_tree_nodes = !draw_tree_nodes;
                            render_data.update_buffers(&display, &rtree, &delaunay, draw_tree_nodes);
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
                            let s = rtree.size() as u32;
                            let seed = [s, s * 35, s + 124, s * 91];
                            let new_points = random_points_with_seed(num, seed);
                            for point in new_points.into_iter() {
                                delaunay.insert(point);
                                rtree.insert(point);
                            }
                            render_data.update_buffers(&display, &rtree, &delaunay, draw_tree_nodes);
                            dirty = true;
                        },
                        _ => (),
                    }
                },
                Event::MouseInput(ElementState::Pressed, MouseButton::Left) => {
                    rtree.insert(last_point);
                    delaunay.insert(last_point);
                    render_data.update_buffers(&display, &rtree, &delaunay, draw_tree_nodes);
                    dirty = true;
                },
                Event::MouseInput(ElementState::Pressed, MouseButton::Right) => {
                    let nn = rtree.nearest_neighbor(&last_point).cloned();
                    if let Some(p) = nn {
                        
                    //     rtree.remove(&p);
                        let handle = delaunay.locate_vertex_with_hint(&p, 0).unwrap().fix();
                    //     delaunay.remove(handle);
                    //     render_data.update_buffers(&display, &rtree, &delaunay, draw_tree_nodes);
                    //     let selection = get_selected_vertices(&rtree, last_point, lookup_mode);
                    //     render_data.update_selection(&display, &selection);
                        if let Some(segment_start) = segment_start_maybe {
                            dirty = true;
                            delaunay.add_constraint(segment_start, handle);
                            render_data.update_buffers(&display, &rtree, &delaunay, draw_tree_nodes);
                            segment_start_maybe = None;
                        } else {
                            segment_start_maybe = Some(handle);
                        }
                    }
                },                    
                Event::MouseMoved(x, y) => {
                    let (w, h) = display.get_framebuffer_dimensions();
                    // Transform x, y into the range [-1 , 1]
                    let y = h as i32 - y;
                    let x = (x as f64 / w as f64) * 2. - 1.;
                    let y = (y as f64 / h as f64) * 2. - 1.;
                    last_point = Point2::new(x, y);
                    let selection = get_selected_vertices(&rtree, last_point, lookup_mode);
                    render_data.update_selection(&display, &selection);
                    dirty = true;
                },
                _ => (),
            }
        }
        if dirty {
            render_data.draw(&display);
        }
    }
}

fn get_selected_vertices(tree: &RTree<Point2<f64>>, point: Point2<f64>,
                         lookup_mode: LookupMode) -> Vec<Point2<f64>> {
    
    const LOOKUP_RADIUS2: f64 = 0.2 * 0.2;
    const N: usize = 10;
    let mut points = Vec::new();
    match lookup_mode {
        LookupMode::Nearest => {
            points.extend(
                tree.nearest_neighbor(&point).iter().cloned());
        },
        LookupMode::InCircle => {
            points.extend(tree.lookup_in_circle(
                &point, &LOOKUP_RADIUS2).iter().cloned());
        },
        LookupMode::NearestN => {
            points.extend(tree.nearest_n_neighbors(
                &point, N));
        },
        LookupMode::CloseN => {
            points.extend(tree.close_neighbor(&point).iter().cloned());
        },
    }
    points
}

fn print_help() {
    println!("H - print this help dialog");
    println!("F - toggle drawing of r-tree nodes");
    println!("M - change lookup mode");
    println!("A - add 10 random points.");
    println!("B - add 100 random points.");
    println!("--------------------------");
    println!("Left click: Add single point.");
    println!("Right click: Delete closest point.");
}

pub fn random_points_in_range<S: SpadeNum + Rand + SampleRange + BaseNum>(range: S, size: usize, seed: [u32; 4]) -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range.clone(), range.clone());
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_points_with_seed<S: SpadeNum + BaseFloat + Rand + SampleRange>(size: usize, seed: [u32; 4]) -> Vec<Point2<S>> {
    random_points_in_range(S::one(), size, seed)
}
