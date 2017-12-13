// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

/*
 * This example is an interactive demo showing the features of spade's delaunay
 * triangulation and R-Tree. Press h for help.
 */


use graphics::{RenderData};
use spade::rtree::{RTree};
use cgmath::{Point2};
use glium::{DisplayBuild};
use glium::glutin::{Event, ElementState, MouseButton};
use glium::glutin::VirtualKeyCode;
use rand::Rng;

#[derive(Clone, Copy)]
pub enum LookupMode {
    Nearest,
    NearestN,
    InCircle,
    CloseN,
}

impl ::std::fmt::Display for LookupMode {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}", match self {
            &LookupMode::Nearest => "Nearest neighbor",
            &LookupMode::NearestN => "Nearest N neighbors",
            &LookupMode::InCircle => "Contained in Circle",
            &LookupMode::CloseN => "Close neighbor",
        })
    }
}

pub fn run() {

    let display = ::glium::glutin::WindowBuilder::new()
        .with_dimensions(800, 800)
        .with_title("Interactive Demo".to_string())
        .build_glium()
        .unwrap();

    let seed = ::rand::thread_rng().gen();
    let initial_points = ::random_points_with_seed(10000, seed);
    let mut rtree = RTree::bulk_load(initial_points);

    let mut render_data = RenderData::new(&display);
    render_data.update_rtree_buffers(&display, &rtree);

    let mut last_point = Point2::new(0., 0.);
    let mut lookup_mode = LookupMode::Nearest;

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
                            render_data.update_rtree_buffers(&display, &rtree);
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
                            let seed = ::rand::thread_rng().gen();
                            let new_points = ::random_points_with_seed(num, seed);
                            for point in new_points.into_iter() {
                                rtree.insert(point);
                            }
                            render_data.update_rtree_buffers(&display, &rtree);
                            dirty = true;
                        },
                        _ => (),
                    }
                },
                Event::MouseInput(ElementState::Pressed, MouseButton::Left) => {
                    rtree.insert(last_point);
                    render_data.update_rtree_buffers(&display, &rtree);
                    dirty = true;
                },
                Event::MouseInput(ElementState::Pressed, MouseButton::Right) => {
                    let nn = rtree.nearest_neighbor(&last_point).cloned();
                    if let Some(p) = nn {
                        rtree.remove(&p);
                        render_data.update_rtree_buffers(&display, &rtree);
                        let selection = get_selected_vertices(&rtree, last_point, lookup_mode);
                        render_data.update_selection(&display, &selection);
                        dirty = true;
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
    println!("M - change lookup mode");
    println!("A - add 10 random points.");
    println!("B - add 100 random points.");
    println!("--------------------------");
    println!("Left click: Add single point.");
    println!("Right click: Delete closest point.");
}
