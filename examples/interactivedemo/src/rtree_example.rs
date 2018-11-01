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

use cgmath::Point2;
use glium::glutin;
use glium::glutin::{ElementState, Event, MouseButton, VirtualKeyCode, WindowEvent};
use graphics::RenderData;
use rand::distributions::Standard;
use rand::Rng;
use spade::rtree::RTree;

#[derive(Clone, Copy)]
pub enum LookupMode {
    Nearest,
    NearestN,
    InCircle,
    CloseN,
}

impl ::std::fmt::Display for LookupMode {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                &LookupMode::Nearest => "Nearest neighbor",
                &LookupMode::NearestN => "Nearest N neighbors",
                &LookupMode::InCircle => "Contained in Circle",
                &LookupMode::CloseN => "Close neighbor",
            }
        )
    }
}

pub fn run() {
    let mut events_loop = glium::glutin::EventsLoop::new();
    let window = ::glium::glutin::WindowBuilder::new()
        .with_dimensions(glutin::dpi::LogicalSize::new(800f64, 800f64))
        .with_title("Interactive Demo".to_string());

    let context = glium::glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();

    let seed = ::rand::thread_rng().sample(Standard);
    let initial_points = ::random_points_with_seed(30, &seed);
    let mut rtree = RTree::bulk_load(initial_points);

    let mut render_data = RenderData::new(&display);
    render_data.update_rtree_buffers(&display, &rtree);

    let mut last_point = Point2::new(0., 0.);
    let mut lookup_mode = LookupMode::Nearest;

    println!("Interactive Demo");
    print_help();
    let mut dirty = false;
    events_loop.run_forever(|event| {
        match event {
            Event::WindowEvent { event, .. } => {
                match event {
                    WindowEvent::Refresh => render_data.draw(&display),
                    WindowEvent::CloseRequested => return glutin::ControlFlow::Break,
                    WindowEvent::MouseInput { state, button, .. }
                        if state == ElementState::Pressed && button == MouseButton::Left =>
                    {
                        rtree.insert(last_point);
                        render_data.update_rtree_buffers(&display, &rtree);
                        dirty = true;
                    }
                    WindowEvent::MouseInput { state, button, .. }
                        if state == ElementState::Pressed && button == MouseButton::Right =>
                    {
                        let nn = rtree.nearest_neighbor(&last_point).cloned();
                        if let Some(p) = nn {
                            rtree.remove(&p);
                            render_data.update_rtree_buffers(&display, &rtree);
                            let selection = get_selected_vertices(&rtree, last_point, lookup_mode);
                            render_data.update_selection(&display, &selection);
                            dirty = true;
                        }
                    }
                    WindowEvent::CursorMoved { position, .. } => {
                        let (w, h) = display.get_framebuffer_dimensions();
                        let x = position.x as i32;
                        let y = position.y as i32;
                        // Transform x, y into the range [-1 , 1]
                        let y = h as i32 - y;
                        let x = (x as f64 / w as f64) * 2. - 1.;
                        let y = (y as f64 / h as f64) * 2. - 1.;
                        last_point = Point2::new(x, y);
                        let selection = get_selected_vertices(&rtree, last_point, lookup_mode);
                        render_data.update_selection(&display, &selection);
                        dirty = true;
                    }
                    WindowEvent::KeyboardInput { input, .. }
                        if input.state == ElementState::Pressed =>
                    {
                        match input.virtual_keycode {
                            Some(VirtualKeyCode::Escape) => return glutin::ControlFlow::Break,
                            Some(VirtualKeyCode::H) => {
                                print_help();
                            }
                            Some(VirtualKeyCode::F) => {
                                render_data.update_rtree_buffers(&display, &rtree);
                                dirty = true;
                            }
                            Some(VirtualKeyCode::M) => {
                                match lookup_mode {
                                    LookupMode::Nearest => lookup_mode = LookupMode::NearestN,
                                    LookupMode::NearestN => lookup_mode = LookupMode::InCircle,
                                    LookupMode::InCircle => lookup_mode = LookupMode::CloseN,
                                    LookupMode::CloseN => lookup_mode = LookupMode::Nearest,
                                }
                                println!("Changed lookup mode to {}", lookup_mode);
                            }
                            Some(VirtualKeyCode::A) | Some(VirtualKeyCode::B) => {
                                // Insert some random points
                                let num = if input.virtual_keycode == Some(VirtualKeyCode::A) {
                                    10usize
                                } else {
                                    100
                                };
                                let seed = ::rand::thread_rng().sample(Standard);
                                let new_points = ::random_points_with_seed(num, &seed);
                                for point in new_points.into_iter() {
                                    rtree.insert(point);
                                }
                                render_data.update_rtree_buffers(&display, &rtree);
                                dirty = true;
                            }
                            _ => (),
                        }
                    }
                    _ => (),
                }
            }
            _ => (),
        }
        if dirty {
            render_data.draw(&display);
        }
        glutin::ControlFlow::Continue
    });
}

fn get_selected_vertices(
    tree: &RTree<Point2<f64>>,
    point: Point2<f64>,
    lookup_mode: LookupMode,
) -> Vec<Point2<f64>> {
    const LOOKUP_RADIUS2: f64 = 0.2 * 0.2;
    const N: usize = 10;
    let mut points = Vec::new();
    match lookup_mode {
        LookupMode::Nearest => {
            points.extend(tree.nearest_neighbor(&point).iter().cloned());
        }
        LookupMode::InCircle => {
            points.extend(
                tree.lookup_in_circle(&point, &LOOKUP_RADIUS2)
                    .iter()
                    .cloned(),
            );
        }
        LookupMode::NearestN => {
            points.extend(tree.nearest_n_neighbors(&point, N));
        }
        LookupMode::CloseN => {
            points.extend(tree.close_neighbor(&point).iter().cloned());
        }
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
