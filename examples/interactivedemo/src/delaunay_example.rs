// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::graphics::RenderData;
use cgmath::Point2;
use glium::glutin;
use glium::glutin::{ElementState, Event, MouseButton, VirtualKeyCode, WindowEvent};
use rand::distributions::Standard;
use rand::Rng;
use spade::delaunay::DelaunayTriangulation;
use spade::kernels::FloatKernel;

pub type ExampleTriangulation = DelaunayTriangulation<Point2<f64>, FloatKernel>;

pub fn run() {
    let mut events_loop = glium::glutin::EventsLoop::new();
    let window = ::glium::glutin::WindowBuilder::new()
        .with_dimensions(glutin::dpi::LogicalSize::new(800f64, 800f64))
        .with_title("Delauany Demo".to_string());

    let context = glium::glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();

    let mut delaunay = DelaunayTriangulation::with_tree_locate();

    let mut render_data = RenderData::new(&display);

    let mut last_point = Point2::new(0., 0.);

    println!("Delaunay Demo");
    print_help();
    let mut dirty = false;
    events_loop.run_forever(|event| {
        match event {
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::Refresh => render_data.draw(&display),
                WindowEvent::CloseRequested => return glutin::ControlFlow::Break,
                WindowEvent::MouseInput { state, button, .. }
                    if state == ElementState::Pressed && button == MouseButton::Left =>
                {
                    delaunay.insert(last_point);
                    render_data.update_delaunay_buffers(&display, &delaunay);
                    dirty = true;
                }
                WindowEvent::MouseInput { state, button, .. }
                    if state == ElementState::Pressed && button == MouseButton::Right =>
                {
                    let nn = delaunay.nearest_neighbor(&last_point).map(|p| p.fix());
                    if let Some(handle) = nn {
                        delaunay.remove(handle);
                        render_data.update_delaunay_buffers(&display, &delaunay);
                        let selection = get_selected_vertices(&delaunay, last_point);
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
                    let selection = get_selected_vertices(&delaunay, last_point);
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
                        Some(VirtualKeyCode::A) | Some(VirtualKeyCode::B) => {
                            // Insert some random points
                            let num = if input.virtual_keycode == Some(VirtualKeyCode::A) {
                                10usize
                            } else {
                                100
                            };
                            let mut rng = ::rand::thread_rng();

                            let seed = rng.sample(Standard);
                            let new_points = crate::random_points_with_seed(num, &seed);
                            for point in new_points.into_iter() {
                                delaunay.insert(point);
                            }
                            render_data.update_delaunay_buffers(&display, &delaunay);
                            dirty = true;
                        }
                        _ => (),
                    }
                }
                _ => (),
            },
            _ => (),
        }
        if dirty {
            render_data.draw(&display);
        }
        glutin::ControlFlow::Continue
    });
}

fn get_selected_vertices(delaunay: &ExampleTriangulation, point: Point2<f64>) -> Vec<Point2<f64>> {
    let mut points = Vec::new();
    points.extend(delaunay.nearest_neighbor(&point).map(|p| (*p).clone()));
    points
}

fn print_help() {
    println!("H - print this help dialog");
    println!("A - add 10 random points.");
    println!("B - add 100 random points.");
    println!("--------------------------");
    println!("Left click: Add single point.");
    println!("Right click: Delete closest point.");
    println!();
}
