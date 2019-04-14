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
use glium::glutin::VirtualKeyCode;
use glium::glutin::{ElementState, Event, MouseButton, WindowEvent};
use rand::distributions::Standard;
use rand::Rng;
use spade::delaunay::ConstrainedDelaunayTriangulation;
use spade::kernels::FloatKernel;

pub type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>, FloatKernel>;

pub fn run() {
    let mut events_loop = glium::glutin::EventsLoop::new();

    let window = ::glium::glutin::WindowBuilder::new()
        .with_dimensions(glutin::dpi::LogicalSize::new(800f64, 800f64))
        .with_title("CDT Demo".to_string());

    let context = glium::glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();
    let mut cdt = Cdt::new();

    let mut render_data = RenderData::new(&display);

    let mut last_point = Point2::new(0., 0.);
    let mut last_handle = None;

    println!("CDT Demo");
    print_help();
    let mut dirty = false;
    events_loop.run_forever(|event| {
        if let Event::WindowEvent { event, .. } = event {
            match event {
                WindowEvent::Refresh => render_data.draw(&display),
                WindowEvent::CloseRequested => return glutin::ControlFlow::Break,
                WindowEvent::MouseInput { state, button, .. }
                    if state == ElementState::Pressed && button == MouseButton::Left =>
                {
                    cdt.insert(last_point);
                    render_data.update_cdt_buffers(&display, &cdt);
                    dirty = true;
                }
                WindowEvent::MouseInput { state, button, .. }
                    if state == ElementState::Pressed && button == MouseButton::Right =>
                {
                    let nn = cdt.nearest_neighbor(&last_point).map(|p| p.fix());
                    if let Some(handle) = nn {
                        if let Some(last) = last_handle {
                            if cdt.can_add_constraint(last, handle) {
                                cdt.add_constraint(last, handle);
                                render_data.update_cdt_buffers(&display, &cdt);
                            }
                            last_handle = None;
                            render_data.update_selection_lines(&display, &[]);
                            dirty = true;
                        } else {
                            last_handle = Some(handle);
                        }
                    }
                }
                WindowEvent::CursorMoved { position, .. } => {
                    let (w, h) = display.get_framebuffer_dimensions();
                    // Transform x, y into the range [-1 , 1]
                    let x = position.x as i32;
                    let y = position.y as i32;
                    let y = h as i32 - y;
                    let x = (f64::from(x) / f64::from(w)) * 2. - 1.;
                    let y = (f64::from(y) / f64::from(h)) * 2. - 1.;
                    last_point = Point2::new(x, y);
                    let selection = get_selected_vertices(&cdt, last_point);
                    render_data.update_selection(&display, &selection);
                    if let Some(last_handle) = last_handle {
                        let highlight_line = vec![*cdt.vertex(last_handle), last_point];
                        render_data.update_selection_lines(&display, &highlight_line);
                    }
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
                                cdt.insert(point);
                            }
                            render_data.update_cdt_buffers(&display, &cdt);
                            dirty = true;
                        }
                        Some(VirtualKeyCode::D) => {
                            let nn = cdt.nearest_neighbor(&last_point).map(|p| p.fix());
                            if let Some(handle) = nn {
                                cdt.remove(handle);
                                render_data.update_cdt_buffers(&display, &cdt);
                                let selection = get_selected_vertices(&cdt, last_point);
                                render_data.update_selection(&display, &selection);
                                render_data.update_selection_lines(&display, &[]);
                                last_handle = None;
                                dirty = true;
                            }
                        }
                        _ => (),
                    }
                }
                _ => (),
            }
        }
        if dirty {
            render_data.draw(&display);
        }
        glutin::ControlFlow::Continue
    });
}

#[allow(clippy::map_clone)]
fn get_selected_vertices(cdt: &Cdt, point: Point2<f64>) -> Vec<Point2<f64>> {
    let mut points = Vec::new();
    points.extend(cdt.nearest_neighbor(&point).map(|p| *p));
    points
}

fn print_help() {
    println!("H - print this help dialog");
    println!("A - add 10 random points.");
    println!("B - add 100 random points.");
    println!("D - delete closest point.");
    println!("--------------------------");
    println!("Left click: Add single point.");
    println!("Right click: start / end adding a constraint.");
    println!();
}
