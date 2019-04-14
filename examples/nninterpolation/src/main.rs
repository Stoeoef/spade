// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

/*
 * This example showcases spade's interpolation features.
 * Press h for help.
 * See ./interpolation.rs for interpolation related code and
 * ./delaunay_creation.rs for code related to the generation
 * of the randomized triangulation.
 *
 * *Note*: This demo uses kiss3d which uses an old version of
 * nalgebra. This nalgebra version is incompatible with spade, that's
 * the reason why we're using cgmath points in the delaunay triangulation
 * and nalgebra points for the rendering in kiss3d. Once kiss3d updates,
 * only nalgebra will be used.
 */
#![warn(clippy::all)]
extern crate cgmath;
extern crate kiss3d;
extern crate nalgebra;
extern crate noise;
extern crate rand;
extern crate spade;

mod constants;
mod delaunay_creation;
mod interpolation;

use nalgebra as na;

use cgmath as cg;
use cgmath::EuclideanSpace;

use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::light::Light;
use kiss3d::resource::Mesh;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;

use std::cell::RefCell;
use std::rc::Rc;

use crate::delaunay_creation::Delaunay;
use crate::interpolation::interpolation_methods::{
    BarycentricInterpolation, FarinC1Interpolation, NaturalNeighborInterpolation,
    SibsonC1Interpolation,
};
use crate::interpolation::{Grid, InterpolationMethod};

struct InterpolationRenderData {
    edges: Vec<(na::Point3<f32>, na::Point3<f32>)>,
    mesh: Rc<RefCell<Mesh>>,
    title: &'static str,
}

pub fn cg_vec_to_na(vec: cg::Vector3<f64>) -> na::Point3<f32> {
    na::Point3::new(vec.x as f32, vec.y as f32, vec.z as f32)
}

impl InterpolationRenderData {
    fn new<I: InterpolationMethod>(delaunay: &Delaunay) -> InterpolationRenderData {
        let grid = Grid::<I>::from_delaunay_interpolation(delaunay);
        let (vertices, indices) = grid.get_triangles();
        let mesh = Mesh::new(vertices, indices, None, None, false);
        InterpolationRenderData {
            edges: grid.get_edges(),
            mesh: Rc::new(RefCell::new(mesh)),
            title: I::title(),
        }
    }
}

fn print_help() {
    println!("Interpolation Demo");
    println!("H - print this help");
    println!("N - show / hide normals");
    println!("G - switch interpolation method");
    println!("D - show / hide delaunay triangulation");
    println!("T - toggle display method of interpolated mesh");
}

#[derive(PartialEq, Eq)]
enum DelaunayVisibility {
    All,
    OnlyLines,
    None,
}

#[derive(PartialEq, Eq)]
enum GridRenderType {
    Lines,
    Polygons,
}

impl DelaunayVisibility {
    fn next(&self) -> DelaunayVisibility {
        use crate::DelaunayVisibility::*;
        match self {
            All => OnlyLines,
            OnlyLines => None,
            None => All,
        }
    }
}

impl GridRenderType {
    fn next(&self) -> GridRenderType {
        use crate::GridRenderType::*;
        match self {
            Lines => Polygons,
            Polygons => Lines,
        }
    }
}

fn main() {
    let mut window = Window::new("Delaunay Demo");
    window.set_light(Light::StickToCamera);

    print_help();

    let mut delaunay_visibility = DelaunayVisibility::All;
    let mut grid_render_type = GridRenderType::Lines;
    let mut show_normals = false;

    // Create delaunay triangulation and its mesh
    let delaunay = crate::delaunay_creation::generate_random_triangulation();
    let delaunay_mesh = create_mesh_from_triangulation(&delaunay);
    let delaunay_mesh = Rc::new(RefCell::new(delaunay_mesh));
    let mut delaunay_node = window.add_mesh(delaunay_mesh.clone(), na::Vector3::new(1.0, 1.0, 1.0));
    delaunay_node.enable_backface_culling(false);
    let delaunay_lines = extract_edges(&delaunay);

    let interpolation_meshes = [
        InterpolationRenderData::new::<BarycentricInterpolation>(&delaunay),
        InterpolationRenderData::new::<NaturalNeighborInterpolation>(&delaunay),
        InterpolationRenderData::new::<SibsonC1Interpolation>(&delaunay),
        InterpolationRenderData::new::<FarinC1Interpolation>(&delaunay),
    ];

    let mut cur_interpolation_mesh_node: Option<SceneNode> = None;
    let mut cur_interpolation_mesh_index = 0;

    let normals = get_normals(&delaunay);

    while window.render() {
        for event in window.events().iter() {
            let mut update_interpolation_mesh = false;
            match event.value {
                WindowEvent::Key(Key::H, Action::Press, _) => print_help(),
                WindowEvent::Key(Key::N, Action::Press, _) => show_normals = !show_normals,
                WindowEvent::Key(Key::G, Action::Press, _) => {
                    cur_interpolation_mesh_index += 1;
                    update_interpolation_mesh = true;
                    if cur_interpolation_mesh_index > interpolation_meshes.len() {
                        cur_interpolation_mesh_index = 0;
                    }
                    if cur_interpolation_mesh_index < interpolation_meshes.len() {
                        println!(
                            "Change interpolation method to {}",
                            interpolation_meshes[cur_interpolation_mesh_index].title
                        );
                    }
                }
                WindowEvent::Key(Key::T, Action::Press, _) => {
                    grid_render_type = grid_render_type.next();
                    update_interpolation_mesh = true;
                }
                WindowEvent::Key(Key::D, Action::Press, _) => {
                    delaunay_visibility = delaunay_visibility.next();
                    if delaunay_visibility == DelaunayVisibility::All {
                        delaunay_node = window
                            .scene_mut()
                            .add_mesh(delaunay_mesh.clone(), na::Vector3::new(1.0, 1.0, 1.));
                        delaunay_node.enable_backface_culling(false);
                    } else {
                        delaunay_node.unlink();
                    }
                }
                _ => {}
            }
            if update_interpolation_mesh {
                if let Some(mut mesh_node) = cur_interpolation_mesh_node {
                    mesh_node.unlink();
                    cur_interpolation_mesh_node = None;
                }
                if cur_interpolation_mesh_index < interpolation_meshes.len()
                    && grid_render_type == GridRenderType::Polygons
                {
                    let mut new_node = window.scene_mut().add_mesh(
                        interpolation_meshes[cur_interpolation_mesh_index]
                            .mesh
                            .clone(),
                        na::Vector3::new(1.0, 1.0, 1.0),
                    );
                    new_node.enable_backface_culling(false);
                    cur_interpolation_mesh_node = Some(new_node);
                }
            }
        }

        if delaunay_visibility == DelaunayVisibility::All
            || delaunay_visibility == DelaunayVisibility::OnlyLines
        {
            let color = na::Point3::new(0.8, 0.5, 0.2);
            for &(from, to) in &delaunay_lines {
                window.draw_line(&from, &to, &color);
            }
        }

        if grid_render_type == GridRenderType::Lines {
            if let Some(mesh) = interpolation_meshes.get(cur_interpolation_mesh_index) {
                let color = na::Point3::new(0.5, 0.8, 0.2);
                for &(from, to) in &mesh.edges {
                    window.draw_line(&from, &to, &color);
                }
            }
        }

        if show_normals {
            let color = na::Point3::new(0.5, 0.5, 1.0);
            for &(from, to) in &normals {
                window.draw_line(&from, &to, &color);
            }
        }
    }
}

fn get_normals(delaunay: &Delaunay) -> Vec<(na::Point3<f32>, na::Point3<f32>)> {
    let mut result = Vec::new();
    for v in delaunay.vertices() {
        let n = v.normal;
        let p = v.position_3d();
        result.push((cg_vec_to_na(p.to_vec()), cg_vec_to_na(p.to_vec() - n * 0.3)));
    }
    result
}

fn extract_edges(delaunay: &Delaunay) -> Vec<(na::Point3<f32>, na::Point3<f32>)> {
    let offset = cg::Vector3::new(0., 0., -0.01);
    let mut lines = Vec::new();
    for edge in delaunay.edges() {
        let from_pos = cg_vec_to_na(edge.from().position_3d().to_vec() + offset);
        let to_pos = cg_vec_to_na(edge.to().position_3d().to_vec() + offset);
        lines.push((from_pos, to_pos));
    }
    lines
}

fn create_mesh_from_triangulation(delaunay: &Delaunay) -> Mesh {
    let mut coords = Vec::new();
    let mut faces = Vec::new();
    for vertex in delaunay.vertices() {
        coords.push(cg_vec_to_na(vertex.position_3d().to_vec()));
    }
    for triangle in delaunay.triangles() {
        let triangle = triangle.as_triangle();
        let h0 = triangle[0].fix();
        let h1 = triangle[1].fix();
        let h2 = triangle[2].fix();
        faces.push(na::Point3::new(h0 as u16, h1 as u16, h2 as u16));
    }
    Mesh::new(coords, faces, None, None, false)
}
