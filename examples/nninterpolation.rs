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

extern crate kiss3d;
extern crate nalgebra;
extern crate rand;
extern crate noise;
extern crate spade;

use nalgebra::{Vector2, Vector3, Point3, Cast, Repeat};
use kiss3d::window::Window;
use kiss3d::light::Light;
use kiss3d::resource::Mesh;

use rand::Rng;
use rand::distributions::{IndependentSample, Range};
use spade::{DelaunayTriangulation, HasPosition};
use std::rc::Rc;
use std::cell::RefCell;

const SAMPLE_REGION: f64 = 3.;
const FREQUENCY: f64 = 1.;
const NUM_POINTS: usize = 100;
const MAX_HEIGHT: f64 = 1.5;
const GRID_SUBDIVISIONS: usize = 50;

struct VectorWithHeight {
    point: Vector2<f64>,
    height: f64,
}

impl HasPosition for VectorWithHeight {
    type Vector = Vector2<f64>;
    fn position(&self) -> Vector2<f64> {
        self.point
    }
}

impl VectorWithHeight {
    fn position_3d(&self) -> Vector3<f64> {
        Vector3::new(self.point.x, self.point.y, self.height)
    }
}

fn main() {
    let mut window = Window::new("Delaunay Demo");

    window.set_light(Light::StickToCamera);
    
    let (mesh, delaunay) = generate_random_mesh();
    let mesh = Rc::new(RefCell::new(mesh));
    let mut node = window.add_mesh(mesh, Repeat::repeat(1.));
    node.enable_backface_culling(false);
    
    let lines = extract_edges(&delaunay);
    let grid = get_interpolated_grid(&delaunay);

    while window.render() {
        let color = Point3::new(0.8, 0.5, 0.2);
        for &(from, to) in &lines {
            window.draw_line(&from.to_point(), &to.to_point(), &color);
        }

        let color = Point3::new(0.5, 0.8, 0.2);
        for &(from, to) in &grid {
            window.draw_line(&from.to_point(), &to.to_point(), &color);
        }

    }
}

fn get_interpolated_grid(
    delaunay: &DelaunayTriangulation<VectorWithHeight>) -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let mut result = Vec::new();
    const GRID_SIZE: f64 = SAMPLE_REGION * 1.1;
    const OFFSET: f64 = -0.1;
    let start = Vector2::new(-GRID_SIZE, -GRID_SIZE);
    let step = GRID_SIZE * 2. / GRID_SUBDIVISIONS as f64;
    for x in 0 .. GRID_SUBDIVISIONS {
        for y in 0 .. GRID_SUBDIVISIONS {
            let pos = Vector2::new(x as f64, y as f64);
            let from = start + pos * step;
            for dim in 0 .. 2 {
                let mut dir = Vector2::new(0., 0.);
                dir[dim] = 1.;
                let to = start + (pos + dir) * step;
                let from_height = delaunay.nn_interpolation(from, |v| v.height);
                let to_height = delaunay.nn_interpolation(to, |v| v.height);
                if let (Some(from_height), Some(to_height)) = (from_height, to_height) {
                    let from = VectorWithHeight { point: from, height: from_height + OFFSET };
                    let to = VectorWithHeight { point: to, height: to_height + OFFSET };
                    result.push((Cast::from(from.position_3d()), Cast::from(to.position_3d())));
                }
            }
        }
    }
    result
}

fn generate_random_mesh() -> (Mesh, DelaunayTriangulation<VectorWithHeight>) {
    let delaunay = generate_random_triangulation();
    (create_mesh_from_triangulation(&delaunay), delaunay)
}

fn extract_edges(delaunay: &DelaunayTriangulation<VectorWithHeight>)
                 -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let offset = Vector3::new(0., 0., -0.01);
    let mut lines = Vec::new();
    for edge in delaunay.subdiv().edges() {
        let from_pos = Cast::from(edge.from_handle().position_3d() + offset);
        let to_pos = Cast::from(edge.to_handle().position_3d() + offset);
        lines.push((from_pos, to_pos));
    }
    lines
}

fn generate_random_triangulation() -> DelaunayTriangulation<VectorWithHeight> {

    let mut rng = ::rand::thread_rng();
    let seed = ::noise::Seed::new(rng.gen());
    let mut delaunay = DelaunayTriangulation::new();

    let range = Range::new(-SAMPLE_REGION, SAMPLE_REGION);
    for _ in 0 .. NUM_POINTS {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        let height = ::noise::open_simplex2(&seed, &[x * FREQUENCY, y * FREQUENCY]) * MAX_HEIGHT;
        delaunay.insert(VectorWithHeight { point: Vector2::new(x, y), height: height });
    }
    
    delaunay
}

fn create_mesh_from_triangulation(delaunay: &DelaunayTriangulation<VectorWithHeight>) -> Mesh {
    let mut coords = Vec::new();
    let mut faces = Vec::new();
    for vertex in delaunay.subdiv().vertices() {
        coords.push(Cast::from(vertex.position_3d().to_point()));
    }
    for triangle in delaunay.triangles() {
        let h0 = triangle.0[0].fix();
        let h1 = triangle.0[1].fix();
        let h2 = triangle.0[2].fix();
        faces.push(Cast::from(Point3::new(h0, h1, h2)));
    }
    return Mesh::new(coords, faces, None, None, false);
}
