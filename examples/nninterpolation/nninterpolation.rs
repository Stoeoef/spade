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
extern crate glfw;

use nalgebra::{Vector2, Vector3, Point3, Cast, Repeat};
use kiss3d::window::Window;
use kiss3d::light::Light;
use kiss3d::resource::Mesh;

use rand::Rng;
use rand::distributions::{IndependentSample, Range};
use spade::{DelaunayTriangulation, HasPosition, TrivialKernel};
use std::rc::Rc;
use std::cell::RefCell;

const SAMPLE_REGION: f64 = 3.5;
const FREQUENCY: f64 = 1.;
const NUM_POINTS: usize = 120;
const MAX_HEIGHT: f64 = 1.5;
const GRID_SUBDIVISIONS: usize = 300;

fn print_help() {
    println!("Interpolation Demo");
    println!("H - print this help");
    println!("N - show / hide normals");
    println!("G - toggle grid");
    println!("D - show / hide delaunay triangulation");
}

#[derive(PartialEq, Eq)]
enum DelaunayVisibility {
    All,
    OnlyLines,
    None,
}

#[derive(PartialEq, Eq)]
enum GridVisibility {
    ShowC0,
    ShowC1,
    ShowFarin,
    None,
}

impl GridVisibility {
    fn next(&self) -> GridVisibility {
        use GridVisibility::*;
        match self {
            &ShowC0 => ShowC1,
            &ShowC1 => ShowFarin,
            &ShowFarin => None,
            &None => ShowC0,
        }
    }

    fn print_description(&self) {
        use GridVisibility::*;
        match self {
            &ShowC0 => println!("Displaying c0 interpolant"),
            &ShowC1 => println!("Displaying sibson's c1 interpolant"),
            &ShowFarin => println!("Displaying farin's c1 interpolant"),
            _ => { },
        }
    }
}
    

impl DelaunayVisibility {
    fn next(&self) -> DelaunayVisibility {
        use DelaunayVisibility::*;
        match self {
            &All => OnlyLines,
            &OnlyLines => None,
            &None => All,
        }
    }
}

struct VectorWithHeight {
    point: Vector2<f64>,
    height: f64,
    gradient: Vector2<f64>,
    normal: Vector3<f64>,
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

    fn new(point: Vector2<f64>, height: f64) -> VectorWithHeight {
        VectorWithHeight {
            point: point,
            height: height,
            gradient: Vector2::new(0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, 0.0),
        }
    }
}


fn main() {
    use glfw::WindowEvent;
    use glfw::Action::Press;
    use glfw::Key;
    let mut window = Window::new("Delaunay Demo");

    window.set_light(Light::StickToCamera);

    let mut grid_visibility = GridVisibility::ShowC0;
    let mut delaunay_visibility = DelaunayVisibility::All;
    let mut show_normals = false;

    let (mesh, delaunay) = generate_random_mesh();
    let mesh = Rc::new(RefCell::new(mesh));
    let mut node = window.add_mesh(mesh.clone(), Repeat::repeat(1.0));
    node.enable_backface_culling(false);

    let lines = extract_edges(&delaunay);
    let grid_c1 = get_interpolated_grid(&delaunay, |d, p| d.nn_interpolation_c1_sibson(
        p, |v| v.height, |v| v.gradient).unwrap());

    let grid_c0 = get_interpolated_grid(&delaunay, |d, p| d.nn_interpolation(
        p, |v| v.height).unwrap());

    let grid_farin = get_interpolated_grid(
        &delaunay, |d, p| d.nn_interpolation_c1_farin(p, |v| v.height, |v| v.gradient).unwrap());


    let normals = get_normals(&delaunay);

    while window.render() {
        for event in window.events().iter() {
            match event.value {
                WindowEvent::Key(Key::H, _, Press, _) => print_help(),
                WindowEvent::Key(Key::N, _, Press, _) => show_normals = !show_normals,
                WindowEvent::Key(Key::G, _, Press, _) => {
                    grid_visibility = grid_visibility.next();
                    grid_visibility.print_description();
                }
                WindowEvent::Key(Key::D, _, Press, _) => {
                    delaunay_visibility = delaunay_visibility.next();
                    if delaunay_visibility == DelaunayVisibility::All {
                        node = window.scene_mut().add_mesh(mesh.clone(), Repeat::repeat(1.0));
                        node.enable_backface_culling(false);
                    } else {
                        node.unlink();
                    }

                }
                _ => { },
            }
        }

        if delaunay_visibility == DelaunayVisibility::All 
            || delaunay_visibility == DelaunayVisibility::OnlyLines
        {
            let color = Point3::new(0.8, 0.5, 0.2);
            for &(from, to) in &lines {
                window.draw_line(&from.to_point(), &to.to_point(), &color);
            }
        }

        let color = Point3::new(0.5, 0.8, 0.2);
        match grid_visibility {
            GridVisibility::ShowC0 => {
                for &(from, to) in &grid_c0 {
                    window.draw_line(&from.to_point(), &to.to_point(), &color);
                }
            },
            GridVisibility::ShowC1 => {
                for &(from, to) in &grid_c1 {
                    window.draw_line(&from.to_point(), &to.to_point(), &color);
                }
            },
            GridVisibility::ShowFarin => {
                for &(from, to) in &grid_farin {
                    window.draw_line(&from.to_point(), &to.to_point(), &color);
                }
            }
            _ => { },
        }

        if show_normals {
            let color = Point3::new(0.5, 0.5, 1.0);
            for &(from, to) in &normals {
                window.draw_line(&from.to_point(), &to.to_point(), &color);
            }
        }
    }
}

fn get_normals(delaunay: &DelaunayTriangulation<VectorWithHeight, TrivialKernel>)
    -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let mut result = Vec::new();
    for v in delaunay.vertices() {
        let n = v.normal;
        let p = v.position_3d();
        result.push((Cast::from(p), Cast::from(p - n * 0.3)));
    }
    result
}

fn get_interpolated_grid<F>(
    delaunay: &DelaunayTriangulation<VectorWithHeight, TrivialKernel>, f: F) 
    -> Vec<(Vector3<f32>, Vector3<f32>)> 
    where F: Fn(&DelaunayTriangulation<VectorWithHeight, TrivialKernel>, &Vector2<f64>) -> f64
{
    let mut result = Vec::new();
    const GRID_SIZE: f64 = SAMPLE_REGION * 1.05;
    const SCALE: f64 = 2.0 * GRID_SIZE / (GRID_SUBDIVISIONS as f64);
    let grid_offset = Vector2::repeat(GRID_SIZE);
    let transform = |v: Vector2<f64>| v * SCALE - grid_offset;

    const OFFSET: f64 = -0.01;
    let mut values = [[0.0; GRID_SUBDIVISIONS + 1]; GRID_SUBDIVISIONS + 1];
    for x in 0 .. GRID_SUBDIVISIONS + 1 {
        for y in 0 .. GRID_SUBDIVISIONS + 1 {
            let pos = transform(Vector2::new(x as f64, y as f64));
            let value = f(delaunay, &pos);
            values[x][y] = value;
        }
    }
    for x in 0 .. GRID_SUBDIVISIONS {
        for y in 0 .. GRID_SUBDIVISIONS {
            let from_val = values[x][y] + OFFSET;
            let from_pos = transform(Vector2::new(x as f64, y as f64));
            let from = VectorWithHeight::new(from_pos, from_val);
            for &(to_x, to_y) in &[(x + 1, y), (x, y + 1)] {
                let to_val = values[to_x][to_y] + OFFSET;
                let to_pos = transform(Vector2::new(to_x as f64, to_y as f64));
                let to = VectorWithHeight::new(to_pos, to_val);
                result.push((Cast::from(from.position_3d()),
                             Cast::from(to.position_3d())));
            }
        }
    }
    result
} 

fn generate_random_mesh() -> (Mesh, DelaunayTriangulation<VectorWithHeight, TrivialKernel>) {
    let mut delaunay = generate_random_triangulation();
    delaunay.estimate_gradients(&(|v| v.height), &(|v, g| v.gradient = g));
    delaunay.estimate_normals(&(|v| v.height), &(|v: &mut VectorWithHeight, n| v.normal = n));
    (create_mesh_from_triangulation(&delaunay), delaunay)
}

fn extract_edges(delaunay: &DelaunayTriangulation<VectorWithHeight, TrivialKernel>)
                 -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let offset = Vector3::new(0., 0., -0.01);
    let mut lines = Vec::new();
    for edge in delaunay.edges() {
        let from_pos = Cast::from(edge.from_handle().position_3d() + offset);
        let to_pos = Cast::from(edge.to_handle().position_3d() + offset);
        lines.push((from_pos, to_pos));
    }
    lines
}

fn generate_random_triangulation() -> DelaunayTriangulation<VectorWithHeight, TrivialKernel> {

    let mut rng = ::rand::thread_rng();
    let seed = ::noise::Seed::new(rng.gen());
    let mut delaunay = DelaunayTriangulation::new();

    let range = Range::new(-SAMPLE_REGION, SAMPLE_REGION);
    for _ in 0 .. NUM_POINTS {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        let height = ::noise::open_simplex2(&seed, &[x * FREQUENCY, y * FREQUENCY]) * MAX_HEIGHT;
        // Try out some other height functions, like those:
        // let height = (x * x + y * y) * 0.3;
        // let height = (x * 3.).sin() + (y - 2.).exp();
        delaunay.insert(VectorWithHeight::new(Vector2::new(x, y), height));
    }
    
    delaunay
}

fn create_mesh_from_triangulation(delaunay: &DelaunayTriangulation<VectorWithHeight, TrivialKernel>) -> Mesh {
    let mut coords = Vec::new();
    let mut faces = Vec::new();
    for vertex in delaunay.vertices() {
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
