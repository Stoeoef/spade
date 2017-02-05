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
use spade::delaunay::{DelaunayTriangulation, FloatDelaunayTriangulation, TriangulationWalkLookup};
use spade::{HasPosition,};
use std::rc::Rc;
use std::cell::RefCell;

type Delaunay = FloatDelaunayTriangulation<
        VectorWithHeight, TriangulationWalkLookup<Vector2<f64>>>;

const SAMPLE_REGION: f64 = 3.5;
const FREQUENCY: f64 = 1.;
const NUM_POINTS: usize = 120;
const MAX_HEIGHT: f64 = 1.5;
const GRID_SUBDIVISIONS: usize = 300;

struct Grid {
    grid: [[f64; GRID_SUBDIVISIONS + 1]; GRID_SUBDIVISIONS + 1],
}

struct InterpolationRenderData {
    edges: Vec<(Vector3<f32>, Vector3<f32>)>,
    mesh: Rc<RefCell<Mesh>>,
}

impl InterpolationRenderData {
    fn from_grid(grid: &Grid) -> InterpolationRenderData {
        let (vertices, indices) = grid.get_triangles();
        let mesh = Mesh::new(vertices, indices, None, None, false);
        InterpolationRenderData {
            edges: grid.get_edges(),
            mesh: Rc::new(RefCell::new(mesh)),
        }
    }
}

impl Grid {
    fn get_edges(&self) -> Vec<(Vector3<f32>, Vector3<f32>)> {
        let mut result = Vec::new();
        for x in 0 .. GRID_SUBDIVISIONS {
            for y in 0 .. GRID_SUBDIVISIONS {
                let from_val = self.grid[x][y] + Self::offset();
                let from_pos = Self::transform(Vector2::new(x as f64, y as f64));
                let from = VectorWithHeight::new(from_pos, from_val);
                for &(to_x, to_y) in &[(x + 1, y), (x, y + 1)] {
                    let to_val = self.grid[to_x][to_y] + Self::offset();
                    let to_pos = Self::transform(Vector2::new(to_x as f64, to_y as f64));
                    let to = VectorWithHeight::new(to_pos, to_val);
                    result.push((Cast::from(from.position_3d()),
                                 Cast::from(to.position_3d())));
                }
            }
        }
        result
    }

    fn get_triangles(&self) -> (Vec<Point3<f32>>, Vec<Point3<u32>>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for x in 0 .. GRID_SUBDIVISIONS + 1 {
            for y in 0 .. GRID_SUBDIVISIONS + 1 {
                let val = self.grid[x][y] + Self::offset();
                let pos = Self::transform(Vector2::new(x as f64, y as f64));
                vertices.push(Point3::new(pos.x as f32, pos.y as f32, val as f32));
            }
        }
        for x in 0 .. GRID_SUBDIVISIONS {
            for y in 0 .. GRID_SUBDIVISIONS {
                let index = |x, y| x * (GRID_SUBDIVISIONS + 1) + y;
                let v00 = index(x, y) as u32;
                let v10 = index(x + 1, y) as u32;
                let v01 = index(x, y + 1) as u32;
                let v11 = index(x + 1, y + 1) as u32;
                indices.push(Point3::new(v00, v10, v11));
                indices.push(Point3::new(v00, v11, v01));
            }
        }
        (vertices, indices)
    }

    fn from_delaunay_interpolation<F>(delaunay: &Delaunay, f: F) -> Grid
        where F: Fn(&Delaunay, &Vector2<f64>) -> f64
    {
        let mut values = [[0.0; GRID_SUBDIVISIONS + 1]; GRID_SUBDIVISIONS + 1];
        for x in 0 .. GRID_SUBDIVISIONS + 1 {
            for y in 0 .. GRID_SUBDIVISIONS + 1 {
                let pos = Self::transform(Vector2::new(x as f64, y as f64));
                let value = f(delaunay, &pos);
                values[x][y] = value;
            }
        }
        Grid {
            grid: values,
        }
    }

    fn offset() -> f64 {
        -0.01
    }

    fn grid_size() -> f64 {
        SAMPLE_REGION * 1.05
    }

    fn scale() -> f64 {
        2.0 * Self::grid_size() / (GRID_SUBDIVISIONS as f64)
    }

    fn grid_offset() -> Vector2<f64> {
        Vector2::repeat(Self::grid_size())
    }

    fn transform(v: Vector2<f64>) -> Vector2<f64> {
        v * Self::scale() - Self::grid_offset()
    }
}

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
    ShowBarycentric,
    None,
}

#[derive(PartialEq, Eq)]
enum GridRenderType {
    Lines,
    Polygons,
}

impl GridVisibility {
    fn next(&self) -> GridVisibility {
        use GridVisibility::*;
        match self {
            &ShowC0 => ShowC1,
            &ShowC1 => ShowFarin,
            &ShowFarin => ShowBarycentric,
            &ShowBarycentric => None,
            &None => ShowC0,
        }
    }

    fn print_description(&self) {
        use GridVisibility::*;
        match self {
            &ShowC0 => println!("Displaying c0 interpolant"),
            &ShowC1 => println!("Displaying sibson's c1 interpolant"),
            &ShowFarin => println!("Displaying farin's c1 interpolant"),
            &ShowBarycentric => println!("Displaying barycentric interpolant"),
            &None => { },
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

impl GridRenderType {
    fn next(&self) -> GridRenderType {
        use GridRenderType::*;
        match self {
            &Lines => Polygons,
            &Polygons => Lines
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
    let mut grid_render_type = GridRenderType::Lines;
    let mut show_normals = false;

    let (mesh, delaunay) = generate_random_mesh();
    let mesh = Rc::new(RefCell::new(mesh));
    let mut delaunay_node = window.add_mesh(mesh.clone(), Repeat::repeat(1.0));
    delaunay_node.enable_backface_culling(false);
    let empty_mesh = Rc::new(RefCell::new(Mesh::new(vec![], vec![], None, None, false)));
    let mut grid_node = window.add_mesh(empty_mesh.clone(), Repeat::repeat(1.0));
    grid_node.unlink();

    let lines = extract_edges(&delaunay);

    let data_barycentric = InterpolationRenderData::from_grid(
        &Grid::from_delaunay_interpolation(&delaunay, |d, p| d.barycentric_interpolation(
            p, |v| v.height).unwrap()));

    let data_c1 = InterpolationRenderData::from_grid(
        &Grid::from_delaunay_interpolation(&delaunay, |d, p| d.nn_interpolation_c1_sibson(
            p, 
            // Check out different smoothness factors
            // 0.5,
            // 2.0,
            // 1.0,
            3.9,
            |v| v.height, |_, v| v.gradient).unwrap()));
    
    let data_c0 = InterpolationRenderData::from_grid(
        &Grid::from_delaunay_interpolation(&delaunay, |d, p| d.nn_interpolation(
            p, |v| v.height).unwrap()));

    let data_farin = InterpolationRenderData::from_grid(
        &Grid::from_delaunay_interpolation(
            &delaunay, |d, p| d.nn_interpolation_c1_farin(p, |v| v.height, |_, v| v.gradient).unwrap()));

    let mut cur_data = None;

    let normals = get_normals(&delaunay);

    while window.render() {
        for event in window.events().iter() {
            match event.value {
                WindowEvent::Key(Key::H, _, Press, _) => print_help(),
                WindowEvent::Key(Key::N, _, Press, _) => show_normals = !show_normals,
                WindowEvent::Key(Key::G, _, Press, _) => {
                    grid_visibility = grid_visibility.next();
                    cur_data = match grid_visibility {
                        GridVisibility::None => None,
                        GridVisibility::ShowC0 => Some(&data_c0),
                        GridVisibility::ShowC1 => Some(&data_c1),
                        GridVisibility::ShowFarin => Some(&data_farin),
                        GridVisibility::ShowBarycentric => Some(&data_barycentric),
                    };
                    if let Some(data) = cur_data {
                        grid_node.unlink();
                        grid_node = window.scene_mut().add_mesh(data.mesh.clone(), Repeat::repeat(1.0));
                        grid_node.enable_backface_culling(false);
                    } else {
                        grid_node.unlink();
                    }
                    grid_visibility.print_description();
                },
                WindowEvent::Key(Key::T, _, Press, _) => {
                    grid_render_type = grid_render_type.next();
                    if grid_render_type == GridRenderType::Polygons {
                        if let Some(data) = cur_data {
                            grid_node.unlink();
                            grid_node = window.scene_mut().add_mesh(data.mesh.clone(), Repeat::repeat(1.0));
                            grid_node.enable_backface_culling(false);
                        }
                    } else {
                        grid_node.unlink();
                    }
                },
                WindowEvent::Key(Key::D, _, Press, _) => {
                    delaunay_visibility = delaunay_visibility.next();
                    if delaunay_visibility == DelaunayVisibility::All {
                        delaunay_node = window.scene_mut().add_mesh(mesh.clone(), Repeat::repeat(1.0));
                        delaunay_node.enable_backface_culling(false);
                    } else {
                        delaunay_node.unlink();
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

        if grid_render_type == GridRenderType::Lines {
            if let Some(data) = cur_data {
                let color = Point3::new(0.5, 0.8, 0.2);
                for &(from, to) in &data.edges {
                    window.draw_line(&from.to_point(), &to.to_point(), &color);
                }
            }
        }

        if show_normals {
            let color = Point3::new(0.5, 0.5, 1.0);
            for &(from, to) in &normals {
                window.draw_line(&from.to_point(), &to.to_point(), &color);
            }
        }
    }
}

fn get_normals(delaunay: &Delaunay) -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let mut result = Vec::new();
    for v in delaunay.vertices() {
        let n = v.normal;
        let p = v.position_3d();
        result.push((Cast::from(p), Cast::from(p - n * 0.3)));
    }
    result
}

fn generate_random_mesh() -> (Mesh, Delaunay) {
    let mut delaunay = generate_random_triangulation();
    delaunay.estimate_gradients(&(|v| v.height), &(|v, g| v.gradient = g));
    delaunay.estimate_normals(&(|v| v.height), &(|v: &mut VectorWithHeight, n| v.normal = n));
    (create_mesh_from_triangulation(&delaunay), delaunay)
}

fn extract_edges(delaunay: &Delaunay)
                 -> Vec<(Vector3<f32>, Vector3<f32>)> {
    let offset = Vector3::new(0., 0., -0.01);
    let mut lines = Vec::new();
    for edge in delaunay.edges() {
        let from_pos = Cast::from(edge.from().position_3d() + offset);
        let to_pos = Cast::from(edge.to().position_3d() + offset);
        lines.push((from_pos, to_pos));
    }
    lines
}

fn generate_random_triangulation() -> Delaunay {

    let mut rng = ::rand::thread_rng();
    let seed = ::noise::Seed::new(rng.gen());
    let mut delaunay = DelaunayTriangulation::with_walk_lookup();

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

fn create_mesh_from_triangulation(delaunay: &Delaunay) -> Mesh {
    let mut coords = Vec::new();
    let mut faces = Vec::new();
    for vertex in delaunay.vertices() {
        coords.push(Cast::from(vertex.position_3d().to_point()));
    }
    for triangle in delaunay.triangles() {
        let triangle = triangle.as_triangle();
        let h0 = triangle[0].fix();
        let h1 = triangle[1].fix();
        let h2 = triangle[2].fix();
        faces.push(Cast::from(Point3::new(h0, h1, h2)));
    }
    return Mesh::new(coords, faces, None, None, false);
}
