use nalgebra::{Vector2, Vector3, Point3, Cast};

use delaunay_creation::{Delaunay, VectorWithHeight};
use constants::*;

// Interpolation Methods ------------------------------
pub trait InterpolationMethod {
    fn interpolate(d: &Delaunay, point: Vector2<f64>) -> f64;
    fn title() -> &'static str;
}

pub mod interpolation_methods {
    use super::InterpolationMethod;
    use ::delaunay_creation::Delaunay;
    use nalgebra::Vector2;

    pub struct BarycentricInterpolation;

    impl InterpolationMethod for BarycentricInterpolation {
        fn interpolate(delaunay: &Delaunay, point: Vector2<f64>) -> f64 {
            delaunay.barycentric_interpolation(&point, |v| v.height).unwrap()
        }

        fn title() -> &'static str {
            "barycentric interpolation"
        }
    }

    pub struct NaturalNeighborInterpolation;

    impl InterpolationMethod for NaturalNeighborInterpolation {
        fn interpolate(delaunay: &Delaunay, point: Vector2<f64>) -> f64 {
            delaunay.nn_interpolation(&point, |v| v.height).unwrap()
        }

        fn title() -> &'static str {
            "natural neighbor interpolation"
        }
    }

    pub struct SibsonC1Interpolation;
    impl InterpolationMethod for SibsonC1Interpolation {
        fn interpolate(delaunay: &Delaunay, point: Vector2<f64>) -> f64 {
            delaunay.nn_interpolation_c1_sibson(
                &point,
                // Check out different smoothness factors
                // 0.5,
                // 2.0,
                1.0,
                // The second function defines the gradient of a point
                |v| v.height, |_, v| v.gradient).unwrap()
        }

        fn title() -> &'static str {
            "sibson's c1 interpolation"
        }
    }

    pub struct FarinC1Interpolation;
    impl InterpolationMethod for FarinC1Interpolation {
        fn interpolate(delaunay: &Delaunay, point: Vector2<f64>) -> f64 {
            delaunay.nn_interpolation_c1_farin(
                &point,
                // The second function defines the gradient of a point
                |v| v.height, |_, v| v.gradient).unwrap()
        }

        fn title() -> &'static str {
            "farin's c1 interpolation"
        }
    }
}

/*
 * Caches interpolated values on a grid and offers methods to 
 * convert these into an edge list or a vertices / indices list
 */
pub struct Grid<I: InterpolationMethod> {
    grid: [[f64; GRID_SUBDIVISIONS + 1]; GRID_SUBDIVISIONS + 1],
    __interpolation: ::std::marker::PhantomData<I>,
}

impl <I: InterpolationMethod> Grid<I> {
    // Returns a list of edges for rendering
    pub fn get_edges(&self) -> Vec<(Vector3<f32>, Vector3<f32>)> {
        let mut result = Vec::new();
        for x in 0 .. GRID_SUBDIVISIONS {
            for y in 0 .. GRID_SUBDIVISIONS {
                let from_val = self.grid[x][y] + OFFSET;
                let from_pos = Self::transform(Vector2::new(x as f64, y as f64));
                let from = VectorWithHeight::new(from_pos, from_val);
                for &(to_x, to_y) in &[(x + 1, y), (x, y + 1)] {
                    let to_val = self.grid[to_x][to_y] + OFFSET;
                    let to_pos = Self::transform(Vector2::new(to_x as f64, to_y as f64));
                    let to = VectorWithHeight::new(to_pos, to_val);
                    result.push((Cast::from(from.position_3d()),
                                 Cast::from(to.position_3d())));
                }
            }
        }
        result
    }

    // Returns a list of vertices and a list of triangle indices that form the
    // grid's mesh.
    pub fn get_triangles(&self) -> (Vec<Point3<f32>>, Vec<Point3<u32>>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for x in 0 .. GRID_SUBDIVISIONS + 1 {
            for y in 0 .. GRID_SUBDIVISIONS + 1 {
                let val = self.grid[x][y] + OFFSET;
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


    // This will do the actual interpolation and store it in the triangulation
    pub fn from_delaunay_interpolation(delaunay: &Delaunay) -> Grid<I>
    {
        let mut values = [[0.0; GRID_SUBDIVISIONS + 1]; GRID_SUBDIVISIONS + 1];
        for x in 0 .. GRID_SUBDIVISIONS + 1 {
            for y in 0 .. GRID_SUBDIVISIONS + 1 {
                let pos = Self::transform(Vector2::new(x as f64, y as f64));
                let value = I::interpolate(delaunay, pos);
                values[x][y] = value;
            }
        }
        Grid {
            grid: values,
            __interpolation: Default::default(),
        }
    }

    fn transform(v: Vector2<f64>) -> Vector2<f64> {
        v * SCALE - GRID_OFFSET
    }
}
