// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

extern crate rand;
extern crate spade;
extern crate cgmath;
extern crate num;

use spade::delaunay::*;
use spade::kernels::*;
use spade::{SpadeNum, HasPosition};
use rand::{XorShiftRng, SeedableRng};
use rand::distributions::{Range, Distribution};
use rand::distributions::range::SampleRange;
use cgmath::{Point2, BaseNum};

type Delaunay = DelaunayTriangulation<PointWithHeightAndGradient, 
                                      FloatKernel, DelaunayWalkLocate>;
struct PointWithHeightAndGradient(Point2<f64>, f64, Point2<f64>);

impl HasPosition for PointWithHeightAndGradient {
    type Point = Point2<f64>;
    fn position(&self) -> Point2<f64> {
        self.0
    }
}

const SIZE: usize = 6000;
const NUM_ITERATIONS: usize = 10000;
const NUM_QUERY_POINTS: usize = 400;

fn main() {

    let vertices = random_points_with_seed_range_and_origin::<f64>(
        20.0, Point2::new(0., 0.), SIZE, b"Look out through");

    let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
    for (index, p) in vertices.iter().enumerate() {
        delaunay.insert(PointWithHeightAndGradient(*p, (index % 10) as f64,
                                                   Point2::new(0.0, 0.0)));
    }

    delaunay.estimate_gradients(&(|v: &PointWithHeightAndGradient| v.1),
                                &(|v: &mut PointWithHeightAndGradient, g| v.2 = g));

    let query_vertices = random_points_with_seed_range_and_origin::<f64>(
        10.0, Point2::new(0., 0.), NUM_QUERY_POINTS, b"the window. Clin");
    
    println!("Starting interpolation benchmark...");
    bench(&query_vertices, &delaunay, "Barycentric interpolation", 
          |d, p| d.barycentric_interpolation(p, |p| p.1).unwrap());
    bench(&query_vertices, &delaunay, "NN interpolation", 
          |d, p| d.nn_interpolation(p, |p| p.1).unwrap());
    bench(&query_vertices, &delaunay, "Sibson's C1 interpolation", 
          |d, p| d.nn_interpolation_c1_sibson(p, 2.0, |p| p.1, |_, p| p.2).unwrap());
    bench(&query_vertices, &delaunay, "Farins's C1 interpolation", 
          |d, p| d.nn_interpolation_c1_farin(p, |p| p.1, |_, p| p.2).unwrap());


    println!("All done!");
}

fn bench<F>(query_vertices: &Vec<Point2<f64>>, delaunay: &Delaunay, title: &str, f: F) where
    F: Fn(&Delaunay, &Point2<f64>) -> f64
{
    let instant = ::std::time::Instant::now();
    for p in query_vertices {
        for _ in 0 .. NUM_ITERATIONS {
            let height = f(delaunay, p);
            // This assertion should prevent most clever compiler optimizations
            assert!(!height.is_nan());
        }
    }
    let elapsed = instant.elapsed();
    let millis = 1000 * elapsed.as_secs() + elapsed.subsec_nanos() as u64 / 1000000;
    println!("{} benchmark: {:?}ms", title, millis);
}

fn random_points_with_seed_range_and_origin<S: SpadeNum + BaseNum + Copy + SampleRange>(
    range: S, origin: Point2<S>, size: usize, seed: &[u8; 16])
    -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed.clone());
    let range = Range::new(-range, range);
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.sample(&mut rng) + origin.x;
        let y = range.sample(&mut rng) + origin.y;
        points.push(Point2::new(x, y));
    }
    points    
}
