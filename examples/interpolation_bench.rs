// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

extern crate cgmath;
extern crate num;
extern crate rand;
extern crate spade;

use cgmath::{BaseNum, Point2};
use rand::distributions::uniform::SampleUniform;
use rand::{Rng, SeedableRng};
use rand_hc::Hc128Rng;
use spade::delaunay::*;
use spade::kernels::*;
use spade::{HasPosition, SpadeNum};

type Delaunay = DelaunayTriangulation<PointWithHeightAndGradient, FloatKernel, DelaunayWalkLocate>;
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

const SEED_2: &[u8; 32] = b"\xfd\x6f\x5a\xfa\x31\xed\xb4\x7e\xa5\xff\x69\x66\xdc\x73\x0a\x83\
    \xee\x98\x68\x9a\x1b\xbe\x41\x38\x15\x6e\x3b\xbf\x61\x0b\x3d\x60";

fn main() {
    let seed = b"\x2e\x4a\x77\xea\x60\x68\xcc\xb3\x42\x58\x10\x7c\x5a\xfe\xbb\x7e\
    \x40\x4b\x77\x4a\x00\xb9\x52\xa4\x41\x37\xe7\x80\xba\x78\x3c\xa1";

    let vertices =
        random_points_with_seed_range_and_origin::<f64>(20.0, Point2::new(0., 0.), SIZE, seed);

    let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
    for (index, p) in vertices.iter().enumerate() {
        delaunay.insert(PointWithHeightAndGradient(
            *p,
            (index % 10) as f64,
            Point2::new(0.0, 0.0),
        ));
    }

    delaunay.estimate_gradients(
        &(|v: &PointWithHeightAndGradient| v.1),
        &(|v: &mut PointWithHeightAndGradient, g| v.2 = g),
    );

    let query_vertices = random_points_with_seed_range_and_origin::<f64>(
        10.0,
        Point2::new(0., 0.),
        NUM_QUERY_POINTS,
        SEED_2,
    );

    println!("Starting interpolation benchmark...");
    bench(
        &query_vertices,
        &delaunay,
        "Barycentric interpolation",
        |d, p| d.barycentric_interpolation(p, |p| p.1).unwrap(),
    );
    bench(&query_vertices, &delaunay, "NN interpolation", |d, p| {
        d.nn_interpolation(p, |p| p.1).unwrap()
    });
    bench(
        &query_vertices,
        &delaunay,
        "Sibson's C1 interpolation",
        |d, p| {
            d.nn_interpolation_c1_sibson(p, 2.0, |p| p.1, |_, p| p.2)
                .unwrap()
        },
    );
    bench(
        &query_vertices,
        &delaunay,
        "Farins's C1 interpolation",
        |d, p| d.nn_interpolation_c1_farin(p, |p| p.1, |_, p| p.2).unwrap(),
    );

    println!("All done!");
}

fn bench<F>(query_vertices: &Vec<Point2<f64>>, delaunay: &Delaunay, title: &str, f: F)
where
    F: Fn(&Delaunay, &Point2<f64>) -> f64,
{
    let instant = ::std::time::Instant::now();
    for p in query_vertices {
        for _ in 0..NUM_ITERATIONS {
            let height = f(delaunay, p);
            // This assertion should prevent most clever compiler optimizations
            assert!(!height.is_nan());
        }
    }
    let elapsed = instant.elapsed();
    let millis = 1000 * elapsed.as_secs() + elapsed.subsec_nanos() as u64 / 1000000;
    println!("{} benchmark: {:?}ms", title, millis);
}

fn random_points_with_seed_range_and_origin<S: SpadeNum + BaseNum + Copy + SampleUniform>(
    range: S,
    origin: Point2<S>,
    size: usize,
    seed: &[u8; 32],
) -> Vec<Point2<S>> {
    let mut rng = Hc128Rng::from_seed(seed.clone());
    let mut points = Vec::new();
    for _ in 0..size {
        let x = rng.gen_range(-range..range) + origin.x;
        let y = rng.gen_range(-range..range) + origin.y;
        points.push(Point2::new(x, y));
    }
    points
}
