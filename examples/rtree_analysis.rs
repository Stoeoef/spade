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
use num::zero;
use rand::distributions::uniform::SampleUniform;
use rand::{Rng, SeedableRng};
use rand_hc::Hc128Rng;
use spade::rtree::RTree;
use spade::SpadeNum;
use std::fs::File;
use std::io::{stdout, Write};
use std::path::Path;
use std::time::{Duration, Instant};

fn main() {
    run_compare_operations_bench();
}

#[inline(never)]
fn blackbox<T: ?Sized>(_: &T) {}

fn measure<F, T>(result: &mut Vec<u32>, points: &[Point2<f32>], mut operation: F)
where
    F: FnMut(Point2<f32>) -> T,
{
    let now = Instant::now();
    for point in points {
        blackbox(&operation(*point));
    }
    let elapsed = now.elapsed();
    let ns = duration_ns(elapsed) as u32;
    result.push(ns / points.len() as u32);
}

fn duration_ns(duration: Duration) -> u32 {
    duration.as_secs() as u32 * 1_000_000_000 + duration.subsec_nanos()
}

fn run_compare_operations_bench() {
    const MAX_VERTICES: usize = 4000000;
    const ITERATIONS: usize = 40000;
    const NUM_STEPS: usize = 100;
    const CHUNK_SIZE: usize = MAX_VERTICES / NUM_STEPS;

    const SEED: &[u8; 32] = b"\xb0\xff\x53\x3e\x91\xa5\xd7\x15\xa5\x77\xb5\x0f\xa6\x4a\x07\xfa\
                         \x90\x00\x41\x78\x22\x8f\x12\x2f\x52\xc4\xd4\x0d\x3d\xd9\x11\x54";
    const SEED2: &[u8; 32] = b"\x2c\x35\x6a\x07\x06\x30\x6a\x36\x88\x0a\x38\xa4\x35\xe5\xeb\xcb\
        \xbe\x40\x8a\x94\xfc\x84\xe7\x84\x2a\x54\x2c\xbf\x84\x5a\xc7\x2c";

    let vertices = random_points_with_seed::<f32>(MAX_VERTICES, SEED);
    let query_points = random_points_with_seed(ITERATIONS, SEED2);
    // let vertices = random_walk_with_seed::<f32>(1.0, MAX_VERTICES, [3, 1, 4, 1]);
    // let query_points = random_walk_with_seed::<f32>(1.0, ITERATIONS, [203, 3013, 9083, 33156]);

    let mut result_file = File::create(&Path::new("rtree_compare_operations.dat")).unwrap();
    println!("Running benchmark...");

    let mut tree = RTree::new();

    let mut insert_times = Vec::new();
    let mut nearest_neighbor_times = Vec::new();
    let mut unsuccsessful_lookup_times = Vec::new();
    let mut succsessful_lookup_times = Vec::new();
    let mut nn_iterator_times = Vec::new();

    for chunk in vertices.chunks(CHUNK_SIZE) {
        print!(".");
        stdout().flush().unwrap();

        measure(&mut insert_times, chunk, |point| tree.insert(point));
        measure(&mut nearest_neighbor_times, &query_points, |point| {
            tree.nearest_neighbor(&point)
        });
        measure(&mut unsuccsessful_lookup_times, &query_points, |point| {
            tree.lookup(&point)
        });
        measure(&mut succsessful_lookup_times, chunk, |point| {
            tree.lookup(&point)
        });
        measure(&mut nn_iterator_times, &query_points, |point| {
            tree.nearest_neighbor_iterator(&point).skip(10).next()
        });
    }

    // Print all measurements to a file
    let mut print_measurements = |description: &str, measurements: &Vec<u32>| {
        write!(result_file, "\"{}\"\n", description).unwrap();
        for (index, time) in measurements.iter().enumerate() {
            let size = index * CHUNK_SIZE;
            write!(result_file, "{} {}\n", size, time).unwrap();
        }
        write!(result_file, "\n\n").unwrap();
    };

    print_measurements("insert", &insert_times);
    print_measurements("nearest_neighbor", &nearest_neighbor_times);
    print_measurements("successful lookup", &succsessful_lookup_times);
    print_measurements("unsuccessful lookup", &unsuccsessful_lookup_times);
    print_measurements("nn_iterator_times", &nn_iterator_times);

    println!("Done!");
}

pub fn random_points_with_seed<S: SpadeNum + BaseNum + Copy + SampleUniform>(
    size: usize,
    seed: &[u8; 32],
) -> Vec<Point2<S>> {
    let mut rng = Hc128Rng::from_seed(seed.clone());
    let mut points = Vec::new();
    for _ in 0..size {
        let x = rng.gen_range(-S::one()..S::one());
        let y = rng.gen_range(-S::one()..S::one());
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_walk_with_seed<S: SpadeNum + SampleUniform + BaseNum>(
    step: S,
    size: usize,
    seed: &[u8; 32],
) -> Vec<Point2<S>> {
    let mut rng = Hc128Rng::from_seed(seed.clone());
    let mut points = Vec::new();
    let mut last = Point2::new(zero(), zero());
    for _ in 0..size {
        let x = rng.gen_range(-step..step);
        let y = rng.gen_range(-step..step);
        last = Point2::new(last.x + x, last.y + y);
        points.push(last);
    }
    points
}
