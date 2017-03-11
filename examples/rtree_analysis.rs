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
extern crate time;
extern crate num;

use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use spade::SpadeNum;
use spade::rtree::RTree;
use time::Duration;
use cgmath::{Point2, BaseNum};
use std::path::Path;
use std::fs::File;
use std::io::{Write, stdout};
use num::zero;

fn main() {
    run_compare_operations_bench();
}

#[inline(never)]
fn blackbox<T: ?Sized>(_: &T) {
}

fn measure<F, T>(result: &mut Vec<i64>, points: &[Point2<f32>], mut operation: F) 
    where F: FnMut(Point2<f32>) -> T {
    let time = Duration::span(|| {
        for point in points {
            blackbox(&operation(*point));
        }
    }).num_nanoseconds().unwrap();
    result.push(time / points.len() as i64);
 }

fn run_compare_operations_bench() {
    const MAX_VERTICES: usize = 4000000;
    const ITERATIONS: usize = 40000;
    const NUM_STEPS: usize = 100;
    const CHUNK_SIZE: usize = MAX_VERTICES / NUM_STEPS;

    // let vertices = random_points_with_seed::<f32>(MAX_VERTICES, [3, 1, 4, 1]);
    // let query_points = random_points_with_seed(ITERATIONS, [2000, 1443, 2448, 99]);
    let vertices = random_walk_with_seed::<f32>(1.0, MAX_VERTICES, [3, 1, 4, 1]);
    let query_points = random_walk_with_seed::<f32>(1.0, ITERATIONS, [203, 3013, 9083, 33156]);

    let mut result_file = File::create(&Path::new("rtree_compare_operations.dat")).unwrap();
    println!("Running benchmark...");

    let mut insert_times = Vec::new();
    let mut nearest_neighbor_times = Vec::new();
    let mut unsuccsessful_lookup_times = Vec::new();
    let mut succsessful_lookup_times = Vec::new();

    let mut tree = RTree::new();

    for chunk in vertices.chunks(CHUNK_SIZE) {
        print!(".");
        stdout().flush().unwrap();

        measure(&mut insert_times, chunk, |point| tree.insert(point));
        measure(&mut nearest_neighbor_times, &query_points, |point| tree.nearest_neighbor(&point));
        measure(&mut unsuccsessful_lookup_times, &query_points, |point| tree.lookup(&point));
        measure(&mut succsessful_lookup_times, chunk, |point| tree.lookup(&point));
    }

    // Print all measurements to a file
    let mut print_measurements = |description: &str, measurements: &Vec<i64>| {
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

    println!("Done!");
}

pub fn random_points_with_seed<S: SpadeNum + BaseNum + Copy + Rand + SampleRange>(
    size: usize, seed: [u32; 4])
    -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-S::one(), S::one());
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points    
}

pub fn random_walk_with_seed<S: SpadeNum + Rand + SampleRange + BaseNum>(step: S, size: usize, seed: [u32; 4]) -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let rand_range = Range::new(-step, step);
    let mut points = Vec::new();
    let mut last = Point2::new(zero(), zero());
    for _ in 0 .. size {
        let x = rand_range.ind_sample(&mut rng);
        let y = rand_range.ind_sample(&mut rng);
        last = Point2::new(last.x + x, last.y + y);
        points.push(last);
    }
    points
}
