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

extern crate rand;
extern crate spade;
extern crate cgmath;
extern crate time;

use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use spade::{RTree, SpadeNum, LookupStructure};
use time::Duration;
use cgmath::{Vector2};
use std::path::Path;
use std::fs::File;
use std::io::{Write, stdout};


fn main() {
    run_compare_operations_bench();
}

#[inline(never)]
fn blackbox<T: ?Sized>(_: &T) {
}

fn measure<F, T>(result: &mut Vec<i64>, points: &[Vector2<f32>], mut operation: F) 
    where F: FnMut(Vector2<f32>) -> T {
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

    let vertices = random_points_with_seed::<f32>(MAX_VERTICES, [3, 1, 4, 1]);
    let query_points = random_points_with_seed(ITERATIONS, [2000, 1443, 2448, 99]);
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

pub fn random_points_with_seed<S: SpadeNum + Copy + Rand + SampleRange>(
    size: usize, seed: [u32; 4])
    -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-S::one(), S::one());
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points    
}
