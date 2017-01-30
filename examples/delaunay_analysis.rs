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
extern crate num;

use spade::{DelaunayTriangulation, DelaunayKernel, TwoDimensional,
            AdaptiveIntKernel, TrivialKernel, FloatKernel, SpadeNum};
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use time::Duration;
use num::Zero;
use cgmath::{Vector2, BaseNum};
use std::path::Path;
use std::fs::File;
use std::io::{Write};

fn bench<V: TwoDimensional, K: DelaunayKernel<V::Scalar>>(vs: &[V], chunk_size: usize, title: &str)
                                                   -> Vec<i64> {
    println!("{}", title);
    let mut delaunay: DelaunayTriangulation<V, V, K> = DelaunayTriangulation::new();
    let mut result = Vec::new();
    let mut sum = 0;
    for chunk in vs.chunks(chunk_size) {
        let time = Duration::span(|| {
            for vertex in chunk.iter() {
                delaunay.insert(vertex.clone());
            }
        });
        sum += time.num_nanoseconds().unwrap();
        result.push(time.num_nanoseconds().unwrap() / chunk_size as i64);
    };
    assert!(delaunay.num_vertices() > vs.len() / 2);
    println!("time / op: {:?}ns", sum / vs.len() as i64);
    result
}

fn main() {

    const SIZE: usize = 400000;
    const CHUNK_SIZE: usize = SIZE / 100;
    const DO_RANDOM_WALK: bool = true;

    let seed = [661311, 350191, 123021, 231261];
    let (vertices_f64, vertices_i64);
    if DO_RANDOM_WALK {
        vertices_f64 = random_walk_with_seed_and_origin::<f64>(0.001, SIZE, seed);
        vertices_i64 = random_walk_with_seed_and_origin::<i64>(3, SIZE, seed);
    } else {
        vertices_f64 = random_points_with_seed_range_and_origin::<f64>(
            20.0, Vector2::new(1e10, -1e10), SIZE, seed);
        vertices_i64 = random_points_in_range::<i64>(10000, SIZE, seed);
    }
    let vertices_large_range = random_points_in_range::<i64>(
        1000000000, SIZE, seed);

    let f64_time = bench::<_, TrivialKernel>(&vertices_f64, CHUNK_SIZE, "f64 benchmark");
    let i64_time = bench::<_, TrivialKernel>(&vertices_i64, CHUNK_SIZE, "i64 benchmark");
    let apt_time = bench::<_, AdaptiveIntKernel>(&vertices_large_range, CHUNK_SIZE, 
                                                 "AdaptiveIntKernel benchmark");
    let floatk_time = bench::<_, FloatKernel>(&vertices_f64, CHUNK_SIZE, "FloatKernel benchmark");


    let mut result_file = File::create(&Path::new("delaunay_analysis.dat")).unwrap();
    let mut print_measurements = |description: &str, measurements: &Vec<i64>| {
        write!(result_file, "\"{}\"\n", description).unwrap();
        for (index, time) in measurements.iter().enumerate() {
            let size = index * CHUNK_SIZE;
            write!(result_file, "{} {}\n", size, time).unwrap();
        }
        write!(result_file, "\n\n").unwrap();
    };

    print_measurements("f64", &f64_time);
    print_measurements("i64", &i64_time);
    print_measurements("Adaptive", &apt_time);
    print_measurements("FloatKernel", &floatk_time);

    println!("Done!");

}

pub fn random_points_in_range<S: SpadeNum + Rand + SampleRange + BaseNum>(range: S, size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range.clone(), range.clone());
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points
}

pub fn random_points_with_seed_range_and_origin<S: SpadeNum + Copy + Rand + SampleRange>(
    range: S, origin: Vector2<S>, size: usize, seed: [u32; 4])
    -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range, range);
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng) + origin.x;
        let y = range.ind_sample(&mut rng) + origin.y;
        points.push(Vector2::new(x, y));
    }
    points    
}

pub fn random_walk_with_seed_and_origin<S: SpadeNum + Rand + SampleRange + BaseNum>(step: S, size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let rand_range = Range::new(-step, step);
    let mut points = Vec::new();
    let mut last = Vector2::zero();
    for _ in 0 .. size {
        let x = rand_range.ind_sample(&mut rng);
        let y = rand_range.ind_sample(&mut rng);
        last = last + Vector2::new(x, y);
        points.push(last);
    }
    points
}
