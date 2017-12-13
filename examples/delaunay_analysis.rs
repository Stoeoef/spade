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
use spade::rtree::RTree;
use spade::{TwoDimensional, SpadeNum};
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use std::time::{Instant, Duration};
use cgmath::{Point2, BaseNum, EuclideanSpace};
use std::path::Path;
use std::fs::File;
use std::io::{Write};

type DWL = DelaunayWalkLocate;
type Delaunay<Scalar, K, L> = DelaunayTriangulation<Point2<Scalar>, K, L>;

struct BenchSetup<'a, V, K, L> where
    V: TwoDimensional + 'a,
    K: DelaunayKernel<V::Scalar>,
    L: DelaunayLocateStructure<V> {

    title: &'static str,
    delaunay: DelaunayTriangulation<V, K, L>,
    vertices: &'a [V],
    chunk_size: usize,
}

trait Benchmark {
    fn do_bench(&mut self) -> BenchResult;
}

struct BenchResult {
    measurements: Vec<u32>,
    title: &'static str,
}

fn duration_ns(duration: Duration) -> u32 {
    duration.as_secs() as u32 * 1_000_000_000 + duration.subsec_nanos()
}

impl <'a, V, K, L> Benchmark for BenchSetup<'a, V, K, L> where
    V: TwoDimensional + 'a,
    K: DelaunayKernel<V::Scalar>,
    L: DelaunayLocateStructure<V> {
    
    fn do_bench(&mut self) -> BenchResult {
        println!("{}", self.title);
        let mut result = Vec::new();
        let mut sum = 0;
        // let mut delaunay = self.delaunay.clone();
        let vertices = self.vertices;
        for chunk in vertices.chunks(self.chunk_size) {
            let now = Instant::now();
            for vertex in chunk.iter() {
                self.delaunay.insert(vertex.clone());
            }
            let elapsed = now.elapsed();
            let ns = duration_ns(elapsed);
            sum += ns;
            result.push(ns / self.chunk_size as u32);
        };
        assert!(self.delaunay.num_vertices() > self.vertices.len() / 2);
        println!("time / op: {:?}ns", sum / self.vertices.len() as u32);
        BenchResult {
            measurements: result,
            title: self.title,
        }
    }
}

fn main() {

    const SIZE: usize = 400000;
    const CHUNK_SIZE: usize = SIZE / 100;

    let seed = [66719311, 3577332191, 127720921, 1790772261];
    let vertices_f64_walk = random_walk_with_seed_and_origin::<f64>(0.001, SIZE, seed);
    let vertices_i64_walk = random_walk_with_seed_and_origin::<i64>(3, SIZE, seed);
    let vertices_f64_uniform = random_points_with_seed_range_and_origin::<f64>(
            20.0, Point2::new(1e10, -1e10), SIZE, seed);
    let vertices_i64_uniform = random_points_in_range::<i64>(10000, SIZE, seed);

    let mut benches: Vec<Box<Benchmark>> = vec![
        // F64 Benchmarks
        Box::new(BenchSetup {
            title: "f64T - uniform - tree lookup",
            delaunay: Delaunay::<f64, TrivialKernel, RTree<_>>::new(),
            vertices: &vertices_f64_uniform,
            chunk_size: CHUNK_SIZE,
        }),

        Box::new(BenchSetup {
            title: "f64T - uniform - walk lookup",
            delaunay: Delaunay::<f64, TrivialKernel, DWL>::new(),
            vertices: &vertices_f64_uniform,
            chunk_size: CHUNK_SIZE,
        }),

        Box::new(BenchSetup {
            title: "f64T - random walk - tree lookup",
            delaunay: Delaunay::<f64, TrivialKernel, RTree<_>>::new(),
            vertices: &vertices_f64_walk,
            chunk_size: CHUNK_SIZE,
        }),
        Box::new(BenchSetup {
            title: "f64T - random walk - walk lookup",
            delaunay: Delaunay::<f64, TrivialKernel, DWL>::new(),
            vertices: &vertices_f64_walk,
            chunk_size: CHUNK_SIZE,
        }),
        // Float Kernel Benchmarks
        // Box::new(BenchSetup {
        //     title: "f64F - uniform - tree lookup",
        //     delaunay: Delaunay::<f64, FloatKernel, RTree<_>>::new(),
        //     vertices: &vertices_f64_uniform,
        //     chunk_size: CHUNK_SIZE,
        // }),

        // Box::new(BenchSetup {
        //     title: "f64F - uniform - walk lookup",
        //     delaunay: Delaunay::<f64, FloatKernel, DWL>::new(),
        //     vertices: &vertices_f64_uniform,
        //     chunk_size: CHUNK_SIZE,
        // }),

        Box::new(BenchSetup {
            title: "f64F - random walk - tree lookup",
            delaunay: Delaunay::<f64, FloatKernel, RTree<_>>::new(),
            vertices: &vertices_f64_walk,
            chunk_size: CHUNK_SIZE,
        }),
        Box::new(BenchSetup {
            title: "f64F - random walk - walk lookup",
            delaunay: Delaunay::<f64, FloatKernel, DWL>::new(),
            vertices: &vertices_f64_walk,
            chunk_size: CHUNK_SIZE,
        }),
        // i64 Benchmarks
        Box::new(BenchSetup {
            title: "i64T - uniform - tree lookup",
            delaunay: Delaunay::<i64, TrivialKernel, RTree<_>>::new(),
            vertices: &vertices_i64_uniform,
            chunk_size: CHUNK_SIZE,
        }),

        // Box::new(BenchSetup {
        //     title: "i64T - uniform - walk lookup",
        //     delaunay: Delaunay::<i64, TrivialKernel, DWL>::new(),
        //     vertices: &vertices_i64_uniform,
        //     chunk_size: CHUNK_SIZE,
        // }),

        Box::new(BenchSetup {
            title: "i64T - random walk - tree lookup",
            delaunay: Delaunay::<i64, TrivialKernel, RTree<_>>::new(),
            vertices: &vertices_i64_walk,
            chunk_size: CHUNK_SIZE,
        }),
        Box::new(BenchSetup {
            title: "i64T - random walk - walk lookup",
            delaunay: Delaunay::<i64, TrivialKernel, DWL>::new(),
            vertices: &vertices_i64_walk,
            chunk_size: CHUNK_SIZE,
        }),

        // Adaptive Int Benchmarks
        // Box::new(BenchSetup {
        //     title: "i64A - uniform - tree lookup",
        //     delaunay: Delaunay::<i64, AdaptiveIntKernel, RTree<_>>::new(),
        //     vertices: &vertices_i64_uniform,
        //     chunk_size: CHUNK_SIZE,
        // }),

        Box::new(BenchSetup {
            title: "i64A - uniform - walk lookup",
            delaunay: Delaunay::<i64, AdaptiveIntKernel, DWL>::new(),
            vertices: &vertices_i64_uniform,
            chunk_size: CHUNK_SIZE,
        }),

        // Box::new(BenchSetup {
        //     title: "i64A - random walk - tree lookup",
        //     delaunay: Delaunay::<i64, AdaptiveIntKernel, RTree<_>>::new(),
        //     vertices: &vertices_i64_walk,
        //     chunk_size: CHUNK_SIZE,
        // }),

        Box::new(BenchSetup {
            title: "i64A - random walk - walk lookup",
            delaunay: Delaunay::<i64, AdaptiveIntKernel, DWL>::new(),
            vertices: &vertices_i64_walk,
            chunk_size: CHUNK_SIZE,
        }),

    ];
    let mut results = Vec::new();
    for config in &mut benches {
        results.push(config.do_bench());
    }

    let mut result_file = File::create(&Path::new("delaunay_analysis.dat")).unwrap();
    let mut print_measurements = |r: &BenchResult| {
        write!(result_file, "\"{}\"\n", r.title).unwrap();
        for (index, time) in r.measurements.iter().enumerate() {
            let size = index * CHUNK_SIZE;
            write!(result_file, "{} {}\n", size, time).unwrap();
        }
        write!(result_file, "\n\n").unwrap();
    };

    for result in &results {
        print_measurements(result);
    }
    println!("Done!");

}

pub fn random_points_in_range<S: SpadeNum + Rand + SampleRange + BaseNum>(range: S, size: usize, seed: [u32; 4]) -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range.clone(), range.clone());
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_points_with_seed_range_and_origin<S: SpadeNum + BaseNum + Copy + Rand + SampleRange>(
    range: S, origin: Point2<S>, size: usize, seed: [u32; 4])
    -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-range, range);
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng) + origin.x;
        let y = range.ind_sample(&mut rng) + origin.y;
        points.push(Point2::new(x, y));
    }
    points    
}

pub fn random_walk_with_seed_and_origin<S: SpadeNum + Rand + SampleRange + BaseNum>(step: S, size: usize, seed: [u32; 4]) -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let rand_range = Range::new(-step, step);
    let mut points = Vec::new();
    let mut last = Point2::origin();
    for _ in 0 .. size {
        let x = rand_range.ind_sample(&mut rng);
        let y = rand_range.ind_sample(&mut rng);
        last = Point2::new(last.x + x, last.y + y);
        points.push(last);
    }
    points
}
