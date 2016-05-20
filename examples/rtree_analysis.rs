extern crate rtree;
extern crate cgmath;
extern crate time;
extern crate rand;

use rtree::{RTree, RTreeOptions};
use time::Duration;
use rand::{SeedableRng, XorShiftRng};
use rand::distributions::{Range, IndependentSample};
use cgmath::Vector2;
use std::cmp::{min, max};
use std::path::Path;
use std::fs::File;
use std::io::Write;


fn main() {
    // let max_sizes = [4, 5, 6, 10, 20, 40, 80, 150];
    // let min_sizes = [0.2f64, 0.3, 0.4, 0.5, 0.6];
    let max_sizes = [150];
    let min_sizes = [0.5];

    const WRITE_TO_FILE: bool = true;

    let mut result_file = File::create(&Path::new("rtree_analysis.dat")).unwrap();
    if WRITE_TO_FILE {
        write!(result_file, "# max_size min_size time\n").unwrap();
    }
    for max_size in max_sizes.iter() {
        for min_size_factor in min_sizes.iter() {
            let min_size = max(2, (min_size_factor * *max_size as f64) as usize);
            let reinsertion_count = max((min_size as f32 * 0.5) as usize, 1);
            let time = if max_size <= &min_size || reinsertion_count > min_size {
                0
            } else {
                let options = RTreeOptions::new()
                    .set_max_size(*max_size)
                    .set_min_size(min_size)
                    .set_reinsertion_count(reinsertion_count);
                println!("Running new benchmark...");
                println!("Options: {:?}", options);
                let time = bench(options).num_milliseconds();
                println!("Time: {:?} ms", time);
                time
            };
            if WRITE_TO_FILE {
                write!(result_file, "{} {} {}\n", max_size, min_size_factor, time).unwrap();
            }
        }
        if WRITE_TO_FILE {
            write!(result_file, "\n").unwrap();
        }
    }
}

fn random_points(size: usize, seed: [u32; 4]) -> Vec<Vector2<f64>> {
    const SIZE: f32 = 1000.;
    // let mut rng = XorShiftRng::from_seed([1, 3, 3, 7]);
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-SIZE as f64 / 2., SIZE as f64 / 2.);
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points
}

fn bench(options: RTreeOptions) -> Duration {
    const NUM_VERTICES: usize = 20000;
    const NUM_ITERATIONS: usize = 10;
    let points = random_points(NUM_VERTICES, [3, 1, 4, 1]);
    let mut min_duration = Duration::max_value();
    for _ in 0 .. NUM_ITERATIONS {
        let duration = Duration::span(|| {
            let mut rtree = RTree::new_with_options(options.clone());
            for point in points.iter().cloned() {
                rtree.insert(point);
            }
        });
        min_duration = min(min_duration, duration);
    }
    min_duration
}
