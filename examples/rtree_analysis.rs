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
use std::io::{Write, stdout};


fn main() {
    // run_insertion_over_time_bench();
    // run_insertion_with_different_parameters_bench();
    run_lookup_with_different_parameters_bench();
}

fn run_insertion_over_time_bench() {
    let max_sizes = [5, 6, 7, 8, 9, 10, 16];
    const NUM_VERTICES: usize = 1000000;
    const CHUNK_SIZE: usize = NUM_VERTICES / 300;
    const NUM_ITERATIONS: usize = 2;

    let mut result_file = File::create(&Path::new("rtree_insertion_over_time.dat")).unwrap();
    write!(result_file, "# max_size num_vertices time\n").unwrap();

    let points = random_points(NUM_VERTICES, [2, 2, 31123, 998]);
    for max_size in max_sizes.iter() {
        let min_size = (*max_size as f32 * 0.45) as usize;
        let reinsertion_count = (min_size as f32 * 0.8) as usize;
        println!("Running new benchmark...");
        println!("max_size: {}", max_size);
        let mut times = [i64::max_value(); NUM_VERTICES / CHUNK_SIZE];
        let mut total = 0;
        for _ in 0 .. NUM_ITERATIONS {
            let mut tree = RTreeOptions::new()
                .set_max_size(*max_size)
                .set_min_size(min_size)
                .set_reinsertion_count(reinsertion_count)
                .build();
            for (i, chunk) in points.chunks(CHUNK_SIZE).enumerate() {
                if chunk.len() < CHUNK_SIZE {
                    continue;
                }
                let time = Duration::span(|| {
                    for point in chunk.iter().cloned() {
                        tree.insert(point);
                    }
                }).num_nanoseconds().unwrap();
                times[i] = min(times[i], time);
                total += time;
            }
            print!(".");
            stdout().flush().unwrap();
        }
        write!(result_file, "\"max_size = {}\"\n", max_size).unwrap();
        for (i, time) in times.iter().enumerate() {
            if *time == i64::max_value() {
                continue;
            }
            write!(result_file, "{} {}\n", i * CHUNK_SIZE, time).unwrap();
        }
        write!(result_file, "\n\n").unwrap();
        println!("\n Total time: {}", total / NUM_ITERATIONS as i64);
    }
}

fn run_lookup_with_different_parameters_bench() {
    // let max_sizes = [6, 8, 10, 20, 40, 80, 150];
    // let reinsertion_counts = [0.1, 0.3, 0.7];
    let max_sizes = [6];
    let reinsertion_counts = [0.3];

    const NUM_START_VERTICES: usize = 200000;
    const NUM_LOOKUP_VERTICES: usize = 400000;

    let points = random_points(NUM_START_VERTICES, [3, 1, 4, 1]);
    let lookup_points = random_points(NUM_LOOKUP_VERTICES, [312, 330, 9, 931]);
    let mut result_file = File::create(
        &Path::new("rtree_lookup_with_different_parameters.dat")).unwrap();

    for max_size in max_sizes.iter().cloned() {
        let min_size = (max_size as f32 * 0.45) as usize;
        for reinsertion_count_factor in reinsertion_counts.iter().cloned() {
            let reinsertion_count = max((max_size as f32 * reinsertion_count_factor) as usize, 1);
            let options = RTreeOptions::new()
                .set_max_size(max_size)
                .set_min_size(min_size)
                .set_reinsertion_count(reinsertion_count);
            println!("Running new benchmark...");
            println!("Options: {:?}", options);
            let mut start_tree = options.build();
            for point in points.iter().cloned() {
                start_tree.insert(point);
            }
            let time = bench(&start_tree, &lookup_points, |t, p| {
                if t.lookup(&p).is_some() { println!("don't optimize me away"); } })
                .num_milliseconds();
            println!("Time: {:?} ms", time);
            write!(result_file, "{} {} {}\n", max_size, reinsertion_count_factor, time).unwrap();
        };
        write!(result_file, "\n").unwrap();
    }
}

fn run_insertion_with_different_parameters_bench() {
    let max_sizes = [4, 5, 6, 7, 10, 20, 40, 80, 150];
    let min_sizes = [0.2f64, 0.4, 0.45, 0.5, 0.6, 0.8];

    const NUM_START_VERTICES: usize = 200000;
    const NUM_INSERTION_VERTICES: usize = 200000;
    let start_vertices = random_points(NUM_START_VERTICES, [3, 1, 4, 1]);
    let insertion_vertices = random_points(NUM_INSERTION_VERTICES, [4556, 99821, 2156126, 22]);

    let mut result_file = File::create(
        &Path::new("rtree_insertion_with_different_parameters.dat")).unwrap();
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
                let mut start_tree = options.build();
                for point in start_vertices.iter().cloned() {
                    start_tree.insert(point);
                }

                let time = bench(&start_tree, &insertion_vertices, 
                                 |t, p| t.insert(p)).num_milliseconds();
                println!("Time: {:?} ms", time);
                time
            };
            write!(result_file, "{} {} {}\n", max_size, min_size_factor, time).unwrap();
        }
        write!(result_file, "\n").unwrap();
    }
}

fn random_points(size: usize, seed: [u32; 4]) -> Vec<Vector2<f64>> {
    const SIZE: f32 = 1000.;

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

fn bench<F>(start_tree: &RTree<Vector2<f64>>, points: &Vec<Vector2<f64>>, f: F) -> Duration 
    where F: Fn(&mut RTree<Vector2<f64>>, Vector2<f64>)
{
    const NUM_ITERATIONS: usize = 5;

    let mut min_duration = Duration::max_value();
    for _ in 0 .. NUM_ITERATIONS {
        let mut clone = start_tree.clone();
        let duration = Duration::span(|| {
            for point in points.iter().cloned() {
                f(&mut clone, point);
            }
        });
        min_duration = min(min_duration, duration);
    }
    min_duration
}
