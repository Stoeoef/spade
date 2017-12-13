extern crate spade;
extern crate cgmath;
extern crate rand;

use std::path::Path;
use std::fs::File;
use std::io::{Write};
use spade::rtree::RTree;
use spade::SpadeNum;
use cgmath::{Point2, BaseNum};
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::range::SampleRange;
use rand::distributions::{Range, IndependentSample};
use std::time::{Instant, Duration};

static CHUNK_SIZES: &[usize] = &[
    100, 200, 400, 700, 1000, 1500,
    5_000, 10_000, 25_000, 50_000, 100_000, 500_000,
    1_000_000, 1_500_000, 2_000_000, 2_500_000];

pub fn main() {
    let points: Vec<Point2<f32>> = random_points_with_seed(*CHUNK_SIZES.last().unwrap(), [32, 21992, 38332, 218882]);

    let mut result_file = File::create(&Path::new("rtree_compare_bulk_operations.dat")).unwrap();
    println!("Running benchmark...");

    let mut times_sequential = Vec::new();
    let mut times_bulk = Vec::new();
    for chunk_size in CHUNK_SIZES {
        println!("Insertion of {} elements..", chunk_size);
        let mut tree = RTree::new();
        let now = Instant::now();
        for p in points[..*chunk_size].iter().cloned() {
            tree.insert(p);
        }
        let elapsed = now.elapsed();
        times_sequential.push((*chunk_size as u32, duration_ms(elapsed)));
        let now = Instant::now();
        let points_cpy = points[..*chunk_size].to_vec();
        RTree::bulk_load(points_cpy);
        let elapsed = now.elapsed();
        times_bulk.push((*chunk_size as u32, duration_ms(elapsed)));
    }

    // Print all measurements to a file
    let mut print_measurements = |description: &str, measurements: &Vec<(u32, f32)>| {
        write!(result_file, "\"{}\"\n", description).unwrap();
        for &(size, time) in measurements {
            write!(result_file, "{} {}\n", size, time).unwrap();
        }
        write!(result_file, "\n\n").unwrap();
    };

    print_measurements("sequential loading", &times_sequential);
    print_measurements("bulk loading", &times_bulk);
}

fn duration_ms(duration: Duration) -> f32 {
    duration.as_secs() as f32 * 1000. + duration.subsec_nanos() as f32 / 1_000_000.
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
