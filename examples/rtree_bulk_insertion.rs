extern crate spade;
extern crate cgmath;
extern crate rand;

use std::path::Path;
use std::fs::File;
use std::io::{Write};
use spade::rtree::{RTree, RTreeOptions};
use spade::SpadeNum;
use cgmath::{Point2, BaseNum};
use rand::{XorShiftRng, SeedableRng};
use rand::distributions::range::SampleRange;
use rand::distributions::{Range, Distribution};
use std::time::{Instant, Duration};

struct TestTimes {
    chunk_size: usize,
    sequential_insertion_time: f32,
    sequential_insertion_query_time: f32,
    bulk_insertion_time: f32,
    bulk_insertion_query_time: f32,
}

static CHUNK_SIZES: &[usize] = &[
    100, 200, 400, 700, 1000, 1500,
    2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000];
    // 5_000, 10_000, 25_000, 50_000, 100_000, 500_000,
    // 1_000_000];

const QUERY_SIZE: usize = 1000;

struct BenchSetup {
    options: RTreeOptions,
    times: Vec<TestTimes>,
    title: &'static str,
}

impl BenchSetup {
    fn new(min_size: usize, max_size: usize, title: &'static str) -> BenchSetup {
        BenchSetup {
            options: RTreeOptions::new()
                .max_size(max_size)
                .min_size(min_size)
                .reinsertion_count(max_size / 3),
            times: Vec::new(),
            title: title,
        }
    }
}

pub fn main() {
    let insertion_points: Vec<Point2<f32>> = random_points_with_seed(*CHUNK_SIZES.last().unwrap(), b"Back to where it".clone());
    let query_points: Vec<Point2<f32>> = random_points_with_seed(QUERY_SIZE, b"back to...began!".clone());

    let mut benchmarks = vec![
        BenchSetup::new(3, 6, "Default options (3-6)"),
        BenchSetup::new(8, 16, "8-16"),
    ];
    for bench in &mut benchmarks {
        run_benchmark(bench, &insertion_points, &query_points);
    }

    let mut result_file = File::create(&Path::new("rtree_compare_bulk_operations.dat")).unwrap();
    for bench in &benchmarks {
        output_benchmark(&mut result_file, bench)
    }
}

fn run_benchmark(benchmark: &mut BenchSetup, insertion_points: &Vec<Point2<f32>>, query_points: &Vec<Point2<f32>>) {
    let &mut BenchSetup {
        ref title,
        ref mut times,
        ref options
    } = benchmark;

    println!("Running benchmark \"{}\"", title);
    for chunk_size in CHUNK_SIZES {
        // Sequential insertion
        let mut tree = RTree::new_with_options(options.clone());
        let now = Instant::now();
        for p in insertion_points[..*chunk_size].iter().cloned() {
            tree.insert(p);
        }
        let elapsed = now.elapsed();
        let sequential_insertion_time = duration_ns(elapsed) / *chunk_size as f32;
        
        let now = Instant::now();
        for query_point in query_points {
            tree.nearest_neighbor(query_point);
        }
        let elapsed = now.elapsed();
        let sequential_insertion_query_time = duration_ns(elapsed) / query_points.len() as f32;

        // Bulk loading
        let now = Instant::now();
        let points_cpy = insertion_points[..*chunk_size].to_vec();
        RTree::bulk_load_with_options(options.clone(), points_cpy);
        let elapsed = now.elapsed();
        let bulk_insertion_time = duration_ns(elapsed) / *chunk_size as f32;

        let now = Instant::now();
        for query_point in query_points {
            tree.nearest_neighbor(query_point);
        }
        let elapsed = now.elapsed();
        let bulk_insertion_query_time = duration_ns(elapsed) / query_points.len() as f32;
        
        let new_times = TestTimes {
            chunk_size: *chunk_size,
            sequential_insertion_time,
            sequential_insertion_query_time,
            bulk_insertion_time,
            bulk_insertion_query_time,
        };
        times.push(new_times);
        print!(".");
        ::std::io::stdout().flush().unwrap();
    }
    println!();
}

fn print_measurements<F>(result_file: &mut ::std::fs::File, description: &str, measurements: &Vec<TestTimes>, f: F) 
    where F: Fn(&TestTimes) -> f32 {
    write!(result_file, "\"{}\"\n", description).unwrap();
    for test_time in measurements {
        write!(result_file, "{} {}\n", test_time.chunk_size, f(test_time)).unwrap();
    }
    write!(result_file, "\n\n").unwrap();
}

fn output_benchmark(result_file: &mut ::std::fs::File, benchmark: &BenchSetup) {
    let &BenchSetup {
        ref times,
        ref title,
        ..
    } = benchmark;
        
    let description = format!("{} - sequential query time", title);
    print_measurements(result_file, &description, times, |time| time.sequential_insertion_query_time);

    let description = format!("{} - sequential insertion time", title);
    print_measurements(result_file, &description, times, |time| time.sequential_insertion_time);

    let description = format!("{} - bulk query time", title);
    print_measurements(result_file, &description, times, |time| time.bulk_insertion_query_time);

    let description = format!("{} - bulk insertion time", title);
    print_measurements(result_file, &description, times, |time| time.bulk_insertion_time);
}

fn duration_ns(duration: Duration) -> f32 {
    duration.as_secs() as f32 * 1_000_000_000. + duration.subsec_nanos() as f32
}

fn random_points_with_seed<S: SpadeNum + BaseNum + Copy + SampleRange>(
    size: usize, seed: [u8; 16])
    -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-S::one(), S::one());
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.sample(&mut rng);
        let y = range.sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points    
}
