extern crate spade;
extern crate cgmath;
extern crate rand;

use spade::rtree::RTree;
use spade::SpadeNum;
use cgmath::{Point2, BaseNum};
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::range::SampleRange;
use rand::distributions::{Range, IndependentSample};
use std::time::{Instant, Duration};

static CHUNK_SIZES: &[usize] = &[5_000, 10_000, 20_000, 40_000, 90_000, 200_000,
                                 400_000, 800_000, 1_600_000, 3_200_000];

pub fn main() {
    let points: Vec<Point2<f32>> = random_points_with_seed(*CHUNK_SIZES.last().unwrap(), [32, 21992, 38332, 218882]);

    for chunk_size in CHUNK_SIZES {
        let mut tree = RTree::new();
        let now = Instant::now();
        for p in points[..*chunk_size].iter().cloned() {
            tree.insert(p);
        }
        let elapsed = now.elapsed();
        println!("Successive insertion: {} elements, {} ms", chunk_size, duration_ms(elapsed));

        let now = Instant::now();
        let points_cpy = points[..*chunk_size].to_vec();
        RTree::bulk_load(points_cpy);
        let elapsed = now.elapsed();
        println!("Bulk load: {} elements, {} ms", chunk_size, duration_ms(elapsed));
    }
}

fn duration_ms(duration: Duration) -> u32 {
    duration.as_secs() as u32 * 1000 + duration.subsec_nanos() / 1_000_000
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
