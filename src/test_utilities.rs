#![allow(missing_docs)]
use crate::Point2;
use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;

pub const SEED: &[u8; 32] = b"wPYxAkIiHcEmSBAxQFoXFrpYToCe1B71";
pub const SEED2: &[u8; 32] = b"14LzG37Y9EHTcmLW8vBDqWwtYsCeVVyF";

use alloc::vec::Vec;

pub fn random_points_in_range(range: f64, size: usize, seed: &[u8; 32]) -> Vec<Point2<f64>> {
    let mut rng = rand::rngs::StdRng::from_seed(*seed);
    let range = Uniform::new(-range, range);
    let mut points = Vec::with_capacity(size);
    for _ in 0..size {
        let x = range.sample(&mut rng);
        let y = range.sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_points_with_seed(size: usize, seed: &[u8; 32]) -> Vec<Point2<f64>> {
    random_points_in_range(1.0, size, seed)
}
