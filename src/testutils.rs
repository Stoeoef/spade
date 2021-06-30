// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(missing_docs)]
use crate::rtree::RTree;
use crate::traits::SpadeNum;
use cgmath::{BaseFloat, BaseNum, Point2};
use rand::distributions::uniform::SampleUniform;
use rand::{Rng, SeedableRng};
use rand_hc::Hc128Rng;

pub fn random_points_in_range<S: SpadeNum + SampleUniform + BaseNum>(
    range: S,
    size: usize,
    seed: &[u8; 32],
) -> Vec<Point2<S>> {
    let mut rng = Hc128Rng::from_seed(*seed);
    let mut points = Vec::with_capacity(size);
    for _ in 0..size {
        let x = rng.gen_range(-range..range);
        let y = rng.gen_range(-range..range);
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_points_with_seed<S: SpadeNum + BaseFloat + SampleUniform>(
    size: usize,
    seed: &[u8; 32],
) -> Vec<Point2<S>> {
    random_points_in_range(S::one(), size, seed)
}

pub fn create_random_tree<S: SpadeNum + BaseFloat + SampleUniform>(
    size: usize,
    seed: &[u8; 32],
) -> (RTree<Point2<S>>, Vec<Point2<S>>) {
    let mut tree = RTree::new();

    let points = random_points_with_seed(size, seed);
    for point in points.iter() {
        tree.insert(*point)
    }
    (tree, points)
}
