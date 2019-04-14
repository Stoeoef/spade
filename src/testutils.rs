// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(missing_docs)]
use cgmath::{BaseFloat, Point2, BaseNum};
use crate::rtree::RTree;
use rand::{XorShiftRng, SeedableRng};
use rand::distributions::{Range, Distribution};
use rand::distributions::range::SampleRange;
use crate::traits::SpadeNum;

pub fn random_points_in_range<S: SpadeNum + SampleRange + BaseNum>(range: S, size: usize, seed: &[u8; 16]) -> Vec<Point2<S>> {
    let mut rng = XorShiftRng::from_seed(*seed);
    let range = Range::new(-range, range);
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.sample(&mut rng);
        let y = range.sample(&mut rng);
        points.push(Point2::new(x, y));
    }
    points
}

pub fn random_points_with_seed<S: SpadeNum + BaseFloat + SampleRange>(size: usize, seed: &[u8; 16]) -> Vec<Point2<S>> {
    random_points_in_range(S::one(), size, seed)
}

pub fn create_random_tree<S: SpadeNum + BaseFloat + SampleRange>(
    size: usize, seed: &[u8; 16]) -> (
    RTree<Point2<S>>, Vec<Point2<S>>) {
    let mut tree = RTree::new();
    
    let points = random_points_with_seed(size, seed);
    for point in points.iter() {
        tree.insert(*point)
    }
    (tree, points)
}
