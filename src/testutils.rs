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

#![allow(missing_docs)]

use cgmath::{BaseFloat, Vector2, BaseNum};
use rtree::RTree;
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use traits::SpadeNum;

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

pub fn random_points_with_seed<S: SpadeNum + BaseFloat + Rand + SampleRange>(size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    random_points_in_range(S::one(), size, seed)
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


pub fn create_random_tree<S: SpadeNum + BaseFloat + Rand + SampleRange>(
    size: usize, seed: [u32; 4]) -> (
    RTree<Vector2<S>>, Vec<Vector2<S>>) {
    let mut tree = RTree::new();
    
    let points = random_points_with_seed(size, seed);
    for point in points.iter() {
        tree.insert(*point)
    }
    (tree, points)
}
