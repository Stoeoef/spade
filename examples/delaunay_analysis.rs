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

extern crate rand;
extern crate spade;
extern crate cgmath;
extern crate time;
extern crate num;

#[macro_use]
extern crate glium;

mod utils;

use rand::XorShiftRng;
use spade::{DelaunayTriangulation, RTree, HasPosition, SpatialObject, BigVec2};
use num::bigint::{ToBigInt, BigInt};
use utils::*;
use time::Duration;
use cgmath::Vector2;

struct PointWithIndex {
    point: Vector2<f64>,
    index: usize,
}

impl HasPosition for PointWithIndex {
    type Vector = Vector2<f64>;

    fn position(&self) -> Vector2<f64> {
        self.point
    }
}

fn main() {
    const MAX_VERTICES: usize = 400000;
    const MAX_VERTICES_SMALL: usize = MAX_VERTICES / 10;

    let seed = [6663111, 351, 23151, 23126];
    let vertices = random_points_with_seed::<f64>(MAX_VERTICES, seed);
    // let mut delaunay = DelaunayTriangulation::new();
    // println!("f64 benchmark");
    // let time = Duration::span(|| {
    //     for vertex in vertices.iter() {
    //         delaunay.insert(*vertex);
    //     }
    // });
    // println!("time: {:?}", time);
    // println!("i64 benchmark");
    let vertices = random_points_with_seed_and_range::<i64>(10000, MAX_VERTICES, seed);
    // let mut delaunay = DelaunayTriangulation::new();
    // let time = Duration::span(|| {
    //     for vertex in vertices.iter() {
    //         delaunay.insert(*vertex);
    //     }
    // });
    // println!("time: {:?}", time);
    
    println!("BigInt benchmark");
    let vertices: Vec<BigVec2<BigInt>> = vertices.iter().take(MAX_VERTICES_SMALL)
        .map(|v| (*v).into()).collect();
    let mut delaunay = DelaunayTriangulation::new();
    let time = Duration::span(|| {
        for vertex in vertices.iter() {
            delaunay.insert(vertex.clone());
        }
    });
    println!("time: {:?}", time);
}
