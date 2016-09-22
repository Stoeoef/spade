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

use spade::{DelaunayTriangulation, DelaunayKernel,
            AdaptiveIntKernel, TrivialKernel, FloatKernel, VectorN};
use spade::testutils::*;
use time::Duration;
use cgmath::Vector2;

fn bench<V: VectorN, K: DelaunayKernel<V::Scalar>>(vs: &[V], title: &str) {
    println!("{}", title);
    let mut delaunay: DelaunayTriangulation<V, K> = DelaunayTriangulation::new();
    let time = Duration::span(|| {
        for vertex in vs.iter() {
            delaunay.insert(vertex.clone());
        }
    });
    assert!(delaunay.num_vertices() > vs.len() / 2);
    println!("time / op: {:?}ns", time.num_nanoseconds().unwrap() / vs.len() as i64);
}

fn main() {

    const SIZE: usize = 400000;

    let seed = [661311, 350191, 123021, 231261];
    let vertices_f64 = random_points_with_seed_range_and_origin::<f64>(
        20.0, Vector2::new(1e10, -1e10), SIZE, seed);
    let vertices_i64 = random_points_in_range::<i64>(10000, SIZE, seed);
    let vertices_large_range = random_points_in_range::<i64>(
        1000000000, SIZE, seed);

    bench::<_, TrivialKernel>(&vertices_f64, "f64 benchmark");
    bench::<_, TrivialKernel>(&vertices_i64, "i64 benchmark");
    bench::<_, AdaptiveIntKernel>(&vertices_large_range, "AdaptiveIntKernel benchmark");
    bench::<_, FloatKernel>(&vertices_f64, "FloatKernel benchmark");
}
