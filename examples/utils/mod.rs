// Copyright 2016 The RTree Developers. For a full listing of the authors,
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

#![allow(dead_code)]
pub mod exampleapplication;

// Rust generates warnings for this crate when running cargo test, allowing some dead
// code seems to be the less obnoxious solution.
use cgmath::{Vector2, Vector3, Array};
use cgmath::conv::*;
use spade::{BoundingRect};

#[derive(Clone, Copy)]
pub struct Vertex {
    pub pos: [f32; 2],
    pub color: [f32; 3],
}

#[cfg(not(test))]
implement_vertex!(Vertex, pos, color);
impl Vertex {
    pub fn new(pos: [f32; 2], color: [f32; 3]) -> Vertex {
        Vertex { pos: pos, color: color }
    }
}


// wrap around if a value overflows
// fn wrap<S: BaseFloat>(val: S) -> S {
//     let two = one::<S>() + one::<S>();
//     if val > one::<S>() {
//         two - val
//     } else if val < -one::<S>() {
//         two + val
//     } else {
//         val
//     }
// }

// pub fn random_walk<S: BaseFloat + Rand + SampleRange>(
//     size: usize, step_size: S, seed: [u32; 4]) -> Vec<Vector2<S>> {
//     let mut rng = XorShiftRng::from_seed(seed);
//     let range = Range::new(-step_size, step_size);
//     let mut points = Vec::new();
//     let mut last = Vector2::from_value(zero::<S>());
//     for _ in 0 .. size {
//         let dx = range.ind_sample(&mut rng);
//         let dy = range.ind_sample(&mut rng);
//         let mut next: Vector2<S> = Vector2::new(dx, dy) + last;
//         next.x = wrap(next.x);
//         next.y = wrap(next.y);
//         last = next.clone();
//         points.push(next);
//     }
//     points
// }


pub fn push_rectangle(vec: &mut Vec<Vertex>, rect: &BoundingRect<Vector2<f32>>, 
                                    color: &Vector3<f32>) {
    let v0: Vector2<_> = rect.lower();
    let v2: Vector2<_> = rect.upper();
    let v1 = Vector2::new(v2.x, v0.y);
    let v3 = Vector2::new(v0.x, v2.y);
    vec.extend([v0, v1, v1, v2, v2, v3, v3, v0].iter().cloned().map(
        |v| Vertex::new(array2(v), array3(color.clone()))));
}

pub fn push_cross(vec: &mut Vec<Vertex>, pos: &Vector2<f32>,
                                color: &Vector3<f32>) {
    let mut delta =  Vector2::from_value(0.015f32);
    let v0 = *pos + delta;
    let v1 = *pos - delta;
    delta.x *= -1.0;
    let v2 = *pos + delta;
    let v3 = *pos - delta;
    vec.extend([v0, v1, v2, v3].iter().cloned().map(
        |v| Vertex::new(array2(v), array3(color.clone()))));
}

pub fn get_color_for_depth(depth: usize) -> Vector3<f32> {
    match depth {
        0 => Vector3::new(1., 1., 1.),
        1 => Vector3::new(1., 0., 1.) * 0.85,
        2 => Vector3::new(0., 1., 1.) * 0.7,
        3 => Vector3::new(0., 0., 1.) * 0.55,
        4 => Vector3::new(1., 0., 0.) * 0.4,
        _ => Vector3::new(0., 1., 1.) * 0.25,
    }
}
