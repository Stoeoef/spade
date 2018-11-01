// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use cgmath::Vector2;

pub const SAMPLE_REGION: f64 = 3.5;
pub const FREQUENCY: f64 = 1.;
pub const NUM_POINTS: usize = 120;
pub const MAX_HEIGHT: f64 = 1.5;

pub const GRID_SUBDIVISIONS: usize = 250;
pub const OFFSET: f64 = -0.01;
pub const GRID_SIZE: f64 = SAMPLE_REGION * 1.05;
pub const SCALE: f64 = 2.0 * GRID_SIZE / (GRID_SUBDIVISIONS as f64);
pub const GRID_OFFSET: Vector2<f64> = Vector2 { x: GRID_SIZE, y: GRID_SIZE };
