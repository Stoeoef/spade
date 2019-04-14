// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use rand::Rng;

use cgmath::{EuclideanSpace, Point2, Point3, Vector3};
use spade::delaunay::{DelaunayTriangulation, DelaunayWalkLocate, FloatDelaunayTriangulation};
use spade::HasPosition;

use crate::constants::*;
use noise::{NoiseFn, Seedable};

pub type Delaunay = FloatDelaunayTriangulation<PointWithHeight, DelaunayWalkLocate>;

pub struct PointWithHeight {
    point: Point2<f64>,
    pub height: f64,
    pub gradient: Point2<f64>,
    // We don't need the normal for interpolation purposes. We store it only for
    // visualization.
    pub normal: Vector3<f64>,
}

impl HasPosition for PointWithHeight {
    type Point = Point2<f64>;
    fn position(&self) -> Point2<f64> {
        self.point
    }
}

impl PointWithHeight {
    pub fn position_3d(&self) -> Point3<f64> {
        Point3::new(self.point.x, self.point.y, self.height)
    }

    pub fn new(point: Point2<f64>, height: f64) -> PointWithHeight {
        PointWithHeight {
            point: point,
            height: height,
            gradient: Point2::new(0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, 0.0),
        }
    }
}

// Triangulation creation and normal estimation
pub fn generate_random_triangulation() -> Delaunay {
    let mut rng = ::rand::thread_rng();
    let mut delaunay = DelaunayTriangulation::with_walk_locate();
    let noise = ::noise::OpenSimplex::new().set_seed(rng.gen());
    for _ in 0..NUM_POINTS {
        let x = rng.gen_range(-SAMPLE_REGION, SAMPLE_REGION);
        let y = rng.gen_range(-SAMPLE_REGION, SAMPLE_REGION);
        let height = noise.get([x * FREQUENCY, y * FREQUENCY]) * MAX_HEIGHT;
        // Try out some other height functions, like those:
        // let height = (x * x + y * y) * 0.3;
        // let height = (x * 3.).sin() + (y - 2.).exp();
        delaunay.insert(PointWithHeight::new(Point2::new(x, y), height));
    }

    // Note that, for interpolation, we only need the gradients. For visualization
    // purposes, the normals are also generated and stored within the vertices
    delaunay.estimate_gradients(&(|v| v.height), &(|v, g| v.gradient = g));
    delaunay.estimate_normals(
        &(|v| v.height),
        &(|v: &mut PointWithHeight, n: Point3<_>| v.normal = n.to_vec()),
    );

    delaunay
}
