#![feature(test)]

extern crate rtree;
extern crate cgmath;
extern crate rand;
extern crate test;

use rtree::{RTree};
use cgmath::{Vector2};
use rand::{XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use test::Bencher;

fn random_points(size: usize, seed: [u32; 4]) -> Vec<Vector2<f64>> {
    const SIZE: f32 = 1000.;
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-SIZE as f64 / 2., SIZE as f64 / 2.);
    let mut points = Vec::with_capacity(size);
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points
}

// #[bench]
// fn delaunay_begin_insertion(b: &mut Bencher) {
//     let mut delaunay = DelaunayTriangulation::new();
//     let mut points = random_points(100, [1, 2, 3, 4]);
//     b.iter(move || {
//         for p in points.drain(..) {
//             delaunay.insert(p);
//         }
//     });
// }

// #[bench]
// fn delaunay_populated_insertion(b: &mut Bencher) {
//     let mut delaunay = DelaunayTriangulation::new();
//     let mut points = random_points(10000, [1, 2, 3, 4]);
//     for p in points.drain(..) {
//         delaunay.insert(p);
//     }

//     let mut points = random_points(100, [1, 3, 3, 7]);
//     b.iter(move || {
//         for p in points.drain(..) {
//             delaunay.insert(p);
//         }
//     });
// }

#[bench]
fn rtree_begin_insertion(b: &mut Bencher) {
    let mut points = random_points(100, [1, 2, 9, 9]);
    b.iter(move || {
        let mut rtree = RTree::new();
        for p in points.drain(..) {
            rtree.insert(p);
        }
    });
}

#[bench]
fn rtree_populated_insertion(b: &mut Bencher) {
    let mut points = random_points(100000, [1, 3, 3, 8]);
    let mut rtree = RTree::new();
    for p in points.drain(..) {
        rtree.insert(p);
    }
    let mut points = random_points(100, [2, 3, 5, 7]);
    b.iter(move || {
        for p in points.drain(..) {
            rtree.insert(p);
        }
    });
}
