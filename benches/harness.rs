//! Benchmark harness and related utilities

use cgmath::*;
use rand::distributions::range::SampleRange;
use rand::distributions::{Distribution, Range};
use rand::*;
use spade::*;

pub fn uniform_points_in_range<S: SpadeNum + SampleRange + BaseNum, R: Rng>(
    range: S,
    size: usize,
    rng: &mut R,
) -> Vec<Point2<S>> {
    let range = Range::new(-range.clone(), range.clone());

    (0..size)
        .map(|_| Point2::new(range.sample(rng), range.sample(rng)))
        .collect()
}

pub fn uniform_i64s<R: Rng>(size: usize, rng: &mut R) -> Vec<Point2<i64>> {
    uniform_points_in_range(1 << 10, size, rng)
}

pub fn uniform_f64s<R: Rng>(size: usize, rng: &mut R) -> Vec<Point2<f64>> {
    uniform_points_in_range(f64::one(), size, rng)
}

#[macro_export]
macro_rules! delaunay {
    ($k: ty, $l: ty) => {
        spade::delaunay::DelaunayTriangulation::<_, $k, $l>::new()
    };
}

#[macro_export]
macro_rules! benchmark_kernels {
    ($s: expr, $mac: ident) => {
        ParameterizedBenchmark::new(
            "uniform_f64/float_kernel/tree_locate",
            $mac!(delaunay!(FloatKernel, RTree<_>), uniform_f64s),
            $s,
        )
        .with_function(
            "uniform_f64/float_kernel/walk_locate",
            $mac!(delaunay!(FloatKernel, DelaunayWalkLocate), uniform_f64s),
        )
        .with_function(
            "uniform_f64/trivial_kernel/tree_locate",
            $mac!(delaunay!(TrivialKernel, RTree<_>), uniform_f64s),
        )
        .with_function(
            "uniform_f64/trivial_kernel/walk_locate",
            $mac!(delaunay!(TrivialKernel, DelaunayWalkLocate), uniform_f64s),
        )
        .with_function(
            "uniform_i64/adaptive_kernel/tree_locate",
            $mac!(delaunay!(AdaptiveIntKernel, RTree<_>), uniform_i64s),
        )
        .with_function(
            "uniform_i64/adaptive_kernel/walk_locate",
            $mac!(
                delaunay!(AdaptiveIntKernel, DelaunayWalkLocate),
                uniform_i64s
            ),
        )
    };
}
