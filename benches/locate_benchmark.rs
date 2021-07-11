use std::iter::FromIterator;

use criterion::{measurement::WallTime, BenchmarkGroup, Criterion};
use rand::distributions::uniform::SampleUniform;
use spade::{DelaunayTriangulation, LastUsedVertexHintGenerator, Point2, SpadeNum, Triangulation};

use crate::benchmark_utilities::{
    random_walk_distribution, uniform_distribution, uniform_f32, uniform_f64, SEED, SEED2,
};

pub fn locate_benchmark(c: &mut Criterion) {
    const RANGE: f64 = 1.0e9;
    const NUM_ELEMENTS: usize = 100_000;

    fn single_locate_benchmark<I: IntoIterator<Item = Point2<S>>, S: SpadeNum + SampleUniform>(
        group: &mut BenchmarkGroup<WallTime>,
        name: &'static str,
        range: S,
        elements: I,
    ) where
        S::Sampler: Copy,
    {
        let triangulation =
            DelaunayTriangulation::<_, (), (), (), LastUsedVertexHintGenerator>::from_iter(
                uniform_distribution(*SEED, range).take(NUM_ELEMENTS),
            );

        let mut elements = elements.into_iter();

        group.bench_function(name, |b| {
            b.iter_with_setup(
                || elements.next().unwrap(),
                |point| triangulation.locate(point),
            )
        });
    }

    let mut group = c.benchmark_group("locate benchmark (uniform)");

    single_locate_benchmark(&mut group, "locate (f64)", RANGE, uniform_f64());
    single_locate_benchmark(&mut group, "locate (f32)", RANGE as f32, uniform_f32());

    group.finish();

    let mut group = c.benchmark_group("locate benchmark (random walk)");

    const STEP_SIZE: f64 = RANGE / (NUM_ELEMENTS as f64);

    single_locate_benchmark(
        &mut group,
        "locate (f64)",
        RANGE,
        random_walk_distribution(STEP_SIZE, *SEED2),
    );
    single_locate_benchmark(
        &mut group,
        "locate (f32)",
        RANGE as f32,
        random_walk_distribution(STEP_SIZE as f32, *SEED2),
    );

    group.finish();
}
