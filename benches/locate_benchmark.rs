use std::{iter::FromIterator, time::Duration};

use criterion::{measurement::WallTime, BenchmarkGroup, Criterion};
use rand::distributions::uniform::SampleUniform;
use spade::{
    DelaunayTriangulation, HierarchyHintGeneratorWithBranchFactor, HintGenerator, Point2, SpadeNum,
    Triangulation,
};

use crate::benchmark_utilities::{uniform_distribution, uniform_f64, SEED2};

pub fn locate_benchmark(c: &mut Criterion) {
    const RANGE: f64 = 1.0e9;
    const NUM_ELEMENTS: usize = 100_000;

    fn single_locate_benchmark<
        H: HintGenerator<S>,
        I: IntoIterator<Item = Point2<S>>,
        S: SpadeNum + SampleUniform,
    >(
        group: &mut BenchmarkGroup<WallTime>,
        name: String,
        range: S,
        elements: I,
    ) where
        S::Sampler: Copy,
    {
        let triangulation = DelaunayTriangulation::<_, (), (), (), H>::from_iter(
            uniform_distribution(*SEED2, range).take(NUM_ELEMENTS),
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
    group
        .warm_up_time(Duration::from_secs(2))
        .measurement_time(Duration::from_secs(4));

    fn single_hierarchy<const BRANCH_FACTOR: u32>(group: &mut BenchmarkGroup<WallTime>) {
        single_locate_benchmark::<HierarchyHintGeneratorWithBranchFactor<f64, BRANCH_FACTOR>, _, _>(
            group,
            format!("locate (hierarchy<{:02}>), f64", BRANCH_FACTOR),
            RANGE,
            uniform_f64(),
        );
    }

    single_hierarchy::<2>(&mut group);
    single_hierarchy::<3>(&mut group);
    single_hierarchy::<4>(&mut group);
    single_hierarchy::<5>(&mut group);
    single_hierarchy::<8>(&mut group);
    single_hierarchy::<10>(&mut group);
    single_hierarchy::<14>(&mut group);
    single_hierarchy::<16>(&mut group);
    single_hierarchy::<18>(&mut group);
    single_hierarchy::<20>(&mut group);
    single_hierarchy::<32>(&mut group);
    single_hierarchy::<40>(&mut group);
    single_hierarchy::<50>(&mut group);
    single_hierarchy::<64>(&mut group);
    single_hierarchy::<70>(&mut group);
    single_hierarchy::<80>(&mut group);
    single_hierarchy::<90>(&mut group);

    group.finish();
}
