use criterion::measurement::WallTime;
use criterion::*;

use rand::distributions::uniform::SampleUniform;

use spade::{
    DelaunayTriangulation, HierarchyHintGeneratorWithBranchFactor, LastUsedVertexHintGenerator,
    Point2, SpadeNum, Triangulation,
};

use crate::benchmark_utilities::*;

type LastUsedVertexTriangulation<S> =
    DelaunayTriangulation<Point2<S>, (), (), (), LastUsedVertexHintGenerator>;

type HierarchyTriangulation<S, const SHIFT: u32> =
    DelaunayTriangulation<Point2<S>, (), (), (), HierarchyHintGeneratorWithBranchFactor<S, SHIFT>>;

pub fn bulk_load_benchmark(c: &mut Criterion) {
    fn single_bulk_load_benchmark<
        T: Triangulation<Vertex = Point2<S>>,
        S: SpadeNum + SampleUniform,
    >(
        group: &mut BenchmarkGroup<WallTime>,
        name: &'static str,
        range: S,
        is_uniform: bool,
        use_bulk_load: bool,
        sizes: &[usize],
    ) where
        S::Sampler: Copy,
    {
        for size in sizes {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::new(name, size), &size, |b, &size| {
                let dist = || uniform_distribution(*SEED, range).take(*size);
                let uniform_dist = || random_walk_distribution(S::one(), *SEED).take(*size);

                if is_uniform {
                    if use_bulk_load {
                        b.iter(|| T::bulk_load(dist().collect()))
                    } else {
                        b.iter(|| T::from_iter(dist()))
                    }
                } else {
                    if use_bulk_load {
                        b.iter(|| T::bulk_load(uniform_dist().collect()))
                    } else {
                        b.iter(|| T::from_iter(uniform_dist()))
                    }
                }
            });
        }
    }

    let mut group = c.benchmark_group("bulk load benchmark (small)");

    let ref small_sizes = [7000, 20_000, 40_000, 100_000];

    // single_bulk_load_benchmark::<HierarchyTriangulation<_, 16>, _>(
    //     &mut group,
    //     "hierarchy16 (f64), uniform insertion",
    //     RANGE,
    //     true,
    //     false,
    //     small_sizes,
    // );

    single_bulk_load_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f64, uniform insertion, bulk)",
        RANGE as f64,
        true,
        true,
        small_sizes,
    );

    // single_bulk_load_benchmark::<HierarchyTriangulation<_, 16>, _>(
    //     &mut group,
    //     "hierarchy16 (f64, local insertion)",
    //     RANGE,
    //     false,
    //     false,
    //     small_sizes,
    // );

    // single_bulk_load_benchmark::<LastUsedVertexTriangulation<_>, _>(
    //     &mut group,
    //     "last used vertex (f64, local insertion)",
    //     RANGE,
    //     false,
    //     false,
    //     small_sizes,
    // );

    group.finish();

    // let mut group = c.benchmark_group("insert benchmark (local insertion)");

    // const STEP_SIZE: f64 = 1.0;
    // let ref big_sizes = [8000, 32_000, 65_000, 85_000, 120_000, 180_000];

    // single_bulk_load_benchmark::<HierarchyTriangulation<_, 16>, _>(
    //     &mut group,
    //     "hierarchy16 (f64)",
    //     RANGE as f64,
    //     true,
    //     big_sizes,
    // );

    // single_bulk_load_benchmark::<LastUsedVertexTriangulation<_>, _>(
    //     &mut group,
    //     "last used vertex (f64)",
    //     STEP_SIZE,
    //     false,
    //     big_sizes,
    // );

    // group.finish();
}
