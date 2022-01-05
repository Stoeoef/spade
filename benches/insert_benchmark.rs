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

pub fn insert_benchmark(c: &mut Criterion) {
    fn single_insert_benchmark<T: Triangulation<Vertex = Point2<S>>, S: SpadeNum + SampleUniform>(
        group: &mut BenchmarkGroup<WallTime>,
        name: &'static str,
        range: S,
        is_uniform: bool,
        sizes: &[usize],
    ) where
        S::Sampler: Copy,
    {
        for size in sizes {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::new(name, size), &size, |b, &size| {
                // Create triangulation to test queries on
                if is_uniform {
                    b.iter(|| {
                        let mut result = T::default();
                        for vertex in uniform_distribution(*SEED, range).take(*size) {
                            result.insert(vertex).unwrap();
                        }
                    });
                } else {
                    b.iter(|| {
                        let mut result = T::default();
                        for vertex in random_walk_distribution(range, *SEED).take(*size) {
                            result.insert(vertex).unwrap();
                        }
                    });
                }
            });
        }
    }

    let mut group = c.benchmark_group("insert benchmark (small)");

    let small_sizes = &[200, 800, 1000, 2000, 3500, 5000, 7000];

    single_insert_benchmark::<HierarchyTriangulation<_, 16>, _>(
        &mut group,
        "hierarchy16 (f64), uniform insertion",
        RANGE,
        true,
        small_sizes,
    );

    single_insert_benchmark::<HierarchyTriangulation<_, 16>, _>(
        &mut group,
        "hierarchy16 (f32, uniform insertion)",
        RANGE as f32,
        true,
        small_sizes,
    );

    single_insert_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f64, uniform insertion)",
        RANGE as f64,
        true,
        small_sizes,
    );

    single_insert_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f32, uniform insertion)",
        RANGE as f32,
        true,
        small_sizes,
    );

    single_insert_benchmark::<HierarchyTriangulation<_, 16>, _>(
        &mut group,
        "hierarchy16 (f64, local insertion)",
        RANGE,
        false,
        small_sizes,
    );

    single_insert_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f64, local insertion)",
        RANGE,
        false,
        small_sizes,
    );

    group.finish();

    let mut group = c.benchmark_group("insert benchmark (local insertion)");

    const STEP_SIZE: f64 = 1.0;
    let big_sizes = &[8000, 32_000, 65_000, 85_000, 120_000, 180_000];

    single_insert_benchmark::<HierarchyTriangulation<_, 16>, _>(
        &mut group,
        "hierarchy16 (f64)",
        RANGE as f64,
        true,
        big_sizes,
    );

    single_insert_benchmark::<HierarchyTriangulation<_, 16>, _>(
        &mut group,
        "hierarchy16 (f32)",
        RANGE as f32,
        true,
        big_sizes,
    );

    single_insert_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f64)",
        STEP_SIZE,
        false,
        big_sizes,
    );

    single_insert_benchmark::<LastUsedVertexTriangulation<_>, _>(
        &mut group,
        "last used vertex (f32)",
        STEP_SIZE as f32,
        false,
        big_sizes,
    );

    group.finish();
}
