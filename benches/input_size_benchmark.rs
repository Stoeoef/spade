use criterion::measurement::{Measurement, WallTime};
use criterion::*;

use rand::distributions::uniform::SampleUniform;

use spade::{
    DelaunayTriangulation, LastUsedVertexHintGenerator, Point2, RTreeHintGenerator, SpadeNum,
    Triangulation,
};

use crate::benchmark_utilities::*;

type RTreeTriangulation<S> = DelaunayTriangulation<Point2<S>, (), (), (), RTreeHintGenerator<S>>;
type LastUsedVertexTriangulation<S> =
    DelaunayTriangulation<Point2<S>, (), (), (), LastUsedVertexHintGenerator>;

pub fn input_size_benchmark(c: &mut Criterion) {
    fn single_locate_benchmark<
        T: Triangulation<Vertex = Point2<S>>,
        I: IntoIterator<Item = Point2<S>>,
        S: SpadeNum + SampleUniform,
    >(
        group: &mut BenchmarkGroup<WallTime>,
        name: &'static str,
        range: S,
        elements: I,
    ) where
        S::Sampler: Copy,
    {
        // Setup iteration parameters
        let sizes = [500, 2000, 8000, 32_000, 65_000, 85_000, 120_000, 180_000];

        let mut iterator = elements.into_iter();

        for size in sizes {
            group.throughput(Throughput::Elements(size as u64));
            group.bench_with_input(BenchmarkId::new(name, size), &size, |b, &size| {
                // Create triangulation to test queries on
                let triangulation = T::from_iter(uniform_distribution(*SEED2, range).take(size));

                b.iter_with_setup(
                    || iterator.next().unwrap(),
                    |pt| {
                        triangulation.locate(pt);
                    },
                );
            });
        }
    }

    let mut group = c.benchmark_group("input size benchmark");
    group.sample_size(10);

    single_locate_benchmark::<RTreeTriangulation<_>, _, _>(
        &mut group,
        "rtree (f64, uniform)",
        RANGE,
        uniform_f64(),
    );
    single_locate_benchmark::<LastUsedVertexTriangulation<_>, _, _>(
        &mut group,
        "last_used (f64, uniform)",
        RANGE,
        uniform_f64(),
    );

    single_locate_benchmark::<RTreeTriangulation<_>, _, _>(
        &mut group,
        "rtree (f32, uniform)",
        RANGE as f32,
        uniform_f32(),
    );
    single_locate_benchmark::<LastUsedVertexTriangulation<_>, _, _>(
        &mut group,
        "last used (f32, uniform)",
        RANGE as f32,
        uniform_f32(),
    );

    group.finish();
}
