use criterion::*;
use spade::Triangulation;
use spade::{DelaunayTriangulation, FloatTriangulation, Point2};

use crate::benchmark_utilities::*;

pub fn interpolation_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("interpolation benchmarks");

    let points = uniform_distribution(*SEED, 1.0)
        .take(50)
        .collect::<Vec<_>>();
    let triangulation = DelaunayTriangulation::<_>::bulk_load(points).unwrap();

    let query_point = Point2::new(0.5, -0.2);

    group.bench_function("nearest neighbor interpolation", |b| {
        b.iter(|| {
            triangulation
                .nearest_neighbor(query_point)
                .unwrap()
                .position()
        });
    });

    let barycentric = triangulation.barycentric();
    group.bench_function("barycentric interpolation", |b| {
        b.iter(|| barycentric.interpolate(|v| v.position().x + v.position().y, query_point));
    });

    let natural_neighbor = &triangulation.natural_neighbor();
    group.bench_function("natural neighbor interpolation (c0)", |b| {
        b.iter(|| natural_neighbor.interpolate(|v| v.position().x + v.position().y, query_point));
    });

    group.bench_function("natural neighbor interpolation (c1)", |b| {
        b.iter(|| {
            natural_neighbor.interpolate_gradient(
                |v| v.position().x + v.position().y,
                |_| [0.0, 0.0],
                0.5,
                query_point,
            )
        });
    });
    group.finish();
}
