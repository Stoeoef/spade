mod bench_utilities;

use std::time::Duration;

use criterion::{
    criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, BenchmarkId, Criterion,
    SamplingMode, Throughput,
};
use delaunay_compare::{cdt_crate, delaunator_crate, spade_crate, DelaunayCrate};

pub fn creation_benchmark(c: &mut Criterion) {
    fn run_single<Crate: DelaunayCrate>(
        group: &mut BenchmarkGroup<WallTime>,
        sizes: &[usize],
        crate_name: &'static str,
    ) {
        for size in sizes {
            group.throughput(Throughput::Elements(*size as u64));

            let name = format!("{} (local insertion)", crate_name);
            let mut delaunay_crate = Crate::default();
            delaunay_crate.init(bench_utilities::walk_f64().take(*size));
            group.bench_with_input(BenchmarkId::new(name, size), &size, |b, _| {
                b.iter(|| delaunay_crate.run_creation())
            });

            let name = format!("{} (uniform)", crate_name);
            let mut delaunay_crate = Crate::default();
            delaunay_crate.init(bench_utilities::uniform_f64().take(*size));
            group.bench_with_input(BenchmarkId::new(name, size), size, |b, _| {
                b.iter(|| delaunay_crate.run_creation())
            });
        }
    }

    fn run_all(mut group: BenchmarkGroup<WallTime>, sizes: &[usize]) {
        run_single::<spade_crate::SpadeCrate>(&mut group, sizes, "spade 2");
        run_single::<spade_crate::SpadeCrateWithHierarchy>(&mut group, sizes, "spade 2 hierarchy");
        run_single::<cdt_crate::CdtCrate>(&mut group, sizes, "cdt");
        run_single::<delaunator_crate::DelaunatorCrate>(&mut group, sizes, "delaunator");
        group.finish();
    }

    let mut group = c.benchmark_group("comparison: creation benchmark (small)");
    let small = [2000, 4000, 6000, 8000, 10_000, 12_000, 14_000];

    group.warm_up_time(Duration::from_secs(1));
    group.sample_size(50);
    group.measurement_time(Duration::from_secs(3));

    run_all(group, &small);

    let mut group = c.benchmark_group("comparison: creation benchmark (big)");
    let big = [50_000, 100_000, 150_000, 200_000, 250_000];
    group.sample_size(10);
    group.sampling_mode(SamplingMode::Flat);

    run_all(group, &big);
}

criterion_group!(benches, creation_benchmark);
criterion_main!(benches);
