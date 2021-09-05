use criterion::*;

mod benchmark_utilities;
mod bulk_load_benchmark;
mod insert_benchmark;
mod locate_benchmark;

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = locate_benchmark::locate_benchmark, insert_benchmark::insert_benchmark, bulk_load_benchmark::bulk_load_benchmark
}

criterion_main!(benches);
