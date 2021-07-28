use criterion::*;

mod benchmark_utilities;
mod insert_benchmark;
mod locate_benchmark;

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = locate_benchmark::locate_benchmark, insert_benchmark::insert_benchmark
}

criterion_main!(benches);
