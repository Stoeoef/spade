use criterion::*;

mod benchmark_utilities;
mod input_size_benchmark;
mod locate_benchmark;

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = locate_benchmark::locate_benchmark, input_size_benchmark::input_size_benchmark
}

criterion_main!(benches);
