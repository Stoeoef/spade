use std::time::Duration;

use criterion::*;

mod benchmark_utilities;
mod bulk_load_benchmark;
mod bulk_load_vs_incremental;
mod interpolation_benchmark;
mod locate_benchmark;

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(50).warm_up_time(Duration::from_secs(1));
    targets =
    bulk_load_benchmark::bulk_load_benchmark,
    bulk_load_vs_incremental::bulk_load_vs_incremental_benchmark,
    interpolation_benchmark::interpolation_benchmark,
    locate_benchmark::locate_benchmark,
}

criterion_main!(benches);
