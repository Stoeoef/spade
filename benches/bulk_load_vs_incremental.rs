use criterion::*;

use crate::benchmark_utilities::*;

pub fn bulk_load_vs_incremental_benchmark(c: &mut Criterion) {
    {
        use CreationMethod::*;
        let mut group = c.benchmark_group("bulk vs incremental loading");
        for creation_method in [BulkLoad, Incremental, JumpAndWalk] {
            let config = CreationBenchConfig {
                creation_method,
                sample_size: SampleSize::SmallSampleSet,
                hint_generator_type: HintGeneratorType::Hierarchy,
                sample_distribution: SampleDistribution::Uniform,
            };
            config.apply(&mut group);
        }
        group.finish();
    };
}
