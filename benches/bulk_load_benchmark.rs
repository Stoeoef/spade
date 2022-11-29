use criterion::*;

use crate::benchmark_utilities::*;

pub fn bulk_load_benchmark(c: &mut Criterion) {
    use HintGeneratorType::*;
    use SampleDistribution::*;

    let mut group = c.benchmark_group("bulk load benchmark (small)");

    for (hint_generator_type, sample_distribution) in [
        (LastUsedVertex, RandomWalk),
        (LastUsedVertex, Uniform),
        (Hierarchy, Uniform),
        (Hierarchy, RandomWalk),
    ] {
        let config = CreationBenchConfig {
            creation_method: CreationMethod::BulkLoad,
            sample_size: SampleSize::SmallSampleSet,
            hint_generator_type,
            sample_distribution,
        };
        config.apply(&mut group);
    }
    group.finish();

    let mut group = c.benchmark_group("bulk load benchmark (big)");

    let config = CreationBenchConfig {
        creation_method: CreationMethod::BulkLoad,
        sample_size: SampleSize::LargeSampleSet,
        hint_generator_type: LastUsedVertex,
        sample_distribution: Uniform,
    };
    config.apply(&mut group);

    group.finish();
}
