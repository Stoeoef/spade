use std::fmt::{self, Display, Formatter};

use criterion::{measurement::WallTime, BenchmarkGroup, BenchmarkId, Throughput};
use rand::{distributions::uniform::SampleUniform, Rng, SeedableRng};
use spade::{
    DelaunayTriangulation, HierarchyHintGenerator, LastUsedVertexHintGenerator, Point2, SpadeNum,
    Triangulation,
};

type LastUsedTriangulation =
    DelaunayTriangulation<Point2<f64>, (), (), (), LastUsedVertexHintGenerator>;
type HierarchyTriangulation =
    DelaunayTriangulation<Point2<f64>, (), (), (), HierarchyHintGenerator<f64>>;

pub const SEED: &[u8; 32] = b"\xfb\xdc\x4e\xa0\x30\xde\x82\xba\x69\x97\x3c\x52\x49\x4d\x00\xca
\x5c\x21\xa3\x8d\x5c\xf2\x34\x4e\x58\x7d\x80\x16\x66\x23\x30";

pub const SEED2: &[u8; 32] = b"\x65\xb8\x8e\x3c\x9f\xd8\xb4\x92\xf2\x2e\xcd\x18\xf7\xbc\x2b\x29
dd\xd2\xf3\xa8\x3b\xf4\xe5\x65\xa0\x10\x60\xb1\xce\x22";

pub const RANGE: f64 = 1.0e9;

pub fn uniform_distribution<S: SpadeNum + SampleUniform>(
    seed: [u8; 32],
    range: S,
) -> impl Iterator<Item = Point2<S>>
where
    S::Sampler: Copy,
{
    let range = rand::distributions::Uniform::new_inclusive(-range, range);
    let mut rng = rand::rngs::StdRng::from_seed(seed);
    std::iter::from_fn(move || Some(Point2::new(rng.sample(range), rng.sample(range))))
}

pub fn uniform_f64() -> impl Iterator<Item = Point2<f64>> {
    uniform_distribution(*SEED, RANGE)
}

pub fn random_walk_distribution<S: SpadeNum + SampleUniform>(
    step_size: S,
    seed: [u8; 32],
) -> impl Iterator<Item = Point2<S>>
where
    S::Sampler: Copy,
{
    let range = rand::distributions::Uniform::new_inclusive(-step_size, step_size);
    let mut last_x = S::zero();
    let mut last_y = S::one();

    let mut rng = rand::rngs::StdRng::from_seed(seed);
    let step_fn = move || {
        last_x = last_x + rng.sample(range);
        last_y = last_y + rng.sample(range);

        Some(Point2::new(last_x, last_y))
    };
    std::iter::from_fn(step_fn)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CreationMethod {
    BulkLoad,
    Incremental,
}

impl Display for CreationMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            CreationMethod::BulkLoad => "bulk loading",
            CreationMethod::Incremental => "incremental loading",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SampleSize {
    LargeSampleSet,
    SmallSampleSet,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HintGeneratorType {
    LastUsedVertex,
    Hierarchy,
}

impl Display for HintGeneratorType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            HintGeneratorType::LastUsedVertex => "last used vertex heuristic",
            HintGeneratorType::Hierarchy => "hierarchy lookup",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SampleDistribution {
    RandomWalk,
    Uniform,
}

impl Display for SampleDistribution {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            SampleDistribution::RandomWalk => "locally clustered",
            SampleDistribution::Uniform => "uniformly distributed",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CreationBenchConfig {
    pub creation_method: CreationMethod,
    pub sample_size: SampleSize,
    pub hint_generator_type: HintGeneratorType,
    pub sample_distribution: SampleDistribution,
}

impl CreationBenchConfig {
    pub fn apply(&self, group: &mut BenchmarkGroup<WallTime>) {
        let sizes = self.get_sizes();
        let name = self.get_name();

        for size in sizes {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::new(name.clone(), size), &size, |b, &size| {
                let data = self.get_distribution().take(*size).collect::<Vec<_>>();
                match self.creation_method {
                    CreationMethod::BulkLoad => b.iter(|| self.bulk_load(data.clone())),
                    CreationMethod::Incremental => b.iter(|| self.sequential_load(data.clone())),
                }
            });
        }
    }

    fn get_sizes(&self) -> &[usize] {
        match self.sample_size {
            SampleSize::LargeSampleSet => &[100_000, 200_000, 400_000, 800_000],
            SampleSize::SmallSampleSet => &[1000, 2000, 3000, 6000, 10_000, 15_000],
        }
    }

    fn get_name(&self) -> String {
        format!(
            "{creation_method} on {hint_generator_type}, {sample_distribution}",
            creation_method = self.creation_method,
            hint_generator_type = self.hint_generator_type,
            sample_distribution = self.sample_distribution
        )
    }

    fn get_distribution(&self) -> Box<dyn Iterator<Item = Point2<f64>>> {
        let range = 1000.0;
        match self.sample_distribution {
            SampleDistribution::RandomWalk => Box::new(random_walk_distribution(range, *SEED)),
            SampleDistribution::Uniform => Box::new(uniform_distribution(*SEED, range)),
        }
    }

    fn bulk_load(&self, data: Vec<Point2<f64>>) {
        match self.hint_generator_type {
            HintGeneratorType::LastUsedVertex => {
                LastUsedTriangulation::bulk_load(data).unwrap();
            }
            HintGeneratorType::Hierarchy => {
                HierarchyTriangulation::bulk_load(data).unwrap();
            }
        };
    }

    fn sequential_load(&self, data: Vec<Point2<f64>>) {
        match self.hint_generator_type {
            HintGeneratorType::LastUsedVertex => {
                let mut triangulation = LastUsedTriangulation::new();
                for point in data {
                    triangulation.insert(point).unwrap();
                }
            }
            HintGeneratorType::Hierarchy => {
                let mut triangulation = HierarchyTriangulation::new();
                for point in data {
                    triangulation.insert(point).unwrap();
                }
            }
        }
    }
}
