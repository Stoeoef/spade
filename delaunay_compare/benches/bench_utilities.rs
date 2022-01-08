use num_traits::Float;
use rand::rngs::StdRng;
use rand::{distributions::uniform::SampleUniform, Rng, SeedableRng};

pub const SEED: &[u8; 32] = b"\xfb\xdc\x4e\xa0\x30\xde\x82\xba\x69\x97\x3c\x52\x49\x4d\x00\xca
\x5c\x21\xa3\x8d\x5c\xf2\x34\x4e\x58\x7d\x80\x16\x66\x23\x30";

pub const RANGE: f64 = 1.0e9;

pub fn uniform_distribution<S>(seed: [u8; 32], range: S) -> impl Iterator<Item = [S; 2]>
where
    S: SampleUniform + Float,
    S::Sampler: Copy,
{
    let range = rand::distributions::Uniform::new_inclusive(-range, range);
    let mut rng = StdRng::from_seed(seed);
    std::iter::from_fn(move || Some([rng.sample(range), rng.sample(range)]))
}

pub fn uniform_f64() -> impl Iterator<Item = [f64; 2]> {
    uniform_distribution(*SEED, RANGE)
}

pub fn walk_f64() -> impl Iterator<Item = [f64; 2]> {
    random_walk_distribution(1.0, *SEED)
}

pub fn random_walk_distribution<S>(step_size: S, seed: [u8; 32]) -> impl Iterator<Item = [S; 2]>
where
    S: SampleUniform + Float,
    S::Sampler: Copy,
{
    let range = rand::distributions::Uniform::new_inclusive(-step_size, step_size);
    let mut last_x = S::zero();
    let mut last_y = S::one();

    let mut rng = StdRng::from_seed(seed);
    let step_fn = move || {
        last_x = last_x + rng.sample(range);
        last_y = last_y + rng.sample(range);

        Some([last_x, last_y])
    };
    std::iter::from_fn(step_fn)
}
