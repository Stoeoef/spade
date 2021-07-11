use rand::{distributions::uniform::SampleUniform, Rng, SeedableRng};
use rand_hc::Hc128Rng;
use spade::{Point2, SpadeNum};

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
    let mut rng = Hc128Rng::from_seed(seed);
    std::iter::from_fn(move || Some(Point2::new(rng.sample(range), rng.sample(range))))
}

pub fn uniform_f64() -> impl Iterator<Item = Point2<f64>> {
    uniform_distribution(*SEED, RANGE)
}

pub fn uniform_f32() -> impl Iterator<Item = Point2<f32>> {
    uniform_distribution(*SEED, RANGE as f32)
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

    let mut rng = Hc128Rng::from_seed(seed);
    let step_fn = move || {
        last_x = last_x + rng.sample(range);
        last_y = last_y + rng.sample(range);

        Some(Point2::new(last_x, last_y))
    };
    std::iter::from_fn(step_fn)
}
