use criterion::*;
use rand::*;
use rand_hc::Hc128Rng;

use spade::{delaunay::*, kernels::*, rtree::*};

const SEED: &[u8; 32] = b"\xfb\xdc\x4e\xa0\x30\xde\x82\xba\x69\x97\x3c\x52\x49\x4d\x00\xca
\x5c\x21\xa3\x8d\x5c\xf2\x34\x4e\x58\x7d\x80\x16\x66\x23\x30";

fn insertion_benchmarks(c: &mut Criterion) {
    // Setup iteration parameters
    const MAX_SIZE: usize = 15;
    let sizes = (0..MAX_SIZE).step_by(2).map(|s| 1 << s);

    macro_rules! insert_block {
        ($e: expr, $data: expr) => {
            move |b, param| {
                let mut rng = Hc128Rng::from_seed(*SEED);
                let data = ($data)(*param, &mut rng);
                b.iter_with_large_drop(|| {
                    let mut t = $e;
                    for d in &data {
                        t.insert(*d);
                    }
                    t
                })
            }
        };
    }
    c.bench("insert", benchmark_kernels!(sizes, insert_block));
}

fn locate_benchmarks(c: &mut Criterion) {
    // Setup iteration parameters
    const MAX_SIZE: usize = 15;
    let sizes = (0..MAX_SIZE).step_by(2).map(|s| 1 << s);

    macro_rules! locate_block {
        ($e: expr, $data: expr) => {
            move |b, param| {
                let mut rng = Hc128Rng::from_seed(*SEED);

                // Create datastructure to test queries on
                let t = {
                    let mut t = $e;
                    let data = ($data)(*param, &mut rng);
                    for d in data {
                        t.insert(d);
                    }
                    t
                };

                use cgmath::Point2;
                b.iter_with_setup(
                    || Point2::new(rng.gen(), rng.gen()),
                    |pt| {
                        t.locate(&pt);
                    },
                );
            }
        };
    }
    c.bench("locate", benchmark_kernels!(sizes, locate_block));
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(25);
    targets = insertion_benchmarks, locate_benchmarks
}
criterion_main!(benches);

mod harness;
use harness::*;
