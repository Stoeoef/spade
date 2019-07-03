use criterion::*;
use spade::{
    delaunay::*,
    kernels::*,
    rtree::*,
};
use rand::*;

const SEED: &'static [u8; 16] = b"rand::thread_rng";

fn insertion_benchmarks(c: &mut Criterion) {
    // Setup iteration parameters
    const MAX_SIZE: usize = 15;
    let sizes = (0..MAX_SIZE).step_by(2).map(|s| 1 << s);

    macro_rules! insert_block {
        ($e: expr, $data: expr) => {
            move |b, param| {
                let mut rng = XorShiftRng::from_seed(*SEED);
                let data = ($data)(*param, &mut rng);
                b.iter_with_large_drop(|| {
                    let mut t = $e;
                    for d in &data { t.insert(*d); }
                    t
                })
            }
        }
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
                let mut rng = XorShiftRng::from_seed(*SEED);

                // Create datastructure to test queries on
                let t = {
                    let mut t = $e;
                    let data = ($data)(*param, &mut rng);
                    for d in data { t.insert(d); }
                    t
                };

                use cgmath::Point2;
                b.iter_with_setup(
                    || { Point2::new(rng.gen(), rng.gen()) },
                    |pt| { t.locate(&pt); }
                );

            }
        }
    }
    c.bench("locate", benchmark_kernels!(sizes, locate_block));

}

criterion_group!{
    name = benches;
    config = Criterion::default().sample_size(25);
    targets = insertion_benchmarks, locate_benchmarks
}
criterion_main!(benches);

mod harness;
use harness::*;
