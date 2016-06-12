use cgmath::{BaseFloat, Vector2, one};
use rtree::RTree;
use rand::{Rand, XorShiftRng, SeedableRng};
use rand::distributions::{Range, IndependentSample};
use rand::distributions::range::SampleRange;
use num::Bounded;
use traits::RTreeNum;

pub fn random_points_with_seed<S: BaseFloat + Rand + SampleRange>(size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
    let mut rng = XorShiftRng::from_seed(seed);
    let range = Range::new(-one::<S>(), one());
    let mut points = Vec::new();
    for _ in 0 .. size {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        points.push(Vector2::new(x, y));
    }
    points
}

pub fn create_random_tree<S: BaseFloat + Rand + SampleRange + RTreeNum>(
    size: usize, seed: [u32; 4]) -> (
    RTree<Vector2<S>>, Vec<Vector2<S>>) {
    let mut tree = RTree::new();
    
    let points = random_points_with_seed(size, seed);
    for point in points.iter() {
        tree.insert(point.clone())
    }
    (tree, points)
}
