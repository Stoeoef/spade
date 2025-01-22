use core::{
    f64,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{
    DelaunayTriangulation, HasPosition, Point2, SpadeNum, Triangulation, TriangulationExt,
};

use super::FixedVertexHandle;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use alloc::vec::Vec;

/// A structure used to speed up common operations on delaunay triangulations like insertion and geometry queries by providing
/// hints on where to start searching for elements.
///
/// Without a hint, these operations run in `O(sqrt(n))` for `n` uniformly distributed vertices. Most time is spent by
/// "walking" to the queried site (e.g. the face that is being inserted to), starting at a random vertex. A hint generator can
/// speed this up by either using heuristics or a spatial data structure to determine where to start walking closer to the target
/// site.
///
/// Hints can also be given manually by using the `...with_hint` methods (e.g.
/// [Triangulation::insert_with_hint])
///
/// Usually, you should not need to implement this trait. Spade currently implements two common hint generators that should
/// fulfill most needs:
///  - A heuristic that uses the last inserted vertex as hint ([LastUsedVertexHintGenerator])
///  - A hint generator based on a hierarchy of triangulations that improves walk time to `O(log(n))`
///     ([HierarchyHintGenerator])
pub trait HintGenerator<S: SpadeNum>: Default {
    /// Returns a vertex handle that should be close to a given position.
    ///
    /// The returned vertex handle may be invalid.
    fn get_hint<TR: Triangulation<Vertex = V>, V: HasPosition<Scalar = S>>(
        &self,
        position: Point2<S>,
        triangulation: &TR,
    ) -> FixedVertexHandle;

    /// Notifies the hint generator that an element was looked up
    fn notify_vertex_lookup(&self, vertex: FixedVertexHandle);
    /// Notifies the hint generator that a new vertex is inserted
    fn notify_vertex_inserted(&mut self, vertex: FixedVertexHandle, vertex_position: Point2<S>);
    /// Notifies the hint generator that a vertex was removed
    fn notify_vertex_removed(
        &mut self,
        swapped_in_point: Option<Point2<S>>,
        vertex: FixedVertexHandle,
        vertex_position: Point2<S>,
    );

    /// Creates a new hint generator initialized to give hints for a specific triangulation
    fn initialize_from_triangulation<TR, V>(triangulation: &TR) -> Self
    where
        TR: Triangulation<Vertex = V>,
        V: HasPosition<Scalar = S>;
}

/// A hint generator that returns the last used vertex as hint.
///
/// This is useful if multiple insertion or locate queries are spatially close instead of randomly distributed.
/// The run time of insertion and locate queries will be bounded by a constant in this case.
///
/// This heuristic requires only a constant additional amount of memory.
///
/// # Example
/// ```
/// use spade::{DelaunayTriangulation, LastUsedVertexHintGenerator, Point2, Triangulation};
///
/// type LastUsedVertexTriangulation =
///     DelaunayTriangulation<Point2<f64>, (), (), (), LastUsedVertexHintGenerator>;
///
/// let mut triangulation = LastUsedVertexTriangulation::new();
/// // Start using the triangulation, e.g. by inserting vertices
/// triangulation.insert(Point2::new(0.0, 0.0));
/// ```
#[derive(Default, Debug)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde")
)]
pub struct LastUsedVertexHintGenerator {
    // Serde does not implement `(De)Serialize` for `AtomicUsize` in no_std environments.
    #[cfg_attr(feature = "serde", serde(skip))]
    index: AtomicUsize,
}

impl Clone for LastUsedVertexHintGenerator {
    fn clone(&self) -> Self {
        Self {
            index: AtomicUsize::new(self.index.load(Ordering::Relaxed)),
        }
    }
}

impl<S: SpadeNum> HintGenerator<S> for LastUsedVertexHintGenerator {
    fn get_hint<TR: Triangulation<Vertex = V>, V: HasPosition<Scalar = S>>(
        &self,
        _: Point2<S>,
        _triangulation: &TR,
    ) -> FixedVertexHandle {
        FixedVertexHandle::new(self.index.load(Ordering::Relaxed))
    }

    fn notify_vertex_lookup(&self, vertex: FixedVertexHandle) {
        self.index.store(vertex.index(), Ordering::Relaxed);
    }

    fn notify_vertex_inserted(&mut self, vertex: FixedVertexHandle, _: Point2<S>) {
        <Self as HintGenerator<S>>::notify_vertex_lookup(self, vertex);
    }

    fn notify_vertex_removed(
        &mut self,
        _swapped_in_point: Option<Point2<S>>,
        vertex: FixedVertexHandle,
        _vertex_position: Point2<S>,
    ) {
        // Use the previous vertex handle as next hint. This should be a good hint if vertices
        // were inserted in local batches.
        let hint = FixedVertexHandle::new(vertex.index().saturating_sub(1));
        <Self as HintGenerator<S>>::notify_vertex_lookup(self, hint);
    }

    fn initialize_from_triangulation<TR, V>(_: &TR) -> Self
    where
        TR: Triangulation,
        V: HasPosition<Scalar = S>,
    {
        Self::default()
    }
}

/// A hint generator based on a hierarchy of triangulations optimized for randomly accessing elements of
/// the triangulation.
///
/// Using this hint generator results in an insertion and lookup performance of O(log(n)) by keeping a
/// few layers of sparsely populated Delaunay Triangulations. These layers can then be quickly traversed
/// before diving deeper into the next, more detailed layer where the search is furthermore refined.
///
/// # Type parameters
///  - `S`: The scalar type used by the triangulation
///
/// # Example
/// ```
/// use spade::{Point2, Triangulation, DelaunayTriangulation, HierarchyHintGenerator};
///
/// pub type HierarchyTriangulation = DelaunayTriangulation<Point2<f64>, (), (), (), HierarchyHintGenerator<f64>>;
///
/// let mut triangulation = HierarchyTriangulation::new();
/// // Start using the triangulation, e.g. by inserting vertices
/// triangulation.insert(Point2::new(0.0, 0.0));
/// ```
pub type HierarchyHintGenerator<S> = HierarchyHintGeneratorWithBranchFactor<S, 16>;

#[derive(Debug, Clone)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde")
)]
#[doc(hidden)]
pub struct HierarchyHintGeneratorWithBranchFactor<S: SpadeNum, const BRANCH_FACTOR: u32> {
    hierarchy: Vec<DelaunayTriangulation<Point2<S>>>,
    num_elements_of_base_triangulation: usize,
}

impl<S: SpadeNum, const BRANCH_FACTOR: u32> Default
    for HierarchyHintGeneratorWithBranchFactor<S, BRANCH_FACTOR>
{
    fn default() -> Self {
        Self {
            hierarchy: Vec::new(),
            num_elements_of_base_triangulation: 0,
        }
    }
}

impl<S: SpadeNum, const BRANCH_FACTOR: u32> HintGenerator<S>
    for HierarchyHintGeneratorWithBranchFactor<S, BRANCH_FACTOR>
{
    fn get_hint<TR: Triangulation<Vertex = V>, V: HasPosition<Scalar = S>>(
        &self,
        position: Point2<S>,
        _triangulation: &TR,
    ) -> FixedVertexHandle {
        let mut nearest = FixedVertexHandle::new(0);
        for layer in self.hierarchy.iter().rev().skip(1) {
            nearest = layer.walk_to_nearest_neighbor(nearest, position).fix();
            let hint_generator: &LastUsedVertexHintGenerator = layer.hint_generator();
            <LastUsedVertexHintGenerator as HintGenerator<S>>::notify_vertex_lookup(
                hint_generator,
                nearest,
            );
            nearest = FixedVertexHandle::new(nearest.index() * BRANCH_FACTOR as usize);
        }
        nearest
    }

    fn notify_vertex_lookup(&self, _: FixedVertexHandle) {}

    fn notify_vertex_inserted(&mut self, vertex: FixedVertexHandle, vertex_position: Point2<S>) {
        self.num_elements_of_base_triangulation += 1;

        // Find first layer to insert into. Insert into all higher layers.
        let mut index = vertex.index() as u32;

        let mut remainder = 0;
        for triangulation in &mut self.hierarchy {
            remainder = index % BRANCH_FACTOR;
            index /= BRANCH_FACTOR;

            if remainder == 0 {
                triangulation.insert(vertex_position).unwrap();
            } else {
                break;
            }
        }

        if remainder == 0 {
            let mut new_layer = DelaunayTriangulation::new();
            let position_of_vertex_0 = self
                .hierarchy
                .first()
                .map(|layer| layer.vertex(FixedVertexHandle::new(0)).position())
                .unwrap_or(vertex_position);
            new_layer.insert(position_of_vertex_0).unwrap();
            self.hierarchy.push(new_layer);
        }
    }

    fn notify_vertex_removed(
        &mut self,
        mut swapped_in_point: Option<Point2<S>>,
        vertex: FixedVertexHandle,
        _vertex_position: Point2<S>,
    ) {
        let index = vertex.index() as u32;

        let mut current_divisor = BRANCH_FACTOR;
        self.num_elements_of_base_triangulation -= 1;
        let mut last_layer_size = self.num_elements_of_base_triangulation;

        for triangulation in &mut self.hierarchy {
            let remainder = index % current_divisor;
            let index_to_remove = index / current_divisor;
            current_divisor *= BRANCH_FACTOR;

            if remainder == 0 {
                // The current handle is part of this layer and must be removed.
                if let Some(swapped_point) = swapped_in_point.as_ref() {
                    if (triangulation.num_vertices() - 1) * (BRANCH_FACTOR as usize)
                        != last_layer_size
                    {
                        // Only insert a new element if the swapped element is not already present
                        // in the layer
                        triangulation.insert(*swapped_point).unwrap();
                    }
                }
                triangulation.remove(FixedVertexHandle::new(index_to_remove as usize));
            }

            let prev_num_vertices = last_layer_size as u32;
            // Divide by BRANCH_FACTOR and round up
            let max_num_vertices = (prev_num_vertices + BRANCH_FACTOR - 1) / BRANCH_FACTOR;
            if triangulation.num_vertices() as u32 > max_num_vertices {
                // The layer contains too many elements. Remove the last.
                let vertex_to_pop = FixedVertexHandle::new(triangulation.num_vertices() - 1);
                swapped_in_point = None;
                triangulation.remove(vertex_to_pop);
            }

            last_layer_size = triangulation.num_vertices();
        }

        if let [.., ref before_last, _] = self.hierarchy.as_slice() {
            if before_last.num_vertices() == 1 {
                // Last layer has become irrelevant
                self.hierarchy.pop();
            }
        }
    }

    fn initialize_from_triangulation<TR, V>(triangulation: &TR) -> Self
    where
        TR: Triangulation<Vertex = V>,
        V: HasPosition<Scalar = S>,
    {
        let mut result = Self::default();
        for vertex in triangulation.vertices() {
            result.notify_vertex_inserted(vertex.fix(), vertex.position());
        }
        result
    }
}

/// Using the Jump-and-Walk algorithm to find points.
#[derive(Default, Debug)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde")
)]
pub struct JumpAndWalkHintGenerator;

impl<S: SpadeNum> HintGenerator<S> for JumpAndWalkHintGenerator {
    fn get_hint<TR: Triangulation<Vertex = V>, V: HasPosition<Scalar = S>>(
        &self,
        position: Point2<S>,
        triangulation: &TR,
    ) -> FixedVertexHandle {
        let mut vertices = triangulation.s().vertices();
        let Some(first) = vertices.next() else {
            panic!("no vertices");
        };
        let len = vertices.len();
        let step = (len / approx_cube_root(len)).max(1);

        // Jump
        let mut best_dist_2 = first.position().distance_2(position);
        let mut best_vert = first.fix();
        for vert in vertices.step_by(step) {
            let dist_2 = vert.position().distance_2(position);
            if dist_2 < best_dist_2 {
                best_dist_2 = dist_2;
                best_vert = vert.fix();
            }
        }

        best_vert
    }

    fn notify_vertex_lookup(&self, _vertex: FixedVertexHandle) {}
    fn notify_vertex_inserted(&mut self, _vertex: FixedVertexHandle, _vertex_position: Point2<S>) {}
    fn notify_vertex_removed(
        &mut self,
        _swapped_in_point: Option<Point2<S>>,
        _vertex: FixedVertexHandle,
        _vertex_position: Point2<S>,
    ) {
    }

    fn initialize_from_triangulation<TR, V>(_triangulation: &TR) -> Self
    where
        TR: Triangulation<Vertex = V>,
        V: HasPosition<Scalar = S>,
    {
        Self
    }
}

fn approx_cube_root(value: usize) -> usize {
    const LEADING_ZEROS_TO_GUESS: [u32; 64] = [
        1905390, 1512309, 1200320, 952695, 756155, 600160, 476348, 378078, 300080, 238174, 189039,
        150040, 119087, 94520, 75020, 59544, 47260, 37510, 29772, 23630, 18755, 14886, 11815, 9378,
        7443, 5908, 4689, 3722, 2954, 2345, 1861, 1477, 1173, 931, 739, 587, 466, 370, 294, 233,
        185, 147, 117, 93, 74, 59, 47, 37, 30, 24, 19, 15, 12, 10, 8, 6, 5, 4, 3, 3, 2, 2, 2, 1,
    ];
    if value < 2 {
        return value;
    }
    let index = value.leading_zeros() as usize;
    let mut guess = LEADING_ZEROS_TO_GUESS[index] as usize;

    // One iteration of Newton's method.
    guess = (value / (guess * guess) + guess * 2) / 3;

    guess
}

#[cfg(test)]
mod test {
    use rand::{prelude::SliceRandom, RngCore, SeedableRng};

    use crate::{
        handles::FixedVertexHandle, test_utilities, DelaunayTriangulation, InsertionError, Point2,
        Triangulation, TriangulationExt,
    };

    use super::approx_cube_root;

    use alloc::vec::Vec;

    const BRANCH_FACTOR: u32 = 3;

    type HierarchyTriangulation = DelaunayTriangulation<
        Point2<f64>,
        (),
        (),
        (),
        super::HierarchyHintGeneratorWithBranchFactor<f64, BRANCH_FACTOR>,
    >;

    #[test]
    fn hierarchy_hint_generator_test() -> Result<(), InsertionError> {
        let vertices = test_utilities::random_points_with_seed(1025, test_utilities::SEED);
        let triangulation = HierarchyTriangulation::bulk_load(vertices)?;

        hierarchy_sanity_check(&triangulation);
        Ok(())
    }

    fn hierarchy_sanity_check(triangulation: &HierarchyTriangulation) {
        for vertex in triangulation.vertices() {
            let position = vertex.position();
            let base_index = vertex.fix().index() as u32;
            let mut power = BRANCH_FACTOR;

            for layer in &triangulation.hint_generator().hierarchy {
                if base_index % power == 0 {
                    let corresponding_vertex =
                        FixedVertexHandle::new((base_index / power) as usize);
                    assert_eq!(layer.vertex(corresponding_vertex).position(), position);
                    power *= BRANCH_FACTOR;
                } else {
                    assert!(layer.locate_vertex(position).is_none());
                }
            }
        }
    }

    #[test]
    fn test_hierarchy_hint_generator_bulk_load_small() -> Result<(), InsertionError> {
        let mut rng = rand::rngs::StdRng::from_seed(*test_utilities::SEED);
        let mut seed_fn = || {
            let mut seed = [0u8; 32];
            rng.fill_bytes(seed.as_mut());
            seed
        };

        for size in 0..5 {
            let vertices = test_utilities::random_points_with_seed(size, &seed_fn());
            let triangulation = HierarchyTriangulation::bulk_load(vertices)?;
            hierarchy_sanity_check(&triangulation);
            triangulation.sanity_check();
        }

        for size in 1..20 {
            let vertices = test_utilities::random_points_with_seed(1 + size * 26, &seed_fn());
            let triangulation = HierarchyTriangulation::bulk_load(vertices)?;
            hierarchy_sanity_check(&triangulation);
        }
        Ok(())
    }

    #[test]
    fn hierarchy_hint_generator_removal_test() -> Result<(), InsertionError> {
        let vertices = test_utilities::random_points_with_seed(300, test_utilities::SEED);
        let mut triangulation = HierarchyTriangulation::bulk_load(vertices)?;

        let mut rng = rand::rngs::StdRng::from_seed(*test_utilities::SEED2);
        while let Some(to_remove) = triangulation
            .fixed_vertices()
            .collect::<Vec<_>>()
            .choose(&mut rng)
        {
            triangulation.remove(*to_remove);
            hierarchy_sanity_check(&triangulation);
        }
        Ok(())
    }

    #[test]
    fn test_approx_cube_root() {
        for i in 0usize.. {
            let Some(i3) = (i * i).checked_mul(i) else {
                break;
            };
            let sut = approx_cube_root(i3);
            // Tolerance of 0.015% of the input.
            let tolerance = i3 / 6859;
            assert!(sut.abs_diff(i) <= tolerance, "{sut} ~= {i}");
        }
    }
}
