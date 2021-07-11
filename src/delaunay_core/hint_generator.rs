use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(feature = "rtree")]
pub use rtree_hint_generator::RTreeHintGenerator;

use crate::{Point2, SpadeNum};

use super::FixedVertexHandle;

#[cfg(feature = "serde")]
use serde_crate::{Deserialize, Serialize};

/// A structure used to speed up common operations on delaunay triangulations like insertion and geometry queries by providing
/// hints on where to start searching for elements.
///
/// Without a hint, these operations run in `O(sqrt(n))` for `n` uniformly distributed vertices. Most time is spent by
/// "walking" to the queried site (e.g. the face that is being inserted to), starting at a random vertex. A hint generator can
/// speed this up by either using heuristics or a spatial data structure to determine where to start walking closer to the target
/// site.
///
/// Hints can also be given manually by using the `...with_hint` methods (e.g.
/// [Triangulation::insert_with_hint](crate::Triangulation::insert_with_hint))
///
/// Usually, you should not need to implement this trait. Spade currently implements two common hint generators that should
/// fulfill most needs:
///  - A heuristic that uses the last inserted vertex as hint ([LastUsedVertexHintGenerator])
///  - A hint generator based on r-trees that improves walk time to `O(log(n))` ([RTreeHintGenerator])
///
pub trait HintGenerator<S: SpadeNum>: Default {
    /// Returns a vertex handle that should be close to a given position.
    ///
    /// The returned vertex handle may be invalid.
    fn get_hint(&self, position: Point2<S>) -> FixedVertexHandle;

    /// Notifies the hint generator that an element was looked up
    fn notify_vertex_lookup(&self, vertex: FixedVertexHandle);
    /// Notifies the hint generator that a new vertex is inserted
    fn notify_vertex_inserted(&mut self, vertex: FixedVertexHandle, vertex_position: Point2<S>);
    /// Notifies the hint generator that a vertex was removed
    fn notify_vertex_removed(&mut self, vertex: FixedVertexHandle, vertex_position: Point2<S>);
}

/// A hint generator that returns the last used vertex as hint.
///
/// This is useful if multiple insertion or locate queries are spatially close instead of randomly distributed.
/// The run time of insertion and locate queries will be bounded by a constant in this case.
///
/// This heuristic takes only a constant additional amount of memory.
/// Since this is a very good out of the box heuristic without any strong disadvantage it is being used by Spade by default.
#[derive(Default, Debug)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
pub struct LastUsedVertexHintGenerator {
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
    fn get_hint(&self, _: Point2<S>) -> FixedVertexHandle {
        FixedVertexHandle::new(self.index.load(Ordering::Relaxed))
    }

    fn notify_vertex_lookup(&self, vertex: FixedVertexHandle) {
        self.index.store(vertex.index(), Ordering::Relaxed);
    }

    fn notify_vertex_inserted(&mut self, vertex: FixedVertexHandle, _: Point2<S>) {
        <Self as HintGenerator<S>>::notify_vertex_lookup(self, vertex);
    }

    fn notify_vertex_removed(&mut self, vertex: FixedVertexHandle, _: Point2<S>) {
        // Use the previous vertex handle as next hint. This should be a good hint if vertices
        // were inserted in local batches.
        let hint = FixedVertexHandle::new(vertex.index().saturating_sub(1));
        <Self as HintGenerator<S>>::notify_vertex_lookup(self, hint);
    }
}

#[cfg(feature = "rtree")]
mod rtree_hint_generator {
    #[cfg(feature = "serde")]
    use serde_crate::{Deserialize, Serialize};

    use rstar::{primitives::PointWithData, Envelope, RTree, RTreeNum, SelectionFunction, AABB};

    use crate::{handles::FixedVertexHandle, HintGenerator, Point2, SpadeNum};

    /// A hint generator based on an underlying r-tree.
    ///
    /// **Requires the `rtree` feature.
    /// An r-tree is a data structure optimized for finding the nearest neighbor in a point set quickly.
    /// This generator comes with a certain cost (both in terms of storage and increased best case run time)
    /// and should only be used for _random access_ queries.
    ///
    /// In this scenario, this heuristic improves the average locate and insertion time from `O(sqrt(n))` to
    /// `O(log(n))` for n inserted vertices. However, the best-case time increases from `O(1)` to `O(log(n))`.
    ///
    /// It requires `O(n)` additional memory.
    ///
    /// # Example
    /// The hint generator is specified as a type parameter of [DelaunayTriangulation](crate::DelaunayTriangulation):
    /// ```
    /// use spade::{DelaunayTriangulation, Triangulation, Point2, RTreeHintGenerator};
    ///
    /// let triangulation = DelaunayTriangulation::<Point2<f32>, (), (), (), RTreeHintGenerator<f32>>::new();
    /// ```
    #[derive(Clone, Debug)]
    #[cfg_attr(
        feature = "serde",
        derive(Serialize, Deserialize),
        serde(crate = "serde_crate")
    )]

    pub struct RTreeHintGenerator<S: RTreeNum> {
        rtree: RTree<PointWithData<FixedVertexHandle, [S; 2]>>,
    }

    impl<S: RTreeNum> Default for RTreeHintGenerator<S> {
        fn default() -> Self {
            Self {
                rtree: RTree::default(),
            }
        }
    }

    impl<S: RTreeNum + SpadeNum> HintGenerator<S> for RTreeHintGenerator<S> {
        fn get_hint(&self, position: Point2<S>) -> FixedVertexHandle {
            let selector = CloseNeighborSelectionFunction {
                target_point: [position.x, position.y],
            };
            self.rtree
                .locate_with_selection_function(selector)
                .map(|entry| entry.data)
                .next()
                .unwrap_or(FixedVertexHandle::new(0))
        }

        fn notify_vertex_lookup(&self, _: FixedVertexHandle) {}

        fn notify_vertex_inserted(
            &mut self,
            vertex: FixedVertexHandle,
            vertex_position: Point2<S>,
        ) {
            self.rtree.insert(PointWithData::new(
                vertex,
                [vertex_position.x, vertex_position.y],
            ));
        }

        fn notify_vertex_removed(&mut self, _: FixedVertexHandle, vertex_position: Point2<S>) {
            self.rtree
                .remove_at_point(&[vertex_position.x, vertex_position.y]);
        }
    }

    struct CloseNeighborSelectionFunction<S> {
        target_point: [S; 2],
    }

    /// Performs a depth first search and stop at the first leaf.
    ///
    /// This will usually not be the closest neighbor but makes up for it by being much faster.
    impl<S: RTreeNum> SelectionFunction<PointWithData<FixedVertexHandle, [S; 2]>>
        for CloseNeighborSelectionFunction<S>
    {
        fn should_unpack_parent(&self, envelope: &AABB<[S; 2]>) -> bool {
            envelope.contains_point(&self.target_point)
        }
    }
}

#[cfg(test)]
mod test {

    #[cfg(feature = "rtree")]
    #[allow(unused)]
    fn make_sure_hint_generators_are_send_and_sync() {
        // This just needs to compile
        fn foo<T: Send + Sync>(_: T) {}

        foo(super::LastUsedVertexHintGenerator::default());
        foo(super::RTreeHintGenerator::<f64>::default());
    }
}
