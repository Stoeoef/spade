use super::{CircularIterator, NextBackFn};
use crate::{
    delaunay_core::DCEL,
    handles::{DirectedEdgeHandle, FixedDirectedEdgeHandle},
};

pub struct HullIterator<'a, V, DE, UE, F> {
    inner_iterator: CircularIterator<'a, V, DE, UE, F, HullNextBackFn>,
}

struct HullNextBackFn;

impl NextBackFn for HullNextBackFn {
    fn next<'a, V, DE, UE, F>(
        edge_handle: DirectedEdgeHandle<'a, V, DE, UE, F>,
    ) -> DirectedEdgeHandle<'a, V, DE, UE, F> {
        edge_handle.next()
    }

    fn next_back<'a, V, DE, UE, F>(
        edge_handle: DirectedEdgeHandle<'a, V, DE, UE, F>,
    ) -> DirectedEdgeHandle<'a, V, DE, UE, F> {
        edge_handle.prev()
    }
}

impl<'a, V, DE, UE, F> HullIterator<'a, V, DE, UE, F> {
    pub(crate) fn new(dcel: &'a DCEL<V, DE, UE, F>) -> Self {
        let outer_face = dcel.outer_face();

        let inner_iterator = if let Some(first_edge) = outer_face.adjacent_edge() {
            CircularIterator::new(first_edge)
        } else {
            CircularIterator::new_empty(DirectedEdgeHandle::new(
                dcel,
                FixedDirectedEdgeHandle::new(0),
            ))
        };

        Self { inner_iterator }
    }
}

impl<'a, V, DE, UE, F> Iterator for HullIterator<'a, V, DE, UE, F> {
    type Item = DirectedEdgeHandle<'a, V, DE, UE, F>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner_iterator.next()
    }
}

impl<'a, V, DE, UE, F> DoubleEndedIterator for HullIterator<'a, V, DE, UE, F> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner_iterator.next_back()
    }
}

#[cfg(test)]
mod test {
    use crate::test_utilities::{random_points_with_seed, SEED};
    use crate::{DelaunayTriangulation, Point2, Triangulation};

    #[test]
    fn test_empty_hull() {
        let mut triangulation = DelaunayTriangulation::<Point2<f64>>::new();
        assert!(triangulation.convex_hull().next().is_none());

        triangulation.insert(Point2::new(0.0, 0.0));
        assert!(triangulation.convex_hull().next().is_none());

        triangulation.insert(Point2::new(0.0, 1.0));
        assert_eq!(triangulation.convex_hull().count(), 2);
    }

    #[test]
    fn test_bigger_triangulation() {
        let vertices = random_points_with_seed(100, &SEED);
        let triangulation: DelaunayTriangulation<_> = vertices.into_iter().collect();

        let convex_hull: Vec<_> = triangulation.convex_hull().collect();
        let mut reversed: Vec<_> = triangulation.convex_hull().rev().collect();

        for ch_edge in &convex_hull {
            assert!(ch_edge.face().is_outer())
        }

        reversed.reverse();
        assert_eq!(convex_hull, reversed);

        assert_eq!(triangulation.convex_hull_size(), convex_hull.len());
    }
}
