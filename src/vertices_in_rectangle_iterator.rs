use std::collections::HashSet;

use crate::{
    edges_in_rectangle_iterator::{is_point_inside, EdgesInRectangleIterator},
    handles::{FixedVertexHandle, VertexHandle},
    HasPosition, Point2, Triangulation,
};

/// Iterator over a rectangular area of a Delaunay triangulation. Returned by
/// [Triangulation::get_vertices_in_rectangle].
///
/// The item type is [VertexHandle].
#[derive(PartialEq, Clone, Debug)]
pub struct VerticesInRectangleIterator<'a, V, DE = (), UE = (), F = ()>
where
    V: HasPosition,
{
    iterator: EdgesInRectangleIterator<'a, V, DE, UE, F>,
    already_visited: HashSet<FixedVertexHandle>,
    pending: Option<VertexHandle<'a, V, DE, UE, F>>,
}

impl<'a, V, DE, UE, F> VerticesInRectangleIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    pub(crate) fn new<T>(
        triangulation: &'a T,
        lower: Point2<V::Scalar>,
        upper: Point2<V::Scalar>,
    ) -> Self
    where
        T: Triangulation<Vertex = V, DirectedEdge = DE, UndirectedEdge = UE, Face = F>,
    {
        Self {
            iterator: EdgesInRectangleIterator::new(triangulation, lower, upper),
            already_visited: HashSet::new(),
            pending: None,
        }
    }
}

impl<'a, V, DE, UE, F> Iterator for VerticesInRectangleIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    type Item = VertexHandle<'a, V, DE, UE, F>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(result) = self.pending.take() {
            return Some(result);
        }
        let lower = self.iterator.lower;
        let upper = self.iterator.upper;

        for edge in self.iterator.by_ref() {
            let [v0, v1] = edge.vertices();
            let v0_valid = self.already_visited.insert(v0.fix())
                && is_point_inside(lower, upper, v0.position());
            let v1_valid = self.already_visited.insert(v1.fix())
                && is_point_inside(lower, upper, v1.position());
            match (v0_valid, v1_valid) {
                (true, false) => return Some(v0),
                (false, true) => return Some(v1),
                (true, true) => {
                    self.pending = Some(v0);
                    return Some(v1);
                }
                (false, false) => continue,
            }
        }
        None
    }
}
