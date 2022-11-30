use std::collections::HashSet;

use crate::{
    edges_in_rectangle_iterator,
    handles::{FixedVertexHandle, VertexHandle},
    HasPosition, Point2, Triangulation,
};

/// Iterator over a rectangular area of a Delaunay triangulation. Returned by
/// [Triangulation::get_vertices_in_rectangle].
///
/// The item type is [VertexHandle].
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct VerticesInRectangleIterator<'a, V, DE = (), UE = (), F = ()>
where
    V: HasPosition,
{
    lower: Point2<V::Scalar>,
    upper: Point2<V::Scalar>,
    todo: Vec<VertexHandle<'a, V, DE, UE, F>>,
    already_visited: HashSet<FixedVertexHandle>,
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
        let todo: Vec<_> = triangulation
            .get_edges_in_rectangle(lower, upper)
            .flat_map(|edge| edge.vertices())
            .find(|point| {
                edges_in_rectangle_iterator::is_point_inside(lower, upper, point.position())
            })
            .into_iter()
            .collect();
        let already_visited = todo.iter().map(|v| v.fix()).collect();
        Self {
            lower,
            upper,
            todo,
            already_visited,
        }
    }
}

impl<'a, V, DE, UE, F> Iterator for VerticesInRectangleIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    type Item = VertexHandle<'a, V, DE, UE, F>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(next) = self.todo.pop() {
            let neighbors = next.out_edges().map(|edge| edge.to());
            for neighbor in neighbors {
                if self.already_visited.insert(neighbor.fix()) {
                    self.todo.push(neighbor);
                }
            }
            if edges_in_rectangle_iterator::is_point_inside(self.lower, self.upper, next.position())
            {
                return Some(next);
            }
        }
        None
    }
}
