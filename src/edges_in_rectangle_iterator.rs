/// An iterator over all edges within a rectangle.
///
/// The tricky part is how to check that an item has already been visited. A naive approach can be
/// implemented with a single `HashSet`. However, this will require O(n) memory
/// (where `n` is the number of returned edges).
///
/// To avoid this, the iterator will first find an initial edge loop (e.g. a face) within the
/// rectangle. Then, this edge loop is expanded with each `next()` call until all edges in the
/// loop lie outside of the rectangle (or no other edge is found).
///
/// This approach makes sure that edges which are completely _within_ the loop can be removed
/// early on as they will not be required for future "already visited" checks. This reduces the
/// maximum required additional memory to some multiple of the longest loop length.
///
/// A directed edge e is expanded by replacing `e` with `e.rev().next()` and `e.rev().prev()` unless
/// they are already a part of the loop. Some additional logic is required to cover outer edges.
///
/// The loops is stored within `remaining_hull_set` and `remaining_hull_vec`. `remaining_hull_set`
/// is the "true" current loop, `remaining_hull_vec` is only used to pop elements from that set.
use std::collections::HashSet;

use crate::{
    delaunay_core::math,
    handles::{DirectedEdgeHandle, UndirectedEdgeHandle},
    HasPosition, Point2, SpadeNum, Triangulation,
};

/// Iterator over a rectangular area of a Delaunay triangulation. Returned by
/// [Triangulation::get_edges_in_rectangle].
///
/// The item type is [UndirectedEdgeHandle].
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct EdgesInRectangleIterator<'a, V, DE = (), UE = (), F = ()>
where
    V: HasPosition,
{
    remaining_hull_set: HashSet<DirectedEdgeHandle<'a, V, DE, UE, F>>,
    remaining_hull_vec: Vec<DirectedEdgeHandle<'a, V, DE, UE, F>>,
    pub(crate) lower: Point2<V::Scalar>,
    pub(crate) upper: Point2<V::Scalar>,
}

impl<'a, V, DE, UE, F> EdgesInRectangleIterator<'a, V, DE, UE, F>
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
        if triangulation.all_vertices_on_line() {
            // Simply pre-filter all edges for degenerate triangulations
            return Self {
                remaining_hull_set: HashSet::new(),
                remaining_hull_vec: triangulation
                    .undirected_edges()
                    .map(|edge| edge.as_directed())
                    .filter(|edge| is_handle_inside(lower, upper, *edge))
                    .collect(),
                lower,
                upper,
            };
        }

        let mid = lower.add(upper).mul(0.5.into());
        let remaining_hull_vec = match triangulation.locate(mid) {
            crate::PositionInTriangulation::OnVertex(vertex) => {
                let out_edge = triangulation.vertex(vertex).out_edge();
                if let Some(out_edge) = out_edge {
                    if let Some(face) = out_edge.face().as_inner() {
                        face.adjacent_edges().into()
                    } else {
                        vec![out_edge, out_edge.rev()]
                    }
                } else {
                    vec![]
                }
            }
            crate::PositionInTriangulation::OnEdge(edge) => {
                let edge = triangulation.directed_edge(edge);
                vec![edge, edge.rev()]
            }
            crate::PositionInTriangulation::OnFace(face) => {
                triangulation.face(face).adjacent_edges().into()
            }
            crate::PositionInTriangulation::OutsideOfConvexHull(edge) => {
                Self::get_intersecting_outside_edge(lower, upper, triangulation.directed_edge(edge))
            }
            crate::PositionInTriangulation::NoTriangulation => vec![],
        };

        Self {
            remaining_hull_set: remaining_hull_vec.iter().copied().collect(),
            remaining_hull_vec,
            lower,
            upper,
        }
    }

    fn get_intersecting_outside_edge(
        lower: Point2<V::Scalar>,
        upper: Point2<V::Scalar>,
        mut edge: DirectedEdgeHandle<'a, V, DE, UE, F>,
    ) -> Vec<DirectedEdgeHandle<'a, V, DE, UE, F>> {
        let v0 = lower.to_f64();
        let v1 = Point2::new(lower.x, upper.y).to_f64();
        let v2 = upper.to_f64();
        let v3 = Point2::new(upper.x, lower.y).to_f64();

        let edges = [(v0, v1), (v1, v2), (v2, v3), (v3, v0)];

        let mut intersections = Vec::with_capacity(4);

        let mut moved_to_next = false;
        let mut moved_to_prev = false;

        loop {
            let e0 = edge.from().position().to_f64();
            let e1 = edge.to().position().to_f64();
            for &(v0, v1) in &edges {
                let (s0, s1) = get_edge_intersections(v0, v1, e0, e1);

                if (0.0..=1.0).contains(&s0) {
                    intersections.push(s1);
                }
            }

            match intersections.as_slice() {
                [] => return vec![],
                &[single] => {
                    if single < 0.0 {
                        edge = edge.prev();
                    } else if single > 1.0 {
                        edge = edge.next();
                    } else {
                        return vec![edge];
                    }
                }
                &[first, second, ..] => {
                    if first < 0.0 && second < 0.0 {
                        edge = edge.prev();
                        moved_to_prev = true;
                    } else if first > 1.0 && second > 1.0 {
                        edge = edge.next();
                        moved_to_next = true;
                    } else {
                        return vec![edge];
                    }
                }
            }
            if moved_to_next && moved_to_prev {
                // Special case: Sometimes, all edges will cause two intersections while
                // no edge crosses the rectangle
                // The algorithm will indefinitely switch directions in this situation
                // We must make sure to prevent this endless loop by returning no intersections.
                return vec![];
            }
            intersections.clear();
        }
    }

    fn is_handle_inside(&self, handle: DirectedEdgeHandle<V, DE, UE, F>) -> bool {
        is_handle_inside(self.lower, self.upper, handle)
    }
}

fn is_handle_inside<V, DE, UE, F>(
    lower: Point2<V::Scalar>,
    upper: Point2<V::Scalar>,
    next_handle: DirectedEdgeHandle<V, DE, UE, F>,
) -> bool
where
    V: HasPosition,
{
    is_inside(
        lower.to_f64(),
        upper.to_f64(),
        next_handle.from().position().to_f64(),
        next_handle.to().position().to_f64(),
    )
}

fn is_inside(lower: Point2<f64>, upper: Point2<f64>, from: Point2<f64>, to: Point2<f64>) -> bool {
    if is_point_inside(lower, upper, from) || is_point_inside(lower, upper, to) {
        return true;
    }

    let v0 = lower.to_f64();
    let v1 = Point2::new(lower.x, upper.y).to_f64();
    let v2 = upper.to_f64();
    let v3 = Point2::new(upper.x, lower.y).to_f64();

    let edges = [(v0, v1), (v1, v2), (v2, v3), (v3, v0)];

    for (v0, v1) in edges {
        let (s0, s1) = get_edge_intersections(v0, v1, from, to);
        if s0.is_infinite() {
            // Edges are colinear
            if !math::side_query(from, to, v0).is_on_line() {
                // Edges are parallel but not overlapping
                continue;
            }

            let p1 = math::project_point(from, to, v0);
            let p2 = math::project_point(from, to, v1);

            if p1.is_on_edge() || p2.is_on_edge() {
                return true;
            }
        } else if (0.0..=1.0).contains(&s0) && (0.0..=1.0).contains(&s1) {
            return true;
        }
    }
    false
}

pub(crate) fn is_point_inside<S: SpadeNum>(
    lower: Point2<S>,
    upper: Point2<S>,
    point: Point2<S>,
) -> bool {
    point.all_component_wise(lower, |a, b| a >= b) && point.all_component_wise(upper, |a, b| a <= b)
}

fn get_edge_intersections(
    v0: Point2<f64>,
    v1: Point2<f64>,
    e0: Point2<f64>,
    e1: Point2<f64>,
) -> (f64, f64) {
    let x4: f64 = v1.x;
    let x3: f64 = v0.x;
    let x2: f64 = e1.x;
    let x1: f64 = e0.x;
    let y4: f64 = v1.y;
    let y3: f64 = v0.y;
    let y2: f64 = e1.y;
    let y1: f64 = e0.y;
    let divisor = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
    let (s0, s1);
    if divisor == 0.0 {
        s0 = f64::INFINITY;
        s1 = f64::INFINITY;
    } else {
        s0 = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / divisor;
        s1 = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / divisor;
    }
    (s0, s1)
}

impl<'a, V, DE, UE, F> Iterator for EdgesInRectangleIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    type Item = UndirectedEdgeHandle<'a, V, DE, UE, F>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(next_edge) = self.remaining_hull_vec.pop() {
            if next_edge.rev().is_outer_edge() {
                if self.is_handle_inside(next_edge) {
                    return Some(next_edge.as_undirected());
                } else {
                    continue;
                }
            }

            if !next_edge.is_outer_edge() {
                if !self.remaining_hull_set.remove(&next_edge) {
                    // next_edge is a leftover that still resided in remaining_hull_vec - ignore it
                    continue;
                }

                if self.remaining_hull_set.remove(&next_edge.rev())
                    && self.is_handle_inside(next_edge)
                {
                    // next_edge is an "inner edge" of the current hull. Make sure both get removed and
                    // return it once.
                    return Some(next_edge.as_undirected());
                }
            }

            if !self.is_handle_inside(next_edge) {
                continue;
            }

            let next_edge = next_edge.rev();

            let next = next_edge.next();
            let prev = next_edge.prev();

            self.remaining_hull_vec.push(prev);
            self.remaining_hull_set.insert(prev);

            self.remaining_hull_vec.push(next);
            self.remaining_hull_set.insert(next);

            return Some(next_edge.as_undirected());
        }

        None
    }
}

#[cfg(test)]
mod test {
    use crate::{
        edges_in_rectangle_iterator::{is_handle_inside, is_inside, is_point_inside},
        test_utilities::random_points_with_seed,
        ConstrainedDelaunayTriangulation, DelaunayTriangulation, InsertionError, Point2,
        Triangulation,
    };

    use super::get_edge_intersections;

    #[test]
    fn test_empty() {
        let d = DelaunayTriangulation::<Point2<_>>::default();
        let edges = d.get_edges_in_rectangle(Point2::new(-1.0, -2.0), Point2::new(2.0, 2.0));
        assert_eq!(edges.count(), 0);
    }

    #[test]
    fn test_single_face() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(3.0, 2.0),
            Point2::new(-2.0, 1.0),
            Point2::new(-3.0, -4.0),
        ];
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        test_rectangle_iterator(&d);
        Ok(())
    }

    #[test]
    fn test_rectangle_center_on_vertex() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(3.0, 2.0),
            Point2::new(0.0, 0.0),
            Point2::new(-2.0, 1.0),
            Point2::new(-3.0, -4.0),
        ];
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        test_rectangle_iterator(&d);
        Ok(())
    }

    #[test]
    fn test_small() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(3.0, 2.0),
            Point2::new(-2.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(1.0, -4.0),
        ];
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        let edges = d.get_edges_in_rectangle(Point2::new(-10.0, -10.0), Point2::new(10.0, 10.0));
        assert_eq!(edges.count(), d.num_undirected_edges());

        Ok(())
    }

    #[test]
    fn test_random() -> Result<(), InsertionError> {
        for size in [4, 52, 122] {
            let vertices = random_points_with_seed(size, crate::test_utilities::SEED);
            let d = DelaunayTriangulation::<_>::bulk_load(vertices.clone())?;
            test_rectangle_iterator(&d);
            let mut c = ConstrainedDelaunayTriangulation::<_>::bulk_load(vertices)?;
            let constraints = random_points_with_seed(size, crate::test_utilities::SEED2);
            for points in constraints.as_slice().chunks_exact(2) {
                let from = c.insert(points[0])?;
                let to = c.insert(points[1])?;
                if c.can_add_constraint(from, to) {
                    c.add_constraint(from, to);
                }
            }
            test_rectangle_iterator(&c);
        }
        Ok(())
    }

    fn test_rectangle_iterator(d: &impl Triangulation<Vertex = Point2<f64>>) {
        let areas = [
            (Point2::new(-2.0, -2.0), Point2::new(-0.1, -0.1)),
            (Point2::new(-10.0, -10.0), Point2::new(10.0, 10.0)),
            (Point2::new(-0.5, -0.5), Point2::new(0.5, 0.5)),
            (Point2::new(-0.1, -10.), Point2::new(0.1, 10.)),
            (Point2::new(-5.0, -0.1), Point2::new(5.0, 0.1)),
            (Point2::new(-0.9, -0.9), Point2::new(0.9, 0.9)),
            (Point2::new(0.0, 0.0), Point2::new(0.0, 0.0)),
            (Point2::new(20.1, 20.1), Point2::new(10.0, 10.0)),
            (Point2::new(-2.0, -2.0), Point2::new(0.0, 0.0)),
            (Point2::new(-2.0, 0.0), Point2::new(0.0, 2.0)),
            (Point2::new(0.0, 0.0), Point2::new(2.0, 2.0)),
            (Point2::new(-2.0, -2.0), Point2::new(-1.0, -1.0)),
            (Point2::new(-2.0, -2.0), Point2::new(-1.0, -0.5)),
            (Point2::new(-2.0, -2.0), Point2::new(-1.0, 0.0)),
        ];

        for (lower, upper) in areas {
            let edges = d.get_edges_in_rectangle(lower, upper);
            let sequential_count = d
                .undirected_edges()
                .filter(|edge| is_handle_inside(lower, upper, edge.as_directed()))
                .count();
            assert_eq!(edges.count(), sequential_count);

            let vertices = d.get_vertices_in_rectangle(lower, upper);
            let sequential_count = d
                .vertices()
                .filter(|vertex| is_point_inside(lower, upper, vertex.position()))
                .count();
            assert_eq!(vertices.count(), sequential_count);
        }
    }

    #[test]
    fn test_medium_triangulation() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(-7., -5.5),
            Point2::new(-4., -6.5),
            Point2::new(-5., -9.),
            Point2::new(-6., 6.),
            Point2::new(-8., -6.),
            Point2::new(3., 3.),
        ];
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        test_rectangle_iterator(&d);
        Ok(())
    }

    #[test]
    fn test_convex_hull_edge_cases() -> Result<(), InsertionError> {
        let coords = [-1.0f64, 0.1, 1.0];
        let mut vertices = Vec::new();
        for x in &coords {
            for y in &coords {
                vertices.push(Point2::new(*x, *y));
            }
        }
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        test_rectangle_iterator(&d);
        Ok(())
    }

    #[test]
    fn test_degenerate() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.5),
            Point2::new(2.0, 1.0),
            Point2::new(4.0, 2.0),
        ];
        let d = DelaunayTriangulation::<_>::bulk_load(vertices)?;
        test_rectangle_iterator(&d);
        Ok(())
    }

    #[test]
    fn test_edge_intersections() {
        let sut = get_edge_intersections;
        let v0 = Point2::new(0.0, 0.0);
        let v1 = Point2::new(2.0, 0.0);

        let (s0, s1) = sut(v0, v1, Point2::new(1.0, 1.0), Point2::new(1.0, -1.0));
        assert_eq!(s0, 0.5);
        assert_eq!(s1, 0.5);

        let (s0, s1) = sut(v0, v1, Point2::new(0.5, 1.0), Point2::new(0.5, -1.0));
        assert_eq!(s0, 0.25);
        assert_eq!(s1, 0.5);

        let (s0, s1) = sut(v0, v1, v0, v1);
        assert!(s0.is_infinite());
        assert!(s1.is_infinite());

        let (s0, s1) = sut(v0, v1, Point2::new(1.0, 0.0), Point2::new(20.0, 0.0));
        assert!(s0.is_infinite());
        assert!(s1.is_infinite());

        let (s0, s1) = sut(v0, v1, Point2::new(1.0, 3.0), Point2::new(20.0, 3.0));
        assert!(s0.is_infinite());
        assert!(s1.is_infinite());
    }

    #[test]
    fn test_is_inside() {
        fn check(from: Point2<f64>, to: Point2<f64>, expected: bool) {
            let lower = Point2::new(-1.0, -1.0);
            let upper = Point2::new(1.0, 1.0);
            assert_eq!(is_inside(lower, upper, from, to), expected);
            assert_eq!(is_inside(lower, upper, to, from), expected);
        }

        check(Point2::new(0.0, 0.0), Point2::new(2.0, 3.0), true);
        check(Point2::new(-2.0, -3.0), Point2::new(2.0, 4.0), true);
        check(Point2::new(-2.0, 0.5), Point2::new(2.0, 0.5), true);

        check(Point2::new(-0.25, -0.25), Point2::new(0.1, 0.1), true);

        check(Point2::new(-0.25, -0.25), Point2::new(1.0, 1.0), true);
        check(Point2::new(-1.0, -1.0), Point2::new(0.5, 1.0), true);

        check(Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0), true);
        check(Point2::new(-1.0, -1.0), Point2::new(2.0, 2.0), true);

        check(Point2::new(1.0, 1.0), Point2::new(3.0, 4.0), true);

        check(Point2::new(-1.0, -1.0), Point2::new(-2.0, 0.0), true);

        check(Point2::new(-1.0, -1.0), Point2::new(-1.0, 1.0), true);
        check(Point2::new(-1.0, -1.0), Point2::new(-1.0, 2.0), true);
        check(Point2::new(-1.0, -2.0), Point2::new(-1.0, 2.0), true);
        check(Point2::new(-2.0, -1.0), Point2::new(-0.5, -1.0), true);

        check(Point2::new(-2.0, 2.0), Point2::new(3.0, -3.0), true);
        check(Point2::new(-0.25, 0.25), Point2::new(-20.0, 20.0), true);

        check(Point2::new(-2.0, 0.0), Point2::new(0.0, -2.0), true);

        check(Point2::new(-2.0, -2.0), Point2::new(-1.0, -3.0), false);
        check(Point2::new(-2.0, 0.0), Point2::new(0.0, -2.1), false);
        check(Point2::new(-1.0, -2.0), Point2::new(1.0, -2.0), false);

        check(Point2::new(-2.0, -1.0), Point2::new(-1.5, -1.0), false);
        check(Point2::new(1.5, -1.0), Point2::new(2.0, -1.0), false);
    }
}
