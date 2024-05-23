use crate::{HasPosition, Point2, Triangulation, TriangulationExt};
use crate::delaunay_core::math;
use crate::handles::{DirectedEdgeHandle, FixedVertexHandle, VertexHandle};

/// An iterator over all intersections of a straight line across the triangulation.
///
/// This iterator can, for example, be used to detect which edges prevent the insertion
/// of a new constraint edge.
///
/// Three different intersection kinds are possible:
///  - The line crosses a different edge
///  - The line goes through a vertex
///  - The line overlaps with an existing edge
///
/// These intersections are defined by the [Intersection] enum. Intersections are always returned in the same
/// order as the line, i.e. the closest intersection to the starting point is returned first.
///
/// # Example
///
/// This triangulation is created by the code below:
///
#[doc = include_str!("../images/line_intersection_iterator_scenario.svg")]
///
/// The line (s -> e) generates 5 intersections (labeled with `0..=4`): two edge intersections,
/// a vertex intersection, an edge overlap and  a final vertex intersection.
///
/// ```
/// use spade::{LineIntersectionIterator, Point2, DelaunayTriangulation};
/// # use spade::InsertionError;
/// # fn main() -> Result<(), InsertionError> {
/// let vertices = vec![
///     Point2::new(-30.0, 20.0),  // v0
///     Point2::new(0.0, -20.0),   // v1
///     Point2::new(0.0, 20.0),    // v2
///     Point2::new(14.0, 0.0),    // v3
///     Point2::new(30.0, 0.0),    // v4
/// ];
///
/// let triangulation = DelaunayTriangulation::<_>::bulk_load_stable(vertices)?;
/// for intersection in LineIntersectionIterator::new(
///     &triangulation,
///     Point2::new(-30.0, 0.0),
///     Point2::new(40.0, 0.0),
/// ) {
///     println!("{:?}", intersection);
/// }
/// # Ok(())
/// # }
/// ```
///
/// Output (simplified):
/// ```text
/// EdgeIntersection(DirectedEdgeHandle v0 -> v1)
/// EdgeIntersection(DirectedEdgeHandle v2 -> v1)
/// VertexIntersection(VertexHandle(3))
/// EdgeOverlap(DirectedEdgeHandle v3 -> v4)
/// VertexIntersection(VertexHandle(4))
/// ```
pub struct LineIntersectionIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
{
    cur_intersection: Option<Intersection<'a, V, DE, UE, F>>,
    line_from: Point2<V::Scalar>,
    line_to: Point2<V::Scalar>,
}

/// An intersection that can occur when moving through a triangulation along a straight line.
///
/// This is used as return type for [LineIntersectionIterator].
#[allow(clippy::enum_variant_names)]
pub enum Intersection<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    /// Indicates that the line is either crossing or touching an existing edge.
    /// The line's destination will always be either on the edge or on its left side (in a right-handed coordinate system).
    EdgeIntersection(DirectedEdgeHandle<'a, V, DE, UE, F>),
    /// Indicates that the line is touching a vertex.
    /// A line beginning or starting on a vertex also generates this intersection. A "line" beginning and starting on the same
    /// vertex will also return this intersection.
    VertexIntersection(VertexHandle<'a, V, DE, UE, F>),
    /// Indicates that a line is (partially) overlapping an existing edge.
    ///
    /// This implies that the line points in the same direction as the edge and that they share a common line segment.
    EdgeOverlap(DirectedEdgeHandle<'a, V, DE, UE, F>),
}

impl<'a, V, DE, UE, F> core::fmt::Debug for Intersection<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        use self::Intersection::*;
        match self {
            EdgeIntersection(handle) => write!(f, "EdgeIntersection({:?})", handle),
            VertexIntersection(handle) => write!(f, "VertexIntersection({:?})", handle),
            EdgeOverlap(handle) => write!(f, "EdgeOverlap({:?})", handle),
        }
    }
}

impl<'a, V, DE, UE, F> PartialEq for Intersection<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    fn eq(&self, other: &Self) -> bool {
        use self::Intersection::*;
        match (self, other) {
            (&EdgeIntersection(h0), &EdgeIntersection(h1)) => h0 == h1,
            (&VertexIntersection(h0), &VertexIntersection(h1)) => h0 == h1,
            (&EdgeOverlap(h0), &EdgeOverlap(h1)) => h0 == h1,
            _ => false,
        }
    }
}

impl<'a, V, DE, UE, F> Copy for Intersection<'a, V, DE, UE, F> where V: HasPosition {}

impl<'a, V, DE, UE, F> Clone for Intersection<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, V, DE, UE, F> Intersection<'a, V, DE, UE, F>
where
    V: HasPosition,
{
    /// Returns the intersected edge if this is an edge intersection or `None` otherwise.
    pub fn as_edge_intersection(&self) -> Option<DirectedEdgeHandle<'a, V, DE, UE, F>> {
        match self {
            Intersection::EdgeIntersection(ref edge) => Some(*edge),
            _ => None,
        }
    }
}

impl<'a, V, DE, UE, F> LineIntersectionIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
{
    /// Creates a new `LineIntersectionIterator` covering an arbitrary line.
    /// See [LineIntersectionIterator] for more information.
    pub fn new<T>(
        delaunay: &'a T,
        line_from: Point2<V::Scalar>,
        line_to: Point2<V::Scalar>,
    ) -> LineIntersectionIterator<'a, V, DE, UE, F>
    where
        T: Triangulation<Vertex = V, DirectedEdge = DE, UndirectedEdge = UE, Face = F>,
    {
        let first_intersection = Self::get_first_intersection(delaunay, line_from, line_to);
        LineIntersectionIterator {
            cur_intersection: first_intersection,
            line_from,
            line_to,
        }
    }

    /// Creates a new line iterator covering the line spanned by two existing vertices.
    ///
    /// Both start and end vertex are part of the iteration result.
    ///
    /// # Example
    /// ```
    /// use spade::{DelaunayTriangulation, Intersection, LineIntersectionIterator, Point2, Triangulation};
    /// # use spade::InsertionError;
    /// # fn main() -> Result<(), InsertionError> {
    ///     let mut triangulation = DelaunayTriangulation::<Point2<f64>>::new();
    ///     let v0 = triangulation.insert(Point2::new(0.0, 0.0))?;
    ///     let v1 = triangulation.insert(Point2::new(1.0, 1.0))?;
    ///
    ///     let expected_edge_overlap = triangulation.get_edge_from_neighbors(v0, v1).unwrap();
    ///     let all_intersections = LineIntersectionIterator::new_from_handles(&triangulation, v0, v1).collect::<Vec<_>>();
    ///     
    ///     let v0 = triangulation.vertex(v0);
    ///     let v1 = triangulation.vertex(v1);
    ///     assert_eq!(all_intersections, vec![
    ///         Intersection::VertexIntersection(v0),
    ///         Intersection::EdgeOverlap(expected_edge_overlap),
    ///         Intersection::VertexIntersection(v1),
    ///     ]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new_from_handles<T>(
        delaunay: &T,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
    ) -> LineIntersectionIterator<V, DE, UE, F>
    where
        T: Triangulation<Vertex = V, DirectedEdge = DE, UndirectedEdge = UE, Face = F>,
    {
        let from = delaunay.vertex(from);
        let line_from = from.position();
        let to = delaunay.vertex(to);
        let line_to = to.position();

        LineIntersectionIterator {
            cur_intersection: Some(Intersection::VertexIntersection(from)),
            line_from,
            line_to,
        }
    }

    fn get_first_intersection<T>(
        delaunay: &'a T,
        line_from: Point2<V::Scalar>,
        line_to: Point2<V::Scalar>,
    ) -> Option<Intersection<'a, V, DE, UE, F>>
    where
        T: Triangulation<Vertex = V, DirectedEdge = DE, UndirectedEdge = UE, Face = F>,
    {
        use crate::PositionInTriangulation::*;

        match delaunay.locate_with_hint_option_core(line_from, None) {
            OutsideOfConvexHull(edge_handle) => {
                let mut edge = delaunay.directed_edge(edge_handle);
                let line_from_query = edge.side_query(line_from);

                loop {
                    if line_from_query.is_on_line() {
                        let dist_from = edge.from().position().distance_2(line_from);
                        let dist_to = edge.to().position().distance_2(line_from);
                        let vertex = if dist_to < dist_from {
                            edge.to()
                        } else {
                            edge.from()
                        };
                        return Some(Intersection::VertexIntersection(vertex));
                    }
                    let line_to_query = edge.side_query(line_to);

                    if line_to_query.is_on_left_side() {
                        return None;
                    }

                    let edge_from = edge.from().position();
                    let edge_to = edge.to().position();
                    let edge_from_query = math::side_query(line_from, line_to, edge_from);
                    let edge_to_query = math::side_query(line_from, line_to, edge_to);

                    match (
                        edge_from_query.is_on_left_side(),
                        edge_to_query.is_on_left_side_or_on_line(),
                    ) {
                        (true, true) => edge = edge.prev(),
                        (false, false) => edge = edge.next(),
                        (false, true) => {
                            if edge_to_query.is_on_line() {
                                return Some(Intersection::VertexIntersection(edge.to()));
                            }
                            if edge_from_query.is_on_line() {
                                return Some(Intersection::VertexIntersection(edge.from()));
                            }
                            return Some(Intersection::EdgeIntersection(edge.rev()));
                        }
                        (true, false) => panic!("Unexpected edge topology. This is a bug."),
                    }
                }
            }
            OnFace(face_handle) => get_first_edge_from_edge_ring(
                delaunay.face(face_handle).adjacent_edges().iter().copied(),
                line_from,
                line_to,
            ),
            OnVertex(vertex_handle) => Some(Intersection::VertexIntersection(
                delaunay.vertex(vertex_handle),
            )),
            OnEdge(edge) => {
                let edge = delaunay.directed_edge(edge);
                let edge_from = edge.from().position();
                let edge_to = edge.to().position();
                let from_query = math::side_query(line_from, line_to, edge_from);
                let to_query = math::side_query(line_from, line_to, edge_to);
                if from_query.is_on_line() && to_query.is_on_line() {
                    let dist_to = edge_to.sub(line_to).length2();
                    let dist_from = edge_from.sub(line_to).length2();
                    if dist_to < dist_from {
                        Some(Intersection::EdgeOverlap(edge))
                    } else {
                        Some(Intersection::EdgeOverlap(edge.rev()))
                    }
                } else {
                    let edge_query = edge.side_query(line_to);
                    if edge_query.is_on_left_side() {
                        Some(Intersection::EdgeIntersection(edge))
                    } else {
                        Some(Intersection::EdgeIntersection(edge.rev()))
                    }
                }
            }
            NoTriangulation => {
                if let Some(next_vertex) = delaunay.vertices().next() {
                    let single_vertex = next_vertex.position();
                    let projection = math::project_point(line_from, line_to, single_vertex);
                    if projection.is_on_edge() {
                        let query = math::side_query(line_from, line_to, single_vertex);
                        if query.is_on_line() {
                            return Some(Intersection::VertexIntersection(next_vertex));
                        }
                    }
                }
                None
            }
        }
    }

    fn get_next(&mut self) -> Option<Intersection<'a, V, DE, UE, F>> {
        use self::Intersection::*;
        match self.cur_intersection {
            Some(EdgeIntersection(cur_edge)) => {
                match trace_direction_out_of_edge(cur_edge, self.line_from, self.line_to) {
                    EdgeOutDirection::ConvexHull => None,
                    EdgeOutDirection::VertexIntersection(vertex) => {
                        Some(VertexIntersection(vertex))
                    }
                    EdgeOutDirection::EdgeIntersection(edge) => Some(EdgeIntersection(edge)),
                    EdgeOutDirection::NoIntersection => None,
                }
            }
            Some(VertexIntersection(vertex)) => {
                if vertex.position() == self.line_to {
                    // Target point was reached - the iteration can finish
                    None
                } else {
                    match trace_direction_out_of_vertex(vertex, self.line_to) {
                        VertexOutDirection::ConvexHull => None,
                        VertexOutDirection::EdgeOverlap(edge) => Some(EdgeOverlap(edge)),
                        VertexOutDirection::EdgeIntersection(edge) => {
                            if edge.side_query(self.line_to).is_on_right_side() {
                                // The target point was skipped over - the iteration can finish
                                None
                            } else {
                                Some(EdgeIntersection(edge))
                            }
                        }
                    }
                }
            }
            Some(EdgeOverlap(edge)) => {
                if self.line_from == self.line_to {
                    None
                } else if math::project_point(self.line_from, self.line_to, edge.to().position())
                    .is_on_edge()
                {
                    Some(VertexIntersection(edge.to()))
                } else {
                    None
                }
            }
            None => None,
        }
    }
}

impl<'a, V, DE, UE, F> Iterator for LineIntersectionIterator<'a, V, DE, UE, F>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
{
    type Item = Intersection<'a, V, DE, UE, F>;

    fn next(&mut self) -> Option<Self::Item> {
        let cur = self.cur_intersection;
        self.cur_intersection = self.get_next();
        cur
    }
}

fn get_first_edge_from_edge_ring<'a, I, V, DE, UE, F>(
    ring: I,
    line_from: Point2<V::Scalar>,
    line_to: Point2<V::Scalar>,
) -> Option<Intersection<'a, V, DE, UE, F>>
where
    I: IntoIterator<Item = DirectedEdgeHandle<'a, V, DE, UE, F>>,
    V: HasPosition,
{
    use self::Intersection::*;
    for edge in ring {
        let cur_from = edge.from().position();
        let cur_to = edge.to().position();

        debug_assert!(math::side_query(cur_from, cur_to, line_from).is_on_left_side_or_on_line());
        if math::intersects_edge_non_collinear(line_from, line_to, cur_from, cur_to) {
            if math::side_query(line_from, line_to, cur_from).is_on_line() {
                return Some(VertexIntersection(edge.from()));
            } else if math::side_query(line_from, line_to, cur_to).is_on_line() {
                return Some(VertexIntersection(edge.to()));
            }
            return Some(EdgeIntersection(edge.rev()));
        }
    }
    None
}

pub(super) enum VertexOutDirection<'a, V, DE, UE, F> {
    ConvexHull,
    EdgeOverlap(DirectedEdgeHandle<'a, V, DE, UE, F>),
    EdgeIntersection(DirectedEdgeHandle<'a, V, DE, UE, F>),
}

pub(super) fn trace_direction_out_of_vertex<V, DE, UE, F>(
    vertex: VertexHandle<V, DE, UE, F>,
    line_to: Point2<V::Scalar>,
) -> VertexOutDirection<V, DE, UE, F>
where
    V: HasPosition,
{
    let mut current_edge = match vertex.out_edge() {
        Some(edge) => edge,
        None => return VertexOutDirection::ConvexHull,
    };

    let mut current_query = current_edge.side_query(line_to);

    let iterate_ccw = current_query.is_on_left_side();

    loop {
        if current_query.is_on_line() && !current_edge.project_point(line_to).is_before_edge() {
            return VertexOutDirection::EdgeOverlap(current_edge);
        }

        let next = if iterate_ccw {
            current_edge.ccw()
        } else {
            current_edge.cw()
        };

        let next_query = next.side_query(line_to);
        if next_query.is_on_line() && !next.project_point(line_to).is_before_edge() {
            return VertexOutDirection::EdgeOverlap(next);
        }

        let face_between_current_and_next = if iterate_ccw {
            current_edge.face()
        } else {
            next.face()
        };
        if face_between_current_and_next.is_outer() {
            return VertexOutDirection::ConvexHull;
        }

        if iterate_ccw == next_query.is_on_right_side() {
            let segment_edge = if iterate_ccw {
                current_edge.next()
            } else {
                current_edge.rev().prev()
            };
            return VertexOutDirection::EdgeIntersection(segment_edge.rev());
        }

        current_query = next_query;
        current_edge = next;
    }
}

pub(super) enum EdgeOutDirection<'a, V, DE, UE, F> {
    ConvexHull,
    VertexIntersection(VertexHandle<'a, V, DE, UE, F>),
    EdgeIntersection(DirectedEdgeHandle<'a, V, DE, UE, F>),
    NoIntersection,
}

pub(super) fn trace_direction_out_of_edge<V, DE, UE, F>(
    edge: DirectedEdgeHandle<V, DE, UE, F>,
    line_from: Point2<V::Scalar>,
    line_to: Point2<V::Scalar>,
) -> EdgeOutDirection<V, DE, UE, F>
where
    V: HasPosition,
{
    debug_assert!(
        edge.side_query(line_to).is_on_left_side_or_on_line(),
        "The target must be on the left side of the current edge"
    );

    let e_prev = if edge.is_outer_edge() {
        // Iteration reached an edge of the convex hull, we're done.
        return EdgeOutDirection::ConvexHull;
    } else {
        edge.prev()
    };

    let o_next = edge.next();

    // Find out which edges of the left face intersect the line
    let e_prev_inter = e_prev.intersects_edge_non_collinear(line_from, line_to);
    let o_next_inter = o_next.intersects_edge_non_collinear(line_from, line_to);

    match (e_prev_inter, o_next_inter) {
        (true, false) => EdgeOutDirection::EdgeIntersection(e_prev.rev()),
        (false, true) => EdgeOutDirection::EdgeIntersection(o_next.rev()),
        (true, true) => {
            // Both edges intersect - this means the line is cutting through a common point
            EdgeOutDirection::VertexIntersection(e_prev.from())
        }
        (false, false) => EdgeOutDirection::NoIntersection,
    }
}

#[cfg(test)]
mod test {
    use self::Intersection::*;
    use super::*;
    use crate::{InsertionError, Point2, Triangulation as _};

    use alloc::{vec, vec::Vec};

    type Triangulation = crate::DelaunayTriangulation<Point2<f64>>;

    fn reverse<'a, V, DE, UE, F>(
        intersection: &Intersection<'a, V, DE, UE, F>,
    ) -> Intersection<'a, V, DE, UE, F>
    where
        V: HasPosition,
    {
        match intersection {
            EdgeIntersection(edge) => EdgeIntersection(edge.rev()),
            VertexIntersection(vertex) => VertexIntersection(*vertex),
            EdgeOverlap(edge) => EdgeOverlap(edge.rev()),
        }
    }

    fn check(
        delaunay: &Triangulation,
        from: Point2<f64>,
        to: Point2<f64>,
        mut expected: Vec<Intersection<Point2<f64>, (), (), ()>>,
    ) {
        let collected: Vec<_> = LineIntersectionIterator::new(delaunay, from, to).collect();
        assert_eq!(collected, expected);
        let mut reversed = Vec::new();
        let rev_collected: Vec<_> = LineIntersectionIterator::new(delaunay, to, from).collect();
        for intersection in rev_collected.iter() {
            reversed.push(reverse(intersection));
        }
        expected.reverse();
        assert_eq!(reversed, expected);
    }

    fn create_test_triangulation() -> Result<
        (
            Triangulation,
            FixedVertexHandle,
            FixedVertexHandle,
            FixedVertexHandle,
            FixedVertexHandle,
        ),
        InsertionError,
    > {
        let v0 = Point2::new(-2.0, -2.0);
        let v1 = Point2::new(2.0, 2.0);
        let v2 = Point2::new(1.0, -1.0);
        let v3 = Point2::new(-1.0, 1.0);

        let mut delaunay = Triangulation::new();
        let v0 = delaunay.insert(v0)?;
        let v1 = delaunay.insert(v1)?;
        let v2 = delaunay.insert(v2)?;
        let v3 = delaunay.insert(v3)?;

        Ok((delaunay, v0, v1, v2, v3))
    }

    #[test]
    fn test_single_line_intersection() -> Result<(), InsertionError> {
        let (delaunay, _, _, v2, v3) = create_test_triangulation()?;
        let from = Point2::new(-0.5, -0.5);
        let to = Point2::new(0.5, 0.5);
        let edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        check(&delaunay, from, to, vec![EdgeIntersection(edge)]);
        Ok(())
    }

    #[test]
    fn test_empty_inner_intersection() -> Result<(), InsertionError> {
        let (delaunay, _, _, _, _) = create_test_triangulation()?;
        let from = Point2::new(-0.5, -0.5);
        let to = Point2::new(-0.25, -0.25);
        assert!(LineIntersectionIterator::new(&delaunay, from, to)
            .next()
            .is_none());

        assert!(LineIntersectionIterator::new(&delaunay, from, from)
            .next()
            .is_none());
        Ok(())
    }

    #[test]
    fn test_between_vertices_intersection() -> Result<(), InsertionError> {
        let (delaunay, v0, v1, v2, v3) = create_test_triangulation()?;
        let from = Point2::new(-2.0, -2.0);
        let to = Point2::new(2.0, 2.0);
        let edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        let first = VertexIntersection(delaunay.vertex(v0));
        let last = VertexIntersection(delaunay.vertex(v1));
        let edges: Vec<_> = LineIntersectionIterator::new(&delaunay, from, to).collect();
        assert_eq!(edges, vec![first, EdgeIntersection(edge), last]);
        Ok(())
    }

    #[test]
    fn test_mixed_intersections() -> Result<(), InsertionError> {
        let (mut delaunay, _, v1, v2, v3) = create_test_triangulation()?;
        let v4 = delaunay.insert(Point2::new(1.0, 1.0))?;
        let from = Point2::new(-1.0, -1.0);
        let to = Point2::new(2.0, 2.0);
        let intersection_edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        let overlap_edge = delaunay.get_edge_from_neighbors(v4, v1).unwrap();
        check(
            &delaunay,
            from,
            to,
            vec![
                EdgeIntersection(intersection_edge),
                VertexIntersection(delaunay.vertex(v4)),
                EdgeOverlap(overlap_edge),
                VertexIntersection(delaunay.vertex(v1)),
            ],
        );
        Ok(())
    }

    #[test]
    fn test_out_of_hull_intersections() -> Result<(), InsertionError> {
        let (ref d, v0, v1, v2, v3) = create_test_triangulation()?;

        let edge20 = d.get_edge_from_neighbors(v2, v0).unwrap();
        let edge20 = EdgeIntersection(edge20);
        let edge30 = d.get_edge_from_neighbors(v3, v0).unwrap();
        let edge30 = EdgeIntersection(edge30);
        let edge12 = d.get_edge_from_neighbors(v1, v2).unwrap();
        let edge12 = EdgeIntersection(edge12);
        let edge32 = d.get_edge_from_neighbors(v3, v2).unwrap();
        let o32 = EdgeOverlap(edge32);
        let edge32 = EdgeIntersection(edge32);

        let v0 = VertexIntersection(d.vertex(v0));
        let v2 = VertexIntersection(d.vertex(v2));
        let v3 = VertexIntersection(d.vertex(v3));

        // No intersection
        let from = Point2::new(-2.0, 1.0);
        let to = Point2::new(-2.0, 0.0);
        check(d, from, to, vec![]);
        // One intersection
        let from = Point2::new(-2., 0.);
        let to = Point2::new(-2., -4.0);
        check(d, from, to, vec![v0]);
        let from = Point2::new(-0.5, -0.5);
        check(d, from, to, vec![edge20]);
        // Two intersections
        let from = Point2::new(-2.0, 0.0);
        let to = Point2::new(0., -2.);
        check(d, from, to, vec![edge30, edge20]);
        let from = Point2::new(-3.0, 3.0);
        let to = Point2::new(3.0, -3.0);
        check(d, from, to, vec![v3, o32, v2]);
        // Three intersections
        let from = Point2::new(-2.0, 0.0);
        let to = Point2::new(2., -1.);
        check(d, from, to, vec![edge30, edge32, edge12]);
        Ok(())
    }

    #[test]
    fn test_on_line_intersection() -> Result<(), InsertionError> {
        let (d, _, v1, v2, v3) = create_test_triangulation()?;

        let edge = d.get_edge_from_neighbors(v2, v3).unwrap();
        let e32 = EdgeIntersection(edge.rev());
        let o23 = EdgeOverlap(edge);
        let o32 = EdgeOverlap(edge.rev());

        let v1 = VertexIntersection(d.vertex(v1));
        let v2 = VertexIntersection(d.vertex(v2));
        let v3 = VertexIntersection(d.vertex(v3));

        let from = Point2::new(0.0, 0.0);
        let to = Point2::new(0.0, 0.0);
        let collected: Vec<_> = LineIntersectionIterator::new(&d, from, to).collect();
        assert!(collected == vec![o23] || collected == vec![o32]);
        let to = Point2::new(0.2, 0.2);
        check(&d, from, to, vec![e32]);
        let to = Point2::new(2.0, 2.0);
        check(&d, from, to, vec![e32, v1]);
        let to = Point2::new(-30.0, 30.0);
        check(&d, from, to, vec![o23, v3]);
        let to = Point2::new(30.0, -30.0);
        check(&d, from, to, vec![o32, v2]);
        let from = Point2::new(-30.0, 30.0);
        check(&d, from, to, vec![v3, o32, v2]);
        Ok(())
    }

    #[test]
    fn test_intersecting_zero_vertices() {
        let delaunay = Triangulation::new();
        let mut iterator = LineIntersectionIterator::new(
            &delaunay,
            Point2::new(0.5, 1.234),
            Point2::new(3.223, 42.0),
        );
        assert!(iterator.next().is_none());
    }

    #[test]
    fn test_intersection_when_passing_equal_line_from_and_line_to() -> Result<(), InsertionError> {
        let (delaunay, v1, _, v3, v4) = create_test_triangulation()?;
        let v1_position = delaunay.s().vertex(v1).position();
        let edge = delaunay.get_edge_from_neighbors(v3, v4).unwrap();
        let v1 = delaunay.s().vertex(v1);

        check(
            &delaunay,
            v1_position,
            v1_position,
            vec![VertexIntersection(v1)],
        );
        let origin = Point2::new(0.0, 0.0);
        let intersections: Vec<_> =
            LineIntersectionIterator::new(&delaunay, origin, origin).collect();

        assert_eq![intersections.len(), 1];
        assert!(
            intersections[0] == EdgeOverlap(edge) || intersections[0] == EdgeOverlap(edge.rev())
        );
        Ok(())
    }

    #[test]
    fn test_intersecting_single_vertex() -> Result<(), InsertionError> {
        let mut delaunay = Triangulation::new();
        let pos = Point2::new(0.5, 0.5);
        let v0 = delaunay.insert(pos)?;
        let v0 = delaunay.vertex(v0);
        let from = Point2::new(1.0, 0.0);
        let to = Point2::new(0.0, 1.0);
        check(&delaunay, from, to, vec![VertexIntersection(v0)]);
        check(&delaunay, pos, pos, vec![VertexIntersection(v0)]);
        let to = Point2::new(1.234, 42.0);
        check(&delaunay, from, to, vec![]);
        Ok(())
    }

    #[test]
    fn test_intersecting_degenerate_triangulation() -> Result<(), InsertionError> {
        let mut d = Triangulation::new();

        let v2 = d.insert(Point2::new(0.0, 0.0))?;
        let v3 = d.insert(Point2::new(1.0, 1.0))?;
        let v1 = d.insert(Point2::new(-2.0, -2.0))?;
        let v0 = d.insert(Point2::new(-3.0, -3.0))?;

        let e01 = d.get_edge_from_neighbors(v0, v1).unwrap();
        let e12 = d.get_edge_from_neighbors(v1, v2).unwrap();
        let e23 = d.get_edge_from_neighbors(v2, v3).unwrap();

        let v0 = VertexIntersection(d.vertex(v0));
        let v1 = VertexIntersection(d.vertex(v1));
        let v2 = VertexIntersection(d.vertex(v2));
        let v3 = VertexIntersection(d.vertex(v3));

        let intersection_23 = EdgeIntersection(e23);
        let e01 = EdgeOverlap(e01);
        let e12 = EdgeOverlap(e12);
        let e23 = EdgeOverlap(e23);

        let from = Point2::new(-4.0, -4.0);
        let to = Point2::new(5.0, 5.0);
        check(&d, from, to, vec![v0, e01, v1, e12, v2, e23, v3]);
        let to = Point2::new(-2.0, -2.0);
        check(&d, from, to, vec![v0, e01, v1]);
        let to = Point2::new(-0.5, -0.5);
        check(&d, from, to, vec![v0, e01, v1, e12]);

        let from = Point2::new(1.0, -0.0);
        let to = Point2::new(0.0, 1.0);
        check(&d, from, to, vec![intersection_23]);
        let to = Point2::new(0.5, 0.5);
        check(&d, from, to, vec![intersection_23]);

        let from = Point2::new(0.5, -0.5);
        let to = Point2::new(-0.5, 0.5);
        check(&d, from, to, vec![v2]);

        let single_vertex = Point2::new(-2.0, -2.0);
        check(&d, single_vertex, single_vertex, vec![v1]);

        Ok(())
    }

    #[test]
    fn test_intersecting_through_point_ending_on_face() -> Result<(), InsertionError> {
        let mut d = Triangulation::new();

        let v0 = d.insert(Point2::new(0.0, 0.0))?;
        let v1 = d.insert(Point2::new(1.0, 1.0))?;
        let v2 = d.insert(Point2::new(-1.0, 1.0))?;
        let v3 = d.insert(Point2::new(0.0, 2.0))?;
        d.insert(Point2::new(-1.0, 3.0))?;
        d.insert(Point2::new(1.0, 3.0))?;

        let e = d.get_edge_from_neighbors(v2, v1).unwrap();

        let v0 = VertexIntersection(d.vertex(v0));
        let e = EdgeIntersection(e);
        let v3 = VertexIntersection(d.vertex(v3));

        let from = Point2::new(0.0, 0.0);
        let to = Point2::new(0.0, 2.5);
        check(&d, from, to, vec![v0, e, v3]);
        Ok(())
    }
}
