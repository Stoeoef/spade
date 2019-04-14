use crate::primitives::SimpleEdge;
use crate::delaunay::*;
use self::delaunay_basic::BasicDelaunaySubdivision;
use crate::traits::{HasPosition2D};
use crate::point_traits::{PointN, TwoDimensional, PointNExtensions};

pub struct LineIntersectionIterator<'a, T, V, E=()> where
    T: BasicDelaunaySubdivision<V, EdgeType=E> + 'a,
    V: HasPosition2D + 'a,
    V::Point: TwoDimensional,
    E: Default + Copy + 'a,
{
    cur_intersection: Option<Intersection<'a, V, E>>,
    line: SimpleEdge<V::Point>,
    delaunay: &'a T,
}

pub enum Intersection<'a, V, E=()> where
    V: 'a, E: 'a
{
    EdgeIntersection(EdgeHandle<'a, V, E>),
    VertexIntersection(VertexHandle<'a, V, E>),
    EdgeOverlap(EdgeHandle<'a, V, E>),
}

impl <'a, V, E> ::std::fmt::Debug for Intersection<'a, V, E> where
    V: 'a,
    E: Default,
{
     fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
         use self::Intersection::*;
         match self {
             &EdgeIntersection(handle) => write!(f, "EdgeIntersection({:?})", handle),
             &VertexIntersection(handle) => write!(f, "VertexIntersection({:?})", handle),
             &EdgeOverlap(handle) => write!(f, "EdgeOverlap({:?})", handle),
         }
    }
}

impl <'a, V, E> PartialEq for Intersection<'a, V, E> {
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

impl <'a, V, E> Copy for Intersection<'a, V, E> {}

impl <'a, V, E> Clone for Intersection<'a, V, E> {
    fn clone(&self) -> Self {
        use self::Intersection::*;
        match self {
            &EdgeIntersection(handle) => EdgeIntersection(handle),
            &VertexIntersection(handle) => VertexIntersection(handle),
            &EdgeOverlap(handle) => EdgeOverlap(handle),
        }
    }
}

impl <'a, T, V, E>  LineIntersectionIterator<'a, T, V, E> where
    T: BasicDelaunaySubdivision<V, EdgeType=E> + 'a,
    V: HasPosition2D + 'a,
    V::Point: TwoDimensional,
    E: Default + Copy,
{
    pub fn new(delaunay: &'a T, from: &V::Point, to: &V::Point) -> LineIntersectionIterator<'a, T, V, E> {
        let line = SimpleEdge::new(from.clone(), to.clone());
        let first_intersection = Self::get_first_intersection(delaunay, &line);
        LineIntersectionIterator {
            cur_intersection: first_intersection,
            line: line,
            delaunay: delaunay,
        }
    }

    pub fn new_from_handles(delaunay: &'a T, from: FixedVertexHandle, to: FixedVertexHandle)
        -> LineIntersectionIterator<'a, T, V, E> {
        let from = delaunay.s().vertex(from);
        let to = delaunay.s().vertex(to);
        let line = SimpleEdge::new(from.position(), to.position());
        LineIntersectionIterator {
            cur_intersection: Some(Intersection::VertexIntersection(from)),
            line: line,
            delaunay: delaunay,
        }
    }

    fn get_first_intersection(delaunay: &'a T, line: &SimpleEdge<V::Point>) -> Option<Intersection<'a, V, E>> {
        use crate::delaunay::PositionInTriangulation::*;
        match delaunay.locate_with_hint_option(&line.from, None) {
            InTriangle(face_handle) => {
                Self::get_first_edge_from_edge_ring(face_handle.adjacent_edges(), line)
            },
            OnPoint(vertex_handle) => {
                Some(Intersection::VertexIntersection(vertex_handle))
            },
            OnEdge(edge) => {
                let edge_from = edge.from().position();
                let edge_to = edge.to().position();
                let from_query = line.side_query::<T::Kernel>(&edge_from);
                let to_query = line.side_query::<T::Kernel>(&edge_to);
                if from_query.is_on_line() && to_query.is_on_line() {
                    let dist_to = edge_to.sub(&line.to).length2();
                    let dist_from = edge_from.sub(&line.to).length2();
                    if dist_to < dist_from {
                        Some(Intersection::EdgeOverlap(edge))
                    } else {
                        Some(Intersection::EdgeOverlap(edge.sym()))
                    }
                } else {
                    Some(Intersection::EdgeIntersection(edge))
                }
            }
            OutsideConvexHull(edge) => {
                let query = T::to_simple_edge(edge).side_query::<T::Kernel>(&line.from);
                if query.is_on_line() {
                    let dist_to = edge.to().position().sub(&line.from).length2();
                    let dist_from = edge.from().position().sub(&line.from).length2();
                    let vertex = if dist_to < dist_from {
                        edge.to()
                    } else {
                        edge.from()
                    };
                    Some(Intersection::VertexIntersection(vertex))
                } else {
                    let edges = delaunay.get_convex_hull_edges_for_point(edge.fix(), &line.from);
                    Self::get_first_edge_from_edge_ring( 
                        edges.iter().map(|e| delaunay.s().edge(*e)), 
                        line)
                }
            },
            NoTriangulationPresent => {
                if delaunay.s().num_vertices() == 0 {
                    None
                } else {
                    let vertex = delaunay.s().vertices().next().unwrap();
                    if line.side_query::<T::Kernel>(&vertex.position()).is_on_line() {
                        if line.is_projection_on_edge(&vertex.position()) {
                            Some(Intersection::VertexIntersection(vertex))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                }
            }
        }
    }

    fn get_first_edge_from_edge_ring<I>(ring: I, line: &SimpleEdge<V::Point>) -> Option<Intersection<'a, V, E>> where
        I: IntoIterator<Item=EdgeHandle<'a, V, E>>,
    {
        use self::Intersection::*;
        for edge in ring {
            let cur_edge = SimpleEdge::new(edge.from().position(),
                                           edge.to().position());

            debug_assert!(cur_edge.side_query::<T::Kernel>(&line.from).is_on_left_side_or_on_line());
            if line.intersects_edge_non_collinear::<T::Kernel>(&cur_edge) {
                if line.side_query::<T::Kernel>(&cur_edge.from).is_on_line() {
                    return Some(VertexIntersection(edge.from()));
                } else if line.side_query::<T::Kernel>(&cur_edge.to).is_on_line() {
                    return Some(VertexIntersection(edge.to()));
                }
                return Some(EdgeIntersection(edge.sym()));
            }
        }
        None
    }

    fn get_next(&mut self) -> Option<Intersection<'a, V, E>> {
        use self::Intersection::*;
        match self.cur_intersection {
            Some(EdgeIntersection(cur_edge)) => {
                if cur_edge.face() == self.delaunay.infinite_face() {
                    // The line has reached the convex hull
                    None
                } else {
                    let ref target = self.line.to;
                    assert!(T::to_simple_edge(cur_edge)
                            .side_query::<T::Kernel>(&target).is_on_left_side_or_on_line(),
                            "The target must be on the left side of the current edge");

                    let e_prev = cur_edge.o_prev();
                    let o_next = cur_edge.o_next();
                    // Find out which edges of the left face intersect the line
                    let e_prev_inter = T::to_simple_edge(e_prev)
                        .intersects_edge_non_collinear::<T::Kernel>(&self.line);
                    let o_next_inter = T::to_simple_edge(o_next)
                        .intersects_edge_non_collinear::<T::Kernel>(&self.line);
                    match (e_prev_inter, o_next_inter) {
                        (true, false) => Some(EdgeIntersection(e_prev.sym())),
                        (false, true) => Some(EdgeIntersection(cur_edge.o_next().sym())),
                        (true, true) => {
                            // Both edges intersect - this mean line is cutting through a common point
                            Some(VertexIntersection(e_prev.from()))
                        },
                        (false, false) => None,
                    }
                }
            },
            Some(VertexIntersection(vertex)) => {
                if vertex.position() == self.line.to {
                    // Target point was reached - the iteration can finish
                    None
                } else {
                    let line = SimpleEdge::new(vertex.position(), self.line.to.clone());
                    for edge in vertex.ccw_out_edges() {
                        let pos = edge.to().position();
                        if line.side_query::<T::Kernel>(&pos).is_on_line() {
                            let p1 = &edge.from().position();
                            let p2 = &edge.to().position();
                            let dir = p2.sub(p1);
                            let s = self.line.to.sub(p1).dot(&dir);
                            if s > ::num::zero::<<V::Point as PointN>::Scalar>() {
                                return Some(EdgeOverlap(edge));
                            }
                        }
                    }
                    let ring_iterator = vertex.ccw_out_edges()
                        .filter_map(|e| {
                            let edge = e.o_next();
                            if edge.face() == self.delaunay.infinite_face() {
                                None
                            } else {
                                Some(edge)
                            }
                        });
                    Self::get_first_edge_from_edge_ring(ring_iterator, &line)
                }
            },
            Some(EdgeOverlap(edge)) => {
                if self.line.from == self.line.to {
                    None
                } else {
                    if self.line.is_projection_on_edge(&edge.to().position()) {
                        Some(VertexIntersection(edge.to()))
                    } else {
                        None
                    }
                }
            },
            None => None
        }
    }
}

impl <'a, T, V, E> Iterator for LineIntersectionIterator<'a, T, V, E> 
    where 
          T: BasicDelaunaySubdivision<V, EdgeType=E> + 'a,
          V: HasPosition2D + 'a,
          V::Point: TwoDimensional,
          E: Default + Copy + 'a,
{
    type Item = Intersection<'a, V, E>;

    fn next(&mut self) -> Option<Self::Item> {
        let cur = self.cur_intersection;
        self.cur_intersection = self.get_next();
        cur
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use cgmath::Point2;
    use crate::kernels::FloatKernel;
    use self::Intersection::*;

    type Triangulation = DelaunayTriangulation<Point2<f64>, FloatKernel>;

    fn reverse<'a, V>(intersection: &Intersection<'a, V>) -> Intersection<'a, V> {
        match intersection {
            &EdgeIntersection(edge) => EdgeIntersection(edge.sym()),
            &VertexIntersection(vertex) => VertexIntersection(vertex),
            &EdgeOverlap(edge) => EdgeOverlap(edge.sym()),
        }
    }

    fn check(delaunay: &Triangulation, 
             from: Point2<f64>,
             to: Point2<f64>,
             mut expected: Vec<Intersection<Point2<f64>>>) {
        let collected: Vec<_> = LineIntersectionIterator::new(delaunay, &from, &to).collect();
        assert_eq!(collected, expected);
        let mut reversed = Vec::new();
        let rev_collected: Vec<_> = LineIntersectionIterator::new(delaunay, &to, &from).collect();
        for intersection in rev_collected.iter() {
            reversed.push(reverse(intersection));
        }
        expected.reverse();
        assert_eq!(reversed, expected);
    }

    fn create_test_triangulation() -> (Triangulation, 
                                       FixedVertexHandle, 
                                       FixedVertexHandle, 
                                       FixedVertexHandle, 
                                       FixedVertexHandle) {
        let mut delaunay = Triangulation::new();
        let v0 = delaunay.insert(Point2::new(-2.0, -2.0));
        let v1 = delaunay.insert(Point2::new(2.0, 2.0));
        let v2 = delaunay.insert(Point2::new(1.0, -1.0));
        let v3 = delaunay.insert(Point2::new(-1.0, 1.0));
        (delaunay, v0, v1, v2, v3)
    }
    
    #[test]
    fn test_single_line_intersection() {
        let (delaunay, _, _, v2, v3) = create_test_triangulation();
        let from = Point2::new(-0.5, -0.5);
        let to = Point2::new(0.5, 0.5);
        let edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        check(&delaunay, from, to, vec![EdgeIntersection(edge)]);
    }

    #[test]
    fn test_empty_inner_intersection() {
        let (delaunay, _, _, _, _) = create_test_triangulation();
        let from = Point2::new(-0.5, -0.5);
        let to = Point2::new(-0.25, -0.25);
        assert!(LineIntersectionIterator::new(&delaunay, &from, &to).next().is_none());
    }

    #[test]
    fn test_between_vertices_intersection() {
        let (delaunay, v0, v1, v2, v3) = create_test_triangulation();
        let from = Point2::new(-2.0, -2.0);
        let to = Point2::new(2.0, 2.0);
        let edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        let first = VertexIntersection(delaunay.vertex(v0));
        let last = VertexIntersection(delaunay.vertex(v1));
        let edges: Vec<_> = LineIntersectionIterator::new(&delaunay, &from, &to).collect();
        assert_eq!(edges, vec![first, EdgeIntersection(edge), last]);
    }

    #[test]
    fn test_mixed_intersections() {
        let (mut delaunay, _, v1, v2, v3) = create_test_triangulation();
        let v4 = delaunay.insert(Point2::new(1.0, 1.0));
        let from = Point2::new(-1.0, -1.0);
        let to = Point2::new(2.0, 2.0);
        let intersection_edge = delaunay.get_edge_from_neighbors(v3, v2).unwrap();
        let overlap_edge = delaunay.get_edge_from_neighbors(v4, v1).unwrap();
        check(&delaunay, from, to, vec![
            EdgeIntersection(intersection_edge),
            VertexIntersection(delaunay.vertex(v4)),
            EdgeOverlap(overlap_edge),
            VertexIntersection(delaunay.vertex(v1))]);
    }

    #[test]
    fn test_out_of_hull_intersections() {
        let (ref d, v0, v1, v2, v3) = create_test_triangulation();

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
    }

    #[test]
    fn test_on_line_intersection() {
        let (d, _, v1, v2, v3) = create_test_triangulation();

        let edge = d.get_edge_from_neighbors(v2, v3).unwrap();
        let e32 = EdgeIntersection(edge.sym());
        let o23 = EdgeOverlap(edge);
        let o32 = EdgeOverlap(edge.sym());

        let v1 = VertexIntersection(d.vertex(v1));
        let v2 = VertexIntersection(d.vertex(v2));
        let v3 = VertexIntersection(d.vertex(v3));

        let from = Point2::new(0.0, 0.0);
        let to = Point2::new(0.0, 0.0);
        let collected: Vec<_> = LineIntersectionIterator::new(&d, &from, &to).collect();
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
    }

    #[test]
    fn test_degenerate_iteration() {
        let mut t = Triangulation::new();
        let v0 = t.insert(Point2::new(0.0, 0.0));

        let from = Point2::new(-1.0, -1.0);
        let to = Point2::new(-0.5, -0.5);
        check(&t, from, to, vec![]);

        let to = Point2::new(1.0, 1.0);
        check(&t, from, to, vec![VertexIntersection(t.vertex(v0))]);

        let v2 = t.insert(Point2::new(2.0, 2.0));
        let v1 = t.insert(Point2::new(1.0, 1.0));
        assert!(t.is_degenerate());
        let to = Point2::new(3.0, 3.0);

        let e01 = t.get_edge_from_neighbors(v0, v1).unwrap();
        let e12 = t.get_edge_from_neighbors(v1, v2).unwrap();

        let o01 = EdgeOverlap(e01);
        let o12 = EdgeOverlap(e12);

        let v0 = VertexIntersection(t.vertex(v0));
        let v1 = VertexIntersection(t.vertex(v1));
        let v2 = VertexIntersection(t.vertex(v2));
        
        check(&t, from, to, vec![v0, o01, v1, o12, v2]);
        let from = Point2::new(1.0, 0.0);
        let to = Point2::new(0.0, 1.0);
        check(&t, from, to, vec![EdgeIntersection(e01)]);

        let from = Point2::new(0.0, 2.0);
        let to = Point2::new(2.0, 0.0);
        check(&t, from, to, vec!(v1));

        let from = Point2::new(1.0, 1.0);
        let to = Point2::new(1.5, 1.5);
        check(&t, from, to, vec![v1, EdgeOverlap(e12)])
    }
}
