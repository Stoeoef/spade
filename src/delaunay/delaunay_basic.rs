use delaunay::*;
use self::dcel::*;
use traits::HasPosition2D;
use smallvec::SmallVec;
use kernels::DelaunayKernel;
use point_traits::{PointN, TwoDimensional};
use primitives::SimpleEdge;
use std::collections::HashSet;

pub trait Subdivision<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional
{
    /// Creates a dynamic vertex handle from a fixed vertex handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V>;

    /// Returns a mutable reference to the vertex data referenced by a
    /// `FixedVertexHandle`.
    fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V;

    /// Creates a dynamic face handle from a fixed face handle.
    /// May panic if the faces was invalidated by a previous vertex
    /// removal.
    fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V>;

    /// Creates a dynamic edge handle from a fixed edge handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V>;

    /// Returns the number of vertices in this triangulation.
    fn num_vertices(&self) -> usize;

    /// Returns the number of faces in this triangulation.
    /// *Note*: This count does include the infinite face.
    fn num_faces(&self) -> usize;

    /// Returns the number of triangles in this triangulation.
    /// As there is always exactly one face not a triangle, this is
    /// `self.num_faces() - 1`.
    fn num_triangles(&self) -> usize;

    /// Returns the number of edges in this triangulation.
    fn num_edges(&self) -> usize;

    /// Returns an iterator over all triangles.
    fn triangles(&self) -> FacesIterator<V>;

    /// Returns an iterator over all edges.
    fn edges(&self) -> EdgesIterator<V>;

    /// Returns an iterator over all vertices.
    fn vertices(&self) -> VerticesIterator<V>;

    fn get_edge_from_vertices(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> Option<EdgeHandle<V>>;
}

pub trait Locateable<V, K>: Subdivision<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    /// Returns information about the location of a point in a triangulation.
    ///
    /// Additionally, a hint can be given to speed up computation. The hint should be a vertex close
    /// to the position that is being looked up.
    fn locate_with_hint
        (&self,
         point: &V::Point,
         hint: FixedVertexHandle)
         -> PositionInTriangulation<VertexHandle<V>, FaceHandle<V>, EdgeHandle<V>>;

    fn locate_vertex_with_hint(&self,
                               point: &V::Point,
                               hint: FixedVertexHandle)
                               -> Option<VertexHandle<V>> {
        match self.locate_with_hint(point, hint) {
            PositionInTriangulation::OnPoint(h) => Some(h),
            _ => None,
        }
    }
}

pub trait HasSubdivision<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    fn s(&self) -> &DCEL<V>;
    fn s_mut(&mut self) -> &mut DCEL<V>;
}

impl<T, V, K> Subdivision<V, K> for T
    where T: HasSubdivision<V, K>,
          V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    /// Creates a dynamic vertex handle from a fixed vertex handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V> {
        self.s().vertex(handle)
    }

    /// Returns a mutable reference to the vertex data referenced by a
    /// `FixedVertexHandle`.
    fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V {
        self.s_mut().vertex_mut(handle)
    }

    /// Creates a dynamic face handle from a fixed face handle.
    /// May panic if the faces was invalidated by a previous vertex
    /// removal.
    fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V> {
        self.s().face(handle)
    }

    /// Creates a dynamic edge handle from a fixed edge handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V> {
        self.s().edge(handle)
    }

    /// Returns the number of vertices in this triangulation.
    fn num_vertices(&self) -> usize {
        self.s().num_vertices()
    }

    /// Returns the number of faces in this triangulation.
    /// *Note*: This count does include the infinite face.
    fn num_faces(&self) -> usize {
        self.s().num_faces()
    }

    /// Returns the number of triangles in this triangulation.
    /// As there is always exactly one face not a triangle, this is
    /// `self.num_faces() - 1`.
    fn num_triangles(&self) -> usize {
        self.s().num_faces() - 1
    }

    /// Returns the number of edges in this triangulation.
    fn num_edges(&self) -> usize {
        self.s().num_edges()
    }

    /// Returns an iterator over all triangles.
    fn triangles(&self) -> FacesIterator<V> {
        let mut result = self.s().faces();
        // Skip the outer face
        result.next();
        result
    }

    /// Returns an iterator over all edges.
    fn edges(&self) -> EdgesIterator<V> {
        self.s().edges()
    }

    /// Returns an iterator over all vertices.
    fn vertices(&self) -> VerticesIterator<V> {
        self.s().vertices()
    }

    fn get_edge_from_vertices(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> Option<EdgeHandle<V>> {
        from_neighbors(self.s(), from, to)
    }
}
impl<T, V, K> Locateable<V, K> for T
    where T: Subdivision<V, K> + BasicDelaunaySubdivision<V, K>,
          V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    fn locate_with_hint
        (&self,
         point: &V::Point,
         hint: FixedVertexHandle)
         -> PositionInTriangulation<VertexHandle<V>, FaceHandle<V>, EdgeHandle<V>> {
        use self::PositionInTriangulation::*;
        match self.locate_with_hint_fixed(point, hint) {
            NoTriangulationPresent => NoTriangulationPresent,
            InTriangle(face) => InTriangle(self.s().face(face)),
            OutsideConvexHull(edge) => OutsideConvexHull(self.s().edge(edge)),
            OnPoint(vertex) => OnPoint(self.s().vertex(vertex)),
            OnEdge(edge) => OnEdge(self.s().edge(edge)),
        }
    }
}

pub trait BasicDelaunaySubdivision<V, K>: HasSubdivision<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{

    fn is_defined_legal(&self, _: FixedEdgeHandle) -> bool {
        false
    }

    fn insert_outside_convex_hull(&mut self,
                                  closest_edge: FixedEdgeHandle,
                                  t: V)
                                  -> FixedVertexHandle {
        let position = t.position();
        let mut ch_edges = self.get_convex_hull_edges_for_point(closest_edge, &position);
        let new_handle = self.s_mut().insert_vertex(t);
        // Make new connections
        let mut last_edge =
            self.s_mut().connect_edge_to_isolated_vertex(*ch_edges.last().unwrap(), new_handle);

        for edge in ch_edges.iter().rev() {
            last_edge = self.s_mut().create_face(last_edge, *edge);
            // Reverse last_edge
            last_edge = self.s_mut()
                .edge(last_edge)
                .sym()
                .fix();
        }
        self.legalize_edges(&mut ch_edges, &position);
        new_handle
    }

    fn insert_into_triangle(&mut self, face: FixedFaceHandle, t: V) -> FixedVertexHandle {
        let position = t.position();

        let new_handle = self.s_mut().insert_vertex(t);
        let (e0, e1, e2) = {
            let face = self.s_mut().face(face);
            let adj = face.adjacent_edge().unwrap();
            (adj.o_prev().fix(), adj.fix(), adj.o_next().fix())
        };

        let mut last_edge = self.s_mut().connect_edge_to_isolated_vertex(e2, new_handle);
        last_edge = self.s_mut()
            .edge(last_edge)
            .sym()
            .fix();
        last_edge = self.s_mut().create_face(e0, last_edge);
        last_edge = self.s_mut()
            .edge(last_edge)
            .sym()
            .fix();
        self.s_mut().create_face(e1, last_edge);
        let mut edges = SmallVec::new();
        edges.push(e0);
        edges.push(e1);
        edges.push(e2);
        self.legalize_edges(&mut edges, &position);
        new_handle
    }

    fn insert_on_edge(&mut self, edge: FixedEdgeHandle, t: V) -> FixedVertexHandle {
        let position = t.position();
        let new_handle = self.s_mut().insert_vertex(t);
        let mut illegal_edges = SmallVec::new();
        let (from, to) = {
            let edge = self.s_mut().edge(edge);
            (edge.from().fix(), edge.to().fix())
        };

        let left_handle_opt = self.get_left_triangle(from, to);
        let right_handle_opt = self.get_right_triangle(from, to);
        let edge_handle = from_neighbors(&self.s_mut(), from, to).unwrap().fix();
        self.s_mut().split_edge(edge_handle, new_handle);
        if let Some(left_handle) = left_handle_opt {
            let edge1 = from_neighbors(&self.s_mut(), to, left_handle).unwrap().fix();
            let edge0 = from_neighbors(&self.s_mut(), left_handle, from).unwrap().fix();
            let edge_mid = from_neighbors(&self.s_mut(), from, new_handle).unwrap().fix();

            self.s_mut().create_face(edge_mid, edge0);
            illegal_edges.push(edge0);
            illegal_edges.push(edge1);
        }
        if let Some(right_handle) = right_handle_opt {
            let edge0 = from_neighbors(&self.s_mut(), from, right_handle).unwrap().fix();
            let edge1 = from_neighbors(&self.s_mut(), right_handle, to).unwrap().fix();
            let edge_mid = from_neighbors(&self.s_mut(), to, new_handle).unwrap().fix();
            self.s_mut().create_face(edge_mid, edge1);
            illegal_edges.push(edge0);
            illegal_edges.push(edge1);
        }
        self.legalize_edges(&mut illegal_edges, &position);
        new_handle
    }



    fn get_convex_hull_edges_for_point(&self,
                                       first_edge: FixedEdgeHandle,
                                       point: &V::Point)
                                       -> SmallVec<[FixedEdgeHandle; 16]>
        where V: HasPosition2D,
              V::Point: TwoDimensional,
              K: DelaunayKernel<<V::Point as PointN>::Scalar>
    {
        let mut result = SmallVec::new();
        let first_edge = self.s().edge(first_edge);
        debug_assert!(K::side_query(&Self::to_simple_edge(first_edge), point).is_on_left_side());

        let mut last_edge = first_edge;
        result.push(last_edge.fix());
        // Follow the first edge in cw and ccw direction
        loop {
            last_edge = last_edge.o_next();
            let query = K::side_query(&Self::to_simple_edge(last_edge), point);
            if query.is_on_left_side() {
                result.push(last_edge.fix());
            } else {
                break;
            }
        }

        last_edge = first_edge;
        loop {
            last_edge = last_edge.o_prev();
            let query = K::side_query(&Self::to_simple_edge(last_edge), point);
            if query.is_on_left_side() {
                result.insert(0, last_edge.fix());
            } else {
                break;
            }
        }
        result
    }

    fn to_simple_edge<'a>(edge: EdgeHandle<'a, V>) -> SimpleEdge<V::Point> {
        let from = (edge.from()).position();
        let to = (edge.to()).position();
        SimpleEdge::new(from, to)
    }

    fn legalize_edges(&mut self, edges: &mut SmallVec<[FixedEdgeHandle; 16]>, position: &V::Point)
        where V: HasPosition2D,
              V::Point: TwoDimensional,
              K: DelaunayKernel<<V::Point as PointN>::Scalar>
    {
        while let Some(e) = edges.pop() {
            if !self.is_ch_edge(e) && !self.is_defined_legal(e) {
                let (v0, v1, v2, e1, e2);
                {
                    let edge = self.s_mut().edge(e);
                    v0 = (*edge.from()).position();
                    v1 = (*edge.to()).position();
                    v2 = (*edge.cw().to()).position();
                    e1 = edge.sym().o_next().fix();
                    e2 = edge.sym().o_prev().fix();
                }
                debug_assert!(K::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(K::is_ordered_ccw(position, &v0, &v1));
                if K::contained_in_circumference(&v1, &v2, &v0, position) {
                    // The edge is illegal
                    self.s_mut().flip_cw(e);
                    edges.push(e1);
                    edges.push(e2);
                }
            }
        }
    }

    fn is_ch_edge(&self, edge: FixedEdgeHandle) -> bool {
        let edge = self.s().edge(edge);
        let sym = edge.sym();
        edge.face().fix() == 0 || sym.face().fix() == 0
    }

    fn get_left_triangle(&self,
                         from: FixedVertexHandle,
                         to: FixedVertexHandle)
                         -> Option<FixedVertexHandle> {
        let edge_handle = from_neighbors(&self.s(), from, to).unwrap();
        let ccw_handle = edge_handle.ccw().to();
        let query = K::side_query(&Self::to_simple_edge(edge_handle),
                                  &(*ccw_handle).position());
        if query.is_on_left_side() {
            debug_assert!(from_neighbors(&self.s(), ccw_handle.fix(), to).is_some());
            Some(ccw_handle.fix())
        } else {
            None
        }
    }

    fn get_right_triangle(&self,
                          from: FixedVertexHandle,
                          to: FixedVertexHandle)
                          -> Option<FixedVertexHandle> {
        self.get_left_triangle(to, from)
    }

    fn locate_with_hint_fixed
        (&self,
         point: &V::Point,
         start: FixedVertexHandle)
         -> PositionInTriangulation<FixedVertexHandle, FixedFaceHandle, FixedEdgeHandle> {
        let mut cur_edge = self.s()
            .vertex(start)
            .out_edge()
            .expect("Cannot start search with an isolated vertex");
        let mut cur_query = K::side_query(&Self::to_simple_edge(cur_edge), point);
        // Invariant: point must not be on the right side of cur_edge
        if cur_query.is_on_right_side() {
            cur_edge = cur_edge.sym();
            cur_query = cur_query.reversed();
        }
        loop {
            assert!(cur_query.is_on_left_side_or_on_line());
            if cur_edge.face().fix() == 0 {
                if cur_query.is_on_line() {
                    cur_edge = cur_edge.sym();
                } else {
                    return PositionInTriangulation::OutsideConvexHull(cur_edge.fix());
                }
            }
            let from_pos = (*cur_edge.from()).position();
            if &from_pos == point {
                return PositionInTriangulation::OnPoint(cur_edge.from().fix());
            }
            // Check if cur_edge.o_next is also on the left side
            let next = cur_edge.o_next();
            let to_pos = (*next.from()).position();
            if &to_pos == point {
                return PositionInTriangulation::OnPoint(cur_edge.to().fix());
            }

            let next_query = K::side_query(&Self::to_simple_edge(next), point);
            if next_query.is_on_right_side_or_on_line() {
                // We continue walking into the face right of next
                cur_edge = next.sym();
                cur_query = next_query.reversed();
            } else {
                // Check if cur_edge.o_prev is also on the left side
                let prev = cur_edge.o_prev();
                let prev_query = K::side_query(&Self::to_simple_edge(prev), point);
                if prev_query.is_on_right_side_or_on_line() {
                    // We continue walking into the face right of prev
                    cur_edge = prev.sym();
                    cur_query = prev_query.reversed();
                } else {
                    if cur_query.is_on_line() {
                        return PositionInTriangulation::OnEdge(cur_edge.fix());
                    }
                    return PositionInTriangulation::InTriangle(cur_edge.face().fix());
                }
            }
        }
    }

    fn fill_hole(&mut self, loop_edges: Vec<FixedEdgeHandle>) {
        let mut border_edges = HashSet::new();
        
        for e in &loop_edges {
            border_edges.insert(*e);
            border_edges.insert(self.s().edge(*e).sym().fix());
        }

        let last_edge = *loop_edges.last().unwrap();

        // Fill the hole
        let mut todo = Vec::new();
        for i in 2 .. loop_edges.len() - 1 {
            let edge = self.s_mut().create_face(last_edge, loop_edges[i]);
            todo.push(edge);
        }
        // Legalize edges
        while let Some(fixed_edge_handle) = todo.pop() {
            let (v0, v1, vl, vr, e1, e2, e3, e4);
            {
                let edge = self.s().edge(fixed_edge_handle);
                v0 = (*edge.from()).position();
                v1 = (*edge.to()).position();
                vl = (*edge.ccw().to()).position();
                vr = (*edge.cw().to()).position();
                e1 = edge.cw().fix();
                e2 = edge.ccw().fix();
                e3 = edge.sym().cw().fix();
                e4 = edge.sym().ccw().fix();
            }
            if !K::contained_in_circumference(&v0, &v1, &vl, &vr) {
                // Flip edge
                self.s_mut().flip_cw(fixed_edge_handle);
                
                for e in &[e1, e2, e3, e4] {
                    if !border_edges.contains(e) {
                        todo.push(*e);
                    }
                }
            }
        }
    }
}
