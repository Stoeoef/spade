use crate::delaunay::*;
use self::dcel::*;
use self::delaunay_locate::VertexEntry;
use crate::traits::HasPosition2D;
use smallvec::SmallVec;
use crate::kernels::DelaunayKernel;
use crate::point_traits::{PointN, PointNExtensions, TwoDimensional};
use crate::primitives::SimpleEdge;
use std::collections::HashSet;

type FixedPosition = PositionInTriangulation<FixedVertexHandle, FixedFaceHandle, FixedEdgeHandle>;
type DynamicPosition<'a, V, E> = PositionInTriangulation<VertexHandle<'a, V, E>,
                                                         FaceHandle<'a, V, E>,
                                                         EdgeHandle<'a, V, E>>;

pub trait Subdivision<V>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
{
    type Kernel: DelaunayKernel<<V::Point as PointN>::Scalar>;
    type EdgeType: Default;

    /// Creates a dynamic vertex handle from a fixed vertex handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V, Self::EdgeType>;

    /// Returns a mutable reference to the vertex data referenced by a
    /// `FixedVertexHandle`.
    fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V;

    /// Creates a dynamic face handle from a fixed face handle.
    /// May panic if the faces was invalidated by a previous vertex
    /// removal.
    fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V, Self::EdgeType>;

    /// Creates a dynamic edge handle from a fixed edge handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V, Self::EdgeType>;

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
    fn triangles(&self) -> FacesIterator<V, Self::EdgeType>;

    /// Returns an iterator over all edges.
    fn edges(&self) -> EdgesIterator<V, Self::EdgeType>;

    /// Returns an iterator over all vertices.
    fn vertices(&self) -> VerticesIterator<V, Self::EdgeType>;
}

pub trait HasSubdivision<V>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
{
    type Kernel: DelaunayKernel<<V::Point as PointN>::Scalar>;
    type EdgeType: Default + Copy;

    fn s(&self) -> &DCEL<V, Self::EdgeType>;
    fn s_mut(&mut self) -> &mut DCEL<V, Self::EdgeType>;
}

impl<T, V> Subdivision<V> for T
    where T: HasSubdivision<V>,
          V: HasPosition2D,
          V::Point: TwoDimensional,
{
    type Kernel = T::Kernel;
    type EdgeType = T::EdgeType;

    /// Creates a dynamic vertex handle from a fixed vertex handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V, Self::EdgeType> {
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
    fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V, Self::EdgeType> {
        self.s().face(handle)
    }

    /// Creates a dynamic edge handle from a fixed edge handle.
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V, Self::EdgeType> {
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
    fn triangles(&self) -> FacesIterator<V, Self::EdgeType> {
        let mut result = self.s().faces();
        // Skip the outer face
        result.next();
        result
    }

    /// Returns an iterator over all edges.
    fn edges(&self) -> EdgesIterator<V, Self::EdgeType> {
        self.s().edges()
    }

    /// Returns an iterator over all vertices.
    fn vertices(&self) -> VerticesIterator<V, Self::EdgeType> {
        self.s().vertices()
    }
}

// This trait implements most common operations on regular delaunay triangulations and cdts.
pub trait BasicDelaunaySubdivision<V>: HasSubdivision<V>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
{
    type LocateStructure: DelaunayLocateStructure<V::Point>;

    fn locate_structure(&self) -> &Self::LocateStructure;
    fn locate_structure_mut(&mut self) -> &mut Self::LocateStructure;

    fn all_points_on_line(&self) -> bool;
    fn set_all_points_on_line(&mut self, new_value: bool);

    fn is_defined_legal(&self, _: FixedEdgeHandle) -> bool {
        false
    }

    fn infinite_face(&self) -> FaceHandle<V, Self::EdgeType> {
        self.s().face(0)
    }

    fn handle_legal_edge_split(&mut self, _edges: &[FixedEdgeHandle; 4]) { }

    fn initial_insertion(&mut self, position: FixedPosition, t: V) -> 
        Result<FixedVertexHandle, FixedVertexHandle> 
    {
        assert!(self.all_points_on_line());
        // Inserts points if no points are present or if all points
        // lie on the same line
        match position {
            PositionInTriangulation::NoTriangulationPresent => {
                // Insert the first or second point
                let new_handle = self.s_mut().insert_vertex(t);
                if self.s().num_vertices() == 2 {
                    let (handle0, handle1) = {
                        let mut iter = self.s().fixed_vertices();
                        (iter.next().unwrap(), iter.next().unwrap())
                    };
                    self.s_mut().connect_two_isolated_vertices(handle0, handle1, 0);
                }
                Ok(new_handle)
            },
            PositionInTriangulation::OnEdge(edge) => {
                let new_vertex = self.s_mut().insert_vertex(t);
                self.s_mut().split_edge(edge, new_vertex);
                Ok(new_vertex)
            },
            PositionInTriangulation::OnPoint(vertex) => {
                self.s_mut().update_vertex(vertex, t);
                Result::Err(vertex)
            },
            PositionInTriangulation::OutsideConvexHull(edge) => {
                let line = Self::to_simple_edge(self.s().edge(edge));
                let pos = t.position();
                let query = line.side_query::<Self::Kernel>(&pos);
                if query.is_on_line() {
                    let prev_edge = {
                        let edge = self.s().edge(edge);
                        let dist_to = pos.sub(&edge.to().position()).length2();
                        let dist_from = pos.sub(&edge.from().position()).length2();
                        if dist_to < dist_from {
                            edge.fix()
                        } else {
                            edge.sym().fix()
                        }
                    };
                    let new_handle = self.s_mut().insert_vertex(t);
                    self.s_mut().connect_edge_to_isolated_vertex(prev_edge, new_handle);
                    Ok(new_handle)
                } else {
                    // The new vertex does not lay on the same line as the 
                    // points before
                    self.set_all_points_on_line(false);
                    let oriented_edge = if query.is_on_right_side() {
                        // insert_outside_convex_hull requires that the new point
                        // is to the left of the given edge
                        self.s().edge(edge).sym().fix()
                    } else {
                        edge
                    };
                    Ok(self.insert_outside_convex_hull(oriented_edge, t))
                        
                }
            },
            PositionInTriangulation::InTriangle(_) => {
                panic!("Impossible control flow path. This is a bug.");
            }
        }
    }

    fn get_default_hint(&self, coord: &V::Point) -> FixedVertexHandle {
        let hint = self.locate_structure().find_close_handle(coord);
        if hint < self.s().num_vertices() {
            hint
        } else {
            0
        }
    }

    fn insert_with_hint_option(&mut self, 
                               t: V,
                               hint: Option<FixedVertexHandle>)
                               -> FixedVertexHandle {
        let pos = t.position();
        let position_in_triangulation = self.locate_with_hint_option_fixed(&pos, hint);
        let insertion_result = if self.all_points_on_line() {
            self.initial_insertion(position_in_triangulation, t)
        } else {
            match position_in_triangulation {
                PositionInTriangulation::OutsideConvexHull(edge) => {
                    Result::Ok(self.insert_outside_convex_hull(edge, t))
                },
                PositionInTriangulation::InTriangle(face) => {
                    Result::Ok(self.insert_into_triangle(face, t))
                },
                PositionInTriangulation::OnEdge(edge) => {
                    if self.is_defined_legal(edge) {
                        // If the edge is defined as legal the resulting edges must
                        // be redefined as legal
                        let (from, to);
                        {
                            let edge = self.s().edge(edge);
                            from = edge.from().fix();
                            to = edge.to().fix();
                        }
                        let new_handle = self.insert_on_edge(edge, t);
                        let handles = {
                            let e1 = self.s().get_edge_from_neighbors(from, new_handle).unwrap();
                            let e2 = self.s().get_edge_from_neighbors(new_handle, to).unwrap();
                            [e1.fix(), e1.sym().fix(), e2.fix(), e2.sym().fix()]
                        };
                        self.handle_legal_edge_split(&handles);
                        Result::Ok(new_handle)
                    } else {
                        Result::Ok(self.insert_on_edge(edge, t))
                    }
                },
                PositionInTriangulation::OnPoint(vertex) => {
                    self.s_mut().update_vertex(vertex, t);
                    Result::Err(vertex)
                },
                PositionInTriangulation::NoTriangulationPresent => {
                    panic!("Impossible control path. This is a bug");
                }
            }
        };
        match insertion_result {
            Result::Ok(new_handle) => {
                self.locate_structure_mut().insert_vertex_entry(VertexEntry {
                    point: pos,
                    handle: new_handle
                });
                new_handle
            },
            Result::Err(update_handle) => {
                update_handle
            }
        }
    }

    fn locate_with_hint_option(&self, point: &V::Point, hint: Option<FixedVertexHandle>) 
                               -> DynamicPosition<V, Self::EdgeType> 
    {
        use self::PositionInTriangulation::*;
        match self.locate_with_hint_option_fixed(point, hint) {
            NoTriangulationPresent => NoTriangulationPresent,
            InTriangle(face) => InTriangle(self.s().face(face)),
            OutsideConvexHull(edge) => OutsideConvexHull(self.s().edge(edge)),
            OnPoint(vertex) => OnPoint(self.s().vertex(vertex)),
            OnEdge(edge) => OnEdge(self.s().edge(edge)),
        }
    }

    fn locate_with_hint_option_fixed(&self, point: &V::Point, hint: Option<FixedVertexHandle>) -> 
        PositionInTriangulation<FixedVertexHandle, FixedFaceHandle, FixedEdgeHandle> {
        if self.all_points_on_line() {
            self.brute_force_locate(point)
        } else {
            let start = hint.unwrap_or_else(|| self.get_default_hint(point));
            self.locate_with_hint_fixed(point, start) 
        }
    }

    fn brute_force_locate(&self, point: &V::Point)
                          -> FixedPosition {
        assert!(self.all_points_on_line());
        for vertex in self.s().vertices() {
            if &vertex.position() == point {
                return PositionInTriangulation::OnPoint(vertex.fix());
            }
        }
        if self.s().num_vertices() <= 1 {
            return PositionInTriangulation::NoTriangulationPresent;
        }

        let mut line = Self::to_simple_edge(self.s().edges().next().unwrap());
        let mut query = line.side_query::<Self::Kernel>(point);
        if query.is_on_right_side() {
            // Always make sure the point is on the left side of the query line
            query = query.reversed();
            ::std::mem::swap(&mut line.from, &mut line.to);
        }
        let dir = line.to.sub(&line.from);
        let mut vertices: Vec<_> = self.s().vertices().map(
            |v| (v.fix(), dir.dot(&(*v).position()))).collect();
        // Sort vertices according to their position on the line
        vertices.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());

        let get_edge = |index: usize| {
            let n1 = vertices[index - 1].0;
            let n2 = vertices[index].0;
            self.s().get_edge_from_neighbors(n1, n2).unwrap().fix()
        };
        let vertex_pos = dir.dot(point);
        match vertices.binary_search_by(|pos| pos.1.partial_cmp(&vertex_pos).unwrap()) {
            Ok(index) => {
                if query.is_on_line() {
                    // This should not be reachable
                    PositionInTriangulation::OnPoint(vertices[index].0)
                } else {
                    PositionInTriangulation::OutsideConvexHull(get_edge(1))
                }
            },
            Err(index) => {
                let len = vertices.len();
                if index == 0 {
                    PositionInTriangulation::OutsideConvexHull(get_edge(1))
                } else if index == len {
                    PositionInTriangulation::OutsideConvexHull(get_edge(len - 1))
                } else {
                    let edge = get_edge(index);
                    if query.is_on_line() {
                        PositionInTriangulation::OnEdge(edge)
                    } else {
                        PositionInTriangulation::OutsideConvexHull(edge)
                    }
                }
            }
        }
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
        let edge_handle = self.s().get_edge_from_neighbors(from, to).unwrap().fix();
        self.s_mut().split_edge(edge_handle, new_handle);
        if let Some(left_handle) = left_handle_opt {
            let edge1 = self.s().get_edge_from_neighbors(to, left_handle).unwrap().fix();
            let edge0 = self.s().get_edge_from_neighbors(left_handle, from).unwrap().fix();
            let edge_mid = self.s().get_edge_from_neighbors(from, new_handle).unwrap().fix();

            self.s_mut().create_face(edge_mid, edge0);
            illegal_edges.push(edge0);
            illegal_edges.push(edge1);
        }
        if let Some(right_handle) = right_handle_opt {
            let edge0 = self.s().get_edge_from_neighbors(from, right_handle).unwrap().fix();
            let edge1 = self.s().get_edge_from_neighbors(right_handle, to).unwrap().fix();
            let edge_mid = self.s().get_edge_from_neighbors(to, new_handle).unwrap().fix();
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
    {
        let mut result = SmallVec::new();
        let first_edge = self.s().edge(first_edge);
        debug_assert!(Self::Kernel::side_query(&Self::to_simple_edge(first_edge), point).is_on_left_side());

        let mut last_edge = first_edge;
        result.push(last_edge.fix());
        // Follow the first edge in cw and ccw direction
        loop {
            last_edge = last_edge.o_next();
            let query = Self::Kernel::side_query(&Self::to_simple_edge(last_edge), point);
            if query.is_on_left_side() {
                result.push(last_edge.fix());
            } else {
                break;
            }
        }

        last_edge = first_edge;
        loop {
            last_edge = last_edge.o_prev();
            let query = Self::Kernel::side_query(&Self::to_simple_edge(last_edge), point);
            if query.is_on_left_side() {
                result.insert(0, last_edge.fix());
            } else {
                break;
            }
        }
        result
    }

    fn to_simple_edge(edge: EdgeHandle<V, Self::EdgeType>) -> SimpleEdge<V::Point> {
        let from = (edge.from()).position();
        let to = (edge.to()).position();
        SimpleEdge::new(from, to)
    }

    fn legalize_edges(&mut self, edges: &mut SmallVec<[FixedEdgeHandle; 16]>, position: &V::Point)
        where V: HasPosition2D,
              V::Point: TwoDimensional,
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
                debug_assert!(Self::Kernel::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(Self::Kernel::is_ordered_ccw(position, &v0, &v1));
                if Self::Kernel::contained_in_circumference(&v1, &v2, &v0, position) {
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
        let edge_handle = self.s().get_edge_from_neighbors(from, to).unwrap();
        let ccw_handle = edge_handle.ccw().to();
        let query = Self::Kernel::side_query(&Self::to_simple_edge(edge_handle),
                                  &(*ccw_handle).position());
        if query.is_on_left_side() {
            debug_assert!(self.s().get_edge_from_neighbors(ccw_handle.fix(), to).is_some());
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
        let mut cur_query = Self::Kernel::side_query(&Self::to_simple_edge(cur_edge), point);
        // Invariant: point must not be on the right side of cur_edge
        if cur_query.is_on_right_side() {
            cur_edge = cur_edge.sym();
            cur_query = cur_query.reversed();
        }
        loop {
            assert!(cur_query.is_on_left_side_or_on_line());
            if cur_edge.face() == self.infinite_face() {
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

            let next_query = Self::Kernel::side_query(&Self::to_simple_edge(next), point);
            if next_query.is_on_right_side_or_on_line() {
                // We continue walking into the face right of next
                cur_edge = next.sym();
                cur_query = next_query.reversed();
            } else {
                // Check if cur_edge.o_prev is also on the left side
                let prev = cur_edge.o_prev();
                let prev_query = Self::Kernel::side_query(&Self::to_simple_edge(prev), point);
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

    fn remove(&mut self, vertex: FixedVertexHandle) -> V {
        let mut neighbors = Vec::new();
        let mut ch_removal = false;
        // Check if we remove from the convex hull
        // At the same time, all neigboring edges are collected
        for edge in self.s().vertex(vertex).ccw_out_edges() {
            if edge.face() == self.infinite_face() {
                ch_removal = true;
                neighbors.clear();
                let next = edge.ccw();
                for edge in next.ccw_iter() {
                    neighbors.push(edge.to().fix());
                }
                break;
            }
            neighbors.push(edge.to().fix());
        }
        ch_removal |= neighbors.is_empty();
        
        let infinite = self.infinite_face().fix();

        let VertexRemovalResult { updated_vertex, data } = 
            self.s_mut().remove_vertex(vertex, Some(infinite));

        // Remove point from locate structure
        let vertex_pos = data.position();
        self.locate_structure_mut().remove_vertex_entry(&VertexEntry::new(vertex_pos, vertex));

        if let Some(updated_vertex) = updated_vertex {
            // Update locate structure if necessary
            let pos = (*self.s().vertex(vertex)).position();
            self.locate_structure_mut().update_vertex_entry(VertexEntry {
                point: pos,
                handle: vertex,
            });
            // Update a vertex entry in the neighbors vector
            for n in &mut neighbors {
                if *n == updated_vertex {
                    *n = vertex;
                    break;
                }
            }
        }
        
        if ch_removal {
            if self.all_points_on_line() {
                if neighbors.len() > 1 {
                    self.repair_edge(&neighbors);
                }
            } else {
                // We removed a vertex from the convex hull
                self.repair_convex_hull(&neighbors);
                if self.s().num_faces() == 1 {
                    self.set_all_points_on_line(true);
                }
            }
        } else {
            // Removed an inner vertex
            let loop_edges: Vec<_> = {
                let first = self.s().get_edge_from_neighbors(
                    neighbors[0], neighbors[1]).unwrap();
                first.o_next_iterator().map(|e| e.fix()).collect()
            };
            self.fill_hole(loop_edges);
        }
        data
    }

    fn repair_edge(&mut self, vertices: &[FixedVertexHandle]) {
        assert_eq!(vertices.len(), 2);
        assert!(self.all_points_on_line());
        let v0 = vertices[0];
        let v1 = vertices[1];
        debug_assert!(self.s().get_edge_from_neighbors(v0, v1).is_none());
        let v0_out = self.s().vertex(v0).out_edge().map(|e| e.sym().fix());
        let v1_out = self.s().vertex(v1).out_edge().map(|e| e.sym().fix());
        match (v0_out, v1_out) {
            (None, None) => self.s_mut().connect_two_isolated_vertices(v0, v1, 0),
            (Some(edge), None) => self.s_mut().connect_edge_to_isolated_vertex(edge, v1),
            (None, Some(edge)) => self.s_mut().connect_edge_to_isolated_vertex(edge, v0),
            (Some(v0_edge), Some(v1_edge)) => {
                let v1_edge = self.s().edge(v1_edge).sym().fix();
                self.s_mut().connect_edge_to_edge(v0_edge, v1_edge)
            }
        };
    }

    fn repair_convex_hull(&mut self, vertices: &[FixedVertexHandle]) {
        // A vertex from the convex hull has been removed. This removal can create
        // multiple 'pockets' in the hull that need to be re-triangulated. 
        // 'vertices' contains all vertices that were adjacent to the removed point
        // in ccw order.
        let mut ch: Vec<FixedVertexHandle> = Vec::new();

        // First, determine the new convex hull. Since the points are ordered
        // counterclockwise, this can be done in O(vertices.len()) .
        // This is similar to calculating an upper convex hull.
        for v in vertices {
            let vertex = self.s().vertex(*v);
            let v_pos = (*vertex).position();
            loop {
                if ch.len() >= 2 {
                    let p0 = (*self.s().vertex(ch[ch.len() - 2])).position();
                    let p1 = (*self.s().vertex(ch[ch.len() - 1])).position();
                    let edge = SimpleEdge::new(p0, p1);
                    if Self::Kernel::side_query(&edge, &v_pos).is_on_left_side() {
                        ch.pop();
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
            ch.push(*v);
        }
        // In the next step, analyze if any pockets were created. Follow the
        // edges of the newly created convex hull and check if the edge also exists
        // in the triangulation. If not, a pocket was found that must be filled.
        for vs in ch.windows(2) {
            let v0 = vs[0];
            let v1 = vs[1];
            if self.s().get_edge_from_neighbors(v0, v1).is_none() {
                // The edge does not exists, get all edges of the pocket
                let mut edges = Vec::new();
                let pos = vertices.iter().position(|v| *v == v0).unwrap();
                {
                    let mut cur_edge = self.s().get_edge_from_neighbors(
                        v0, vertices[pos + 1]).unwrap();
                    loop {
                        edges.push(cur_edge.fix());
                        cur_edge = cur_edge.o_next();
                        if cur_edge.from().fix() == v1 {
                            // We have reached the convex hull again, this
                            // closes the pocket.
                            break;
                        }
                    }
                }
                let first = self.s().edge(*edges.first().unwrap()).fix();
                let last = self.s().edge(*edges.last().unwrap()).fix();
                edges.push(self.s_mut().create_face(last, first));
                // Fill the pocket
                self.fill_hole(edges);
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
        for edge in &loop_edges[2 .. loop_edges.len() - 1] {
            let edge = self.s_mut().create_face(last_edge, *edge);
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
            if !Self::Kernel::contained_in_circumference(&v0, &v1, &vl, &vr) {
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

    #[cfg(test)]
    fn sanity_check(&self) {
        self.s().sanity_check();
        for face in self.s().faces().skip(1) {
            let triangle = face.as_triangle();
            assert!(Self::Kernel::is_ordered_ccw(
                &(*triangle[0]).position(),
                &(*triangle[1]).position(),
                &(*triangle[2]).position()));
        }
        if self.all_points_on_line() {
            assert_eq!(self.s().num_faces(), 1);
            assert_eq!(self.s().num_edges() as i32, 0.max(self.s().num_vertices() as i32 - 1));
            for edge in self.s().edges() {
                assert_eq!(edge.face(), self.infinite_face());
            }
        } else {
            for vertex in self.s().vertices() {
                assert!(vertex.out_edge().is_some());
            }
            for edge in self.s().edges() {
                assert_ne!(edge.face(), edge.sym().face());
            }
        }
    }
}
