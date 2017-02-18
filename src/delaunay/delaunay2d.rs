// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module offers a delaunay triangulation and various iterators for inspecting
//! its structure. The main data structure is `DelaunayTriangulation`.

use num::{One, Float, Zero, one, zero};
use traits::{SpatialObject, HasPosition2D, SpadeFloat, HasPosition};
use vector_traits::{VectorN, VectorNExtensions, TwoDimensional, ThreeDimensional};
use kernels::{DelaunayKernel, TrivialKernel, FloatKernel};
use primitives::{SimpleEdge, SimpleTriangle};
use std::marker::PhantomData;
use std::collections::{HashSet};

use delaunay::*;

pub type FloatDelaunayTriangulation<T, L> = DelaunayTriangulation<T, FloatKernel, L>;
pub type IntDelaunayTriangulation<T, L> = DelaunayTriangulation<T, TrivialKernel, L>;

/// Yields information about a point's position in triangulation.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum PositionInTriangulation<V: Copy, F: Copy, E: Copy> {
    /// The point is contained in a triangle. The three given edge handles traverse the triangle in
    /// counterclockwise direction
    InTriangle(F),
    /// The point is outside the convex hull. The two given handles mark an edge on the convex
    /// hull that is close to the point.
    OutsideConvexHull(E),
    /// An object with this position has already been inserted. Its handle is given.
    OnPoint(V),
    /// The point lies on the edge between two vertices. The handles of these vertices are given.
    OnEdge(E),
    /// There is no valid triangulation yet, thus, all points inserted so far lie on a line.
    NoTriangulationPresent,
}

fn calculate_convex_polygon_double_area<V>(
    poly: &Vec<V>) -> V::Scalar 
    where V: TwoDimensional {
    let mut sum = zero();
    for vertices in poly[1 .. ].windows(2) {
        let triangle = SimpleTriangle::new(poly[0].clone(), vertices[0].clone(), vertices[1].clone());
        sum = sum + triangle.double_area();
    }
    sum
}

/// A 2D Delaunay triangulation.
/// 
/// A delaunay triangulation is a special triangulation of a set of points that fulfills some
/// suitable properties for geometric operations like interpolation.
///
/// Objects that are inserted into the triangulation have to implement the `HasPosition2D` trait.
/// The trait is implemented for all types that implement `TwoDimensional`, like `Vector2` from
/// the `cgmath` and `nalgebra` package or `[S; 2]` for `S: SpadeNum`.
/// 
/// A straightforward delaunay triangulation implementation will suffer from precision problems:
/// various geometric queries can fail if imprecise calculations (like native `f32` / `f64` operations)
/// are used. Those failures can yield to incorrect results or panics at runtime.
/// To prevent those crashes, Spade offers a few "calculation kernels" that may fit the 
/// individual needs of an application.
///
/// # Example
///
/// ```
/// extern crate nalgebra;
/// extern crate spade;
///
/// use nalgebra::Vector2;
/// use spade::delaunay::FloatDelaunayTriangulation;
///
/// # fn main() {
///   let mut delaunay = FloatDelaunayTriangulation::with_walk_lookup();
///   delaunay.insert(Vector2::new(0.0, 1.0));
///   delaunay.insert(Vector2::new(0.0, -1.0));
///   delaunay.insert(Vector2::new(1.0, 0.0));
///   delaunay.insert(Vector2::new(-1.0, 0.0));
///   delaunay.insert(Vector2::new(0.0, 0.0));
///   for face in delaunay.triangles() {
///     let triangle = face.as_triangle();
///     println!("found triangle: {:?} -> {:?} -> {:?}", *triangle[0], *triangle[1], *triangle[2]);
///   }
///   for edge in delaunay.edges() {
///     println!("found an edge: {:?} -> {:?}", *edge.from(), *edge.to());
///   }
/// # }
/// ```
///
/// # Iterating
/// A triangulation has three elements - vertices, edges and triangles - that can be iterated over.
/// Use `vertices()` `edges()` and `triangles()` to call appropriate, non-mutating iterators.
///
/// # Mutating
/// Adding vertices is always a legal operation. Removing vertices is not possible.
/// Mutation of a vertex is possible with `lookup_mut(..)`.
///
/// While `VertexHandle`s are intended to be used for short lived 
/// iteration purposes, `FixedVertexHandle`s can be made persistent and have 
/// a `'static` lifetime. A `VertexHandle` can be transformed to its fixed
/// variant by calling `VertexHandle::fix()`.
/// Use `handle` or `handle_mut` to resolve a `FixedVertexHandle`.
/// 
/// A vertex position must not be changed after the vertex has been inserted.
///
/// # Type Parameters
/// `DelaunayTriangulation` has three type parameters: `V`, `K` and `L`.
/// `V: HasPosition2D` defines the delaunay's vertex type. 
/// `K: DelaunayKernel` defines the triangulations calculation kernel.
/// For more information, see `kernels::DelaunayKernel`.
/// `L` Defines the delaunay lookup structure.
/// For more information, see `DelaunayLookupStructure`.
pub struct DelaunayTriangulation<V, K, L = RTreeDelaunayLookup<<V as HasPosition>::Vector>>
    where V: HasPosition2D,
          V::Vector: TwoDimensional,
          L: DelaunayLookupStructure<V::Vector>,
{
    __kernel: PhantomData<K>,
    s: DCEL<V>,
    all_points_on_line: bool,
    lookup: L,
}

impl<V, K, L> Clone for DelaunayTriangulation<V, K, L> 
    where V: HasPosition2D + Clone,
          V::Vector: TwoDimensional,
          L: DelaunayLookupStructure<V::Vector> {

    fn clone(&self) -> DelaunayTriangulation<V, K, L> {
        DelaunayTriangulation {
          __kernel: Default::default(),
            s: self.s.clone(),
            all_points_on_line: self.all_points_on_line,
            lookup: self.lookup.clone(),
        }
    }
}

impl <V, K> DelaunayTriangulation<V, K>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          V::Vector: TwoDimensional,
{

    pub fn with_tree_lookup() -> DelaunayTriangulation<V, K> {
        DelaunayTriangulation::new()
    }

    pub fn with_walk_lookup() -> DelaunayTriangulation<V, K, TriangulationWalkLookup<V::Vector>> {
        DelaunayTriangulation::new()
    }
}

impl <V, K, L> DelaunayTriangulation<V, K, L>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          V::Vector: TwoDimensional,
          L: DelaunayLookupStructure<V::Vector>
{

    /// Creates a new Delaunay triangulation.
    ///
    /// Using this method directly can be a bit cumbersome due to type annotations, consider using
    /// The short hand definitions `FloatDelaunayTriangulation` or `IntDelaunayTriangulation` in
    /// combination with the methods `with_walk_lookup` and `with_tree_lookup`.
    /// If you want to have full control over a triangulations kernel and lookup structure, you
    /// can use this method like this:
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// use spade::delaunay::{DelaunayTriangulation, TriangulationWalkLookup};
    /// use spade::kernels::TrivialKernel;
    /// # fn main() {
    /// let mut triangulation = DelaunayTriangulation::<_, TrivialKernel, TriangulationWalkLookup<_>>::new();
    /// # triangulation.insert(nalgebra::Vector2::new(332f64, 123f64));
    /// # }
    /// ```
    /// 
    /// Usually, the omitted types (the triangulation's vertex type) can be inferred from a call
    /// to `insert`.
    pub fn new() -> DelaunayTriangulation<V, K, L> {
        DelaunayTriangulation {
            __kernel: Default::default(),
            s: DCEL::new(),
            all_points_on_line: true,
            lookup: Default::default(),
        }
    }

    /// Creates a dynamic vertex handle from a fixed handle.
    /// May panic if the handle was not obtained from this triangulation.
    pub fn handle(&self, handle: FixedVertexHandle) -> VertexHandle<V> {
        self.s.vertex(handle)
    }

    /// Returns the number of vertices in this triangulation.
    pub fn num_vertices(&self) -> usize {
        self.s.num_vertices()
    }

    /// Returns the number of faces in this triangulation.
    pub fn num_faces(&self) -> usize {
        self.s.num_faces()
    }

    /// Returns the number of finite faces in this triangulation.
    /// As there is always exactly one face, `self.num_faces() - 1`.
    pub fn num_finite_faces(&self) -> usize {
        self.s.num_faces() - 1
    }

    /// Returns the number of edges in this triangulation.
    pub fn num_edges(&self) -> usize {
        self.s.num_edges()
    }

    /// Returns an iterator over all triangles.
    pub fn triangles(&self) -> FacesIterator<V> {
        let mut result = self.s.faces();
        // Skip the outer face
        result.next();
        result
    }

    /// Returns an iterator over all edges.
    pub fn edges(&self) -> EdgesIterator<V> {
        self.s.edges()
    }

    /// Returns an iterator over all vertices.
    pub fn vertices(&self) -> VerticesIterator<V> {
        self.s.vertices()
    }

    pub fn infinite_face(&self) -> FaceHandle<V> {
        self.s.face(0)
    }

    pub fn is_degenerate(&self) -> bool {
        self.all_points_on_line
    }

    fn initial_insertion(&mut self, t: V) -> Result<FixedVertexHandle, FixedVertexHandle> {
        assert!(self.all_points_on_line);
        // Inserts points if no points are present or if all points
        // lie on the same line
        let new_pos = t.position();
        for vertex in self.s.fixed_vertices() {
            let pos = (*self.s.vertex(vertex)).position();
            if pos == new_pos {
                self.s.update_vertex(vertex, t);
                return Result::Err(vertex);
            }
        }

        if self.s.num_vertices() <= 1 {
            return Result::Ok(self.s.insert_vertex(t));
        }

        // Check if the new point is on the same line as all points in the
        // triangulation
        let from = (*self.s.vertex(0)).position();
        let to = (*self.s.vertex(1)).position();
        let edge = SimpleEdge::new(from.clone(), to.clone());
        if K::side_query(&edge, &new_pos).is_on_line() {
            return Result::Ok(self.s.insert_vertex(t));
        }
        // The point does not lie on the same line as all other points.
        // Start creating a triangulation
        let dir = to.sub(&from);
        let mut vertices: Vec<_> = self.s.vertices().map(
            |v| (v.fix(), dir.dot(&(*v).position()))).collect();
        // Sort vertices according to their position on the line
        vertices.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());

        // Create line
        let is_ccw = K::is_ordered_ccw(&new_pos, &from, &to);
        let mut last_edge = self.s.connect_two_isolated_vertices(vertices[0].0, vertices[1].0, 0);
        let mut edges = vec![last_edge];
        for v in vertices.iter().skip(2) {
            let edge = self.s.connect_edge_to_isolated_vertex(last_edge, v.0);
            edges.push(self.s.edge(edge).fix());
            last_edge = edge;
        }
        if is_ccw {
            edges.reverse();
        }
        let new_vertex = self.s.insert_vertex(t);
        // Connect all points on the line to the new vertex
        let mut last_edge = *edges.first().unwrap();
        if !is_ccw {
            last_edge = self.s.edge(last_edge).sym().fix();
        }
        last_edge = self.s.connect_edge_to_isolated_vertex(last_edge, new_vertex);
        for e in edges {
            let e = if !is_ccw {
                self.s.edge(e).sym().fix()
            } else {
                e
            };
            last_edge = self.s.create_face(last_edge, e);
            last_edge = self.s.edge(last_edge).sym().fix();
        }
        self.all_points_on_line = false;
        Result::Ok(new_vertex)
    }

    fn get_convex_hull_edges_for_point(&self, first_edge: FixedEdgeHandle,
                                       point: &V::Vector) -> Vec<FixedEdgeHandle> 
    {
        let mut result = Vec::new();
        let first_edge = self.s.edge(first_edge);
        debug_assert!(K::side_query(&to_simple_edge(&first_edge), point).is_on_left_side());

        let mut last_edge = first_edge;
        result.push(last_edge.fix());
        // Follow the first edge in cw and ccw direction
        loop {
            last_edge = last_edge.o_next();
            let query = K::side_query(&to_simple_edge(&last_edge), point);
            if query.is_on_left_side() {
                result.push(last_edge.fix());
            } else {
                break;
            }
        }

        last_edge = first_edge;
        loop {
            last_edge = last_edge.o_prev();
            let query = K::side_query(&to_simple_edge(&last_edge), point);
            if query.is_on_left_side() {
                result.insert(0, last_edge.fix());
            } else {
                break;
            }
        }
        result
    }

    fn insert_outside_convex_hull(
        &mut self, closest_edge: FixedEdgeHandle, t: V) -> FixedVertexHandle {
        let position = t.position();
        let ch_edges = self.get_convex_hull_edges_for_point(closest_edge, &position);
        let new_handle = self.s.insert_vertex(t);
        // Make new connections
        let mut last_edge = self.s.connect_edge_to_isolated_vertex(
            *ch_edges.last().unwrap(), new_handle);

        for edge in ch_edges.iter().rev() {
            last_edge = self.s.create_face(last_edge, *edge);
            // Reverse last_edge
            last_edge = self.s.edge(last_edge).sym().fix();
        }
        self.legalize_edges(ch_edges, new_handle);
        new_handle
    }
    
    fn insert_into_triangle(&mut self, face: FixedFaceHandle, t: V) -> FixedVertexHandle {
        let new_handle = self.s.insert_vertex(t);
        let (e0, e1, e2) = {
            let face = self.s.face(face);
            let adj = face.adjacent_edge().unwrap();
            (adj.o_prev().fix(), adj.fix(), adj.o_next().fix())
        };

        let mut last_edge = self.s.connect_edge_to_isolated_vertex(e2, new_handle);
        last_edge = self.s.edge(last_edge).sym().fix();
        last_edge = self.s.create_face(e0, last_edge);
        last_edge = self.s.edge(last_edge).sym().fix();
        self.s.create_face(e1, last_edge);
        self.legalize_edges(vec![e0, e1, e2], new_handle);
        new_handle
    }

    fn is_ch_edge(&self, edge: FixedEdgeHandle) -> bool {
        let edge = self.s.edge(edge);
        let sym = edge.sym();
        let inf = self.infinite_face();
        edge.face() == inf || sym.face() == inf
    }

    fn get_left_triangle(&self, edge: (FixedVertexHandle, FixedVertexHandle)) 
                         -> Option<FixedVertexHandle> {
        let edge_handle = EdgeHandle::from_neighbors(&self.s, edge.0, edge.1).unwrap();
        let ccw_handle = edge_handle.ccw().to();
        let query = K::side_query(&to_simple_edge(&edge_handle), &(*ccw_handle).position());
        if query.is_on_left_side() {
            debug_assert!(EdgeHandle::from_neighbors(&self.s, ccw_handle.fix(), edge.1).is_some());
            Some(ccw_handle.fix())
        } else {
            None
        }
    }

    fn get_right_triangle(&self, edge: (FixedVertexHandle, FixedVertexHandle))
        -> Option<FixedVertexHandle> {
        self.get_left_triangle((edge.1, edge.0))
    }

    fn insert_on_edge(&mut self, edge: FixedEdgeHandle, t: V) -> FixedVertexHandle {
        let new_handle = self.s.insert_vertex(t);
        let mut illegal_edges = Vec::new();
        let (from, to) = {
            let edge = self.s.edge(edge);
            (edge.from().fix(), edge.to().fix())
        };

        let left_handle_opt = self.get_left_triangle((from, to));
        let right_handle_opt = self.get_right_triangle((from, to));
        let edge_handle = EdgeHandle::from_neighbors(&self.s, from, to).unwrap().fix();
        self.s.split_edge(edge_handle, new_handle);
        if let Some(left_handle) = left_handle_opt {
            let edge1 = EdgeHandle::from_neighbors(&self.s, to, left_handle).unwrap().fix();
            let edge0 = EdgeHandle::from_neighbors(&self.s, left_handle, from).unwrap().fix();
            let edge_mid = EdgeHandle::from_neighbors(&self.s, from, new_handle).unwrap().fix();

            self.s.create_face(edge_mid, edge0);
            illegal_edges.push(edge0);
            illegal_edges.push(edge1);
        }
        if let Some(right_handle) = right_handle_opt {
            let edge0 = EdgeHandle::from_neighbors(&self.s, from, right_handle).unwrap().fix();
            let edge1 = EdgeHandle::from_neighbors(&self.s, right_handle, to).unwrap().fix();
            let edge_mid = EdgeHandle::from_neighbors(&self.s, to, new_handle).unwrap().fix();
            self.s.create_face(edge_mid, edge1);
            illegal_edges.push(edge0);
            illegal_edges.push(edge1);

        }
        self.legalize_edges(illegal_edges, new_handle);
        new_handle
    }

    /// Returns information about the location of a point in a triangulation.
    pub fn locate(
        &self, point: &V::Vector) -> PositionInTriangulation<VertexHandle<V>, FaceHandle<V>, EdgeHandle<V>> {
        use self::PositionInTriangulation::*;

        match self.locate_with_hint_option_fixed(point, None) {
            NoTriangulationPresent => NoTriangulationPresent,
            InTriangle(face) => InTriangle(self.s.face(face)),
            OutsideConvexHull(edge) => OutsideConvexHull(self.s.edge(edge)),
            OnPoint(vertex) => OnPoint(self.s.vertex(vertex)),
            OnEdge(edge) => OnEdge(self.s.edge(edge)),
        }
    }

    fn locate_with_hint_option_fixed(&self, point: &V::Vector, hint: Option<FixedVertexHandle>) -> 
        PositionInTriangulation<FixedVertexHandle, FixedFaceHandle, FixedEdgeHandle> {
        if self.all_points_on_line {
            // TODO: We might want to check if the point is on the line or already contained
            PositionInTriangulation::NoTriangulationPresent
        } else {
            let start = hint.unwrap_or_else(|| self.lookup.find_close_handle(point).unwrap_or(0));
            self.locate_with_hint_fixed(point, start) 
        }
    }

    fn locate_with_hint_fixed(&self, point: &V::Vector, start: FixedVertexHandle) -> 
        PositionInTriangulation<FixedVertexHandle, FixedFaceHandle, FixedEdgeHandle> 
    {
        let mut cur_edge = self.s.vertex(start).out_edge().expect("Cannot start search with an isolated vertex");
        let mut cur_query = K::side_query(&to_simple_edge(&cur_edge), point);
        // Invariant: point must not be on the right side of cur_edge
        if cur_query.is_on_right_side() {
            cur_edge = cur_edge.sym();
            cur_query = cur_query.sym();
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
       
            let next_query = K::side_query(&to_simple_edge(&next), point);
            if next_query.is_on_right_side_or_on_line() {
                // We continue walking into the face right of next
                cur_edge = next.sym();
                cur_query = next_query.sym();
            } else {
                // Check if cur_edge.o_prev is also on the left side
                let prev = cur_edge.o_prev();
                let prev_query = K::side_query(&to_simple_edge(&prev), point);
                if prev_query.is_on_right_side_or_on_line() {
                    // We continue walking into the face right of prev
                    cur_edge = prev.sym();
                    cur_query = prev_query.sym();
                } else {
                    if cur_query.is_on_line() {
                        return PositionInTriangulation::OnEdge(cur_edge.fix());
                    }
                    return PositionInTriangulation::InTriangle(cur_edge.face().fix());
                }
            }
        }
    }

    /// Inserts a new vertex into the triangulation
    /// This operations runs in O(log n) on average, n denotes the number of vertices contained
    /// in the triangulation. If the point has already been contained in the
    /// triangulation, the old vertex is overwritten.
    ///
    /// Returns a handle to the new vertex. Use this handle with
    /// `DelaunayTriangulation::handle(..)` to refer to the vertex.
    pub fn insert(&mut self, t: V) -> FixedVertexHandle {
        self.insert_with_hint_option(t, None)
    }

    pub fn insert_with_hint(&mut self, t: V, hint: FixedVertexHandle) -> FixedVertexHandle {
        self.insert_with_hint_option(t, Some(hint))
    }


    fn insert_with_hint_option(&mut self, t: V, hint: Option<FixedVertexHandle>) -> FixedVertexHandle {
        let pos = t.position();
        let position_in_triangulation = self.locate_with_hint_option_fixed(&pos, hint);
        let insertion_result = match position_in_triangulation {
            PositionInTriangulation::OutsideConvexHull(edge) => {
                Result::Ok(self.insert_outside_convex_hull(edge, t))
            },
            PositionInTriangulation::InTriangle(face) => {
                Result::Ok(self.insert_into_triangle(face, t))
            },
            PositionInTriangulation::OnEdge(edge) => {
                Result::Ok(self.insert_on_edge(edge, t))
            },
            PositionInTriangulation::OnPoint(vertex) => {
                self.s.update_vertex(vertex, t);
                Result::Err(vertex)
            },
            PositionInTriangulation::NoTriangulationPresent => {
                self.initial_insertion(t)
            }
        };
        match insertion_result {
            Result::Ok(new_handle) => {
                self.lookup.insert_vertex_entry(VertexEntry {
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

    fn legalize_edges(&mut self, mut edges: Vec<FixedEdgeHandle>, new_vertex: FixedVertexHandle) {
        let position = (*self.s.vertex(new_vertex)).position();
        while let Some(e) = edges.pop() {
            if !self.is_ch_edge(e) {
                let (v0, v1, v2, e1, e2);
                {
                    let edge = self.s.edge(e);
                    v0 = (*edge.from()).position();
                    v1 = (*edge.to()).position();
                    v2 = (*edge.cw().to()).position();
                    e1 = edge.sym().o_next().fix();
                    e2 = edge.sym().o_prev().fix();
                }
                debug_assert!(K::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(K::is_ordered_ccw(&position, &v0, &v1));
                if K::contained_in_circumference(&v1, &v2, &v0, &position) {
                    // The edge is illegal
                    self.s.flip_cw(e);
                    edges.push(e1);
                    edges.push(e2);
                }
            }
        }
    }

    /// Attempts to remove a vertex from the triangulation.
    pub fn locate_and_remove(&mut self, point: &V::Vector) -> Option<V> {
        use self::PositionInTriangulation::*;
        match self.locate_with_hint_option_fixed(point, None) {
            OnPoint(handle) => Some(self.remove(handle)),
            _ => None,
        }
    }


    /// Removes a vertex from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    /// *Note*: This operation will invalidate a single vertex handle.
    /// Do not rely on any fixed vertex handle created before the removal
    /// to be valid.
    pub fn remove(&mut self, vertex: FixedVertexHandle) -> V {
        let mut neighbors = Vec::new();
        let mut ch_removal = false;
        for edge in self.handle(vertex).ccw_out_edges() {
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
        
        let infinite = self.infinite_face().fix();
        let VertexRemovalResult { updated_vertex, data } = 
            self.s.remove_vertex(vertex, Some(infinite));
        if let Some(updated_vertex) = updated_vertex {
            // Update rtree if necessary
            let pos = (*self.handle(vertex)).position();
            self.lookup.update_vertex_entry(VertexEntry {
                point: pos,
                handle: vertex,
            });
            for mut n in &mut neighbors {
                if *n == updated_vertex {
                    *n = vertex;
                    break;
                }
            }
        }
        
        // Remove deleted point from internal rtree
        let vertex_pos = data.position();
        self.lookup.remove_vertex_entry(&vertex_pos);

        if !self.all_points_on_line {
            if ch_removal {
                // We removed a vertex from the convex hull
                self.repair_convex_hull(&neighbors);
                if self.s.num_faces() == 1 {
                    self.make_degenerate(); 
                }
            } else {
                let loop_edges: Vec<_> = {
                    let first = EdgeHandle::from_neighbors(
                        &self.s, neighbors[0], neighbors[1]).unwrap();
                    first.o_next_iterator().map(|e| e.fix()).collect()
                };
                self.fill_hole(loop_edges);
            }
    }
        data
    }

    fn fill_hole(&mut self, loop_edges: Vec<FixedEdgeHandle>) {
        let mut border_edges = HashSet::new();
        
        for e in &loop_edges {
            border_edges.insert(*e);
            border_edges.insert(self.s.edge(*e).sym().fix());
        }

        let last_edge = *loop_edges.last().unwrap();

        // Fill the hole
        let mut todo = Vec::new();
        for i in 2 .. loop_edges.len() - 1 {
            let edge = self.s.create_face(last_edge, loop_edges[i]);
            todo.push(edge);
        }
        // Legalize edges
        while let Some(fixed_edge_handle) = todo.pop() {
            let (v0, v1, vl, vr, e1, e2, e3, e4);
            {
                let edge = self.s.edge(fixed_edge_handle);
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
                self.s.flip_cw(fixed_edge_handle);
                
                for e in &[e1, e2, e3, e4] {
                    if !border_edges.contains(e) {
                        todo.push(*e);
                    }
                }
            }
        }
    }

    fn repair_convex_hull(&mut self, vertices: &Vec<FixedVertexHandle>) {
        let mut ch: Vec<FixedVertexHandle> = Vec::new();
        for v in vertices {
            let vertex = self.s.vertex(*v);
            let v_pos = (*vertex).position();
            loop {
                if ch.len() >= 2 {
                    let p0 = (*self.s.vertex(ch[ch.len() - 2])).position();
                    let p1 = (*self.s.vertex(ch[ch.len() - 1])).position();
                    let edge = SimpleEdge::new(p0, p1);
                    if K::side_query(&edge, &v_pos).is_on_left_side() {
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
        for vs in ch.windows(2) {
            let v0 = vs[0];
            let v1 = vs[1];
            if EdgeHandle::from_neighbors(&self.s, v0, v1).is_none() {
                let mut edges = Vec::new();
                let pos = vertices.iter().position(|v| *v == v0).unwrap();
                {
                    let mut cur_edge = EdgeHandle::from_neighbors(
                        &self.s, v0, vertices[pos + 1]).unwrap();
                    loop {
                        edges.push(cur_edge.fix());
                        cur_edge = cur_edge.o_next();
                        if cur_edge.from().fix() == v1 {
                            break;
                        }
                    }
                }
                let first = self.s.edge(*edges.first().unwrap()).fix();
                let last = self.s.edge(*edges.last().unwrap()).fix();
                edges.push(self.s.create_face(last, first));
                self.fill_hole(edges);
            }
        }
    }

    fn make_degenerate(&mut self) {
        // Assume all points lie on a line.
        self.s.clear_edges_and_faces();
        self.all_points_on_line = true;
    }

    #[cfg(test)]
    fn sanity_check(&self) {
        self.s.sanity_check();
        for face in self.triangles() {
            let triangle = face.as_triangle();
            assert!(K::is_ordered_ccw(
                &(*triangle[0]).position(),
                &(*triangle[1]).position(),
                &(*triangle[2]).position()));
        }
        if self.is_degenerate() {
            assert_eq!(self.num_finite_faces(), 0);
            assert_eq!(self.num_edges(), 0);
        } else {
            for vertex in self.vertices() {
                assert!(vertex.out_edge().is_some());
            }
        }
        for edge in self.edges() {
            assert!(edge.face() != edge.sym().face());
        }
    }
}

impl <V, K> DelaunayTriangulation<V, K, RTreeDelaunayLookup<V::Vector>>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          V::Vector: TwoDimensional {

    /// Checks if the triangulation contains an object with a given coordinate.
    pub fn lookup(&self, point: &V::Vector) -> Option<VertexHandle<V>> {
        let handle = self.lookup.lookup(point);
        handle.map(|h| self.s.vertex(h.handle))
    }

    /// Removes a vertex at a given position from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    /// *Note*: This operation will invalidate a single vertex handle.
    /// Do not rely on any fixed vertex handle created before the removal
    /// to be valid.
    pub fn lookup_and_remove(&mut self, point: &V::Vector) -> Option<V> {
        let handle = self.lookup(point).map(|h| h.fix());
        handle.map(|h| self.remove(h))
    }

    // /// Returns all vertices contained in a rectangle.
    // pub fn lookup_in_rect(&self, rect: &BoundingRect<V::Vector>) -> Vec<VertexHandle<V>> {
    //     let fixed_handles = self.points.lookup_in_rectangle(rect);
    //     fixed_handles.iter().map(|entry| self.s.vertex(entry.handle)).collect()
    // }

    // /// Returns all vertices contained in a circle.
    // pub fn lookup_in_circle(&self, center: &V::Vector, 
    //                         radius2: &<V::Vector as VectorN>::Scalar) -> Vec<VertexHandle<V>> {
    //     let fixed_handles = self.points.lookup_in_circle(center, radius2);
    //     fixed_handles.iter().map(|entry| self.s.vertex(entry.handle)).collect()
    // }
}

impl <V, K, L> DelaunayTriangulation<V, K, L>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          L: DelaunayLookupStructure<V::Vector>,
          V::Vector: TwoDimensional {

    pub fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V> {
        self.s.vertex(handle)
    }    

    /// Returns a mutable reference to the vertex data referenced by a 
    /// `FixedVertexHandle`.
    pub fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V {
        self.s.vertex_mut(handle)
    }

    // /// Checks if the triangulation contains an object and returns a mutable
    // /// reference to it. Note that this method will return a mutable reference, while
    // /// `lookup(..)` returns a vertex handle.
    // pub fn lookup_mut(&mut self, point: &V::Vector) -> Option<&mut V> {
    //     let handle = self.points.lookup(point).map(|e| e.handle);
    //     if let Some(handle) = handle {
    //         Some(self.handle_mut(handle))
    //     } else {
    //         None
    //     }
    // }
}

impl <V, K, L> DelaunayTriangulation<V, K, L> 
    where V: HasPosition2D, <V::Vector as VectorN>::Scalar: SpadeFloat,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar> ,
          L: DelaunayLookupStructure<V::Vector>,
          V::Vector: TwoDimensional
{
    pub fn barycentric_interpolation<F> (&self, point: &V::Vector, f: F) 
                                         -> Option<<V::Vector as VectorN>::Scalar> 
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar {
        let vertices = match self.locate(point) {
            PositionInTriangulation::NoTriangulationPresent => return None,
            PositionInTriangulation::OnPoint(v) => vec![v],
            PositionInTriangulation::OnEdge(e) => vec![e.from(), e.to()],
            PositionInTriangulation::InTriangle(f) => {
                let vs = f.as_triangle();
                vec![vs[0], vs[1], vs[2]]
            }
            PositionInTriangulation::OutsideConvexHull(e) => vec![e.from(), e.to()],
        };
        if vertices.len() == 1 {
            Some(f(&*vertices[0]))
        } else if vertices.len() == 2 {
            let p0 = vertices[0].position();
            let p1 = vertices[1].position();
            let one = <<V::Vector as VectorN>::Scalar>::one();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(zero(), edge.project_point(point), one);
            let w0 = one - w1;
            Some(w1 * f(&*vertices[1]) + w0 * f(&*vertices[0]))
        } else {
            let triangle = ::primitives::SimpleTriangle::new(
                vertices[0].position(),
                vertices[1].position(),
                vertices[2].position());
            let b_coords = triangle.barycentric_interpolation(point);
            let w0 = f(&*vertices[0]);
            let w1 = f(&*vertices[1]);
            let w2 = f(&*vertices[2]);
            Some(w0 * b_coords[0] + w1 * b_coords[1] + w2 * b_coords[2])
        }
    }
    

    /// Performs a natural neighbor interpolation for a given position.
    ///
    /// Returns `None` if the triangulation has no triangles yet.
    /// Points outside of the convex hull will be interpolated as well.
    /// This operation runs in O(log n) for n inserted vertices.
    ///
    /// # Example
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    ///
    /// # use nalgebra::{Vector2};
    /// use spade::HasPosition;
    /// use spade::delaunay::{FloatDelaunayTriangulation};
    ///
    /// struct PointWithHeight {
    ///   point: Vector2<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f64>;
    ///     fn position(&self) -> Vector2<f64> {
    ///       self.point
    ///     }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_lookup();
    ///   delaunay.insert(PointWithHeight { point: Vector2::new(0.0, 0.0), height: 5. });
    ///   delaunay.insert(PointWithHeight { point: Vector2::new(1.0, 0.0), height: 0. });
    ///   delaunay.insert(PointWithHeight { point: Vector2::new(0.0, 1.0), height: 0. });
    ///   let lookup = Vector2::new(0.2, 0.2);
    ///   // Interpolate the points height
    ///   let interpolated = delaunay.nn_interpolation(&lookup, |p| p.height).unwrap();
    ///   // and insert it afterwards.
    ///   delaunay.insert(PointWithHeight { point: lookup, height: interpolated });
    ///   // Data points themselves will always yield their own height
    ///   assert_eq!(delaunay.nn_interpolation(&Vector2::new(0.0, 0.0), |p| p.height),
    ///              Some(5.0));
    /// }
    /// ```
    pub fn nn_interpolation<R, F: Fn(&V) -> R>(
        &self, point: &V::Vector, f: F) -> Option<R>
        where R: ::std::ops::Mul<<V::Vector as VectorN>::Scalar, Output=R> + ::std::ops::Add<Output=R> 
    + ::std::ops::Sub<R> {
        
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);

        let mut sum = None;
        for (index, fixed_handle) in nns.iter().enumerate() {
            let new = f(&*self.s.vertex(*fixed_handle)) * ws[index];
            sum = match sum {
                Some(val) => Some(val + new),
                None => Some(new),
            }
        }
        sum
    }

    fn get_weights(&self, nns: &Vec<FixedVertexHandle>, 
                   point: &V::Vector) -> Vec<<V::Vector as VectorN>::Scalar> {

        if nns.len() == 1 {
            return vec![one()];
        }

        if nns.len() == 2 {
            let p0 = (*self.s.vertex(nns[0])).position();
            let p1 = (*self.s.vertex(nns[1])).position();
            let one = <<V::Vector as VectorN>::Scalar>::one();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(zero(), edge.project_point(point), one);
            let w0 = one - w1;
            return vec![w0, w1];
        }
        let mut result = Vec::new();
        let len = nns.len();

        let mut point_cell = Vec::new();
        for (index, cur) in nns.iter().enumerate() {
            let cur_pos = (*self.s.vertex(*cur)).position();
            let next = (*self.s.vertex(nns[(index + 1) % len])).position();
            let triangle = SimpleTriangle::new(next, cur_pos, point.clone());
            point_cell.push(triangle.circumcenter());
        }

        let mut areas = Vec::new();

        for (index, cur) in nns.iter().enumerate() {
            let cur_pos = (*self.s.vertex(*cur)).position();
            let prev = nns[((index + len) - 1) % len];
            let next = nns[(index + 1) % len];
            let mut ccw_edge = EdgeHandle::from_neighbors(&self.s, *cur, prev).unwrap();
            let mut polygon = Vec::new();
            polygon.push(point_cell[((index + len) - 1) % len].clone());
            loop {
                if ccw_edge.to().fix() == next {
                    break;
                }
                let cw_edge = ccw_edge.cw();
                let triangle = SimpleTriangle::new((*ccw_edge.to()).position().clone(), 
                                                   (*cw_edge.to()).position().clone(),
                                                   cur_pos.clone());
                polygon.push(triangle.circumcenter());
                ccw_edge = cw_edge;
            }
            polygon.push(point_cell[index].clone());
            areas.push(calculate_convex_polygon_double_area(&polygon));
        }
        let mut sum = zero();
        for area in &areas {
            sum += *area;
        }
        for i in 0 .. len {
            // Normalize weights
            result.push(areas[i] / sum);
        }
        result
    }

    fn get_natural_neighbors(&self, position: &V::Vector) -> Vec<FixedVertexHandle> {
        match self.locate_with_hint_option_fixed(position, None) {
            PositionInTriangulation::InTriangle(face) => {
                let vs = self.s.face(face).as_triangle();
                self.inspect_flips(vec![(vs[0].fix(), vs[2].fix()), (vs[2].fix(), vs[1].fix()), 
                                        (vs[1].fix(), vs[0].fix())], position)
            },
            PositionInTriangulation::OnEdge(edge) => {
                let (h0, h1) = {
                    let edge = self.s.edge(edge);
                    (edge.from().fix(), edge.to().fix())
                };
                let left_opt = self.get_left_triangle((h0, h1));
                let right_opt = self.get_right_triangle((h0, h1));
                if let (Some(left_handle), Some(right_handle)) = (left_opt, right_opt) {
                    let mut direct_neighbors = Vec::new();
                    direct_neighbors.push((h0, left_handle));
                    direct_neighbors.push((left_handle, h1));
                    direct_neighbors.push((h1, right_handle));
                    direct_neighbors.push((right_handle, h0));
                    self.inspect_flips(direct_neighbors, position)
                } else {
                    vec![h0, h1]
                }
            },
            PositionInTriangulation::OutsideConvexHull(edge) => {
                // Get the closest edge on the convex hull
                let edges = self.get_convex_hull_edges_for_point(
                    edge, position);
                let mut min_dist = <V::Vector as VectorN>::Scalar::max_value();
                let mut min_edge = 0;
                for cur_edge in edges {
                    let edge = self.s.edge(cur_edge);
                    let simple = to_simple_edge(&edge);
                    let new_dist = simple.distance2(position);
                    if new_dist < min_dist {
                        min_edge = cur_edge;
                        min_dist = new_dist;
                    }
                }
                let min_edge = self.s.edge(min_edge);
                let from_pos = (*min_edge.from()).position();
                let to_pos = (*min_edge.to()).position();
                let simple = to_simple_edge(&min_edge);
                if simple.is_projection_on_edge(position) {
                    vec![min_edge.from().fix(), min_edge.to().fix()]
                } else {
                    // Return the closer point
                    if from_pos.distance2(position) < to_pos.distance2(position) {
                        vec![min_edge.from().fix()]
                    } else {
                        vec![min_edge.to().fix()]
                    }
                }
            },
            PositionInTriangulation::OnPoint(fixed_handle) => {
                vec![fixed_handle]
            },
            _ => Vec::new()
        }
    }

    fn inspect_flips(&self, mut edges: Vec<(FixedVertexHandle, FixedVertexHandle)>, position: &V::Vector) -> Vec<FixedVertexHandle> {
        // TODO: edges should store edge references instead of tuples to accomodate
        // the underlying dcel structure better.
        let mut result = Vec::new();
        while let Some((h0, h1)) = edges.pop() {
            if let Some(h2) = self.get_left_triangle((h0, h1)) {
                let v0 = (*self.s.vertex(h0)).position();
                let v1 = (*self.s.vertex(h1)).position();
                let v2 = (*self.s.vertex(h2)).position();
                debug_assert!(!K::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(K::is_ordered_ccw(position, &v1, &v0));
                if K::contained_in_circumference(&v2, &v1, &v0, position) {
                    // The would be illegal
                    // Add edges in ccw order
                    edges.push((h0, h2));
                    edges.push((h2, h1));
                } else {
                    // Edge is legal
                    result.push(h0);
                }
            } else {
                // Edge is part of the convex hull
                result.push(h0);
            }
        }
        result
    }

    /// Estimates a normal value for a given vertex.
    ///
    /// This assumes that the triangulation models some kind of height field, given by the
    /// function `f`.
    /// The normal is the weighted and normalized average of the normals of all triangles
    /// adjacent to the given vertex.
    pub fn estimate_normal<F, RV>(&self, v: FixedVertexHandle, f: &F) -> RV
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar,
              RV: ThreeDimensional<Scalar=<V::Vector as VectorN>::Scalar> {
        let mut v_pos = RV::new();
        let neighbor_positions: Vec<_> = {
            let handle = self.handle(v);
            let v_2d = (*handle).position();
            *v_pos.nth_mut(0) = *v_2d.nth(0);
            *v_pos.nth_mut(1) = *v_2d.nth(1);
            *v_pos.nth_mut(2)
                = f(&*handle);

            handle.ccw_out_edges().map(
                |e| {
                    let pos = (*e.to()).position();
                    let mut result = RV::new();
                    *result.nth_mut(0) = *pos.nth(0);
                    *result.nth_mut(1) = *pos.nth(1);
                    *result.nth_mut(2) = f((&*e.to()));
                    result
                }).collect()
        };
        let mut final_normal = RV::new();
        for index in 0 .. neighbor_positions.len() {
            let p0 = neighbor_positions[index].clone();
            let p1 = neighbor_positions[
                (index + 1) % neighbor_positions.len()].clone();
            let d0 = v_pos.sub(&p0);
            let d1 = v_pos.sub(&p1);
            let normal = d0.cross(&d1);
            if *normal.nth(2) > zero() {
                final_normal = final_normal.add(&normal);
            }

        }
        // Normalize
        final_normal.div(final_normal.length2().sqrt())
    }

    /// Estimates and returns the gradient for a single vertex in this triangulation.
    pub fn estimate_gradient<F>(&self, v: FixedVertexHandle, f: &F) -> V::Vector 
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar {
        use cgmath::Vector3;
        let normal = self.estimate_normal::<_, Vector3<_>>(v, f);
        // Calculate gradient from normal
        let mut gradient = V::Vector::new();
        *gradient.nth_mut(0) = normal.x;
        *gradient.nth_mut(1) = normal.y;
        let g2 = gradient.length2();
        if g2 != zero() {
            let one = <V::Vector as VectorN>::Scalar::one();
            let d2 = one - normal.z * normal.z;
            let a2 = d2 / (one - d2);
            gradient = gradient.mul((a2 / g2).sqrt());
        }
      
        gradient
    }

    /// Interpolates a data point on this triangulation according to Sibson's c1 interpolant.
    ///
    /// The interpolation given by `nn_interpolation` is not differentiable at the triangulation's
    /// data points. Sibson introduced another interpolation scheme that takes the gradient of each
    /// data point into account and offers an interpolation that is differentiable (c1) at the data
    /// points.
    /// The interpolation needs to know the gradients of the points natural neighbors, though.
    /// Spade can estimate them automatically, see `estimate_gradient` and `estimate_gradients`.
    /// The value that should be interpolated is given by `f`, the gradient of a vertex must
    /// be given by `g`.
    ///
    /// # flatness
    /// An additional flatness factor determines how flat the triangulation will be around
    /// the datapoints. A flatness factor of 0.5 is the factor used in sibson's original interpolant.
    /// A flatness of 0.0 or lower is nearly identical to sibson's original interpolant
    /// (`nn_interpolation(..)`). A factor of (exactly) 1.0 should yield best performance since
    /// an exponentiation can be omitted.
    /// # Example
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// 
    /// # use nalgebra::{Vector2, Vector3};
    /// use spade::delaunay::FloatDelaunayTriangulation;
    /// use spade::HasPosition;
    ///
    /// struct PointWithHeight {
    ///   point: Vector2<f64>,
    ///   gradient: Vector2<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f64>;
    ///     fn position(&self) -> Vector2<f64> {
    ///       self.point
    ///     }
    /// }
    ///
    /// impl PointWithHeight {
    ///   fn new(point: Vector2<f64>, height: f64) -> PointWithHeight {
    ///     // Initialize the gradient to any value since it will be overwritten
    ///     PointWithHeight { point: point, height: height, gradient: Vector2::new(0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_lookup();
    ///   // Insert some points here... (skipped)
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(0.0, 0.0), 5.));
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(1.0, 0.0), 0.));
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(0.0, 1.0), 2.));
    ///   // Estimate all gradients and store them:
    ///   delaunay.estimate_gradients(&(|v: &PointWithHeight| v.height),
    ///                               &(|v: &mut PointWithHeight, g| v.gradient = g));
    ///   
    ///   // Now we can use the gradients for interpolation, flatness is set to 2.0:
    ///   let interpolated = delaunay.nn_interpolation_c1_sibson(
    ///      &Vector2::new(0.5, 0.2), 2.0, |v| v.height, |_, v| v.gradient);
    ///   println!("interpolated: {}", interpolated.unwrap());
    /// }
    /// ```
    pub fn nn_interpolation_c1_sibson<F, G> (&self, point: &V::Vector, 
                                             flatness: <V::Vector as VectorN>::Scalar,
                                             f: F, g: G) 
                                             -> Option<<V::Vector as VectorN>::Scalar> 
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar,
              G: Fn(&Self, &VertexHandle<V>) -> V::Vector {
        
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        if ws.len() == 1 {
            return Some(f(&*self.handle(*nns.first().unwrap())));
        }
        let mut sum_c0 = zero();
        let mut sum_c1 = zero();
        let mut sum_c1_weights = zero();
        let mut alpha = <V::Vector as VectorN>::Scalar::zero();
        let mut beta = <V::Vector as VectorN>::Scalar::zero();
        for (index, fixed_handle) in nns.iter().enumerate() {
            let handle = self.s.vertex(*fixed_handle);
            let pos_i = (*handle).position();
            let h_i = f(&*handle);
            let diff = pos_i.sub(point);
            let r_i2 = diff.length2();
            let r_i = r_i2.powf(flatness);
            let c1_weight_i = ws[index] / r_i;
            let grad_i = g(&self, &handle);
            let zeta_i = h_i + diff.dot(&grad_i);
            alpha += c1_weight_i * r_i;
            beta += c1_weight_i * r_i2;
            sum_c1_weights += c1_weight_i;
            sum_c1 += zeta_i * c1_weight_i;
            sum_c0 += h_i * ws[index];
        }
        alpha /= sum_c1_weights;
        sum_c1 /= sum_c1_weights;
        let result = (alpha * sum_c0 + beta * sum_c1) / (alpha + beta);
        Some(result)
    }


    /// Interpolates a data point on this triangulation using Farin's c1 interpolant.
    ///
    /// This method is used in the same way as `nn_interpolation_c1_sibson`. The resulting
    /// interpolant is very similar to sibson's c1 interpolant but uses a different algorithm.
    pub fn nn_interpolation_c1_farin<F, G>(&self, point: &V::Vector, f: F, g: G) 
                                           -> Option<<V::Vector as VectorN>::Scalar>
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar,
              G: Fn(&Self, &VertexHandle<V>) -> V::Vector  {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        let handles: Vec<_> = nns.iter().map(|v| self.handle(*v)).collect();

        let one = <V::Vector as VectorN>::Scalar::one();
        let two = one + one;
        let three = one + one + one;
        let four = two + two;
        let six = three * two;

        let get_edge_height = |i2: usize, i1: usize| {
            // Calculates an edge control point.
            // The edge is given by handles[i2] -> handles[i1]
            let p1 = (*handles[i1]).position();
            let p2 = (*handles[i2]).position();
            let diff = p1 .sub(&p2);
            let grad = g(&self, &handles[i2]);
            f(&*handles[i2]) - grad.dot(&diff) / three
        };

        let mut result = <V::Vector as VectorN>::Scalar::zero();

        for first_index in 0 .. nns.len() {
            for second_index in first_index .. nns.len() {
                for third_index in second_index .. nns.len() {
                    // Determine control point "height"
                    let control_height;
                    let norm;
                    if first_index == second_index && second_index == third_index {
                        // Control point is one of the original data points
                        control_height = f(&*handles[first_index]);
                        norm = six;
                    } else {
                        if first_index == second_index || first_index == third_index ||
                            second_index == third_index 
                        {
                            // Point lies on an edge of the bezier simplex
                            let (i2, i1) = if first_index == second_index {
                                (first_index, third_index)
                            } else if first_index == third_index {
                                (first_index, second_index)
                            } else {
                                (second_index, first_index)
                            };
                            control_height = get_edge_height(i2, i1);
                            norm = two;
                        } else {
                            // We have an inner control point, first != second != third
                            // Get all 6 edge control points of the triangle spanned by
                            // first_index, second_index and third_index
                            let cur_handles = [first_index, second_index, third_index];
                            const EDGE_POINTS: [(usize, usize); 6] = [(0, 1), (0, 2), (1, 0), 
                                                                      (1, 2), (2, 0), (2, 1)];
                            let mut edge_contrib = <V::Vector as VectorN>::Scalar::zero();
                            for &(ref i1, ref i2) in &EDGE_POINTS {
                                let i2 = cur_handles[*i2];
                                let i1 = cur_handles[*i1];
                                // TODO: We could cache the edge points in some way instead of
                                // calculating them anew for each inner point...
                                edge_contrib += get_edge_height(i2, i1);
                            }
                            edge_contrib /= four;
                            
                            let inner_contrib = (f(&*handles[first_index]) +
                                                 f(&*handles[second_index]) + 
                                                 f(&*handles[third_index])) / six;
                            control_height = edge_contrib - inner_contrib;
                            norm = one;
                        }
                    }
                    // Add control height to result, weight it with the appropriate natural
                    // neighbor coordinates.
                    result += six * ws[first_index] * ws[second_index] * ws[third_index] 
                        * control_height / norm;
                }
            }
        }
        Some(result)
    }
}

impl <V, K, L> DelaunayTriangulation<V, K, L> 
    where V: HasPosition2D, <V::Vector as VectorN>::Scalar: SpadeFloat,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          L: DelaunayLookupStructure<V::Vector>,
          V::Vector: TwoDimensional,
{

    /// Estimates a normal for each vertex in the triangulation.
    ///
    /// `f` must yield a "height" value for each vertex,
    /// `g` is a callback function that can be used to store the calculated normals.
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// 
    /// # use nalgebra::{Vector2, Vector3};
    /// use spade::delaunay::FloatDelaunayTriangulation;
    /// use spade::HasPosition;
    ///
    /// struct PointWithHeight {
    ///   point: Vector2<f64>,
    ///   normal: Vector3<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f64>;
    ///     fn position(&self) -> Vector2<f64> {
    ///       self.point
    ///     }
    /// }
    /// impl PointWithHeight {
    ///   fn new(point: Vector2<f64>, height: f64) -> PointWithHeight {
    ///     PointWithHeight { point: point, height: height, normal: Vector3::new(0., 0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_lookup();
    ///   // Insert some points here... (skipped)
    ///   # delaunay.insert(PointWithHeight { point: Vector2::new(0.0, 0.0), height: 5., normal: Vector3::new(0., 0., 0.)});
    ///   // Then, estimate all normals at once:
    ///   delaunay.estimate_normals(&(|v: &PointWithHeight| v.height), &(|v: &mut PointWithHeight, n| v.normal = n));
    ///   
    ///   // And print them
    ///   for vertex in delaunay.vertices() {
    ///      println!("vertex: {:?}, normal: {:?}", vertex.position(), vertex.normal);
    ///   }
    /// }
    /// ```
    pub fn estimate_normals<F, G, RV>(&mut self, f: &F, g: G)
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar,
              G: Fn(&mut V, RV),
              RV: ThreeDimensional<Scalar=<V::Vector as VectorN>::Scalar> {
        for v in 0 .. self.num_vertices() {
            let normal = self.estimate_normal::<F, RV>(v, f);
            g(self.s.vertex_mut(v), normal);
        }
    }
    /// Estimates gradients for all vertices in this triangulation.
    ///
    /// `f` yields the value for which a gradient should be calculated,
    /// `g` can be used to store the gradient. See `estimate_normals` for a similar example.
    ///
    /// Internally, the normal for each vertex is calculated first. Then, an appropriate
    /// gradient is calculated.
    pub fn estimate_gradients<F, G>(&mut self, f: &F, g: &G) 
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar,
              G: Fn(&mut V, V::Vector) {
        for v in 0 .. self.num_vertices() {
            let gradient = self.estimate_gradient::<F>(v, f);
            g(self.s.vertex_mut(v), gradient);
        }
    }
}

fn to_simple_edge<'a, V>(edge: &EdgeHandle<'a, V>) -> SimpleEdge<V::Vector> 
    where V: HasPosition + 'a,
{
    let from = (*edge.from()).position();
    let to = (*edge.to()).position();
    SimpleEdge::new(from, to)
}

#[cfg(test)]
mod test {
    use super::{FloatDelaunayTriangulation, IntDelaunayTriangulation};
    use cgmath::{Vector2};
    use testutils::*;
    use rand::{SeedableRng, XorShiftRng, Rng};
    use rand::distributions::{Range, IndependentSample};
    use traits::{HasPosition};

    #[test]
    fn test_inserting_one_point() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        assert_eq!(d.num_vertices(), 0);
        d.insert(Vector2::new(0f64, 0f64));
        assert_eq!(d.num_vertices(), 1);
        d.sanity_check();
    }

    #[test]
    fn test_inserting_three_points() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0f64, 0f64));
        d.insert(Vector2::new(1f64, 0f64));
        d.insert(Vector2::new(0f64, 1f64));
        assert_eq!(d.num_vertices(), 3);
        d.sanity_check();
    }

    #[test]
    fn test_inserting_three_points_cw() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0f64, 0f64));
        d.insert(Vector2::new(0f64, 1f64));
        d.insert(Vector2::new(1f64, 0f64));
        assert_eq!(d.num_vertices(), 3);
        d.sanity_check();
    }

    #[test]
    fn test_inserting_four_points() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0f64, 0f64));
        d.insert(Vector2::new(1f64, 0f64));
        d.insert(Vector2::new(0f64, 1f64));
        d.insert(Vector2::new(1f64, 1f64));
        assert_eq!(d.num_vertices(), 4);
        d.sanity_check();
    }

    #[test]
    fn test_iterate_faces() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(-1f64, -1f64));
        d.insert(Vector2::new(1f64, 1f64));
        d.insert(Vector2::new(-1f64, 1f64));
        d.insert(Vector2::new(1f64, -1f64));
        assert_eq!(d.triangles().count(), 2);
        const SIZE: usize = 1000;
        let points = random_points_with_seed::<f64>(SIZE, [2, 3, 112, 2000]);
        for p in points {
            d.insert(p);
        }
        assert_eq!(d.triangles().count(), SIZE * 2 + 2);
        for _ in 0 .. SIZE / 2 {
            d.remove(5);
        }
        assert_eq!(d.triangles().count(), SIZE + 2);
        for t in d.triangles() {
            let adj = t.adjacent_edge().expect("Triangle must have an adjacent edge");
            let next = adj.o_next();
            let prev = adj.o_prev();
            // Check if the edge is really adjacent to a triangle
            assert_eq!(next.o_next(), prev);
            assert_eq!(prev.o_prev(), next);
            assert_eq!(next.to(), prev.from());
        }

    }
    
    #[test]
    fn test_insert_points() {
        // Just check if this won't crash
        const SIZE: usize = 10000;
        let mut points = random_points_with_seed::<f64>(SIZE, [1, 3, 3, 7]);
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
        d.sanity_check();
    }

    #[test]
    fn test_insert_integral_points() {
        const SIZE: usize = 10000;
        let mut points = random_points_in_range(1000i64, SIZE, [100934, 1235, 701, 12355]);
        let mut d = IntDelaunayTriangulation::with_tree_lookup();
        for p in points.drain(..) {
            d.insert(p);
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_outside_convex_hull() {
        const NUM: usize = 100;
        let mut rng = XorShiftRng::from_seed([94, 62, 2010, 2016]);
        let range = Range::new(0., ::std::f64::consts::PI);
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for _ in 0 .. NUM {
            let ang = range.ind_sample(&mut rng);
            let vec = Vector2::new(ang.sin(), ang.cos()) * 100.;
            d.insert(vec);
        }
        assert_eq!(d.num_vertices(), NUM);
        d.sanity_check();
    }


    #[test]
    fn test_insert_same_point_small() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        let points = vec![Vector2::new(0.0, 0.0),
                          Vector2::new(0.0, 1.0),
                          Vector2::new(1.0, 0.0),
                          Vector2::new(0.2, 0.1)];
        for p in &points {
            d.insert(*p);
        }
        for p in &points {
            d.insert(*p);
        }
        assert_eq!(d.num_vertices(), points.len());
        d.sanity_check();
    }

    #[test]
    fn test_insert_same_point() {
        const SIZE: usize = 300;
        let mut points = random_points_with_seed::<f64>(SIZE, [2, 123, 43, 7]);
        let mut d = FloatDelaunayTriangulation::with_walk_lookup();
        for p in &points {
            d.insert(*p);
        }
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
        d.sanity_check();
    }

    #[test]
    fn test_insert_point_on_ch_edge() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0., 0f64));
        d.insert(Vector2::new(1., 0.));
        d.insert(Vector2::new(0., 1.));

        d.insert(Vector2::new(0., 0.4));
        d.sanity_check();
    }        

    #[test]
    fn test_insert_on_edges() {
        // Just check if this won't crash
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0., 0f64));
        d.insert(Vector2::new(1., 0.));
        d.insert(Vector2::new(0., 1.));
        d.insert(Vector2::new(1., 1.));
        d.insert(Vector2::new(0.5, 0.5));
        d.insert(Vector2::new(0.2, 0.2));
        d.insert(Vector2::new(0., 0.4));
        d.insert(Vector2::new(1., 0.5));
        d.insert(Vector2::new(0.5, 1.));
        d.insert(Vector2::new(0.7, 0.));
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0., 1f64));
        for i in -50 .. 50 {
            d.insert(Vector2::new(i as f64, 0.));
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line_2() {
        use super::PositionInTriangulation;
        // This test inserts the line first
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();

        for i in -50 .. 50 {
            d.insert(Vector2::new(i as f64, 0.));
        }
        
        assert_eq!(d.locate(&Vector2::new(10., 12.)),
                   PositionInTriangulation::NoTriangulationPresent);

        for i in -10 .. 10 {
            d.insert(Vector2::new(i as f64, 0.5 * (i as f64)));
        }
        d.sanity_check();
    }


    #[test]
    fn test_insert_points_on_grid() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0., 1f64));
        d.insert(Vector2::new(0.0, 0.0));
        d.insert(Vector2::new(1.0, 0.0));
        for y in 0 .. 20 {
            for x in 0 .. 7 {
                d.insert(Vector2::new(x as f64, y as f64));
            }
        }
        d.sanity_check();
    }

    #[test]
    fn crash_test1() {
        let points = [Vector2::new(-0.47000003, -0.5525),
                      Vector2::new(-0.45499998, -0.055000007),
                      Vector2::new(0.049999952, -0.52),
                      Vector2::new(-0.10310739, -0.37901995),
                      Vector2::new(-0.29053342, -0.20643954),
                      Vector2::new(-0.19144729, -0.42079023)];
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in points.iter().cloned() {
            d.insert(point);
        }
        d.sanity_check();
    }

    #[test]
    fn crash_test2() {
        let points = [
            Vector2 { x: 536f64, y: -1024f64 }, 
            Vector2 { x: 248.00000000000003, y: -1072f64 }, 
            Vector2 { x: 464f64, y: -1036f64 }, 
            Vector2 { x: 616f64, y: -1004f64 }
        ];
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in points.iter() {
            d.insert(*point);
        }
        d.sanity_check();
    }


    struct PointWithHeight {
        point: Vector2<f64>,
        height: f64,
    }

    impl HasPosition for PointWithHeight {
        type Vector = Vector2<f64>;
        fn position(&self) -> Vector2<f64> {
            self.point
        }
    }

    impl PointWithHeight {
        fn new(x: f64, y: f64, height: f64) -> PointWithHeight {
            PointWithHeight {
                point: Vector2::new(x, y),
                height: height
            }
        }
    }

    #[test]
    fn test_natural_neighbor_interpolation() {
        let mut points = vec![
            PointWithHeight::new(0.0, 0.0, 0.0),
            PointWithHeight::new(1.0, 0.0, 0.0),
            PointWithHeight::new(0.0, 1.0, 0.0),
            PointWithHeight::new(1.0, 1.0, 0.0),
            PointWithHeight::new(2.0, 0.0, 0.0),
            PointWithHeight::new(2.0, 1.0, 0.0),
            PointWithHeight::new(3.0, 0.0, 1.0),
            PointWithHeight::new(3.0, 1.0, 1.0),
            PointWithHeight::new(4.0, 0.0, 1.0),
            PointWithHeight::new(4.0, 1.0, 1.0),
        ];
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in points.drain(..) {
            d.insert(point);
        }
        assert_eq!(d.nn_interpolation(&Vector2::new(0.5, 0.5), |p| p.height), Some(0.0));
        assert_eq!(d.nn_interpolation(&Vector2::new(0.2, 0.8), |p| p.height), Some(0.0));
        assert_eq!(d.nn_interpolation(&Vector2::new(3.5, 1.), |p| p.height), Some(1.0));
        assert_eq!(d.nn_interpolation(&Vector2::new(-20., 0.2), |p| p.height), Some(0.0));
        let height = d.nn_interpolation(&Vector2::new(3.2, 0.9), |p| p.height).unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        let height = d.nn_interpolation(&Vector2::new(3.5, 0.5), |p| p.height).unwrap();
        assert!((height - 1.0).abs() < 0.00001);
    }

    #[test]
    fn test_insert_points_with_increasing_distance() {
        use cgmath::InnerSpace;
        let mut points = random_points_with_seed::<f64>(1000, [2, 23, 493, 712]);
        points.sort_by(|p1, p2| p1.magnitude2().partial_cmp(&p2.magnitude2()).unwrap());
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in &points {
            d.insert(*point);
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_grid_with_increasing_distance() {
        // This test inserts points on a grid with increasing distance 
        // from (0., 0.)
        use cgmath::InnerSpace;
        let mut points = Vec::new();
        const SIZE: i64 = 7;
        for x in -SIZE .. SIZE {
            for y in -SIZE .. SIZE {
                let point = Vector2::new(x as f64, y as f64);
                points.push(point);
            }
        }
        points.sort_by(|p1, p2| p1.magnitude2().partial_cmp(&p2.magnitude2()).unwrap());
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in &points {
            d.insert(*point);
        }
        d.sanity_check();
    }

    #[test]
    fn test_remove_in_triangle() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(-1.0, 0.0f64));
        d.insert(Vector2::new(1.0, 0.0f64));
        d.insert(Vector2::new(0.0, 1.0f64));
        let to_remove = d.insert(Vector2::new(0.0, 0.5));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 3);
        // Reinsert the last point, just to see if a crash occurs
        d.insert(Vector2::new(0.0, 0.5));
        d.sanity_check();
    }

    #[test]
    fn test_remove_in_quad() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0.0, 0.0f64));
        d.insert(Vector2::new(1.0, 0.0f64));
        d.insert(Vector2::new(0.0, 1.0f64));
        d.insert(Vector2::new(1.0, 1.0f64));
        let to_remove = d.insert(Vector2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        let to_remove = d.insert(Vector2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        d.insert(Vector2::new(0.5, 0.6));
        d.sanity_check();
    }

    #[test]
    fn seven_point_test() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0.0, 0.0f64));
        d.insert(Vector2::new(1.0, 0.0f64));
        d.insert(Vector2::new(1.0, 1.0f64));
        d.insert(Vector2::new(0.5, 1.5));
        d.insert(Vector2::new(0.0, 1.0f64));
        d.insert(Vector2::new(0.55, 0.5));
        let last = d.insert(Vector2::new(0.4, 0.9));
        d.remove(last);
    }

    #[test]
    fn test_remove_inner() {
        use ::rand::{SeedableRng, Rng};

        let mut points = random_points_with_seed::<f64>(1000, [22, 231, 493, 712]);
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in &points {
            d.insert(*point);
        }
        // Insert an outer quad since we don't want to remove vertices from
        // the convex hull.
        d.insert(Vector2::new(-2.0, -2.0));
        d.insert(Vector2::new(-2.0, 2.0));
        d.insert(Vector2::new(2.0, -2.0));
        d.insert(Vector2::new(2.0, 2.0));
        // Now remove all inner points
        let mut rng = ::rand::XorShiftRng::from_seed([10, 10, 20, 1203031]);
        rng.shuffle(&mut points);
        assert_eq!(d.num_vertices(), 1004);
        for point in &points {
            let handle = d.lookup(point).unwrap().fix();
            d.remove(handle);
        }
        assert_eq!(d.num_vertices(), 4);
        d.sanity_check();
    }
    
    #[test]
    fn test_remove_outer() {
        use cgmath::InnerSpace;
        let mut points = random_points_with_seed::<f64>(1000, [1022, 35611, 2493, 7212]);
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in &points {
            d.insert(*point);
        }
        points.sort_by(|p1, p2| p1.magnitude2().partial_cmp(&p2.magnitude2()).unwrap());
        for point in points[3 ..].iter().rev() {
            let handle = d.lookup(point).unwrap().fix();
            d.remove(handle);
        }
        d.sanity_check();
    }

    #[test]
    fn test_removal_and_insertion() {
        let points = random_points_with_seed::<f64>(1000, [10221, 325611, 20493, 72212]);
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        for point in &points {
            d.insert(*point);
        }
        let mut rng = XorShiftRng::from_seed([92211, 88213, 3992, 2991]);
        for _ in 0 .. 1000 {
            if rng.gen() {
                // Insert new random point
                d.insert(rng.gen());
            } else {
                // Remove random point
                let range = Range::new(0, d.num_vertices());
                let handle = range.ind_sample(&mut rng);
                d.remove(handle);
            }
        }
        d.sanity_check();
    }

    #[test]
    fn test_remove_until_degenerate() {
        let mut d = FloatDelaunayTriangulation::with_tree_lookup();
        d.insert(Vector2::new(0., 0f64));
        d.insert(Vector2::new(1., 0.));
        d.insert(Vector2::new(0., 1.));
        d.insert(Vector2::new(0., 0.5));
        d.insert(Vector2::new(0., 0.25));        
        d.insert(Vector2::new(0., 0.75));        
        assert_eq!(d.num_finite_faces(), 4);
        assert!(d.lookup_and_remove(&Vector2::new(1., 0.)).is_some());
        d.sanity_check();
        assert_eq!(d.num_faces(), 1);
        assert!(d.is_degenerate());
        while d.num_vertices() != 0 {
            d.remove(0);
        }
        assert!(d.is_degenerate());
        d.sanity_check();
        d.insert(Vector2::new(0.5, 0.5));
        d.insert(Vector2::new(0.2, 0.5));
        d.insert(Vector2::new(1.5, 0.0));
        d.sanity_check();
    }
}
