// Copyright 2016 The Spade Developers. For a full listing of the authors,
// refer to the Cargo.toml file at the top-level directory of this distribution.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! This module offers a delaunay triangulation and various iterators for inspecting
//! its structure. The main data structure is `DelaunayTriangulation`.
use num::{One, Float, Zero, one, zero};
use traits::{SpatialObject, HasPosition2D, SpadeFloat, HasPosition};
use vector_traits::{VectorN, TwoDimensional, ThreeDimensional};
use kernels::{DelaunayKernel, TrivialKernel};
use rtree::RTree;
use boundingvolume::BoundingRect;
use primitives::{SimpleEdge, SimpleTriangle};
use std::marker::PhantomData;

use planarsubdivision::{PlanarSubdivision, FixedVertexHandle, EdgeHandle, AllEdgesIterator,
                        AllVerticesIterator, VertexHandle, SectorInfo};

/// Yields information about a point's position in triangulation.
#[derive(Debug)]
pub enum PositionInTriangulation<H> {
    /// The point is contained in a triangle. The three given handles are ordered counterclockwise.
    InTriangle([H; 3]),
    /// The point is outside the convex hull. The two given handles mark an edge on the convex
    /// hull that is close to the point.
    OutsideConvexHull(H, H),
    /// An object with this position has already been inserted. Its handle is given.
    OnPoint(H),
    /// The point lies on the edge between two objects. The handles of these objects are given.
    OnEdge(H, H),
    /// There is no valid triangulation yet, thus, all points inserted so far lie on a line.
    NoTriangulationPresent,
}

/// Represents a triangle in a delaunay triangulation.
pub struct DelaunayTriangle<'a, V, K>(pub [VertexHandle<'a, V, K>; 3]) where V: HasPosition2D + 'a, K: 'a, V::Vector: TwoDimensional;

/// Iterator over all triangles in a delaunay triangulation.
pub struct DelaunayTriangleIterator<'a, V, K> 
    where V: HasPosition2D + 'a, 
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar> + 'a,
          V::Vector: TwoDimensional
{
    delaunay: &'a DelaunayTriangulation<V, K>,
    current_index: FixedVertexHandle,
    current_neighbors: Vec<(FixedVertexHandle, FixedVertexHandle)>,
}

impl <'a, V, K> Iterator for DelaunayTriangleIterator<'a, V, K>
    where V: HasPosition2D + 'a, 
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          V::Vector: TwoDimensional {
    type Item = DelaunayTriangle<'a, V, K>;

    fn next(&mut self) -> Option<DelaunayTriangle<'a, V, K>> {
        let ref s = self.delaunay.s;
        loop {
            if let Some((h1, h2)) = self.current_neighbors.pop() {
                let h0 = self.current_index - 1;
                // Check if these neighbors form a triangle
                let h0 = s.handle(h0);
                let h1 = s.handle(h1);
                let h2 = s.handle(h2);
                let v0 = h0.position();
                let v1 = h1.position();
                let v2 = h2.position();
                let edge = SimpleEdge::new(v0, v1);
                if K::side_query(&edge, &v2).is_on_left_side() {
                    return Some(DelaunayTriangle([h0, h1, h2]))
                } else {
                    continue;
                }
            } else {
                let h0 = self.current_index;
                if h0 >= self.delaunay.num_vertices() {
                    return None
                }
                let neighbors = s.handle(h0).fixed_neighbors();
                for i in 0 .. neighbors.len() {
                    let h1 = neighbors[i];
                    if h0 < h1 {
                        let h2 = neighbors[(i + 1) % neighbors.len()];
                        if h0 < h2 {
                            self.current_neighbors.push((h1, h2));
                        }
                    }
                }
                self.current_index += 1;
            }
        }
    }
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

/// An entry of the delaunay triangulation's r-tree.
#[derive(Clone, Debug)]
struct PointEntry<V> {
    point: V,
    handle: FixedVertexHandle,
}

impl<V> HasPosition for PointEntry<V> where V: VectorN {
    type Vector = V;
    fn position(&self) -> V {
        self.point.clone()
    }
}

/// A 2D Delaunay triangulation.
/// 
/// A delaunay triangulation is a special triangulation of a set of points that fulfills some
/// suitable properties for geometric operations like interpolation.
///
/// This triangulation works with `Vector2`-vectors from the `cgmath` and `nalgebra` package.
/// Objects that are inserted into the triangulation have to implement the `HasPosition2D` trait.
/// 
/// Implementing delaunay triangulations is all about precision: the various geometric queries
/// can fail if imprecise calculations are used (like native `f32` / `f64` operations), 
/// resulting in crashes at run time. To prevent those crashes, Spade offers a few "calculation
/// kernels" that can fit the individual needs of an application. Refer to the following table
/// and the documentation of each kernel for more information:
///
/// |  |  Vector types: | When to use: | Properties: |
/// |-------------------|-------------------------------------------|------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
/// | `TrivialKernel` | recommended: f64, i64 also supported: f32, i32 | For i64 coordinates in the range of ± 100000. For f64 coordinates if performance is your  main concern. | Fastest performance. Can crash due to rounding (float types) and overflow (int types) issues. |
/// | `FloatKernel` | f64 | Recommended for f64 coordinates. | Still pretty fast. Uses adaptive precise arithmetic to prevent (theoretically) all crashes. |
/// | `AdaptiveIntKernel` | i64 | For i64 coordinates with large value range. | Slower than `FloatKernel`. Consider casting to floats and using those instead. |
/// # Creation
/// A triangulation can be created with `<kernel_name>::new_triangulation()`, e.g.
///
/// ```
/// # extern crate nalgebra;
/// # extern crate spade;
/// use spade::{DelaunayKernel, FloatKernel};
/// # fn main() {
/// let mut triangulation = FloatKernel::new_triangulation();
/// // Once you insert the first piece of data, type inference is able to fully determine the
/// // triangulation's type - in this case, it contains nalgebra::Vector2<f64>.
/// triangulation.insert(nalgebra::Vector2::new(332f64, 123f64));
/// # }
/// ```
///
/// Also, `Default` is implemented, creating a triangulation with the trivial kernel.
///
/// # Example
///
/// ```
/// extern crate nalgebra;
/// extern crate spade;
///
/// use nalgebra::{Vector2};
/// use spade::{DelaunayTriangulation, DelaunayTriangle};
///
/// # fn main() {
///   let mut delaunay = DelaunayTriangulation::default();
///   delaunay.insert(Vector2::new(0.0, 1.0));
///   delaunay.insert(Vector2::new(0.0, -1.0));
///   delaunay.insert(Vector2::new(1.0, 0.0));
///   delaunay.insert(Vector2::new(-1.0, 0.0));
///   delaunay.insert(Vector2::new(0.0, 0.0));
///   for DelaunayTriangle(vertices) in delaunay.triangles() {
///     println!("found triangle: {:?} -> {:?} -> {:?}", *vertices[0], *vertices[1], *vertices[2]);
///   }
///   for edge in delaunay.edges() {
///     println!("found an edge: {:?} -> {:?}", *edge.from_handle(), *edge.to_handle());
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
pub struct DelaunayTriangulation<V: HasPosition2D, K>
    where V::Vector: TwoDimensional
{
    __kernel: PhantomData<K>,
    s: PlanarSubdivision<V, K>,
    points: RTree<PointEntry<V::Vector>>,
    all_points_on_line: bool,
}

impl<V: HasPosition2D + Clone, K> Clone for DelaunayTriangulation<V, K>  where V::Vector: TwoDimensional {
    fn clone(&self) -> DelaunayTriangulation<V, K> {
        DelaunayTriangulation {
          __kernel: Default::default(),
            s: self.s.clone(),
            points: self.points.clone(),
            all_points_on_line: self.all_points_on_line,
        }
    }
}

impl <V> Default for DelaunayTriangulation<V, TrivialKernel> 
    where V: HasPosition2D,
          V::Vector: TwoDimensional {
    fn default() -> DelaunayTriangulation<V, TrivialKernel> {
        DelaunayTriangulation::new()
    }
}

impl <V, K> DelaunayTriangulation<V, K>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar>,
          V::Vector: TwoDimensional,
{
    /// Creates a new Delaunay triangulation.
    ///
    /// Using this method directly can be a bit cumbersome due to type annotations, consider using
    /// `Default::default()` or `<kernel name>::new_triangulation()` for a `DelaunayKernel` if
    /// possible. Otherwise, you may need to specify the kernel like this:
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// use spade::{DelaunayTriangulation, FloatKernel};
    /// # fn main() {
    /// let mut triangulation = DelaunayTriangulation::<_, FloatKernel>::new();
    /// 
    /// # triangulation.insert(nalgebra::Vector2::new(332f64, 123f64));
    /// # }
    /// ```
    /// 
    /// Usually, the omitted type (the triangulation's vertex type) can be inferred after
    /// `insert` is called.
    pub fn new() -> DelaunayTriangulation<V, K> {
        DelaunayTriangulation {
            __kernel: Default::default(),
            s: PlanarSubdivision::new(),
            points: RTree::new(),
            all_points_on_line: true,
        }
    }

    /// Creates a dynamic vertex handle from a fixed handle.
    /// May panic if the handle was not obtained from this triangulation.
    pub fn handle(&self, handle: FixedVertexHandle) -> VertexHandle<V, K> {
        self.s.handle(handle)
    }

    /// Returns a mutable reference to the vertex data referenced by a 
    /// `FixedVertexHandle`. May panic if the handle was not obtained from this
    /// triangulation.
    pub fn handle_mut(&mut self, handle: FixedVertexHandle) -> &mut V {
        self.s.mut_data(handle)
    }

    /// Checks if the triangulation contains an object with a given coordinate.
    pub fn lookup(&self, point: &V::Vector) -> Option<VertexHandle<V, K>> {
        let handle = self.points.lookup(point);
        handle.map(|h| self.s.handle(h.handle))
    }

    /// Checks if the triangulation contains an object and returns a mutable
    /// reference to it. Note that this will return a reference, while
    /// `lookup(..)` returns a vertex handle.
    pub fn lookup_mut(&mut self, point: &V::Vector) -> Option<&mut V> {
        let handle = self.points.lookup(point);
        if let Some(entry) = handle {
            Some(self.s.mut_data(entry.handle))

        } else {
            None
        }
    }

    /// Returns all vertices contained in a rectangle.
    pub fn lookup_in_rect(&self, rect: &BoundingRect<V::Vector>) -> Vec<VertexHandle<V, K>> {
        let fixed_handles = self.points.lookup_in_rectangle(rect);
        fixed_handles.iter().map(|entry| self.s.handle(entry.handle)).collect()
    }

    /// Returns all vertices contained in a circle.
    pub fn lookup_in_circle(&self, center: &V::Vector, 
                            radius2: &<V::Vector as VectorN>::Scalar) -> Vec<VertexHandle<V, K>> {
        let fixed_handles = self.points.lookup_in_circle(center, radius2);
        fixed_handles.iter().map(|entry| self.s.handle(entry.handle)).collect()
    }

    /// Returns the number of vertices in this triangulation.
    pub fn num_vertices(&self) -> usize {
        self.s.num_vertices()
    }

    /// Returns an iterator over all triangles.
    pub fn triangles(&self) -> DelaunayTriangleIterator<V, K> {
        DelaunayTriangleIterator {
            delaunay: &self,
            current_index: 0,
            current_neighbors: Vec::new(),
        }
    }

    /// Returns an iterator over all edges.
    pub fn edges(&self) -> AllEdgesIterator<V, K> {
        self.s.edges()
    }

    /// Returns an iterator over all vertices.
    pub fn vertices(&self) -> AllVerticesIterator<V, K> {
        self.s.vertices()
    }

    fn initial_insertion(&mut self, t: V) -> Result<FixedVertexHandle, FixedVertexHandle> {
        assert!(self.all_points_on_line);
        // Inserts points if no points are present or if all points
        // lie on the same line
        let new_pos = t.position();
        if let Some(entry) = self.points.lookup(&new_pos) {
            self.s.update_vertex(entry.handle, t);
            return Result::Err(entry.handle);
        }

        if self.points.size() <= 1 {
            return Result::Ok(self.s.insert_vertex(t));
        }

        // Check if the new point is on the same line as all points in the
        // triangulation
        let from = self.s.handle(0).position();
        let to = self.s.handle(1).position();
        let edge = SimpleEdge::new(from.clone(), to.clone());
        if K::side_query(&edge, &new_pos).is_on_line() {
            return Result::Ok(self.s.insert_vertex(t));
        }
        // The point does not lie on the same line as all other points.
        // Start creating a triangulation
        let dir = to.clone() - from.clone();
        let mut vertices: Vec<_> = self.s.vertices().map(
            |v| (v.fix(), dir.dot(&v.position()))).collect();
        // Sort vertices according to their position on the line
        vertices.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());

        // Create line
        for vs in vertices.windows(2) {
            self.s.connect(vs[0].0, vs[1].0);
        }

        let new_vertex = self.s.insert_vertex(t);
        // Connect all points on the line to the new vertex
        for &(v, _) in &vertices {
            self.s.connect(v, new_vertex);
        }
        self.all_points_on_line = false;
        Result::Ok(new_vertex)
    }

    fn get_convex_hull_edges_for_point(&self, first_edge: (FixedVertexHandle, FixedVertexHandle),
                                       point: &V::Vector) -> Vec<FixedVertexHandle> 
    {
        let mut result = Vec::new();
 
        let first_edge = EdgeHandle::from_neighbors(&self.s, first_edge.0, first_edge.1).unwrap();
        debug_assert!(K::side_query(&first_edge.to_simple_edge(), point).is_on_right_side());

        let mut last_edge = first_edge.clone();
        result.push(last_edge.from_handle().fix());
        result.push(last_edge.to_handle().fix());
        // Follow the first edge in cw and ccw direction
        loop {
            let next_edge = last_edge.rev().ccw();
            let query = K::side_query(&next_edge.to_simple_edge(), point);
            if query.is_on_right_side() {
                result.push(next_edge.to_handle().fix());
            } else {
                break;
            }
            last_edge = next_edge;
        }

        last_edge = first_edge.clone();
        loop {
            let next_edge = last_edge.cw().rev();
            let query = K::side_query(&next_edge.to_simple_edge(), (point));
            if query.is_on_right_side() {
                result.insert(0, next_edge.from_handle().fix());
            } else {
                break;
            }
            last_edge = next_edge;
        }
        result
    }

    fn insert_outside_convex_hull(
        &mut self, closest_edge: (FixedVertexHandle, FixedVertexHandle), t: V) -> FixedVertexHandle {
        let position = t.position();
        let handles = self.get_convex_hull_edges_for_point(closest_edge, &position);
        let new_handle = self.s.insert_vertex(t);
        // Make new connections
        let mut illegal_edges = Vec::with_capacity(handles.len() - 1);
        for cur_handle in &handles {
            self.s.connect(*cur_handle, new_handle);
        }
        for from_to in handles.windows(2) {
            illegal_edges.push((from_to[0], from_to[1]));
        }
        self.legalize_edges(illegal_edges, new_handle);
        new_handle
    }
    
    fn insert_into_triangle(&mut self, vertices: [FixedVertexHandle; 3], t: V) -> FixedVertexHandle {
        let new_handle = self.s.insert_vertex(t);
        
        let first_vertex = vertices[0];
        let second_vertex = vertices[1];
        let third_vertex = vertices[2];

        self.s.connect(first_vertex, new_handle);
        self.s.connect(second_vertex, new_handle);
        self.s.connect(third_vertex, new_handle);
        let illegal_edges = vec![(second_vertex, first_vertex),
                                 (third_vertex, second_vertex),
                                 (first_vertex, third_vertex)];
        self.legalize_edges(illegal_edges, new_handle);
        new_handle
    }

    fn get_left_triangle(&self, edge: (FixedVertexHandle, FixedVertexHandle)) 
                         -> Option<FixedVertexHandle> {
        let edge_handle = EdgeHandle::from_neighbors(&self.s, edge.0, edge.1).unwrap();
        let ccw_handle = edge_handle.ccw().to_handle();
        let query = K::side_query(&edge_handle.to_simple_edge(), &ccw_handle.position());
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

    fn insert_on_edge(&mut self, edge: (FixedVertexHandle, FixedVertexHandle), t: V) -> FixedVertexHandle {
        let new_handle = self.s.insert_vertex(t);
        let mut illegal_edges = Vec::new();

        let left_handle_opt = self.get_left_triangle(edge);
        let right_handle_opt = self.get_right_triangle(edge);
        assert!(self.s.disconnect(edge.0, edge.1));
        if let Some(left_handle) = left_handle_opt {
            self.s.connect(left_handle, new_handle);
            illegal_edges.push((left_handle, edge.1));
            illegal_edges.push((edge.0, left_handle));
        }
        if let Some(right_handle) = right_handle_opt {
            self.s.connect(right_handle, new_handle);
            illegal_edges.push((right_handle, edge.0));
            illegal_edges.push((edge.1, right_handle));
        }

        self.s.connect(edge.0, new_handle);
        self.s.connect(edge.1, new_handle);
 
        self.legalize_edges(illegal_edges, new_handle);

        new_handle
    }

    /// Returns information about the location of a point in a triangulation.
    pub fn get_position_in_triangulation(
        &self, point: &V::Vector) -> PositionInTriangulation<VertexHandle<V, K>> {
        use PositionInTriangulation::*;

        match self.get_position_in_triangulation_fixed(point) {
            NoTriangulationPresent => NoTriangulationPresent,
            InTriangle(fs) => InTriangle(
                [self.s.handle(fs[0]), self.s.handle(fs[1]), self.s.handle(fs[2])]),
            OutsideConvexHull(h1, h2) => OutsideConvexHull(self.s.handle(h1), self.s.handle(h2)),
            OnPoint(h) => OnPoint(self.s.handle(h)),
            OnEdge(h1, h2) => OnEdge(self.s.handle(h1), self.s.handle(h2)),
        }
    }

    fn get_position_in_triangulation_fixed(&self, point: &V::Vector) -> PositionInTriangulation<FixedVertexHandle> {
        if self.all_points_on_line {
            // TODO: We might want to check if the point is on the line or already contained
            PositionInTriangulation::NoTriangulationPresent
        } else {
            let start = self.points.close_neighbor(point).unwrap().handle;
            self.get_position_in_triangulation_from_start_point_fixed(start, point) 
        }
    }

    fn get_position_in_triangulation_from_start_point_fixed(
        &self, start: FixedVertexHandle, point: &V::Vector) -> PositionInTriangulation<FixedVertexHandle> {
        let mut cur_handle = self.s.handle(start);
        loop {
            let from_pos = cur_handle.position();
            let sector = match cur_handle.sector_info(point) {
                SectorInfo::NoSector => panic!("Found isolated point. This is a bug."),
                SectorInfo::InSector(sector) => sector
            };
            let neighbors = cur_handle.fixed_neighbors();
            let cw_handle = self.s.handle(neighbors[sector - 1]);
            let ccw_handle = self.s.handle(neighbors[sector % neighbors.len()]);
            let cw_pos = cw_handle.position();
            let ccw_pos = ccw_handle.position();
            let edge = SimpleEdge::new(cw_pos.clone(), ccw_pos.clone());
            if cw_pos == *point {
                return PositionInTriangulation::OnPoint(cw_handle.fix());
            }
            if ccw_pos == *point {
                return PositionInTriangulation::OnPoint(ccw_handle.fix());
            }
            if from_pos == *point {
                return PositionInTriangulation::OnPoint(cur_handle.fix());
            }
            if K::point_on_edge(
                &SimpleEdge::new(from_pos.clone(), cw_pos.clone()), point) {
                return PositionInTriangulation::OnEdge(cw_handle.fix(), cur_handle.fix());
            }
            if K::point_on_edge(
                &SimpleEdge::new(from_pos.clone(), ccw_pos.clone()), point) {
                return PositionInTriangulation::OnEdge(cur_handle.fix(), ccw_handle.fix());
            }

            // Check if the segment is a reflex angle (> 180°)
            let query = K::side_query(&edge, &from_pos);
            if query.is_on_right_side_or_on_line() {
                // The segment forms a reflex angle and is part of the convex hull
                // In some rare cases, we must keep marching towards the point,
                // otherwise we will return an edge e of the convex hull with
                // point being on the right side of e.
                let cw_edge = SimpleEdge::new(from_pos.clone(), cw_pos.clone());
                let ccw_edge = SimpleEdge::new(from_pos.clone(), ccw_pos.clone());

                let cw_edge_feasible = K::side_query(&cw_edge, point).is_on_left_side();
                let ccw_edge_feasible = K::side_query(&ccw_edge, point).is_on_right_side();
                match (cw_edge_feasible, ccw_edge_feasible) {
                    (_, true) => {
                        // Both edges are part of the convex hull, return any of them
                        return PositionInTriangulation::OutsideConvexHull(
                            cur_handle.fix(), ccw_handle.fix());
                    },
                    (true, false) => return PositionInTriangulation::OutsideConvexHull(
                        cw_handle.fix(), cur_handle.fix()),
                    (false, false) => { 
                        // Continue walking
                        let distance_cw = cw_pos.distance2(point);
                        let distance_ccw = ccw_pos.distance2(point);
                        cur_handle = if distance_cw < distance_ccw { cw_handle } else { ccw_handle };
                        continue;
                    },
                }
            }
            // Check if point is contained within the triangle formed by this segment
            if K::side_query(&edge, point).is_on_left_side() {
                // Point lies inside of a triangle
                debug_assert!(EdgeHandle::from_neighbors(&self.s, cw_handle.fix(), ccw_handle.fix()).is_some());
                return PositionInTriangulation::InTriangle([
                    cw_handle.fix(), ccw_handle.fix(), cur_handle.fix()]);
            }
            // We didn't reach the point yet - continue walking
            let distance_cw = cw_pos.distance2(point);
            let distance_ccw = ccw_pos.distance2(point);
            cur_handle = if distance_cw < distance_ccw { cw_handle } else { ccw_handle };
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
        let pos = t.position();
        let insertion_result = match self.get_position_in_triangulation_fixed(&pos) {
            PositionInTriangulation::OutsideConvexHull(v1, v2) => {
                Result::Ok(self.insert_outside_convex_hull((v1, v2), t))
            },
            PositionInTriangulation::InTriangle(vertices) => {
                Result::Ok(self.insert_into_triangle(vertices, t))
            },
            PositionInTriangulation::OnEdge(v1, v2) => {
                Result::Ok(self.insert_on_edge((v1, v2), t))
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
                self.points.insert(PointEntry { point: pos, handle: new_handle });
                new_handle
            },
            Result::Err(update_handle) => {
                update_handle
            }
        }
    }

    fn legalize_edges(&mut self, mut edges: Vec<(FixedVertexHandle, FixedVertexHandle)>, new_vertex: FixedVertexHandle) {
        let position = self.s.handle(new_vertex).position();
         while let Some((h0, h1)) = edges.pop() {
            if let Some(h2) = self.get_left_triangle((h0, h1)) {
                let v0 = self.s.handle(h0).position();
                let v1 = self.s.handle(h1).position();
                let v2 = self.s.handle(h2).position();
                debug_assert!(!K::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(K::is_ordered_ccw(&position, &v1, &v0));
                if K::contained_in_circumference(&v2, &v1, &v0, &position) {
                    // The edge is illegal
                    let edge = EdgeHandle::from_neighbors(&self.s, h0, h1).unwrap().fix();
                    self.s.flip_edge_cw(&edge);
                    edges.push((h2, h1));
                    edges.push((h0, h2));
                    debug_assert!(EdgeHandle::from_neighbors(&self.s, h0, h2).is_some());
                    debug_assert!(EdgeHandle::from_neighbors(&self.s, h1, h2).is_some());
                }
            }
        }
    }
}

impl <V, K> DelaunayTriangulation<V, K> 
    where V: HasPosition2D, <V::Vector as VectorN>::Scalar: SpadeFloat,
          K: DelaunayKernel<<V::Vector as VectorN>::Scalar> ,
          V::Vector: TwoDimensional
{

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
    /// # use spade::{DelaunayTriangulation, DelaunayTriangle, HasPosition};
    ///
    /// struct PointWithHeight {
    ///   point: Vector2<f32>,
    ///   height: f32,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f32>;
    ///     fn position(&self) -> Vector2<f32> {
    ///       self.point
    ///     }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = DelaunayTriangulation::default();
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
            let handle = self.s.handle(*fixed_handle);
            let new = f(&*handle) * ws[index];
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
            let p0 = self.s.handle(nns[0]).position();
            let p1 = self.s.handle(nns[1]).position();
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
            let cur_pos = self.s.handle(*cur).position();
            let next = self.s.handle(nns[(index + 1) % len]).position();
            let triangle = SimpleTriangle::new(next, cur_pos, point.clone());
            point_cell.push(triangle.circumcenter());
        }

        let mut areas = Vec::new();

        for (index, cur) in nns.iter().enumerate() {
            let cur_pos = self.s.handle(*cur).position();
            let prev = nns[((index + len) - 1) % len];
            let next = nns[(index + 1) % len];
            let mut ccw_edge = EdgeHandle::from_neighbors(&self.s, *cur, prev).unwrap();
            let mut polygon = Vec::new();
            polygon.push(point_cell[((index + len) - 1) % len].clone());
            loop {
                if ccw_edge.to_handle().fix() == next {
                    break;
                }
                let cw_edge = ccw_edge.cw();
                let triangle = SimpleTriangle::new(ccw_edge.to_handle().position().clone(), 
                                                   cw_edge.to_handle().position().clone(), cur_pos.clone());
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
        match self.get_position_in_triangulation_fixed(position) {
            PositionInTriangulation::InTriangle(vs) => {
                self.inspect_flips(vec![(vs[0], vs[2]), (vs[2], vs[1]), (vs[1], vs[0])], position)
            },
            PositionInTriangulation::OnEdge(h0, h1) => {
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
            PositionInTriangulation::OutsideConvexHull(h0, h1) => {
                // Get the closest edge on the convex hull
                let edges = self.get_convex_hull_edges_for_point(
                    (h0, h1), position);
                let mut min_dist = <V::Vector as VectorN>::Scalar::max_value();
                let mut min_from = 0;
                let mut min_to = 0;
                for cur_edge in edges.windows(2) {
                    let p0 = self.handle(cur_edge[0]).position();
                    let p1 = self.handle(cur_edge[1]).position();
                    let new_dist = SimpleEdge::new(p0, p1).distance2(position);
                    if new_dist < min_dist {
                        min_from = cur_edge[0];
                        min_to = cur_edge[1];
                        min_dist = new_dist;
                    }
                }
                let p0 = self.handle(min_from).position();
                let p1 = self.handle(min_to).position();
                if SimpleEdge::new(p0.clone(), p1.clone())
                    .is_projection_on_edge(position) {
                    vec![min_from, min_to]
                } else {
                    // Return the closer point
                    if p0.distance2(position) < p1.distance2(position) {
                        vec![min_from]
                    } else {
                        vec![min_to]
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
        let mut result = Vec::new();
        while let Some((h0, h1)) = edges.pop() {
            if let Some(h2) = self.get_left_triangle((h0, h1)) {
                let v0 = self.s.handle(h0).position();
                let v1 = self.s.handle(h1).position();
                let v2 = self.s.handle(h2).position();
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
    /// # use spade::{DelaunayTriangulation, DelaunayTriangle, HasPosition};
    ///
    /// struct PointWithHeight {
    ///   point: Vector2<f32>,
    ///   normal: Vector3<f32>,
    ///   height: f32,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f32>;
    ///     fn position(&self) -> Vector2<f32> {
    ///       self.point
    ///     }
    /// }
    /// impl PointWithHeight {
    ///   fn new(point: Vector2<f32>, height: f32) -> PointWithHeight {
    ///     PointWithHeight { point: point, height: height, normal: Vector3::new(0., 0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = DelaunayTriangulation::default();
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
            g(self.s.mut_data(v), normal);
        }
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
            let v_2d = handle.position();
            v_pos[0] = v_2d[0];
            v_pos[1] = v_2d[1];
            v_pos[2] = f(&*handle);

            handle.neighbors().map(
                |n| {
                    let pos = n.position();
                    let mut result = RV::new();
                    result[0] = pos[0];
                    result[1] = pos[1];
                    result[2] = f(&*n);
                    result
                }).collect()
        };
        let mut final_normal = RV::new();
        for index in 0 .. neighbor_positions.len() {
            let p0 = neighbor_positions[index].clone();
            let p1 = neighbor_positions[
                (index + 1) % neighbor_positions.len()].clone();
            let d0 = v_pos.clone() - p0;
            let d1 = v_pos.clone() - p1;
            let normal = d0.cross(&d1);
            if normal[2] > zero() {
                final_normal = final_normal + normal;
            }

        }
        // Normalize
        final_normal.clone() / final_normal.length2().sqrt()
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
            let gradient = self.estimate_gradient::<F, G>(v, f);
            g(self.s.mut_data(v), gradient);
        }
    }

    /// Estimates and returns the gradient for a single vertex in this triangulation.
    pub fn estimate_gradient<F, G>(&self, v: FixedVertexHandle, f: &F) -> V::Vector 
        where F: Fn(&V) -> <V::Vector as VectorN>::Scalar {
        use cgmath::Vector3;
        let normal = self.estimate_normal::<_, Vector3<_>>(v, f);
        // Calculate gradient from normal
        let mut gradient = V::Vector::new();
        gradient[0] = normal.x;
        gradient[1] = normal.y;
        let g2 = gradient.length2();
        if g2 != zero() {
            let one = <V::Vector as VectorN>::Scalar::one();
            let d2 = one - normal.z * normal.z;
            let a2 = d2 / (one - d2);
            gradient = gradient * (a2 / g2).sqrt();
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
    /// # Example
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// 
    /// # use nalgebra::{Vector2, Vector3};
    /// # use spade::{DelaunayTriangulation, DelaunayTriangle, HasPosition};
    /// struct PointWithHeight {
    ///   point: Vector2<f32>,
    ///   gradient: Vector2<f32>,
    ///   height: f32,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Vector = Vector2<f32>;
    ///     fn position(&self) -> Vector2<f32> {
    ///       self.point
    ///     }
    /// }
    ///
    /// impl PointWithHeight {
    ///   fn new(point: Vector2<f32>, height: f32) -> PointWithHeight {
    ///     // Initialize the gradient to any value since it will be overwritten
    ///     PointWithHeight { point: point, height: height, gradient: Vector2::new(0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = DelaunayTriangulation::default();
    ///   // Insert some points here...
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(0.0, 0.0), 5.));
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(1.0, 0.0), 0.));
    ///   # delaunay.insert(PointWithHeight::new(Vector2::new(0.0, 1.0), 2.));
    ///   // Estimate all gradients and store them:
    ///   delaunay.estimate_gradients(&(|v: &PointWithHeight| v.height),
    ///                               &(|v: &mut PointWithHeight, g| v.gradient = g));
    ///   
    ///   // Now we can use the gradients for interpolation:
    ///  let interpolated = delaunay.nn_interpolation_c1_sibson(
    ///      &Vector2::new(0.5, 0.2), |v| v.height, |v| v.gradient);
    ///  println!("interpolated: {}", interpolated.unwrap());
    /// }
    /// ```
    pub fn nn_interpolation_c1_sibson<F: Fn(&V) -> <V::Vector as VectorN>::Scalar, G: Fn(&V) -> V::Vector> (
        &self, point: &V::Vector, f: F, g: G) -> Option<<V::Vector as VectorN>::Scalar> {
        
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        let mut sum_c0 = zero();
        let mut sum_c1 = zero();
        let mut sum_c1_weights = zero();
        let mut alpha = <V::Vector as VectorN>::Scalar::zero();
        let mut beta = <V::Vector as VectorN>::Scalar::zero();
        for (index, fixed_handle) in nns.iter().enumerate() {
            let handle = self.s.handle(*fixed_handle);
            let pos_i = handle.position();
            let h_i = f(&*handle);
            let diff = pos_i - point.clone();
            let r_i = diff.length2().sqrt();
            let c1_weight_i = ws[index] / r_i;
            let grad_i = g(&*handle);
            let zeta_i = h_i + diff.dot(&grad_i);
            alpha += c1_weight_i * r_i;
            beta += c1_weight_i * r_i * r_i;
            sum_c1_weights += c1_weight_i;
            sum_c1 += zeta_i * c1_weight_i;
            sum_c0 += h_i * ws[index];

        }
        alpha /= sum_c1_weights;
        sum_c1 /= sum_c1_weights;
        let result = (alpha * sum_c0 + beta * sum_c1) / (alpha + beta);
        Some(result)
    }
}

#[cfg(test)]
mod test {
    use super::{DelaunayTriangulation};
    use cgmath::{Vector2};
    use testutils::*;
    use rand::{SeedableRng, XorShiftRng};
    use rand::distributions::{Range, IndependentSample};
    use traits::{HasPosition};
    use kernels::{FloatKernel, DelaunayKernel};

    #[test]
    fn test_inserting_one_point() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        assert_eq!(d.num_vertices(), 0);
        d.insert(Vector2::new(0f32, 0f32));
        assert_eq!(d.num_vertices(), 1);
    }

    #[test]
    fn test_inserting_two_points() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(0f32, 1f32));
        assert_eq!(d.num_vertices(), 2);
    }

    #[test]
    fn test_inserting_three_points() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(1f32, 0f32));
        d.insert(Vector2::new(0f32, 1f32));
        assert_eq!(d.num_vertices(), 3);
    }

    #[test]
    fn test_iterate_faces() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(1f32, 0f32));
        d.insert(Vector2::new(1f32, 1f32));
        d.insert(Vector2::new(2f32, 1f32));
        assert_eq!(d.triangles().count(), 2);
    }
    
    #[test]
    fn test_insert_points() {
        // Just check if this won't crash
        const SIZE: usize = 50000;
        let mut points = random_points_with_seed::<f64>(SIZE, [1, 3, 3, 7]);
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
    }

    #[test]
    fn test_insert_integral_points() {
        const SIZE: usize = 100000;
        let mut points = random_points_in_range(1000i64, SIZE, [100934, 1235, 701, 12355]);
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for p in points.drain(..) {
            d.insert(p);
        }
    }

    #[test]
    fn test_insert_outside_convex_hull() {
        const NUM: usize = 100;
        let mut rng = XorShiftRng::from_seed([94, 62, 2010, 2016]);
        let range = Range::new(0., ::std::f64::consts::PI);
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for _ in 0 .. NUM {
            let ang = range.ind_sample(&mut rng);
            let vec = Vector2::new(ang.sin(), ang.cos()) * 100.;
            d.insert(vec);
        }
        assert_eq!(d.num_vertices(), NUM);
    }

    #[test]
    fn test_insert_same_point() {
        const SIZE: usize = 30;
        let mut points = random_points_with_seed::<f64>(SIZE, [2, 123, 43, 7]);
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for p in &points {
            d.insert(*p);
        }
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
    }

    #[test]
    fn test_insert_on_edges() {
        // Just check if this won't crash
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0., 0f32));
        d.insert(Vector2::new(1., 0.));
        d.insert(Vector2::new(0., 1.));
        d.insert(Vector2::new(1., 1.));
        d.insert(Vector2::new(0.5, 0.5));
        d.insert(Vector2::new(0.2, 0.2));
        d.insert(Vector2::new(0., 0.4));
        d.insert(Vector2::new(1., 0.5));
        d.insert(Vector2::new(0.5, 1.));
        d.insert(Vector2::new(0.7, 0.));
    }

    #[test]
    fn test_insert_points_on_line() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0., 1f32));
        for i in -50 .. 50 {
            d.insert(Vector2::new(i as f32, 0.));
        }
    }

    #[test]
    fn test_insert_points_on_line_2() {
        // This test inserts the line first
        let mut d: DelaunayTriangulation<_, _> = Default::default();

        for i in -50 .. 50 {
            d.insert(Vector2::new(i as f32, 0.));
        }
        
        for i in -10 .. 10 {
            d.insert(Vector2::new(i as f32, 0.5 * (i as f32)));
        }
    }


    #[test]
    fn test_insert_points_on_grid() {
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        d.insert(Vector2::new(0., 1f32));
        d.insert(Vector2::new(0.0, 0.0));
        d.insert(Vector2::new(1.0, 0.0));
        for y in 0 .. 20 {
            for x in 0 .. 7 {
                d.insert(Vector2::new(x as f32, y as f32));
            }
        }
    }

    #[test]
    fn crash_test1() {
        let points = [Vector2::new(-0.47000003, -0.5525),
                      Vector2::new(-0.45499998, -0.055000007),
                      Vector2::new(0.049999952, -0.52),
                      Vector2::new(-0.10310739, -0.37901995),
                      Vector2::new(-0.29053342, -0.20643954),
                      Vector2::new(-0.19144729, -0.42079023)];
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for point in points.iter().cloned() {
            d.insert(point);
        }
    }

    #[test]
    fn crash_test2() {
        let points = [
            Vector2 { x: 536f64, y: -1024f64 }, 
            Vector2 { x: 248.00000000000003, y: -1072f64 }, 
            Vector2 { x: 464f64, y: -1036f64 }, 
            Vector2 { x: 616f64, y: -1004f64 }
        ];
        let mut d = FloatKernel::new_triangulation();
        for point in points.iter() {
            d.insert(*point);
        }
    }


    struct PointWithHeight {
        point: Vector2<f32>,
        height: f32,
    }

    impl HasPosition for PointWithHeight {
        type Vector = Vector2<f32>;
        fn position(&self) -> Vector2<f32> {
            self.point
        }
    }

    impl PointWithHeight {
        fn new(x: f32, y: f32, height: f32) -> PointWithHeight {
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
        let mut d: DelaunayTriangulation<_, _> = Default::default();
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
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for point in &points {
            d.insert(*point);
        }
    }

    #[test]
    fn test_insert_points_on_grid_with_increasing_distance() {
        // This test inserts points on a grid in increasing order from (0., 0.)
        use cgmath::InnerSpace;
        let mut points = Vec::new();
        const SIZE: i32 = 7;
        for x in -SIZE .. SIZE {
            for y in -SIZE .. SIZE {
                let point = Vector2::new(x as f64, y as f64);
                points.push(point);
            }
        }
        points.sort_by(|p1, p2| p1.magnitude2().partial_cmp(&p2.magnitude2()).unwrap());
        let mut d: DelaunayTriangulation<_, _> = Default::default();
        for point in &points {
            d.insert(*point);
        }
    }
}
