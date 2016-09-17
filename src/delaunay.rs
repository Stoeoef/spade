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

use nalgebra::{Matrix3, Matrix4};
use num::{one, zero};
use traits::{SpatialObject, HasPosition, RTreeFloat, VectorN, RTreeNum};
use rtree::RTree;
use boundingvolume::BoundingRect;
use primitives::{SimpleEdge, SimpleTriangle};

use planarsubdivision::{PlanarSubdivision, FixedVertexHandle, EdgeHandle,
                        VertexHandle, contained_in_circle_segment};

#[derive(Debug)]
pub enum PositionInTriangulation {
    InTriangle([FixedVertexHandle; 3]),
    OutsideConvexHull(FixedVertexHandle, FixedVertexHandle),
    OnPoint(FixedVertexHandle),
    OnEdge(FixedVertexHandle, FixedVertexHandle),
    NoTriangulationPresent,
}

pub struct DelaunayTriangle<'a, V>(pub [VertexHandle<'a, V>; 3]) where V: HasPosition + 'a;

pub struct DelaunayTriangleIterator<'a, V> where V: HasPosition + 'a {
    delaunay: &'a DelaunayTriangulation<V>,
    current_index: FixedVertexHandle,
    current_neighbors: Vec<(FixedVertexHandle, FixedVertexHandle)>,
}

impl <'a, V> Iterator for DelaunayTriangleIterator<'a, V> where V: HasPosition + 'a {
    type Item = DelaunayTriangle<'a, V>;

    fn next(&mut self) -> Option<DelaunayTriangle<'a, V>> {
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
                if edge.side_query(&v2).is_on_left_side() {
                    return Some(DelaunayTriangle([h0, h1, h2]))
                } else {
                    continue;
                }
            } else {
                let h0 = self.current_index;
                if h0 >= self.delaunay.num_vertices() {
                    return None
                }
                let neighbors = s.handle(h0).fixed_neighbors_vec().clone();
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
    poly: &Vec<V>) -> V::Scalar where V: VectorN {
    let mut sum = zero();
    for vertices in poly[1 .. ].windows(2) {
        let triangle = SimpleTriangle::new(poly[0].clone(), vertices[0].clone(), vertices[1].clone());
        sum = sum + triangle.double_area();
    }
    sum
}

fn determinant_mat3<S>(m: &Matrix3<S>) -> S where S: RTreeNum {
    m[(0, 0)].clone() * m[(1, 1)].clone() * m[(2, 2)].clone() 
        + m[(0, 1)].clone() * m[(1, 2)].clone() * m[(2, 0)].clone()
        + m[(0, 2)].clone() * m[(1, 0)].clone() * m[(2, 1)].clone()
        - m[(0, 2)].clone() * m[(1, 1)].clone() * m[(2, 0)].clone()
        - m[(0, 1)].clone() * m[(1, 0)].clone() * m[(2, 2)].clone()
        - m[(0, 0)].clone() * m[(1, 2)].clone() * m[(2, 1)].clone()
}

fn determinant_mat4<S>(m: &Matrix4<S>) -> S where S: RTreeNum {
   let m0 = Matrix3::new(m[(1, 1)].clone(), m[(2, 1)].clone(), m[(3, 1)].clone(),
                         m[(1, 2)].clone(), m[(2, 2)].clone(), m[(3, 2)].clone(),
                         m[(1, 3)].clone(), m[(2, 3)].clone(), m[(3, 3)].clone());
    let m1 = Matrix3::new(m[(0, 1)].clone(), m[(2, 1)].clone(), m[(3, 1)].clone(),
                          m[(0, 2)].clone(), m[(2, 2)].clone(), m[(3, 2)].clone(),
                          m[(0, 3)].clone(), m[(2, 3)].clone(), m[(3, 3)].clone());
    let m2 = Matrix3::new(m[(0, 1)].clone(), m[(1, 1)].clone(), m[(3, 1)].clone(),
                          m[(0, 2)].clone(), m[(1, 2)].clone(), m[(3, 2)].clone(),
                          m[(0, 3)].clone(), m[(1, 3)].clone(), m[(3, 3)].clone());
    let m3 = Matrix3::new(m[(0, 1)].clone(), m[(1, 1)].clone(), m[(2, 1)].clone(),
                          m[(0, 2)].clone(), m[(1, 2)].clone(), m[(2, 2)].clone(),
                          m[(0, 3)].clone(), m[(1, 3)].clone(), m[(2, 3)].clone());
    
    m[(0, 0)].clone() * determinant_mat3(&m0) -
        m[(1, 0)].clone() * determinant_mat3(&m1) +
        m[(2, 0)].clone() * determinant_mat3(&m2) -
        m[(3, 0)].clone() * determinant_mat3(&m3)
}

fn contained_in_circumference<V>(v1: &V, v2: &V, v3: &V, p: &V) -> bool where V: VectorN {
                          
    let matrix = Matrix4::new(v1[0].clone(), v2[0].clone(), v3[0].clone(), p[0].clone(),
                              v1[1].clone(), v2[1].clone(), v3[1].clone(), p[1].clone(),
                              
                              v1[0].clone() * v1[0].clone() + v1[1].clone() * v1[1].clone(),
                              v2[0].clone() * v2[0].clone() + v2[1].clone() * v2[1].clone(),
                              v3[0].clone() * v3[0].clone() + v3[1].clone() * v3[1].clone(),
                              p[0].clone() * p[0].clone() + p[1].clone() * p[1].clone(),
                              
                              one(), one(), one(), one());

    determinant_mat4(&matrix) < zero()
}

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

pub struct DelaunayTriangulation<V: HasPosition> {
    pub s: PlanarSubdivision<V>,
    points: RTree<PointEntry<V::Vector>>,
    all_points_on_line: bool,
}

impl <V: HasPosition> Default for DelaunayTriangulation<V> {
    fn default() -> DelaunayTriangulation<V> {
        DelaunayTriangulation::new()
    }
}

impl <V: HasPosition> DelaunayTriangulation<V> {
    pub fn new() -> DelaunayTriangulation<V> {
        DelaunayTriangulation {
            s: PlanarSubdivision::new(),
            points: RTree::new(),
            all_points_on_line: true,
        }
    }
    
    pub fn lookup(&self, point: &V::Vector) -> Option<VertexHandle<V>> {
        let handle = self.points.lookup(point);
        handle.map(|h| self.s.handle(h.handle))
    }

    pub fn lookup_in_rect(&self, rect: &BoundingRect<V::Vector>) -> Vec<VertexHandle<V>> {
        let fixed_handles = self.points.lookup_in_rectangle(rect);
        fixed_handles.iter().map(|entry| self.s.handle(entry.handle)).collect()
    }

    pub fn lookup_in_circle(&self, center: &V::Vector, 
                            radius2: &<V::Vector as VectorN>::Scalar) -> Vec<VertexHandle<V>> {
        let fixed_handles = self.points.lookup_in_circle(center, radius2);
        fixed_handles.iter().map(|entry| self.s.handle(entry.handle)).collect()
    }

    pub fn num_vertices(&self) -> usize {
        self.s.num_vertices()
    }

    pub fn triangles(&self) -> DelaunayTriangleIterator<V> {
        DelaunayTriangleIterator {
            delaunay: &self,
            current_index: 0,
            current_neighbors: Vec::new(),
        }
    }

    pub fn subdiv(&self) -> &PlanarSubdivision<V> {
        &self.s
    }

    fn initial_insertion(&mut self, t: V) -> Option<FixedVertexHandle> {
        assert!(self.all_points_on_line);
        // Inserts points if no points are present or if all points
        // lie on the same line
        let new_pos = t.position();
        if let Some(entry) = self.points.lookup(&new_pos) {
            self.s.update_vertex(entry.handle, t);
            return None;
        }

        if self.points.size() <= 1 {
            return Some(self.s.insert_vertex(t));
        }

        // Check if the new point is on the same line as all points in the
        // triangulation
        let from = self.s.handle(0).position();
        let to = self.s.handle(1).position();
        let edge = SimpleEdge::new(from.clone(), to.clone());
        if edge.side_query(&new_pos).is_on_line() {
            return Some(self.s.insert_vertex(t));
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
        Some(new_vertex)
    }

    fn get_convex_hull_edges_for_point(&self, first_edge: (FixedVertexHandle, FixedVertexHandle),
                                       point: &V::Vector) -> Vec<FixedVertexHandle> 
    {
        let mut result = Vec::new();
 
        let first_edge = EdgeHandle::from_neighbors(&self.s, first_edge.0, first_edge.1).unwrap();
        debug_assert!(first_edge.to_simple_edge().side_query(point).is_on_right_side());

        let mut last_edge = first_edge.clone();
        result.push(last_edge.from_handle().fix());
        result.push(last_edge.to_handle().fix());
        // Follow the first edge in cw and ccw direction
        loop {
            let next_edge = last_edge.rev().ccw();
            let query = next_edge.to_simple_edge().side_query(point);
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
            let query = next_edge.to_simple_edge().side_query(point);
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
        let query = edge_handle.to_simple_edge().side_query(&ccw_handle.position());
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

    #[inline(never)]
    pub fn get_position_in_triangulation(&self, point: &V::Vector) -> PositionInTriangulation {
        if self.all_points_on_line {
            PositionInTriangulation::NoTriangulationPresent
        } else {
            let start = self.points.close_neighbor(point).unwrap().handle;
            self.get_position_in_triangulation_from_start_point(start, point) 
        }
    }

    fn get_position_in_triangulation_from_start_point(
        &self, start: FixedVertexHandle, point: &V::Vector) -> PositionInTriangulation {
        let mut cur_handle = self.s.handle(start);
        loop {
            let from_pos = cur_handle.position();
            let (mut sector, cw_handle, ccw_handle);
            {
                let neighbors = cur_handle.fixed_neighbors_vec();
                sector = neighbors.len();
                for i in 1 .. neighbors.len() {
                    let cw_pos = self.s.handle(neighbors[i - 1]).position();
                    let ccw_pos = self.s.handle(neighbors[i]).position();
                    if contained_in_circle_segment(
                        &from_pos, &cw_pos, &ccw_pos, point) {
                        sector = i;
                        break;
                    }
                }
                cw_handle = self.s.handle(neighbors[sector - 1]).clone();
                ccw_handle = self.s.handle(neighbors[sector % neighbors.len()]).clone();
            }
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
            if SimpleEdge::new(from_pos.clone(), cw_pos.clone()).is_point_on_edge(point) {
                return PositionInTriangulation::OnEdge(cw_handle.fix(), cur_handle.fix());
            }
            if SimpleEdge::new(from_pos.clone(), ccw_pos.clone()).is_point_on_edge(point) {
                return PositionInTriangulation::OnEdge(cur_handle.fix(), ccw_handle.fix());
            }

            // Check if the segment is a reflex angle (> 180°)
            let query = edge.side_query(&from_pos);
            if query.is_on_right_side_or_on_line() {
                // The segment forms a reflex angle and is part of the convex hull
                // In some rare cases, we must keep marching towards the point,
                // otherwise we will return an edge e of the convex hull with
                // point being on the right side of e.
                let cw_edge = SimpleEdge::new(from_pos.clone(), cw_pos.clone());
                let ccw_edge = SimpleEdge::new(from_pos.clone(), ccw_pos.clone());

                let cw_edge_feasible = cw_edge.side_query(point).is_on_left_side();
                let ccw_edge_feasible = ccw_edge.side_query(point).is_on_right_side();
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
            if edge.side_query(point).is_on_left_side() {
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

    pub fn insert(&mut self, t: V) {
        let pos = t.position();
        let new_handle_opt = match self.get_position_in_triangulation(&pos) {
            PositionInTriangulation::OutsideConvexHull(v1, v2) => {
                Some(self.insert_outside_convex_hull((v1, v2), t))
            },
            PositionInTriangulation::InTriangle(vertices) => {
                Some(self.insert_into_triangle(vertices, t))
            },
            PositionInTriangulation::OnEdge(v1, v2) => {
                Some(self.insert_on_edge((v1, v2), t))
            },
            PositionInTriangulation::OnPoint(vertex) => {
                self.s.update_vertex(vertex, t);
                None
            },
            PositionInTriangulation::NoTriangulationPresent => {
                self.initial_insertion(t)
            }                
        };
        if let Some(new_handle) = new_handle_opt {
            self.points.insert(PointEntry { point: pos, handle: new_handle });
        }
    }

    fn legalize_edges(&mut self, mut edges: Vec<(FixedVertexHandle, FixedVertexHandle)>, new_vertex: FixedVertexHandle) {
        let position = self.s.handle(new_vertex).position();
         while let Some((h0, h1)) = edges.pop() {
            if let Some(h2) = self.get_left_triangle((h0, h1)) {
                let v0 = self.s.handle(h0).position();
                let v1 = self.s.handle(h1).position();
                let v2 = self.s.handle(h2).position();
                debug_assert!(!SimpleTriangle::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(SimpleTriangle::is_ordered_ccw(&position, &v1, &v0));
                if contained_in_circumference(&v2, &v1, &v0, &position) {
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

impl <V> DelaunayTriangulation<V> where V: HasPosition, <V::Vector as VectorN>::Scalar: RTreeFloat {

    pub fn nn_interpolation<F: Fn(&V) -> <V::Vector as VectorN>::Scalar>(
        &self, point: &V::Vector, f: F) -> Option<<V::Vector as VectorN>::Scalar> {

        let interpolate_on_edge = |h0, h1| {
            let h0 = self.s.handle(h0);
            let h1 = self.s.handle(h1);
            let edge = SimpleEdge::new(h0.position(), h1.position());
            let w1 = ::clamp::clamp(zero(), edge.project_point(point), one());
            let (f0, f1) = (f(&*h0), f(&*h1));
            f0 * (one::<<V::Vector as VectorN>::Scalar>() - w1) + f1 * w1
        };

        let nns = match self.get_position_in_triangulation(point) {
            PositionInTriangulation::InTriangle(vs) => {
                self.inspect_flips(vec![(vs[0], vs[2]), (vs[2], vs[1]), (vs[1], vs[0])], point)
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
                    self.inspect_flips(direct_neighbors, point)
                } else {
                    return Some(interpolate_on_edge(h0, h1))
                }
            },
            PositionInTriangulation::OutsideConvexHull(h0, h1) => {
                return Some(interpolate_on_edge(h0, h1));
            },
            PositionInTriangulation::OnPoint(fixed_handle) => {
                let handle = self.s.handle(fixed_handle);
                return Some(f(&*handle));
            },
            _ => return None,
        };
        let ws = self.get_weights(&nns, point);
        let mut sum = zero();
        for (index, fixed_handle) in nns.iter().enumerate() {
            let handle = self.s.handle(*fixed_handle);
            sum += ws[index] * f(&*handle);
        }
        Some(sum)
    }
    
    fn get_weights(&self, nns: &Vec<FixedVertexHandle>, 
                   point: &V::Vector) -> Vec<<V::Vector as VectorN>::Scalar> {
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

    fn inspect_flips(&self, mut edges: Vec<(FixedVertexHandle, FixedVertexHandle)>, position: &V::Vector) -> Vec<FixedVertexHandle> {
        let mut result = Vec::new();
        while let Some((h0, h1)) = edges.pop() {
            if let Some(h2) = self.get_left_triangle((h0, h1)) {
                let v0 = self.s.handle(h0).position();
                let v1 = self.s.handle(h1).position();
                let v2 = self.s.handle(h2).position();
                debug_assert!(!SimpleTriangle::is_ordered_ccw(&v2, &v1, &v0));
                debug_assert!(SimpleTriangle::is_ordered_ccw(position, &v1, &v0));
                if contained_in_circumference(&v2, &v1, &v0, position) {
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
}

#[cfg(test)]
mod test {
    use super::{DelaunayTriangulation, contained_in_circumference};
    use cgmath::{Vector2, Array};
    use rand::{Rand, SeedableRng, XorShiftRng};
    use rand::distributions::{Range, IndependentSample};
    use rand::distributions::range::SampleRange;
    use traits::{HasPosition, RTreeNum};

    #[test]
    fn test_inserting_one_point() {
        let mut d = DelaunayTriangulation::new();
        assert_eq!(d.num_vertices(), 0);
        d.insert(Vector2::new(0f32, 0f32));
        assert_eq!(d.num_vertices(), 1);
    }

    #[test]
    fn test_inserting_two_points() {
        let mut d = DelaunayTriangulation::new();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(0f32, 1f32));
        assert_eq!(d.num_vertices(), 2);
    }

    #[test]
    fn test_inserting_three_points() {
        let mut d = DelaunayTriangulation::new();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(1f32, 0f32));
        d.insert(Vector2::new(0f32, 1f32));
        assert_eq!(d.num_vertices(), 3);
    }

    #[test]
    fn test_iterate_faces() {
        let mut d = DelaunayTriangulation::new();
        d.insert(Vector2::new(0f32, 0f32));
        d.insert(Vector2::new(1f32, 0f32));
        d.insert(Vector2::new(1f32, 1f32));
        d.insert(Vector2::new(2f32, 1f32));
        assert_eq!(d.triangles().count(), 2);
    }
    
    // #[test]
    // fn test_triangle_entry_partial_eq() {
    //     let v0 = Vector2::new(0.0f32, 0.0);
    //     let v1 = Vector2::new(1.0, 0.0);
    //     let v2 = Vector2::new(0.0, 2.0);
    //     let v3 = Vector2::new(1.0, 1.0);
    //     let t1 = TriangleEntry::new([v0, v1, v2], [0, 1, 2]);
    //     let t2 = TriangleEntry::new([v2, v0, v1], [2, 0, 1]);
    //     assert!(t1 == t2);
    //     let t3 = TriangleEntry::new([v1, v2, v0], [1, 2, 0]);
    //     assert!(t1 == t3);
    //     let t4 = TriangleEntry::new([v0, v1, v3], [0, 1, 2]);
    //     assert!(t1 != t4);
    //     let t5 = TriangleEntry::new([v2, v0, v1], [3, 0, 1]);
    //     assert!(t1 != t5);
    // }

    #[test]
    fn test_contained_in_circumference() {
        let (a1, a2, a3) = (1f32, 2f32, 3f32);
        let offset = Vector2::new(0.5, 0.7);
        let v1 = Vector2::new(a1.sin(), a1.cos()) * 2. + offset;
        let v2 = Vector2::new(a2.sin(), a2.cos()) * 2. + offset;
        let v3 = Vector2::new(a3.sin(), a3.cos()) * 2. + offset;
        assert!(contained_in_circumference(&v1, &v2, &v3, &offset));
        let shrunk = (v1 - offset) * 0.9 + offset;
        assert!(contained_in_circumference(&v1, &v2, &v3, &shrunk));
        let expanded = (v1 - offset) * 1.1 + offset;
        assert!(!contained_in_circumference(&v1, &v2, &v3, &expanded));
        assert!(!contained_in_circumference(&v1, &v2, &v3, &(Vector2::from_value(2.0) + offset)));
        assert!(contained_in_circumference(&Vector2::new(0f32, -1f32), &Vector2::new(0f32, 0f32),
                                           &Vector2::new(-1f32, 0f32), &Vector2::new(0f32, 1f32)));
    }

    fn random_points_in_range<S: RTreeNum + Rand + SampleRange>(range: S, size: usize, seed: [u32; 4]) -> Vec<Vector2<S>> {
        // let mut rng = XorShiftRng::from_seed([1, 3, 3, 7]);
        let mut rng = XorShiftRng::from_seed(seed);
        let range = Range::new(-range.clone(), range.clone());
        let mut points = Vec::with_capacity(size);
        for _ in 0 .. size {
            let x = range.ind_sample(&mut rng);
            let y = range.ind_sample(&mut rng);
            points.push(Vector2::new(x, y));
        }
        points
    }

    fn random_points(size: usize, seed: [u32; 4]) -> Vec<Vector2<f64>> {
        random_points_in_range(500., size, seed)
    }

    #[test]
    fn test_insert_points() {
        // Just check if this won't crash
        const SIZE: usize = 200000;
        let mut points = random_points(SIZE, [1, 3, 3, 7]);
        let mut delaunay = DelaunayTriangulation::new();
        for p in points.drain(..) {
            delaunay.insert(p);
        }
        assert_eq!(delaunay.num_vertices(), SIZE);
    }

    #[test]
    fn test_insert_integral_points() {
        const SIZE: usize = 100000;
        let mut points = random_points_in_range(1000i64, SIZE, [100934, 1235, 701, 12355]);
        let mut delaunay = DelaunayTriangulation::new();
        for p in points.drain(..) {
            delaunay.insert(p);
        }
    }

    #[test]
    fn test_insert_outside_convex_hull() {
        const NUM: usize = 100;
        let mut rng = XorShiftRng::from_seed([94, 62, 2010, 2016]);
        let range = Range::new(0., ::std::f64::consts::PI);
        let mut delaunay = DelaunayTriangulation::new();
        for _ in 0 .. NUM {
            let ang = range.ind_sample(&mut rng);
            let vec = Vector2::new(ang.sin(), ang.cos()) * 100.;
            delaunay.insert(vec);
        }
        assert_eq!(delaunay.num_vertices(), NUM);
    }

    #[test]
    fn test_insert_same_point() {
        const SIZE: usize = 30;
        let mut points = random_points(SIZE, [2, 123, 43, 7]);
        let mut delaunay = DelaunayTriangulation::new();
        for p in &points {
            delaunay.insert(*p);
        }
        for p in points.drain(..) {
            delaunay.insert(p);
        }
        assert_eq!(delaunay.num_vertices(), SIZE);
    }

    #[test]
    fn test_insert_on_edges() {
        // Just check if this won't crash
        let mut delaunay = DelaunayTriangulation::new();
        delaunay.insert(Vector2::new(0., 0f32));
        delaunay.insert(Vector2::new(1., 0.));
        delaunay.insert(Vector2::new(0., 1.));
        delaunay.insert(Vector2::new(1., 1.));
        delaunay.insert(Vector2::new(0.5, 0.5));
        delaunay.insert(Vector2::new(0.2, 0.2));
        delaunay.insert(Vector2::new(0., 0.4));
        delaunay.insert(Vector2::new(1., 0.5));
        delaunay.insert(Vector2::new(0.5, 1.));
        delaunay.insert(Vector2::new(0.7, 0.));
    }

    #[test]
    fn test_insert_points_on_line() {
        let mut delaunay = DelaunayTriangulation::new();
        delaunay.insert(Vector2::new(0., 1f32));
        for i in -50 .. 50 {
            delaunay.insert(Vector2::new(i as f32, 0.));
        }
    }

    #[test]
    fn test_insert_points_on_line_2() {
        // This test inserts the line first
        let mut delaunay = DelaunayTriangulation::new();

        for i in -50 .. 50 {
            delaunay.insert(Vector2::new(i as f32, 0.));
        }
        
        for i in -10 .. 10 {
            delaunay.insert(Vector2::new(i as f32, 0.5 * (i as f32)));
        }
    }


    #[test]
    fn test_insert_points_on_grid() {
        let mut delaunay = DelaunayTriangulation::new();
        delaunay.insert(Vector2::new(0., 1f32));
        delaunay.insert(Vector2::new(0.0, 0.0));
        delaunay.insert(Vector2::new(1.0, 0.0));
        for y in 0 .. 20 {
            for x in 0 .. 7 {
                delaunay.insert(Vector2::new(x as f32, y as f32));
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
        let mut delaunay = DelaunayTriangulation::new();
        for point in points.iter().cloned() {
            delaunay.insert(point);
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
        let mut delaunay = DelaunayTriangulation::new();
        for point in points.drain(..) {
            delaunay.insert(point);
        }
        assert_eq!(delaunay.nn_interpolation(&Vector2::new(0.5, 0.5), |p| p.height), Some(0.0));
        assert_eq!(delaunay.nn_interpolation(&Vector2::new(0.2, 0.8), |p| p.height), Some(0.0));
        assert_eq!(delaunay.nn_interpolation(&Vector2::new(3.5, 1.), |p| p.height), Some(1.0));
        assert_eq!(delaunay.nn_interpolation(&Vector2::new(-20., 0.2), |p| p.height), Some(0.0));
        let height = delaunay.nn_interpolation(&Vector2::new(3.2, 0.9), |p| p.height).unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        let height = delaunay.nn_interpolation(&Vector2::new(3.5, 0.5), |p| p.height).unwrap();
        assert!((height - 1.0).abs() < 0.00001);
    }

    #[test]
    fn test_insert_points_with_increasing_distance() {
        use cgmath::InnerSpace;
        let mut points = random_points(1000, [2, 23, 493, 712]);
        points.sort_by(|p1, p2| p1.magnitude2().partial_cmp(&p2.magnitude2()).unwrap());
        let mut delaunay = DelaunayTriangulation::new();
        for point in &points {
            delaunay.insert(*point);
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
        let mut delaunay = DelaunayTriangulation::new();
        for point in &points {
            delaunay.insert(*point);
        }
    }
}
