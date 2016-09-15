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

use traits::{HasPosition, RTreeFloat, VectorN};
use primitives::SimpleEdge;
use std::default::Default;
use std::ops::Deref;

pub struct PlanarSubdivision<V> where V: HasPosition {
    vertices: Vec<VertexEntry<V>>,
}

pub fn contained_in_circle_segment<V: VectorN>(
    origin: V, cw_edge: V, ccw_edge: V, query_point: V) -> bool 
    where V::Scalar: RTreeFloat
{
    let cw_edge = SimpleEdge::new(origin, cw_edge);
    let ccw_edge = SimpleEdge::new(origin, ccw_edge);

    let cw_info = cw_edge.side_query(query_point);
    let ccw_info = ccw_edge.side_query(query_point);

    let is_cw_left = cw_info.is_on_left_side_or_on_line();
    let is_ccw_right = ccw_info.is_on_right_side_or_on_line();
    // Check if segment forms an angler sharper than 180 deg
    if ccw_edge.side_query(cw_edge.to).is_on_left_side() {
        // More than 180 deg
        is_cw_left || is_ccw_right
    } else {
        // Less than 180 deg
        is_cw_left && is_ccw_right        
    }
}

impl<V: HasPosition> PlanarSubdivision<V> where <V::Vector as VectorN>::Scalar: RTreeFloat {

    pub fn new() -> PlanarSubdivision<V> {
        Default::default()
    }

    pub fn insert_vertex<'a> (&'a mut self, vertex: V) 
                              -> FixedVertexHandle
    {
        let index = self.vertices.len();
        let new_entry = VertexEntry::new(vertex);
        self.vertices.push(new_entry);
        index
    }

    pub fn connect(&mut self, from: FixedVertexHandle, to: FixedVertexHandle)
        -> FixedEdgeHandle
    {
        assert!(from != to);
        let mut to_index = 0;
        for i in 0 .. 2 {
            let (from, to) = if i == 0 { (to, from) } else { (from, to) };
            let from_pos = self.handle(from).position();
            let to_pos = self.handle(to).position();
            // Find sector that contains 'to'
            // We must keep the out edges sorted ccw
            let mut sector;
            {
                let ref neighbors: Vec<_> = self.entry(from).neighbors;
                sector = neighbors.len();
                for i in 1 .. neighbors.len() {
                    let cw_pos = self.handle(neighbors[i - 1]).position();
                    let ccw_pos = self.handle(neighbors[i]).position();
                    if contained_in_circle_segment(from_pos, cw_pos, 
                                                   ccw_pos, to_pos) {
                        sector = i;
                        break;
                    }
                }
            }
            if i == 0 {
                to_index = sector;
            }
            debug_assert!({
                let ref ns = self.entry(from).neighbors;
                if ns.len() >= 2 {
                    let cw_pos = self.handle(ns[sector - 1]).position();
                    let ccw_pos = self.handle(ns[sector % ns.len()]).position();
                    contained_in_circle_segment(from_pos, cw_pos, ccw_pos, to_pos)
                } else {
                    true
                }
            });
            self.mut_entry(from).neighbors.insert(sector, to);
        }
        FixedEdgeHandle::new(from.clone(), to_index)
    }
    
    pub fn disconnect(&mut self, v1: FixedVertexHandle, v2: FixedVertexHandle)
        -> bool {
            let v1_index;
            let v2_index;
            if let Some(edge) = EdgeHandle::from_neighbors(self, v1, v2) {
                v1_index = edge.to_index;
                v2_index = edge.rev().to_index;
            } else {
                return false
            }
            self.mut_entry(v1).neighbors.remove(v1_index);
            self.mut_entry(v2).neighbors.remove(v2_index);
            true
        }

    fn entry(&self, handle: FixedVertexHandle)
             -> &VertexEntry<V> {
        &self.vertices[handle.clone()]
    }

    fn entry_option(&self, handle: FixedVertexHandle)
                    -> Option<&VertexEntry<V>> 
    {
        self.vertices.get(handle)
    }

    pub fn update_vertex(&mut self, handle: FixedVertexHandle,
                         vertex: V) {
        assert_eq!(vertex.position(), self.handle(handle).position());
        self.mut_entry(handle).data = vertex;
    }

    fn mut_entry(&mut self, handle: FixedVertexHandle)
             -> &mut VertexEntry<V> {
        &mut self.vertices[handle.clone()]
    }

    pub fn flip_edge_cw(&mut self, fixed_handle: &FixedEdgeHandle) {
        // Get handles
        let (edge, rev, new1, new2) = {
            let edge_handle = self.edge_handle(fixed_handle);
            let rev = edge_handle.rev();
            let new1 = edge_handle.ccw().to_handle().fix();
            let new2 = rev.ccw().to_handle().fix();
            (edge_handle.fix(), rev.fix(), new1, new2)
        };
        self.mut_entry(edge.from_handle).neighbors.remove(edge.to_index);
        self.mut_entry(rev.from_handle).neighbors.remove(rev.to_index);
        self.connect(new1, new2);
    }

    pub fn handle(&self, fixed_handle: FixedVertexHandle)
        -> VertexHandle<V> {
        VertexHandle::new(&self, fixed_handle.clone())
    }

    pub fn edge_handle(&self, fixed_handle: &FixedEdgeHandle) 
    -> EdgeHandle<V> {
        EdgeHandle::new(&self, fixed_handle.from_handle,
                        fixed_handle.to_index)
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn fixed_vertices(&self) -> Box<Iterator<Item=FixedVertexHandle>> {
        Box::new(0 .. self.vertices.len())
    }

    pub fn vertices<'a>(&'a self) -> AllVerticesIterator<'a, V> {
        AllVerticesIterator { subdiv: self, cur_vertex: 0 }
    }

    pub fn edges<'a>(&'a self) -> AllEdgesIterator<'a, V> {
        AllEdgesIterator::new(self)
    }
}

pub struct AllVerticesIterator<'a, V: HasPosition + 'a> {
    subdiv: &'a PlanarSubdivision<V>,
    cur_vertex: FixedVertexHandle,
}

impl <'a, V: HasPosition> Iterator for AllVerticesIterator<'a, V> 
    where <V::Vector as VectorN>::Scalar: RTreeFloat {
    type Item = VertexHandle<'a, V>;

    fn next(&mut self) -> Option<VertexHandle<'a, V>> {
        if self.cur_vertex >= self.subdiv.num_vertices() {
            None
        } else {
            self.cur_vertex += 1;
            Some(self.subdiv.handle((self.cur_vertex - 1)))
        }
    }
}

pub struct AllEdgesIterator<'a, V: HasPosition + 'a> {
    subdiv: &'a PlanarSubdivision<V>,
    cur_vertex: FixedVertexHandle,
    cur_vertex_index: usize,
}

impl <'a, V: HasPosition> AllEdgesIterator<'a, V> 
    where <V::Vector as VectorN>::Scalar: RTreeFloat
{
    pub fn new(subdiv: &'a PlanarSubdivision<V>) -> Self {
        AllEdgesIterator {
            subdiv: subdiv,
            cur_vertex: 0,
            cur_vertex_index: 0,
        }
    }
}

impl <'a, V: HasPosition + 'a> Iterator for AllEdgesIterator<'a, V> 
    where <V::Vector as VectorN>::Scalar: RTreeFloat
{
    type Item = EdgeHandle<'a, V>;
    
    fn next(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some(cur_entry) = self.subdiv.entry_option(self.cur_vertex) {
            if self.cur_vertex_index >= cur_entry.neighbors.len() {
                self.cur_vertex += 1;
                self.cur_vertex_index = 0;
                // TODO: This can cause a stack overflow for many isolated vertices
                self.next()
            } else {
                let result = EdgeHandle::new(
                    self.subdiv, self.cur_vertex, self.cur_vertex_index);
                self.cur_vertex_index += 1;
                if result.to_handle().fix() > self.cur_vertex {
                    Some(result)
                } else {
                    self.next()
                }
            }
        } else {
            None
        }
    }
}

impl <V: HasPosition> Default for PlanarSubdivision<V> {
    fn default() -> PlanarSubdivision<V> {
        PlanarSubdivision {
            vertices: Vec::new()
        }
    }
}

struct VertexEntry<V> where V: HasPosition {
    data: V,
    neighbors: Vec<FixedVertexHandle>,
}

impl <V: HasPosition> VertexEntry<V> {
    fn new(data: V) -> VertexEntry<V> {
        VertexEntry {
            data: data,
            neighbors: Vec::with_capacity(6),
        }
    }
}

pub struct VertexHandle<'a, V> where V: HasPosition + 'a, <V::Vector as VectorN>::Scalar: RTreeFloat {
    subdiv: &'a PlanarSubdivision<V>,
    fixed: FixedVertexHandle,
}

impl <'a, V> Clone for VertexHandle<'a, V> where V: HasPosition + 'a,
<V::Vector as VectorN>::Scalar: RTreeFloat {
    fn clone(&self) -> VertexHandle<'a, V> {
        VertexHandle::new(self.subdiv, self.fixed)
    }
}

struct CircularIterator<T, F> where 
    for <'a> F: Fn(&'a T) -> T, T: PartialEq {
    cur_handle: Option<T>,
    start: T,
    is_done: bool,
    f: F,
}

impl <T, F> CircularIterator<T, F> where
    for <'a> F: Fn(&T) -> T, T: PartialEq + Clone {

    fn new(start: T, f: F) -> CircularIterator<T, F> {
        CircularIterator {
            cur_handle: None,
            start: start,
            is_done: false,
            f: f,
        }
    }
}

impl <T, F> Iterator for CircularIterator<T, F> where
    for <'a> F: Fn(&T) -> T, T: PartialEq + Clone
{
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.is_done {
            None
        } else {
            let next = {
                let ref cur_handle = self.cur_handle.clone().unwrap_or_else(
                    || self.start.clone());
                (self.f)(cur_handle)
            };
            self.is_done |= next == self.start;
            self.cur_handle = Some(next);
            self.cur_handle.clone()
        }
    }
}

pub struct EdgeHandle<'a, V> where V: HasPosition + 'a {
    subdiv: &'a PlanarSubdivision<V>,
    from_handle: FixedVertexHandle,
    to_index: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FixedEdgeHandle {
    from_handle: FixedVertexHandle,
    to_index: usize
}

impl FixedEdgeHandle {
    fn new(from_handle: FixedVertexHandle, to_index: usize)
        -> FixedEdgeHandle
    {
        FixedEdgeHandle { from_handle: from_handle,
                          to_index: to_index }
    }
}

impl <'a, V> Clone for EdgeHandle<'a, V> where V: HasPosition + 'a, 
<V::Vector as VectorN>::Scalar: RTreeFloat {
    fn clone(&self) -> EdgeHandle<'a, V> {
        EdgeHandle::new(self.subdiv, self.from_handle, self.to_index)
    }
}

impl <'a, V> PartialEq for EdgeHandle<'a, V> where V: HasPosition + 'a {
    fn eq(&self, rhs: &EdgeHandle<'a, V>) -> bool {
        self.from_handle == rhs.from_handle && self.to_index == rhs.to_index
    }
}

impl <'a, V: HasPosition> EdgeHandle<'a, V> where <V::Vector as VectorN>::Scalar: RTreeFloat {

    pub fn from_neighbors(
        subdiv: &'a PlanarSubdivision<V>, from: FixedVertexHandle,
        to: FixedVertexHandle) -> Option<EdgeHandle<'a, V>> {
        let from_handle = subdiv.handle(from);
        if let Some(index) = from_handle.fixed_neighbors_vec().iter()
            .position(|e| *e == to)
        {
            Some(EdgeHandle::new(subdiv, from, index))
        } else {
            None
        }
    }

    fn new(subdiv: &'a PlanarSubdivision<V>, from_handle: FixedVertexHandle,
           to_index: usize) -> EdgeHandle<'a, V> {
        EdgeHandle {
            subdiv: subdiv,
            from_handle: from_handle,
            to_index: to_index,
        }
    }

    pub fn from_handle(&self) -> VertexHandle<'a, V> {
        self.subdiv.handle(self.from_handle)
    }

    pub fn to_handle(&self) -> VertexHandle<'a, V> {
        self.subdiv.handle(self.subdiv.entry(self.from_handle).neighbors[self.to_index])
    }

    pub fn ccw_edges(&self) -> Box<Iterator<Item=EdgeHandle<'a, V>> + 'a> {
        let clone: EdgeHandle<'a, V> = self.clone();
        Box::new(CircularIterator::new(clone, |h| h.ccw()))
    }

    pub fn fix(&self) -> FixedEdgeHandle {
        FixedEdgeHandle {
            from_handle: self.from_handle,
            to_index: self.to_index,
        }
    }


    pub fn ccw(&self) -> EdgeHandle<'a, V> {
        let from_handle = self.from_handle();
        let ccw_index = (self.to_index + 1) % from_handle.num_neighbors();
        EdgeHandle::new(self.subdiv, self.from_handle, ccw_index)
    }

    pub fn cw(&self) -> EdgeHandle<'a, V> {
        let from_handle = self.from_handle();
        let cw_index = if self.to_index == 0 {
            from_handle.num_neighbors() - 1
        } else {
            self.to_index - 1
        };
        EdgeHandle::new(self.subdiv, self.from_handle, cw_index)
    }

    pub fn rev(&self) -> EdgeHandle<'a, V> {
        let to_handle = self.to_handle().fix();
        EdgeHandle::from_neighbors(
            &self.subdiv, to_handle, self.from_handle)
            .expect("Expected symmetrical edge. This is a bug.")
    }

    pub fn to_simple_edge(&self) -> SimpleEdge<V::Vector> {
        SimpleEdge::new(self.from_handle().position(),
                        self.to_handle().position())
    }
}

impl <'a, V: HasPosition> VertexHandle<'a, V> where <V::Vector as VectorN>::Scalar: RTreeFloat {
    fn new(subdiv: &'a PlanarSubdivision<V>,
           fixed: FixedVertexHandle) -> VertexHandle<'a, V> {
        VertexHandle {
            subdiv: subdiv,
            fixed: fixed,
        }
    }

    pub fn fixed_neighbors(&self) -> Box<Iterator<Item=FixedVertexHandle> + 'a> {
        let entry: &'a VertexEntry<V> = self.subdiv.entry(self.fixed);
        Box::new(entry.neighbors.iter().cloned())
    }

    pub fn fixed_neighbors_vec(&self) -> &Vec<FixedVertexHandle> {
        &self.subdiv.entry(self.fixed).neighbors
    }

    pub fn neighbors(&self) -> Box<Iterator<Item=VertexHandle<'a, V>> + 'a> {
        let entry: &'a VertexEntry<V> = self.subdiv.entry(self.fixed);
        let subdiv = self.subdiv;
        Box::new(entry.neighbors.iter().map(move |h| subdiv.handle(*h)))
    }

    pub fn out_edges(&self) -> Box<Iterator<Item=EdgeHandle<'a, V>> + 'a> {
        let num_neighbors = self.num_neighbors();
        let subdiv = self.subdiv;
        let fixed = self.fixed;
        Box::new((0 .. num_neighbors).map(move |index| EdgeHandle::new(
            subdiv, fixed, index)))
    }

    fn num_neighbors(&self) -> usize {
        let entry: &'a VertexEntry<V> = self.subdiv.entry(self.fixed);
        entry.neighbors.len()
    }

    pub fn fix(&self) -> FixedVertexHandle {
        self.fixed
    }
}

impl <'a, V> Deref for VertexHandle<'a, V> where V: HasPosition,
<V::Vector as VectorN>::Scalar: RTreeFloat {
    type Target = V;

    fn deref(&self) -> &V {
        &self.subdiv.entry(self.fixed).data
    }
}

pub type FixedVertexHandle = usize;

#[cfg(test)]
mod test {
    use super::{PlanarSubdivision, FixedVertexHandle, EdgeHandle, CircularIterator,
                contained_in_circle_segment};
    use cgmath::Vector2;

    fn create_subdiv_with_triangle() -> (PlanarSubdivision<Vector2<f32>>,
                                         FixedVertexHandle, FixedVertexHandle, 
                                         FixedVertexHandle) {
        let mut subdiv = PlanarSubdivision::new();
        let h1 = subdiv.insert_vertex(Vector2::new(0.0f32, 0.0));
        let h2 = subdiv.insert_vertex(Vector2::new(1.0f32, 0.0));
        let h3 = subdiv.insert_vertex(Vector2::new(0.0f32, 1.0));
        subdiv.connect(h1, h2);
        subdiv.connect(h2, h3);
        subdiv.connect(h3, h1);
        (subdiv, h1, h2, h3)
    }
    
    #[test]
    fn test_insert_vertex() {
        let mut s = PlanarSubdivision::new();
        assert_eq!(s.num_vertices(), 0);
        let fixed = s.insert_vertex(Vector2::new(0.0f32, 0.0));
        {
            let handle = s.handle(fixed);
            assert_eq!(handle.fixed_neighbors().count(), 0);
            assert_eq!(s.num_vertices(), 1);
        }
        s.insert_vertex(Vector2::new(0.0f32, 0.0));
        assert_eq!(s.num_vertices(), 2);
    }

    #[test]
    fn test_connect_vertices_simple() {
        let mut s = PlanarSubdivision::new();
        let vec1 = Vector2::new(0.0f32, 0.0);
        let vec2 = Vector2::new(0.0f32, 1.0);
        let fixed1 = s.insert_vertex(vec1);
        let fixed2 = s.insert_vertex(vec2);
        s.connect(fixed1, fixed2);
        
        let (handle1, handle2) = (s.handle(fixed1), s.handle(fixed2));
        assert_eq!(handle1.fixed_neighbors().collect::<Vec<_>>(),
                   vec![fixed2]);
        assert_eq!(handle2.fixed_neighbors().collect::<Vec<_>>(),
                   vec![fixed1]);
    }

    #[test]
    fn test_disconnect() {
        let mut s = PlanarSubdivision::new();
        let vec1 = Vector2::new(0.0f32, 0.0);
        let vec2 = Vector2::new(0.0f32, 1.0);
        let fixed1 = s.insert_vertex(vec1);
        let fixed2 = s.insert_vertex(vec2);
        assert!(!s.disconnect(fixed1, fixed2));
        s.connect(fixed1, fixed2);
        assert!(s.disconnect(fixed1, fixed2));
        let (handle1, handle2) = (s.handle(fixed1), s.handle(fixed2));
        assert_eq!(handle1.fixed_neighbors().collect::<Vec<_>>(), Vec::new());
        assert_eq!(handle2.fixed_neighbors().collect::<Vec<_>>(), Vec::new());
        
    }

    #[derive(Clone, PartialEq)]
    struct ListAndIndex {
        list: Vec<i32>,
        index: usize,
    }

    impl ListAndIndex {
        fn next(&self) -> ListAndIndex {
            let next_index = (self.index + 1) % self.list.len();
            ListAndIndex {
                list: self.list.clone(),
                index: next_index,
            }
        }
    }

    #[test]
    fn test_circular_iterator() {
        let l = ListAndIndex {
            list: vec![0i32, 1, 2, 3],
            index: 1,
        };

        let all: Vec<_> = CircularIterator::new(l, |list| list.next()).collect();
        assert_eq!(all.len(), 4);
        let mut last_value: Option<ListAndIndex> = None;
        for l in all.iter() {
            if let Some(ref value) = last_value {
                assert_eq!(l.index, (value.index + 1) % 4);
            }
            last_value = Some(l.clone());
        }
    }

    #[test]
    fn test_connect_vertices_retains_ccw_ordering() {
        // After insertion, edges must still be ordered ccw
        let mut s = PlanarSubdivision::new();
        let v0 = s.insert_vertex(Vector2::new(0.0f32, 0.0));
        let vs = vec![Vector2::new(0.0f32, 1.0),
                      Vector2::new(1.0f32, 0.0),
                      Vector2::new(0.0f32, -1.0),
                      Vector2::new(-1.0f32, 0.0)];

        let vs: Vec<_> = vs.into_iter().map(|v| s.insert_vertex(v)).collect();
        for v in vs.iter() {
            s.connect(v0, *v);
        }
        let neighbors: Vec<_> = s.handle(v0).fixed_neighbors().collect();
        let positions: Vec<_> = vs.iter()
            .map(|v| neighbors.iter().position(|n| n == v).unwrap()).collect();
        assert!((positions[0] + 1) % 4 == positions[3]);
        assert!((positions[3] + 1) % 4 == positions[2]);
        assert!((positions[2] + 1) % 4 == positions[1]);
        assert!((positions[1] + 1) % 4 == positions[0]);

    }

    #[test]
    fn test_out_edges_count() {
        let (mut subdiv, f1, f2, f3) = create_subdiv_with_triangle();
        {
            let h1 = subdiv.handle(f1);
            let h2 = subdiv.handle(f2);
            let h3 = subdiv.handle(f3);
            assert_eq!(h1.out_edges().count(), 2);
            assert_eq!(h2.out_edges().count(), 2);
            assert_eq!(h3.out_edges().count(), 2);
        }
        let f4 = subdiv.insert_vertex(Vector2::new(1., 1.));
        {
            let h4 = subdiv.handle(f4);
            assert_eq!(h4.out_edges().count(), 0);
        }
        subdiv.connect(f1, f4);
        let h1 = subdiv.handle(f1);
        let h4 = subdiv.handle(f4);
        assert_eq!(h1.out_edges().count(), 3);
        assert_eq!(h4.out_edges().count(), 1);
    }

    #[test]
    fn test_edge_from_neighbors() {
        let mut s = PlanarSubdivision::new();
        let f1 = s.insert_vertex(Vector2::new(0f32, 0f32));
        let f2 = s.insert_vertex(Vector2::new(1f32, 0f32));
        let f3 = s.insert_vertex(Vector2::new(0f32, 1f32));
        s.connect(f1, f2);
        s.connect(f2, f3);
        assert!(EdgeHandle::from_neighbors(&s, f1, f3).is_none());
        assert!(EdgeHandle::from_neighbors(&s, f3, f1).is_none());
        assert!(EdgeHandle::from_neighbors(&s, f1, f2).is_some());
        assert!(EdgeHandle::from_neighbors(&s, f2, f1).is_some());
        assert!(EdgeHandle::from_neighbors(&s, f1, f2).is_some());
        assert!(EdgeHandle::from_neighbors(&s, f3, f2).is_some());
        assert!(EdgeHandle::from_neighbors(&s, f2, f3).is_some());

        let edge = EdgeHandle::from_neighbors(&s, f2, f3).unwrap();
        assert_eq!(edge.from_handle().fix(), f2);
        assert_eq!(edge.to_handle().fix(), f3);
    }

    #[test]
    fn test_rev_edge() {
        let (subdiv, f1, _, _) = create_subdiv_with_triangle();
        let h1 = subdiv.handle(f1);
        assert_eq!(h1.out_edges().count(), 2);
        for edge in h1.out_edges() {
            let rev_edge = edge.rev();
            assert_eq!(rev_edge.from_handle().fix(),
                       edge.to_handle().fix());
            assert_eq!(rev_edge.to_handle().fix(),
                       edge.from_handle().fix());
        }
    }

    #[test]
    fn test_flip_edge() {
        let mut subdiv = PlanarSubdivision::new();
        let h0 = subdiv.insert_vertex(Vector2::new(-1f32, -1f32));
        let h1 = subdiv.insert_vertex(Vector2::new(-1f32, 1f32));
        let h2 = subdiv.insert_vertex(Vector2::new(1f32, 1f32));
        let h3 = subdiv.insert_vertex(Vector2::new(1f32, -1f32));
        subdiv.connect(h0, h1);
        subdiv.connect(h1, h2);
        subdiv.connect(h2, h3);
        subdiv.connect(h3, h0);

        subdiv.connect(h0, h2);
        let edge = EdgeHandle::from_neighbors(&subdiv, h0, h2).unwrap().fix();
        subdiv.flip_edge_cw(&edge);
        let flipped = EdgeHandle::from_neighbors(&subdiv, h1, h3);
        assert!(flipped.is_some());
    }

    #[test]
    fn test_contained_in_circle_segment() {
        let v0 = Vector2::new(0f32, 0f32);
        let v1 = Vector2::new(10f32, 0f32);
        let v2 = Vector2::new(0f32, 2f32);
        let v3 = Vector2::new(-1.8f32, 0f32);
        let v4 = Vector2::new(0f32, -15f32);
        
        let t1 = Vector2::new(1f32, 1f32);
        let t2 = Vector2::new(1f32, -1f32);
        let t3 = Vector2::new(-1f32, 1f32);
        let t4 = Vector2::new(-1f32, -1f32);

        assert!(!contained_in_circle_segment(v0, v2, v3, t1));
        assert!(!contained_in_circle_segment(v0, v2, v3, t2));
        assert!(contained_in_circle_segment(v0, v2, v3, t3));
        assert!(!contained_in_circle_segment(v0, v2, v3, t4));


        assert!(contained_in_circle_segment(v0, v1, v3, t1));
        assert!(!contained_in_circle_segment(v0, v1, v3, t2));
        assert!(contained_in_circle_segment(v0, v1, v3, t3));
        assert!(!contained_in_circle_segment(v0, v1, v3, t4));

        assert!(contained_in_circle_segment(v0, v1, v4, t1));
        assert!(!contained_in_circle_segment(v0, v1, v4, t2));
        assert!(contained_in_circle_segment(v0, v1, v4, t3));
        assert!(contained_in_circle_segment(v0, v1, v4, t4));
    }

    #[test]
    fn test_all_edges_iterator() {
        let mut subdiv = PlanarSubdivision::new();
        let h0 = subdiv.insert_vertex(Vector2::new(-1f32, -1f32));
        let h1 = subdiv.insert_vertex(Vector2::new(-1f32, 1f32));
        let h2 = subdiv.insert_vertex(Vector2::new(1f32, 1f32));
        let h3 = subdiv.insert_vertex(Vector2::new(1f32, -1f32));
        subdiv.connect(h0, h1);
        subdiv.connect(h1, h2);
        subdiv.connect(h2, h3);
        subdiv.connect(h3, h0);
        subdiv.connect(h0, h2);

        for i in 0 .. 20 {
            subdiv.insert_vertex(Vector2::new(i as f32, 10.));
        }
        assert_eq!(subdiv.edges().count(), 5);
    }
}

