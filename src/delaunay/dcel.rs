// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

/// Handle to a vertex.
///
/// This handle is "fixed", meaning it is intended to be used for
/// mutation (e.g., removing a vertex) or storage (e.g., storing
/// references to vertices for later usage).
pub type FixedVertexHandle = usize;
/// Handle to an edge.
///
/// This handle is "fixed", meaning it is intended to be used
/// for storage. Note that removal operations will invalidate
/// edge handles.
pub type FixedEdgeHandle = usize;
/// Handle to a face.
///
/// This handle is "fixed", meaning it is intended to be used
/// for storage. Note that removal operations will invalidate
/// face handles.
pub type FixedFaceHandle = usize;

pub struct VertexRemovalResult<V> {
    pub updated_vertex: Option<FixedVertexHandle>,
    pub data: V,
}

#[derive(Debug, Copy, Clone)]
struct FaceEntry {
    adjacent_edge: Option<FixedEdgeHandle>,
}

#[derive(Debug, Copy, Clone)]
struct VertexEntry<V> {
    data: V,
    out_edge: Option<FixedEdgeHandle>,
}

impl <V> VertexEntry<V> {
    fn new(data: V) -> VertexEntry<V> {
        VertexEntry {
            data: data,
            out_edge: None,
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
struct HalfEdgeEntry {
    next: FixedEdgeHandle,
    prev: FixedEdgeHandle,
    twin: FixedEdgeHandle,
    origin: FixedVertexHandle,
    face: FixedFaceHandle,
}
    

#[derive(Clone)]
pub struct DCEL<V> {
    vertices: Vec<VertexEntry<V>>,
    faces: Vec<FaceEntry>,
    edges: Vec<HalfEdgeEntry>,
}

impl <V> DCEL<V> {

    pub fn new() -> DCEL<V> {
        DCEL {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: vec![FaceEntry { adjacent_edge: None } ],
        }
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len() / 2
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    pub fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V> {
        VertexHandle::new(self, handle)
    }

    pub fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V> {
        EdgeHandle::new(self, handle)
    }

    pub fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V> {
        FaceHandle::new(self, handle)
    }

    pub fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V {
        &mut self.vertices[handle].data
    }

    pub fn insert_vertex(&mut self, vertex: V) -> FixedVertexHandle {
        self.vertices.push(VertexEntry::new(vertex));
        self.vertices.len() - 1
    }

    pub fn connect_two_isolated_vertices(&mut self, v0: FixedVertexHandle, v1: FixedVertexHandle,
                                         face: FixedFaceHandle) -> FixedEdgeHandle {
        assert!(self.vertices[v0].out_edge.is_none(), "v0 is not isolated");
        assert!(self.vertices[v1].out_edge.is_none(), "v1 is not isolated");
        assert!(self.faces[face].adjacent_edge.is_none(), "face must not contain any adjacent edges");
        let edge_index = self.edges.len();
        let twin_index = edge_index + 1;
        let edge = HalfEdgeEntry {
            next: twin_index,
            prev: twin_index,
            twin: twin_index,
            origin: v0,
            face: face,
        };
        self.edges.push(edge);

        let twin = HalfEdgeEntry {
            next: edge_index,
            prev: edge_index,
            twin: edge_index,
            origin: v1,
            face: face,
        };
        self.edges.push(twin);

        self.vertices[v0].out_edge = Some(edge_index);
        self.vertices[v1].out_edge = Some(twin_index);

        self.faces[face].adjacent_edge = Some(edge_index);

        edge_index
    }

    pub fn connect_edge_to_isolated_vertex(&mut self, prev_handle: FixedEdgeHandle,
                                           vertex: FixedVertexHandle) -> FixedEdgeHandle {
        assert!(self.vertices[vertex].out_edge.is_none(), "Given vertex is not isolated");
        let prev = self.edges[prev_handle];

        let edge_index = self.edges.len();
        let twin_index = edge_index + 1;
        let edge = HalfEdgeEntry {
            next: twin_index,
            prev: prev_handle,
            twin: twin_index,
            origin: self.edges[prev.twin].origin,
            face: prev.face,
        };
        self.edges.push(edge);

        let twin = HalfEdgeEntry {
            next: prev.next,
            prev: edge_index,
            twin: edge_index,
            origin: vertex,
            face: prev.face,
        };
        self.edges.push(twin);

        self.edges[prev_handle].next = edge_index;
        self.edges[prev.next].prev = twin_index;
        
        self.vertices[vertex].out_edge = Some(twin_index);
        edge_index
    }

    pub fn clear_edges_and_faces(&mut self) {
        self.edges.clear();
        self.faces.clear();
        self.faces.push(FaceEntry { adjacent_edge: None });
        for vertex in &mut self.vertices {
            vertex.out_edge = None;
        }
    }

    pub fn remove_vertex(&mut self, vertex_handle: FixedVertexHandle,
                         remaining_face: Option<FixedFaceHandle>) -> VertexRemovalResult<V> {
        while let Some(out_edge) = self.vertices[vertex_handle].out_edge {
            self.remove_edge(out_edge, remaining_face);
        }
        let data = self.vertices.swap_remove(vertex_handle).data;
        let updated_vertex = if self.vertices.len() == vertex_handle {
            None
        } else {
            // Update origin of all out edges
            let to_update: Vec<_> = self.vertex(vertex_handle).ccw_out_edges().map(
                |e| e.fix()).collect();
            for e in to_update {
                self.edges[e].origin = vertex_handle;
            }
            Some(self.vertices.len())
        };

        VertexRemovalResult {
            updated_vertex: updated_vertex,
            data: data,
        }
    }

    fn connect_edge_to_edge(&mut self, prev_edge_handle: FixedEdgeHandle,
                            next_edge_handle: FixedEdgeHandle) -> FixedEdgeHandle {
        let edge_index = self.edges.len();
        let twin_index = edge_index + 1;
        let next_edge = self.edges[next_edge_handle];
        let prev_edge = self.edges[prev_edge_handle];
        let edge = HalfEdgeEntry {
            next: next_edge_handle,
            prev: prev_edge_handle,
            twin: twin_index,
            origin: self.edges[prev_edge.twin].origin,
            face: next_edge.face,
        };
        self.edges.push(edge);

        let twin = HalfEdgeEntry {
            next: prev_edge.next,
            prev: next_edge.prev,
            twin: edge_index,
            origin: next_edge.origin,
            face: next_edge.face,
        };
        self.edges.push(twin);
        
        self.edges[next_edge_handle].prev = edge_index;
        self.edges[prev_edge_handle].next = edge_index;
        self.edges[next_edge.prev].next = twin_index;
        self.edges[prev_edge.next].prev = twin_index;

        edge_index
    }

    pub fn split_edge(&mut self, edge_handle: FixedEdgeHandle, 
                      split_vertex: FixedVertexHandle) -> FixedEdgeHandle {
        assert!(self.vertices[split_vertex].out_edge.is_none(),
                "Given vertex must be isolated");
        let edge = self.edges[edge_handle];
        let twin = self.edges[edge.twin];
        
        let new_edge_index = self.edges.len();
        let new_twin_index = new_edge_index + 1;
        let new_edge = HalfEdgeEntry {
            next: edge.next,
            prev: edge_handle,
            twin: new_twin_index,
            origin: split_vertex,
            face: edge.face,
        };

        let new_twin = HalfEdgeEntry {
            next: edge.twin,
            prev: twin.prev,
            twin: new_edge_index,
            origin: twin.origin,
            face: twin.face,
        };

        self.edges[edge.next].prev = new_edge_index;
        self.edges[twin.prev].next = new_twin_index;

        self.edges[edge_handle].next = new_edge_index;
        self.edges[edge.twin].prev = new_twin_index;

        self.edges[edge.twin].origin = split_vertex;
        self.vertices[twin.origin].out_edge = Some(new_twin_index);
        self.vertices[split_vertex].out_edge = Some(new_edge_index);

        self.edges.push(new_edge);
        self.edges.push(new_twin);
        new_edge_index
    }

    pub fn remove_edge(&mut self, edge_handle: FixedEdgeHandle, 
                       remaining_face: Option<FixedFaceHandle>) {
        let edge = self.edges[edge_handle];

        let twin = self.edges[edge.twin];

        self.edges[edge.prev].next = twin.next;
        self.edges[twin.next].prev = edge.prev;
        self.edges[edge.next].prev = twin.prev;
        self.edges[twin.prev].next = edge.next;

        let (to_remove, to_keep) = if remaining_face == Some(twin.face) {
            (edge, twin)
        } else {
            (twin, edge)
        };
        
        if edge.prev == edge.twin && edge.next == edge.twin {
            // We remove an isolated edge
            self.faces[to_keep.face].adjacent_edge = None;
        } else {
            let new_adjacent_edge = if edge.prev != edge.twin {
                edge.prev
            } else {
                edge.next
            };
            self.faces[to_keep.face].adjacent_edge = Some(new_adjacent_edge);
            self.edges[new_adjacent_edge].face = to_keep.face;
        }

        if edge.prev == edge.twin {
            self.vertices[edge.origin].out_edge = None;
        } else {
            self.vertices[edge.origin].out_edge = Some(twin.next);
        }

        if edge.next == edge.twin {
            self.vertices[twin.origin].out_edge = None;
        } else {
            self.vertices[twin.origin].out_edge = Some(edge.next);
        }
        
        // We must remove the larger index first to prevent the other edge
        // from being updated
        if edge_handle > edge.twin {
            self.swap_out_edge(edge_handle);
            self.swap_out_edge(edge.twin);
        } else {
            self.swap_out_edge(edge.twin);
            self.swap_out_edge(edge_handle);
        }
        if edge.face != twin.face {
            let neighs: Vec<_> = self.face(to_keep.face).adjacent_edges().map(|e| e.fix()).collect();
            for n in neighs {
                self.edges[n].face = to_keep.face
            }
            self.remove_face(to_remove.face);
        }
    }

    fn remove_face(&mut self, face: FixedFaceHandle) {
        self.faces.swap_remove(face);
        if self.faces.len() > face {
            let neighs: Vec<_> = self.face(face).adjacent_edges().map(|e| e.fix()).collect();
            for n in neighs {
                self.edges[n].face = face;
            }
        }
    }

    fn swap_out_edge(&mut self, edge_handle: FixedEdgeHandle) {
        self.edges.swap_remove(edge_handle);
        if self.edges.len() > edge_handle {
            // Update edge index
            let old_handle = self.edges.len();
            let edge = self.edges[edge_handle];
            self.edges[edge.next].prev = edge_handle;
            self.edges[edge.prev].next = edge_handle;
            self.edges[edge.twin].twin = edge_handle;

            if self.vertices[edge.origin].out_edge == Some(old_handle) {
                self.vertices[edge.origin].out_edge = Some(edge_handle);
            }
            self.faces[edge.face].adjacent_edge = Some(edge_handle);

        }
    }

    pub fn create_face(&mut self, 
                       prev_edge_handle: FixedEdgeHandle, 
                       next_edge_handle: FixedEdgeHandle) -> FixedEdgeHandle {
        let edge_index = self.connect_edge_to_edge(prev_edge_handle, next_edge_handle);
        
        let new_face = self.num_faces();

        self.faces.push(FaceEntry {
            adjacent_edge: Some(edge_index)
        });

        // Set the face to the left of the new edge
        let mut cur_edge = edge_index;

        loop {
            self.edges[cur_edge].face = new_face;
            cur_edge = self.edges[cur_edge].next;
            if cur_edge == edge_index {
                break;
            }
        }
        let twin = self.edges[edge_index].twin;
        self.faces[self.edges[twin].face].adjacent_edge = Some(twin);
        edge_index
    }

    pub fn update_vertex(&mut self, handle: FixedVertexHandle, data: V) {
        self.vertices[handle].data = data;
    }

    pub fn edges(&self) -> EdgesIterator<V> {
        EdgesIterator::new(&self)
    }

    pub fn vertices(&self) -> VerticesIterator<V> {
        VerticesIterator::new(&self)
    }

    pub fn fixed_vertices(&self) -> FixedVerticesIterator {
        (0 .. self.num_vertices())
    }

    pub fn faces(&self) -> FacesIterator<V> {
        FacesIterator::new(&self)
    }

    pub fn flip_cw(&mut self, e: FixedEdgeHandle) {
        let en = self.edges[e].next;
        let ep = self.edges[e].prev;
        let t = self.edges[e].twin;
        let tn = self.edges[t].next;
        let tp = self.edges[t].prev;

        self.edges[en].next = e;
        self.edges[en].prev = tp;
        self.edges[e].next = tp;
        self.edges[e].prev = en;
        self.edges[tp].next = en;
        self.edges[tp].prev = e;
        
        self.edges[tn].next = t;
        self.edges[tn].prev = ep;
        self.edges[t].next = ep;
        self.edges[t].prev = tn;
        self.edges[ep].next = tn;
        self.edges[ep].prev = t;

        self.vertices[self.edges[e].origin].out_edge = Some(tn);
        self.vertices[self.edges[t].origin].out_edge = Some(en);

        self.edges[e].origin = self.edges[ep].origin;
        self.edges[t].origin = self.edges[tp].origin;

        self.faces[self.edges[e].face].adjacent_edge = Some(e);
        self.faces[self.edges[t].face].adjacent_edge = Some(t);

        self.edges[tp].face = self.edges[e].face;
        self.edges[ep].face = self.edges[t].face;
    }

    #[doc(hidden)]
    #[allow(dead_code)]
    pub fn print(&self) {
        for (index, edge) in self.edges.iter().enumerate() {
            println!("edge {}: {:#?}", index, edge);
        }
        for (index, vertex) in self.vertices.iter().enumerate() {
            println!("vertex {}: {:?}", index, vertex.out_edge);
        }
        for (index, face) in self.faces.iter().enumerate() {
            println!("face {}: {:?}", index, face);
        }
    }

    #[cfg(test)]
    pub fn sanity_check(&self) {
        for (index, face) in self.faces.iter().enumerate() {
            if let Some(adj) = face.adjacent_edge {
                assert_eq!(self.edges[adj].face, index);
            }
        }
        for (index, vertex) in self.vertices.iter().enumerate() {
            if let Some(out_edge) = vertex.out_edge {
                assert_eq!(self.edges[out_edge].origin, index);
            }
        }
        for edge in self.edges() {
            assert_eq!(edge, edge.o_next().o_prev());
            assert_eq!(edge, edge.o_prev().o_next());
            assert_eq!(edge, edge.sym().sym());
        }
    }
}

/// An iterator that iterates over the edges adjacent to a face.
///
/// The iterator will traverse the edges in oriented order.
/// This order is counterclockwise for right handed coordinate systems
/// or clockwise for left handed systems.
pub struct ONextIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    cur_until: Option<(FixedEdgeHandle, FixedEdgeHandle)>
}

impl <'a, V> ONextIterator<'a, V> where V: 'a {
    fn new_empty(dcel: &'a DCEL<V>) -> ONextIterator<'a, V> {
        ONextIterator {
            dcel: dcel,
            cur_until: None,
        }
    }

    fn new(dcel: &'a DCEL<V>, edge: FixedEdgeHandle) -> ONextIterator<'a, V> {
        let edge = dcel.edge(edge);
        ONextIterator {
            dcel: dcel,
            cur_until: Some((edge.fix(), edge.o_prev().fix())),
        }
    }
}

impl <'a, V> Iterator for ONextIterator<'a, V> where V: 'a {
    type Item = EdgeHandle<'a, V>;

    fn next(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some((cur, until)) = self.cur_until {
            let cur_handle = self.dcel.edge(cur);
            if cur == until {
                self.cur_until = None;
            } else {
                let new_cur = cur_handle.o_next().fix();
                self.cur_until = Some((new_cur, until));
            }
            Some(cur_handle)
        } else {
            None
        }
    }
}

impl <'a, V> DoubleEndedIterator for ONextIterator<'a, V> where V: 'a {
    fn next_back(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some((cur, until)) = self.cur_until {
            let until_handle = self.dcel.edge(until);
            if cur == until {
                self.cur_until = None;
            } else {
                let new_until = until_handle.o_prev().fix();
                self.cur_until = Some((cur, new_until));
            }
            Some(until_handle)
        } else {
            None
        }
    }
}

/// An iterator that iterates over the outgoing edges from a vertex.
///
/// The edges will be iterated in counterclockwise order. Note that
/// this assumes that you use a right handed coordinate system,
/// otherwise the sense of orientation is inverted.
pub struct CCWIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    cur_until: Option<(FixedEdgeHandle, FixedEdgeHandle)>,
}

impl <'a, V> CCWIterator<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>, vertex: FixedVertexHandle) -> CCWIterator<'a, V> {
        let cur_until = if let Some(edge) = dcel.vertex(vertex).out_edge() {
            Some((edge.ccw().fix(), edge.fix()))
        } else {
            None
        };
        CCWIterator {
            dcel: dcel,
            cur_until: cur_until,
        }
    }

    fn from_edge(dcel: &'a DCEL<V>, edge: FixedEdgeHandle) -> CCWIterator<'a, V> {
        let edge = dcel.edge(edge);
        CCWIterator {
            dcel: dcel,
            cur_until: Some((edge.fix(), edge.cw().fix())),
        }
    }
}

impl <'a, V> Iterator for CCWIterator<'a, V> where V: 'a {
    type Item = EdgeHandle<'a, V>;

    fn next(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some((cur, until)) = self.cur_until {
            let cur_handle = self.dcel.edge(cur);
            if cur == until {
                self.cur_until = None;
            } else {
                let new_cur = cur_handle.ccw().fix();
                self.cur_until = Some((new_cur, until));
            }
            Some(cur_handle)
        } else {
            None
        }
    }
}

impl <'a, V> DoubleEndedIterator for CCWIterator<'a, V> where V: 'a {
    fn next_back(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some((cur, until)) = self.cur_until {
            let until_handle = self.dcel.edge(until);
            if cur == until {
                self.cur_until = None;
            } else {
                let new_until = until_handle.cw().fix();
                self.cur_until = Some((cur, new_until));
            }
            Some(until_handle)
        } else {
            None
        }
    }
}


pub struct FacesIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    current: FixedFaceHandle,
}

impl <'a, V> FacesIterator<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>) -> FacesIterator<'a, V> {
        FacesIterator {
            dcel: dcel,
            current: 0,
        }
    }
}

impl <'a, V> Iterator for FacesIterator<'a, V> where V: 'a {
    type Item = FaceHandle<'a, V>;
    
    fn next(&mut self) -> Option<FaceHandle<'a, V>> {
        if self.current < self.dcel.num_faces() {
            let result = FaceHandle::new(self.dcel, self.current);
            self.current += 1;
            Some(result)
        } else {
            None
        }
    }
}

type FixedVerticesIterator = ::std::ops::Range<usize>;

pub struct VerticesIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    current: FixedVertexHandle,
}

impl <'a, V> VerticesIterator<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>) -> VerticesIterator<'a, V> {
        VerticesIterator {
            dcel: dcel,
            current: 0,
        }
    }
}

impl <'a, V> Iterator for VerticesIterator<'a, V> where V: 'a {
    type Item = VertexHandle<'a, V>;
    
    fn next(&mut self) -> Option<VertexHandle<'a, V>> {
        if self.current < self.dcel.num_vertices() {
            let result = VertexHandle::new(self.dcel, self.current);
            self.current += 1;
            Some(result)
        } else {
            None
        }
    }
}

pub struct EdgesIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    current: FixedEdgeHandle,
}

impl <'a, V> EdgesIterator<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>) -> Self {
        EdgesIterator {
            dcel: dcel,
            current: 0
        }
    }
}

impl <'a, V> Iterator for EdgesIterator<'a, V> {
    type Item = EdgeHandle<'a, V>;

    fn next(&mut self) -> Option<EdgeHandle<'a, V>> {
        if let Some(edge) = self.dcel.edges.get(self.current) {
            let twin = edge.twin;
            self.current += 1;
            if self.current - 1 < twin {
                Some(EdgeHandle::new(self.dcel, self.current - 1))
            } else {
                self.next()
            }
        } else {
            None
        }
    }
}

/// A handle to a directed edge.
///
/// This handle has methods that can be used to
/// retrieve the edges neighboring edges, faces
/// and vertices.
pub struct EdgeHandle<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    handle: FixedEdgeHandle,
}

/// A handle to a vertex.
///
/// This handle has methods that can be used to retrieve
/// the vertice's outgoing edges.
pub struct VertexHandle<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    handle: FixedVertexHandle,
}

/// A handle to a face.
///
/// This handle has methods that can be used to get all
/// edges adjacent to this face.
pub struct FaceHandle<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    handle: FixedFaceHandle,
}

impl <'a, V> ::std::fmt::Debug for VertexHandle<'a, V> where V: 'a {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "VertexHandle({:?})", self.handle)
    }
}

impl <'a, V> PartialEq for VertexHandle<'a, V> where V: 'a {
    fn eq(&self, other: &VertexHandle<'a, V>) -> bool {
        self.handle == other.handle
    }
}

impl <'a, V> Copy for VertexHandle<'a, V> where V: 'a { }

impl <'a, V> VertexHandle<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>, handle: FixedVertexHandle) -> VertexHandle<'a, V> {
        VertexHandle {
            dcel: dcel,
            handle: handle,
        }
    }

    /// Returns an outgoing edge.
    ///
    /// If the vertex has multiple outgoing edges, any of them is returned.
    pub fn out_edge(&self) -> Option<EdgeHandle<'a, V>> {
        self.dcel.vertices[self.handle].out_edge.map(|e| self.dcel.edge(e))
    }

    /// Returns all outgoing edges in counter clockwise order.
    ///
    /// Note that this assumes that you use a right handed coordinate system,
    /// otherwise the sense of orientation is inverted.
    pub fn ccw_out_edges(&self) -> CCWIterator<'a, V> {
        CCWIterator::new(self.dcel, self.handle)
    }

    /// Creates a fixed vertex handle from this dynamic handle.
    pub fn fix(&self) -> FixedVertexHandle {
        self.handle
    }
}

impl <'a, V> Clone for VertexHandle<'a, V> where V: 'a {
    fn clone(&self) -> Self {
        VertexHandle::new(self.dcel, self.handle)
    }
}

impl <'a, V> ::std::ops::Deref for VertexHandle<'a, V> {
    type Target = V;
    
    fn deref(&self) -> &V {
        &self.dcel.vertices[self.handle].data
    }
}

impl <'a, V> Copy for EdgeHandle<'a, V> where V: 'a { }

impl <'a, V> Clone for EdgeHandle<'a, V> where V: 'a {
    fn clone(&self) -> Self {
        EdgeHandle::new(self.dcel, self.handle)
    }
}

impl <'a, V> PartialEq for EdgeHandle<'a, V> where V: 'a {
    fn eq(&self, other: &EdgeHandle<'a, V>) -> bool {
        self.handle == other.handle
    }
}

impl <'a, V> ::std::fmt::Debug for EdgeHandle<'a, V> where V: 'a {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "EdgeHandle({:?})", self.handle)
    }
}

pub fn from_neighbors<'a, V> (
    dcel: &'a DCEL<V>,
    from: FixedVertexHandle,
    to: FixedVertexHandle) -> Option<EdgeHandle<'a, V>> {
    let vertex = dcel.vertex(from);
    for edge in vertex.ccw_out_edges() {
        if edge.to().fix() == to {
            return Some(edge);
        }
    }
    None
}

impl <'a, V> EdgeHandle<'a, V> where V: 'a {

    fn new(dcel: &'a DCEL<V>, handle: FixedEdgeHandle) -> Self {
        EdgeHandle {
            dcel: dcel,
            handle: handle,
        }
    }

    /// Creates a fixed vertex handle from this dynamic handle.
    pub fn fix(&self) -> FixedEdgeHandle {
        self.handle
    }

    /// Returns the edge's source vertex.
    pub fn from(&self) -> VertexHandle<'a, V> {
        let edge = self.dcel.edges[self.handle];
        VertexHandle::new(self.dcel, edge.origin)
    }

    /// Returns the oriented next edge.
    ///
    /// The oriented next edge shares the same face as this edge.
    /// When traversing the face's edges in oriented order,
    /// this edge is the predecessor of the oriented next edge.
    /// "Oriented" means counterclockwise for right handed 
    /// coordinate systems.
    pub fn o_next(&self) -> EdgeHandle<'a, V> {
        EdgeHandle::new(self.dcel, self.dcel.edges[self.handle].next)
    }

    /// Returns the oriented previous edge.
    ///
    /// The oriented previous edge shares the same face as this edge.
    /// When traversing the face's edges in oriented order,
    /// this edge is the successor of the oriented previous edge.
    /// "Oriented" means counterclockwise for right handed
    /// coordinate systems.
    pub fn o_prev(&self) -> EdgeHandle<'a, V> {
        EdgeHandle::new(self.dcel, self.dcel.edges[self.handle].prev)
    }

    /// Returns an iterator over all edges sharing the same face
    /// as this edge.
    ///
    /// The face's edges will be traversed in oriented order.
    /// This order is counterclockwise for right handed coordinate
    /// systems or clockwise for left handed systems.
    pub fn o_next_iterator(&self) -> ONextIterator<'a, V> {
        ONextIterator::new(self.dcel, self.handle)
    }

    /// Returns the edges destination vertex.
    pub fn to(&self) -> VertexHandle<'a, V> {
        self.sym().from()
    }
    
    /// Returns the face located to the left of this edge.
    pub fn face(&self) -> FaceHandle<'a, V> {
        self.dcel.face(self.dcel.edges[self.handle].face)
    }

    /// Returns this edge's mirror edge.
    pub fn sym(&self) -> EdgeHandle<'a, V> {
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[self.handle].twin,
        }
    }

    /// Returns the next edge in clockwise direction.
    ///
    /// Note that this assumes that you use a right handed coordinate system,
    /// otherwise the sense of orientation is inverted.
    pub fn cw(&self) -> EdgeHandle<'a, V> {
        let twin = self.sym().handle;
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[twin].next,
        }
    }

    /// Returns the next edge in counter clockwise direction.
    ///
    /// Note that this assumes that you use a right handed coordinate system,
    /// otherwise the sense of orientation is inverted.
    pub fn ccw(&self) -> EdgeHandle<'a, V> {
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[self.handle].prev,
        }.sym()
    }

    /// Returns an iterator over all edges in counter clockwise
    /// order.
    ///
    /// Note that this assumes that you use a right handed coordinate system,
    /// otherwise the sense of orientation is inverted.
    pub fn ccw_iter(&self) -> CCWIterator<'a, V> {
        CCWIterator::from_edge(self.dcel, self.handle)
    }
}

impl <'a, V> Copy for FaceHandle<'a, V> where V: 'a { }

impl <'a, V> Clone for FaceHandle<'a, V> where V: 'a {
    fn clone(&self) -> Self {
        FaceHandle::new(self.dcel, self.handle)
    }
}

impl <'a, V> PartialEq for FaceHandle<'a, V> where V: 'a {
    fn eq(&self, other: &Self) -> bool {
        self.handle == other.handle
    }
}

impl <'a, V> ::std::fmt::Debug for FaceHandle<'a, V> where V: 'a {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "FaceHandle({:?})", self.handle)
    }
}

impl <'a, V> FaceHandle<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>, handle: FixedFaceHandle) -> FaceHandle<'a, V> {
        FaceHandle {
            dcel: dcel,
            handle: handle,
        }
    }

    /// Tries to interpret this face as a triangle, returning its 3 vertices.
    ///
    /// # Panic
    /// This method will panic if the face does not form a triangle.
    pub fn as_triangle(&self) -> [VertexHandle<'a, V>; 3] {
        let adjacent = self.dcel.faces[self.handle].adjacent_edge
            .expect("Face has no adjacent edge");
        let edge = self.dcel.edge(adjacent);
        let prev = edge.o_prev();
        debug_assert!(prev.o_prev() == edge.o_next(),
                      "Face does not form a triangle");
        [prev.from(), edge.from(), edge.to()]
    }

    /// Returns an edge that is adjacent to this face.
    ///
    /// If this face has multiple adjacent edges, any of them is returned.
    pub fn adjacent_edge(&self) -> Option<EdgeHandle<'a, V>> {
        self.dcel.faces[self.handle].adjacent_edge.map(
            |e| EdgeHandle::new(self.dcel, e))
    }

    /// Returns an iterator that iterates over all adjacent edges.
    ///
    /// The edges are traversed in oriented order.
    /// This order will be counterclockwise for right handed coordinate
    /// system or clockwise for left handed systems.
    pub fn adjacent_edges(&self) -> ONextIterator<'a, V> {
        if let Some(adj) = self.dcel.faces[self.handle].adjacent_edge {
            ONextIterator::new(self.dcel, adj)
        } else {
            ONextIterator::new_empty(self.dcel)
        }
    }

    /// Creates a fixed vertex handle from this dynamic vertex handle.
    pub fn fix(&self) -> FixedFaceHandle {
        self.handle
    }
}


#[cfg(test)]
mod test {
    use super::{HalfEdgeEntry, DCEL};

    #[test]
    fn test_create_triangle() {
        let mut dcel = DCEL::new();
        let v0 = dcel.insert_vertex(());
        let v1 = dcel.insert_vertex(());
        let v2 = dcel.insert_vertex(());
        let e01 = dcel.connect_two_isolated_vertices(v0, v1, 0);
        let e12 = dcel.connect_edge_to_isolated_vertex(e01, v2);
        let e20 = dcel.create_face(e12, e01);
        let t01 = dcel.edges[e01].twin;
        let t12 = dcel.edges[e12].twin;
        let t20 = dcel.edges[e20].twin;
        assert_eq!(dcel.edges[e01], 
                   HalfEdgeEntry {
                       next: e12,
                       prev: e20,
                       twin: t01,
                       origin: 0,
                       face: 1,
                   });
        assert_eq!(dcel.edges[e12], 
                   HalfEdgeEntry {
                       next: e20,
                       prev: e01,
                       twin: t12,
                       origin: 1,
                       face: 1,
                   });
        assert_eq!(dcel.edges[e20], 
                   HalfEdgeEntry {
                       next: e01,
                       prev: e12,
                       twin: t20,
                       origin: 2,
                       face: 1,
                   });
        assert_eq!(dcel.edges[t01].face, 0);
        assert_eq!(dcel.edges[t12].face, 0);
        assert_eq!(dcel.edges[t20].face, 0);
    }

    #[test]
    fn test_flip() {
        let mut dcel = DCEL::new();
        let v0 = dcel.insert_vertex(());
        let v1 = dcel.insert_vertex(());
        let v2 = dcel.insert_vertex(());
        let v3 = dcel.insert_vertex(());

        let e01 = dcel.connect_two_isolated_vertices(v0, v1, 0);
        let e12 = dcel.connect_edge_to_isolated_vertex(e01, v2);
        let e23 = dcel.connect_edge_to_isolated_vertex(e12, v3);
        let e30 = dcel.create_face(e23, e01);
        let e_flip = dcel.create_face(e30, e23);
        assert_eq!(dcel.edges[e_flip], 
                   HalfEdgeEntry {
                       next: e23,
                       prev: e30,
                       twin: dcel.edges[e_flip].twin,
                       origin: 0,
                       face: 2,
                   });
        dcel.flip_cw(e_flip);
        let twin = dcel.edges[e_flip].twin;
        assert_eq!(dcel.edges[e_flip],
                   HalfEdgeEntry {
                       next: e12,
                       prev: e23,
                       twin: twin,
                       origin: 3,
                       face: 2,
                   });
        assert_eq!(dcel.edges[twin],
                   HalfEdgeEntry {
                       next: e30,
                       prev: e01,
                       twin: e_flip,
                       origin: 1,
                       face: 1,
                   });
    }

    #[test]
    fn test_split() {
        let mut dcel = DCEL::new();
        let v0 = dcel.insert_vertex(());
        let v1 = dcel.insert_vertex(());
        let v2 = dcel.insert_vertex(());
        let v3 = dcel.insert_vertex(());
        let e01 = dcel.connect_two_isolated_vertices(v0, v1, 0);
        let t01 = dcel.edge(e01).sym().fix();
        let e12 = dcel.connect_edge_to_isolated_vertex(e01, v2);
        let t12 = dcel.edge(e12).sym().fix();
        let e20 = dcel.create_face(e12, e01);
        let t20 = dcel.edge(e20).sym().fix();

        let e_split = dcel.split_edge(e20, v3);
        let t_split = dcel.edge(e_split).sym().fix();
        assert_eq!(dcel.edges[e20], 
                   HalfEdgeEntry {
                       next: e_split,
                       prev: e12,
                       twin: t20,
                       origin: v2,
                       face: 1,
                   });
        assert_eq!(dcel.edges[e_split], 
                   HalfEdgeEntry {
                       next: e01,
                       prev: e20,
                       twin: t_split,
                       origin: v3,
                       face: 1,
                   });
        assert_eq!(dcel.edges[t_split],
                   HalfEdgeEntry {
                       next: t20,
                       prev: t01,
                       origin: v0,
                       twin: e_split,
                       face: 0,
                   });
        assert_eq!(dcel.edges[t20],
                   HalfEdgeEntry {
                       next: t12,
                       prev: t_split,
                       origin: v3,
                       twin: e20,
                       face: 0,
                   });
        assert_eq!(dcel.edges[t01].next, t_split);
        assert_eq!(dcel.edges[e01].prev, e_split);
        assert_eq!(dcel.edges[t12].prev, t20);
        assert_eq!(dcel.edges[e12].next, e20);
        assert!(dcel.vertices[v3].out_edge == Some(e_split)
                || dcel.vertices[v3].out_edge == Some(t20));
        dcel.sanity_check();

    }

    #[test]
    fn test_cw_ccw() {
        let mut dcel = DCEL::new();
        let v0 = dcel.insert_vertex(());
        let v1 = dcel.insert_vertex(());
        let v2 = dcel.insert_vertex(());
        let v3 = dcel.insert_vertex(());

        let e01 = dcel.connect_two_isolated_vertices(v0, v1, 0);
        let e12 = dcel.connect_edge_to_isolated_vertex(e01, v2);
        let e23 = dcel.connect_edge_to_isolated_vertex(e12, v3);
        let e30 = dcel.create_face(e23, e01);
        let e02 = dcel.create_face(e30, e23);

        let e02 = dcel.edge(e02);
        assert_eq!(e02.cw().fix(), e01);
        assert_eq!(e02.ccw().fix(), dcel.edges[e30].twin);
    }

    #[test]
    fn pentagon_test() {
        let mut dcel = DCEL::new();
        let mut v = Vec::new();
        for _ in 0 .. 5 {
            v.push(dcel.insert_vertex(()));
        }

        let e01 = dcel.connect_two_isolated_vertices(v[0], v[1], 0);
        let e12 = dcel.connect_edge_to_isolated_vertex(e01, v[2]);
        let e23 = dcel.connect_edge_to_isolated_vertex(e12, v[3]);
        let e34 = dcel.connect_edge_to_isolated_vertex(e23, v[4]);
        let e40 = dcel.create_face(e34, e01);

        let e02 = dcel.create_face(e40, e23);
        let e03 = dcel.create_face(e40, e34);
        let entry = dcel.edges[e02];
        assert_eq!(entry.next, e23);
        assert_eq!(entry.prev, dcel.edges[e03].twin);
        assert_eq!(entry.origin, v[0]);
    }

    #[test]
    fn test_ccw_iterator() {
        let mut dcel = DCEL::new();
        let mut vs = Vec::new();
        let central = dcel.insert_vertex(());
        assert_eq!(dcel.vertex(central).ccw_out_edges().next(), None);

        for _ in 0 .. 5 {
            vs.push(dcel.insert_vertex(()));
        }
        let mut last_edge = dcel.connect_two_isolated_vertices(central, vs[0], 0);
        last_edge = dcel.edge(last_edge).sym().fix();
        for i in 1 .. 5 {
            last_edge = dcel.connect_edge_to_isolated_vertex(last_edge, vs[i]);
            last_edge = dcel.edge(last_edge).sym().fix();
        }
        
        let out_edge = dcel.vertex(central).out_edge().unwrap();
        let mut neighs: Vec<_> = out_edge.ccw_iter().map(|e| e.to().fix()).collect();
        assert_eq!(neighs.len(), 5);
        for i in 0 .. 5 {
            let first = neighs[i];
            let second = neighs[(i + 1) % 5];
            assert_eq!(first - 1, second % 5);
        }
        let revs: Vec<_> = out_edge.ccw_iter().rev().map(|e| e.to().fix()).collect();
        neighs.reverse();
        assert_eq!(neighs, revs);
    }

    #[test]
    fn test_o_next_iterator() {
        let mut dcel = DCEL::new();
        let mut vs = Vec::new();
        for _ in 0 .. 5 {
            vs.push(dcel.insert_vertex(()));
        }
        
        let mut last_edge = dcel.connect_two_isolated_vertices(vs[0], vs[1], 0);
        let mut edges = vec![last_edge];
        for i in 2 .. 5 {
            last_edge = dcel.connect_edge_to_isolated_vertex(last_edge, vs[i]);
            edges.push(last_edge);
        }
        edges.push(dcel.connect_edge_to_edge(last_edge, vs[0]));
        
        let mut iterated: Vec<_> = dcel.edge(edges[0]).o_next_iterator().map(|e| e.fix()).collect();
        assert_eq!(iterated, edges);

        let rev: Vec<_> = dcel.edge(edges[0]).o_next_iterator().rev().map(|e| e.fix()).collect();
        iterated.reverse();
        assert_eq!(iterated, rev);
    }
}
