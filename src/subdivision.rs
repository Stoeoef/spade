// Copyright 2016 The RTree Developers. For a full listing of the authors,
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

use num::{Float};
use traits::PointObject;
use cgmath::{Vector2, EuclideanVector, Angle, Rad};

pub struct EdgeIterator<'a, V> where V: 'a {
    start: FixedEdgeHandle,
    current: FixedEdgeHandle,
    subdivision: &'a Subdivision2<V>,
    finished: bool,
    f: Box<Fn(&Subdivision2<V>, &FixedEdgeHandle) -> FixedEdgeHandle>,
}

impl <'a, V> EdgeIterator<'a, V> {
    fn new(subdivision: &'a Subdivision2<V>, start: FixedEdgeHandle,
           f: Box<Fn(&Subdivision2<V>, &FixedEdgeHandle) -> FixedEdgeHandle>) -> EdgeIterator<'a, V> {
        EdgeIterator {
            subdivision: subdivision,
            current: start,
            start: start,
            finished: false,
            f: f,
        }
    }
}

impl <'a, V> Iterator for EdgeIterator<'a, V> {
    type Item = FixedEdgeHandle;
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let next = (self.f)(self.subdivision, &self.current);
        let old_current = ::std::mem::replace(&mut self.current, next);
        if self.current == self.start {
            self.finished = true;
        }
        Some(old_current)
    }
}

#[derive(Clone, Debug)]
pub struct FixedFaceHandle(FixedEdgeHandle);

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct FixedEdgeHandle {
    e: usize,
    r: usize,
}

impl FixedEdgeHandle {
    fn rot(&self) -> FixedEdgeHandle {
        FixedEdgeHandle { e: self.e, r: (self.r + 1) % 4 }
    }

    fn rot_inv(&self) -> FixedEdgeHandle {
        FixedEdgeHandle { e: self.e, r: (self.r + 3) % 4 }
    }

    pub fn sym(&self) -> FixedEdgeHandle {
        FixedEdgeHandle { e: self.e, r: (self.r + 2) % 4 }
    }

    pub fn left_face(&self) -> FixedFaceHandle {
        FixedFaceHandle(self.rot_inv())
    }

    pub fn right_face(&self) -> FixedFaceHandle {
        FixedFaceHandle(self.rot())
    }
}

pub struct AllEdgesIterator {
    num_edges: usize,
    cur: usize,
}

impl Iterator for AllEdgesIterator {
    type Item = FixedEdgeHandle;
    
    fn next(&mut self) -> Option<FixedEdgeHandle> {
        if self.cur >= self.num_edges {
            None
        } else {
            self.cur += 1;
            Some(FixedEdgeHandle { e: self.cur - 1, r: 0 })
        }
    }
}

impl AllEdgesIterator {
    fn new(num_edges: usize) -> AllEdgesIterator {
        AllEdgesIterator {
            num_edges: num_edges,
            cur: 0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct FixedVertexHandle(FixedEdgeHandle);

pub struct Subdivision2<V> {
    edges: Vec<[(usize, FixedEdgeHandle); 4]>,
    vertices: Vec<V>,
    num_faces: usize,
}

const UNKOWN_FACE_INDEX: usize = 0;

impl <V> Subdivision2<V> {

    pub fn new() -> Subdivision2<V> {
        Subdivision2 {
            edges: Vec::new(),
            vertices: Vec::new(),
            num_faces: 1,
        }
    }

    pub fn add_new_edge(&mut self, from: V, to: V) -> FixedEdgeHandle {
        let from_index = self.vertices.len();
        self.vertices.push(from);
        let to_index = from_index + 1;
        self.vertices.push(to);
        self.add_new_edge_with_indices(from_index, to_index)
    }

    fn add_new_edge_with_indices(&mut self, from_index: usize, 
                                 to_index: usize) -> FixedEdgeHandle {
        let edge_index = self.edges.len();
        let new_edge = [(from_index, FixedEdgeHandle { e: edge_index, r: 0 }),
                        (UNKOWN_FACE_INDEX, FixedEdgeHandle { e: edge_index, r: 3 }),
                        (to_index, FixedEdgeHandle { e: edge_index, r: 2 }),
                        (UNKOWN_FACE_INDEX, FixedEdgeHandle { e: edge_index, r: 1 })];
        self.edges.push(new_edge);
        FixedEdgeHandle { e: edge_index, r: 0 }
    }

    pub fn edges(&self) -> AllEdgesIterator {
        AllEdgesIterator::new(self.num_edges())
    }

    pub fn num_faces(&self) -> usize {
        self.num_faces
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
    
    pub fn connect_to_new(&mut self, from: &FixedEdgeHandle,
                          to: V) -> FixedEdgeHandle {
        let from_index = self.get_from_index(&from);
        let to_index = self.vertices.len();
        self.vertices.push(to);
        let edge = self.add_new_edge_with_indices(from_index, to_index);
        self.splice(&edge, from);
        edge
    }

    pub fn connect(&mut self, edge_a: &FixedEdgeHandle,
                   edge_b: &FixedEdgeHandle) -> FixedEdgeHandle {
        let to_a = self.get_to_index(edge_a);
        let from_b = self.get_from_index(edge_b);
        let edge = self.add_new_edge_with_indices(to_a, from_b);
        let create_face = {
            // Check if a new face is created
            let until = self.get_from_index(edge_a);
            let mut create_face = false;
            for cur in self.face_edges(&edge_b.left_face()) {
                if self.get_to_index(&cur) == until {
                    create_face = true;
                    break;
                }
            }
            create_face
        };
        let lnext = self.left_face_next(edge_a);
        self.splice(&edge, &lnext);
        self.splice(&edge.sym(), edge_b);
        if create_face {
            let mut cur = edge_b.clone();
            let until = cur;
            loop {
                cur = self.left_face_next(&cur);
                if cur == until {
                    break;
                }
            }
            self.num_faces += 1;
        }
        edge
    }

    pub fn deref_vertex(&self, handle: &FixedVertexHandle) -> &V {
        &self.vertices[self.get_from_index(&handle.0)]
    }

    pub fn eq_faces(&self, lhs: &FixedFaceHandle, rhs: &FixedFaceHandle) -> bool {
        self.get_from_index(&lhs.0) == self.get_from_index(&rhs.0)
    }

    fn eq_vertices(&self, lhs: &FixedVertexHandle, rhs: &FixedVertexHandle) -> bool {
        self.get_from_index(&lhs.0) == self.get_from_index(&rhs.0)
    }

    pub fn face_edges(&self, face: &FixedFaceHandle) -> EdgeIterator<V> 
    {
        EdgeIterator::new(self, face.0.rot(), Box::new(|s, e| s.left_face_next(e)))
    }

    pub fn out_edges(&self, vertex: &FixedVertexHandle) -> EdgeIterator<V> {
            EdgeIterator::new(self, vertex.0, Box::new(|s, e| s.ccw(e)))
    }

    pub fn vertex(&self, index: usize) -> &V {
        &self.vertices[index]
    }

    pub fn from_vertex_handle(&self, edge: &FixedEdgeHandle) -> FixedVertexHandle {
        FixedVertexHandle(edge.clone())
    }

    pub fn to_vertex_handle(&self, edge: &FixedEdgeHandle) -> FixedVertexHandle {
        FixedVertexHandle(edge.sym())
    }

    pub fn from_vertex(&self, edge: &FixedEdgeHandle) -> &V {
        self.deref_vertex(&self.from_vertex_handle(edge))
    }

    pub fn to_vertex(&self, edge: &FixedEdgeHandle) -> &V {
        self.deref_vertex(&self.to_vertex_handle(edge))
    }

    fn get_from_index(&self, edge: &FixedEdgeHandle) -> usize {
        self.edges[edge.e][edge.r].0
    }

    // fn get_from_index_mut(&mut self, edge: &FixedEdgeHandle) -> &mut usize {
    //     &mut self.edges[edge.e][edge.r].0
    // }

    fn get_to_index(&self, edge: &FixedEdgeHandle) -> usize {
        self.get_from_index(&edge.sym())
    }

    pub fn left_face_next(&self, edge: &FixedEdgeHandle) -> FixedEdgeHandle {
        self.ccw(&edge.rot_inv()).rot()
    }

    fn right_face_next(&self, edge: &FixedEdgeHandle) -> FixedEdgeHandle {
        self.ccw(&edge.rot()).rot_inv()
    }

    fn r_prev(&self, edge: &FixedEdgeHandle) -> FixedEdgeHandle {
        self.ccw(&edge.sym())
    }

    pub fn ccw(&self, edge: &FixedEdgeHandle) -> FixedEdgeHandle {
        self.edges[edge.e][edge.r].1
    }

    pub fn cw(&self, edge: &FixedEdgeHandle) -> FixedEdgeHandle {
        self.edges[edge.e][(edge.r + 1) % 4].1.rot()
    }

    fn ccw_mut(&mut self, edge: &FixedEdgeHandle) -> &mut FixedEdgeHandle {
        &mut self.edges[edge.e][edge.r].1
    }

    fn splice(&mut self, a: &FixedEdgeHandle, b: &FixedEdgeHandle) {
        let onext_a = self.ccw(a);
        let onext_b = self.ccw(b);

        let dual_onext_a = self.ccw(a).rot();
        let dual_onext_a_onext = self.ccw(&dual_onext_a);
        let dual_onext_b = self.ccw(b).rot();
        let dual_onext_b_onext = self.ccw(&dual_onext_b);

        *self.ccw_mut(a) = onext_b;
        *self.ccw_mut(b) = onext_a;
        *self.ccw_mut(&dual_onext_a) = dual_onext_b_onext;
        *self.ccw_mut(&dual_onext_b) = dual_onext_a_onext;
    }
}

impl <V> Subdivision2<V> where V: PointObject {
    pub fn connect_to_new_point(&mut self, from: &FixedVertexHandle, to: V) -> FixedEdgeHandle {
        let from_edge = self.insertion_edge(from, &to.position());
        self.connect_to_new(&from_edge, to)
    }

    pub fn connect_points(&mut self, from: &FixedVertexHandle, to: &FixedVertexHandle) -> FixedEdgeHandle {
        let from_vertex = self.deref_vertex(from).position();
        let to_vertex = self.deref_vertex(to).position();
        let from_edge = self.insertion_edge(from, &to_vertex).sym();
        let to_edge = self.insertion_edge(to, &from_vertex);
        self.connect(&from_edge, &to_edge)
    }

    pub fn insertion_edge(&self, vertex: &FixedVertexHandle, point: &Vector2<V::Scalar>) -> FixedEdgeHandle {
        let vertex_position = self.deref_vertex(vertex).position();
        let norm = (point - vertex_position).normalize();
        let point_angle = Angle::atan2(norm.x, norm.y);
        let mut ccw_edge_angle = Rad::new(Float::infinity());
        let mut ccw_edge = vertex.0;
        let mut smallest_angle = Rad::new(Float::infinity());
        let mut smallest_angle_edge = vertex.0;
        for edge in self.out_edges(&vertex) {
            let to = self.deref_vertex(&self.to_vertex_handle(&edge)).position();
            let new_edge_dir = (to - vertex_position).normalize();
            let new_angle = Angle::atan2(new_edge_dir.x, new_edge_dir.y);
            if new_angle < smallest_angle {
                smallest_angle = new_angle;
                smallest_angle_edge = edge;
            }
            if new_angle > point_angle && new_angle < ccw_edge_angle {
                ccw_edge_angle = new_angle;
                ccw_edge = edge;
            }
        }
        if ccw_edge_angle != Rad::new(Float::infinity()) {
            ccw_edge
        } else {
            smallest_angle_edge
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Subdivision2};

    #[test]
    fn test_make_new_edge() {
        let mut s = Subdivision2::<()>::new();
        let edge = s.add_new_edge_with_indices(0, 1);
        assert!(!s.eq_vertices(&s.from_vertex_handle(&edge), 
                               &s.to_vertex_handle(&edge)));
        assert!(s.eq_faces(&edge.left_face(), &edge.right_face()));
        assert!(s.left_face_next(&edge) == edge.sym());
        assert!(s.right_face_next(&edge) == edge.sym());
        assert!(s.ccw(&edge) == edge);
        assert!(s.cw(&edge) == edge);
    }
    
    #[test]
    fn test_connect_to_new() {
        let mut s = Subdivision2::<u32>::new();
        let edge1 = s.add_new_edge_with_indices(0, 1);
        let edge2 = s.connect_to_new(&edge1.sym(), 3);
        
        assert!(s.left_face_next(&edge1) == edge2);
        assert!(s.right_face_next(&edge1) == edge1.sym());
        assert!(s.ccw(&edge1) == edge1);
        assert!(s.ccw(&edge2) == edge1.sym());
        assert!(s.left_face_next(&edge2) == edge2.sym());
        assert!(s.right_face_next(&edge2) == edge1);
    }

    #[test]
    fn test_splice() {
        let mut s = Subdivision2::<usize>::new();
        let edge1 = s.add_new_edge(0, 1);
        let edge2 = s.add_new_edge(1, 2);
        let l_next = s.left_face_next(&edge1);
        s.splice(&l_next, &edge2);
        assert!(s.left_face_next(&edge1) == edge2);
        assert!(s.right_face_next(&edge1) == edge1.sym());
        assert!(s.ccw(&edge1) == edge1);
        assert!(s.ccw(&edge2) == edge1.sym());
        assert!(s.left_face_next(&edge2) == edge2.sym());
        assert!(s.right_face_next(&edge2) == edge1);
    }

    #[test]
    fn test_polygon_adds_face() {
        for size in 3 .. 20 {
            let mut s = Subdivision2::<usize>::new();
            let mut edges = Vec::new();
            let mut last = s.add_new_edge(0, 1);
            edges.push(last.clone());
            for i in 2 .. size {
                last = s.connect_to_new(&last.sym(), i);
                edges.push(last.clone());
            }
            let from = s.to_vertex_handle(&last);
            let to = s.from_vertex_handle(&edges.first().unwrap());
            s.connect(&last, edges.first().unwrap());

            assert_eq!(s.num_faces(), 2);

            let new_face = edges.first().unwrap().left_face();
            let mut num_edges = 0;
            for edge in s.face_edges(&new_face) {
                assert!(s.eq_faces(&new_face, &edge.left_face()));
                num_edges += 1;
            }
            assert_eq!(num_edges, size);
        }
    }

    #[test]
    fn test_two_triangle_properties() {
        /* 
         * Creates 2 triangles: edge1 -> edge2 -> edge3
         * and edge1 -> edge4 -> edge5
         */
        let mut s = Subdivision2::<usize>::new();
        let edge1 = s.add_new_edge(0, 1);
        let edge2 = s.connect_to_new(&edge1.sym(), 2);
        let edge3 = s.connect(&edge2, &edge1);
        let edge4 = s.connect_to_new(&edge2, 3);
        let edge5 = s.connect(&edge4, &edge1);
        assert_eq!(s.num_faces(), 3);
        assert_eq!(s.face_edges(&edge1.left_face()).count(), 3);
        assert_eq!(s.face_edges(&edge1.right_face()).count(), 3);
        assert_eq!(s.face_edges(&edge2.left_face()).count(), 4);

        assert_eq!(s.left_face_next(&edge1), edge4);
        assert_eq!(s.left_face_next(&edge4), edge5);

        assert_eq!(s.right_face_next(&edge1), edge3);
        assert_eq!(s.right_face_next(&edge3), edge2);
        
    }

    #[test]
    fn test_two_triangle_properties2() {
        /* 
         * Creates 2 triangles: edge1 -> edge2 -> edge3
         * and edge1 -> edge4 -> edge5,
         * The second triangle is to the right of the first triangle
         */
        let mut s = Subdivision2::<usize>::new();
        let edge1 = s.add_new_edge(0, 1);
        let edge2 = s.connect_to_new(&edge1.sym(), 2);
        let edge3 = s.connect(&edge2, &edge1);
        let edge4 = s.connect_to_new(&edge1.sym(), 3);
        let edge5 = s.connect(&edge4, &edge3.sym());
        assert_eq!(s.num_faces(), 3);

        let mut cur = edge3.sym();

        assert_eq!(s.face_edges(&edge1.left_face()).count(), 3);
        assert_eq!(s.face_edges(&edge1.right_face()).count(), 3);

        assert_eq!(s.left_face_next(&edge1), edge2);
        assert_eq!(s.left_face_next(&edge2), edge3);

        assert_eq!(s.right_face_next(&edge1), edge5);
        assert_eq!(s.right_face_next(&edge5), edge4);
    }

    #[test]
    fn test_line_doesnt_add_face() {
        // Creating a line should not add a face
        let mut s = Subdivision2::<usize>::new();
        let edge1 = s.add_new_edge(0, 1);
        let edge2 = s.add_new_edge(2, 3);
        s.connect(&edge1, &edge2);
        assert_eq!(s.num_faces(), 1);
    }

    #[test]
    fn test_insert_point_into_triangle() {
        let mut s = Subdivision2::<Vector<f32>>::new();
        let edge1 = s.add_new_edge(Vector2::new(-1., -1.), Vector2::new(1., -1.));
        let edge2 = s.connect_to_new(&edge1.sym(), 2);
        let edge3 = s.connect(&edge2, &edge1);

        let inner1 = s.connect_to_new(&edge1, 3);
        let inner2 = s.connect(&edge2, &inner1.sym());
        let inner3 = s.connect(&edge3, &inner1.sym());

        assert_eq!(s.num_faces(), 4);
        for edge in &[&edge1, &edge2, &edge3] {
            print!("edge.faces")
            assert_eq!(s.face_edges(&edge.left_face()).count(), 3);
            assert_eq!(s.face_edges(&edge.right_face()).count(), 3);            
        }
    }


}
