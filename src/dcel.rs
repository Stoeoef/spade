type FixedVertexHandle = usize;
type HalfEdgeHandle = usize;
type FaceHandle = usize;

#[derive(Debug)]
pub struct VertexEntry<V> {
    data: V,
    out_edge: Option<HalfEdgeHandle>,
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
    next: HalfEdgeHandle,
    prev: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    origin: FixedVertexHandle,
    face: FaceHandle,
}
    

pub struct DCEL<V> {
    vertices: Vec<VertexEntry<V>>,
    edges: Vec<HalfEdgeEntry>,
    num_faces: usize,
}

impl <V> DCEL<V> {

    pub fn new() -> DCEL<V> {
        DCEL {
            vertices: Vec::new(),
            edges: Vec::new(),
            num_faces: 1,
        }
    }

    pub fn insert_vertex(&mut self, vertex: V) -> FixedVertexHandle {
        self.vertices.push(VertexEntry::new(vertex));
        self.vertices.len() - 1
    }

    pub fn connect(&mut self, v0: FixedVertexHandle, v1: FixedVertexHandle, shared_face: FaceHandle) -> HalfEdgeHandle {
        let edge_index = self.edges.len();
        let twin_index = edge_index + 1;
        let edge = HalfEdgeEntry {
            next: 0,
            prev: 0,
            twin: twin_index,
            origin: v0,
            face: shared_face,
        };
        self.edges.push(edge);

        let twin = HalfEdgeEntry {
            next: 0,
            prev: 0,
            twin: edge_index,
            origin: v1,
            face: shared_face,
        };
        self.edges.push(twin);

        self.connect_to_out_edge(edge_index, twin_index, shared_face);
        self.connect_to_out_edge(twin_index, edge_index, shared_face);

        self.vertices[v0].out_edge = Some(edge_index);
        self.vertices[v1].out_edge = Some(twin_index);
        edge_index
    }

    fn connect_to_out_edge(&mut self, edge: HalfEdgeHandle, twin: HalfEdgeHandle, shared_face: FaceHandle) {
        if let Some(out_edge) = self.vertices[self.edges[edge].origin].out_edge {
            // Find an edge adjacent to shared_face
            let mut cur_edge = out_edge;
            loop {
                let cur_twin = self.edges[cur_edge].twin;
                if self.edges[cur_twin].face == shared_face {
                    self.edges[twin].next = self.edges[cur_twin].next;
                    self.edges[cur_twin].next = edge;
                    self.edges[cur_edge].prev = twin;
                    self.edges[edge].prev = self.edges[cur_edge].twin;
                    break;
                }
                cur_edge = self.edges[cur_twin].next;
                if cur_edge == out_edge {
                    panic!("The given shared face is not adjacent to any existing edge of v0");
                }
            }
        } else {
            self.edges[twin].next = edge;
            self.edges[edge].prev = twin;
        }
    }

    pub fn create_face(&mut self, prev_edge_handle: HalfEdgeHandle, next_edge_handle: HalfEdgeHandle) -> HalfEdgeHandle {
        let edge_index = self.edges.len();
        let twin_index = edge_index + 1;
        let new_face = self.num_faces;
        self.num_faces += 1;
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


        // Set the left face of twin to the new face
        let mut cur_edge = edge_index;
        loop {
            self.edges[cur_edge].face = new_face;
            cur_edge = self.edges[cur_edge].next;
            if cur_edge == edge_index {
                break;
            }
        }
        edge_index
    }

    pub fn edge(&self, handle: HalfEdgeHandle) -> EdgeHandle<V> {
        EdgeHandle::new(self, handle)
    }

    pub fn edges(&self) -> EdgesIterator<V> {
        EdgesIterator::new(&self)
    }

    pub fn flip_cw(&mut self, e: HalfEdgeHandle) {
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
    }
}

pub struct EdgesIterator<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    current: HalfEdgeHandle,
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

#[derive(Copy, Clone)]
pub struct EdgeHandle<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    handle: HalfEdgeHandle,
}

pub struct VertexHandle<'a, V> where V: 'a {
    dcel: &'a DCEL<V>,
    handle: FixedVertexHandle,
}

impl <'a, V> ::std::ops::Deref for VertexHandle<'a, V> {
    type Target = V;
    
    fn deref(&self) -> &V {
        &self.dcel.vertices[self.handle].data
    }
}

impl <'a, V> EdgeHandle<'a, V> where V: 'a {

    fn new(dcel: &'a DCEL<V>, handle: HalfEdgeHandle) -> Self {
        EdgeHandle {
            dcel: dcel,
            handle: handle,
        }
    }

    pub fn fix(&self) -> HalfEdgeHandle {
        self.handle
    }

    pub fn from(&self) -> VertexHandle<'a, V> {
        let edge = self.dcel.edges[self.handle];
        VertexHandle::new(self.dcel, edge.origin)
    }

    pub fn to(&self) -> VertexHandle<'a, V> {
        self.sym().from()
    }

    pub fn sym(&self) -> EdgeHandle<'a, V> {
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[self.handle].twin,
        }
    }

    pub fn cw(&self) -> EdgeHandle<'a, V> {
        let twin = self.sym().handle;
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[twin].next,
        }
    }

    pub fn ccw(&self) -> EdgeHandle<'a, V> {
        EdgeHandle {
            dcel: self.dcel,
            handle: self.dcel.edges[self.handle].prev,
        }.sym()
    }
}


impl <'a, V> VertexHandle<'a, V> where V: 'a {
    fn new(dcel: &'a DCEL<V>, handle: FixedVertexHandle) -> VertexHandle<'a, V> {
        VertexHandle {
            dcel: dcel,
            handle: handle,
        }
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
        let e01 = dcel.connect(v0, v1, 0);
        let e12 = dcel.connect(v1, v2, 0);
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

        let e01 = dcel.connect(v0, v1, 0);
        let e12 = dcel.connect(v1, v2, 0);
        let e23 = dcel.connect(v2, v3, 0);
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
        assert_eq!(dcel.edges[e_flip],
                   HalfEdgeEntry {
                       next: e12,
                       prev: e23,
                       twin: dcel.edges[e_flip].twin,
                       origin: 3,
                       face: 2,
                   });
        assert_eq!(dcel.edges[dcel.edges[e_flip].twin],
                   HalfEdgeEntry {
                       next: e30,
                       prev: e01,
                       twin: e_flip,
                       origin: 1,
                       face: 1,
                   });
        
    }

    #[test]
    fn test_cw_ccw() {
        let mut dcel = DCEL::new();
        let v0 = dcel.insert_vertex(());
        let v1 = dcel.insert_vertex(());
        let v2 = dcel.insert_vertex(());
        let v3 = dcel.insert_vertex(());

        let e01 = dcel.connect(v0, v1, 0);
        let e12 = dcel.connect(v1, v2, 0);
        let e23 = dcel.connect(v2, v3, 0);
        let e30 = dcel.create_face(e12, e01);
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

        let e01 = dcel.connect(v[0], v[1], 0);
        dcel.connect(v[1], v[2], 0);
        let e23 = dcel.connect(v[2], v[3], 0);
        let e34 = dcel.connect(v[3], v[4], 0);
        let e40 = dcel.create_face(e34, e01);

        let e02 = dcel.create_face(e40, e23);
        let e03 = dcel.create_face(e40, e34);
        let entry = dcel.edges[e02];
        assert_eq!(entry.next, e23);
        assert_eq!(entry.prev, dcel.edges[e03].twin);
        assert_eq!(entry.origin, v[0]);
    }
}
