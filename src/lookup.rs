use traits::{HasPosition, SpatialObject};
use vector_traits::VectorN;
use rtree::RTree;
use dcel::FixedVertexHandle;

pub type RTreeDelaunayLookup<V> = RTree<VertexEntry<V>>;

pub trait SpatialLookupStructure<T: SpatialObject> : LookupStructure<T> {
    fn nearest_neighbor(&self, point: &T::Vector) -> Option<&T>;
}

pub trait LookupStructure<T: SpatialObject> {
    
    fn insert(&mut self, object: T);
    fn lookup_and_remove(&mut self, query_point: &T::Vector) -> Option<T>;
    fn lookup(&self, point: &T::Vector) -> Option<&T>;
    fn lookup_mut(&mut self, point: &T::Vector) -> Option<&mut T>;
}

pub trait DelaunayLookupStructure<T: VectorN> : Default + Clone {

    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>);
    fn update_vertex_entry(&mut self, new_entry: VertexEntry<T>);
    fn remove_vertex_entry(&mut self, to_remove: &T);

    fn find_close_handle(&self, point: &T) -> Option<FixedVertexHandle>;
}

/// An entry of the delaunay triangulation's internal r-tree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VertexEntry<V> where V: VectorN {
    pub point: V,
    pub handle: FixedVertexHandle,
}

impl<V> HasPosition for VertexEntry<V> 
    where V: VectorN {
    type Vector = V;
    fn position(&self) -> V {
        self.point.clone()
    }
}

#[derive(Clone)]
pub struct TriangulationWalkLookup<T: VectorN> {
    last: Option<VertexEntry<T>>,
}

impl <T: VectorN> Default for TriangulationWalkLookup<T> {
    fn default() -> Self {
        TriangulationWalkLookup {
            last: None,
        }
    }
}

impl <T: VectorN> DelaunayLookupStructure<T> for TriangulationWalkLookup<T> {

    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>) {
        self.last = Some(entry);
    }

    fn update_vertex_entry(&mut self, new_entry: VertexEntry<T>) {
        if self.last.as_ref().map(|e| &e.point) == Some(&new_entry.point) {
            self.last = Some(new_entry);
        }
    }

    fn remove_vertex_entry(&mut self, to_remove: &T) {
        if self.last.as_ref().map(|e| &e.point) == Some(to_remove) {
            self.last = None;
        }
    }

     fn find_close_handle(&self, _: &T) -> Option<FixedVertexHandle> {
         self.last.as_ref().map(|e| e.handle)
     }
}

impl <T: VectorN> DelaunayLookupStructure<T> for RTree<VertexEntry<T>> {
    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>) {
        self.insert(entry);
    }

    fn update_vertex_entry(&mut self, new_entry: VertexEntry<T>) {
        let old_pos = new_entry.point.clone();
        *self.lookup_mut(&old_pos).unwrap() = new_entry;
    }

    fn remove_vertex_entry(&mut self, to_remove: &T) {
        assert!(self.lookup_and_remove(to_remove).is_some());
    }

    fn find_close_handle(&self, point: &T) -> Option<FixedVertexHandle> {
        self.close_neighbor(point).map(|e| e.handle)
    }
}
