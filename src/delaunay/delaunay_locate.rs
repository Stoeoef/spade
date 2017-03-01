// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use traits::HasPosition;
use point_traits::PointN;
use rtree::RTree;
use delaunay::FixedVertexHandle;

/// Locate strategy that uses an r-tree to locate points in O(log(n)) time.
pub type RTreeDelaunayLocate<V> = RTree<VertexEntry<V>>;

/// Locate strategy for delaunay triangulations.
/// 
/// Many operations of a delaunay triangulation, like insertion, require to find
/// the triangle that contains a certain point. For larger triangulations, this step
/// tends to take most time, since it might be necessary to traverse large parts of
/// the triangulation until the desired triangle is found. To mitigate this, spade
/// offers different methods on how point are located. Each method implements this 
/// trait. It is recommended to choose from one of the given implementations.
pub trait DelaunayLocateStructure<T: PointN> : Default + Clone {
    
    /// This method is called when a new vertex entry has been inserted.
    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>);
    /// This method is called when a vertex has been updated.
    fn update_vertex_entry(&mut self, new_entry: VertexEntry<T>);
    /// This method is callend when a vertex has been removed.
    fn remove_vertex_entry(&mut self, to_remove: &T);
    /// Returns, if possible, a vertex handle that is close to the given point.
    fn find_close_handle(&self, point: &T) -> Option<FixedVertexHandle>;
}

/// An entry of the delaunay triangulation's internal r-tree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VertexEntry<V> where V: PointN {
    pub point: V,
    pub handle: FixedVertexHandle,
}

impl<V> HasPosition for VertexEntry<V> 
    where V: PointN {
    type Point = V;
    fn position(&self) -> V {
        self.point.clone()
    }
}

/// Locate strategy that walks along the edges of a triangulation until the target point
/// is found. This approach takes O(sqrt(n)) time on average, with a worst case of O(n) and
/// a best case of O(1).
#[derive(Clone)]
pub struct TriangulationWalkLocate<T: PointN> {
    last: Option<VertexEntry<T>>,
}

impl <T: PointN> Default for TriangulationWalkLocate<T> {
    fn default() -> Self {
        TriangulationWalkLocate {
            last: None,
        }
    }
}

impl <T: PointN> DelaunayLocateStructure<T> for TriangulationWalkLocate<T> {

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

impl <T: PointN> DelaunayLocateStructure<T> for RTree<VertexEntry<T>> {
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
