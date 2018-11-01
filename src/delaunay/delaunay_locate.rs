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
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Clone)]
#[deprecated(since="1.0.1", note="Replaced by DelaunayWalkLocate")]
#[allow(deprecated)]
#[allow(missing_docs)]
pub struct TriangulationWalkLocate<T: PointN> {
    marker: ::std::marker::PhantomData<T>,
    locate: DelaunayWalkLocate,
}

#[allow(deprecated)]
impl <T: PointN> Default for TriangulationWalkLocate<T> {
    fn default() -> Self {
        TriangulationWalkLocate {
            marker: Default::default(),
            locate: Default::default(),
        }
    }
}

#[allow(deprecated)]
impl <T: PointN> DelaunayLocateStructure<T> for TriangulationWalkLocate<T> {
    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>) {
        self.locate.insert_vertex_entry(entry);
    }

    fn update_vertex_entry(&mut self, new_entry: VertexEntry<T>) {
        self.locate.update_vertex_entry(new_entry);
    }

    fn remove_vertex_entry(&mut self, to_remove: &VertexEntry<T>) {
        self.locate.remove_vertex_entry(to_remove);
    }

    fn find_close_handle(&self, point: &T) -> FixedVertexHandle {
        self.locate.find_close_handle(point)
    }

    fn new_query_result(&self, entry: FixedVertexHandle) {
        DelaunayLocateStructure::<T>::new_query_result(&self.locate, entry);
    }
}

/// Locate strategy that uses an r-tree to locate points in O(log(n)) time.
#[deprecated(since="1.0.1", note="Renamed to DelaunayTreeLocate")]
pub type RTreeDelaunayLocate<V> = RTree<VertexEntry<V>>;

/// Locate strategy that uses an r-tree to locate points in O(log(n)) time.
pub type DelaunayTreeLocate<V> = RTree<VertexEntry<V>>;

/// Locate strategy for Delaunay triangulations.
/// 
/// Many operations of a Delaunay triangulation, like insertion, require to find
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
    fn remove_vertex_entry(&mut self, to_remove: &VertexEntry<T>);
    /// Returns, if possible, a vertex handle that is close to the given point.
    fn find_close_handle(&self, point: &T) -> FixedVertexHandle;
    /// Notifies the locate structure about the result of a query.
    fn new_query_result(&self, entry: FixedVertexHandle);
}

/// An entry of the Delaunay triangulation's internal r-tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
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

impl <V: PointN> VertexEntry<V> {
    pub fn new(point: V, handle: FixedVertexHandle) -> VertexEntry<V> {
        VertexEntry {
            point: point,
            handle: handle,
        }
    }
}

/// Locate strategy that walks along the edges of a triangulation until the target point
/// is found.
///
/// This approach takes O(sqrt(n)) time on average, with a worst case of O(n) and
/// a best case of O(1).
/// This strategy works especially well if subsequent queries like insertion, interpolation
/// or locate queries, are performed close to each other, as the result of the last query
/// operation will be used as hint for the next operation.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
pub struct DelaunayWalkLocate {
    #[cfg_attr(feature = "serde_serialize", serde(skip))]
    last: Arc<AtomicUsize>,
}

impl Default for DelaunayWalkLocate {
    fn default() -> Self {
        DelaunayWalkLocate {
            last: Arc::new(AtomicUsize::new(0)),
        }
    }
}

impl <T: PointN>  DelaunayLocateStructure<T> for DelaunayWalkLocate {

    fn insert_vertex_entry(&mut self, entry: VertexEntry<T>) {
        self.last.store(entry.handle, Ordering::Relaxed);
    }

    fn update_vertex_entry(&mut self, _: VertexEntry<T>) {
    }

    fn remove_vertex_entry(&mut self, _: &VertexEntry<T>) {
    }

     fn find_close_handle(&self, _: &T) -> FixedVertexHandle {
         self.last.load(Ordering::Relaxed)
     }

    fn new_query_result(&self, entry: FixedVertexHandle) {
        self.last.store(entry, Ordering::Relaxed);
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

    fn remove_vertex_entry(&mut self, to_remove: &VertexEntry<T>) {
        assert!(self.lookup_and_remove(&to_remove.point).is_some());
    }

    fn find_close_handle(&self, point: &T) -> FixedVertexHandle {
        self.close_neighbor(point).map(|e| e.handle).unwrap()
    }

    fn new_query_result(&self, _: FixedVertexHandle) {
    }
}
