// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of an n-dimensional r*-tree.

use misc::min_inline;
use std::sync::Arc;
use traits::{SpatialObject};
use ::TwoDimensional;
use point_traits::{PointN, PointNExtensions};
use num::{zero};
use boundingrect::BoundingRect;
use std::iter::Once;

#[cfg(feature = "serde_serialize")]
use serde::{Serialize, Deserialize};

#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
pub struct RTreeOptions {
    max_size: usize,
    min_size: usize,
    reinsertion_count: usize,
}

impl Default for RTreeOptions {
    fn default() -> RTreeOptions {
        RTreeOptions::new()
    }
}

#[doc(hidden)]
impl RTreeOptions {
    pub fn new() -> Self {
        RTreeOptions {
            max_size: 6,
            min_size: 3,
            reinsertion_count: 2,
        }
    }

    pub fn max_size(mut self, max_size: usize) -> Self {
        assert!(max_size > self.min_size, "Max size must be larger than min size");
        self.max_size = max_size;
        self
    }

    pub fn min_size(mut self, min_size: usize) -> Self {
        assert!(self.max_size > min_size, "Min size must be smaller than max size");
        self.min_size = min_size;
        self
    }

    pub fn reinsertion_count(mut self, reinsertion_count: usize) -> Self {
        assert!(0 < reinsertion_count, "Reinsertion cannot be zero");
        assert!(self.max_size > reinsertion_count, "Reinsertion count must be smaller than max size");
        self.reinsertion_count = reinsertion_count;
        self
    }

    pub fn build<T: SpatialObject>(self) -> RTree<T> {
        RTree::new_with_options(self)
    }
}

/// Iterates over all entries in an r-tree.
/// Returned by `RTree::iter()`
pub struct RTreeIterator<'a, T> 
    where T: SpatialObject + 'a {
    data: &'a DirectoryNodeData<T>,
    cur_index: usize, 
    cur_iterator: Option<Box<RTreeNodeIterator<'a, T>>>,
}

#[allow(missing_docs)]
pub enum RTreeNodeIterator<'a, T> 
    where T: SpatialObject + 'a {
    LeafIterator(Once<&'a T>),
    DirectoryNodeIterator(RTreeIterator<'a, T>),
}

impl <'a, T> RTreeIterator<'a, T> 
    where T: SpatialObject {
    fn new(data: &'a DirectoryNodeData<T>) -> RTreeIterator<'a, T> {
        RTreeIterator {
            data: data,
            cur_index: 0,
            cur_iterator: data.children.first().map(
                |child| Box::new(RTreeNodeIterator::new(child))),
        }
    }
}

impl <'a, T> Iterator for RTreeIterator<'a, T>
    where T: SpatialObject {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        if let Some(cur_iterator) = self.cur_iterator.as_mut() {
            if let Some(next) = cur_iterator.next() {
                // Child iterator can still iterate
                Some(next)
            } else {
                loop {
                    // Change to the next child
                    self.cur_index += 1;
                    if let Some(child_node) = self.data.children.get(self.cur_index) {
                        // Set a new iterator...
                        *cur_iterator = Box::new(RTreeNodeIterator::new(child_node));
                        // ... and call it
                        let next = cur_iterator.next();
                        if next.is_some() {
                            return next;
                        }
                    } else {
                        // We've iterated through all of our children
                        return None;
                    }
                }
            }
        } else {
            None
        }
    }
}

impl <'a, T> RTreeNodeIterator<'a, T>
    where T: SpatialObject {

    fn new(node: &'a RTreeNode<T>) -> RTreeNodeIterator<'a, T> {
        use self::RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match node {
            &RTreeNode::Leaf(ref b) => LeafIterator(::std::iter::once(b)),
            &RTreeNode::DirectoryNode(ref data) => 
                DirectoryNodeIterator(RTreeIterator::new(data)),
        }
    }
}

impl <'a, T> Iterator for RTreeNodeIterator<'a, T>
    where T: SpatialObject {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        use self::RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match self {
            &mut LeafIterator(ref mut once) => once.next(),
            &mut DirectoryNodeIterator(ref mut iter) => iter.next(),
        }
    }
}

/// An iterator yielding the elements of an `RTree` ordered by their distance to a query point.
/// 
/// This `struct` is created by the `nearest_neighbor_iter` method on `RTree`
pub struct NearestNeighborIterator<'a, T>
    where T: SpatialObject + 'a
{
    nodes: ::std::collections::binary_heap::BinaryHeap<RTreeNodeDistanceWrapper<'a, T>>,
    query_point: T::Point,
}

struct RTreeNodeDistanceWrapper<'a, T>
    where T: SpatialObject + 'a
{
    node: &'a RTreeNode<T>,
    distance: <T::Point as PointN>::Scalar
}

impl <'a, T> PartialEq for RTreeNodeDistanceWrapper<'a, T>
    where T: SpatialObject + 'a
{
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance
    }
}

impl <'a, T> PartialOrd for RTreeNodeDistanceWrapper<'a, T> 
    where T: SpatialObject + 'a
{
    fn partial_cmp(&self, other: &Self) -> Option<::std::cmp::Ordering> {
        // Inverse comparison creates a min heap
        other.distance.partial_cmp(&self.distance)
    }
}

impl <'a, T> Eq for RTreeNodeDistanceWrapper<'a, T>
    where T: SpatialObject + 'a
{

}

impl <'a, T> Ord for RTreeNodeDistanceWrapper<'a, T> where T: SpatialObject + 'a {
    fn cmp(&self, other: &Self) -> ::std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl <'a, T> NearestNeighborIterator<'a, T>
    where T: SpatialObject + 'a
{
    fn new(root: &'a DirectoryNodeData<T>, query_point: T::Point) -> Self {
        let mut result = NearestNeighborIterator {
            nodes: Default::default(),
            query_point: query_point,
        };
        result.extend_heap(&root.children);
        result
    }

    fn extend_heap(&mut self, children: &'a [RTreeNode<T>]) {
        let query_point = self.query_point.clone();
        self.nodes.extend(children.iter().map(|child| {
            let distance = match child {
                RTreeNode::DirectoryNode(ref data) => data.mbr().min_dist2(&query_point),
                RTreeNode::Leaf(ref t) => t.distance2(&query_point),
            };

            RTreeNodeDistanceWrapper {
                node: child,
                distance: distance
            }
        }));
    }
}

impl <'a, T> Iterator for NearestNeighborIterator<'a, T>
    where T: SpatialObject + 'a
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(current) = self.nodes.pop() {
            match current {
                RTreeNodeDistanceWrapper { node: RTreeNode::DirectoryNode(ref data), .. } => {
                    self.extend_heap(&data.children);
                },
                RTreeNodeDistanceWrapper { node: RTreeNode::Leaf(ref t), ..} => {
                    return Some(t);
                },
            }
        }
        None
    }
}

impl <T> DirectoryNodeData<T>
    where T: SpatialObject + Clone,
          T::Point: TwoDimensional {

    fn bulk_load(options: Arc<RTreeOptions>, mut elements: &mut[T]) -> DirectoryNodeData<T> {
        let m = options.max_size;
        if elements.len() <= m {
            // Reached leaf level
            let elements: Vec<_> = elements.into_iter().map(|e| RTreeNode::Leaf(e.clone())).collect();
            return DirectoryNodeData::new_parent(elements, 1, options);
        }


        let depth = (elements.len() as f32).log(m as f32).ceil() as usize;
        let n_subtree = (m as f32).powi(depth as i32 - 1);
        let remaining_clusters = (elements.len() as f32 / n_subtree).ceil() as usize;

        let num_vertical_slices = (remaining_clusters as f32).sqrt().ceil() as usize;
        let vertical_slice_num_elements = (elements.len() + num_vertical_slices - 1) / num_vertical_slices;
        let mut children = Vec::with_capacity(m + 1);
        create_clusters(&mut elements, vertical_slice_num_elements, 0);
        
        
        let num_clusters_per_slice = (remaining_clusters + num_vertical_slices - 1) / num_vertical_slices;
        for mut slice in elements.chunks_mut(vertical_slice_num_elements) {
            let cluster_num_elements = 
                (slice.len() + num_clusters_per_slice - 1) / num_clusters_per_slice;

            create_clusters(&mut slice, cluster_num_elements, 1);

            for cluster in slice.chunks_mut(cluster_num_elements) {
                let child = DirectoryNodeData::bulk_load(options.clone(), cluster);
                children.push(RTreeNode::DirectoryNode(child));
            }
        }
        DirectoryNodeData::new_parent(children, depth, options)
    }
}

#[doc(hidden)]
impl <T> DirectoryNodeData<T>
    where T: SpatialObject {
    pub fn children(&self) -> &Vec<RTreeNode<T>> {
        &self.children
    }

    pub fn depth(&self) -> usize {
        self.depth
    }

    pub fn mbr(&self) -> BoundingRect<T::Point> {
        self.bounding_box.clone().unwrap()
    }

    fn new(depth: usize, options: Arc<RTreeOptions>) -> DirectoryNodeData<T> {
        DirectoryNodeData {
            bounding_box: None,
            children: Vec::with_capacity(options.max_size + 1),
            options: options,
            depth: depth,
        }
    }

    #[inline]
    fn new_parent(children: Vec<RTreeNode<T>>, depth: usize, options: Arc<RTreeOptions>
                  ) -> DirectoryNodeData<T> {
        let mut result = DirectoryNodeData {
            bounding_box: None,
            children: children,
            depth: depth,
            options: options
        };
        result.update_mbr();
        result
    }

    #[inline]
    fn update_mbr(&mut self) {
        if let Some(first) = self.children.first() {
            let mut new_mbr = first.mbr();
            for child in &self.children[1 .. ] {
                new_mbr.add_rect(&child.mbr());
            }
            self.bounding_box = Some(new_mbr);
        } else {
            self.bounding_box = None;
        }
    }

    #[inline]
    fn update_mbr_with_element(&mut self, element_bb: &BoundingRect<T::Point>) {
        if let Some(ref mut bb) = self.bounding_box {
            bb.add_rect(element_bb);
        }  else {
            self.bounding_box = Some(element_bb.clone());
        }
    }

    fn insert(&mut self, t: RTreeNode<T>, state: &mut InsertionState) -> InsertionResult<T> {
        // Adjust own mbr - the element will most likely become a child of this node
        self.update_mbr_with_element(&t.mbr());
        if t.depth() + 1 == self.depth {
            // Force insertion into this node
            self.add_children(vec![t]);
            return self.resolve_overflow(state);
        }
        let expand = {
            let follow = self.choose_subtree(&t);
            follow.insert(t, state)
        };
        match expand {
            InsertionResult::Split(child) => {
                // Insert into own list
                self.add_children(vec![child]);
                self.resolve_overflow(state)
            },
            result @ InsertionResult::Reinsert(_) => {
                // Reinsertion can shrink the mbr
                self.update_mbr();
                result
            },
            complete => complete,
        }
    }

    fn resolve_overflow(&mut self, state: &mut InsertionState) -> InsertionResult<T> {
        if self.children.len() > self.options.max_size {
            if state.did_reinsert(self.depth) {
                // We did already reinsert on that level - split this node
                let offsplit = self.split();
                InsertionResult::Split(offsplit)
            } else {
                // We didn't attempt to reinsert yet - give it a try
                state.mark_reinsertion(self.depth);
                let reinsertion_nodes = self.reinsert();
                InsertionResult::Reinsert(reinsertion_nodes)
            }
        } else {
            InsertionResult::Complete
        }
    }

    fn split(&mut self) -> RTreeNode<T> {
        let axis = self.get_split_axis();
        assert!(self.children.len() >= 2);
        // Sort along axis
        self.children.sort_by(|l, r| l.mbr().lower().nth(axis).partial_cmp(&r.mbr().lower().nth(axis)).unwrap());
        let mut best = (zero(), zero());
        let mut best_index = self.options.min_size;

        for k in self.options.min_size .. self.children.len() - self.options.min_size + 1 {
            let mut first_mbr = self.children[k - 1].mbr();
            let mut second_mbr = self.children[k].mbr();
            let (l, r) = self.children.split_at(k);
            for child in l {
                first_mbr.add_rect(&child.mbr());
            }
            for child in r {
                second_mbr.add_rect(&child.mbr());
            }

            let overlap_value = first_mbr.intersect(&second_mbr).area();
            let area_value = first_mbr.area() + second_mbr.area();
            let new_best = (overlap_value, area_value);
            if new_best < best || k == self.options.min_size{
                best = new_best;
                best_index = k;
            }
        }
        let offsplit = self.children.split_off(best_index);
        let result = RTreeNode::DirectoryNode(DirectoryNodeData::new_parent(offsplit, self.depth,
                                                                            self.options.clone()));
        self.update_mbr();
        result
    }

    fn reinsert(&mut self) -> Vec<RTreeNode<T>> {
        let center = self.mbr().center();
        // Sort with increasing order so we can use Vec::split_off
        self.children.sort_by(|l, r| {
            let l_center = l.mbr().center();
            let r_center = r.mbr().center();
            l_center.sub(&center).length2().partial_cmp(&(r_center.sub(&center)).length2()).unwrap()
        });
        let num_children = self.children.len();
        let result = self.children.split_off(num_children - self.options.reinsertion_count);
        self.update_mbr();
        result
    }

    fn get_split_axis(&mut self) -> usize {
        let mut best_goodness = zero();
        let mut best_axis = 0;
        for axis in 0 .. T::Point::dimensions() {
            // Sort children along the current axis
            self.children.sort_by(|l, r| l.mbr().lower().nth(axis)
                                  .partial_cmp(&r.mbr().lower().nth(axis)).unwrap());
            for k in self.options.min_size .. self.children.len() - self.options.min_size + 1 {
                let mut first_mbr = self.children[k - 1].mbr();
                let mut second_mbr = self.children[k].mbr();
                let (l, r) = self.children.split_at(k);
                for child in l {
                    first_mbr.add_rect(&child.mbr());
                }
                for child in r {
                    second_mbr.add_rect(&child.mbr());
                }

                let margin_value = first_mbr.half_margin() + second_mbr.half_margin();
                if best_goodness > margin_value || axis == 0 {
                    best_axis = axis;
                    best_goodness = margin_value;
                }
            }
        }
        best_axis
    }

    fn choose_subtree(&mut self, node: &RTreeNode<T>) -> &mut DirectoryNodeData<T> {
        assert!(self.depth >= 2, "Cannot choose subtree on this level");
        let insertion_mbr = node.mbr();
        let mut inclusion_count = 0;
        let mut min_area = zero();
        let mut min_index = 0;
        let mut first = true;
        for (index, child) in self.children.iter().enumerate() {
            let mbr = child.mbr();
            if mbr.contains_rect(&insertion_mbr) {
                inclusion_count += 1;
                let area = mbr.area();
                if area < min_area || first {
                    min_area = area;
                    min_index = index;
                    first = false;
                }
            }
        }
        if inclusion_count == 0 {
            // No inclusion found, subtree depends on overlap and area increase
            let all_leaves = self.depth <= 2;
            let mut min = (zero(), zero(), zero());

            for (index, child1) in self.children.iter().enumerate() {
                let mbr = child1.mbr();
                let mut new_mbr = mbr.clone();
                new_mbr.add_rect(&insertion_mbr);
                let overlap_increase = if all_leaves {
                    // Calculate minimal overlap increase
                    let mut overlap: <T::Point as PointN>::Scalar = zero();
                    let mut new_overlap: <T::Point as PointN>::Scalar = zero();
                    for child2 in self.children.iter() {
                        if child1 as *const _ != child2 as *const _ {
                            let child_mbr = child2.mbr();
                            overlap = overlap.clone() + mbr.intersect(&child_mbr).area();
                            new_overlap = new_overlap.clone() + new_mbr.intersect(&child_mbr).area();
                        }
                    }
                    let overlap_increase = new_overlap - overlap;
                    overlap_increase
                } else {
                    // Don't calculate overlap increase if not all children are leaves
                    zero()
                };
                // Calculate area increase and area
                let area = new_mbr.area();
                let area_increase = area.clone() - mbr.area();
                let new_min = (overlap_increase, area_increase, area);
                if new_min < min || index == 0 {
                    min = new_min;
                    min_index = index;
                }
            }
        }
        if let RTreeNode::DirectoryNode(ref mut data) = self.children[min_index] {
            data
        } else {
            panic!("There must not be leaves on this depth")
        }
    }

    fn add_children(&mut self, mut new_children: Vec<RTreeNode<T>>) {
        if let &mut Some(ref mut bb) = &mut self.bounding_box {
            for child in &new_children {
                bb.add_rect(&child.mbr());
            }
            self.children.append(&mut new_children);
            return;
        } 
        if let Some(first) = new_children.first() {
            let mut bb = first.mbr();
            for child in &new_children[1 ..] {
                bb.add_rect(&child.mbr());
            }
            self.bounding_box = Some(bb);
        }
        self.children.append(&mut new_children);
    }

    fn close_neighbor(&self, point: &T::Point) -> Option<&T> {
        if self.children.is_empty() {
            return None;
        }
        let mut follow = self;
        loop {
            let mut min_min_dist = zero();
            let mut new_follow = follow.children.first().unwrap();
            let mut first = true;
            for child in follow.children.iter() {
                let min_dist = child.mbr().min_dist2(point);
                if min_dist < min_min_dist || first {
                    new_follow = child;
                    min_min_dist = min_dist;
                    first = false;
                }
            }
            match new_follow {
                &RTreeNode::DirectoryNode(ref data) => {
                    follow = data;
                },
                &RTreeNode::Leaf(ref t) => {
                    return Some(t)
                }
            }
        }
    }

    fn extend_heap<'a>(
            heap: &mut ::std::collections::BinaryHeap<RTreeNodeDistanceWrapper<'a, T>>,
            children: &'a [RTreeNode<T>],
            query_point: &T::Point,
            prune_distance_option: &mut Option<<T::Point as PointN>::Scalar>) {
            for child in children {
                let distance = match child {
                    RTreeNode::DirectoryNode(ref data) => data.mbr().min_dist2(query_point),
                    RTreeNode::Leaf(ref t) => t.distance2(query_point)
                };
                let do_prune;
                if let Some(prune_distance) = prune_distance_option.clone() {
                    if distance > prune_distance {
                        do_prune = true;
                    } else {
                        let new_prune_distance = child.mbr().min_max_dist2(query_point);
                        *prune_distance_option = Some(min_inline(prune_distance, new_prune_distance));
                        do_prune = false;
                    }
                } else {
                    *prune_distance_option = Some(child.mbr().min_max_dist2(query_point));
                    do_prune = false;
                }
                if !do_prune {
                    heap.push(RTreeNodeDistanceWrapper {
                        distance: distance,
                        node: child,
                    });
                }
            }
        }

    fn nearest_neighbor(&self, point: &T::Point) -> Option<&T> {
        
        let mut smallest_min_max = None;
        let mut heap = ::std::collections::binary_heap::BinaryHeap::new();
        Self::extend_heap(&mut heap, &self.children, point, &mut smallest_min_max);
        while let Some(current) = heap.pop() {
            match current {
                RTreeNodeDistanceWrapper { node: RTreeNode::DirectoryNode(ref data), .. } => {
                    Self::extend_heap(&mut heap, &data.children, point, &mut smallest_min_max);
                },
                RTreeNodeDistanceWrapper { node: RTreeNode::Leaf(ref t), ..} => {
                    return Some(t);
                },
            }
        }
        None
    }

    fn nearest_neighbors<'a>(&'a self, point: &T::Point,
                         mut nearest_distance: Option<<T::Point as PointN>::Scalar>, 
                             result: &mut Vec<&'a T>) -> Option<<T::Point as PointN>::Scalar> {
        // Calculate smallest minmax-distance
        let mut smallest_min_max: <T::Point as PointN>::Scalar = zero();
        let mut first = true;
        for child in self.children.iter() {
            let new_min = child.mbr().min_max_dist2(point);
            if first {
                first = false;
                smallest_min_max = new_min;
            } else {
                smallest_min_max = min_inline(smallest_min_max, new_min);
            }
        }
        let mut sorted: Vec<_> = self.children.iter().collect();
        sorted.sort_by(|l, r| l.mbr().min_dist2(point).partial_cmp(
            &r.mbr().min_dist2(point)).unwrap());
        for child in sorted.iter() {
            let min_dist = child.mbr().min_dist2(point);
            if min_dist > smallest_min_max || nearest_distance.clone().map(|d| min_dist > d).unwrap_or(false) {
                // Prune this element
                continue;
            }
            match child.nearest_neighbors(point, nearest_distance.clone(), result) {
                Some(nearest) => {
                    nearest_distance = Some(nearest);
                },
                None => {}
            }
        }
        nearest_distance
    }

    fn lookup_and_remove(&mut self, point: &T::Point) -> Option<T> {
        let contains = self.bounding_box.as_ref().map(|bb | bb.contains_point(point)).unwrap_or(false);
        if contains {
            let mut children = ::std::mem::replace(&mut self.children, 
                                                   Vec::new());
            let mut result = None;
            for child in children.drain(..) {
                match child {
                    RTreeNode::DirectoryNode(mut data) => {
                        result = data.lookup_and_remove(point).or(result);
                        if !data.children.is_empty() {
                            // Don't add a node if it has become empty
                            self.children.push(RTreeNode::DirectoryNode(data));
                        }
                    },
                    RTreeNode::Leaf(b) => {
                        if b.contains(point) {
                            result = Some(b);
                        } else {
                            self.children.push(RTreeNode::Leaf(b))
                        }
                    }
                }
            }
            if result.is_some() {
                // Update the mbr if we did remove an element
                self.update_mbr();
            }
            result
        } else {
            None
        }
    }
    
    fn lookup(&self, point: &T::Point) -> Option<&T> {
        let mut todo_list = Vec::with_capacity(40);
        todo_list.push(self);
        while let Some(next) = todo_list.pop() {
            if next.mbr().contains_point(point) {
                for child in next.children.iter() {
                    match child {
                        &RTreeNode::DirectoryNode(ref data) => {
                            todo_list.push(data);
                        },
                        &RTreeNode::Leaf(ref obj) => {
                            if obj.contains(point) {
                                return Some(obj);
                            }
                        },
                    }
                }
            }
        }
        None
    }

    fn lookup_in_circle<'b>(&'b self, result: &mut Vec<&'b T>, origin: &T::Point,
                            radius2: &<T::Point as PointN>::Scalar)
    {
        // Only look at children whose mbr intersects the circle
        for child in self.children.iter().filter(|c| {
            let min_dist2 = c.mbr().min_dist2(origin);
            min_dist2 <= *radius2
        }) {
            match child {
                &RTreeNode::DirectoryNode(ref data) =>
                    data.lookup_in_circle(result, origin, radius2),
                &RTreeNode::Leaf(ref t) => {
                    if t.distance2(origin) < *radius2 {
                        result.push(t);
                    }
                },
            }
        }
    }

    fn lookup_in_rectangle<'b>(&'b self, result: &mut Vec<&'b T>,
                               query_rect: &BoundingRect<T::Point>) {
        for child in self.children.iter().filter(|c| c.mbr().intersects(query_rect)) {
            match child {
                &RTreeNode::DirectoryNode(ref data) => data.lookup_in_rectangle(result, query_rect),
                &RTreeNode::Leaf(ref t) => {
                    if t.mbr().intersects(query_rect) {
                        result.push(t);
                    }
                }
            }
        }
    }
}

impl <T> DirectoryNodeData<T>
    where T: SpatialObject {
    fn lookup_mut(&mut self, point: &T::Point) -> Option<&mut T> {
        let mut todo_list = Vec::with_capacity(40);
        todo_list.push(self);
        while let Some(next) = todo_list.pop() {
            if next.mbr().contains_point(point) {
                for child in next.children.iter_mut() {
                    match child {
                        &mut RTreeNode::DirectoryNode(ref mut data) => {
                            todo_list.push(data);
                        },
                        &mut RTreeNode::Leaf(ref mut obj) => {
                            if (*obj).contains(point) {
                                return Some(obj);
                            }
                        },
                    }
                }
            }
        }
        None
    }
}

#[doc(hidden)]
impl <T> DirectoryNodeData<T>
    where T: SpatialObject + PartialEq {

    pub fn remove(&mut self, to_remove: &T) -> bool {
        let contains = self.bounding_box.as_ref().map(
            |bb| bb.contains_rect(&to_remove.mbr())).unwrap_or(false);
        if contains {
            let mut result = false;
            let mut remove_index = None;
            for (index, child) in self.children.iter_mut().enumerate() {
                match child {
                    &mut RTreeNode::DirectoryNode(ref mut data) => {
                        if data.remove(to_remove) {
                            result = true;
                            if data.children.is_empty() {
                                // Mark this child for removal as it has become empty
                                remove_index = Some(index);
                            }
                            break;
                        }
                    },
                    &mut RTreeNode::Leaf(ref t) => {
                        if t == to_remove {
                            remove_index = Some(index);
                            result = true;
                            break;
                        }
                    }
                }
            }
            if let Some(to_remove) = remove_index {
                self.children.remove(to_remove);
                self.update_mbr();
            }
            result
        } else {
            false
        }
    }

    fn contains(&self, obj: &T) -> bool {
        let contains = self.bounding_box.as_ref().map(
            |bb| bb.contains_rect(&obj.mbr())).unwrap_or(false);
        if contains {
            for child in self.children.iter() {
                match child {
                    &RTreeNode::DirectoryNode(ref data) => {
                        if data.contains(obj) {
                            return true;
                        }
                    },
                    &RTreeNode::Leaf(ref t) => {
                        if t == obj {
                            return true
                        }
                    }
                }
            }
        }
        false
    }
}

enum InsertionResult<T>
    where T: SpatialObject {
    Complete,
    Split(RTreeNode<T>),
    Reinsert(Vec<RTreeNode<T>>),
}

struct InsertionState {
 reinsertions: Vec<bool>,
}

impl InsertionState {
    fn new(max_depth: usize) -> InsertionState {
        let mut reinsertions = Vec::with_capacity(max_depth + 1);
        reinsertions.resize(max_depth, false);
        InsertionState {
            reinsertions: reinsertions,
        }
    }

    fn did_reinsert(&self, depth: usize) -> bool {
        self.reinsertions[depth]
    }

    fn mark_reinsertion(&mut self, depth: usize) {
        self.reinsertions[depth] = true;
    }
}

#[doc(hidden)]
impl <T> RTreeNode<T>
    where T: SpatialObject {
    pub fn depth(&self) -> usize {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.depth,
            _ => 0
        }
    }

    pub fn mbr(&self) -> BoundingRect<T::Point> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.bounding_box.clone().unwrap(),
            &RTreeNode::Leaf(ref t) => t.mbr(),
        }
    }

    fn nearest_neighbors<'a>(&'a self, point: &T::Point, 
                             nearest_distance: Option<<T::Point as PointN>::Scalar>,
                             result: &mut Vec<&'a T>) -> Option<<T::Point as PointN>::Scalar> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.nearest_neighbors(point, nearest_distance, result),
            &RTreeNode::Leaf(ref t) => {
                let distance = t.distance2(point);
                match nearest_distance {
                    Some(nearest) => {                
                        if distance <= nearest {
                            if distance < nearest {
                                // We've found a new minimum element, remove all other neighbors found so far
                                result.clear();
                            }
                            result.push(t);
                            Some(distance)
                        } else {
                            // This object is not among the nearest neigbors
                            None
                        }
                    },
                    None => {
                        result.push(t);
                        Some(distance)
                    }
                }
            }
        }
    }
}

#[doc(hidden)]
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde_serialize", serde(bound(serialize="T: Serialize, T::Point: Serialize", deserialize="T: Deserialize<'de>, T::Point: Deserialize<'de>")))]
pub struct DirectoryNodeData<T>
    where T: SpatialObject {
    bounding_box: Option<BoundingRect<T::Point>>,
    children: Vec<RTreeNode<T>>,
    depth: usize,
    options: Arc<RTreeOptions>,
}

#[doc(hidden)]
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde_serialize", serde(bound(serialize="T: Serialize, T::Point: Serialize", deserialize="T: Deserialize<'de>, T::Point: Deserialize<'de>")))]
pub enum RTreeNode<T>
    where T: SpatialObject {
    Leaf(T),
    DirectoryNode(DirectoryNodeData<T>),
}

/// A rust implementation of n dimensional r*-trees
///
/// [R-trees](https://en.wikipedia.org/wiki/R-tree) provide efficient nearest-neighbor searches for
/// many objects. [R*-trees](https://en.wikipedia.org/wiki/R*_tree) (&quot;R-Star-Trees&quot;) 
/// are a common variant of r-trees and use more advanced heuristics to improve query performance.
/// Instead of linear time complexity, r-trees yield logarithmic complexity for look-up operations
/// and nearest neighbor queries. Inserting into an r-tree runs in O(log(n)) time on average.
/// Also, a bulk loading algorithm is implemented for faster and higher quality tree creation.
/// Some simple geometric primitives that can be inserted into an r-tree can be found in 
/// the `primitives` module. If your object is not among those, consider implementing the
/// `SpatialObject` trait.
/// 
///
/// ```
/// # extern crate nalgebra;
/// # extern crate spade;
///
/// use nalgebra::{Point4};
/// use spade::rtree::RTree;
///
/// # fn main() {
///   let mut tree = RTree::new();
///   tree.insert(Point4::new(13i32, 10, 10, 37));
/// # }
/// ```
/// # Basic Example
///
/// ```
/// extern crate cgmath; // Alternatively: use nalgebra or [f32; 2]
/// extern crate spade;
///
/// use spade::rtree::RTree;
/// use cgmath::Point2;
///
/// fn main() {
/// let mut rtree = RTree::new();
/// // Insert two points
/// rtree.insert(Point2::new(0.5, 0.5f32));
/// rtree.insert(Point2::new(1.0, 1.0f32));
///
/// if rtree.lookup(&Point2::new(0.5, 0.5)).is_some() {
///   println!("Found point at [0.5, 0.5]/");
/// }
/// 
/// let nearest = rtree.nearest_neighbor(&Point2::new(1.5, 1.5)).unwrap();
/// println!("nearest neighbor at [1.5, 1.5]: {:?}", nearest);
///
/// // Iterate over all elements
/// for point in rtree.iter() {
///   println!("Found point: {:?}", point);
/// }
/// }
/// ```
#[derive(Clone)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde_serialize", serde(bound(serialize="T: Serialize, T::Point: Serialize", deserialize="T: Deserialize<'de>, T::Point: Deserialize<'de>")))]
pub struct RTree<T> where T: SpatialObject {
    root: DirectoryNodeData<T>,
    size: usize,
}

impl <T> ::std::fmt::Debug for RTree<T> where T: SpatialObject + ::std::fmt::Debug {

    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> Result<(), ::std::fmt::Error> {
        let mut iter = self.iter();
        write!(f, "RTree {{")?;
        if let Some(next) = iter.next() {
            next.fmt(f)?
        }
        for t in iter {
            write!(f, ", {:?}", t)?;
        }
        write!(f, "}}")
    }
}

impl<T> Default for RTree<T> where T: SpatialObject {
    fn default() -> RTree<T> {
        RTree::new()
    }
}

impl<T> RTree<T> 
    where T: SpatialObject {
    /// Creates an empty r*-tree.
    pub fn new() -> RTree<T> {
        RTree::new_with_options(Default::default())
    }

    /// Returns the trees minimal bounding box.
    pub fn mbr(&self) -> Option<BoundingRect<T::Point>> {
        self.root.bounding_box.clone()
    }

    #[doc(hidden)]
    pub fn new_with_options(options: RTreeOptions) -> RTree<T> {
        let options = Arc::new(options);
        RTree {
            root: DirectoryNodeData::new(1, options),
            size: 0,
        }
    }

    /// Returns the number of elements contained in this r-tree.
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns an iterator over all contained elements.
    pub fn iter(&self) -> RTreeIterator<T> {
        RTreeIterator::new(&self.root)
    }
    
    #[doc(hidden)]
    pub fn root(&self) -> &DirectoryNodeData<T> {
        // This access is only needed for one of the examples
        &self.root
    }

    /// Returns the nearest neighbor.
    ///
    /// Returns `None` if the tree is empty.
    pub fn nearest_neighbor(&self, query_point: &T::Point) -> Option<&T> {
        let result = self.root.nearest_neighbor(query_point);
        if result.is_none() && self.size > 0 {
            // In some cases, the nearest element can be pruned by its own MinMax distance
            // Use the iterator for these cases (which doesn't prune)
            self.nearest_neighbor_iterator(query_point).next()
        } else {
            result
        }
    }

    /// Returns an object close to a given point. This operation is faster than
    /// `nearest_neighbor` but will not neccessarily yield the real nearest neighbor.
    pub fn close_neighbor(&self, point: &T::Point) -> Option<&T> {
        self.root.close_neighbor(point)
    }

    /// Returns the nearest neighbors of a given point.
    ///
    /// All returned values will have the exact same distance from the given query point.
    /// Returns an empty `Vec` if the tree is empty.
    pub fn nearest_neighbors(&self, query_point: &T::Point) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.nearest_neighbors(query_point, None, &mut result);
        }
        result
    }

    /// Returns the nearest n neighbors.
    pub fn nearest_n_neighbors(&self, query_point: &T::Point, n: usize) -> Vec<&T> {
        self.nearest_neighbor_iterator(query_point).take(n).collect()        
    }

   /// Returns an iterator over the nearest neighbors of a point.
    pub fn nearest_neighbor_iterator(&self, query_point: &T::Point) -> NearestNeighborIterator<T> {
        NearestNeighborIterator::new(&self.root, query_point.clone())
    }
    


    /// Returns all objects (partially) contained in a rectangle
    pub fn lookup_in_rectangle(&self, query_rect: &BoundingRect<T::Point>) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.lookup_in_rectangle(&mut result, query_rect);
        }
        result
    }

    /// Returns all objects (partially) contained in a circle.
    ///
    /// Note that `radius2` is the circle's squared radius, not the actual radius.
    /// An object is contained if a part of it lies within the circle.
    pub fn lookup_in_circle(&self, circle_origin: &T::Point, 
                            radius2: &<T::Point as PointN>::Scalar) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.lookup_in_circle(&mut result, circle_origin.into(), radius2);
        }
        result
    }
}

impl <T> RTree<T>
    where T: SpatialObject + Clone,
          T::Point: TwoDimensional {

    /// Creates a new rtree with some initial elements.
    ///
    /// This method should run faster than inserting all elements sequentially.
    /// Also, the resulting rtree should have better quality in terms of
    /// query performance. This is an implementation of the 
    /// [OMT algorithm](http://ftp.informatik.rwth-aachen.de/Publications/CEUR-WS/Vol-74/files/FORUM_18.pdf).
    ///
    /// The current implementation only works for two dimensional data.
    pub fn bulk_load(elements: Vec<T>) -> RTree<T> {
        Self::bulk_load_with_options(Default::default(), elements)
    }

    #[doc(hidden)]
    pub fn bulk_load_with_options(options: RTreeOptions, mut elements: Vec<T>) -> RTree<T> {
        let options = Arc::new(options);
        RTree {
            root: DirectoryNodeData::bulk_load(options, &mut elements),
            size: elements.len(),
        }
    }
}

impl<T> RTree<T> 
    where T: SpatialObject {
    /// Searches for an element at a given position.
    ///
    /// If `query_point` is contained by one object in the tree, this object will be returned.
    /// If multiple objects contain the point, only one of them will be returned.
    pub fn lookup(&self, query_point: &T::Point) -> Option<&T> {
        if self.size > 0 {
            self.root.lookup(query_point)
        } else {
            None
        }
    }

    /// Searches for an element at a given position and returns a mutable
    /// reference.
    /// If `query_point` is contained by multiple objects in the tree,
    /// one of them will be returned.
    /// *Do not change the object's minimal bounding box*.
    pub fn lookup_mut(&mut self, query_point: &T::Point) -> Option<&mut T> {
        if self.size > 0 {
            self.root.lookup_mut(query_point)
        } else {
            None
        }
    }

    /// Inserts a new element into the tree.
    ///
    /// This will require `O(log(n))` operations on average, where n is the number of
    /// elements contained in the tree.
    pub fn insert(&mut self, t: T) {
        let mut state = InsertionState::new(self.root.depth + 1);
        let mut insertion_stack = vec![RTreeNode::Leaf(t)];
        loop {
            if let Some(next) = insertion_stack.pop() {
                match self.root.insert(next, &mut state) {
                    InsertionResult::Split(node) => {
                        // The root node was split, create a new root and increase depth
                        let new_depth = self.root.depth + 1;
                        let options = self.root.options.clone();
                        let old_root = ::std::mem::replace(
                            &mut self.root, DirectoryNodeData::new(
                                new_depth, options));
                        self.root.add_children(vec![RTreeNode::DirectoryNode(old_root), node]);
                    },
                    InsertionResult::Reinsert(nodes) => {
                        // Schedule elements for reinsertion
                        insertion_stack.extend(nodes);
                    },
                    _ => {},
                }
            } else {
                break;
            }
        }
        self.size += 1;
    }

    /// Searches for an element and removes it.
    ///
    /// If the given point is contained by one object in the tree, this object is being removed
    /// and returned. If the point is contained by multiple objects, only one of them is removed and
    /// returned.
    pub fn lookup_and_remove(&mut self, query_point: &T::Point) -> Option<T> {
        if self.size > 0 {
            let result = self.root.lookup_and_remove(query_point);
            if result.is_some() {
                if self.root.children.is_empty() {
                    self.root.depth = 1;
                }
                self.size -= 1;
            }
            result
        } else {
            None
        }
    }
}

impl <T> RTree<T>
    where T: SpatialObject + PartialEq {

    /// Removes an object from the tree.
    ///
    /// Locates and removes an object from the tree, returning
    /// `true` if the element could be removed.
    /// If multiple object's are equal to `to_remove`, only one
    /// will be deleted.
    pub fn remove(&mut self, obj: &T) -> bool {
        if self.size == 0 {
            return false;
        }
        let result = self.root.remove(obj);
        if self.root.children.is_empty() {
            self.root.depth = 1;
        }
        if result {
            self.size -= 1;
        }
        result
    }

    /// Returns `true` if a given object is contained in this tree.
    pub fn contains(&self, obj: &T) -> bool {
        self.root.contains(obj)
    }
}

#[inline]
fn create_clusters<T: SpatialObject>(array: &mut [T], cluster_size: usize, dimension: usize) {
    let comp = |l: &T, r: &T| {
        let l_mbr = l.mbr();
        let r_mbr = r.mbr();
        l_mbr.lower().nth(dimension).partial_cmp(&r_mbr.lower().nth(dimension)).unwrap()
    };

    let mut cur = 0;
    while cur <= array.len() {
        ::pdqselect::select_by(&mut array[cur..], cluster_size, &comp);
        cur += cluster_size;
    }
}

#[cfg(test)]
mod test {
    use super::{RTree};
    use boundingrect::BoundingRect;
    use primitives::{SimpleTriangle, SimpleEdge};
    use cgmath::{Point2, InnerSpace};
    use num::Float;
    use traits::SpatialObject;
    use testutils::*;

    #[test]
    fn test_tree_with_integral_points() {
        // This test should compile
        let mut tree = RTree::new();
        tree.insert(Point2::new(13, 37));
        assert!(tree.lookup(&Point2::new(13, 37)).is_some())
    }

    #[test]
    fn test_tree_with_array_points() {
        // This test should compile
        let mut tree = RTree::<[i32; 3]>::new();
        tree.insert([13i32, 37, 12]);
        assert!(tree.lookup(&[13, 37, 12]).is_some())
    }

    #[test]
    fn test_nearest_neighbor_empty() {
        let tree: RTree<Point2<f32>> = RTree::new();
        assert!(tree.nearest_neighbor(&Point2::new(10., 20f32)).is_none());
    }

    #[test]
    fn test_nearest_neighbor() {
        let (tree, points) = create_random_tree::<f32>(1000, b"nearest-neighbor");
        let sample_points = random_points_with_seed(100, b"aoppllb,1 0 s0,p");
        for sample_point in &sample_points {
            let mut nearest = None;
            let mut closest_dist = Float::infinity();
            for point in &points {
                let new_dist = (point - sample_point).magnitude2();
                if new_dist < closest_dist {
                    closest_dist = new_dist;
                    nearest = Some(point);
                }
            }
            assert!(nearest == tree.nearest_neighbor(sample_point));
        }
    }

    #[test]
    fn test_nearest_neighbor_lines() {
        let points = random_points_with_seed::<f32>(200, b"ThirstyChicken?=");
        let lines: Vec<_> = points.chunks(2)
                                .map(|xs| SimpleEdge::new(xs[0], xs[1]))
                                .collect();
        let tree = RTree::bulk_load(lines.clone());
        let sample_points = random_points_with_seed(1000, b"S)nug( life .rat");
        for sample_point in &sample_points {
            let mut nearest = None;
            let mut closest_dist = Float::infinity();
            for line in &lines {
                let new_dist = line.distance2(sample_point);
                if new_dist < closest_dist {
                    closest_dist = new_dist;
                    nearest = Some(line);
                }
            }
            assert_eq!(nearest, tree.nearest_neighbor(sample_point));
        }
    }

    #[test]
    fn test_nearest_neighbor_iterator() {

        let (tree, mut points) = create_random_tree::<f32>(100, b"switzerlandchees");
        let sample_points = random_points_with_seed(10, b"15yummydelicious");
        for sample_point in &sample_points {
            points.sort_by(|l, r| l.distance2(sample_point).partial_cmp(&r.distance2(sample_point)).unwrap());
            let collected: Vec<_> = tree.nearest_neighbor_iterator(sample_point).cloned().collect();
            assert_eq!(points, collected);
        }
    }

    #[test]
    fn test_nearest_neighbor_iterator_empty() {
        let tree: RTree<[f32; 2]> = RTree::new();
        assert_eq!(tree.nearest_neighbor_iterator(&[0.0, 0.0]).next(), None);
    }

    #[test]
    fn test_lookup_in_circle() {
        let (tree, points) = create_random_tree::<f32>(1000, b"never turn your ");
        let sample_points = random_points_with_seed(100, b"back on friends~");
        const RADIUS: f32 = 20.;
        for sample_point in &sample_points {
            let mut expected = Vec::new();
            for point in &points {
                let new_dist = (point - sample_point).magnitude2();
                if new_dist < RADIUS * RADIUS {
                    expected.push(point);
                }
            }
            let points = tree.lookup_in_circle(sample_point, &(RADIUS * RADIUS));
            assert_eq!(points.len(), expected.len());
            for p in &points {
                assert!(expected.contains(p));
            }
            for p in &expected {
                assert!(points.contains(p));
            }
        }
    }

    #[test]
    fn test_lookup_in_rect() {
        use cgmath::{EuclideanSpace, Vector2};

        let (tree, points) = create_random_tree::<f32>(1000, b"the heir of the ");
        let sample_points = random_points_with_seed(100, b"highlord you bet");
        const SIZE: f32 = 20.;
        for sample_point in &sample_points {
            let sample_rect = BoundingRect::from_corners(
                sample_point, &Point2::from_vec(sample_point.to_vec() + Vector2::new(SIZE, SIZE)));
            let mut expected = Vec::new();
            for point in &points {
                if sample_rect.contains_point(point) {
                    expected.push(point);
                }
            }
            let points = tree.lookup_in_rectangle(&sample_rect);
            assert_eq!(points.len(), expected.len());
            for p in &points {
                assert!(expected.contains(p));
            }
            for p in &expected {
                assert!(points.contains(p));
            }
        }
    }

    #[test]
    fn test_nearest_neighbors() {
        let mut tree = RTree::new();
        assert!(tree.nearest_neighbors(&Point2::new(1, 0)).is_empty());
        tree.insert(Point2::new(1, 0));
        tree.insert(Point2::new(0, 1));
        tree.insert(Point2::new(-1, 0));
        tree.insert(Point2::new(0, -1));
        tree.insert(Point2::new(3, 0));
        tree.insert(Point2::new(2, 1));
        tree.insert(Point2::new(2, -1));
        assert_eq!(tree.nearest_neighbors(&Point2::new(0, 0)).len(), 4);
        assert_eq!(tree.nearest_neighbors(&Point2::new(1, 0)).len(), 1);
        assert_eq!(tree.nearest_neighbors(&Point2::new(2, 0)).len(), 4);
    }

    #[test]
    fn test_lookup() {
        let (mut tree, mut points) = create_random_tree::<f32>(10000, b"The enemy of min");
        let sample_points = random_points_with_seed(1000, b"Isn't he of your");
        for sample_point in &sample_points {
            assert!(tree.lookup(sample_point).is_none());
            assert!(tree.lookup_mut(sample_point).is_none());
        }
        for point in points.iter_mut() {
            assert!(tree.lookup(point) == Some(point));
            assert!(tree.lookup_mut(point) == Some(point));
        }
    }

    #[test]
    fn test_lookup_and_remove() {
        let (mut tree, points) = create_random_tree::<f32>(10000, b"kind? And finall");
        let sample_points = random_points_with_seed(1000, b"you may follow m");
        for sample_point in &sample_points {
            assert!(!tree.lookup_and_remove(sample_point).is_some());
        }

        // Test if all points are still there
        for point in &points {
            assert_eq!(tree.lookup(point), Some(point));
        }
        // Now remove all points
        for point in &points {
            assert_eq!(tree.lookup_and_remove(point).as_ref(), Some(point));
        }
        assert!(tree.root.children.is_empty());
        tree.insert(Point2::new(20., 10.));
    }

    #[test]
    fn test_remove() {
        let random_points = random_points_with_seed(300, b"Farewell, he sai");
        let mut triangles = Vec::with_capacity(100);
        for ps in random_points.chunks(3) {
            triangles.push(SimpleTriangle::new(ps[0], ps[1], ps[2]));
        }

        for ps in random_points.chunks(3) {
            // Insert every triangle twice
            triangles.push(SimpleTriangle::new(ps[2], ps[1], ps[0]));
        }

        let mut tree = RTree::new();
        for _ in 0 .. 2 {
            for triangle in triangles.iter().cloned() {
                tree.insert(triangle);
            }

            // Try to remove a triangle that is not contained
            let triangle = SimpleTriangle::new(Point2::new(0.0, 0.0), 
                                               Point2::new(1.0, 0.0), 
                                               Point2::new(1.0, 1.0));
            assert!(!tree.remove(&triangle));
            let mut size = 200usize;
            for triangle in &triangles {
                assert!(tree.remove(triangle));
                size -= 1;
                assert_eq!(tree.size(), size);
            }
        }
    }

    #[test]
    fn test_remove_line() {
        let mut tree = RTree::new();
        let edge = SimpleEdge::new([0f32, 0.], [1., 1.]);
        tree.insert(edge.clone());
        tree.insert(SimpleEdge::new([3., 4.], [0., 2.]));
        tree.insert(SimpleEdge::new([-3., -4.], [0., 2.]));
        tree.remove(&edge);
        assert_eq!(tree.size(), 2);
    }


    #[test]
    fn test_iteration() {
        let (tree, reference_points) = create_random_tree::<f32>(100, b"Nightfall~quietl");
        // Check if the set of reference points and the points given by
        // iteration are equal
        assert_eq!(tree.iter().count(), 100);
        let points: Vec<_> = tree.iter().map(|v| v.clone()).collect();
        for p in points.iter() {
            assert!(reference_points.contains(p));
        }
        for p in reference_points.iter() {
            assert!(points.contains(p));
        }
    }
    
    #[test]
    fn test_higher_dimensions() {
        use nalgebra::Point4;
        use rand::{XorShiftRng, SeedableRng};
        use rand::distributions::{Standard, Distribution};

        let mut tree: RTree<Point4<f32>> = RTree::new();
        let mut rng = XorShiftRng::from_seed(b"crept in and cha".clone());
        let mut entries = Vec::new();
        for _ in 0 .. 500 {
            let (x, y, z, w) = Standard.sample(&mut rng);
            let entry = Point4::new(x, y, z, w);
            entries.push(entry);
            tree.insert(entry);
        }

        for entry in &entries {
            assert!(tree.lookup(entry).is_some());
            assert_eq!(tree.nearest_neighbor(entry), Some(entry))
        }
    }

    #[test]
    fn test_bulk_load() {
        const MAX_POINTS: usize = 200;
        let points = random_points_with_seed::<f32>(
            MAX_POINTS, b"Immortal lands l");

        // Check if no number of points crashes
        let mut tree = RTree::new();
        for num in 0 .. MAX_POINTS {
            tree = RTree::bulk_load(points[..num + 1].to_vec());
        }
        assert_eq!(tree.iter().count(), MAX_POINTS);
        // Check if all points have been inserted
        for p in &points {
            assert!(tree.lookup(p).is_some());
        }
    }

    #[cfg(feature="serde")]
    #[test]
    fn test_serialization() {
        use serde_json;
        const SIZE: usize = 1000;
        let points = random_points_with_seed::<f32>(SIZE, b"No sign of life ");
        let tree = RTree::bulk_load(points.clone());
        let json = serde_json::to_string(&tree).expect("Serializing tree failed");
        let parsed: RTree<Point2<f32>> = serde_json::from_str(&json).expect("Deserializing tree failed");
        assert_eq!(parsed.size(), SIZE);
        for point in &points {
            assert!(parsed.contains(point));
        }
    }
}
