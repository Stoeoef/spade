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

use cgmath::{Vector2, InnerSpace, PartialOrd};
use cgmath::conv::array2;
use std::sync::Arc;
use traits::{SpatialObject};
use num::{Bounded, zero};
use boundingvolume::BoundingRect;
use std::iter::Once;
use misc::length2;

#[doc(hidden)]
#[derive(Eq, PartialEq, Clone, Debug)]
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
    pub fn new() -> RTreeOptions {
        RTreeOptions {
            max_size: 6,
            min_size: 3,
            reinsertion_count: 2,
        }
    }

    pub fn set_max_size(mut self, max_size: usize) -> RTreeOptions {
        assert!(max_size > self.min_size);
        self.max_size = max_size;
        self
    }

    pub fn set_min_size(mut self, min_size: usize) -> RTreeOptions {
        assert!(self.max_size > min_size);
        self.min_size = min_size;
        self
    }

    pub fn set_reinsertion_count(mut self, reinsertion_count: usize) -> RTreeOptions {
        assert!(0 < reinsertion_count && self.max_size > reinsertion_count);
        self.reinsertion_count = reinsertion_count;
        self
    }

    pub fn build<T: SpatialObject>(self) -> RTree<T> {
        RTree::new_with_options(self)
    }
}


pub struct RTreeIterator<'a, T> where T: SpatialObject + 'a {
    data: &'a DirectoryNodeData<T>,
    cur_index: usize, 
    cur_iterator: Option<Box<RTreeNodeIterator<'a, T>>>,
}

#[doc(hidden)]
pub enum RTreeNodeIterator<'a, T> where T: SpatialObject + 'a {
    LeafIterator(Once<&'a T>),
    DirectoryNodeIterator(RTreeIterator<'a, T>),
}

impl <'a, T> RTreeIterator<'a, T> where T: SpatialObject {
    fn new(data: &'a DirectoryNodeData<T>) -> RTreeIterator<'a, T> {
        RTreeIterator {
            data: data,
            cur_index: 0,
            cur_iterator: data.children.first().map(
                |child| Box::new(RTreeNodeIterator::new(child))),
        }
    }
}

impl <'a, T> Iterator for RTreeIterator<'a, T> where T: SpatialObject {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        if let Some(mut cur_iterator) = self.cur_iterator.as_mut() {
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

impl <'a, T> RTreeNodeIterator<'a, T> where T: SpatialObject {

    fn new(node: &'a RTreeNode<T>) -> RTreeNodeIterator<'a, T> {
        use RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match node {
            &RTreeNode::Leaf(ref t) => LeafIterator(::std::iter::once(t)),
            &RTreeNode::DirectoryNode(ref data) => 
                DirectoryNodeIterator(RTreeIterator::new(data)),
        }
    }
}

impl <'a, T> Iterator for RTreeNodeIterator<'a, T> where T: SpatialObject {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        use RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match self {
            &mut LeafIterator(ref mut once) => once.next(),
            &mut DirectoryNodeIterator(ref mut iter) => iter.next(),
        }
    }
}

impl <T> DirectoryNodeData<T> where T: SpatialObject {

    pub fn children(&self) -> &Vec<RTreeNode<T>> {
        &self.children
    }

    pub fn depth(&self) -> usize {
        self.depth
    }

    pub fn mbr(&self) -> BoundingRect<T::Scalar> {
        self.bounding_box.clone()
    }

    fn new(depth: usize, options: Arc<RTreeOptions>) -> DirectoryNodeData<T> {
        DirectoryNodeData {
            bounding_box: BoundingRect::new(),
            children: Box::new(Vec::with_capacity(options.max_size + 1)),
            options: options,
            depth: depth,
        }
    }

    fn new_parent(mut children: Box<Vec<RTreeNode<T>>>, depth: usize, options: Arc<RTreeOptions>
                  ) -> DirectoryNodeData<T> {
        let missing = options.max_size + 1 - children.len();
        children.reserve_exact(missing);
        let mut result = DirectoryNodeData {
            bounding_box: BoundingRect::new(),
            children: children,
            depth: depth,
            options: options
        };
        result.update_mbr();
        result
    }

    #[inline]
    fn update_mbr(&mut self) {
        self.bounding_box = BoundingRect::new();
        for child in self.children.iter() {
            self.bounding_box.add_rect(&child.mbr());
        }
    }

    #[inline]
    fn update_mbr_with_element(&mut self, element_bb: &BoundingRect<T::Scalar>) {
        self.bounding_box.add_rect(element_bb);
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
            let mut follow = self.choose_subtree(&t);
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
        self.children.sort_by(|l, r| l.mbr().lower()[axis].partial_cmp(&r.mbr().lower()[axis]).unwrap());
        let mut best = (Bounded::max_value(), Bounded::max_value());
        let mut best_index = self.options.min_size;

        for k in self.options.min_size .. self.children.len() - self.options.min_size + 1 {
            let mut first_mbr = BoundingRect::new();
            let mut second_mbr = BoundingRect::new();
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
            if new_best < best {
                best = new_best;
                best_index = k;
            }
        }
        let offsplit = Box::new(self.children.split_off(best_index));
        let result = RTreeNode::DirectoryNode(DirectoryNodeData::new_parent(offsplit, self.depth,
                                                                            self.options.clone()));
        self.update_mbr();
        result
    }

    fn reinsert(&mut self) -> Vec<RTreeNode<T>> {
        let center = self.bounding_box.center();
        // Sort with increasing order so we can use Vec::split_off
        self.children.sort_by(|l, r| {
            let l_center = l.mbr().center();
            let r_center = r.mbr().center();
            length2(&(l_center - center)).partial_cmp(
                &length2(&(r_center - center))).unwrap()
        });
        let num_children = self.children.len();
        let result = self.children.split_off(num_children - self.options.reinsertion_count);
        self.update_mbr();
        result
    }

    fn get_split_axis(&mut self) -> usize {
        let mut best_goodness = Bounded::max_value();
        let mut best_axis = 0;
        for axis in 0 .. 2usize {
            // Sort children along the current axis
            self.children.sort_by(|l, r| l.mbr().lower()[axis]
                                  .partial_cmp(&r.mbr().lower()[axis]).unwrap());
            for k in self.options.min_size .. self.children.len() - self.options.min_size + 1 {
                let mut first_mbr = BoundingRect::new();
                let mut second_mbr = BoundingRect::new();
                let (l, r) = self.children.split_at(k);
                for child in l {
                    first_mbr.add_rect(&child.mbr());
                }
                for child in r {
                    second_mbr.add_rect(&child.mbr());
                }

                let margin_value = first_mbr.half_margin() + second_mbr.half_margin();
                if best_goodness > margin_value {
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
        let mut min_area = Bounded::max_value();
        let mut min_index = 0;
        for (index, child) in self.children.iter().enumerate() {
            let mbr = child.mbr();
            if mbr.contains_rect(&insertion_mbr) {
                inclusion_count += 1;
                let area = mbr.area();
                if area < min_area {
                    min_area = area;
                    min_index = index;
                }
            }
        }
        if inclusion_count == 0 {
            // No inclusion found, subtree depends on overlap and area increase
            let all_leaves = self.depth <= 2;

            let mut min = (Bounded::max_value(), 
                           Bounded::max_value(), Bounded::max_value());

            for (index, child1) in self.children.iter().enumerate() {
                let mbr = child1.mbr();
                let mut new_mbr = mbr.clone();
                new_mbr.add_rect(&insertion_mbr);
                let overlap_increase = if all_leaves {
                    // Calculate minimal overlap increase
                    let mut overlap: T::Scalar = zero();
                    let mut new_overlap: T::Scalar = zero();
                    for child2 in self.children.iter() {
                        if child1 as *const RTreeNode<T> != child2 as *const RTreeNode<T> {
                            let child_mbr = child2.mbr();
                            overlap = overlap + mbr.intersect(&child_mbr).area();
                            new_overlap = new_overlap +  new_mbr.intersect(&child_mbr).area();
                        }
                    }
                    new_overlap - overlap
                } else {
                    // Don't calculate overlap increase if not all children are leaves
                    zero()
                };
                // Calculate area increase and area
                let area_increase = new_mbr.area() - mbr.area();
                let area = new_mbr.area();
                let new_min = (overlap_increase, area_increase, area);
                if new_min <= min {
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
        self.children.append(&mut new_children);
    }

    fn nearest_neighbor(&self, point: &Vector2<T::Scalar>, 
                        mut nearest_distance: T::Scalar) -> Option<&T> {
        let mut nearest = None;
        // Calculate smallest minmax-distance
        let mut smallest_min_max: T::Scalar = Bounded::max_value();
        for child in self.children.iter() {
            smallest_min_max = smallest_min_max.partial_min(child.mbr().min_max_dist2(point));
        }
        let mut sorted: Vec<_> = self.children.iter().collect();
        sorted.sort_by(|l, r| l.mbr().min_dist2(point).partial_cmp(
            &r.mbr().min_dist2(point)).unwrap());
        for child in sorted.iter() {
            let min_dist = child.mbr().min_dist2(point);
            if min_dist > smallest_min_max || min_dist > nearest_distance {
                // Prune this element
                continue;
            }
            match child.nearest_neighbor(point, nearest_distance) {
                Some(t) => {
                    nearest_distance = t.distance(array2(point.clone())).into();
                    nearest = Some(t);
                },
                None => {}
            }
        }
        nearest
    }

    fn nearest_n_neighbors<'a>(&'a self, point: &Vector2<T::Scalar>, n: usize, result: &mut Vec<&'a T>) {

        // let mut sorted: Vec<_> = self.children.iter().collect();
        // sorted.sort_by(|l, r| l.mbr().min_dist2(point).partial_cmp(
        //     &r.mbr().min_dist2(point)).unwrap());
        for child in self.children.iter() {
        // for child in sorted.iter().cloned() {
            let min_dist = child.mbr().min_dist2(point);
            if result.len() == n && min_dist >= result.last().unwrap().distance(array2(point.clone())) {
                // Prune this element
                continue;
            }
            match child {
                &RTreeNode::DirectoryNode(ref data) => {
                    data.nearest_n_neighbors(point, n, result);
                },
                &RTreeNode::Leaf(ref t) => {
                    let distance = t.distance(array2(point.clone()));
                    if result.len() != n || distance < result.last().unwrap().distance(array2(point.clone())) {
                        if result.len() == n {
                            result.pop();
                        }
                        let index = match result.binary_search_by(|e| e.distance(array2(point.clone())).partial_cmp(
                            &distance).unwrap()) {
                            Ok(index) => index,
                            Err(index) => index,
                        };
                        result.insert(index, t);
                    }
                }
            }
        }
    }

    fn lookup_and_remove(&mut self, point: &Vector2<T::Scalar>) -> Option<T> {
        let contains = self.bounding_box.contains_point(&array2(*point));
        if contains {
            let mut children = ::std::mem::replace(&mut self.children, 
                                                   Box::new(Vec::new()));
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
                    RTreeNode::Leaf(t) => {
                        if t.contains(array2(point.clone())) {
                            result = Some(t);
                        } else {
                            self.children.push(RTreeNode::Leaf(t))
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
    
    fn lookup(&self, point: &Vector2<T::Scalar>) -> Option<&T> {
        if self.bounding_box.contains_point(&array2(*point)) {
            for child in self.children.iter() {
                match child {
                    &RTreeNode::DirectoryNode(ref data) => {
                        let result = data.lookup(point);
                        if result.is_some() {
                            return result;
                        }
                    },
                    &RTreeNode::Leaf(ref t) => {
                        if t.contains(array2(point.clone())) {
                            return Some(t);
                        }
                    }
                }
            }
        }
        None
    }

    fn lookup_in_circle<'b>(&'b self, result: &mut Vec<&'b T>, origin: &Vector2<T::Scalar>,
                        radius2: &T::Scalar)
    {
        // Only look at children whose mbr intersects the circle
        for child in self.children.iter().filter(|c| {
            let min_dist2 = c.mbr().min_dist2(origin);
            min_dist2 < *radius2
        }) {
            match child {
                &RTreeNode::DirectoryNode(ref data) =>
                    data.lookup_in_circle(result, origin, radius2),
                &RTreeNode::Leaf(ref t) => {
                    if t.distance(array2(origin.clone())) < *radius2 {
                        result.push(t);
                    }
                },
            }
        }
    }
}

impl <T> DirectoryNodeData<T> where T: SpatialObject + PartialEq {

    pub fn remove(&mut self, to_remove: &T) -> bool {
        let contains = self.bounding_box.contains_rect(&to_remove.mbr());
        if contains {
            let mut children = ::std::mem::replace(&mut self.children, 
                                                   Box::new(Vec::new()));
            let mut result = false;
            for child in children.drain(..) {
                match child {
                    RTreeNode::DirectoryNode(mut data) => {
                        result = data.remove(to_remove) | result;
                        if !data.children.is_empty() {
                            // Don't add a node if it has become empty
                            self.children.push(RTreeNode::DirectoryNode(data));
                        }
                    },
                    RTreeNode::Leaf(t) => {
                        if t == *to_remove {
                            result = true;
                        } else {
                            self.children.push(RTreeNode::Leaf(t))
                        }
                    }
                }
            }
            if result {
                // Update the mbr if we did remove an element
                self.update_mbr();
            }
            result
        } else {
            false
        }
    }
}

enum InsertionResult<T> where T: SpatialObject {
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

impl <T> RTreeNode<T> where T: SpatialObject {

    pub fn depth(&self) -> usize {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.depth,
            _ => 0
        }
    }

    pub fn mbr(&self) -> BoundingRect<T::Scalar> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.bounding_box.clone(),
            &RTreeNode::Leaf(ref t) => t.mbr(),
        }
    }

    fn nearest_neighbor(&self, point: &Vector2<T::Scalar>, nearest_distance: T::Scalar) 
                        -> Option<&T> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.nearest_neighbor(point, nearest_distance),
            &RTreeNode::Leaf(ref t) => {
                let distance = t.distance(array2(point.clone()));
                if distance < nearest_distance {
                    Some(t)
                } else {
                    None
                }
            }
        }
    }
}

#[doc(hidden)]
#[derive(Clone)]
pub struct DirectoryNodeData<T> where T: SpatialObject {
    bounding_box: BoundingRect<T::Scalar>,
    children: Box<Vec<RTreeNode<T>>>,
    depth: usize,
    options: Arc<RTreeOptions>,
}

#[doc(hidden)]
#[derive(Clone)]
pub enum RTreeNode<T> where T: SpatialObject {
    Leaf(T),
    DirectoryNode(DirectoryNodeData<T>),
}

/// A pure rust r*-tree implementation.
///
/// Note that the structure does use the [r-star-tree heuristics](https://en.wikipedia.org/wiki/R*_tree), even though it is called
/// `RTree`.
/// Allows to insert various objects into the tree and to make some common queries.
#[derive(Clone)]
pub struct RTree<T> where T: SpatialObject {

    root: DirectoryNodeData<T>,
    size: usize,
}

impl<T> Default for RTree<T> where T: SpatialObject {
    fn default() -> RTree<T> {
        RTree::new()
    }
}

impl<T> RTree<T> where T: SpatialObject {
    /// Creates an empty r*-tree.
    pub fn new() -> RTree<T> {
        RTree::new_with_options(Default::default())
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
    
    /// Inserts a new element into the tree.
    ///
    /// This will require `O(log(n))` operations on average, where n is the number of elements contained
    /// in the tree.
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
    pub fn lookup_and_remove(&mut self, query_point: [T::Scalar; 2]) -> Option<T> {
        let result = self.root.lookup_and_remove(&query_point.into());
        if result.is_some() {
            self.size -= 1;
        }
        result
    }

    #[doc(hidden)]
    pub fn root(&self) -> &DirectoryNodeData<T> {
        &self.root
    }

    /// Searches for an element at a given position.
    ///
    /// If `query_point` is contained by one object in the tree, this object will be returned.
    /// If multiple objects contain the point, only one of them will be returned.
    pub fn lookup(&self, point: [T::Scalar; 2]) -> Option<&T> {
        self.root.lookup(&point.into())
    }
    
    /// Returns the nearest neighbor.
    ///
    /// Returns `None` if the tree is empty.
    pub fn nearest_neighbor(&self, query_point: [T::Scalar; 2]) -> Option<&T> {
        self.root.nearest_neighbor(&query_point.into(), Bounded::max_value())
    }

    /// Returns the nearest n neighbors.
    pub fn nearest_n_neighbors(&self, query_point: [T::Scalar; 2], n: usize) -> Vec<&T> {
        let mut result = Vec::new();
        self.root.nearest_n_neighbors(&query_point.into(), n, &mut result);
        result
    }

    /// Returns all objects (partially) contained in a circle.
    ///
    /// Note that `radius2` is the circle's squared radius, not the actual radius.
    /// An object is contained if a part of it lies within the circle.
    pub fn lookup_in_circle(&self, circle_origin: [T::Scalar; 2], radius2: T::Scalar) -> Vec<&T> {
        let mut result = Vec::new();
        self.root.lookup_in_circle(&mut result, &circle_origin.into(), &radius2);
        result
    }
    
}

impl <T> RTree<T> where T: SpatialObject + PartialEq {

    /// Removes an object from the tree.
    ///
    /// Locates and removes an object from the tree, returning
    /// `true` if the element could be removed.
    /// If multiple object's are equal to `to_remove`, only one
    /// will be deleted.
    pub fn remove(&mut self, obj: &T) -> bool {
        let result = self.root.remove(obj);
        if result {
            self.size -= 1;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::{RTree};
    use primitives::SimpleTriangle;
    use cgmath::{Vector2, InnerSpace};
    use cgmath::conv::{array2};
    use num::Float;
    use testutils::*;

    #[test]
    fn test_tree_with_integral_vectors() {
        // This test should compile
        let mut tree = RTree::new();
        tree.insert(Vector2::new(13, 37));
        assert!(tree.lookup([13, 37]).is_some())
    }

    #[test]
    fn test_nearest_neighbor() {
        let (tree, points) = create_random_tree::<f32>(1000, [10, 233, 588812, 411112]);
        let sample_points = random_points_with_seed(100, [66, 123, 12345, 112]);
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
            assert!(nearest == tree.nearest_neighbor(array2(sample_point.clone())));
        }
    }

    #[test]
    fn test_lookup() {
        let (tree, points) = create_random_tree::<f32>(10000, [9, 8, 7, 6]);
        let sample_points = random_points_with_seed(1000, [2, 1, 0, 3]);
        for sample_point in &sample_points {
            assert!(!tree.lookup(array2(sample_point.clone())).is_some());
        }
        for point in &points {
            assert!(tree.lookup(array2(point.clone())) == Some(point));
        }
    }

    #[test]
    fn test_lookup_and_remove() {
        let (mut tree, points) = create_random_tree::<f32>(10000, [3141, 59, 26, 53]);
        let sample_points = random_points_with_seed(1000, [2, 3, 0, 22991]);
        for sample_point in &sample_points {
            assert!(!tree.lookup_and_remove(array2(sample_point.clone())).is_some());
        }

        // Test if all points are still there
        for point in &points {
            assert_eq!(tree.lookup(array2(point.clone())), Some(point));
        }
        // Now remove all points
        for point in &points {
            assert_eq!(tree.lookup_and_remove(array2(point.clone())).as_ref(), Some(point));
        }
        assert!(tree.root.children.is_empty());
    }

    #[test]
    fn test_remove() {
        let random_points = random_points_with_seed(300, [911, 110, 123, 454]);
        let mut triangles = Vec::with_capacity(100);
        for ps in random_points.chunks(3) {
            let ps = [ps[0], ps[1], ps[2]];
            triangles.push(SimpleTriangle::new(ps));
        }
        let mut tree = RTree::new();
        for triangle in triangles.iter().cloned() {
            tree.insert(triangle);
        }
        // Try to remove a triangle that is not contained
        let triangle = SimpleTriangle::new([Vector2::new(0.0, 0.0), 
                                            Vector2::new(1.0, 0.0), 
                                            Vector2::new(1.0, 1.0)]);
        assert!(!tree.remove(&triangle));
        let mut size = 100usize;
        for triangle in triangles.iter() {
            assert!(tree.remove(triangle));
            size -= 1;
            assert_eq!(tree.size(), size);
        }
    }

    #[test]
    fn test_iteration() {
        let (tree, reference_points) = create_random_tree::<f32>(100, [1, 10, 100, 1000]);
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
}
