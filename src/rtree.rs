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

use misc::min_inline;
use std::sync::Arc;
use traits::{SpatialObject};
use vector_traits::{VectorN, VectorNExtensions};
use num::{zero};
use boundingvolume::BoundingRect;
use std::iter::Once;
use std::borrow::{Borrow, BorrowMut};

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

    pub fn build<T: SpatialObject, B: Borrow<T>>(self) -> RTree<T, B> {
        RTree::new_with_options(self)
    }
}

#[doc(hidden)]
pub struct NearestNeighborIterator<'a, T, B>
    where T: SpatialObject + 'a,
          B: Borrow<T> + 'a {
    next_items: Vec<&'a T>,
    min_dist: <T::Vector as VectorN>::Scalar,
    tree: &'a RTree<T, B>,
    query_point: T::Vector,
}

impl<'a, T, B> NearestNeighborIterator<'a, T, B>
    where T: SpatialObject + 'a, 
          B: Borrow<T>
{
    fn new(tree: &'a RTree<T, B>, query_point: T::Vector) -> NearestNeighborIterator<'a, T, B> {
        NearestNeighborIterator {
            next_items: Vec::new(),
            min_dist: zero(),
            tree: tree,
            query_point: query_point,
        }
    }
}

impl <'a, T, B> Iterator for NearestNeighborIterator<'a, T, B> 
    where T: SpatialObject + 'a,
          B: Borrow<T>
{
    type Item = &'a T;
    
    fn next(&mut self) -> Option<&'a T> {
        match self.next_items.pop() {
            Some(t) => Some(t),
            None => {
                // Get new nearest elements
                let min_dist_opt = self.tree.root.nearest_neighbors_with_min_dist(
                    &self.query_point, None, self.min_dist.clone(), &mut self.next_items);
                if let Some(min_dist) = min_dist_opt {
                    self.min_dist = min_dist;
                }
                self.next_items.pop()
            },
        }
    }
}

/// Iterates over all entries in an r-tree.
/// Returned by `RTree::iter()`
pub struct RTreeIterator<'a, T, B> 
    where T: SpatialObject + 'a,
          B: Borrow<T> + 'a
    
{
    data: &'a DirectoryNodeData<T, B>,
    cur_index: usize, 
    cur_iterator: Option<Box<RTreeNodeIterator<'a, T, B>>>,
}

#[allow(missing_docs)]
pub enum RTreeNodeIterator<'a, T, B> 
    where T: SpatialObject + 'a,
          B: Borrow<T> + 'a
{
    LeafIterator(Once<&'a T>),
    DirectoryNodeIterator(RTreeIterator<'a, T, B>),
}

impl <'a, T, B> RTreeIterator<'a, T, B> 
    where T: SpatialObject,
          B: Borrow<T>
{
    fn new(data: &'a DirectoryNodeData<T, B>) -> RTreeIterator<'a, T, B> {
        RTreeIterator {
            data: data,
            cur_index: 0,
            cur_iterator: data.children.first().map(
                |child| Box::new(RTreeNodeIterator::new(child))),
        }
    }
}

impl <'a, T, B> Iterator for RTreeIterator<'a, T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
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

impl <'a, T, B> RTreeNodeIterator<'a, T, B>
    where T: SpatialObject,
          B: Borrow<T>
{

    fn new(node: &'a RTreeNode<T, B>) -> RTreeNodeIterator<'a, T, B> {
        use RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match node {
            &RTreeNode::Leaf(ref b) => LeafIterator(::std::iter::once(b.borrow())),
            &RTreeNode::DirectoryNode(ref data) => 
                DirectoryNodeIterator(RTreeIterator::new(data)),
        }
    }
}

impl <'a, T, B> Iterator for RTreeNodeIterator<'a, T, B>
    where T: SpatialObject,
          B: Borrow<T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        use RTreeNodeIterator::{LeafIterator, DirectoryNodeIterator};
        match self {
            &mut LeafIterator(ref mut once) => once.next(),
            &mut DirectoryNodeIterator(ref mut iter) => iter.next(),
        }
    }
}

#[doc(hidden)]
impl <T, B> DirectoryNodeData<T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
    pub fn children(&self) -> &Vec<RTreeNode<T,B>> {
        &self.children
    }

    pub fn depth(&self) -> usize {
        self.depth
    }

    pub fn mbr(&self) -> BoundingRect<T::Vector> {
        self.bounding_box.clone().unwrap()
    }

    fn new(depth: usize, options: Arc<RTreeOptions>) -> DirectoryNodeData<T, B> {
        DirectoryNodeData {
            bounding_box: None,
            children: Box::new(Vec::with_capacity(options.max_size + 1)),
            options: options,
            depth: depth,
        }
    }

    fn new_parent(mut children: Box<Vec<RTreeNode<T, B>>>, depth: usize, options: Arc<RTreeOptions>
                  ) -> DirectoryNodeData<T, B> {
        let missing = options.max_size + 1 - children.len();
        children.reserve_exact(missing);
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
    fn update_mbr_with_element(&mut self, element_bb: &BoundingRect<T::Vector>) {
        if let Some(ref mut bb) = self.bounding_box {
            bb.add_rect(element_bb);
        }  else {
            self.bounding_box = Some(element_bb.clone());
        }
    }

    fn insert(&mut self, t: RTreeNode<T, B>, state: &mut InsertionState) -> InsertionResult<T, B> {
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

    fn resolve_overflow(&mut self, state: &mut InsertionState) -> InsertionResult<T, B> {
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

    #[inline(never)]
    fn split(&mut self) -> RTreeNode<T, B> {
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
        let offsplit = Box::new(self.children.split_off(best_index));
        let result = RTreeNode::DirectoryNode(DirectoryNodeData::new_parent(offsplit, self.depth,
                                                                            self.options.clone()));
        self.update_mbr();
        result
    }

    #[inline(never)]
    fn reinsert(&mut self) -> Vec<RTreeNode<T, B>> {
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
        for axis in 0 .. T::Vector::dimensions() {
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

    fn choose_subtree(&mut self, node: &RTreeNode<T, B>) -> &mut DirectoryNodeData<T, B> {
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
                    let mut overlap: <T::Vector as VectorN>::Scalar = zero();
                    let mut new_overlap: <T::Vector as VectorN>::Scalar = zero();
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

    fn add_children(&mut self, mut new_children: Vec<RTreeNode<T, B>>) {
        if let &mut Some(ref mut bb) = &mut self.bounding_box {
            for child in &new_children {
                bb.add_rect(&child.mbr());
            }
            self.children.append(&mut new_children);
            return;
        } 
        if let Some(first) = new_children.first() {
            let mut bb = first.mbr();
            for child in new_children.iter().skip(1) {
                bb.add_rect(&child.mbr());
            }
            self.bounding_box = Some(bb);
        }
        self.children.append(&mut new_children);
    }

    fn close_neighbor(&self, point: &T::Vector) -> Option<&T> {
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
                    return Some(t.borrow())
                }
            }
        }
    }

    fn nearest_neighbor(&self, point: &T::Vector, 
                        mut nearest_distance: Option<<T::Vector as VectorN>::Scalar>) -> Option<&T> {
        let mut nearest = None;
        // Calculate smallest minmax-distance
        let mut smallest_min_max: <T::Vector as VectorN>::Scalar = zero();
        let mut first = true;
        for child in self.children.iter() {
            let new_min = child.mbr().min_max_dist2(point);
            smallest_min_max = if first {
                first = false;
                new_min
            } else {
                min_inline(smallest_min_max, new_min)
            };
        }
        let mut sorted = Vec::with_capacity(self.children.len());
        for child in self.children.iter() {
            let min_dist = child.mbr().min_dist2(point);
            if min_dist <= smallest_min_max {
                sorted.push((child, min_dist));
            }
        }
        sorted.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());

        for &(child, ref min_dist) in sorted.iter() {
            if nearest_distance.clone().map(|d| min_dist.clone() > d).unwrap_or(false) {
                // Prune this element
                break;
            }
            match child.nearest_neighbor(point, nearest_distance.clone()) {
                Some(t) => {
                    nearest_distance = Some(t.distance2(point));
                    nearest = Some(t);
                },
                None => {}
            }
        }
        nearest
    }

    fn nearest_neighbors<'a>(&'a self, point: &T::Vector,
                         mut nearest_distance: Option<<T::Vector as VectorN>::Scalar>, 
                             result: &mut Vec<&'a T>) -> Option<<T::Vector as VectorN>::Scalar> {
        // Calculate smallest minmax-distance
        let mut smallest_min_max: <T::Vector as VectorN>::Scalar = zero();
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

    fn nearest_neighbors_with_min_dist<'a>(
        &'a self, point: &T::Vector,
        mut nearest_distance: Option<<T::Vector as VectorN>::Scalar>, 
        min_dist: <T::Vector as VectorN>::Scalar, 
        result: &mut Vec<&'a T>) 
        -> Option<<T::Vector as VectorN>::Scalar> 
    {
        // Calculate smallest minmax-distance
        let mut smallest_min_max: <T::Vector as VectorN>::Scalar = zero();
        for child in self.children.iter() {
            let mbr = child.mbr();
            let min_point = mbr.min_point(point);
            for dim in 0 .. T::Vector::dimensions() {
                for mut p in [mbr.lower(), mbr.upper()].iter_mut() {
                    *p.nth_mut(dim) = min_point.nth(dim).clone();
                    let distance = p.sub(point).length2();
                    if distance < smallest_min_max || dim == 0 {
                        smallest_min_max = distance;
                    }
                }
            }
        }
        let mut sorted: Vec<_> = self.children.iter().collect();
        sorted.sort_by(|l, r| l.mbr().min_dist2(point).partial_cmp(
            &r.mbr().min_dist2(point)).unwrap());
        for child in sorted.iter() {
            let cur_min_dist = child.mbr().min_dist2(point);
            let max_dist = child.mbr().max_dist2(point);
            if max_dist < min_dist.clone() ||
                cur_min_dist > smallest_min_max || 
                nearest_distance.clone().map(|d| cur_min_dist > d).unwrap_or(false) {
                // Prune this element
                continue;
            }
            match child.nearest_neighbors_with_min_dist(point, nearest_distance.clone(), min_dist.clone(), result) {
                Some(nearest) => {
                    nearest_distance = Some(nearest);
                },
                None => {}
            }
        }
        nearest_distance
    }

    fn nearest_n_neighbors<'a>(&'a self, point: &T::Vector, n: usize, result: &mut Vec<&'a T>) {

        for child in self.children.iter() {
            let min_dist = child.mbr().min_dist2(point);
            if result.len() == n && min_dist >= result.last().unwrap().distance2(point) {
                // Prune this element
                continue;
            }
            match child {
                &RTreeNode::DirectoryNode(ref data) => {
                    data.nearest_n_neighbors(point, n, result);
                },
                &RTreeNode::Leaf(ref b) => {
                    let distance = b.borrow().distance2(point);
                    if result.len() != n || distance < result.last().unwrap().distance2(point) {
                        if result.len() == n {
                            result.pop();
                        }
                        let index = match result.binary_search_by(|e| e.distance2(point).partial_cmp(
                            &distance).unwrap()) {
                            Ok(index) => index,
                            Err(index) => index,
                        };
                        result.insert(index, b.borrow());
                    }
                }
            }
        }
    }

    fn lookup_and_remove(&mut self, point: &T::Vector) -> Option<B> {
        let contains = self.bounding_box.as_ref().map(|bb | bb.contains_point(point)).unwrap_or(false);
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
                    RTreeNode::Leaf(b) => {
                        if b.borrow().contains(point) {
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
    
    fn lookup(&self, point: &T::Vector) -> Option<&T> {
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
                            if obj.borrow().contains(point) {
                                return Some(obj.borrow());
                            }
                        },
                    }
                }
            }
        }
        None
    }

    fn lookup_in_circle<'b>(&'b self, result: &mut Vec<&'b T>, origin: &T::Vector,
                            radius2: &<T::Vector as VectorN>::Scalar)
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
                    if t.borrow().distance2(origin) < *radius2 {
                        result.push(t.borrow());
                    }
                },
            }
        }
    }

    fn lookup_in_rectangle<'b>(&'b self, result: &mut Vec<&'b T>,
                               query_rect: &BoundingRect<T::Vector>) {
        for child in self.children.iter().filter(|c| c.mbr().intersects(query_rect)) {
            match child {
                &RTreeNode::DirectoryNode(ref data) => data.lookup_in_rectangle(result, query_rect),
                &RTreeNode::Leaf(ref t) => {
                    if t.borrow().mbr().intersects(query_rect) {
                        result.push(t.borrow());
                    }
                }
            }
        }
    }
}

impl <T, B> DirectoryNodeData<T, B>
    where T: SpatialObject,
          B: BorrowMut<T>
{
    fn lookup_mut(&mut self, point: &T::Vector) -> Option<&mut T> {
        let mut todo_list = Vec::with_capacity(40);
        todo_list.push(self);
        while let Some(mut next) = todo_list.pop() {
            if next.mbr().contains_point(point) {
                for child in next.children.iter_mut() {
                    match child {
                        &mut RTreeNode::DirectoryNode(ref mut data) => {
                            todo_list.push(data);
                        },
                        &mut RTreeNode::Leaf(ref mut obj) => {
                            if (*obj).borrow().contains(point) {
                                return Some(obj.borrow_mut());
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
impl <T, B> DirectoryNodeData<T, B>
    where T: SpatialObject + PartialEq,
          B: Borrow<T>
{

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
                        if t.borrow() == to_remove {
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
                        if t.borrow() == obj {
                            return true
                        }
                    }
                }
            }
        }
        false
    }
}

enum InsertionResult<T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
    Complete,
    Split(RTreeNode<T, B>),
    Reinsert(Vec<RTreeNode<T, B>>),
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
impl <T, B> RTreeNode<T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
    pub fn depth(&self) -> usize {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.depth,
            _ => 0
        }
    }

    pub fn mbr(&self) -> BoundingRect<T::Vector> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.bounding_box.clone().unwrap(),
            &RTreeNode::Leaf(ref t) => t.borrow().mbr(),
        }
    }

    fn nearest_neighbor(&self, point: &T::Vector, nearest_distance: Option<<T::Vector as VectorN>::Scalar>) 
                        -> Option<&T> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.nearest_neighbor(point, nearest_distance),
            &RTreeNode::Leaf(ref t) => {
                let distance = t.borrow().distance2(point);
                if nearest_distance.map(|d| distance < d).unwrap_or(true) {
                    Some(t.borrow())
                } else {
                    None
                }
            }
        }
    }

    fn nearest_neighbors<'a>(&'a self, point: &T::Vector, 
                             nearest_distance: Option<<T::Vector as VectorN>::Scalar>,
                             result: &mut Vec<&'a T>) -> Option<<T::Vector as VectorN>::Scalar> {
        match self {
            &RTreeNode::DirectoryNode(ref data) => data.nearest_neighbors(point, nearest_distance, result),
            &RTreeNode::Leaf(ref t) => {
                let distance = t.borrow().distance2(point);
                match nearest_distance {
                    Some(nearest) => {                
                        if distance <= nearest {
                            if distance < nearest {
                                // We've found a new minimum element, remove all other neighbors found so far
                                result.clear();
                            }
                            result.push(t.borrow());
                            Some(distance)
                        } else {
                            // This object is not among the nearest neigbors
                            None
                        }
                    },
                    None => {
                        result.push(t.borrow());
                        Some(distance)
                    }
                }
            }
        }
    }

    fn nearest_neighbors_with_min_dist<'a>(
        &'a self, point: &T::Vector, nearest_distance: Option<<T::Vector as VectorN>::Scalar>,
        min_dist: <T::Vector as VectorN>::Scalar, result: &mut Vec<&'a T>)
        -> Option<<T::Vector as VectorN>::Scalar> 
    {
        match self {
            &RTreeNode::DirectoryNode(ref data) => 
                data.nearest_neighbors_with_min_dist(point, nearest_distance, min_dist, result),
            &RTreeNode::Leaf(ref t) => {
                let distance = t.borrow().distance2(point);
                if distance <= min_dist {
                    return None;
                }
                match nearest_distance {
                    Some(nearest) => {                
                        if distance <= nearest {
                            if distance < nearest {
                                // We've found a new minimum element, remove all other neighbors found so far
                                result.clear();
                            }
                            result.push(t.borrow());
                            Some(distance)
                        } else {
                            // This object is not among the nearest neigbors
                            None
                        }
                    },
                    None => {
                        result.push(t.borrow());
                        Some(distance)
                    }
                }
            }
        }
    }
}

#[doc(hidden)]
#[derive(Clone)]
pub struct DirectoryNodeData<T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
    bounding_box: Option<BoundingRect<T::Vector>>,
    children: Box<Vec<RTreeNode<T, B>>>,
    depth: usize,
    options: Arc<RTreeOptions>,
}

#[doc(hidden)]
#[derive(Clone)]
pub enum RTreeNode<T, B>
    where T: SpatialObject,
          B: Borrow<T>
{
    Leaf(B),
    DirectoryNode(DirectoryNodeData<T, B>),
}


/// A rust implementation of n dimensional r*-trees
///
/// [R-trees](https://en.wikipedia.org/wiki/R-tree) provide efficient nearest-neighbor searches for
/// many objects. [R*-trees](https://en.wikipedia.org/wiki/R*_tree) (&quot;R-Star-Trees&quot;) 
/// are a common variant of r-trees and use more advanced heuristics to improve query performance. This
/// struct implements r*-trees, despite its name.
/// Instead of linear time complexity, r-trees yield logarithmic complexity
/// for look-up operations and nearest neighbor queries. Inserting into an r-tree runs in O(log(n)) time on average.
/// Some simple simple geometric primitives that can be inserted into an r-tree can be found in 
/// the `primitives` module. If your object is not among those, consider
/// implementing the `SpatialObject` trait.
/// 
/// Note that the `rtree`-structures work with fixed arrays of size 2, 3 or 4 or
/// with the vector types provided by the `nalgebra` and `cgmath` packages.
/// Both integral and floating point scalar types are supported.
///
/// ```
/// # extern crate nalgebra;
/// # extern crate spade;
///
/// # use nalgebra::{Vector4};
/// # use spade::RTree;
///
/// # fn main() {
///   let mut tree = RTree::new();
///   tree.insert(Vector4::new(13i32, 10, 10, 37));
/// # }
/// ```
/// # Basic Example
///
/// ```
/// extern crate cgmath; // Alternatively: use nalgebra or [f32; 2]
/// extern crate spade;
///
/// use spade::RTree;
/// use cgmath::Vector2;
///
/// fn main() {
/// let mut rtree = RTree::new();
/// // Insert two points
/// rtree.insert(Vector2::new(0.5, 0.5f32));
/// rtree.insert(Vector2::new(1.0, 1.0f32));
///
/// if rtree.lookup(&Vector2::new(0.5, 0.5)).is_some() {
///   println!("We'fe found a point at [0.5, 0.5]/");
/// }
/// 
/// let nearest = rtree.nearest_neighbor(&Vector2::new(1.5, 1.5)).unwrap();
/// println!("nearest neighbor at [1.5, 1.5]: {:?}", nearest);
///
/// // Iterate over all elements
/// for point in rtree.iter() {
///   println!("Found point: {:?}", point);
/// }
/// }
/// ```
/// # Type parameters
/// `RTree` has two type parameters `T` and `B`. `T` can be any
/// type that implements `SpatialObject`. If you do not want to
/// insert types with borrow semantics, `B` and `T` are equal.
///
/// ## Inserting pointer like types
/// It is possible to insert pointer like types (like references, boxes, arcs...)
/// that implement `Borrow<T>`. You may need to specify the r-tree's type parameters
/// more precisely as seen in this example:
///
/// ```
/// # extern crate nalgebra;
/// # extern crate spade;
///
/// # use nalgebra::{Vector2};
/// # use spade::RTree;
///
/// # fn main() {
///   let mut tree: RTree<Vector2<_>, Box<_>> = RTree::new();
///   tree.insert(Box::new(Vector2::new(20i32, 10)));
/// # }
/// ```

#[derive(Clone)]
pub struct RTree<T, B = T> where T: SpatialObject, B: Borrow<T> {

    root: DirectoryNodeData<T, B>,
    size: usize,
}

impl<T, B> Default for RTree<T, B> where T: SpatialObject, B: Borrow<T> {
    fn default() -> RTree<T, B> {
        RTree::new()
    }
}

impl<T, B> RTree<T, B> 
    where T: SpatialObject,
          B: Borrow<T>
{
    /// Creates an empty r*-tree.
    pub fn new() -> RTree<T, B> {
        RTree::new_with_options(Default::default())
    }

    /// Returns the trees minimal bounding box.
    pub fn mbr(&self) -> Option<BoundingRect<T::Vector>> {
        self.root.bounding_box.clone()
    }

    #[doc(hidden)]
    pub fn new_with_options(options: RTreeOptions) -> RTree<T, B> {
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
    pub fn iter(&self) -> RTreeIterator<T, B> {
        RTreeIterator::new(&self.root)
    }
    
    /// Inserts a new element into the tree.
    ///
    /// This will require `O(log(n))` operations on average, where n is the number of elements contained
    /// in the tree.
    pub fn insert(&mut self, t: B) {
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

    #[doc(hidden)]
    pub fn root(&self) -> &DirectoryNodeData<T, B> {
        // This access is only needed for one of the examples
        &self.root
    }

    /// Searches for an element at a given position.
    ///
    /// If `query_point` is contained by one object in the tree, this object will be returned.
    /// If multiple objects contain the point, only one of them will be returned.
    pub fn lookup(&self, query_point: &T::Vector) -> Option<&T> {
        if self.size > 0 {
            self.root.lookup(query_point)
        } else {
            None
        }
    }

    /// Returns the nearest neighbor.
    ///
    /// Returns `None` if the tree is empty.
    pub fn nearest_neighbor(&self, query_point: &T::Vector) -> Option<&T> {
        if self.size > 0 {
            self.root.nearest_neighbor(query_point, None)
        } else {
            None
        }
    }

    /// Returns an object close to a given point. This operation is faster than
    /// `nearest_neighbor` but will not neccessarily yield the real nearest neighbor.
    #[inline(never)]
    pub fn close_neighbor(&self, query_point: &T::Vector) -> Option<&T> {
        if self.size > 0 {
            self.root.close_neighbor(query_point)
        } else {
            None
        }
    }

    /// Returns the nearest neighbors of a given point.
    ///
    /// All returned values will have the exact same distance from the given query point.
    /// Returns an empty vector if the tree is empty.
    pub fn nearest_neighbors(&self, query_point: &T::Vector) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.nearest_neighbors(query_point, None, &mut result);
        }
        result
    }

    /// Returns the nearest n neighbors.
    pub fn nearest_n_neighbors(&self, query_point: &T::Vector, n: usize) -> Vec<&T> {
        // let iter= NearestNeighborIterator::new(self, query_point);
        // Iterator::collect(iter.take(n))

        let mut result = Vec::new();
        if self.size > 0 {
            self.root.nearest_n_neighbors(query_point, n, &mut result);
        }
        result
    }

    /// Returns all objects (partially) contained in a rectangle
    pub fn lookup_in_rectangle(&self, query_rect: &BoundingRect<T::Vector>) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.lookup_in_rectangle(&mut result, query_rect);
        }
        result
    }

    // TODO: This doesn't work yet.
    #[doc(hidden)]
    /// Returns an iterator that iterates over the nearest elements next to a
    /// query point.
    pub fn nearest_neighbor_iterator(
        &self, query_point: T::Vector) -> NearestNeighborIterator<T, B> {
        NearestNeighborIterator::new(self, query_point)
    }

    /// Returns all objects (partially) contained in a circle.
    ///
    /// Note that `radius2` is the circle's squared radius, not the actual radius.
    /// An object is contained if a part of it lies within the circle.
    pub fn lookup_in_circle(&self, circle_origin: &T::Vector, 
                            radius2: &<T::Vector as VectorN>::Scalar) -> Vec<&T> {
        let mut result = Vec::new();
        if self.size > 0 {
            self.root.lookup_in_circle(&mut result, circle_origin.into(), radius2);
        }
        result
    }

    /// Searches for an element and removes it.
    ///
    /// If the given point is contained by one object in the tree, this object is being removed
    /// and returned. If the point is contained by multiple objects, only one of them is removed and
    /// returned.
    pub fn lookup_and_remove(&mut self, query_point: &T::Vector) -> Option<B> {
        if self.size > 0 {
            let result = self.root.lookup_and_remove(query_point);
            if result.is_some() {
                self.size -= 1;
            }
            result
        } else {
            None
        }
    }
}

impl<T, B> RTree<T, B> 
    where T: SpatialObject,
          B: BorrowMut<T>
{
    /// Searches for an element at a given position and returns a mutable
    /// reference.
    /// If `query_point` is contained by multiple objects in the tree,
    /// one of them will be returned.
    /// *Do not change the object's minimal bounding box*.
    pub fn lookup_mut(&mut self, query_point: &T::Vector) -> Option<&mut T> {
        if self.size > 0 {
            self.root.lookup_mut(query_point)
        } else {
            None
        }
    }
}

impl <T, B> RTree<T, B>
    where T: SpatialObject + PartialEq,
          B: Borrow<T>
{

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

#[cfg(test)]
mod test {
    use super::{RTree};
    use boundingvolume::BoundingRect;
    use primitives::SimpleTriangle;
    use cgmath::{Vector2, InnerSpace};
    use num::Float;
    use testutils::*;

    #[test]
    fn test_tree_with_integral_vectors() {
        // This test should compile
        let mut tree = RTree::new();
        tree.insert(Vector2::new(13, 37));
        assert!(tree.lookup(&Vector2::new(13, 37)).is_some())
    }

    #[test]
    fn test_tree_with_array_vectors() {
        // This test should compile
        let mut tree = RTree::<[i32; 3], _>::new();
        tree.insert([13i32, 37, 12]);
        assert!(tree.lookup(&[13, 37, 12]).is_some())
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
            assert!(nearest == tree.nearest_neighbor(sample_point));
        }
    }

    #[test]
    fn test_lookup_in_circle() {
        let (tree, points) = create_random_tree::<f32>(1000, [10, 233, 588812, 411112]);
        let sample_points = random_points_with_seed(100, [66, 123, 12345, 112]);
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
        let (tree, points) = create_random_tree::<f32>(1000, [10, 233, 588812, 411112]);
        let sample_points = random_points_with_seed(100, [66, 123, 12345, 112]);
        const SIZE: f32 = 20.;
        for sample_point in &sample_points {
            let sample_rect = BoundingRect::from_corners(
                sample_point, &(*sample_point + Vector2::new(SIZE, SIZE)));
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
        assert!(tree.nearest_neighbors(&Vector2::new(1, 0)).is_empty());
        tree.insert(Vector2::new(1, 0));
        tree.insert(Vector2::new(0, 1));
        tree.insert(Vector2::new(-1, 0));
        tree.insert(Vector2::new(0, -1));
        tree.insert(Vector2::new(3, 0));
        tree.insert(Vector2::new(2, 1));
        tree.insert(Vector2::new(2, -1));
        assert_eq!(tree.nearest_neighbors(&Vector2::new(0, 0)).len(), 4);
        assert_eq!(tree.nearest_neighbors(&Vector2::new(1, 0)).len(), 1);
        assert_eq!(tree.nearest_neighbors(&Vector2::new(2, 0)).len(), 4);
    }

    #[test]
    fn test_lookup() {
        let (mut tree, mut points) = create_random_tree::<f32>(10000, [9, 8, 7, 6]);
        let sample_points = random_points_with_seed(1000, [2, 1, 0, 3]);
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
        let (mut tree, points) = create_random_tree::<f32>(10000, [3141, 59, 26, 53]);
        let sample_points = random_points_with_seed(1000, [2, 3, 0, 22991]);
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
    }

    #[test]
    fn test_remove() {
        let random_points = random_points_with_seed(300, [911, 110, 123, 454]);
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
            let triangle = SimpleTriangle::new(Vector2::new(0.0, 0.0), 
                                               Vector2::new(1.0, 0.0), 
                                               Vector2::new(1.0, 1.0));
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
    
    // TODO: NN iteration doesn't work yet
    #[ignore]
    #[test]
    fn test_nearest_neighbor_iteration() {
        use traits::SpatialObject;
        let (tree, mut points) = create_random_tree::<f32>(100, [9, 8, 7, 6]);
        let query_point = Vector2::new(0., 0.);
        let sorted: Vec<_> = tree.nearest_neighbor_iterator(query_point).cloned().collect();
        points.sort_by(|l, r| l.distance2(&query_point).partial_cmp(&r.distance2(&query_point)).unwrap());
        assert_eq!(sorted, points);
    }


    #[test]
    fn test_higher_dimensions() {
        use cgmath::Vector4;
        use rand::{XorShiftRng, SeedableRng, Rng};
        let mut tree: RTree<Vector4<f32>, Vector4<f32>> = RTree::new();
        let mut rng = XorShiftRng::from_seed([1, 2, 3, 1992]);
        let mut entries = Vec::new();
        for _ in 0 .. 1000 {
            let (x, y, z, w) = (rng.next_f32(), rng.next_f32(), rng.next_f32(), rng.next_f32());
            let entry = Vector4::new(x, y, z, w);
            entries.push(entry);
            tree.insert(entry);
        }

        for entry in &entries {
            assert!(tree.lookup(entry).is_some());
            assert_eq!(tree.nearest_neighbor(entry), Some(entry))
        }
    }

    #[test]
    fn test_insert_borrow_type() {
        // Just check if this compiles
        let mut tree: RTree<Vector2<_>, Box<_>> = RTree::new();
        tree.insert(Box::new(Vector2::new(20i32, 10)));
    }
}

