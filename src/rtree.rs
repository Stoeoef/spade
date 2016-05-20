use cgmath::{Vector2, EuclideanVector};
use std::sync::Arc;
use traits::{SpatialObject};
use num::{Float, zero};
use boundingvolume::BoundingRect;
use std::iter::Once;

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
        assert!(self.min_size > reinsertion_count);
        self.reinsertion_count = reinsertion_count;
        self
    }

    pub fn build<'a, T: SpatialObject>(self) -> RTree<T> {
        RTree::new_with_options(self)
    }
}


pub struct RTreeIterator<'a, T> where T: SpatialObject + 'a {
    data: &'a DirectoryNodeData<T>,
    cur_index: usize, 
    cur_iterator: Option<Box<RTreeNodeIterator<'a, T>>>,
}

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
                let offsplit = self.split();
                InsertionResult::Split(offsplit)
            } else {
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
        let mut best = (Float::infinity(), Float::infinity());
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
            (l_center - center).length2().partial_cmp(&(r_center - center).length2()).unwrap()
        });
        let num_children = self.children.len();
        let result = self.children.split_off(num_children - self.options.reinsertion_count);
        self.update_mbr();
        result
    }

    fn get_split_axis(&mut self) -> usize {
        let mut best_goodness = Float::infinity();
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
        assert!(self.depth >= 2, "Cannot choose subtree if all children are leaves");
        let insertion_mbr = node.mbr();
        let mut inclusion_count = 0;
        let mut min_area = Float::infinity();
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

            let mut min = (Float::infinity(), Float::infinity(), Float::infinity());

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
        let mut smallest_min_max: T::Scalar = Float::infinity();
        for child in self.children.iter() {
            smallest_min_max = smallest_min_max.min(child.mbr().min_max_dist(point));
        }
        let mut sorted = Vec::with_capacity(self.children.len());
        sorted.extend(self.children.iter().map(|e| e));
        sorted.sort_by(|l, r| l.mbr().min_dist(point).partial_cmp(&r.mbr().min_dist(point))
                       .unwrap());
        for child in sorted.iter() {
            let min_dist = child.mbr().min_dist(point);
            if min_dist > smallest_min_max || min_dist > nearest_distance {
                // Prune this element
                continue;
            }
            match child.nearest_neighbor(point, nearest_distance) {
                Some(t) => {
                    nearest_distance = t.distance2(point);
                    nearest = Some(t);
                },
                None => {}
            }
        }
        nearest
    }

    pub fn lookup_and_remove(&mut self, point: &Vector2<T::Scalar>) -> Option<T> {
        let contains = self.bounding_box.contains_point(point);
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
                        if t.contains(point) {
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
    
    pub fn lookup(&self, point: &Vector2<T::Scalar>) -> Option<&T> {
        if self.bounding_box.contains_point(point) {
            for child in self.children.iter() {
                match child {
                    &RTreeNode::DirectoryNode(ref data) => {
                        let result = data.lookup(point);
                        if result.is_some() {
                            return result;
                        }
                    },
                    &RTreeNode::Leaf(ref t) => {
                        if t.contains(point) {
                            return Some(t);
                        }
                    }
                }
            }
        }
        None
    }
}

impl <'a, T> DirectoryNodeData<T> where T: SpatialObject + PartialEq {

    pub fn remove(&mut self, obj: &T) -> bool {
        let contains = self.bounding_box.contains_rect(&obj.mbr());
        if contains {
            let mut children = ::std::mem::replace(&mut self.children, 
                                                   Box::new(Vec::new()));
            let mut result = false;
            for child in children.drain(..) {
                match child {
                    RTreeNode::DirectoryNode(mut data) => {
                        result = data.remove(obj) | result;
                        if !data.children.is_empty() {
                            // Don't add a node if it has become empty
                            self.children.push(RTreeNode::DirectoryNode(data));
                        }
                    },
                    RTreeNode::Leaf(t) => {
                        if t == *obj {
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

    fn depth(&self) -> usize {
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
                let distance = t.distance2(point);
                if distance < nearest_distance {
                    Some(t)
                } else {
                    None
                }
            }
        }
    }
}

struct DirectoryNodeData<T> where T: SpatialObject {
    bounding_box: BoundingRect<T::Scalar>,
    children: Box<Vec<RTreeNode<T>>>,
    depth: usize,
    options: Arc<RTreeOptions>,
}

enum RTreeNode<T> where T: SpatialObject {
    Leaf(T),
    DirectoryNode(DirectoryNodeData<T>),
}

pub struct RTree<T> where T: SpatialObject {
    root: DirectoryNodeData<T>,
    size: usize,
}

impl<'a, T> RTree<T> where T: SpatialObject {
    pub fn new() -> RTree<T> {
        RTree::new_with_options(Default::default())
    }

    pub fn new_with_options(options: RTreeOptions) -> RTree<T> {
        let options = Arc::new(options);
        RTree {
            root: DirectoryNodeData::new(1, options),
            size: 0,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn iter(&self) -> RTreeIterator<T> {
        RTreeIterator::new(&self.root)
    }

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

    pub fn lookup_and_remove(&mut self, point: &Vector2<T::Scalar>) -> Option<T> {
        let result = self.root.lookup_and_remove(point);
        if result.is_some() {
            self.size -= 1;
        }
        result
    }

    pub fn lookup(&self, point: &Vector2<T::Scalar>) -> Option<&T> {
        self.root.lookup(point)
    }
    
    pub fn nearest_neighbor(&self, point: &Vector2<T::Scalar>) -> Option<&T> {
        self.root.nearest_neighbor(point, Float::infinity())
    }
}

impl <'a, T> RTree<T> where T: SpatialObject + PartialEq {
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
    use cgmath::{Vector2, EuclideanVector};
    use num::Float;
    use rand::{XorShiftRng, SeedableRng};
    use rand::distributions::{Range, IndependentSample};

    fn random_points(size: u32) -> Vec<Vector2<f64>> {
        random_points_with_seed(size, [2, 31, 244, 22])
    }

    fn random_points_with_seed(size: u32, seed: [u32; 4]) -> Vec<Vector2<f64>> {
        const SIZE: f32 = 1000.;
        let mut rng = XorShiftRng::from_seed(seed);
        let range = Range::new(-SIZE as f64 / 2., SIZE as f64 / 2.);
        let mut points = Vec::new();
        for _ in 0 .. size {
            let x = range.ind_sample(&mut rng);
            let y = range.ind_sample(&mut rng);
            points.push(Vector2::new(x, y));
        }
        points
    }

    fn create_random_tree<'a>(size: u32) -> (RTree<Vector2<f64>>, Vec<Vector2<f64>>) {
        let mut tree: RTree<Vector2<f64>> = RTree::new();
        
        let points = random_points(size);
        for point in points.iter() {
            tree.insert(point.clone())
        }
        (tree, points)
    }
        
    #[test]
    fn test_nearest_neighbor() {
        let (tree, points) = create_random_tree(1000);
        let sample_points = random_points(100);
        for sample_point in &sample_points {
            let mut nearest = None;
            let mut closest_dist = Float::infinity();
            for point in &points {
                let new_dist = (point - sample_point).length2();
                if new_dist < closest_dist {
                    closest_dist = new_dist;
                    nearest = Some(point);
                }
            }
            assert!(nearest == tree.nearest_neighbor(&sample_point));
        }
    }

    #[test]
    fn test_lookup() {
        let (tree, points) = create_random_tree(10000);
        let sample_points = random_points_with_seed(1000, [2, 1, 0, 3]);
        for sample_point in &sample_points {
            assert!(!tree.lookup(sample_point).is_some());
        }
        for point in &points {
            assert!(tree.lookup(point) == Some(point));
        }
    }

    #[test]
    fn test_lookup_and_remove() {
        let (mut tree, points) = create_random_tree(10000);
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
        let random_points = random_points(300);
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
        let (tree, reference_points) = create_random_tree(100);
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
