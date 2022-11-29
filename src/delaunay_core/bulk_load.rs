use std::cmp::{Ordering, Reverse};

use crate::{HasPosition, InsertionError, Point2, Triangulation, TriangulationExt};

use super::{dcel_operations, FixedDirectedEdgeHandle, FixedUndirectedEdgeHandle};

/// An `f64` wrapper implementing `Ord` and `Eq`.
///
/// This is only used as part of bulk loading.
/// All input coordinates are checked with `validate_coordinate` before they are used, hence
/// `Ord` and `Eq` should always be well defined.
#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
struct FloatOrd(f64);

#[allow(clippy::derive_ord_xor_partial_ord)]
impl Ord for FloatOrd {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Eq for FloatOrd {}

/// Implements a circle-sweep bulk loading algorithm for efficient initialization of Delaunay
/// triangulations.
///
/// The algorithm is motivated by:
///
/// A faster circle-sweep Delaunay triangulation algorithm
/// Ahmad Biniaz, Gholamhossein Dastghaibyfard
/// Advances in Engineering Software,
/// Volume 43, Issue 1,
/// 2012,
/// https://doi.org/10.1016/j.advengsoft.2011.09.003
///
/// Or alternatively: http://cglab.ca/~biniaz/papers/Sweep%20Circle.pdf
///
/// # Overview
///
/// The major reason for the algorithm's good performance lies in an efficient lookup structure
/// for finding *hull edges* at a certain *angle*.
/// "angle" always refers to the angle of a vertex to a center point which is calculated first.
/// The lookup structure is implemented by the `Hull` struct. It has a `get` and `insert` method
/// which can quickly find and update the edges of the hull at a given angle.
///
/// The algorithm is roughly compromised of these steps:
///
///  1. Calculate the median position of all vertices. We call this position `initial_center`.
///  2. Sort all vertices along their distance to this center.
///  3. Build a seed triangulation by inserting vertices (beginning with the closest vertex) into an
///     empty triangulation. Stop once the triangulation has at least one inner face.
///  4. Calculate the final center. The final center is some point inside the seed triangulation (e.g.
///     the average its vertices)
///  5. Initiate the `Hull` lookup structure with the seed triangulation.
///  6. Insert all remaining vertices, beginning with the vertex closest to `initial_center`.
///     This can be done efficiently as the edge "closest" to the new vertex can be identified quickly
///     with `Hull.get`. After each insertion, the hull is partially patched to be more convex
///  7. After all vertices have been inserted: The hull is not necessarily convex. Fill any "hole"
///     in the hull by a process comparable to the graham scan algorithm.
///
/// # Some details
///
/// "angle" does not refer to an actual angle in radians but rather to an approximation that doesn't
/// require trigonometry for calculation. See method `pseudo_angle` for more information.
///
/// In rare cases, step 6 is not able to insert a vertex properly. It will be skipped and inserted
/// regularly at the end (slow path). This may happen especially for very skewed triangulations
/// and might be a good point for investigation if some point sets takes surprisingly long to load.
pub fn bulk_load<V, T>(mut elements: Vec<V>) -> Result<T, InsertionError>
where
    V: HasPosition,
    T: Triangulation<Vertex = V>,
{
    if elements.is_empty() {
        return Ok(T::new());
    }

    let mut point_sum = Point2::<f64>::new(0.0, 0.0);

    for element in &elements {
        crate::validate_vertex(element)?;
        let position = element.position();

        point_sum = point_sum.add(position.to_f64());
    }

    // Set the initial center to the average of all positions. This should be a good choice for most triangulations.
    //
    // The research paper uses a different approach by taking the center of the points' bounding box.
    // However, this position might be far off the center off mass if the triangulation has just a few outliers.
    // This could lead to a very uneven angle distribution as nearly all points are might be in a very small angle segment
    // around the center. This degrades the hull-structure's lookup and insertion performance.
    // For this reason, taking the average appears to be a safer option as most vertices should be distributed around the
    // initial center.
    let initial_center = point_sum.mul(1.0 / (elements.len() as f64));

    let mut result = T::with_capacity(elements.len(), elements.len() * 3, elements.len() * 2);

    // Sort by distance, smallest values last. This allows to pop values depending on their distance.
    elements.sort_unstable_by_key(|e| {
        Reverse(FloatOrd(initial_center.distance_2(e.position().to_f64())))
    });

    while let Some(next) = elements.pop() {
        result.insert(next)?;
        if !result.all_vertices_on_line() && result.num_vertices() >= 4 {
            // We'll need 4 vertices to calculate a center position with good precision.
            // Otherwise, dividing by 3.0 can introduce precision loss and errors.
            break;
        }
    }

    if elements.is_empty() {
        return Ok(result);
    }

    // Get new center that is guaranteed to be within the convex hull
    //
    let center_positions = || {
        result
            .vertices()
            .rev()
            .take(4)
            .map(|v| v.position().to_f64())
    };

    let sum_x = center_positions().map(|p| p.x).sum();
    let sum_y = center_positions().map(|p| p.y).sum();

    // Note that we don't re-sort the elements according to their distance to the newest center. This doesn't seem to
    // be required for the algorithms performance, probably due to the `center` being close to `initial_center`.
    // As of now, it's a unclear how to construct point sets that result in a `center` being farther off
    // `initial center` and what the impact of this would be.
    let center = Point2::new(sum_x, sum_y).mul(0.25);

    let mut hull = loop {
        if let Some(hull) = Hull::from_triangulation(&result, center) {
            break hull;
        }

        // The hull cannot be constructed in some rare cases for very degenerate
        // triangulations. Just insert another vertex and try again. Usually hull generation should succeed eventually.
        if let Some(next) = elements.pop() {
            result.insert(next).unwrap();
        } else {
            return Ok(result);
        }
    };

    let mut buffer = Vec::new();
    let mut skipped_elements = Vec::new();
    while let Some(next) = elements.pop() {
        skipped_elements.extend(
            single_bulk_insertion_step(&mut result, center, &mut hull, next, &mut buffer).err(),
        );
    }

    if cfg!(any(fuzzing, test)) {
        hull_sanity_check(&result, &hull);
    }

    fix_convexity(&mut result);

    for element in skipped_elements {
        result.insert(element)?;
    }

    Ok(result)
}

#[inline(never)] // Prevent inlining for better profiling data
fn single_bulk_insertion_step<TR, T>(
    result: &mut TR,
    center: Point2<f64>,
    hull: &mut Hull,
    element: T,
    buffer_for_edge_legalization: &mut Vec<FixedUndirectedEdgeHandle>,
) -> Result<(), T>
where
    T: HasPosition,
    TR: Triangulation<Vertex = T>,
{
    let next_position = element.position();
    let current_angle = pseudo_angle(next_position.to_f64(), center);

    let edge_hint = hull.get(current_angle);

    let edge = result.directed_edge(edge_hint);
    if edge.side_query(next_position).is_on_right_side_or_on_line() {
        // The position is, for some reason, not on the left side of the edge. This can e.g. happen
        // if the vertices have the same angle. The safest way to include these elements appears to
        // skip them and insert them individually at the end (albeit that's very slow)
        return Err(element);
    }

    assert!(edge.is_outer_edge());

    let edge = edge.fix();

    let new_vertex =
        super::dcel_operations::create_new_face_adjacent_to_edge(result.s_mut(), edge, element);
    let ccw_walk_start = result.directed_edge(edge).prev().rev().fix();
    let cw_walk_start = result.directed_edge(edge).next().rev().fix();

    // Check if the edge that was just connected requires legalization
    result.legalize_edge(edge);

    // At this stage the new vertex was successfully inserted. However, insertions like this will end
    // up in a strongly *star shaped* triangulation instead of a nice nearly-convex blob of faces.
    //
    // To fix this, the algorithm proceeds by connecting some of the adjacent edges and forming new
    // faces. A faces is only created if all of its inner angles are less than 90 degrees. This
    // tends to be a good heuristic that doesn't create too many skewed triangles which would need
    // to be fixed later. Refer to the motivating research paper (see method `bulk_load`) for
    // more information.
    //
    // Before:
    //
    // outer face
    //
    //       v <--- the new vertex
    //      /\
    //     /  \     +---- an edge that should potentially not be adjacent to the outer face
    //    /    \    v
    //   x0----x1--------x2
    //
    // After:
    // *if* the angle between v->x1 and x1->x2 is smaller than 90°, the edge x2->v and its new
    // adjacent face is created.
    let mut current_edge = ccw_walk_start;
    loop {
        let handle = result.directed_edge(current_edge);
        let prev = handle.prev();
        let handle = handle.fix();

        let point_projection =
            super::math::project_point(next_position, prev.to().position(), prev.from().position());
        current_edge = prev.fix();

        // `!point_projection.is_after_edge` is used to identify if the new face's angle will be less
        // than 90°.
        if !point_projection.is_behind_edge() && prev.side_query(next_position).is_on_left_side() {
            let new_edge = dcel_operations::create_single_face_between_edge_and_next(
                result.s_mut(),
                current_edge,
            );

            buffer_for_edge_legalization.clear();
            buffer_for_edge_legalization.push(handle.as_undirected());
            buffer_for_edge_legalization.push(current_edge.as_undirected());
            result.legalize_edges_after_removal(buffer_for_edge_legalization, |_| false);

            current_edge = new_edge;
        } else {
            break;
        }
    }

    let mut current_edge = cw_walk_start;
    // Same as before: Create faces if they will have inner angles less than 90 degrees. This loop
    // goes in the other direction (clockwise).
    loop {
        let handle = result.directed_edge(current_edge);
        let next = handle.next();
        let handle = handle.fix();

        let point_projection =
            super::math::project_point(next_position, next.from().position(), next.to().position());

        let next_fix = next.fix();
        if !point_projection.is_behind_edge() && next.side_query(next_position).is_on_left_side() {
            let new_edge = dcel_operations::create_single_face_between_edge_and_next(
                result.s_mut(),
                current_edge,
            );

            buffer_for_edge_legalization.clear();
            buffer_for_edge_legalization.push(handle.as_undirected());
            buffer_for_edge_legalization.push(next_fix.as_undirected());
            result.legalize_edges_after_removal(buffer_for_edge_legalization, |_| false);

            current_edge = new_edge;
        } else {
            break;
        }
    }

    let new_vertex = result.vertex(new_vertex);
    let outgoing_ch_edge = new_vertex.out_edges().find(|edge| edge.is_outer_edge());

    // Fix the hull
    if let Some(second_edge) = outgoing_ch_edge {
        let first_edge = second_edge.prev();

        let first_angle = pseudo_angle(first_edge.from().position().to_f64(), center);
        let second_angle = pseudo_angle(second_edge.to().position().to_f64(), center);

        hull.insert(
            first_angle,
            current_angle,
            second_angle,
            first_edge.fix(),
            second_edge.fix(),
        );
    }
    Ok(())
}

/// Makes the outer hull convex. Similar to a graham scan.
fn fix_convexity<TR>(triangulation: &mut TR)
where
    TR: Triangulation,
{
    let mut edges_to_validate = Vec::with_capacity(2);
    let mut convex_edges: Vec<FixedDirectedEdgeHandle> = Vec::with_capacity(64);

    let mut current_fixed = triangulation.outer_face().adjacent_edge().unwrap().fix();

    loop {
        let current_handle = triangulation.directed_edge(current_fixed);
        let next_handle = current_handle.next().fix();
        convex_edges.push(current_fixed);
        current_fixed = next_handle;
        while let &[.., edge1_fixed, edge2_fixed] = &*convex_edges {
            let edge1 = triangulation.directed_edge(edge1_fixed);
            let edge2 = triangulation.directed_edge(edge2_fixed);

            let target_position = edge2.to().position();
            // Check if the new edge would violate the convex hull property by turning left
            // The convex hull must only contain right turns
            if edge1.side_query(target_position).is_on_left_side() {
                // Violation detected. It is resolved by inserting a new edge
                edges_to_validate.push(edge1.fix().as_undirected());
                edges_to_validate.push(edge2.fix().as_undirected());

                let new_edge = dcel_operations::create_single_face_between_edge_and_next(
                    triangulation.s_mut(),
                    edge1_fixed,
                );

                convex_edges.pop();
                convex_edges.pop();
                convex_edges.push(new_edge);

                triangulation.legalize_edges_after_removal(&mut edges_to_validate, |_| false);
            } else {
                break;
            }
        }

        if Some(&current_fixed) == convex_edges.get(1) {
            break;
        }
    }
}

struct Segment {
    from: FloatOrd,
    to: FloatOrd,
}

impl Segment {
    fn new(from: FloatOrd, to: FloatOrd) -> Self {
        assert_ne!(from, to);
        Self { from, to }
    }

    /// Returns `true` if this segment does not contain the angle 0.0.
    ///
    /// Pseudo angles wrap back to 0.0 after a full rotation.
    fn is_non_wrapping_segment(&self) -> bool {
        self.from < self.to
    }

    fn contains_angle(&self, angle: FloatOrd) -> bool {
        if self.is_non_wrapping_segment() {
            self.from <= angle && angle < self.to
        } else {
            self.from <= angle || angle < self.to
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Node {
    /// Pseudo-angle of this hull entry
    angle: FloatOrd,

    /// An edge leaving at this hull entry.
    edge: FixedDirectedEdgeHandle,

    /// Neighbors (indexes into the hull)
    left: usize,
    right: usize,
}

/// Implements an efficient angle-to-edge lookup for edges of the hull of a triangulation.
///
/// Refer to `bulk_load` (in `bulk_load.rs`) for more background on how this structure is being used.
///
/// It implements an efficient mapping of (pseudo-)angles to edges. To do so, it stores all inserted
/// edges in a linked list backed by a vec. Finding an edge belonging to a given angle can always
/// be done by iterating through this list until the target angle is found.
/// The entries are stored in a consistent order (either clockwise or counter clockwise)
///
/// This naive sequential search is very slow as it needs to traverse half of the list on average.
/// To speed things up, the space space of valid angles (the half open interval [0, 1) )
/// is partitioned into `n` equally sized buckets.
/// For each bucket, `Hull` stores a reference to the list entry with the *biggest angle* that
/// still belongs into that bucket. A sequential search will begin at this bucket and has to traverse
/// only a few elements before finding the target angle.
/// Since the number of buckets is re-adjusted depending on the number of hull entries, this mapping
/// will now be in O(1) for reasonably evenly distributed triangulations.
#[derive(Debug)]
pub struct Hull {
    buckets: Vec<usize>,
    data: Vec<Node>,

    /// Unused indices in data which might be reclaimed later
    empty: Vec<usize>,
}

impl Hull {
    pub fn from_triangulation<T>(triangulation: &T, center: Point2<f64>) -> Option<Self>
    where
        T: Triangulation,
    {
        assert!(!triangulation.all_vertices_on_line());

        let hull_size = triangulation.convex_hull_size();
        let mut data = Vec::with_capacity(hull_size);

        let mut prev_index = hull_size - 1;

        for (current_index, edge) in triangulation.convex_hull().enumerate() {
            let angle_from = pseudo_angle(edge.from().position().to_f64(), center);
            let angle_to = pseudo_angle(edge.to().position().to_f64(), center);

            if angle_from == angle_to || angle_from.0.is_nan() || angle_to.0.is_nan() {
                // Should only be possible for very degenerate triangulations
                return None;
            }

            let next_index = (current_index + 1) % hull_size;

            data.push(Node {
                angle: angle_from,
                edge: edge.fix(),
                left: prev_index,
                right: next_index,
            });
            prev_index = current_index;
        }
        let mut result = Self {
            buckets: Vec::new(),
            data,
            empty: Vec::new(),
        };

        const INITIAL_NUMBER_OF_BUCKETS: usize = 8;
        result.initialize_buckets(INITIAL_NUMBER_OF_BUCKETS);

        Some(result)
    }

    fn initialize_buckets(&mut self, target_size: usize) {
        self.buckets.clear();
        self.buckets.reserve(target_size);

        const INVALID: usize = usize::MAX;
        self.buckets
            .extend(std::iter::repeat(INVALID).take(target_size));

        let (first_index, current_node) = self
            .data
            .iter()
            .enumerate()
            .find(|(index, _)| !self.empty.contains(index))
            .unwrap();

        let first_index = first_index;
        let mut current_index = first_index;
        let first_bucket = self.ceiled_bucket(current_node.angle);
        self.buckets[first_bucket] = current_index;

        loop {
            let current_node = self.data[current_index];
            let segment = self.segment(&current_node);
            let start_bucket = self.ceiled_bucket(segment.from);
            let end_bucket = self.ceiled_bucket(segment.to);

            self.update_bucket_segment(start_bucket, end_bucket, current_index);

            current_index = current_node.right;

            if current_index == first_index {
                break;
            }
        }
    }

    /// Updates the hull after the insertion of a vertex.
    ///
    /// This method should be called after a vertex `v` has been inserted into the outer face of the
    /// triangulation under construction.
    ///
    /// Such a vertex is guaranteed to have two outgoing edges that are adjacent to the convex hull,
    /// e.g. `e0 -> v -> e1`
    ///
    /// In this scenarios, the parameters should be set as follows:
    /// * `left_angle`: `pseudo_angle(e0.from())`
    /// * `middle_angle`: `pseudo_angle(v.position())`
    /// * `right_angle`: `pseudo_angle(e1.to())`
    /// * `left_edge`: `e0.fix()`
    /// * `right_edge`: `e1.fix()`
    ///
    /// Note that `left_angle` and `right_angle` must already be present in the hull. Otherwise,
    /// calling this method will result in an endless loop.
    fn insert(
        &mut self,
        left_angle: FloatOrd,
        middle_angle: FloatOrd,
        mut right_angle: FloatOrd,
        left_edge: FixedDirectedEdgeHandle,
        mut right_edge: FixedDirectedEdgeHandle,
    ) {
        let left_bucket = self.floored_bucket(left_angle);

        let mut left_index = self.buckets[left_bucket];

        loop {
            let current_node = self.data[left_index];
            if current_node.angle == left_angle {
                break;
            }
            left_index = current_node.right;
        }

        let mut right_index;
        if left_angle == right_angle {
            right_index = left_index;
        } else {
            right_index = self.data[left_index].right;
            loop {
                let current_node = self.data[right_index];
                if current_node.angle == right_angle {
                    break;
                }

                if cfg!(any(fuzzing, test)) {
                    assert!(!self.empty.contains(&right_index));
                }

                // Remove current_node - it is completely overlapped by the new segment
                self.empty.push(right_index);
                self.data[current_node.left].right = current_node.right;
                self.data[current_node.right].left = current_node.left;
                right_index = current_node.right;
            }
        }

        let new_index = self.get_next_index();

        if left_angle == middle_angle {
            self.empty.push(left_index);
            left_index = self.data[left_index].left;
        } else {
            self.data[left_index].edge = left_edge;
        }

        if right_angle == middle_angle {
            if left_angle != right_angle {
                self.empty.push(right_index);
            }
            right_edge = self.data[right_index].edge;
            right_index = self.data[right_index].right;

            right_angle = self.data[right_index].angle;
        }

        // Stich the vertex between left_index and right_index
        self.data[left_index].right = new_index;
        self.data[right_index].left = new_index;

        let new_node = Node {
            angle: middle_angle,
            edge: right_edge,
            left: left_index,
            right: right_index,
        };

        self.push_or_update_node(new_node, new_index);

        // Update bucket entries appropriately
        let left_bucket = self.ceiled_bucket(left_angle);
        let middle_bucket = self.ceiled_bucket(middle_angle);
        let right_bucket = self.ceiled_bucket(right_angle);

        self.update_bucket_segment(left_bucket, middle_bucket, left_index);
        self.update_bucket_segment(middle_bucket, right_bucket, new_index);

        self.adjust_bucket_size_if_necessary();
    }

    fn get_next_index(&mut self) -> usize {
        self.empty.pop().unwrap_or(self.data.len())
    }

    fn update_bucket_segment(&mut self, left_bucket: usize, right_bucket: usize, new_value: usize) {
        if left_bucket <= right_bucket {
            for current_bucket in &mut self.buckets[left_bucket..right_bucket] {
                *current_bucket = new_value;
            }
        } else {
            // Wrap buckets
            for current_bucket in &mut self.buckets[left_bucket..] {
                *current_bucket = new_value;
            }
            for current_bucket in &mut self.buckets[..right_bucket] {
                *current_bucket = new_value;
            }
        }
    }

    fn push_or_update_node(&mut self, node: Node, index: usize) {
        if let Some(existing_node) = self.data.get_mut(index) {
            *existing_node = node;
        } else {
            assert_eq!(self.data.len(), index);
            self.data.push(node);
        }
    }

    /// Gets an edge of the hull which covers a given input angle.
    ///
    /// An edge is considered to cover an input angle if the input angle is contained in the angle
    /// segment spanned by `pseudo_angle(edge.from()) .. pseudo_angle(edge.from())`
    fn get(&self, angle: FloatOrd) -> FixedDirectedEdgeHandle {
        let mut current_handle = self.buckets[self.floored_bucket(angle)];
        loop {
            let current_node = self.data[current_handle];
            let left_angle = current_node.angle;
            let next_angle = self.data[current_node.right].angle;

            if Segment::new(left_angle, next_angle).contains_angle(angle) {
                return current_node.edge;
            }

            current_handle = current_node.right;
        }
    }

    /// Looks up what bucket a given pseudo-angle will fall into.
    fn floored_bucket(&self, angle: FloatOrd) -> usize {
        ((angle.0 * (self.buckets.len()) as f64).floor() as usize) % self.buckets.len()
    }

    fn ceiled_bucket(&self, angle: FloatOrd) -> usize {
        ((angle.0 * (self.buckets.len()) as f64).ceil() as usize) % self.buckets.len()
    }

    fn segment(&self, node: &Node) -> Segment {
        Segment::new(node.angle, self.data[node.right].angle)
    }

    fn adjust_bucket_size_if_necessary(&mut self) {
        let size = self.data.len() - self.empty.len();
        let num_buckets = self.buckets.len();

        const MIN_NUMBER_OF_BUCKETS: usize = 16;
        if num_buckets * 2 < size {
            // Too few buckets - increase bucket count
            self.initialize_buckets(num_buckets * 2);
        } else if num_buckets > size * 4 && num_buckets > MIN_NUMBER_OF_BUCKETS {
            let new_size = num_buckets / 4;
            if new_size >= MIN_NUMBER_OF_BUCKETS {
                // Too many buckets - shrink
                self.initialize_buckets(new_size);
            }
        }
    }
}

/// Returns a pseudo-angle in the 0-1 range, without expensive trigonometry functions
///
/// The angle has the following shape:
/// ```text
///              0.25
///               ^ y
///               |
///               |
///   0           |           x
///   <-----------o-----------> 0.5
///   1           |
///               |
///               |
///               v
///              0.75
/// ```
#[inline]
fn pseudo_angle(a: Point2<f64>, center: Point2<f64>) -> FloatOrd {
    let diff = a.sub(center);

    let p = diff.x / (diff.x.abs() + diff.y.abs());

    FloatOrd(1.0 - (if diff.y > 0.0 { 3.0 - p } else { 1.0 + p }) * 0.25)
}

fn hull_sanity_check(triangulation: &impl Triangulation, hull: &Hull) {
    let non_empty_nodes: Vec<_> = hull
        .data
        .iter()
        .enumerate()
        .filter(|(index, _)| !hull.empty.contains(index))
        .collect();

    for (index, node) in &non_empty_nodes {
        let left_node = hull.data[node.left];
        let right_node = hull.data[node.right];

        let edge = triangulation.directed_edge(node.edge);
        assert!(edge.is_outer_edge());

        assert!(!hull.empty.contains(&node.left));
        assert!(!hull.empty.contains(&node.right));

        assert_eq!(left_node.right, *index);
        assert_eq!(right_node.left, *index);
    }

    for (bucket_index, bucket_node) in hull.buckets.iter().enumerate() {
        assert!(!hull.empty.contains(bucket_node));
        let bucket_start_angle = FloatOrd(bucket_index as f64 / hull.buckets.len() as f64);

        for (node_index, node) in &non_empty_nodes {
            let segment = hull.segment(node);

            if segment.contains_angle(bucket_start_angle) {
                // Make sure the bucket refers to the node with the smallest angle in the same bucket
                assert_eq!(node_index, bucket_node);
            }
        }
    }
}

#[cfg(test)]
mod test {
    use float_next_after::NextAfter;
    use rand::{seq::SliceRandom, SeedableRng};

    use crate::test_utilities::{random_points_with_seed, SEED2};

    use crate::{DelaunayTriangulation, InsertionError, Point2, Triangulation, TriangulationExt};

    use super::Hull;

    #[test]
    fn test_bulk_load_with_small_number_of_vertices() -> Result<(), InsertionError> {
        for size in 0..10 {
            let triangulation =
                DelaunayTriangulation::<_>::bulk_load(random_points_with_seed(size, SEED2))?;

            assert_eq!(triangulation.num_vertices(), size);
            triangulation.sanity_check();
        }
        Ok(())
    }

    #[test]
    fn test_bulk_load_on_grid() -> Result<(), InsertionError> {
        // Inserts vertices on whole integer coordinates. This tends provokes special situations,
        // e.g. points being inserted exactly on a line.
        let mut rng = rand::rngs::StdRng::from_seed(*SEED2);
        const TEST_REPETITIONS: usize = 30;
        const GRID_SIZE: usize = 20;

        for _ in 0..TEST_REPETITIONS {
            let mut vertices = Vec::with_capacity(GRID_SIZE * GRID_SIZE);
            for x in 0..GRID_SIZE {
                for y in 0..GRID_SIZE {
                    vertices.push(Point2::new(x as f64, y as f64));
                }
            }

            vertices.shuffle(&mut rng);
            let triangulation = DelaunayTriangulation::<_>::bulk_load(vertices)?;
            assert_eq!(triangulation.num_vertices(), GRID_SIZE * GRID_SIZE);
            triangulation.sanity_check();
        }
        Ok(())
    }

    #[test]
    fn test_bulk_load_on_epsilon_grid() -> Result<(), InsertionError> {
        // Inserts vertices on a grid spaced a part by the smallest possible f64 step
        let mut rng = rand::rngs::StdRng::from_seed(*SEED2);
        const TEST_REPETITIONS: usize = 30;
        const GRID_SIZE: usize = 20;

        // Contains The first GRID_SIZE f64 values that are >= 0.0
        let mut possible_f64: Vec<_> = Vec::with_capacity(GRID_SIZE);
        let mut current_float = crate::MIN_ALLOWED_VALUE;
        for _ in 0..GRID_SIZE / 2 {
            possible_f64.push(current_float);
            possible_f64.push(-current_float);
            current_float = current_float.next_after(f64::INFINITY);
        }

        for _ in 0..TEST_REPETITIONS {
            let mut vertices = Vec::with_capacity(GRID_SIZE * GRID_SIZE);
            for x in 0..GRID_SIZE {
                for y in 0..GRID_SIZE {
                    vertices.push(Point2::new(possible_f64[x], possible_f64[y]));
                }
            }

            vertices.shuffle(&mut rng);
            let triangulation = DelaunayTriangulation::<_>::bulk_load(vertices)?;
            assert_eq!(triangulation.num_vertices(), GRID_SIZE * GRID_SIZE);
            triangulation.sanity_check();
        }
        Ok(())
    }

    #[test]
    fn test_bulk_load() -> Result<(), InsertionError> {
        const SIZE: usize = 9000;
        let mut vertices = random_points_with_seed(SIZE, SEED2);

        vertices.push(Point2::new(4.0, 4.0));
        vertices.push(Point2::new(4.0, -4.0));
        vertices.push(Point2::new(-4.0, 4.0));
        vertices.push(Point2::new(-4.0, -4.0));

        vertices.push(Point2::new(5.0, 5.0));
        vertices.push(Point2::new(5.0, -5.0));
        vertices.push(Point2::new(-5.0, 5.0));
        vertices.push(Point2::new(-5.0, -5.0));

        vertices.push(Point2::new(6.0, 6.0));
        vertices.push(Point2::new(6.0, -6.0));
        vertices.push(Point2::new(-6.0, 6.0));
        vertices.push(Point2::new(-6.0, -6.0));

        let num_vertices = vertices.len();

        let triangulation = DelaunayTriangulation::<Point2<f64>>::bulk_load(vertices)?;
        triangulation.sanity_check();
        assert_eq!(triangulation.num_vertices(), num_vertices);
        Ok(())
    }

    #[test]
    fn test_hull() -> Result<(), InsertionError> {
        let mut triangulation = DelaunayTriangulation::<_>::new();
        triangulation.insert(Point2::new(1.0, 1.0))?; // Angle: 0.375
        triangulation.insert(Point2::new(1.0, -1.0))?; // Angle: 0.125
        triangulation.insert(Point2::new(-1.0, 1.0))?; // Angle: 0.625
        triangulation.insert(Point2::new(-1.0, -1.0))?; // Angle: 0.875

        let mut hull = Hull::from_triangulation(&triangulation, Point2::new(0.0, 0.0)).unwrap();
        super::hull_sanity_check(&triangulation, &hull);

        let center = Point2::new(0.0, 0.0);
        let additional_elements = [
            Point2::new(0.4, 2.0),
            Point2::new(-0.4, 3.0),
            Point2::new(-0.4, -4.0),
            Point2::new(3.0, 5.0),
        ];

        for (index, element) in additional_elements.iter().enumerate() {
            super::single_bulk_insertion_step(
                &mut triangulation,
                center,
                &mut hull,
                *element,
                &mut Vec::new(),
            )
            .unwrap();
            if index != 0 {
                super::hull_sanity_check(&triangulation, &hull)
            }
        }
        Ok(())
    }
}
