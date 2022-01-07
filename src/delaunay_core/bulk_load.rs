use std::cmp::Ordering;

use crate::{HasPosition, InsertionError, Point2, Triangulation, TriangulationExt};

use super::{dcel_operations, FixedDirectedEdgeHandle, FixedUndirectedEdgeHandle};

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
struct FloatOrd(f64);

impl Eq for FloatOrd {}

pub fn bulk_load<V, T>(elements: Vec<V>) -> Result<T, InsertionError>
where
    V: HasPosition,
    T: Triangulation<Vertex = V>,
{
    if elements.is_empty() {
        return Ok(T::new());
    }

    let mut min = elements[0].position();
    let mut max = elements[0].position();

    for element in &elements {
        crate::validate_vertex(element)?;
        let position = element.position();

        min = min.min(position);
        max = max.max(position);
    }

    let center = min.add(max).mul(0.5f32.into()).to_f64();

    let sorting_fn = |a: &V, b: &V| -> Ordering {
        let (a, b) = (a.position().to_f64(), b.position().to_f64());
        center
            .distance2(b)
            .partial_cmp(&center.distance2(a.position()))
            .expect("Invalid point position - possible NAN detected")
    };

    let mut result = T::with_capacity(elements.len(), elements.len() * 3, elements.len() * 2);

    // Sort by distance, smallest values last. This allows to pop values depending on their distance.
    let mut heap = binary_heap_plus::BinaryHeap::from_vec_cmp(elements, sorting_fn);

    while let Some(next) = heap.pop() {
        result.insert(next)?;
        if !result.all_vertices_on_line() && result.num_vertices() >= 4 {
            // We'll need 4 vertices to calculate a center position with good precision.
            // Otherwise, dividing by 3.0 can introduce precision loss and errors.
            break;
        }
    }

    if heap.is_empty() {
        return Ok(result);
    }

    // Get new center that is guaranteed to be within the convex hull
    let center_positions = || result.vertices().take(4).map(|v| v.position().to_f64());

    let sum_x = center_positions().map(|p| p.x).sum();
    let sum_y = center_positions().map(|p| p.y).sum();
    let center = Point2::new(sum_x, sum_y).mul(0.25);

    // Sort new elements in increasing distance compared to the new center
    let mut elements = heap.into_vec();
    elements.sort_unstable_by(sorting_fn);

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

    result.legalize_edge(edge);

    let mut current_edge = ccw_walk_start;
    loop {
        let handle = result.directed_edge(current_edge);
        let prev = handle.prev();
        let handle = handle.fix();

        let point_projection =
            super::math::project_point(next_position, prev.to().position(), prev.from().position());
        current_edge = prev.fix();

        if !point_projection.is_after_edge() && prev.side_query(next_position).is_on_left_side() {
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
    loop {
        let handle = result.directed_edge(current_edge);
        let next = handle.next();
        let handle = handle.fix();

        let point_projection =
            super::math::project_point(next_position, next.from().position(), next.to().position());

        let next_fix = next.fix();
        if !point_projection.is_after_edge() && next.side_query(next_position).is_on_left_side() {
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

    // Fix hull
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
    /// Pseudo-angle of the point
    angle: FloatOrd,

    /// An edge leaving at this hull entry.
    edge: FixedDirectedEdgeHandle,

    /// Neighbors (indexes into the hull)
    left: usize,
    right: usize,
}

/// The Hull stores a set of points which form a left-to-right order
///
/// Each point is associated with an EdgeIndex into a half-edge data structure,
/// but the Hull does not concern itself with such things.
///
/// The Hull supports one kind of lookup: for a point P, find the point Q with
/// the highest X value that is below P.  When projecting P towards the
/// sweep line, it will intersect the edge beginning at Q; this edge is the one
/// which should be split.
///
/// In addition, the Hull stores a random-access map from PointIndex to
/// HullIndex (if present), for fast lookups without hash traversal.
#[derive(Debug)]
pub struct Hull {
    buckets: Vec<usize>,
    data: Vec<Node>,

    /// Empty elements which can be reused later.
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

            if angle_from == angle_to {
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
            .find(|(index, _)| !self.empty.contains(&*index))
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

/// Returns a pseudo-angle in the 0-1 range, without expensive trig functions
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
