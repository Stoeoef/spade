use crate::{delaunay_core::Dcel, intersection_iterator::LineIntersectionIterator};
use crate::{handles::*, intersection_iterator::Intersection};
use crate::{
    DelaunayTriangulation, HasPosition, HintGenerator, InsertionError, LastUsedVertexHintGenerator,
    Point2, Triangulation, TriangulationExt,
};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use alloc::{vec, vec::Vec};

/// Undirected edge type of a [ConstrainedDelaunayTriangulation] (CDT).
///
/// CDTs need to store if an undirected edge is a constrained edge. To do so, CDTs don't use
/// the configured undirected edge type directly but wrap it into `CdtEdge<UE>` first.
///
/// This type will only be relevant if the triangulation's undirected edge type is being
/// overwritten.
///
/// # Type parameters
/// UE: The user configurable undirected edge type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde")
)]
pub struct CdtEdge<UE>(bool, UE);

impl<UE> CdtEdge<UE> {
    /// Returns `true` if this edge is a constraint edge.
    pub fn is_constraint_edge(&self) -> bool {
        self.0
    }

    fn make_constraint_edge(&mut self) {
        assert!(!self.is_constraint_edge());
        self.0 = true;
    }

    /// Returns the wrapped undirected edge data type.
    pub fn data(&self) -> &UE {
        &self.1
    }

    /// Returns the wrapped undirected edge data type.
    pub fn data_mut(&mut self) -> &mut UE {
        &mut self.1
    }
}

impl<UE: Default> Default for CdtEdge<UE> {
    fn default() -> Self {
        CdtEdge(false, UE::default())
    }
}

impl<UE> AsRef<UE> for CdtEdge<UE> {
    fn as_ref(&self) -> &UE {
        self.data()
    }
}

impl<UE> AsMut<UE> for CdtEdge<UE> {
    fn as_mut(&mut self) -> &mut UE {
        self.data_mut()
    }
}

/// A two dimensional
/// [constrained Delaunay triangulation](https://en.wikipedia.org/wiki/Constrained_Delaunay_triangulation).
///
/// A constrained Delaunay triangulation (CDT) is a triangulation that
/// can contain _constraint edges_. These edges will always be present
/// in the resulting triangulation.
///
#[doc = include_str!("../images/cdt.svg")]
///
/// *Left: A CDT with 4 constraint edges. Right: The same triangulation
/// without constraint edges*
///
///
/// The resulting triangulation
/// does not necessarily fulfill the Delaunay property.
///
/// This implementation currently supports only _weakly intersecting_
/// constraints, thus, constraint edges are allowed to touch at
/// their start or end point but are not allowed to intersect at
/// any interior point.
///
/// The constrained triangulation shares most of the implementation of
/// the usual Delaunay triangulation, refer to `DelaunayTriangulation`
/// for more information about type parameters, iteration, performance
/// and more examples.

///
/// # Example
///
/// ```
/// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
/// # fn try_main() -> Result<(), spade::InsertionError> {
/// let mut cdt = ConstrainedDelaunayTriangulation::<Point2<_>>::new();
/// let v0 = cdt.insert(Point2::new(0f64, 0.0))?;
/// let v1 = cdt.insert(Point2::new(1.0, 0.0))?;
/// cdt.add_constraint(v0, v1);
/// // Alternatively, consider using this shorthand
/// cdt.add_constraint_edge(Point2::new(1.0, 1.0), Point2::new(1.0, 0.0))?;
/// println!("Number of constraints: {}", cdt.num_constraints()); // 2 constraints
/// // Constraints are bidirectional!
/// assert!(cdt.exists_constraint(v1, v0));
/// assert!(cdt.exists_constraint(v0, v1));
/// // Check if a new constraint could be added
/// let from = Point2::new(1.0, -2.0);
/// let to = Point2::new(1.0, 0.0);
/// if !cdt.intersects_constraint(from, to) {
///     // No intersections, the edge can be added
///     cdt.add_constraint_edge(from, to)?;
/// }
/// # Ok(()) }
/// # fn main() { try_main().unwrap() }
/// ```
///
/// # See also
/// Refer to [Triangulation] for most implemented methods on this type.
/// Refer to [DelaunayTriangulation](crate::DelaunayTriangulation) for general
/// information about using Delaunay triangulations.
#[doc(alias = "CDT")]
#[derive(Clone)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde")
)]
pub struct ConstrainedDelaunayTriangulation<
    V,
    DE = (),
    UE = (),
    F = (),
    L = LastUsedVertexHintGenerator,
> where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    s: Dcel<V, DE, CdtEdge<UE>, F>,
    num_constraints: usize,
    lookup: L,
}

impl<V, DE, UE, F, L> Default for ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    fn default() -> Self {
        ConstrainedDelaunayTriangulation {
            s: Default::default(),
            num_constraints: 0,
            lookup: Default::default(),
        }
    }
}

impl<V, DE, UE, F, L> Triangulation for ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    type Vertex = V;
    type DirectedEdge = DE;
    type UndirectedEdge = CdtEdge<UE>;
    type Face = F;
    type HintGenerator = L;

    fn s(&self) -> &Dcel<V, DE, CdtEdge<UE>, F> {
        &self.s
    }

    fn s_mut(&mut self) -> &mut Dcel<V, DE, CdtEdge<UE>, F> {
        &mut self.s
    }

    fn is_defined_legal(&self, edge: FixedUndirectedEdgeHandle) -> bool {
        self.is_constraint_edge(edge)
    }

    fn handle_legal_edge_split(&mut self, handles: [FixedUndirectedEdgeHandle; 2]) {
        self.num_constraints += 1;
        for handle in &handles {
            if !self.is_constraint_edge(*handle) {
                self.s
                    .undirected_edge_data_mut(*handle)
                    .make_constraint_edge();
            }
        }
    }

    fn hint_generator(&self) -> &Self::HintGenerator {
        &self.lookup
    }

    fn hint_generator_mut(&mut self) -> &mut Self::HintGenerator {
        &mut self.lookup
    }

    fn clear(&mut self) {
        self.num_constraints = 0;
        self.s_mut().clear();
        let new_hint_generator = HintGenerator::initialize_from_triangulation(self);
        *self.hint_generator_mut() = new_hint_generator;
    }
}

impl<V, DE, UE, F, L> From<DelaunayTriangulation<V, DE, UE, F, L>>
    for ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    fn from(value: DelaunayTriangulation<V, DE, UE, F, L>) -> Self {
        let dcel = value.dcel;
        let s = dcel.map_undirected_edges(|edge| CdtEdge(false, edge));
        let lookup = value.hint_generator;

        ConstrainedDelaunayTriangulation {
            s,
            num_constraints: 0,
            lookup,
        }
    }
}

impl<V, DE, UE, F, L> ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    /// Removes a vertex from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles.
    pub fn remove(&mut self, vertex: FixedVertexHandle) -> V {
        let num_removed_constraints = self
            .s
            .vertex(vertex)
            .out_edges()
            .map(|edge| edge.is_constraint_edge())
            .filter(|b| *b)
            .count();
        self.num_constraints -= num_removed_constraints;
        self.remove_and_notify(vertex)
    }

    /// Returns the number of constraint edges.
    pub fn num_constraints(&self) -> usize {
        self.num_constraints
    }

    /// Returns `true` if a given edge is a constraint edge.
    pub fn is_constraint_edge(&self, edge: FixedUndirectedEdgeHandle) -> bool {
        self.s.undirected_edge_data(edge).is_constraint_edge()
    }

    /// Checks if two vertices are connected by a constraint edge.
    pub fn exists_constraint(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        self.get_edge_from_neighbors(from, to)
            .map(|e| e.is_constraint_edge())
            .unwrap_or(false)
    }

    /// Checks if a constraint edge can be added.
    ///
    /// Returns `false` if the line from `from` to `to` intersects another
    /// constraint edge.
    pub fn can_add_constraint(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        let line_intersection_iterator = LineIntersectionIterator::new_from_handles(self, from, to);
        !self.contains_any_constraint_edge(line_intersection_iterator)
    }

    /// Checks if a line intersects a constraint edge.
    ///
    /// Returns `true` if the edge from `from` to `to` intersects a
    /// constraint edge.
    pub fn intersects_constraint(
        &self,
        line_from: Point2<V::Scalar>,
        line_to: Point2<V::Scalar>,
    ) -> bool {
        let line_intersection_iterator = LineIntersectionIterator::new(self, line_from, line_to);
        self.contains_any_constraint_edge(line_intersection_iterator)
    }

    fn contains_any_constraint_edge(
        &self,
        mut line_intersection_iterator: LineIntersectionIterator<V, DE, CdtEdge<UE>, F>,
    ) -> bool {
        line_intersection_iterator.any(|intersection| match intersection {
            Intersection::EdgeIntersection(edge) => edge.is_constraint_edge(),
            _ => false,
        })
    }

    /// Creates a several constraint edges by taking and connecting vertices from an iterator.
    ///
    /// Every two sequential vertices in the input iterator will be connected by a constraint edge.
    /// If `closed` is set to true, the first and last vertex will also be connected.
    ///
    /// # Special cases:
    ///  - Does nothing if input iterator is empty
    ///  - Only inserts the single vertex if the input iterator contains exactly one element
    ///
    /// # Example
    /// ```
    /// # fn main() -> Result<(), spade::InsertionError> {
    /// use spade::{ConstrainedDelaunayTriangulation, Point2};
    ///
    /// const NUM_VERTICES: usize = 51;
    ///
    /// let mut cdt = ConstrainedDelaunayTriangulation::<_>::default();
    ///
    /// // Iterates through vertices on a circle
    /// let vertices = (0..NUM_VERTICES).map(|i| {
    ///     let angle = std::f64::consts::PI * 2.0 * i as f64 / NUM_VERTICES as f64;
    ///     let (sin, cos) = angle.sin_cos();
    ///     Point2::new(sin, cos)
    /// });
    ///
    /// cdt.add_constraint_edges(vertices, true)?;
    /// # Ok(()) }
    /// ```
    ///
    /// # Panics
    /// Panics if any of the generated constraints intersects with any other constraint edge.
    pub fn add_constraint_edges(
        &mut self,
        vertices: impl IntoIterator<Item = V>,
        closed: bool,
    ) -> Result<(), InsertionError> {
        let mut iter = vertices.into_iter();
        if let Some(first) = iter.next() {
            let first_handle = self.insert(first)?;
            let mut previous_handle = first_handle;
            let mut current_handle = first_handle;
            for current in iter {
                current_handle = self.insert(current)?;
                self.add_constraint(previous_handle, current_handle);
                previous_handle = current_handle;
            }

            if closed && current_handle != first_handle {
                self.add_constraint(current_handle, first_handle);
            }
        }

        Ok(())
    }

    /// Insert two points and creates a constraint between them.
    ///
    /// Returns `true` if at least one constraint edge was added.
    ///
    /// # Panics
    /// Panics if the new constraint edge intersects with an existing
    /// constraint edge. Use [can_add_constraint](Self::can_add_constraint) to check.
    pub fn add_constraint_edge(&mut self, from: V, to: V) -> Result<bool, InsertionError> {
        let from_handle = self.insert(from)?;
        let to_handle = self.insert(to)?;
        Ok(self.add_constraint(from_handle, to_handle))
    }

    /// Adds a constraint edge between to vertices.
    ///
    /// Returns `true` if at least one constraint edge was added.
    /// Note that the given constraint might be splitted into smaller edges
    /// if a vertex in the triangulation lies exactly on the constraint edge.
    /// Thus, `cdt.exists_constraint(from, to)` is not necessarily `true`
    /// after a call to this function.
    ///
    /// Returns false and does nothing if `from == to`.
    ///
    /// # Panics
    /// Panics if the new constraint edge intersects an existing
    /// constraint edge.
    pub fn add_constraint(&mut self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        use super::intersection_iterator::{
            trace_direction_out_of_edge, trace_direction_out_of_vertex, EdgeOutDirection,
            VertexOutDirection,
        };

        if from == to {
            return false;
        }

        let mut cur_from = from;
        let line_from = self.s().vertex(cur_from).position();
        let line_to = self.vertex(to).position();
        let mut result = false;

        'conflict_regions_loop: loop {
            if cur_from == to {
                return result;
            }

            let first_intersection =
                trace_direction_out_of_vertex(self.s().vertex(cur_from), line_to);

            let first_edge = match first_intersection {
                // Panic: A constraint edge should never be able to cross the convex hull
                VertexOutDirection::ConvexHull => panic!("Should not be reachable. This is a bug."),
                VertexOutDirection::EdgeOverlap(edge) => {
                    cur_from = edge.to().fix();
                    let edge = edge.fix().as_undirected();
                    result |= self.make_constraint_edge(edge);
                    continue;
                }
                VertexOutDirection::EdgeIntersection(edge) => edge,
            };

            let mut border_loop = alloc::collections::VecDeque::new();

            border_loop.push_back(first_edge.rev().next().fix());
            border_loop.push_front(first_edge.rev().prev().fix());

            let mut required_rotation = 2;

            let mut current_edge = first_edge;

            let mut faces_to_remove: Vec<FixedFaceHandle<InnerTag>> =
                vec![first_edge.rev().face().as_inner().unwrap().fix()];
            let mut edges_to_remove: Vec<FixedUndirectedEdgeHandle> =
                vec![first_edge.as_undirected().fix()];

            let check_for_constraint_intersection = |edge| {
                assert!(
                    !self.is_constraint_edge(edge),
                    "Error - constraint edges must not intersect each other"
                );
            };

            check_for_constraint_intersection(current_edge.as_undirected().fix());

            loop {
                match trace_direction_out_of_edge(current_edge, line_from, line_to) {
                    EdgeOutDirection::NoIntersection | EdgeOutDirection::ConvexHull => {
                        panic!("Constraint edge does not end in a vertex. This is a bug.")
                    }
                    EdgeOutDirection::EdgeIntersection(edge) => {
                        check_for_constraint_intersection(edge.as_undirected().fix());

                        let next = edge.rev().next();

                        if next == current_edge {
                            border_loop.push_front(edge.rev().prev().fix());
                        } else {
                            required_rotation += 1;
                            assert_eq!(next.next(), current_edge);
                            border_loop.push_back(next.fix());
                        }
                        let edge_fixed = edge.fix();
                        edges_to_remove.push(edge_fixed.as_undirected());
                        faces_to_remove.push(edge.rev().face().as_inner().unwrap().fix());
                        current_edge = edge;
                        continue;
                    }
                    EdgeOutDirection::VertexIntersection(vertex) => {
                        let vertex_fixed = vertex.fix();
                        border_loop.push_front(current_edge.next().fix());
                        border_loop.push_back(current_edge.prev().fix());
                        faces_to_remove.push(current_edge.face().as_inner().unwrap().fix());

                        // We arrived at another vertex and need to close the constraint region
                        // We'll implicitly add an edge from cur_from to vertex that will be
                        // marked as constraint edge later

                        border_loop.rotate_right(required_rotation);
                        let border_loop_vec: Vec<_> = border_loop.into();

                        // The last edge of border_loop_vec must be part of the added constraint
                        // edge. Otherwise remesh_edge_ring will not create an edge between
                        // cur_from and vertex
                        assert_eq!(
                            self.directed_edge(*border_loop_vec.last().unwrap())
                                .from()
                                .fix(),
                            cur_from
                        );

                        // Remesh border loop
                        let mut isolation_result =
                            super::delaunay_core::dcel_operations::remesh_edge_ring(
                                self.s_mut(),
                                border_loop_vec,
                                edges_to_remove,
                                faces_to_remove,
                            );

                        // Mark constraint edge
                        let target_vertex = self.s().vertex(vertex_fixed);
                        let constraint_edge_index = isolation_result
                            .new_edges
                            .iter()
                            .copied()
                            .position(|edge| {
                                self.s()
                                    .undirected_edge(edge)
                                    .vertices()
                                    .contains(&target_vertex)
                            })
                            .unwrap();

                        let constraint_edge = isolation_result.new_edges[constraint_edge_index];
                        debug_assert!(self
                            .s()
                            .undirected_edge(constraint_edge)
                            .vertices()
                            .iter()
                            .any(|v| v.fix() == cur_from));

                        let constraint_edge = self.s().undirected_edge(constraint_edge).fix();

                        result |= self.make_constraint_edge(constraint_edge);

                        // Don't mark the new constraint edge for legalization
                        isolation_result
                            .new_edges
                            .swap_remove(constraint_edge_index);

                        // All flipped edges need to legalized.
                        while let Some(edge) = isolation_result.new_edges.pop() {
                            let edge_handle = self.directed_edge(edge.as_directed());
                            let e2 = edge_handle.prev();
                            let e4 = edge_handle.rev().prev();
                            let left = e2.from().position();
                            let right = e4.from().position();
                            let from = edge_handle.from().position();
                            let to = edge_handle.to().position();

                            if crate::delaunay_core::math::contained_in_circumference(
                                from, to, left, right,
                            ) {
                                let mut push_if_flip_candidate =
                                    |edge: FixedUndirectedEdgeHandle| {
                                        if isolation_result.is_new_edge(edge)
                                            && edge != constraint_edge
                                        {
                                            isolation_result.new_edges.push(edge);
                                        }
                                    };

                                let e1 = edge_handle.next();
                                let e3 = edge_handle.rev().next();

                                push_if_flip_candidate(e1.fix().as_undirected());
                                push_if_flip_candidate(e2.fix().as_undirected());
                                push_if_flip_candidate(e3.fix().as_undirected());
                                push_if_flip_candidate(e4.fix().as_undirected());

                                crate::delaunay_core::dcel_operations::flip_cw(self.s_mut(), edge);
                            }
                        }
                        crate::delaunay_core::dcel_operations::cleanup_isolated_vertex(
                            self.s_mut(),
                            &mut isolation_result,
                        );

                        cur_from = vertex_fixed;

                        continue 'conflict_regions_loop;
                    }
                }
            }
        }
    }

    fn make_constraint_edge(&mut self, edge: FixedUndirectedEdgeHandle) -> bool {
        if !self.is_constraint_edge(edge) {
            self.s.undirected_edge_data_mut(edge).make_constraint_edge();
            self.num_constraints += 1;
            true
        } else {
            false
        }
    }

    #[cfg(test)]
    pub fn cdt_sanity_check(&self) {
        let num_undirected_edges = self
            .s
            .undirected_edges()
            .filter(|e| e.is_constraint_edge())
            .count();

        assert_eq!(num_undirected_edges, self.num_constraints());
        self.basic_sanity_check();
    }
}

#[cfg(test)]
mod test {
    use super::ConstrainedDelaunayTriangulation;
    use crate::test_utilities::*;
    use crate::{DelaunayTriangulation, InsertionError, Point2, Triangulation};
    use rand::distributions::{Distribution, Uniform};
    use rand::{Rng, SeedableRng};

    type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;
    type Delaunay = DelaunayTriangulation<Point2<f64>>;

    #[test]
    fn test_add_same_from_and_to_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(2.0, 2.0))?;
        cdt.insert(Point2::new(1.0, 2.0))?;

        assert!(!cdt.add_constraint(v0, v0));

        let new_point = Point2::new(3.1, 2.0);
        assert!(!cdt.add_constraint_edge(new_point, new_point)?);

        assert_eq!(0, cdt.num_constraints());
        assert_eq!(4, cdt.num_vertices());

        cdt.cdt_sanity_check();

        Ok(())
    }

    use alloc::{vec, vec::Vec};

    #[test]
    fn test_add_single_simple_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(2.0, 2.0))?;
        let v2 = cdt.insert(Point2::new(1.0, 0.5))?;
        let v3 = cdt.insert(Point2::new(0.5, 1.0))?;
        assert!(cdt.get_edge_from_neighbors(v0, v1).is_none());
        assert!(cdt.get_edge_from_neighbors(v2, v3).is_some());

        assert!(cdt.add_constraint(v1, v0));
        assert!(!cdt.add_constraint(v0, v1));
        let edge = cdt
            .get_edge_from_neighbors(v0, v1)
            .expect("Expected constraint edge")
            .as_undirected()
            .fix();
        assert!(cdt.get_edge_from_neighbors(v2, v3).is_none());
        assert!(cdt.is_constraint_edge(edge));
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_existing_edge_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(2.0, 2.0))?;
        let v2 = cdt.insert(Point2::new(1.0, 0.0))?;
        assert!(cdt.add_constraint(v0, v1));
        assert!(cdt.add_constraint(v0, v2));
        assert!(cdt.add_constraint(v1, v2));
        for edge in cdt.undirected_edges() {
            assert!(cdt.is_constraint_edge(edge.fix()));
        }
        assert!(!cdt.add_constraint(v1, v0));
        assert!(!cdt.add_constraint(v1, v2));
        assert_eq!(cdt.num_constraints, 3);
        Ok(())
    }

    #[test]
    fn test_mid_overlapping_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.5))?;
        let v1 = cdt.insert(Point2::new(2.0, 0.5))?;
        let v2 = cdt.insert(Point2::new(3.0, 0.5))?;
        let v3 = cdt.insert(Point2::new(5.0, 0.5))?;
        cdt.insert(Point2::new(1.0, 1.0))?;
        cdt.insert(Point2::new(1.0, 0.0))?;
        cdt.insert(Point2::new(3.0, 1.0))?;
        cdt.insert(Point2::new(3.0, 0.0))?;
        assert!(cdt.get_edge_from_neighbors(v1, v2).is_some());
        let mut copy = cdt.clone();
        assert!(cdt.add_constraint(v0, v3));
        assert_eq!(cdt.num_constraints(), 3);

        copy.add_constraint(v2, v3);
        assert_eq!(copy.num_constraints(), 1);
        copy.add_constraint(v0, v3);
        assert_eq!(copy.num_constraints(), 3);
        Ok(())
    }

    #[test]
    fn test_add_single_complex_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(1.0, 0.0))?;
        cdt.insert(Point2::new(0.0, 1.0))?;
        cdt.insert(Point2::new(2.0, 1.0))?;

        let v1 = cdt.insert(Point2::new(2.0, 2.0))?;
        assert!(cdt.get_edge_from_neighbors(v0, v1).is_none());
        cdt.add_constraint(v0, v1);
        cdt.cdt_sanity_check();
        let edge = cdt
            .get_edge_from_neighbors(v0, v1)
            .expect("Expected constraint edge")
            .fix()
            .as_undirected();
        assert!(cdt.is_constraint_edge(edge));
        Ok(())
    }

    #[test]
    fn test_add_single_constraint() -> Result<(), InsertionError> {
        let points = random_points_with_seed(1000, SEED);
        let mut cdt = Cdt::new();
        assert_eq!(cdt.num_constraints(), 0);
        let mut handles = Vec::new();
        cdt.cdt_sanity_check();
        for point in points.into_iter() {
            handles.push(cdt.insert(point)?);
        }
        cdt.add_constraint(handles[40], handles[200]);
        assert_eq!(cdt.num_constraints(), 1);
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_add_border_constraint() -> Result<(), InsertionError> {
        let points = random_points_with_seed(1000, SEED);
        let mut cdt = Cdt::new();
        let mut max_y = -::core::f64::MAX;
        for point in points {
            max_y = max_y.max(point.y);
            cdt.insert(point)?;
        }
        let v0 = cdt.insert(Point2::new(-20., max_y + 10.))?;
        let v1 = cdt.insert(Point2::new(20., max_y + 10.))?;
        cdt.add_constraint(v0, v1);
        assert_eq!(cdt.num_constraints(), 1);
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_add_multiple_constraints_overlapping() -> Result<(), InsertionError> {
        test_add_multiple_constraints(true)
    }

    #[test]
    fn test_add_multiple_constraints_non_overlapping() -> Result<(), InsertionError> {
        test_add_multiple_constraints(false)
    }

    fn test_add_multiple_constraints(overlapping: bool) -> Result<(), InsertionError> {
        const RANGE: f64 = 10.;
        let seed = if overlapping { SEED } else { SEED2 };
        let points = random_points_in_range(RANGE, 1000, seed);
        let mut cdt = Cdt::new();
        for point in points {
            cdt.insert(point)?;
        }
        let seed = if overlapping { SEED } else { SEED2 };
        let delaunay_points = random_points_in_range(RANGE * 0.9, 80, seed);
        // Use a delaunay triangulation to "generate" non intersecting constraint edges
        let mut d = Delaunay::new();
        for p in delaunay_points {
            d.insert(p)?;
        }
        let mut used_vertices = ::hashbrown::HashSet::new();

        let mut inserted_constraints = Vec::new();
        for v in d.vertices() {
            // Insert only edges that do not touch at the end points if
            // overlapping is false
            if overlapping || used_vertices.insert(v.fix()) {
                let out_edge = v.out_edge().unwrap();
                let to = out_edge.to();

                used_vertices.insert(to.fix());

                let h0 = cdt.insert(v.position())?;
                let h1 = cdt.insert(to.position())?;

                if cdt.add_constraint(h0, h1) {
                    inserted_constraints.push((h0, h1));
                }
                cdt.cdt_sanity_check();

                assert_eq!(cdt.num_constraints(), inserted_constraints.len());
            }
        }
        // Check if all constraints still exists
        for (from, to) in inserted_constraints {
            assert!(cdt.exists_constraint(from, to));
        }
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn crash_case() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        cdt.insert(Point2::new(-8.403036273981348, -0.2248814041797189))?;
        cdt.insert(Point2::new(-8.353215494321136, 0.6088667888877364))?;
        cdt.insert(Point2::new(-7.811923439447166, -0.20003314976217013))?;
        cdt.insert(Point2::new(-7.710431174668773, 0.40691184742787456))?;

        let v0 = cdt.insert(Point2::new(-8.907731924022768, 1.7433952434737847))?;
        let v1 = cdt.insert(Point2::new(-7.899415172394501, -1.4867902598716558))?;
        cdt.cdt_sanity_check();
        cdt.add_constraint(v0, v1);
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_split_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(1.0, 0.0))?;
        cdt.insert(Point2::new(0.0, 1.0))?;
        let v0 = cdt.insert(Point2::new(0.0, 0.5))?;
        let v_last = cdt.insert(Point2::new(1.0, 0.5))?;
        cdt.add_constraint(v0, v_last);
        assert_eq!(cdt.num_constraints(), 1);
        // These points split an existing constraint
        let v1 = cdt.insert(Point2::new(0.25, 0.5))?;
        assert_eq!(cdt.num_constraints(), 2);
        let v2 = cdt.insert(Point2::new(0.75, 0.5))?;
        assert_eq!(cdt.num_constraints(), 3);
        assert!(cdt.exists_constraint(v0, v1));
        assert!(cdt.exists_constraint(v1, v2));
        assert!(cdt.exists_constraint(v2, v_last));
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_simple_retriangulation() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(1.0, 0.25))?;
        cdt.insert(Point2::new(1.0, -0.25))?;
        let v3 = cdt.insert(Point2::new(2.0, 0.75))?;
        let v4 = cdt.insert(Point2::new(2.5, -0.3))?;
        cdt.insert(Point2::new(2.75, 0.75))?;
        cdt.insert(Point2::new(3.0, 0.75))?;
        cdt.insert(Point2::new(4.0, 0.25))?;
        cdt.insert(Point2::new(4.0, -0.25))?;
        let v7 = cdt.insert(Point2::new(5.0, 0.0))?;
        assert!(cdt.get_edge_from_neighbors(v3, v4).is_some());
        cdt.add_constraint(v0, v7);
        assert!(cdt.get_edge_from_neighbors(v0, v7).is_some());
        assert!(cdt.get_edge_from_neighbors(v3, v4).is_none());

        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_add_constraint_over_point() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(1.0, 0.0))?;
        let v2 = cdt.insert(Point2::new(2.0, 0.0))?;
        cdt.insert(Point2::new(0.0, 1.0))?;
        cdt.add_constraint(v0, v2);
        assert_eq!(cdt.num_constraints(), 2);
        assert!(cdt.exists_constraint(v0, v1));
        assert!(cdt.exists_constraint(v1, v2));
        cdt.cdt_sanity_check();
        Ok(())
    }

    fn test_cdt() -> Result<Cdt, InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(1.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(0.0, 1.0))?;
        cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(1.0, 1.0))?;
        cdt.add_constraint(v0, v1);
        Ok(cdt)
    }

    #[test]
    fn test_check_intersects_constraint_edge() -> Result<(), InsertionError> {
        let cdt = test_cdt()?;
        let from = Point2::new(0.2, 0.2);
        let to = Point2::new(0.6, 0.7);
        assert!(cdt.intersects_constraint(from, to));
        assert!(cdt.intersects_constraint(to, from));
        let to = Point2::new(-0.5, 0.2);
        assert!(!cdt.intersects_constraint(from, to));
        let from = Point2::new(0.5, 0.5);
        assert!(cdt.intersects_constraint(from, to));
        assert!(cdt.intersects_constraint(to, from));
        Ok(())
    }

    #[test]
    fn test_add_constraint_degenerate() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(0.0, 1.0))?;
        assert!(cdt.add_constraint(v0, v1));
        assert!(!cdt.add_constraint(v1, v0));
        assert_eq!(cdt.num_constraints(), 1);
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(0.0, 2.0))?;
        cdt.insert(Point2::new(0.0, 1.0))?;
        assert!(cdt.add_constraint(v0, v1));
        assert_eq!(cdt.num_constraints(), 2);
        Ok(())
    }

    fn random_points_on_line<R>(
        range: i64,
        num_points: usize,
        rng: &mut R,
        line_dir: Point2<f64>,
    ) -> Vec<Point2<f64>>
    where
        R: Rng,
    {
        let mut result = Vec::with_capacity(num_points);
        let range = Uniform::new(-range, range);
        for _ in 0..num_points {
            let factor = range.sample(rng);
            result.push(line_dir.mul(factor as f64));
        }
        result
    }

    #[test]
    fn fuzz_test_on_line() -> Result<(), InsertionError> {
        // Generates points on a single line and randomly connects
        // them with constraints.
        let seed = SEED;
        const RANGE: i64 = 10000;
        const NUM_POINTS: usize = 1000;
        let mut rng = rand::rngs::StdRng::from_seed(*seed);
        let points = random_points_on_line(RANGE, NUM_POINTS, &mut rng, Point2::new(1.0, 0.0));
        let mut cdt = ConstrainedDelaunayTriangulation::<_>::new();
        for ps in points.chunks(2) {
            let from = ps[0];
            let to = ps[1];
            let from = cdt.insert(from)?;
            let to = cdt.insert(to)?;
            let should_add_constraint: bool = rng.gen();
            if from != to && should_add_constraint {
                cdt.add_constraint(from, to);
            }

            cdt.cdt_sanity_check();
        }
        Ok(())
    }

    #[test]
    fn fuzz_test_on_grid() -> Result<(), InsertionError> {
        use rand::seq::SliceRandom;
        // Generates points on a grid and randomly connects
        // them with non intersecting constraints
        let seed = SEED;
        let mut points = Vec::with_capacity((RANGE * RANGE) as usize);
        const RANGE: i64 = 30;
        const NUM_CONSTRAINTS: usize = 2000;
        for x in -RANGE..RANGE {
            for y in -RANGE..RANGE {
                points.push(Point2::new(x as f64, y as f64));
            }
        }
        let mut rng = rand::rngs::StdRng::from_seed(*seed);
        points.shuffle(&mut rng);
        let mut cdt = Cdt::new();
        for p in points {
            cdt.insert(p)?;
        }
        let range = Uniform::new(-RANGE, RANGE);
        let directions_and_offset = [
            (Point2::new(1.0, 0.0), Point2::new(0.0, 1.0)),
            (Point2::new(0.0, 1.0), Point2::new(1.0, 0.0)),
            (Point2::new(1.0, 1.0), Point2::new(0.0, 0.0)),
        ];
        for _ in 0..NUM_CONSTRAINTS {
            let &(direction, offset) = directions_and_offset.choose(&mut rng).unwrap();
            let factor1 = range.sample(&mut rng);
            let factor2 = range.sample(&mut rng);
            let p1 = offset.add(direction.mul(factor1 as f64));
            let p2 = offset.add(direction.mul(factor2 as f64));
            if p1 != p2 {
                cdt.add_constraint_edge(p1, p2)?;
            }
        }
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_panic_when_intersecting_a_constraint_edge() {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0)).unwrap();
        let v1 = cdt.insert(Point2::new(1.0, 0.0)).unwrap();
        cdt.add_constraint(v0, v1);
        cdt.add_constraint(v0, v1);
        cdt.add_constraint_edge(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0))
            .unwrap();
        cdt.add_constraint_edge(Point2::new(0.5, 0.5), Point2::new(0.5, -0.5))
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn test_panic_when_intersecting_a_complex_constraint_edge() {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.5, 2.0)).unwrap();
        cdt.insert(Point2::new(0.0, 1.5)).unwrap();
        cdt.insert(Point2::new(1.0, 1.5)).unwrap();
        cdt.add_constraint_edge(Point2::new(0.0, 0.5), Point2::new(1.0, 0.5))
            .unwrap();
        let v1 = cdt.insert(Point2::new(0.5, 0.0)).unwrap();

        cdt.add_constraint(v0, v1);
    }

    #[test]
    fn test_cdt_remove_degenerate() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        let v1 = cdt.insert(Point2::new(1.0, 0.0))?;
        let v2 = cdt.insert(Point2::new(0.0, 1.0))?;
        cdt.add_constraint(v0, v1);
        cdt.add_constraint(v1, v2);
        cdt.add_constraint(v2, v0);
        assert_eq!(cdt.num_constraints(), 3);
        cdt.remove(v1);
        assert_eq!(cdt.num_constraints(), 1);
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_crash_scenario() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        for point in get_points().iter().cloned() {
            cdt.insert(point)?;
        }

        let from = cdt.insert(Point2::new(3.2348222581121586, -8.136734693290444))?;
        cdt.cdt_sanity_check();
        let to = cdt.insert(Point2::new(-8.839844309691154, -8.930685085211245))?;
        cdt.cdt_sanity_check();

        cdt.add_constraint(from, to);
        cdt.cdt_sanity_check();

        Ok(())
    }

    fn get_points() -> Vec<Point2<f64>> {
        vec![
            Point2::new(-3.947938514986289, -8.016680534876258),
            Point2::new(-4.016029045366132, -9.680855465455608),
            Point2::new(-4.46653326962287, -8.462568264351527),
            Point2::new(-7.033691993749462, -8.88072731817851),
            Point2::new(-6.058360215097096, -8.644637388990939),
        ]
    }

    #[test]
    fn test_add_constraint_edges() -> Result<(), InsertionError> {
        for is_closed in [true, false] {
            let mut cdt = Cdt::new();

            const NUM_VERTICES: usize = 51;
            let vertices = (0..NUM_VERTICES).map(|i| {
                let angle = std::f64::consts::PI * 2.0 * i as f64 / NUM_VERTICES as f64;
                let (sin, cos) = angle.sin_cos();
                Point2::new(sin, cos)
            });

            cdt.add_constraint_edges(vertices, is_closed)?;

            if is_closed {
                assert_eq!(NUM_VERTICES, cdt.num_constraints());
            } else {
                assert_eq!(NUM_VERTICES - 1, cdt.num_constraints());
            }

            cdt.cdt_sanity_check();
        }

        Ok(())
    }

    #[test]
    fn test_add_constraint_edges_empty() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();

        cdt.add_constraint_edges(std::iter::empty(), false)?;
        cdt.add_constraint_edges(std::iter::empty(), true)?;

        assert_eq!(cdt.num_vertices(), 0);
        assert_eq!(cdt.num_constraints(), 0);

        Ok(())
    }

    #[test]
    fn test_add_constraint_edges_single() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();

        cdt.add_constraint_edges([Point2::new(1.0, 1.0)], true)?;
        cdt.add_constraint_edges([Point2::new(2.0, 3.0)], false)?;

        assert_eq!(cdt.num_vertices(), 2);
        assert_eq!(cdt.num_constraints(), 0);

        Ok(())
    }

    #[test]
    fn test_add_constraint_edges_duplicate() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let point = Point2::new(0.0, 1.0);
        cdt.add_constraint_edges([point, point], true)?;
        cdt.add_constraint_edges([point, point], false)?;
        cdt.add_constraint_edges([point, point, point], true)?;
        cdt.add_constraint_edges([point, point, point], false)?;

        assert_eq!(cdt.num_vertices(), 1);
        assert_eq!(cdt.num_constraints(), 0);

        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_clear() -> Result<(), InsertionError> {
        let mut cdt = test_cdt()?;
        cdt.clear();

        assert_eq!(cdt.num_constraints(), 0);
        assert_eq!(cdt.num_all_faces(), 1);
        assert_eq!(cdt.num_vertices(), 0);
        assert_eq!(cdt.num_directed_edges(), 0);
        Ok(())
    }
}
