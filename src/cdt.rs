use alloc::vec::Vec;

use num_traits::{Float, NumCast};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::cdt::ConflictRegionEnd::{EdgeOverlap, Existing};
use crate::delaunay_core::dcel_operations::flip_cw;
use crate::delaunay_core::{bulk_load_cdt, bulk_load_stable};
use crate::{
    delaunay_core::Dcel, intersection_iterator::LineIntersectionIterator, PositionInTriangulation,
    SpadeNum,
};
use crate::{handles::*, intersection_iterator::Intersection};
use crate::{
    DelaunayTriangulation, HasPosition, HintGenerator, InsertionError, LastUsedVertexHintGenerator,
    Point2, Triangulation, TriangulationExt,
};

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

    fn unmake_constraint_edge(&mut self) {
        assert!(self.is_constraint_edge());
        self.0 = false;
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

/// A two-dimensional
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
/// Refer to [DelaunayTriangulation](DelaunayTriangulation) for general
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
    dcel: Dcel<V, DE, CdtEdge<UE>, F>,
    num_constraints: usize,
    hint_generator: L,
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
            dcel: Default::default(),
            num_constraints: 0,
            hint_generator: Default::default(),
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
        &self.dcel
    }

    fn s_mut(&mut self) -> &mut Dcel<V, DE, CdtEdge<UE>, F> {
        &mut self.dcel
    }

    fn is_defined_legal(&self, edge: FixedUndirectedEdgeHandle) -> bool {
        self.is_constraint_edge(edge)
    }

    fn handle_legal_edge_split(&mut self, handles: [FixedDirectedEdgeHandle; 2]) {
        self.num_constraints += 1;
        for handle in handles.iter().map(|e| e.as_undirected()) {
            if !self.is_constraint_edge(handle) {
                self.dcel
                    .undirected_edge_data_mut(handle)
                    .make_constraint_edge();
            }
        }
    }

    fn hint_generator(&self) -> &Self::HintGenerator {
        &self.hint_generator
    }

    fn hint_generator_mut(&mut self) -> &mut Self::HintGenerator {
        &mut self.hint_generator
    }

    fn from_parts(
        dcel: Dcel<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>,
        hint_generator: Self::HintGenerator,
        num_constraints: usize,
    ) -> Self {
        Self {
            dcel,
            num_constraints,
            hint_generator,
        }
    }

    fn into_parts(
        self,
    ) -> (
        Dcel<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>,
        Self::HintGenerator,
        usize,
    ) {
        (self.dcel, self.hint_generator, self.num_constraints)
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
            dcel: s,
            num_constraints: 0,
            hint_generator: lookup,
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
    /// Efficient bulk loading of a constraint delaunay triangulation, including both vertices and constraint edges.
    ///
    /// The edges are given as pairs of vertex indices.
    ///
    /// Note that the vertex order is not preserved by this function - iterating through all vertices will not result in
    /// the same sequence as the input vertices. Use [ConstrainedDelaunayTriangulation::bulk_load_cdt_stable] for a
    /// slower but order preserving variant.
    ///
    /// Input vertices may have the same position. However, only one vertex for each position will be kept. Edges
    /// that go to a discarded vertex are rerouted and still inserted.
    /// It is arbitrary which duplicated vertex remains.
    ///
    /// # Example
    /// ```
    /// # fn main() -> Result<(), spade::InsertionError> {
    /// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    /// let mut vertices = vec![
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(3.0, -3.0),
    ///     Point2::new(-1.0, -2.0),
    ///     Point2::new(-4.0, -5.0),
    /// ];
    /// let mut edges = vec![[0, 1], [1, 2], [2, 3], [3, 4]];
    /// let cdt = ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt(vertices.clone(), edges)?;
    ///
    /// assert_eq!(cdt.num_vertices(), 5);
    /// assert_eq!(cdt.num_constraints(), 4);
    /// // The order will usually change
    /// assert_ne!(cdt.vertices().map(|v| v.position()).collect::<Vec<_>>(), vertices);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if any constraint edges overlap. Panics if the edges contain an invalid index (out of range).
    pub fn bulk_load_cdt(vertices: Vec<V>, edges: Vec<[usize; 2]>) -> Result<Self, InsertionError> {
        let mut result = bulk_load_cdt(vertices, edges)?;
        *result.hint_generator_mut() = L::initialize_from_triangulation(&result);
        Ok(result)
    }

    /// Stable bulk load variant that preserves the input vertex order
    ///
    /// The resulting vertex set will be equal to the input vertex set if their positions are all distinct.
    /// See [ConstrainedDelaunayTriangulation::bulk_load_cdt] for additional details like panic behavior and duplicate
    /// handling.
    ///
    /// # Example
    /// ```
    /// # fn main() -> Result<(), spade::InsertionError> {
    /// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    /// let mut vertices = vec![
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(3.0, -3.0),
    ///     Point2::new(-1.0, -2.0),
    ///     Point2::new(-4.0, -5.0),
    /// ];
    /// let mut edges = vec![[0, 1], [1, 2], [2, 3], [3, 4]];
    /// let cdt = ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt_stable(vertices.clone(), edges)?;
    ///
    /// // The ordered will be preserved:
    /// assert_eq!(cdt.vertices().map(|v| v.position()).collect::<Vec<_>>(), vertices);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// It is fine to include vertex positions multiple times. The resulting order will be the same as if
    /// the duplicates were removed prior to insertion. However, it is unclear *which* duplicates are
    /// removed - e.g. do not assume that always the first duplicated vertex remains.
    ///
    /// ```
    /// # fn main() -> Result<(), spade::InsertionError> {
    /// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    /// let mut vertices = vec![
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 2.0), // Duplicate
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(3.0, -3.0),
    ///     Point2::new(3.0, -3.0), // Duplicate
    ///     Point2::new(-4.0, -5.0),
    /// ];
    /// let mut edges = vec![[0, 1], [2, 3], [4, 5]];
    /// let cdt = ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt_stable(vertices.clone(), edges)?;
    ///
    /// // The choice of deduplicated vertices is arbitrary. In this example, dedup[1] and dedup[2] could
    /// // have been swapped
    /// let dedup = [
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(3.0, -3.0),
    ///     Point2::new(-4.0, -5.0),
    /// ];
    /// assert_eq!(cdt.vertices().map(|v| v.position()).collect::<Vec<_>>(), dedup);
    /// # Ok(())
    /// # }
    /// ```
    pub fn bulk_load_cdt_stable(
        vertices: Vec<V>,
        edges: Vec<[usize; 2]>,
    ) -> Result<Self, InsertionError> {
        let mut result: Self =
            bulk_load_stable(move |vertices| bulk_load_cdt(vertices, edges), vertices)?;
        *result.hint_generator_mut() = L::initialize_from_triangulation(&result);
        Ok(result)
    }

    /// Removes a vertex from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles.
    pub fn remove(&mut self, vertex: FixedVertexHandle) -> V {
        let num_removed_constraints = self
            .dcel
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
        self.dcel.undirected_edge_data(edge).is_constraint_edge()
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
    ///
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
    ///
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
    /// Note that the given constraint might be split into smaller edges
    /// if a vertex in the triangulation lies exactly on the constraint edge.
    /// Thus, `cdt.exists_constraint(from, to)` is not necessarily `true`
    /// after a call to this function.
    ///
    /// Returns false and does nothing if `from == to`.
    ///
    /// # Panics
    ///
    /// Panics if the new constraint edge intersects an existing
    /// constraint edge. Use [Self::try_add_constraint] or [Self::add_constraint_and_split] to work
    /// around that.
    pub fn add_constraint(&mut self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        let initial_num_constraints = self.num_constraints();
        self.try_add_constraint_inner(from, to, |_| panic!("Constraint edges must not intersect."));

        self.num_constraints != initial_num_constraints
    }

    /// Takes a conflict region (expressed as a list of intersecting edges) rotates edges to create
    /// a new constraint edge. Then, the rotated edges (except the new constraint edge)
    /// are legalized to restore the Delaunay property.
    ///
    /// Usually, this step is described as "delete all conflicting edges, then re-triangulate the
    /// hole". Spade avoids the removal of edges by _rotating_ (flipping) them into place instead.
    /// The final constraint edge is created implicitly.
    /// This works as long as the intersecting edges are ordered "along the constraint edge", i.e.
    /// the intersection points increase in distance from the constraint edge origin.
    ///
    /// # Example
    ///
    /// The input conflict region might look like this (assuming the target constraint edge goes
    /// from v0 to v1):
    ///
    /// ```text
    ///     v__________v
    ///   / |        / |\
    ///  /  |      /   | \
    /// v0  |e0  /e1 e2| v1
    ///  \  |  /       | /
    ///   \ |/         |/
    ///     v_________ v
    /// ```
    ///
    /// `conflict_edges` would be set to `vec![e0, e1, e2]` in this case, `target_vertex` would be
    /// `v1`.
    ///
    /// Now, flipping these edges _in this order_ will implicitly create the desired edge:
    ///
    /// After flipping the result looks like this with all edges going out of `v0`:
    ///
    /// ```text
    ///     v_________v
    ///   /     __---  \
    ///  / __---        \
    /// v0--------------v1  
    ///  \ --___        /
    ///   \     --___  /
    ///     v---------v
    ///```
    ///
    /// Now, the new edges can be legalized as usual.
    ///
    /// Returns a handle to the new constraint edge (pointing toward `target_vertex`).
    fn resolve_conflict_region(
        &mut self,
        conflict_edges: Vec<FixedDirectedEdgeHandle>,
        target_vertex: FixedVertexHandle,
    ) -> Option<FixedDirectedEdgeHandle> {
        let first = conflict_edges.first()?;

        let mut temporary_constraint_edges = Vec::new();

        let first = self.directed_edge(*first);

        // These refer to the two edges that go out of the constraint edge origin initially.
        // They are used below but need to be defined declared here to appease the borrow checker.
        let first_border_edge = first.rev().prev().fix();
        let last_border_edge = first.rev().next().fix();

        // Flip all conflict edges in the input order - see function comment.
        for edge in &conflict_edges {
            flip_cw(self.s_mut(), edge.as_undirected());
        }

        // Small optimization: For the legalization, the algorithm doesn't need to look at edges
        // outside the conflict region. They are known to be already legal.
        // To do so, we will make the border edges that encompass the conflict region into temporary
        // constraint edges. The legalization will then skip them. This is undone later,
        let mut make_temporary_edge = |cdt: &mut Self, edge: FixedUndirectedEdgeHandle| {
            // Exclude edges that are already a constraint - those should remain constraint edges
            // and not be undone later!
            if !cdt.undirected_edge(edge).is_constraint_edge() {
                temporary_constraint_edges.push(edge);
                cdt.undirected_edge_data_mut(edge).make_constraint_edge();
            }
        };

        make_temporary_edge(self, first_border_edge.as_undirected());
        make_temporary_edge(self, last_border_edge.as_undirected());

        let mut current = first_border_edge;

        let mut result = None;

        // Loops around all border edges and adds them to the temporary constraint edge list.
        // `first_border_edge` and `last_border_edge` refer to the two border edges that are
        // initially going out of the constraint edge start (the two left most edges in the first
        // ascii drawing of the function comment).
        while current != last_border_edge.rev() {
            let handle = self.directed_edge(current);
            let fixed = handle.fix();
            let next = handle.next().fix().as_undirected();

            current = handle.ccw().fix();
            if target_vertex == handle.to().fix() {
                // This loop also finds the implicitly created constraint edge and makes it an
                // official constraint edge!
                self.make_constraint_edge(fixed.as_undirected());
                result = Some(fixed);
            }
            make_temporary_edge(self, next);
        }

        self.legalize_edges_after_removal(
            &mut conflict_edges
                .into_iter()
                .map(|edge| edge.as_undirected())
                .collect(),
            |_| false,
        );

        // Undo the previously made temporary constraint edges
        for edge in temporary_constraint_edges {
            self.undirected_edge_data_mut(edge).0 = false;
        }

        result
    }

    /// Returns all constraint edges that would prevent creating a new constraint between two points.
    ///
    /// # See also
    ///
    /// See also [Self::get_conflicting_edges_between_vertices]
    pub fn get_conflicting_edges_between_points(
        &self,
        from: Point2<<V as HasPosition>::Scalar>,
        to: Point2<<V as HasPosition>::Scalar>,
    ) -> impl Iterator<Item = DirectedEdgeHandle<V, DE, CdtEdge<UE>, F>> {
        LineIntersectionIterator::new(self, from, to)
            .flat_map(|intersection| intersection.as_edge_intersection())
            .filter(|e| e.is_constraint_edge())
    }

    /// Returns all constraint edges that would prevent inserting a new constraint connecting two existing
    /// vertices.
    ///
    /// # See also
    ///
    /// See also [Self::get_conflicting_edges_between_points]
    pub fn get_conflicting_edges_between_vertices(
        &self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
    ) -> impl Iterator<Item = DirectedEdgeHandle<V, DE, CdtEdge<UE>, F>> {
        LineIntersectionIterator::new_from_handles(self, from, to)
            .flat_map(|intersection| intersection.as_edge_intersection())
            .filter(|e| e.is_constraint_edge())
    }

    fn make_constraint_edge(&mut self, edge: FixedUndirectedEdgeHandle) -> bool {
        if !self.is_constraint_edge(edge) {
            self.dcel
                .undirected_edge_data_mut(edge)
                .make_constraint_edge();
            self.num_constraints += 1;
            true
        } else {
            false
        }
    }

    #[cfg(any(test, fuzzing))]
    pub fn cdt_sanity_check(&self) {
        self.cdt_sanity_check_with_params(true);
    }

    #[cfg(any(test, fuzzing))]
    pub fn cdt_sanity_check_with_params(&self, check_convexity: bool) {
        let num_undirected_edges = self
            .dcel
            .undirected_edges()
            .filter(|e| e.is_constraint_edge())
            .count();

        assert_eq!(num_undirected_edges, self.num_constraints());

        if self.num_constraints() == 0 && check_convexity {
            self.sanity_check();
        } else {
            self.basic_sanity_check(check_convexity);
        }
    }

    /// Removes a constraint edge.
    ///
    /// Does nothing and returns `false` if the given edge is not a constraint edge.
    /// Otherwise, the edge is unmarked and the Delaunay property is restored in its vicinity.
    pub fn remove_constraint_edge(&mut self, edge: FixedUndirectedEdgeHandle) -> bool {
        if self.is_constraint_edge(edge) {
            self.dcel
                .undirected_edge_data_mut(edge)
                .unmake_constraint_edge();
            self.num_constraints -= 1;
            self.legalize_edge(edge.as_directed(), true);
            true
        } else {
            false
        }
    }

    /// Attempts to add a constraint edge. Leaves the triangulation unchanged if the new edge would
    /// intersect an existing constraint edge.
    ///
    /// Returns all constraint edges that connect `from` and `to`. This includes any constraint
    /// edge that was already present.
    /// Multiple edges are returned if the line from `from` to `to` intersects an existing vertex.
    /// Returns an empty list if the new constraint would intersect any existing constraint or if
    /// `from == to`.
    ///
    /// # Example
    ///
    /// ```
    /// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    /// # fn try_main() -> Result<(), spade::InsertionError> {
    /// let mut cdt = ConstrainedDelaunayTriangulation::<Point2<_>>::new();
    /// let v0 = cdt.insert(Point2::new(-1.0, 0.0))?;
    /// let v1 = cdt.insert(Point2::new(1.0, 0.0))?;
    /// let v2 = cdt.insert(Point2::new(0.0, 1.0))?;
    /// let v3 = cdt.insert(Point2::new(0.0, -1.0))?;
    /// let first_constraints = cdt.try_add_constraint(v2, v3);
    /// let second_constraints = cdt.try_add_constraint(v0, v1);
    ///
    /// // The first constraint edge can be added as there are no intersecting constraint edges
    /// assert_eq!(first_constraints.len(), 1);
    /// let edge = cdt.directed_edge(first_constraints[0]);
    /// assert_eq!(edge.from().fix(), v2);
    /// assert_eq!(edge.to().fix(), v3);
    ///
    /// // The second edge should not be created as it intersects the first edge.
    /// assert!(second_constraints.is_empty());
    ///
    /// // Consider comparing this to the number of constraints prior to calling
    /// // `try_add_constraint` to check if any new constraint edge was created.
    /// assert_eq!(cdt.num_constraints(), 1);
    /// # Ok(()) }
    /// # fn main() { try_main().unwrap() }
    /// ```
    pub fn try_add_constraint(
        &mut self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
    ) -> Vec<FixedDirectedEdgeHandle> {
        self.try_add_constraint_inner(from, to, |_| ConflictResolution::Cancel)
    }

    fn try_add_constraint_inner<R>(
        &mut self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
        mut conflict_resolver: R,
    ) -> Vec<FixedDirectedEdgeHandle>
    where
        R: FnMut(DirectedEdgeHandle<V, DE, CdtEdge<UE>, F>) -> ConflictResolution<V>,
    {
        // Constraint edges are added with a two-pass approach:
        // - First, identify potential constraint edge intersections (conflicts). This must be done
        //   beforehand in case that the caller chooses to `ConflictResolution::Cancel` the
        //   operation - no mutation should have happened at this stage.
        let (initial_conflict_regions, all_regions_intact) =
            self.get_conflict_resolutions(from, to, &mut conflict_resolver);
        // - Second, apply the conflict resolutions, e.g. by inserting new split points and by
        //   rotating non-constraint edges that intersect the new constraint edge (see function
        //   `resolve_conflict_region`).

        if all_regions_intact {
            self.resolve_conflict_groups(to, initial_conflict_regions)
        } else {
            self.add_splitting_constraint_edge_fallback(initial_conflict_regions)
        }
    }

    /// Fallback routine to add splitting constraints that cannot be added via the fast path.
    ///
    /// This routine simply adds all split vertices first and then adds any missing constraints.
    /// This avoids edge cases that can arise when the split point for a constraint
    /// intersection lies "too far" off the conflict edge due to imprecise calculations.
    fn add_splitting_constraint_edge_fallback(
        &mut self,
        initial_conflict_regions: Vec<InitialConflictRegion<V>>,
    ) -> Vec<FixedDirectedEdgeHandle> {
        let mut vertices_to_connect = Vec::new();
        let mut temporarily_removed = Vec::new();

        // Phase 1: Add all pending split vertices directly.
        for region in initial_conflict_regions {
            let group_end_vertex = match region.group_end {
                Existing(v) => v,
                ConflictRegionEnd::ConstraintEdgeSplit(new_vertex, edge) => {
                    let new_handle = match new_vertex {
                        Ok(new_vertex) => self
                            .insert(new_vertex)
                            .expect("Failed to insert vertex as expected. This is a bug in spade."),
                        Err(handle) => handle,
                    };

                    let [old_from, old_to] = self.directed_edge(edge).vertices().map(|v| v.fix());
                    // The conflict edge can prevent the forced insertion to the split vertex.
                    // It will be removed temporarily
                    self.remove_constraint_edge(edge.as_undirected());

                    // Re-add the temporarily removed edge later as if it was split by the new
                    // vertex
                    temporarily_removed.push([old_from, new_handle]);
                    temporarily_removed.push([new_handle, old_to]);
                    new_handle
                }
                EdgeOverlap(edge) => self.directed_edge(edge).to().fix(),
            };

            vertices_to_connect.push(group_end_vertex);
        }

        let mut result = Vec::new();
        let mut last_vertex = None;

        // Phase 2: Add all constraint edges
        for vertex in vertices_to_connect {
            if let Some(last_vertex) = last_vertex {
                let new_edges = self.try_add_constraint(last_vertex, vertex);
                // try_add_constraint should always succeed as any conflicting edge should have been
                // removed temporarily
                assert_ne!(new_edges, Vec::new());
                result.extend(new_edges);
            }

            last_vertex = Some(vertex);
        }

        for [from, to] in temporarily_removed {
            self.try_add_constraint(from, to);
        }

        result
    }

    fn get_conflict_resolutions<R>(
        &mut self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
        conflict_resolver: &mut R,
    ) -> (Vec<InitialConflictRegion<V>>, bool)
    where
        R: FnMut(DirectedEdgeHandle<V, DE, CdtEdge<UE>, F>) -> ConflictResolution<V>,
    {
        let mut all_regions_intact = true;
        let mut conflict_groups = Vec::new();
        let mut current_group = Vec::new();
        let mut ignore_next_vertex = false;
        for intersection in LineIntersectionIterator::new_from_handles(self, from, to) {
            match intersection {
                Intersection::EdgeIntersection(edge) => {
                    if !edge.is_constraint_edge() {
                        current_group.push(edge.fix());
                        continue;
                    }

                    // The new constraint intersects an existing constraint edge. Start conflict
                    // resolution.
                    match conflict_resolver(edge) {
                        ConflictResolution::Cancel => {
                            return (Vec::new(), true);
                        }
                        ConflictResolution::Split(new_vertex) => {
                            let position = new_vertex.position();
                            let (overlap_vertex, is_valid) =
                                self.verify_split_position(edge, position);

                            // A region is considered to be intact if the split vertex lies
                            // within the region and not outside or on its border.
                            all_regions_intact &= is_valid;

                            let conflict_edges = core::mem::take(&mut current_group);

                            // overlap_vertex.is_some() indicates that the split position
                            // overlaps an existing vertex. This can happen due to rounding
                            // errors and needs some special handling
                            ignore_next_vertex = overlap_vertex.is_some();

                            let group_end_vertex =
                                overlap_vertex.map(|h| Err(h)).unwrap_or(Ok(new_vertex));

                            let group_end = ConflictRegionEnd::ConstraintEdgeSplit(
                                group_end_vertex,
                                edge.fix(),
                            );

                            conflict_groups.push(InitialConflictRegion {
                                conflict_edges,
                                group_end,
                            });
                        }
                    }
                }
                Intersection::VertexIntersection(v) => {
                    if ignore_next_vertex {
                        ignore_next_vertex = false;
                        continue;
                    }
                    let group_end = Existing(v.fix());
                    let conflict_edges = core::mem::take(&mut current_group);
                    conflict_groups.push(InitialConflictRegion {
                        conflict_edges,
                        group_end,
                    });
                }
                Intersection::EdgeOverlap(edge) => {
                    conflict_groups.push(InitialConflictRegion {
                        conflict_edges: Vec::new(),
                        group_end: EdgeOverlap(edge.fix()),
                    });
                    // The next intersection is going to be edge.to(). It would be incorrect to
                    // create a conflict region from that vertex as that region is already handled
                    // by the GroupEnd::EdgeOverlap cases
                    ignore_next_vertex = true;
                }
            }
        }

        (conflict_groups, all_regions_intact)
    }

    fn verify_split_position(
        &self,
        conflict_edge: DirectedEdgeHandle<V, DE, CdtEdge<UE>, F>,
        split_position: Point2<<V as HasPosition>::Scalar>,
    ) -> (Option<FixedVertexHandle>, bool) {
        // Not every split vertex may lead to a conflict region that will properly contain the
        // split vertex. This can happen as not all split positions can be represented precisely.
        //
        // Instead, these vertices will be handled by a slower fallback routine.
        //
        // A split position is considered to be valid if it lies either *on* the edge that was split
        // or *within any of the edges neighboring faces*.
        match self.locate_with_hint(split_position, conflict_edge.from().fix()) {
            PositionInTriangulation::OnEdge(real_edge) => {
                let is_valid = real_edge.as_undirected() == conflict_edge.fix().as_undirected();
                (None, is_valid)
            }
            PositionInTriangulation::OnFace(face) => {
                let face = face.adjust_inner_outer();
                let is_valid =
                    face == conflict_edge.face().fix() || face == conflict_edge.rev().face().fix();
                (None, is_valid)
            }
            PositionInTriangulation::OutsideOfConvexHull(_) => {
                let is_valid = conflict_edge.is_part_of_convex_hull();
                (None, is_valid)
            }
            PositionInTriangulation::OnVertex(v) => (Some(v), false),
            PositionInTriangulation::NoTriangulation => unreachable!(),
        }
    }

    fn resolve_conflict_groups(
        &mut self,
        final_vertex: FixedVertexHandle,
        conflict_groups: Vec<InitialConflictRegion<V>>,
    ) -> Vec<FixedDirectedEdgeHandle> {
        let mut constraint_edges = Vec::new();
        let mut last_vertex = None;

        for InitialConflictRegion {
            conflict_edges,
            group_end,
        } in conflict_groups
        {
            let mut last_edge = None;
            let target_vertex = match group_end {
                Existing(v) => v,
                ConflictRegionEnd::ConstraintEdgeSplit(v, conflict_edge) => {
                    let v = v.expect(
                        "Expected a new vertex for insertion. \
                    An existing vertex should be handled by the fallback routine. \
                    This is a bug in spade.",
                    );
                    let (new_vertex, [e0, e1]) = self.insert_on_edge(conflict_edge, v);
                    let e1_handle = self.directed_edge(e1);
                    // edge_in / edge_out refer to the edge going into / out of the newly split off
                    // vertex.
                    let edge_out = e1_handle.ccw();
                    let edge_in = e1_handle.cw();

                    if Some(edge_in.to().fix()) == last_vertex {
                        constraint_edges.push(edge_in.fix().rev());
                    }

                    if edge_out.to().fix() == final_vertex {
                        // The edge reaches the target vertex - we're done! However, this can
                        // sometimes omit to make the last edge a constraint. This special case
                        // fixes that issue.
                        last_edge = Some(edge_out.fix());
                    }

                    self.handle_legal_edge_split([e0, e1]);
                    new_vertex
                }
                EdgeOverlap(edge) => {
                    constraint_edges.push(edge);
                    last_vertex = Some(self.directed_edge(edge).to().fix());
                    // No need to resolve conflict regions - there are no conflicting edges in the
                    // GroupEnd::EdgeOverlap case
                    continue;
                }
            };

            constraint_edges.extend(self.resolve_conflict_region(conflict_edges, target_vertex));
            constraint_edges.extend(last_edge);

            last_vertex = Some(target_vertex);
        }

        for edge in &constraint_edges {
            self.make_constraint_edge(edge.as_undirected());
        }

        constraint_edges
    }
}

impl<V, DE, UE, F, L> ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    V::Scalar: Float,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    /// Adds a constraint to the triangulation. Splits any existing constraint edge that would
    /// intersect the new constraint edge.
    ///
    /// The `vertex_constructor` closure is used to convert the position of the intersection into
    /// a vertex. The returned vertex must have exactly the same position as the argument of the
    /// closure.
    ///
    /// Returns all constraint edges that connect `from` and `to`. This includes any constraint
    /// edge that was already present.
    /// Multiple edges are returned if the line from `from` to `to` intersects any existing vertex
    /// or any existing constraint edge.
    /// Returns an empty list if `from == to`.
    ///
    /// # Image example
    ///
    /// This is an input CDT with 3 constraints:
    ///
    #[doc = include_str!("../images/add_constraint_and_split_initial.svg")]
    ///
    /// Calling `add_constraint_and_split(v0, v1, ...)` will result in this CDT:
    ///
    #[doc = include_str!("../images/add_constraint_and_split_added.svg")]
    ///
    /// # Code example
    ///
    ///```
    /// use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    /// # fn try_main() -> Result<(), spade::InsertionError> {
    /// use spade::handles::FixedVertexHandle;
    /// let mut cdt = ConstrainedDelaunayTriangulation::<Point2<_>>::new();
    /// let v0 = cdt.insert(Point2::new(-1.0, 0.0))?;
    /// let v1 = cdt.insert(Point2::new(1.0, 0.0))?;
    /// let v2 = cdt.insert(Point2::new(0.0, 1.0))?;
    /// let v3 = cdt.insert(Point2::new(0.0, -1.0))?;
    /// cdt.add_constraint(v2, v3);
    ///
    /// // Should create a new split vertex at the origin
    /// let second_constraints = cdt.add_constraint_and_split(v0, v1, |v| v);
    ///
    /// // Expect one additional point introduced by splitting the first constraint edge:
    /// assert_eq!(cdt.num_vertices(), 5);
    ///
    /// let v4 = FixedVertexHandle::from_index(4); // Newly created
    ///
    /// // Expect 4 constraints as the first constraint was split:
    /// assert_eq!(cdt.num_constraints(), 4);
    ///
    /// // The second edge should consist of two edges, v0 -> v4 and v4 -> v1
    /// assert_eq!(second_constraints.len(), 2);
    ///
    /// let [e0, e1] = [second_constraints[0], second_constraints[1]];
    /// let [e0, e1] = [e0, e1].map(|e| cdt.directed_edge(e));
    ///
    /// assert_eq!(e0.from().fix(), v0);
    /// assert_eq!(e0.to().fix(), v4);
    /// assert_eq!(e1.from().fix(), v4);
    /// assert_eq!(e1.to().fix(), v1);
    ///
    /// # Ok(()) }
    /// # fn main() { try_main().unwrap() }
    /// ```
    ///
    /// # Precision warning
    ///
    /// Intersection points may not _exactly_ lie on the line between `from` and `to`, either due to
    /// loss of precision or as the exact value may not be representable with the underlying
    /// floating point number.
    ///
    /// Thus, iterating a `LineIntersectionIterator::new_from_handles(&cdt, from, to)` will often
    /// not return only `Intersection::EdgeOverlap` as would be expected. Instead, use the returned
    /// `Vec` to identify the edges that form the new constraint.
    /// The absolute deviation from the correct position should be minimal, especially when using
    /// `f64` coordinates as storage type.
    pub fn add_constraint_and_split<C>(
        &mut self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
        vertex_constructor: C,
    ) -> Vec<FixedDirectedEdgeHandle>
    where
        C: Fn(Point2<<V as HasPosition>::Scalar>) -> V,
    {
        let from_pos = self.vertex(from).position();
        let to_pos = self.vertex(to).position();

        self.try_add_constraint_inner(from, to, |edge| {
            let [p0, p1] = edge.positions();
            let line_intersection = get_edge_intersections(p0, p1, from_pos, to_pos);
            let new_vertex = vertex_constructor(line_intersection);
            assert_eq!(new_vertex.position(), line_intersection);
            ConflictResolution::Split(new_vertex)
        })
    }
}

/// Describes all possible ways in which conflict regions which are created while adding a
/// constraint edge may end.
enum ConflictRegionEnd<V> {
    /// Conflict group ends with an existing vertex
    Existing(FixedVertexHandle),
    /// Special case of "Existing" - the constraint edge overlaps any existing edge which implies
    /// that the conflict group also ends on an existing vertex.
    /// However, it makes sense to handle this specially to prevent having to look up the overlapped
    /// edge later.
    EdgeOverlap(FixedDirectedEdgeHandle),
    /// The conflict region ends in a vertex that splits an existing constraint edge. Usually, this
    /// vertex is constructed anew and given by the `Ok` case.
    /// In rare cases, the split vertex may be an existing vertex that does not lie exactly on the
    /// line due to rounding issues. This is indicated by the `Err` case. The constraint edge that
    /// should be split is the second field.
    ConstraintEdgeSplit(Result<V, FixedVertexHandle>, FixedDirectedEdgeHandle),
}

/// Represents a conflict region that does not yet fully exist as a vertex may be missing. This can
/// happen if adding a constraint edge should split any intersecting existing edge.
/// This will eventually be turned into a "real" conflict group (described as a list of edges) by
/// inserting the missing vertex.
struct InitialConflictRegion<V> {
    conflict_edges: Vec<FixedDirectedEdgeHandle>,
    group_end: ConflictRegionEnd<V>,
}

enum ConflictResolution<V> {
    Cancel,
    Split(V),
}

pub fn get_edge_intersections<S: SpadeNum + Float>(
    p1: Point2<S>,
    p2: Point2<S>,
    p3: Point2<S>,
    p4: Point2<S>,
) -> Point2<S> {
    let p1 = p1.to_f64();
    let p2 = p2.to_f64();
    let p3 = p3.to_f64();
    let p4 = p4.to_f64();

    let a1 = p2.y - p1.y;
    let b1 = p1.x - p2.x;
    let c1 = a1 * p1.x + b1 * p1.y;

    let a2 = p4.y - p3.y;
    let b2 = p3.x - p4.x;
    let c2 = a2 * p3.x + b2 * p3.y;

    let determinant = a1 * b2 - a2 * b1;

    let x: f64;
    let y: f64;
    if determinant == 0.0 {
        x = f64::infinity();
        y = f64::infinity();
    } else {
        x = (b2 * c1 - b1 * c2) / determinant;
        y = (a1 * c2 - a2 * c1) / determinant;
    }

    [x, y]
        .map(|s| <S as NumCast>::from(s).unwrap_or_else(|| (s as f32).into()))
        .into()
}

#[cfg(test)]
mod test {
    use alloc::{vec, vec::Vec};

    use rand::distributions::{Distribution, Uniform};
    use rand::{Rng, SeedableRng};

    use crate::delaunay_core::{FixedDirectedEdgeHandle, TriangulationExt};
    use crate::handles::FixedVertexHandle;
    use crate::test_utilities::*;
    use crate::{DelaunayTriangulation, InsertionError, Point2, Triangulation};

    use super::ConstrainedDelaunayTriangulation;

    type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;
    type Delaunay = DelaunayTriangulation<Point2<f64>>;

    #[test]
    fn test_into() -> Result<(), InsertionError> {
        let points = random_points_with_seed(100, SEED);
        let delaunay = DelaunayTriangulation::<_>::bulk_load(points.clone())?;

        let cdt = Cdt::from(delaunay.clone());

        assert_eq!(delaunay.num_vertices(), cdt.num_vertices());
        assert_eq!(delaunay.num_directed_edges(), cdt.num_directed_edges());
        assert_eq!(cdt.num_constraints, 0);

        Ok(())
    }

    #[test]
    fn test_add_same_from_and_to_constraint() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0))?;
        cdt.insert(Point2::new(2.0, 2.0))?;
        cdt.insert(Point2::new(1.0, 2.0))?;

        assert!(!cdt.add_constraint(v0, v0));
        assert!(cdt.try_add_constraint(v0, v0).is_empty());

        let new_point = Point2::new(3.1, 2.0);
        assert!(!cdt.add_constraint_edge(new_point, new_point)?);

        assert_eq!(0, cdt.num_constraints());
        assert_eq!(4, cdt.num_vertices());

        cdt.cdt_sanity_check();

        Ok(())
    }

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
        let mut max_y = -f64::MAX;
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
        // Use a delaunay triangulation to "generate" non-intersecting constraint edges
        let mut d = Delaunay::new();
        for p in delaunay_points {
            d.insert(p)?;
        }
        let mut used_vertices = hashbrown::HashSet::new();

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
        // them with non-intersecting constraints
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
                let angle = core::f64::consts::PI * 2.0 * i as f64 / NUM_VERTICES as f64;
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

        cdt.add_constraint_edges(core::iter::empty(), false)?;
        cdt.add_constraint_edges(core::iter::empty(), true)?;

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

    #[test]
    fn test_cdt_edge_split_degenerate() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        cdt.add_constraint_edge(Point2::new(-10.0, -10.0), Point2::new(20.0, -10.0))?;
        cdt.insert(Point2::new(0.0, -10.0))?;

        assert_eq!(cdt.num_constraints(), 2);

        Ok(())
    }

    #[test]
    fn infinite_loop_bug() -> Result<(), InsertionError> {
        // See https://github.com/Stoeoef/spade/issues/98
        let mut triangulation = Cdt::default();

        let start = Point2::new(-21.296192169189453, 9.872323036193848);
        let edges = [
            (
                Point2::new(-20.926544189453125, 16.53529167175293),
                Point2::new(-27.772645950317383, 4.197676658630371),
            ),
            (
                Point2::new(-20.03745460510254, 12.93730354309082),
                Point2::new(-20.930097579956055, 11.93786907196045),
            ),
            (
                Point2::new(-15.576859474182129, 8.772907257080078),
                Point2::new(-22.373262405395508, 12.348699569702148),
            ),
            (
                Point2::new(-10.038422584533691, 5.663522243499756),
                Point2::new(-16.382625579833984, 9.09498119354248),
            ),
            (
                Point2::new(0.0, 0.0),
                Point2::new(-13.11422061920166, 7.30709171295166),
            ),
            (
                Point2::new(-19.230497360229492, -3.7645812034606934),
                Point2::new(-7.411926746368408, 3.486957311630249),
            ),
            (
                Point2::new(-25.072885513305664, -9.239323616027832),
                Point2::new(-19.462360382080078, -1.621320366859436),
            ),
            (
                Point2::new(-32.41080856323242, -13.72575855255127),
                Point2::new(-22.58626365661621, -2.076631784439087),
            ),
            (
                Point2::new(-32.41080856323242, -13.72575855255127),
                Point2::new(-25.57504653930664, -4.952820301055908),
            ),
            (
                Point2::new(-33.08932113647461, 0.31093916296958923),
                Point2::new(-25.955543518066406, 0.18878456950187683),
            ),
        ];

        for (p1, p2) in edges {
            let p1 = triangulation.insert(p1)?;
            let p2 = triangulation.insert(p2)?;
            assert!(triangulation.can_add_constraint(p1, p2));
            triangulation.add_constraint(p1, p2);
        }

        triangulation.insert(start)?;
        Ok(())
    }

    #[test]
    pub fn infinite_loop_2() -> Result<(), InsertionError> {
        let lines = [
            [
                Point2::new(0.9296344883099084, 0.03071359966930065),
                Point2::new(0.26031306872107085, 0.34491289915959455),
            ],
            [
                Point2::new(0.7384289920396423, 0.4981747664368982),
                Point2::new(0.06543525273452533, 0.34139896206401854),
            ],
            [
                Point2::new(0.9535295221136963, 0.9114305148801416),
                Point2::new(0.8306091165247367, 0.08959389670590667),
            ],
        ];

        let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();

        for [a, b] in lines {
            let a = cdt.insert(a)?;
            let b = cdt.insert(b)?;

            cdt.add_constraint_and_split(a, b, |v| v);
        }

        // This insertion used to fail as the position could not be located
        cdt.insert(Point2::new(0.5138795569454557, 0.3186272217036502))?;
        Ok(())
    }

    fn get_cdt_for_try_add_constraint() -> Result<Cdt, InsertionError> {
        let vertices = vec![
            Point2::new(0.0, -10.0),
            Point2::new(76.0, 0.0),
            Point2::new(20.0, 20.0),
            Point2::new(20.0, -30.0),
            Point2::new(45.0, 25.0),
            Point2::new(32.0, -35.0),
            Point2::new(60.0, 20.0),
            Point2::new(60.0, -30.0),
            Point2::new(50.0, -34.0),
        ];

        Cdt::bulk_load_cdt_stable(vertices, vec![[3, 2], [5, 4], [7, 6]])
    }

    #[test]
    fn test_single_split() -> Result<(), InsertionError> {
        let vertices = vec![
            Point2::new(-1.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, -1.0),
            Point2::new(0.0, 1.0),
        ];

        let mut cdt = Cdt::bulk_load_cdt_stable(vertices, vec![[2, 3]])?;

        let initial_num_vertices = cdt.num_vertices();
        let from = FixedVertexHandle::from_index(0);
        let to = FixedVertexHandle::from_index(1);

        let edges = cdt.add_constraint_and_split(from, to, |v| v);

        assert_eq!(cdt.num_vertices(), initial_num_vertices + 1);
        assert_eq!(edges.len(), 2);
        check_returned_edges(&mut cdt, &edges, from, to);

        Ok(())
    }

    #[test]
    fn test_multiple_splits() -> Result<(), InsertionError> {
        let mut cdt = get_cdt_for_try_add_constraint()?;

        let initial_num_vertices = cdt.num_vertices();
        let from = FixedVertexHandle::from_index(0);
        let to = FixedVertexHandle::from_index(1);

        let edges = cdt.add_constraint_and_split(from, to, |v| v);

        // 3 new points should be added as the constraint intersects all 3 existing edges
        assert_eq!(cdt.num_vertices(), initial_num_vertices + 3);
        assert_eq!(edges.len(), 4);
        check_returned_edges(&mut cdt, &edges, from, to);

        Ok(())
    }

    #[test]
    fn test_try_add_constraint() -> Result<(), InsertionError> {
        let mut cdt = get_cdt_for_try_add_constraint()?;

        let initial_num_vertices = cdt.num_vertices();
        let initial_num_constraints = cdt.num_constraints();
        let from = FixedVertexHandle::from_index(0);
        let to = FixedVertexHandle::from_index(1);

        // Is expected to fail (return an empty list)
        let edges = cdt.try_add_constraint(from, to);
        assert_eq!(edges, Vec::new());
        assert_eq!(cdt.num_vertices(), initial_num_vertices);
        assert_eq!(cdt.num_constraints(), initial_num_constraints);

        let from = FixedVertexHandle::from_index(2);
        let to = FixedVertexHandle::from_index(3);

        // Try to add on top of an existing edge
        let edges = cdt.try_add_constraint(from, to);
        assert_eq!(edges.len(), 1);

        Ok(())
    }

    #[test]
    fn test_remove_constraint_edge() -> Result<(), InsertionError> {
        let mut cdt = get_cdt_for_try_add_constraint()?;
        for edge in cdt.fixed_undirected_edges() {
            cdt.remove_constraint_edge(edge);
        }
        assert_eq!(cdt.num_constraints, 0);
        cdt.sanity_check();

        let added_edges = cdt.try_add_constraint(
            FixedVertexHandle::from_index(0),
            FixedVertexHandle::from_index(1),
        );
        assert_eq!(added_edges.len(), 1);

        assert!(cdt.remove_constraint_edge(added_edges.first().unwrap().as_undirected()));
        assert_eq!(cdt.num_constraints, 0);
        cdt.sanity_check();

        Ok(())
    }

    #[test]
    fn edge_intersection_precision_test() -> Result<(), InsertionError> {
        let edges = [
            [
                Point2::new(17.064112, -17.96008),
                Point2::new(16.249594, -17.145563),
            ],
            [
                Point2::new(-25.290726, -24.435482),
                Point2::new(-5.6608872, -24.435482),
            ],
            [
                Point2::new(17.878626, -18.774595),
                Point2::new(15.435078, -16.331045),
            ],
        ];
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f32>> =
            ConstrainedDelaunayTriangulation::new();

        for edge in edges.iter() {
            let point_a = cdt.insert(edge[0])?;
            let point_b = cdt.insert(edge[1])?;

            // The intersection calculation of the last edge is susceptible to floating point
            // inaccuracies. Spade has a fallback routine that is more costly but should handle
            // these more robustly. This test is set up to trigger this routine.
            cdt.add_constraint_and_split(point_a, point_b, |v| v);
            cdt.cdt_sanity_check();
        }

        assert_eq!(cdt.num_vertices(), 7);

        // Gather all constraint edges as [from, to] index tuples
        let mut constraint_edges = cdt
            .undirected_edges()
            .filter(|e| e.is_constraint_edge())
            .map(|e| e.vertices().map(|v| v.index()))
            .collect::<Vec<_>>();

        // Normalize to make comparison order-independent
        for edge_pair in &mut constraint_edges {
            edge_pair.sort();
        }
        constraint_edges.sort();

        // Manually checked for correctness...
        assert_eq!(
            constraint_edges,
            vec![[0, 6], [1, 6], [2, 3], [4, 6], [5, 6]]
        );

        Ok(())
    }

    #[test]
    fn edge_intersection_precision_test_2() -> Result<(), InsertionError> {
        let edges = [
            [
                Point2 {
                    x: 18.69314193725586,
                    y: 19.589109420776367,
                },
                Point2 {
                    x: 18.69314193725586,
                    y: 20.40362548828125,
                },
            ],
            [
                Point2 {
                    x: 19.507659912109375,
                    y: 20.40362548828125,
                },
                Point2 {
                    x: 17.878625869750977,
                    y: 18.774595260620117,
                },
            ],
            [
                Point2 {
                    x: 20.322175979614258,
                    y: 21.218143463134766,
                },
                Point2 {
                    x: 15.435077667236328,
                    y: 16.331045150756836,
                },
            ],
        ];
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
            ConstrainedDelaunayTriangulation::new();
        for edge in edges {
            let point_a = cdt.insert(edge[0])?;
            let point_b = cdt.insert(edge[1])?;
            cdt.cdt_sanity_check();
            cdt.add_constraint_and_split(point_a, point_b, |v| v);
            cdt.cdt_sanity_check();
        }

        Ok(())
    }

    #[test]
    fn edge_intersection_precision_test_3() -> Result<(), InsertionError> {
        let edges = [
            [
                Point2 {
                    x: -11.673287,
                    y: -28.37192,
                },
                Point2 {
                    x: -16.214716,
                    y: -43.81278,
                },
            ],
            [
                Point2 {
                    x: 7.4022045,
                    y: -51.355137,
                },
                Point2 {
                    x: -13.92232,
                    y: -36.01863,
                },
            ],
        ];

        // `f32` is important. This makes the intersection of the two edges coincide with an
        // existing vertex, triggering an edge case.
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f32>> =
            ConstrainedDelaunayTriangulation::new();
        let mut returned_constraint_edge_counts = Vec::new();
        for edge in edges {
            let point_a = cdt.insert(edge[0])?;
            let point_b = cdt.insert(edge[1])?;
            returned_constraint_edge_counts
                .push(cdt.add_constraint_and_split(point_a, point_b, |v| v).len());
            cdt.cdt_sanity_check();
        }

        // Usually, 4 constraints should be present. However, due to the overlap of the intersection
        // point, the second call to `add_constraint_and_split` does not add 2 constraint edges.
        // See issue #113 for more information
        assert_eq!(cdt.num_constraints, 3);
        assert_eq!(returned_constraint_edge_counts, vec![1, 1]);

        Ok(())
    }

    fn check_returned_edges(
        cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>,
        edges: &[FixedDirectedEdgeHandle],
        first_vertex: FixedVertexHandle,
        last_vertex: FixedVertexHandle,
    ) {
        cdt.cdt_sanity_check();

        let last = edges.last().expect("Edges cannot be empty");
        let last = cdt.directed_edge(*last);

        let mut current_from = first_vertex;

        for edge in edges {
            let edge = cdt.directed_edge(*edge);
            assert_eq!(edge.from().fix(), current_from);
            current_from = edge.to().fix();
        }

        assert_eq!(last.to().fix(), last_vertex);
    }
}
