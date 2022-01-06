use crate::delaunay_core::iterators::HullIterator;
use crate::delaunay_core::InnerOuterMarker;
use crate::iterators::*;
use crate::HintGenerator;
use crate::{delaunay_core::Dcel, handles::*};
use crate::{HasPosition, InsertionError, Point2, PositionInTriangulation, TriangulationExt};

/// Defines common operations on triangulations.
///
/// These operations are both available for
/// [ConstrainedDelaunayTriangulations](crate::ConstrainedDelaunayTriangulation) as well as
/// regular [DelaunayTriangulations](crate::DelaunayTriangulation).
pub trait Triangulation: Default {
    /// The triangulation's vertex type.
    type Vertex: HasPosition;
    /// The triangulation's edge type. Any new edge is created by using the `Default` trait.
    type DirectedEdge: Default;
    /// The triangulation's undirected edge type. Any new edge is created by using the `Default` trait.
    type UndirectedEdge: Default;
    /// The triangulation's face type. Any new face is created by using the `Default` trait.
    type Face: Default;

    /// The hint generator used by the triangulation. See [HintGenerator] for more information.
    type HintGenerator: HintGenerator<<Self::Vertex as HasPosition>::Scalar>;

    #[doc(hidden)]
    fn s(&self) -> &Dcel<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>;

    #[doc(hidden)]
    fn s_mut(
        &mut self,
    ) -> &mut Dcel<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>;

    #[doc(hidden)]
    fn is_defined_legal(&self, _: FixedUndirectedEdgeHandle) -> bool {
        false
    }

    #[doc(hidden)]
    fn handle_legal_edge_split(&mut self, _: [FixedUndirectedEdgeHandle; 2]) {}

    #[doc(hidden)]
    fn hint_generator(&self) -> &Self::HintGenerator;

    #[doc(hidden)]
    fn hint_generator_mut(&mut self) -> &mut Self::HintGenerator;

    /// Creates a new triangulation.
    ///
    /// A newly created triangulation contains no vertices, no edges and the single
    /// outer face.
    ///
    /// # Example
    /// ```
    /// use spade::{DelaunayTriangulation, Point2, Triangulation};
    ///
    /// let mut triangulation: DelaunayTriangulation<Point2<f64>> = DelaunayTriangulation::new();
    /// // An empty triangulation has no vertices and one face
    /// assert_eq!(triangulation.num_vertices(), 0);
    /// assert_eq!(triangulation.num_all_faces(), 1);
    /// assert_eq!(triangulation.num_inner_faces(), 0); // This count excludes the outer face
    /// assert_eq!(triangulation.num_directed_edges(), 0);
    ///
    /// triangulation.insert(Point2::new(0.0, 1.0));
    /// assert_eq!(triangulation.num_vertices(), 1);
    /// assert_eq!(triangulation.num_inner_faces(), 0);
    ///
    /// triangulation.insert(Point2::new(1.0, 1.0));
    /// // Two vertices define the first edge
    /// assert_eq!(triangulation.num_undirected_edges(), 1);
    ///
    /// triangulation.insert(Point2::new(1.0, 0.0));
    /// assert_eq!(triangulation.num_vertices(), 3);
    // // The third point will generate the first inner face!
    /// assert_eq!(triangulation.num_inner_faces(), 1);
    /// ```
    fn new() -> Self {
        Self::default()
    }

    /// Creates a new triangulation and pre-allocates some space for vertices, edges and faces
    fn with_capacity(num_vertices: usize, num_undirected_edges: usize, num_faces: usize) -> Self {
        let mut result = Self::new();
        result
            .s_mut()
            .reserve_capacity(num_vertices, num_undirected_edges, num_faces);
        result
    }

    /// Creates a new triangulation populated with some vertices.
    ///
    /// This will usually be more efficient than inserting the elements sequentially by calling
    /// [insert](Triangulation::insert).
    ///
    /// Returns an [InsertionError] if any input coordinate is invalid. This method should never fail
    /// if all vertices were successfully checked with [crate::validate_vertex].
    ///
    /// # Runtime
    /// This method has a run time of `O(n)` but will run near linearly in practice.
    /// The runtime can be as worse as `O(n²)` if the inputs are very degenerate, e.g.
    /// if all input vertices lie on the same line.
    fn bulk_load(elements: Vec<Self::Vertex>) -> Result<Self, InsertionError> {
        let mut result: Self = crate::delaunay_core::bulk_load(elements)?;
        let hint_generator = Self::HintGenerator::initialize_from_triangulation(&result);
        *result.hint_generator_mut() = hint_generator;
        Ok(result)
    }

    /// Converts a fixed vertex handle to a reference vertex handle.
    ///
    /// *See also the [handles](crate::handles) module for more information.*
    fn vertex(
        &self,
        handle: FixedVertexHandle,
    ) -> VertexHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        self.s().vertex(handle)
    }

    /// Returns a mutable reference to the associated data of a vertex.
    fn vertex_data_mut(&mut self, handle: FixedVertexHandle) -> &mut Self::Vertex {
        self.s_mut().vertex_data_mut(handle)
    }

    /// Converts a fixed face handle to a reference face handle.
    ///
    /// *See also the [handles](crate::handles) module for more information.*
    fn face<InnerOuter: InnerOuterMarker>(
        &self,
        handle: FixedFaceHandle<InnerOuter>,
    ) -> FaceHandle<InnerOuter, Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        self.s().face(handle)
    }

    /// Returns a reference handle to the single outer face of this triangulation.
    fn outer_face(
        &self,
    ) -> FaceHandle<
        PossiblyOuterTag,
        Self::Vertex,
        Self::DirectedEdge,
        Self::UndirectedEdge,
        Self::Face,
    > {
        self.s().outer_face()
    }

    /// Converts a fixed directed edge handle handle to a reference directed edge handle.
    ///
    /// *See also the [handles](crate::handles) module for more information.*
    fn directed_edge(
        &self,
        handle: FixedDirectedEdgeHandle,
    ) -> DirectedEdgeHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        DirectedEdgeHandle::new(self.s(), handle)
    }

    /// Converts a fixed undirected edge handle to a reference undirected edge handle.
    ///
    /// *See also the [handles](crate::handles) module for more information.*
    fn undirected_edge(
        &self,
        handle: FixedUndirectedEdgeHandle,
    ) -> UndirectedEdgeHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        UndirectedEdgeHandle::new(self.s(), handle)
    }

    /// Returns a mutable reference ot the associated data of an undirected edge.
    fn undirected_edge_data_mut(
        &mut self,
        handle: FixedUndirectedEdgeHandle,
    ) -> &mut Self::UndirectedEdge {
        self.s_mut().undirected_edge_data_mut(handle)
    }

    /// Returns the number of all faces, including the single outer face, of this
    /// triangulation.
    ///
    /// This is always equal to `triangulation.num_inner_faces() + 1`.
    fn num_all_faces(&self) -> usize {
        self.s().num_faces()
    }

    /// Returns the number of inner faces in this triangulation.
    fn num_inner_faces(&self) -> usize {
        self.s().num_faces() - 1
    }

    /// Returns the number of undirected edges in this triangulation.
    fn num_undirected_edges(&self) -> usize {
        self.s().num_undirected_edges()
    }

    /// Returns the number of directed edges in this triangulation.
    fn num_directed_edges(&self) -> usize {
        self.s().num_directed_edges()
    }

    /// Returns the number of edges of the convex hull.
    ///
    /// *See also [convex_hull](Self::convex_hull)*
    ///
    /// # Complexity
    /// This method does not need to iterate through the convex hull and has a complexity of O(1)
    fn convex_hull_size(&self) -> usize {
        if self.all_vertices_on_line() {
            self.num_directed_edges()
        } else {
            let num_inner_edges = self.num_inner_faces() * 3;
            self.num_directed_edges() - num_inner_edges
        }
    }

    /// An iterator visiting all directed edges.
    ///
    /// The iterator type is [DirectedEdgeHandle].
    fn directed_edges(
        &self,
    ) -> DirectedEdgeIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        self.s().directed_edges()
    }

    /// An iterator over all undirected edges.
    ///
    /// The iterator type is [UndirectedEdgeHandle]
    fn undirected_edges(
        &self,
    ) -> UndirectedEdgeIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        self.s().undirected_edges()
    }

    /// Returns the number vertices in this triangulation.
    fn num_vertices(&self) -> usize {
        self.s().num_vertices()
    }

    /// An iterator visiting all vertices.
    ///
    /// The iterator type is [VertexHandle](handles/type.VertexHandle.html).
    fn vertices(
        &self,
    ) -> VertexIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        self.s().vertices()
    }

    /// An iterator visiting all vertices.
    ///
    /// The iterator type is [FixedVertexHandle](handles/struct.FixedVertexHandle.html).
    fn fixed_vertices(&self) -> FixedVertexIterator {
        self.s().fixed_vertices()
    }

    /// An iterator visiting all faces.
    ///
    /// The first returned face is the outer face, all other faces will be inner faces.
    /// The iterator type is [FaceHandle<PossiblyOuterTag, ...>](FaceHandle).
    /// See also [inner_faces()](Self::inner_faces()).
    fn all_faces(
        &self,
    ) -> FaceIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        self.s().faces()
    }

    /// An iterator visiting all inner faces.
    //
    /// The iterator type is [FaceHandle<InnerTag, ...>](FaceHandle).
    fn inner_faces(
        &self,
    ) -> InnerFaceIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        self.s().inner_faces()
    }

    /// An iterator visiting all faces of the Voronoi diagram.
    ///
    /// The iterator type is [VoronoiFace]
    fn voronoi_faces(
        &self,
    ) -> VoronoiFaceIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>
    {
        VoronoiFaceIterator::new(self.s())
    }

    /// An iterator visiting all directed voronoi edges.
    ///
    /// The iterator type is (DirectedVoronoiEdge)[crate::handles::DirectedVoronoiEdge]
    fn directed_voronoi_edges(
        &self,
    ) -> DirectedVoronoiEdgeIterator<
        Self::Vertex,
        Self::DirectedEdge,
        Self::UndirectedEdge,
        Self::Face,
    > {
        DirectedVoronoiEdgeIterator::new(self.s())
    }

    /// An iterator visiting all undirected voronoi edges.
    ///
    /// The iterator type is (UndirectedVoronoiEdge)[crate::handles::UndirectedVoronoiEdge]
    fn undirected_voronoi_edges(
        &self,
    ) -> UndirectedVoronoiEdgeIterator<
        Self::Vertex,
        Self::DirectedEdge,
        Self::UndirectedEdge,
        Self::Face,
    > {
        UndirectedVoronoiEdgeIterator::new(self.s())
    }

    /// Returns information about the location of a point in a triangulation.
    fn locate(
        &self,
        point: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> PositionInTriangulation {
        let hint = self.hint_generator().get_hint(point);
        self.locate_with_hint_option_core(point, Some(hint))
    }

    /// Locates a vertex at a given position.
    ///
    /// Returns `None` if the point could not be found.
    #[allow(clippy::type_complexity)]
    fn locate_vertex(
        &self,
        point: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> Option<VertexHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>>
    {
        match self.locate(point) {
            PositionInTriangulation::OnVertex(vertex) => Some(self.vertex(vertex)),
            _ => None,
        }
    }

    /// Returns an edge between two vertices.
    ///
    /// If the edge does not exist, `None` is returned.
    /// This operation runs in `O(n)` time, where `n` is
    /// the degree of `from`.
    #[allow(clippy::type_complexity)]
    fn get_edge_from_neighbors(
        &self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
    ) -> Option<
        DirectedEdgeHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>,
    > {
        self.s().get_edge_from_neighbors(from, to)
    }

    /// Returns information about the location of a point in a triangulation.
    ///
    /// Additionally, a hint can be given to speed up computation.
    /// The hint should be a vertex close to the position that
    /// is being looked up.
    fn locate_with_hint(
        &self,
        point: Point2<<Self::Vertex as HasPosition>::Scalar>,
        hint: FixedVertexHandle,
    ) -> PositionInTriangulation {
        self.locate_with_hint_option_core(point, Some(hint))
    }

    /// Inserts a new vertex into the triangulation.
    ///
    /// A hint can be given to speed up the process.
    /// The hint should be a handle of a vertex close to the new vertex.
    /// in this case the insertion time can be reduced to O(1) on average
    /// if the hint is close. If the hint is randomized, running time will
    /// be O(sqrt(n)) on average with an O(n) worst case.
    fn insert_with_hint(
        &mut self,
        t: Self::Vertex,
        hint: FixedVertexHandle,
    ) -> Result<FixedVertexHandle, InsertionError> {
        self.insert_with_hint_option(t, Some(hint))
    }

    /// Attempts to remove a vertex from the triangulation.
    ///
    /// Returns the removed vertex data if it could be found.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles
    /// upon successful removal.
    fn locate_and_remove(
        &mut self,
        point: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> Option<Self::Vertex> {
        match self.locate_with_hint_option_core(point, None) {
            PositionInTriangulation::OnVertex(handle) => Some(self.remove_and_notify(handle)),
            _ => None,
        }
    }

    /// Removes a vertex from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles.
    fn remove(&mut self, vertex: FixedVertexHandle) -> Self::Vertex {
        self.remove_and_notify(vertex)
    }

    /// Inserts a new vertex into the triangulation.
    ///
    /// This operation runs in O(log(n)) on average when using a tree
    /// lookup to back up the triangulation, or in O(sqrt(n)) when using
    /// a walk lookup. n denotes the number of vertices, the given
    /// running times assume that input data is given uniformly randomly
    /// distributed. If the point has already been contained in the
    /// triangulation, the old vertex is overwritten.
    ///
    /// Returns either a handle to the new vertex or an error if the vertex could not be inserted.
    /// The triangulation will remain unchanged if an error ocurred.
    ///
    /// Use [vertex](Self::vertex) to retrieve more information about the inserted vertex.
    fn insert(&mut self, vertex: Self::Vertex) -> Result<FixedVertexHandle, InsertionError> {
        self.insert_with_hint_option(vertex, None)
    }

    /// An iterator visiting all undirected edges.
    ///
    /// The iterator type is [FixedUndirectedEdgeHandle](handles/struct.FixedUndirectedEdgeHandle.html).
    fn fixed_undirected_edges(&self) -> FixedUndirectedEdgeIterator {
        FixedUndirectedEdgeIterator::new(self.num_undirected_edges())
    }

    /// An iterator visiting all directed edges.
    ///
    /// The iterator type is [FixedDirectedEdgeHandle](handles/struct.FixedDirectedEdgeHandle.html).
    fn fixed_directed_edges(&self) -> FixedDirectedEdgeIterator {
        FixedDirectedEdgeIterator::new(self.num_directed_edges())
    }

    /// An iterator visiting all faces.
    ///
    /// The first returned element is the outer face. All other
    /// faces are inner faces.
    ///
    /// The iterator type is [FixedFaceHandle<PossiblyOuterTag, ...>](FixedFaceHandle).
    fn fixed_all_faces(&self) -> FixedFaceIterator {
        FixedFaceIterator::new(self.num_all_faces())
    }

    /// An iterator visiting all inner faces of the triangulation.
    ///
    /// The iterator type is [FixedFaceHandle<InnerTag, ...>](FixedFaceHandle).
    fn fixed_inner_faces(&self) -> FixedInnerFaceIterator {
        let mut result = FixedInnerFaceIterator::new(self.num_all_faces());
        result.next();
        result
    }

    /// Returns `true` if all vertices lie on a single line.
    ///
    /// This is always the case for triangulations with 0, 1 or two vertices.
    fn all_vertices_on_line(&self) -> bool {
        self.num_all_faces() == 1
    }

    /// Returns an iterator over all convex hull edges.
    ///
    /// The edges are returned in clockwise order as seen from any point in the triangulation.
    ///
    #[doc = include_str!("../images/convex_hull_scenario.svg")]
    ///
    /// *A triangulation with its convex hull being highlighted. `e0` .. `e5` denote the returned
    /// edges in clockwise order.*
    ///
    /// *See also [convex_hull_size](Self::convex_hull_size)*
    fn convex_hull(
        &self,
    ) -> HullIterator<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        {
            HullIterator::new(self.s())
        }
    }

    /// Returns a mutable reference to the associated data of a face.
    fn face_data_mut<InnerOuter: InnerOuterMarker>(
        &mut self,
        handle: FixedFaceHandle<InnerOuter>,
    ) -> &mut Self::Face {
        self.s_mut().face_data_mut(handle)
    }

    /// Returns a mutable reference to the associated data of a directed edge.
    fn directed_edge_data_mut(
        &mut self,
        handle: FixedDirectedEdgeHandle,
    ) -> &mut Self::DirectedEdge {
        self.s_mut().directed_edge_data_mut(handle)
    }
}
