use std::iter::FromIterator;

use crate::delaunay_core::iterators::HullIterator;
use crate::delaunay_core::InnerOuterMarker;
use crate::HintGenerator;
use crate::{delaunay_core::dcel_operations, iterators::*};
use crate::{delaunay_core::DCEL, handles::*};
use crate::{
    delaunay_core::{dcel_operations::IsolateVertexResult, math},
    HasPosition, Point2, PositionInTriangulation,
};
pub enum PositionWhenAllVerticesOnLine {
    OnEdge(FixedDirectedEdgeHandle),
    OnVertex(FixedVertexHandle),
    NotOnLine(FixedDirectedEdgeHandle),
    ExtendingLine(FixedVertexHandle),
}

pub enum InsertionResult {
    NewlyInserted(FixedVertexHandle),
    Updated(FixedVertexHandle),
}

pub struct RemovalResult<V> {
    pub removed_vertex: V,
    pub swapped_in_vertex: Option<FixedVertexHandle>,
}

/// Defines common operations on triangulations.
///
/// These operations are both available for
/// [ConstrainedDelaunayTriangulations](crate::ConstrainedDelaunayTriangulation) as well as
/// regular [DelaunayTriangulations](crate::DelaunayTriangulation).
pub trait Triangulation: Default + FromIterator<Self::Vertex> {
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
    fn s(&self) -> &DCEL<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>;

    #[doc(hidden)]
    fn s_mut(
        &mut self,
    ) -> &mut DCEL<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face>;

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

    /// Creates a new delaunay triangulation.
    ///
    /// A newly created Delaunay triangulation contains no vertices, no edges and the single
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
    /// triangulation.insert(Point2::new(0.0, 1.0));
    /// assert_eq!(triangulation.num_vertices(), 1);
    /// assert_eq!(triangulation.num_inner_faces(), 0);
    /// triangulation.insert(Point2::new(1.0, 1.0));
    /// // Two vertices define the first edge
    /// assert_eq!(triangulation.num_undirected_edges(), 1);
    /// triangulation.insert(Point2::new(1.0, 0.0));
    /// assert_eq!(triangulation.num_vertices(), 3);
    // // The third point will generate the first inner face!
    /// assert_eq!(triangulation.num_inner_faces(), 1);
    /// ```
    fn new() -> Self {
        Self::default()
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
    fn insert_with_hint(&mut self, t: Self::Vertex, hint: FixedVertexHandle) -> FixedVertexHandle {
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
    /// Returns a handle to the new vertex. Use this handle with
    /// `ConstrainedDelaunayTriangulation::vertex(..)` to refer to it.
    fn insert(&mut self, vertex: Self::Vertex) -> FixedVertexHandle {
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
    fn fixed_inner_faces(&self) -> FixedInnerFacesIterator {
        let mut result = FixedInnerFacesIterator::new(self.num_all_faces());
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
    /// The edges are returned in clockwise order as seen from any point _in_ the triangulation.
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

impl<T> TriangulationExt for T where T: Triangulation + ?Sized {}
pub trait TriangulationExt: Triangulation {
    fn insert_with_hint_option(
        &mut self,
        t: Self::Vertex,
        hint: Option<FixedVertexHandle>,
    ) -> FixedVertexHandle {
        let position = t.position();
        let result = self.insert_with_hint_option_impl(t, hint);

        self.hint_generator_mut()
            .notify_vertex_inserted(result, position);
        result
    }

    fn insert_with_hint_option_impl(
        &mut self,
        t: Self::Vertex,
        hint: Option<FixedVertexHandle>,
    ) -> FixedVertexHandle {
        use PositionInTriangulation::*;

        let insertion_result = match self.num_vertices() {
            0 => return dcel_operations::insert_first_vertex(self.s_mut(), t),
            1 => self.insert_second_vertex(t),
            _ => {
                let pos = t.position();

                if self.all_vertices_on_line() {
                    let location = self.locate_when_all_vertices_on_line(pos);
                    self.insert_when_all_vertices_on_line(location, t)
                } else {
                    let position_in_triangulation = self.locate_with_hint_option_core(pos, hint);
                    match position_in_triangulation {
                        OutsideOfConvexHull(edge) => InsertionResult::NewlyInserted(
                            self.insert_outside_of_convex_hull(edge, t),
                        ),
                        OnFace(face) => {
                            InsertionResult::NewlyInserted(self.insert_into_face(face, t))
                        }
                        OnEdge(edge) => {
                            if self.is_defined_legal(edge.as_undirected()) {
                                // If the edge is defined as legal the resulting edges must
                                // be redefined as legal
                                let handle = self.directed_edge(edge);
                                let from = handle.from().fix();
                                let to = handle.to().fix();

                                let new_handle = self.insert_on_edge(edge, t);
                                let e1 = self.get_edge_from_neighbors(from, new_handle).unwrap();
                                let e2 = self.get_edge_from_neighbors(new_handle, to).unwrap();
                                let handles = [e1.fix().as_undirected(), e2.fix().as_undirected()];

                                self.handle_legal_edge_split(handles);
                                InsertionResult::NewlyInserted(new_handle)
                            } else {
                                InsertionResult::NewlyInserted(self.insert_on_edge(edge, t))
                            }
                        }
                        OnVertex(vertex) => {
                            self.s_mut().update_vertex(vertex, t);
                            InsertionResult::Updated(vertex)
                        }
                        NoTriangulation => panic!("Error during vertex lookup. This is a bug."),
                    }
                }
            }
        };

        match insertion_result {
            InsertionResult::NewlyInserted(new_handle) => new_handle,
            InsertionResult::Updated(update_handle) => update_handle,
        }
    }

    fn locate_when_all_vertices_on_line(
        &mut self,
        position: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> PositionWhenAllVerticesOnLine {
        let mut edge = self
            .directed_edges()
            .next()
            .expect("Must not be called on a triangulations without any edge");
        let query = edge.side_query(position);
        if query.is_on_left_side() {
            return PositionWhenAllVerticesOnLine::NotOnLine(edge.fix());
        } else if query.is_on_right_side() {
            return PositionWhenAllVerticesOnLine::NotOnLine(edge.rev().fix());
        }

        assert!(query.is_on_line());
        // The target position is neither right nor left of the line. Walk along the line until
        // we walk over the position or reach one of the line's end vertices.

        let mut projection = edge.project_point(position);
        if projection.is_before_edge() {
            edge = edge.rev();
            // Make sure that the target vertex comes always after the edge. This makes sure we
            // only need to walk forwards by calling edge.next()
            projection = projection.reversed();
        }
        loop {
            if edge.from().position() == position {
                return PositionWhenAllVerticesOnLine::OnVertex(edge.from().fix());
            } else if edge.to().position() == position {
                return PositionWhenAllVerticesOnLine::OnVertex(edge.to().fix());
            }

            if projection.is_on_edge() {
                return PositionWhenAllVerticesOnLine::OnEdge(edge.fix());
            }

            edge = edge.next();
            projection = edge.project_point(position);
            if projection.is_before_edge() {
                // the target position came after the edge, now it comes before. The only
                // explanation is that we have reached on of the lines end vertices
                return PositionWhenAllVerticesOnLine::ExtendingLine(edge.from().fix());
            }
        }
    }

    fn insert_second_vertex(&mut self, vertex: Self::Vertex) -> InsertionResult {
        assert_eq!(self.num_vertices(), 1);

        let first_vertex = FixedVertexHandle::new(0);
        if self.vertex(first_vertex).position() == vertex.position() {
            self.s_mut().update_vertex(first_vertex, vertex);
            return InsertionResult::Updated(first_vertex);
        }

        let second_vertex = dcel_operations::insert_second_vertex(self.s_mut(), vertex);

        InsertionResult::NewlyInserted(second_vertex)
    }

    fn insert_when_all_vertices_on_line(
        &mut self,
        location: PositionWhenAllVerticesOnLine,
        new_vertex: Self::Vertex,
    ) -> InsertionResult {
        match location {
            PositionWhenAllVerticesOnLine::OnEdge(edge) => {
                let result = dcel_operations::split_edge_when_all_vertices_on_line(
                    self.s_mut(),
                    edge,
                    new_vertex,
                );
                InsertionResult::NewlyInserted(result)
            }
            PositionWhenAllVerticesOnLine::OnVertex(vertex) => InsertionResult::Updated(vertex),
            PositionWhenAllVerticesOnLine::NotOnLine(edge) => {
                let result = self.insert_outside_of_convex_hull(edge, new_vertex);
                InsertionResult::NewlyInserted(result)
            }
            PositionWhenAllVerticesOnLine::ExtendingLine(end_vertex) => {
                let result = dcel_operations::extend_line(self.s_mut(), end_vertex, new_vertex);
                InsertionResult::NewlyInserted(result)
            }
        }
    }

    fn locate_with_hint_option_core(
        &self,
        point: Point2<<Self::Vertex as HasPosition>::Scalar>,
        hint: Option<FixedVertexHandle>,
    ) -> PositionInTriangulation {
        let start = hint.unwrap_or_else(|| self.hint_generator().get_hint(point));
        self.locate_with_hint_fixed_core(point, start)
    }

    fn insert_outside_of_convex_hull(
        &mut self,
        convex_hull_edge: FixedDirectedEdgeHandle,
        new_vertex: Self::Vertex,
    ) -> FixedVertexHandle {
        let position = new_vertex.position();

        let edges_to_connect: Vec<_> = self
            .get_vertex_facing_edges(convex_hull_edge, new_vertex.position())
            .into();

        let result =
            dcel_operations::connect_edge_strip(self.s_mut(), &edges_to_connect, new_vertex);
        self.legalize_edges_after_insertion(edges_to_connect, position);
        result
    }

    fn get_vertex_facing_edges(
        &self,
        start_edge: FixedDirectedEdgeHandle,
        position: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> std::collections::VecDeque<FixedDirectedEdgeHandle> {
        let mut result = std::collections::VecDeque::with_capacity(8);
        let mut current_edge_forward = self.directed_edge(start_edge);

        debug_assert!(current_edge_forward.side_query(position).is_on_left_side());

        result.push_back(current_edge_forward.fix());

        loop {
            current_edge_forward = current_edge_forward.next();
            let cur_query = current_edge_forward.side_query(position);
            if cur_query.is_on_left_side() {
                result.push_back(current_edge_forward.fix());
            } else {
                break;
            }
        }

        let mut current_edge_backward = self.directed_edge(start_edge);
        loop {
            current_edge_backward = current_edge_backward.prev();
            let cur_query = current_edge_backward.side_query(position);
            if cur_query.is_on_left_side() {
                result.push_front(current_edge_backward.fix());
            } else {
                break;
            }
        }
        result
    }

    fn insert_into_face(
        &mut self,
        face: FixedFaceHandle<InnerTag>,
        t: Self::Vertex,
    ) -> FixedVertexHandle {
        let new_handle = dcel_operations::insert_into_triangle(self.s_mut(), t, face);
        self.legalize_vertex(new_handle);
        new_handle
    }

    fn insert_on_edge(
        &mut self,
        edge: FixedDirectedEdgeHandle,
        new_vertex: Self::Vertex,
    ) -> FixedVertexHandle {
        let edge_handle = self.directed_edge(edge);
        let new_vertex = if edge_handle.is_outer_edge() {
            dcel_operations::split_half_edge(self.s_mut(), edge.rev(), new_vertex)
        } else if edge_handle.rev().is_outer_edge() {
            dcel_operations::split_half_edge(self.s_mut(), edge, new_vertex)
        } else {
            dcel_operations::split_edge(self.s_mut(), edge, new_vertex)
        };
        self.legalize_vertex(new_vertex);
        new_vertex
    }

    fn legalize_vertex(&mut self, new_handle: FixedVertexHandle) {
        let position = self.vertex(new_handle).position();
        let edges = self
            .vertex(new_handle)
            .out_edges()
            .filter(|e| !e.is_outer_edge())
            .map(|edge| edge.next().fix())
            .collect();
        self.legalize_edges_after_insertion(edges, position);
    }

    /// The Delaunay property refers to the property that no point lies inside
    /// the circumcircle of any of the triangulation's triangles. Adding a
    /// new point into the triangulations may violate this property, this method
    /// "repairs" it by strategically flipping edges until the property
    /// holds again. Every flip produces more "illegal" edges that may have to
    /// be flipped. However, since any edge is flipped at most once, this
    /// algorithm is known to terminate.
    ///
    /// The given position is the position of a new point may have
    /// invalidated the Delaunay property, the point must be on the left side
    /// of the given edge.
    ///
    /// "Flipping an edge" refers to switching to the other diagonal in a
    /// four sided polygon.
    fn legalize_edges_after_insertion(
        &mut self,
        mut edges: Vec<FixedDirectedEdgeHandle>,
        position: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) {
        while let Some(e) = edges.pop() {
            if !self.is_defined_legal(e.as_undirected()) {
                let edge = self.directed_edge(e);

                //         v2------edge.from()
                //          |     / |
                //          |    /  |
                //          |   /<-edge that might be flipped ("edge")
                //          |  /    |
                //          | V     |
                //  edge.to()-------position of newly inserted point
                //
                //
                // If the edge is flipped, the new quad will look like this:
                //
                //      New illegal edge
                //              |
                //              V
                //         v2-------+
                //          | \     |
                // New      |  \    |
                // illegal->|   \   |
                // edge     |    \  |
                //          |     \ |
                //          +-------position of newly inserted point
                let v2 = edge.rev().opposite_position();
                if let Some(v2) = v2 {
                    let v1 = edge.to().position();
                    let v0 = edge.from().position();
                    debug_assert!(math::is_ordered_ccw(v2, v1, v0));
                    let should_flip = math::contained_in_circumference(v2, v1, v0, position);

                    if should_flip {
                        let e1 = edge.rev().next().fix();
                        let e2 = edge.rev().prev().fix();

                        dcel_operations::flip_cw(self.s_mut(), e.as_undirected());
                        edges.push(e1);
                        edges.push(e2);
                    }
                }
            }
        }
    }

    fn validate_vertex_handle(&self, handle: FixedVertexHandle) -> FixedVertexHandle {
        if handle.index() < self.num_vertices() {
            handle
        } else {
            FixedVertexHandle::new(0)
        }
    }

    fn walk_to_nearest_neighbor(
        &self,
        start: FixedVertexHandle,
        position: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> VertexHandle<Self::Vertex, Self::DirectedEdge, Self::UndirectedEdge, Self::Face> {
        let start_position = self.vertex(start).position();

        if start_position == position {
            return self.vertex(start);
        }

        let mut current_minimal_distance = position.distance2(start_position);
        let mut current_minimum_vertex = self.vertex(start);

        while let Some((next_minimum_index, next_minimal_distance)) = current_minimum_vertex
            .out_edges()
            .filter_map(|out_edge| {
                let next_candidate = out_edge.to();
                let new_distance = next_candidate.position().distance2(position);
                if new_distance < current_minimal_distance {
                    Some((next_candidate, new_distance))
                } else {
                    None
                }
            })
            .next()
        {
            current_minimal_distance = next_minimal_distance;
            current_minimum_vertex = next_minimum_index;
        }

        current_minimum_vertex
    }

    /// "Walks" through the triangulation until it finds the target point.
    fn locate_with_hint_fixed_core(
        &self,
        target_position: Point2<<Self::Vertex as HasPosition>::Scalar>,
        start: FixedVertexHandle,
    ) -> PositionInTriangulation {
        if self.num_vertices() < 2 {
            return match self.vertices().next() {
                Some(single_vertex) if single_vertex.position() == target_position => {
                    PositionInTriangulation::OnVertex(single_vertex.fix())
                }
                _ => PositionInTriangulation::NoTriangulation,
            };
        }

        let start = self.validate_vertex_handle(start);

        //let closest_vertex = self.walk_to_nearest_neighbor(start, target_position);
        let closest_vertex = self.vertex(start);

        let out_edge = closest_vertex
            .out_edge()
            .expect("No out edge found. This is a bug.");

        let mut query = out_edge.side_query(target_position);
        let mut edge = if query.is_on_right_side() {
            out_edge.rev()
        } else {
            out_edge
        };

        // Invariant: position is on the left side or on the line of edge.
        loop {
            let prev = edge.prev();
            if !edge.is_outer_edge() {
                let prev_query = prev.side_query(target_position);
                if prev_query.is_on_right_side() {
                    edge = prev.rev();
                    query = edge.side_query(target_position);
                    continue;
                }

                let next = edge.next();
                let next_query = next.side_query(target_position);

                if next_query.is_on_right_side() {
                    edge = next.rev();
                    query = edge.side_query(target_position);
                    continue;
                }

                self.hint_generator()
                    .notify_vertex_lookup(edge.from().fix());

                // Point must be in triangle or on its lines
                return match (
                    query.is_on_line(),
                    next_query.is_on_line(),
                    prev_query.is_on_line(),
                ) {
                    // Point lies on no line and must be inside the face
                    (false, false, false) => {
                        PositionInTriangulation::OnFace(next.face().fix().adjust_inner_outer())
                    }
                    // Point lies on exactly one line
                    (false, false, true) => PositionInTriangulation::OnEdge(prev.fix()),
                    (false, true, false) => PositionInTriangulation::OnEdge(next.fix()),
                    (true, false, false) => PositionInTriangulation::OnEdge(edge.fix()),
                    // Point lies on exactly two lines. Since the edges cannot be collinear,
                    // the point lies on the intersection
                    (false, true, true) => PositionInTriangulation::OnVertex(prev.from().fix()),
                    (true, false, true) => PositionInTriangulation::OnVertex(edge.from().fix()),
                    (true, true, false) => PositionInTriangulation::OnVertex(next.from().fix()),
                    (true, true, true) => panic!("Invalid triangle. This is a bug"),
                };
            } else {
                // Edge is part of the convex hull
                if query.is_on_line() {
                    let rev = edge.rev();
                    if rev.is_outer_edge() {
                        // Both edge and its rev are part of the convex hull,
                        // hence all points must lie on a line
                        return self.locate_if_all_points_on_line(rev, target_position);
                    }
                    // Triangulation is not degenerated. Just continue with the
                    // inner triangle
                    edge = rev;
                    continue;
                } else {
                    self.hint_generator()
                        .notify_vertex_lookup(edge.from().fix());
                    return PositionInTriangulation::OutsideOfConvexHull(edge.fix());
                }
            }
        }
    }

    fn locate_if_all_points_on_line(
        &self,
        start_edge: DirectedEdgeHandle<
            Self::Vertex,
            Self::DirectedEdge,
            Self::UndirectedEdge,
            Self::Face,
        >,
        position: Point2<<Self::Vertex as HasPosition>::Scalar>,
    ) -> PositionInTriangulation {
        let mut current_edge = start_edge;

        let current_projection = start_edge.project_point(position);

        if current_projection.is_before_edge() {
            current_edge = current_edge.rev();
        }

        loop {
            if current_edge.from().position() == position {
                return PositionInTriangulation::OnVertex(current_edge.from().fix());
            }
            if current_edge.to().position() == position {
                return PositionInTriangulation::OnVertex(current_edge.to().fix());
            }

            let current_projection = current_edge.project_point(position);

            if current_projection.is_after_edge() {
                if current_edge.next() == current_edge.rev() {
                    return PositionInTriangulation::OutsideOfConvexHull(current_edge.fix());
                }
                current_edge = current_edge.next();
            } else {
                // 0.0 <= current_projection <= 1.0
                return PositionInTriangulation::OnEdge(current_edge.fix());
            }
        }
    }

    fn remove_and_notify(&mut self, vertex_to_remove: FixedVertexHandle) -> Self::Vertex {
        let position = self.vertex(vertex_to_remove).position();
        let removal_result = self.remove_core(vertex_to_remove);

        let swapped_in_point = removal_result
            .swapped_in_vertex
            .map(|_| self.vertex(vertex_to_remove).position());

        self.hint_generator_mut().notify_vertex_removed(
            swapped_in_point,
            vertex_to_remove,
            position,
        );

        return removal_result.removed_vertex;
    }

    fn remove_core(&mut self, vertex_to_remove: FixedVertexHandle) -> RemovalResult<Self::Vertex> {
        if self.num_all_faces() <= 1 {
            return dcel_operations::remove_when_degenerate(self.s_mut(), vertex_to_remove);
        }
        let vertex = self.vertex(vertex_to_remove);

        let mut border_loop = Vec::with_capacity(10);
        let mut convex_hull_edge = None;
        for edge in vertex.out_edges().rev() {
            if edge.is_outer_edge() {
                convex_hull_edge = Some(edge.fix());
                break;
            }
            let border_edge = edge.next();
            border_loop.push(border_edge.fix());
        }

        if let Some(convex_hull_edge) = convex_hull_edge {
            let mut isolation_result = self.isolate_convex_hull_vertex(convex_hull_edge);

            dcel_operations::cleanup_isolated_vertex(self.s_mut(), &mut isolation_result);
            dcel_operations::swap_remove_vertex(self.s_mut(), vertex_to_remove)
        } else {
            let mut isolation_result = dcel_operations::isolate_vertex_and_fill_hole(
                self.s_mut(),
                border_loop,
                vertex_to_remove,
            );
            // Not exactly elegant. IsolateVertexResult should maybe be split into two parts
            let new_edges = std::mem::replace(&mut isolation_result.new_edges, Vec::new());
            self.legalize_edges_after_removal(new_edges, |edge| {
                !isolation_result.is_new_edge(edge)
            });
            dcel_operations::cleanup_isolated_vertex(self.s_mut(), &mut isolation_result);
            dcel_operations::swap_remove_vertex(self.s_mut(), vertex_to_remove)
        }
    }

    fn isolate_convex_hull_vertex(
        &mut self,
        convex_hull_out_edge: FixedDirectedEdgeHandle,
    ) -> IsolateVertexResult {
        let mut edges_to_validate = Vec::with_capacity(10);
        let mut convex_edges: Vec<FixedDirectedEdgeHandle> = Vec::with_capacity(10);

        let loop_end = self.directed_edge(convex_hull_out_edge);
        let loop_start = loop_end.ccw().fix();
        let mut current = loop_start;
        let loop_end_next = loop_end.next().fix();
        let loop_end = loop_end.fix();

        loop {
            let current_handle = self.directed_edge(current);
            current = current_handle.ccw().fix();
            let edge = current_handle.next().fix();
            convex_edges.push(edge);

            loop {
                if let &[.., edge1, edge2] = &*convex_edges {
                    let edge1 = self.directed_edge(edge1);
                    let edge2 = self.directed_edge(edge2);

                    let target_position = edge2.to().position();
                    // Check if the new edge would violate the convex hull property by going
                    // "to the left". The convex hull must only contain curves to the left
                    if edge1.side_query(target_position).is_on_left_side() {
                        // Violation detected. It is resolved by flipping the edge that connects
                        // the most recently added edge (edge2) with the point that was removed
                        let edge_to_flip = edge2.prev().fix().rev();
                        dcel_operations::flip_cw(self.s_mut(), edge_to_flip.as_undirected());
                        convex_edges.pop();
                        convex_edges.pop();
                        convex_edges.push(edge_to_flip);
                        edges_to_validate.push(edge_to_flip.as_undirected());
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }

            if current == loop_end {
                break;
            }
        }

        convex_edges.push(loop_end_next);
        let result = dcel_operations::disconnect_edge_strip(self.s_mut(), convex_edges);
        self.legalize_edges_after_removal(edges_to_validate, |_| false);
        result
    }

    /// After a vertex removal, the created hole is stitched together by connecting
    /// all vertices to a single vertex at the border (triangle fan).
    /// Since these new edges can violate the Delaunay property, it must be restored.
    ///
    /// The algorithm works like this:
    ///  - Add all new edges to an "invalid" list
    ///  - While the invalid list is not empty: Determine if flipping the top edge is
    ///    required to restore the Delaunay property locally.
    ///  - If the edge was flipped: Determine the flip polygon. A flip refers
    ///    to switching the diagonal in a four sided polygon which defines the
    ///    flip polygon.
    ///  - Add all edges of the flip polygon to the invalid list if they were
    ///    newly created. Otherwise, the edge is part of the border loop surrounding
    ///    the hole created after the vertex removal. These are known to be valid and
    ///    need not to be checked
    ///
    /// For more details, refer to
    /// Olivier Devillers. Vertex Removal in Two Dimensional Delaunay Triangulation:
    /// Speed-up by Low Degrees Optimization.
    /// https://doi.org/10.1016/j.comgeo.2010.10.001
    ///
    /// Note that the described low degrees optimization is not yet part of this library.
    fn legalize_edges_after_removal<F>(
        &mut self,
        mut edges_to_validate: Vec<FixedUndirectedEdgeHandle>,
        edge_must_not_be_flipped_predicate: F,
    ) where
        F: Fn(FixedUndirectedEdgeHandle) -> bool,
    {
        while let Some(next_edge) = edges_to_validate.pop() {
            // left----to
            //  |     ^ |
            //  |    /  |
            //  |   /<-edge that might be flipped ("next_edge")
            //  |  /    |
            //  | /     |
            // from----right
            //
            // left, from, right and to define the flip polygon.
            if self.is_defined_legal(next_edge) || edge_must_not_be_flipped_predicate(next_edge) {
                continue;
            }

            let edge = self.directed_edge(next_edge.as_directed());
            let e2 = edge.prev();
            let e4 = edge.rev().prev();

            let from = edge.from().position();
            let to = edge.to().position();
            let left = edge.opposite_position();
            let right = edge.rev().opposite_position();

            let should_flip = match (left, right) {
                (Some(left), Some(right)) => {
                    math::contained_in_circumference(from, to, left, right)
                }
                // Handle special cases when evaluating edges next to the convex hull
                (None, Some(right)) => math::is_ordered_ccw(right, from, to),
                (Some(left), None) => math::is_ordered_ccw(left, to, from),
                (None, None) => {
                    panic!("Unexpected geometry. This is a bug in spade.")
                }
            };

            if should_flip {
                let e1 = edge.next();
                let e3 = edge.rev().next();

                edges_to_validate.push(e1.fix().as_undirected());
                edges_to_validate.push(e2.fix().as_undirected());
                edges_to_validate.push(e3.fix().as_undirected());
                edges_to_validate.push(e4.fix().as_undirected());

                let fixed = edge.fix();
                dcel_operations::flip_cw(self.s_mut(), fixed.as_undirected());
            }
        }
    }

    #[cfg(test)]
    fn sanity_check(&self) {
        self.s().sanity_check();
        let all_vertices_on_line = self.s().num_faces() <= 1;

        for face in self.s().inner_faces() {
            let triangle = face.vertices();
            // Check that all vertices are stored in ccw orientation
            assert!(math::side_query(
                triangle[0].position(),
                triangle[1].position(),
                triangle[2].position()
            )
            .is_on_left_side(),);
        }

        for edge in self.s().directed_edges() {
            if all_vertices_on_line {
                assert_eq!(edge.face().fix(), dcel_operations::OUTER_FACE_HANDLE)
            } else {
                assert_ne!(edge.face(), edge.rev().face());
            }
            assert_ne!(edge.from(), edge.to());
        }

        if all_vertices_on_line {
            if self.s().num_vertices() > 1 {
                assert_eq!(self.s().num_undirected_edges(), self.s().num_vertices() - 1);
            } else {
                assert_eq!(self.s().num_undirected_edges(), 0);
            }
            assert_eq!(self.s().num_faces(), 1);
        } else {
            let num_inner_edges = self
                .s()
                .directed_edges()
                .filter(|e| !e.face().is_outer())
                .count();

            let num_inner_faces = self.s().num_faces() - 1;
            assert_eq!(num_inner_faces * 3, num_inner_edges);
        }
    }
}
