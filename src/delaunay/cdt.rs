use delaunay::*;
use self::delaunay_basic::{BasicDelaunaySubdivision, HasSubdivision};
use self::dcel::*;
use self::line_intersection_iterator::*;
use traits::{HasPosition, HasPosition2D};
use std::marker::PhantomData;
use point_traits::{PointN, TwoDimensional};
use kernels::{DelaunayKernel, FloatKernel};
use primitives::SimpleEdge;

/// Type shorthand for a constrained Delaunay triangulation using
/// the precise `FloatKernel`.
pub type FloatCDT<T, L> = ConstrainedDelaunayTriangulation<T, FloatKernel, L>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdtEdge(bool);

impl CdtEdge {
    fn is_constraint_edge(&self) -> bool {
        self.0
    }

    fn make_constraint_edge(&mut self) {
        assert!(!self.is_constraint_edge());
        self.0 = true;
    }
}


impl Default for CdtEdge {
    fn default() -> Self {
        CdtEdge(false)
    }
}

/// A two dimensional constrained Delaunay triangulation.
///
/// A constrained Delaunay triangulation is a triangulation that
/// can contain _constraint edges_. These edges will be present
/// in the resulting triangulation, the resulting triangulation
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
/// use spade::delaunay::FloatCDT;
/// let mut cdt = FloatCDT::with_walk_locate();
/// let v0 = cdt.insert([0.0, 0.0f32]);
/// let v1 = cdt.insert([1.0, 0.0]);
/// cdt.add_constraint(v0, v1);
/// // Alternatively, consider using this shorthand
/// cdt.add_constraint_edge([1.0, 1.0], [1.0, 0.0]);
/// // This should print "2"
/// println!("Number of constraints: {}", cdt.num_constraints());
/// // Constraints are bidirectional!
/// assert!(cdt.exists_constraint(v1, v0));
/// assert!(cdt.exists_constraint(v0, v1));
/// // Check if a new constraint could be added
/// let from = [1.0, -2.0];
/// let to = [1.0, 0.0];
/// if !cdt.intersects_constraint(&from, &to) {
///   // No intersections, the edge can be added
///   cdt.add_constraint_edge(from, to);
/// }
/// ```
#[derive(Clone)]
pub struct ConstrainedDelaunayTriangulation<V, K, L = DelaunayTreeLocate<<V as HasPosition>::Point>>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>,
          L: DelaunayLocateStructure<V::Point>,
{
    s: DCEL<V, CdtEdge>,
    locate_structure: L,
    all_points_on_line: bool,
    num_constraints: usize,
    __kernel: PhantomData<K>,
}

#[derive(Debug)]
enum ConflictRegion {
    ExistingEdge(FixedEdgeHandle),
    Region {
        left_hull: Vec<FixedEdgeHandle>,
        right_hull: Vec<FixedEdgeHandle>,
        conflicts: Vec<FixedEdgeHandle>,
    }
}

impl <V, K> ConstrainedDelaunayTriangulation<V, K>
    where V: HasPosition2D,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>,
          V::Point: TwoDimensional,
{
    /// Shorthand constructor for a triangulation that is backed up
    /// by an r-tree for log(n) insertion and locate time on average.
    pub fn with_tree_locate() -> ConstrainedDelaunayTriangulation<V, K> {
        ConstrainedDelaunayTriangulation::new()
    }

    /// Shorthand constructor for a triangulation that uses the
    /// `DelaunayWalkLocate` strategy for insertion and point location
    /// queries. This yields O(sqrt(n)) insertion time on average for
    /// randomly generated vertices.
    pub fn with_walk_locate() -> 
        ConstrainedDelaunayTriangulation<V, K, DelaunayWalkLocate> 
    {
        ConstrainedDelaunayTriangulation::new()
    }
}

impl<V, K, L> BasicDelaunaySubdivision<V> for ConstrainedDelaunayTriangulation<V, K, L>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>,
          L: DelaunayLocateStructure<V::Point>

{
    type LocateStructure = L;

    fn locate_structure(&self) -> &Self::LocateStructure {
        &self.locate_structure
    }

    fn locate_structure_mut(&mut self) -> &mut Self::LocateStructure {
        &mut self.locate_structure
    }

    fn all_points_on_line(&self) -> bool {
        self.all_points_on_line
    }

    fn set_all_points_on_line(&mut self, new_value: bool) {
        self.all_points_on_line = new_value;
    }

    fn is_defined_legal(&self, edge: FixedEdgeHandle) -> bool {
        self.s.edge_data(edge).is_constraint_edge()
    }

    fn handle_legal_edge_split(&mut self, handles: &[FixedEdgeHandle; 4]) {
        self.num_constraints += 1;
        for h in handles {
            if !self.is_constraint_edge(*h) {
                self.s.edge_data_mut(*h).make_constraint_edge();
            }
        }
    }
}

impl<V, K, L> HasSubdivision<V> for ConstrainedDelaunayTriangulation<V, K, L>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>,
          L: DelaunayLocateStructure<V::Point>,
{
    type Kernel = K;
    type EdgeType = CdtEdge;

    fn s(&self) -> &DCEL<V, CdtEdge> {
        &self.s
    }

    fn s_mut(&mut self) -> &mut DCEL<V, CdtEdge> {
        &mut self.s
    }
}

impl<V, K, L> ConstrainedDelaunayTriangulation<V, K, L>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>,
          L: DelaunayLocateStructure<V::Point>,
{

    /// Creates a new constrained Delaunay triangulation.
    pub fn new() -> ConstrainedDelaunayTriangulation<V, K, L> {
        ConstrainedDelaunayTriangulation {
            s: DCEL::new_with_edge(),
            all_points_on_line: true,
            locate_structure: Default::default(),
            __kernel: Default::default(),
            num_constraints: 0,
        }
    }

    /// Creates a dynamic vertex handle from a fixed vertex handle.
    ///
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    pub fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V, CdtEdge> {
        self.s.vertex(handle)
    }

    /// Returns a mutable reference to the vertex data referenced by a 
    /// `FixedVertexHandle`.
    pub fn vertex_mut(&mut self, handle: FixedVertexHandle) -> &mut V {
        self.s.vertex_mut(handle)
    }

    /// Creates a dynamic face handle from a fixed face handle.
    ///
    /// May panic if the faces was invalidated by a previous vertex
    /// removal.
    pub fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V, CdtEdge> {
        self.s.face(handle)
    }

    /// Creates a dynamic edge handle from a fixed edge handle.
    ///
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    pub fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V, CdtEdge> {
        self.s.edge(handle)
    }

    /// Returns the number of vertices in this triangulation.
    pub fn num_vertices(&self) -> usize {
        self.s.num_vertices()
    }

    /// Returns the number of faces in this triangulation.
    ///
    /// This count does include the infinite face.
    pub fn num_faces(&self) -> usize {
        self.s.num_faces()
    }

    /// Returns the number of triangles in this triangulation.
    ///
    /// As there is always exactly one face not being a triangle, 
    /// this is equivalent to `self.num_faces() - 1`.
    pub fn num_triangles(&self) -> usize {
        self.s.num_faces() - 1
    }

    /// Returns the number of edges in this triangulation.
    pub fn num_edges(&self) -> usize {
        self.s.num_edges()
    }

    /// Returns an iterator over all triangles.
    pub fn triangles(&self) -> FacesIterator<V, CdtEdge> {
        let mut result = self.s.faces();
        // Skip the outer face
        result.next();
        result
    }

    /// Returns an iterator over all edges.
    pub fn edges(&self) -> EdgesIterator<V, CdtEdge> {
        self.s.edges()
    }

    /// Returns an iterator over all vertices.
    pub fn vertices(&self) -> VerticesIterator<V, CdtEdge> {
        self.s.vertices()
    }

    /// Returns a handle to the infinite face.
    pub fn infinite_face(&self) -> FaceHandle<V, CdtEdge> {
        self.s.face(0)
    }

    /// Returns `true` if the triangulation is degenerate
    /// 
    /// A triangulation is degenerate if all vertices of the
    /// triangulation lie on one line.
    pub fn is_degenerate(&self) -> bool {
        self.all_points_on_line
    }

    /// Returns information about the location of a point in a triangulation.
    pub fn locate(
        &self, point: &V::Point) -> PositionInTriangulation<VertexHandle<V, CdtEdge>, FaceHandle<V, CdtEdge>, EdgeHandle<V, CdtEdge>> {
        self.locate_with_hint_option(point, None)
    }

    /// Locates a vertex at a given position.
    ///
    /// Returns `None` if the point could not be found.
    pub fn locate_vertex(&self, point: &V::Point) -> Option<VertexHandle<V, CdtEdge>> {
        match self.locate(point) {
            PositionInTriangulation::OnPoint(vertex) => Some(vertex),
            _ => None
        }
    }

    /// Returns an edge between two vertices.
    ///
    /// If the edge does not exist, `None` is returned.
    /// This operation runs in `O(n)` time, where `n` is
    /// the degree of `from`.
    pub fn get_edge_from_neighbors(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> Option<EdgeHandle<V, CdtEdge>> {
        self.s.get_edge_from_neighbors(from, to)
    }

    /// Returns information about the location of a point in a triangulation.
    ///
    /// Additionally, a hint can be given to speed up computation. 
    /// The hint should be a vertex close to the position that
    /// is being looked up.
    pub fn locate_with_hint(&self, point: &V::Point, hint: FixedVertexHandle) -> PositionInTriangulation<VertexHandle<V, CdtEdge>, FaceHandle<V, CdtEdge>, EdgeHandle<V, CdtEdge>> {
        self.locate_with_hint_option(point, Some(hint))
    }

    /// Inserts a new vertex into the triangulation.
    ///
    /// A hint can be given to speed up the process.
    /// The hint should be a handle of a vertex close to the new vertex. This
    /// method is recommended in combination with `DelaunayWalkLocate`, 
    /// in this case the insertion time can be reduced to O(1) on average
    /// if the hint is close. If the hint is randomized, running time will
    /// be O(sqrt(n)) on average with an O(n) worst case.
    pub fn insert_with_hint(&mut self, t: V, hint: FixedVertexHandle) -> FixedVertexHandle {
        self.insert_with_hint_option(t, Some(hint))
    }



    /// Attempts to remove a vertex from the triangulation.
    ///
    /// Returns the removed vertex data if it could be found.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles 
    /// upon successful removal.
    pub fn locate_and_remove(&mut self, point: &V::Point) -> Option<V> {
        use self::PositionInTriangulation::*;
        match self.locate_with_hint_option_fixed(point, None) {
            OnPoint(handle) => Some(self.remove(handle)),
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
    pub fn remove(&mut self, vertex: FixedVertexHandle) -> V {
        let num_removed_constraints = self.s.vertex(vertex)
            .ccw_out_edges()
            .map(|edge| self.is_constraint_edge(edge.fix()))
            .filter(|b| *b)
            .count();
        self.num_constraints -= num_removed_constraints;
        BasicDelaunaySubdivision::remove(self, vertex)
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
    pub fn insert(&mut self, vertex: V) -> FixedVertexHandle {
        let handle = self.insert_with_hint_option(vertex, None);
        handle
    }


    /// Returns the number of constraint edges.
    pub fn num_constraints(&self) -> usize {
        self.num_constraints
    }

    /// Returns `true` if a given edge is a constraint edge.
    pub fn is_constraint_edge(&self, edge: FixedEdgeHandle) -> bool {
        self.s.edge_data(edge).is_constraint_edge()
    }

    /// Checks if two vertices are connected by a constraint edge.
    pub fn exists_constraint(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        self.get_edge_from_neighbors(from, to)
            .map(|e| self.is_constraint_edge(e.fix()))
            .unwrap_or(false)
    }

    /// Checks if a constraint edge can be added.
    ///
    /// Returns `false` if the line from `from` to `to` intersects another
    /// constraint edge.
    pub fn can_add_constraint(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        from != to && !self.intersects_any(
            LineIntersectionIterator::new_from_handles(self, from, to))
    }

    /// Checks if a line intersects a constraint edge.
    ///
    /// Returns `true` if the edge from `from` to `to` intersects a 
    /// constraint edge.
    pub fn intersects_constraint(&self, from: &V::Point, to: &V::Point) -> bool {
        self.intersects_any(LineIntersectionIterator::new(self, from, to))
    }

    fn intersects_any(&self, mut iter: LineIntersectionIterator<Self, V, CdtEdge>) -> bool {
        iter.any(|e| {
            if let Intersection::EdgeIntersection(edge) = e {
                self.is_constraint_edge(edge.fix())
            } else {
                false
            }
        })
    }
    
    /// Insert two points and creates a constraint between them.
    ///
    /// Returns `true` if at least one constraint edge was added.
    ///
    /// # Panics
    /// Panics if the new constraint edge intersects with an existing 
    /// constraint edge.
    pub fn add_constraint_edge(&mut self, from: V, to: V) -> bool {
        let from_handle = self.insert(from);
        let to_handle = self.insert(to);
        self.add_constraint(from_handle, to_handle)
    }

    /// Adds a constraint edge between to vertices.
    ///
    /// Returns `true` if at least one constraint edge was added.
    /// Note that the given constraint might be splitted into smaller edges
    /// if a vertex in the triangulation lies exactly on the constraint edge.
    /// Thus, `cdt.exists_constraint(from, to)` is not necessarily `true`
    /// after a call to this function.
    ///
    /// # Panics
    /// Panics if the new constraint edge intersects an existing
    /// constraint edge.
    pub fn add_constraint(&mut self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        assert!(from != to, "Constraint begin must be different from constraint end.");
        // Find edges that cross the constrained edge
        let mut cur_from = from;
        let mut result = false;
        while let Some(region) = self.get_next_conflict_region(cur_from, to) {
            cur_from = match &region {
                &ConflictRegion::ExistingEdge(ref edge) => self.edge(*edge).to().fix(),
                &ConflictRegion::Region { ref right_hull, .. } => 
                    self.edge(*right_hull.last().unwrap()).to().fix(),
            };
            result |= self.resolve_conflict_region(region);
        }
        result
    }

    fn resolve_conflict_region(&mut self, region: ConflictRegion) -> bool {
        match region {
            ConflictRegion::ExistingEdge(edge) => {
                let (edge, sym) = {
                    let edge = self.edge(edge);
                    (edge.fix(), edge.sym().fix())
                };
                if !self.is_constraint_edge(edge) {
                    self.s.edge_data_mut(edge).make_constraint_edge();
                    self.s.edge_data_mut(sym).make_constraint_edge();
                    self.num_constraints += 1;
                    true
                } else {
                    false
                }
            },
            ConflictRegion::Region {
                left_hull,
                right_hull,
                mut conflicts,
            } => {
                assert!(!left_hull.is_empty());
                assert!(!right_hull.is_empty());
                // Assert that no conflict edge is a constraint edge
                assert!(conflicts.iter().all(|e| !self.is_constraint_edge(*e)),
                        "Constraint intersects another constraint edge");
                // Remove edges from highest to lowest value to
                // prevent any index change of a conflicting edge
                conflicts.sort_unstable();
                // Delete edges
                let mut left_vertices = Vec::new();
                for edge in &left_hull {
                    let edge = self.edge(*edge);
                    left_vertices.push((edge.from().fix(), edge.to().fix()));
                }

                let mut right_vertices = Vec::new();
                for edge in &right_hull {
                    let edge = self.edge(*edge);
                    right_vertices.push((edge.from().fix(), edge.to().fix()));
                }

                // Remove conflict edges
                for handle in conflicts.iter().rev().cloned() {
                    self.s.remove_edge(handle, None);
                }

                let prev_edge = right_vertices.last().unwrap();
                let prev_edge = self.s.get_edge_from_neighbors(prev_edge.0, prev_edge.1).unwrap().fix();
                let next_edge = right_vertices.first().unwrap();
                let next_edge = self.s.get_edge_from_neighbors(next_edge.0, next_edge.1).unwrap().fix();

                let constraint_edge = self.s.create_face(prev_edge, next_edge);
                let constraint_sym = self.s.edge(constraint_edge).sym().fix();

                // Create new constraint edge
                self.s.edge_data_mut(constraint_edge).make_constraint_edge();
                self.s.edge_data_mut(constraint_sym).make_constraint_edge();

                // Retriangulate the areas
                let mut edges = Vec::new();
                for &(v0, v1) in &right_vertices {
                    let edge = self.s.get_edge_from_neighbors(v0, v1).unwrap().fix();
                    edges.push(edge);
                }
                edges.push(constraint_edge);
                self.fill_hole(edges);
                let mut edges = Vec::new();
                for &(v0, v1) in left_vertices.iter().rev() {
                    let edge = self.s.get_edge_from_neighbors(v0, v1).unwrap().fix();
                    edges.push(edge);
                }
                edges.push(constraint_sym);
                self.fill_hole(edges);
                self.num_constraints += 1;
                true
            }
        }
    }

    fn get_next_conflict_region(&self,
                         v0: FixedVertexHandle,
                         v1: FixedVertexHandle)
                           -> Option<ConflictRegion> {
        use self::Intersection::*;
        if v0 == v1 {
            return None;
        }
        let mut line_iterator = LineIntersectionIterator::new_from_handles(self, v0, v1);

        let v0 = self.vertex(v0);
        let v1 = self.vertex(v1);

        let mut right_hull = Vec::new();
        let mut left_hull = Vec::new();
        let mut intersecting_edges = Vec::new();

        while let Some(intersection) = line_iterator.next() {
            match intersection {
                Intersection::EdgeIntersection(edge) => {
                    let simple_edge = SimpleEdge::new(
                        edge.from().position(),
                        edge.to().position());
                    assert!(simple_edge.side_query::<K>(&v1.position()).is_on_left_side());
                    intersecting_edges.push(edge.fix());
                },
                Intersection::EdgeOverlap(edge) => {
                    return Some(ConflictRegion::ExistingEdge(edge.fix()));
                },
                VertexIntersection(vertex) => {
                    if vertex != v0 {
                        break;
                    }
                },
            }
        }

        // Create left and right hull from intersection list
        let first_edge = self.edge(intersecting_edges[0]).sym();
        left_hull.push(first_edge.o_next().fix());
        right_hull.push(first_edge.o_prev().fix());
        let mut last_intersection = None;
        for edge in &intersecting_edges {
            let edge = self.edge(*edge);
            if let Some(last_intersection) = last_intersection {
                if edge.sym().o_prev() == last_intersection {
                    left_hull.push(edge.sym().o_next().fix());
                } else {
                    right_hull.push(edge.sym().o_prev().fix());
                }
            }
            last_intersection = Some(edge);
        }
        let last_edge = last_intersection.unwrap();
        left_hull.push(last_edge.o_prev().fix());
        right_hull.push(last_edge.o_next().fix());

        Some(ConflictRegion::Region {
            left_hull: left_hull,
            right_hull: right_hull,
            conflicts: intersecting_edges,
        })
    }

    #[cfg(test)]
    fn cdt_sanity_check(&self) {
        let mut count = 0;
        for edge in self.s.edges() {
            if self.is_constraint_edge(edge.fix()) {
                count += 1;
                assert!(self.is_constraint_edge(edge.sym().fix()));
            }
        }
        assert_eq!(count, self.num_constraints());
        self.sanity_check();
    }
}

impl <V, K> ConstrainedDelaunayTriangulation<V, K, DelaunayTreeLocate<V::Point>>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar> {

    /// Locates the nearest neighbor for a given point.
    pub fn nearest_neighbor(&self, 
                            position: &V::Point)
                            -> Option<VertexHandle<V, CdtEdge>> {
        let entry = self.locate_structure().nearest_neighbor(position);
        entry.map(|e| self.vertex(e.handle))
    }
}

#[cfg(test)]
mod test {
    use testutils::*;
    use super::ConstrainedDelaunayTriangulation;
    use super::{DelaunayTriangulation, DelaunayWalkLocate};
    use traits::HasPosition;
    use kernels::{AdaptiveIntKernel, FloatKernel};
    use cgmath::{Point2, Vector2, EuclideanSpace};
    use rand::{XorShiftRng, SeedableRng, Rng};
    use rand::distributions::{Range, IndependentSample};
    use super::delaunay_basic::BasicDelaunaySubdivision;

    type CDT = ConstrainedDelaunayTriangulation<Point2<f64>, FloatKernel>;
    type Delaunay = DelaunayTriangulation<Point2<f64>, FloatKernel, DelaunayWalkLocate>;


    #[test]
    fn test_add_single_simple_constraint() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(2.0, 2.0));
        let v2 = cdt.insert(Point2::new(1.0, 0.5));
        let v3 = cdt.insert(Point2::new(0.5, 1.0));
        assert!(cdt.get_edge_from_neighbors(v0, v1).is_none());
        assert!(cdt.get_edge_from_neighbors(v2, v3).is_some());

        assert!(cdt.add_constraint(v1, v0));
        assert!(!cdt.add_constraint(v0, v1));
        let edge = cdt.get_edge_from_neighbors(v0, v1).expect("Expected constraint edge").fix();
        assert!(cdt.get_edge_from_neighbors(v2, v3).is_none());
        assert!(cdt.is_constraint_edge(edge));
        cdt.cdt_sanity_check();
    }

    #[test]
    fn test_existing_edge_constraint() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(2.0, 2.0));
        let v2 = cdt.insert(Point2::new(1.0, 0.0));
        assert!(cdt.add_constraint(v0, v1));
        assert!(cdt.add_constraint(v0, v2));
        assert!(cdt.add_constraint(v1, v2));
        for edge in cdt.edges() {
            assert!(cdt.is_constraint_edge(edge.fix()));
        }
        assert!(!cdt.add_constraint(v1, v0));
        assert!(!cdt.add_constraint(v1, v2));
        assert_eq!(cdt.num_constraints, 3);
    }

    #[test]
    fn test_mid_overlapping_constraint() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.5));
        let v1 = cdt.insert(Point2::new(2.0, 0.5));
        let v2 = cdt.insert(Point2::new(3.0, 0.5));
        let v3 = cdt.insert(Point2::new(5.0, 0.5));
        cdt.insert(Point2::new(1.0, 1.0));
        cdt.insert(Point2::new(1.0, 0.0));
        cdt.insert(Point2::new(3.0, 1.0));
        cdt.insert(Point2::new(3.0, 0.0));
        assert!(cdt.get_edge_from_neighbors(v1, v2).is_some());
        let mut copy = cdt.clone();
        assert!(cdt.add_constraint(v0, v3));
        assert_eq!(cdt.num_constraints(), 3);
        copy.add_constraint(v2, v3);
        assert_eq!(copy.num_constraints(), 1);
        copy.add_constraint(v0, v3);
        assert_eq!(copy.num_constraints(), 3);
    }

    #[test]
    fn test_add_single_complex_constraint() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        cdt.insert(Point2::new(1.0, 0.0));
        cdt.insert(Point2::new(0.0, 1.0));
        cdt.insert(Point2::new(2.0, 1.0));
        // cdt.insert(Point2::new(1.0, 2.0));
        let v1 = cdt.insert(Point2::new(2.0, 2.0));
        assert!(cdt.get_edge_from_neighbors(v0, v1).is_none());
        cdt.add_constraint(v0, v1);
        let edge = cdt.get_edge_from_neighbors(v0, v1).expect("Expected constraint edge").fix();
        assert!(cdt.is_constraint_edge(edge));
        
    }

    #[test]
    fn test_add_single_constraint() {
        let points = random_points_with_seed::<f64>(1000, [3123, 23311, 2212, 3939291]);
        let mut cdt = CDT::new();
        assert_eq!(cdt.num_constraints(), 0);
        for point in points {
            cdt.insert(point);
        }
        cdt.add_constraint(40, 200);
        assert_eq!(cdt.num_constraints(), 1);
        cdt.cdt_sanity_check();
    }

    #[test]
    fn test_add_border_constraint() {
        let points = random_points_with_seed::<f64>(1000, [3201293, 2863311, 2899762, 356429991]);
        let mut cdt = CDT::new();
        let mut max_y = -::std::f64::MAX;
        for point in points {
            max_y = max_y.max(point.y);
            cdt.insert(point);
        }
        let v0 = cdt.insert(Point2::new(-20., max_y + 10.));
        let v1 = cdt.insert(Point2::new(20., max_y + 10.));
        cdt.add_constraint(v0, v1);
        assert_eq!(cdt.num_constraints(), 1);
        cdt.cdt_sanity_check();
    }

    #[test]
    fn test_add_multiple_constraints_overlapping() {
        test_add_multiple_constraints(true);
    }

    #[test]
    fn test_add_multiple_constraints_non_overlapping() {
        test_add_multiple_constraints(false);
    }

    fn test_add_multiple_constraints(overlapping: bool) {
        const RANGE: f64 = 10.;
        let seed = if overlapping { 
            [3201293, 2863311, 2899762, 356429991] 
        } else {
            [266260112, 1351, 2095, 933321]
        };
        let points =
            random_points_in_range::<f64>(RANGE, 1000, seed);
        let mut cdt = CDT::new();
        for point in points {
            cdt.insert(point);
        }
        let seed = if overlapping {
            [212231, 332919, 380283, 211399]
        } else {
            [91231, 3391089, 3883, 21199]
        };
        let delaunay_points =
            random_points_in_range::<f64>(RANGE * 0.9, 80, seed);
        // Use a delaunay triangulation to "generate" non intersecting constraint edges
        let mut d = Delaunay::new();
        for p in delaunay_points {
            d.insert(p);
        }
        let mut used_vertices = ::std::collections::HashSet::new();
        let mut inserted_constraints = Vec::new();
        for v in d.vertices() {
            // Insert only edges that do not touch at the end points if 
            // overlapping is false
            if overlapping || used_vertices.insert(v.fix()) {
                let out_edge = v.out_edge().unwrap();
                let to = out_edge.to();
                used_vertices.insert(to.fix());
                let h0 = cdt.insert(v.position());
                let h1 = cdt.insert(to.position());
                if cdt.add_constraint(h0, h1) {
                    inserted_constraints.push((h0, h1));
                }
                assert_eq!(cdt.num_constraints(), inserted_constraints.len());
            }
        }
        // Check if all constraints still exists
        for (from, to) in inserted_constraints {
            assert!(cdt.exists_constraint(from, to));
        }
        cdt.cdt_sanity_check();
    }

    #[test]
    fn test_split_constraint() {
        let mut cdt = CDT::new();
        cdt.insert(Point2::new(0.0, 0.0));
        cdt.insert(Point2::new(1.0, 0.0));
        cdt.insert(Point2::new(0.0, 1.0));
        let v0 = cdt.insert(Point2::new(0.0, 0.5));
        let v_last = cdt.insert(Point2::new(1.0, 0.5));
        cdt.add_constraint(v0, v_last);
        assert_eq!(cdt.num_constraints(), 1);
        // These points split an existing constraint
        let v1 = cdt.insert(Point2::new(0.25, 0.5));
        assert_eq!(cdt.num_constraints(), 2);
        let v2 = cdt.insert(Point2::new(0.75, 0.5));
        assert_eq!(cdt.num_constraints(), 3);
        assert!(cdt.exists_constraint(v0, v1));
        assert!(cdt.exists_constraint(v1, v2));
        assert!(cdt.exists_constraint(v2, v_last));
        cdt.cdt_sanity_check();
    }

    #[test]
    fn test_add_constraint_over_point() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(1.0, 0.0));
        let v2 = cdt.insert(Point2::new(2.0, 0.0));
        cdt.insert(Point2::new(0.0, 1.0));
        cdt.add_constraint(v0, v2);
        assert_eq!(cdt.num_constraints(), 2);
        assert!(cdt.exists_constraint(v0, v1));
        assert!(cdt.exists_constraint(v1, v2));
        cdt.cdt_sanity_check();
    }

    fn test_cdt() -> CDT {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(1.0, 0.0));
        let v1 = cdt.insert(Point2::new(0.0, 1.0));
        cdt.insert(Point2::new(0.0, 0.0));
        cdt.insert(Point2::new(1.0, 1.0));
        cdt.add_constraint(v0, v1);
        cdt
    }

    #[test]
    fn test_check_intersects_constraint_edge() {
        let cdt = test_cdt();
        let from = Point2::new(0.2, 0.2);
        let to = Point2::new(0.6, 0.7);
        assert!(cdt.intersects_constraint(&from, &to));
        assert!(cdt.intersects_constraint(&to, &from));
        let to = Point2::new(-0.5, 0.2);
        assert!(!cdt.intersects_constraint(&from, &to));
        let from = Point2::new(0.5, 0.5);
        assert!(cdt.intersects_constraint(&from, &to));
        assert!(cdt.intersects_constraint(&to, &from));
    }

    #[test]
    fn test_add_constraint_degenerate() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(0.0, 1.0));
        assert!(cdt.add_constraint(v0, v1));
        assert!(!cdt.add_constraint(v1, v0));
        assert_eq!(cdt.num_constraints(), 1);
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(0.0, 2.0));
        cdt.insert(Point2::new(0.0, 1.0));
        assert!(cdt.add_constraint(v0, v1));
        assert_eq!(cdt.num_constraints(), 2);
    }

    fn random_points_on_line<R>(range: i64, 
                                num_points: usize,
                                rng: &mut R, 
                                line_dir: Vector2<i64>) -> Vec<Point2<i64>> 
        where R: Rng
    {
        let mut result = Vec::with_capacity(num_points);
        let range = Range::new(-range, range);
        for _ in 0 .. num_points {
            let factor = range.ind_sample(rng);
            result.push(Point2::from_vec(line_dir * factor));
        }
        result
    }

    #[test]
    fn fuzz_test_on_line() {
        // Generates points on a single line and randomly connects
        // them with constraints.
        let seed = [441217, 22531, 4999902, 292791];
        const RANGE: i64 = 10000;
        const NUM_POINTS: usize = 2000;
        let mut rng = XorShiftRng::from_seed(seed);
        let points = random_points_on_line(RANGE, NUM_POINTS, &mut rng, Vector2::new(1, 1));
        let mut cdt = ConstrainedDelaunayTriangulation::<_, AdaptiveIntKernel>::with_walk_locate();
        for ps in points.chunks(2) {
            let from = ps[0];
            let to = ps[1];
            let from = cdt.insert(from);
            let to = cdt.insert(to);
            if from != to && rng.gen() {
                cdt.add_constraint(from, to);
            }
        }
        cdt.sanity_check();
        assert!(cdt.is_degenerate());
    }

    #[test]
    fn fuzz_test_on_grid() {
        // Generates points on a grid and randomly connects
        // them with non intersecting constraints
        let seed = [5587, 568731, 955432, 215512];
        let mut points = Vec::with_capacity((RANGE * RANGE) as usize);
        const RANGE: i64 = 30;
        const NUM_CONSTRAINTS: usize = 2000;
        for x in -RANGE .. RANGE {
            for y in -RANGE .. RANGE {
                points.push(Point2::new(x, y));
            }
        }
        let mut rng = XorShiftRng::from_seed(seed);
        rng.shuffle(&mut points);
        let mut cdt = ConstrainedDelaunayTriangulation::<_, AdaptiveIntKernel>::with_walk_locate();
        for p in points {
            cdt.insert(p);
        }
        let range = Range::new(-RANGE, RANGE);
        let directions_and_offset = [(Vector2::new(1, 0), Point2::new(0, 1)),
                                     (Vector2::new(0, 1), Point2::new(1, 0)),
                                     (Vector2::new(1, 1), Point2::new(0, 0))];
        for _ in 0 .. NUM_CONSTRAINTS {
            let &(direction, offset) = rng.choose(&directions_and_offset).unwrap();
            let factor1 = range.ind_sample(&mut rng);
            let factor2 = range.ind_sample(&mut rng);
            let p1 = offset + direction * factor1;
            let p2 = offset + direction * factor2;
            if p1 != p2 {
                cdt.add_constraint_edge(p1, p2);
            }
        }
        cdt.sanity_check();
    }

    #[test]
    fn test_cdt_remove_degenerate() {
        let mut cdt = CDT::new();
        let v0 = cdt.insert(Point2::new(0.0, 0.0));
        let v1 = cdt.insert(Point2::new(1.0, 0.0));
        let v2 = cdt.insert(Point2::new(0.0, 1.0));
        cdt.add_constraint(v0, v1);
        cdt.add_constraint(v1, v2);
        cdt.add_constraint(v2, v0);
        assert_eq!(cdt.num_constraints(), 3);
        assert!(!cdt.is_degenerate());
        cdt.remove(v1);
        assert_eq!(cdt.num_constraints(), 1);
        assert!(cdt.is_degenerate());
    }
}
