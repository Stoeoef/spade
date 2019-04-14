// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::kernels::{DelaunayKernel, FloatKernel, TrivialKernel};
use crate::point_traits::{PointN, PointNExtensions, ThreeDimensional, TwoDimensional};
use crate::primitives::{SimpleEdge, SimpleTriangle};
use crate::traits::{HasPosition, HasPosition2D, SpadeFloat, SpatialObject};
use num::{one, zero, Float, One, Zero};
use smallvec::{smallvec, SmallVec};
use std::marker::PhantomData;

use self::dcel::*;
use self::delaunay_basic::{BasicDelaunaySubdivision, HasSubdivision};
use self::delaunay_locate::*;
use crate::delaunay::*;

/// Type shorthand for a Delaunay triangulation with `f64` coordinates that uses `FloatKernel`
/// for geometric calculations.
pub type FloatDelaunayTriangulation<T, L> = DelaunayTriangulation<T, FloatKernel, L>;
/// Type shorthand for a Delaunay triangulation with `i64` or `i32` coordinates that uses
/// the trivial kernel for geometric calculations.
pub type IntDelaunayTriangulation<T, L> = DelaunayTriangulation<T, TrivialKernel, L>;

/// Stores information about a point's position in triangulation.
///
/// Used as a return type of `DelaunayTriangulation::locate(..)`.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum PositionInTriangulation<V: Copy, F: Copy, E: Copy> {
    /// The point is contained in a triangle.
    InTriangle(F),
    /// The point is outside the convex hull. The given edge is an edge that
    /// is close to the queried position.
    OutsideConvexHull(E),
    /// A vertex with this position has already been inserted. Its handle is given.
    OnPoint(V),
    /// The point lies on an edge.
    OnEdge(E),
    /// There is no valid triangulation yet, thus, less than two points where
    /// inserted.
    NoTriangulationPresent,
}

/// A two dimensional Delaunay triangulation.
///
/// A Delaunay triangulation is a special triangulation of a set of points that fulfills some
/// suitable properties for geometric operations like interpolation.
/// There is also an [own chapter](https://stoeoef.gitbooks.io/spade-user-manual/content/delaunay-triangulation.html) in the user guide covering spade's triangulation.
///
/// Objects that are inserted into the triangulation have to implement the `HasPosition2D` trait.
/// The trait is implemented for all types that implement `TwoDimensional`, like `Point2` from
/// the `cgmath` and `nalgebra` package or `[S; 2]` for `S: SpadeNum`.
///
/// A straightforward Delaunay triangulation implementation will suffer from precision problems:
/// various geometric queries can fail if imprecise calculations (like native `f32` / `f64` operations)
/// are used. Those failures can yield to incorrect results or panics at runtime.
/// To prevent those crashes, Spade offers a few "calculation kernels" that may fit the
/// individual needs of an application. See `spade::kernels` for more information.
///
/// # Example
///
/// ```
/// extern crate nalgebra;
/// extern crate spade;
///
/// use nalgebra::Point2;
/// use spade::delaunay::FloatDelaunayTriangulation;
///
/// # fn main() {
///   let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
///   delaunay.insert(Point2::new(0.0, 1.0));
///   delaunay.insert(Point2::new(0.0, -1.0));
///   delaunay.insert(Point2::new(1.0, 0.0));
///   delaunay.insert(Point2::new(-1.0, 0.0));
///   delaunay.insert(Point2::new(0.0, 0.0));
///   for face in delaunay.triangles() {
///     let triangle = face.as_triangle();
///     println!("Found triangle: {:?} -> {:?} -> {:?}", *triangle[0], *triangle[1], *triangle[2]);
///   }
///   for edge in delaunay.edges() {
///     println!("Found an edge: {:?} -> {:?}", *edge.from(), *edge.to());
///   }
/// # }
/// ```
///
/// # Iterating
/// A triangulation consists of three elements - vertices, edges and triangles - that can be iterated over.
/// Use `vertices()` `edges()` and `triangles()` to call appropriate, non-mutating iterators.
///
/// # Infinite face and convex hull
/// Every triangulation is surrounded by the _infinite face_. This face can be retrieved by calling
/// `infinite_face()`. Iterating the adjacent edges of the infinite face will yield the triangulation's
/// convex hull. See `FaceHandle::adjacent_edges()`.
///
/// # Mutating
/// Vertices can be added and removed. Mutation of vertex data is possible with `lookup_mut(..)`,
/// although this mutation must not alter the vertex position.
///
/// # Interpolation
/// Vertices can store various user data, the points could for example represent meteorological samples
/// measuring temperature or air moisture. These values can be interpolated smoothly within the
/// triangulation's convex hull using a variety of interpolation methods. There are currently four supported
/// interpolation methods:
///  - `barycentric_interpolation(..)`
///  - `nn_interpolation(..)`
///  - `nn_interpolation_c1_sibson(..)`
///  - `nn_interpolation_c1_farin(..)`
///
/// # Type parameters
/// `DelaunayTriangulation` has three type parameters: `V`, `K` and `L`.
/// `V: HasPosition2D` defines the triangulation's vertex type.
/// `K: DelaunayKernel` defines the triangulations calculation kernel.
/// For more information, see `spade::kernels`.
/// `L` Defines the locate structure.
/// For more information, see `DelaunayLocateStructure`.
///
/// # Performance
/// Performance of insertion, interpolation and other queries heavily relies on
/// - the locate structure being used
/// - the closeness of subsequent queries
/// - the kernel being used
/// Depending on the use case, it can vary between O(1) to O(sqrt(n)) on average.
/// The guide has an [own chapter](https://stoeoef.gitbooks.io/spade-user-manual/content/triangulation-performance.html)
/// about performance.
///
/// ## Auto Hinting
/// Since version 1.1, spade uses the result of the last query as hint for the next query when
/// using `DelaunayWalkLocate` as locate strategy. As a consequence, subsequent
/// queries - like insertion, interpolation or nearest neighbor queries - will
/// run in O(1) if the query locations are close to each other.
#[derive(Debug)]
#[cfg_attr(feature = "serde_serialize", derive(Serialize, Deserialize))]
pub struct DelaunayTriangulation<V, K, L = DelaunayTreeLocate<<V as HasPosition>::Point>>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
{
    __kernel: PhantomData<*const K>,
    s: DCEL<V>,
    all_points_on_line: bool,
    locate_structure: L,
}

impl<V, K, L> BasicDelaunaySubdivision<V> for DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
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
}

impl<V, K, L> HasSubdivision<V> for DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
{
    type Kernel = K;
    type EdgeType = ();

    fn s(&self) -> &DCEL<V> {
        &self.s
    }

    fn s_mut(&mut self) -> &mut DCEL<V> {
        &mut self.s
    }
}

impl<V, K, L> Clone for DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D + Clone,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
{
    fn clone(&self) -> DelaunayTriangulation<V, K, L> {
        DelaunayTriangulation {
            __kernel: Default::default(),
            s: self.s.clone(),
            all_points_on_line: self.all_points_on_line,
            locate_structure: self.locate_structure.clone(),
        }
    }
}

impl<V, K, L> Default for DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
{
    fn default() -> Self {
        DelaunayTriangulation::new()
    }
}

impl<V, K> DelaunayTriangulation<V, K>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
{
    /// Shorthand constructor for a Delaunay triangulation that is backed up by an r-tree for
    /// log(n) insertion and locate time on average.
    pub fn with_tree_locate() -> DelaunayTriangulation<V, K> {
        DelaunayTriangulation::new()
    }

    /// Shorthand constructor for a Delaunay triangulation that uses the
    /// `DelaunayWalkLocate` strategy for insertion and point location
    /// queries. This yields O(sqrt(n)) insertion time on average for
    /// randomly generated vertices.
    pub fn with_walk_locate() -> DelaunayTriangulation<V, K, DelaunayWalkLocate> {
        DelaunayTriangulation::new()
    }
}

impl<V, K, L> DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
    L: DelaunayLocateStructure<V::Point>,
{
    /// Creates a new Delaunay triangulation.
    ///
    /// Using this method directly can be a bit cumbersome due to type annotations, consider using
    /// The short hand definitions `FloatDelaunayTriangulation` or `IntDelaunayTriangulation` in
    /// combination with the methods `with_walk_locate` and `with_tree_locate`.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    /// use spade::delaunay::{DelaunayTriangulation, DelaunayWalkLocate};
    /// use spade::kernels::FloatKernel;
    /// # fn main() {
    /// let mut triangulation = DelaunayTriangulation::<_, FloatKernel, DelaunayWalkLocate>::new();
    /// # triangulation.insert(nalgebra::Point2::new(332f64, 123f64));
    /// # }
    /// ```
    ///
    /// Usually, the omitted types (the triangulation's vertex type) can be inferred from a call
    /// to `insert`.
    pub fn new() -> DelaunayTriangulation<V, K, L> {
        DelaunayTriangulation {
            __kernel: Default::default(),
            s: DCEL::new(),
            all_points_on_line: true,
            locate_structure: Default::default(),
        }
    }

    /// Creates a dynamic vertex handle from a fixed vertex handle.
    ///
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    pub fn vertex(&self, handle: FixedVertexHandle) -> VertexHandle<V> {
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
    pub fn face(&self, handle: FixedFaceHandle) -> FaceHandle<V> {
        self.s.face(handle)
    }

    /// Creates a dynamic edge handle from a fixed edge handle.
    ///
    /// May panic if the handle was invalidated by a previous vertex
    /// removal.
    pub fn edge(&self, handle: FixedEdgeHandle) -> EdgeHandle<V> {
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
    pub fn triangles(&self) -> FacesIterator<V> {
        let mut result = self.s.faces();
        // Skip the outer face
        result.next();
        result
    }

    /// Returns an iterator over all undirected edges.
    pub fn edges(&self) -> EdgesIterator<V> {
        self.s.edges()
    }

    /// Returns an iterator over all vertices.
    pub fn vertices(&self) -> VerticesIterator<V> {
        self.s.vertices()
    }

    /// Returns a handle to the infinite face.
    pub fn infinite_face(&self) -> FaceHandle<V> {
        self.s.face(0)
    }

    /// Returns `true` if the triangulation is degenerate
    ///
    /// A triangulation is degenerate if all vertices of the
    /// triangulation lie on one line.
    pub fn is_degenerate(&self) -> bool {
        self.all_points_on_line
    }

    /// Returns an edge between two vertices.
    ///
    /// If the edge does not exist, `None` is returned.
    /// This operation runs in `O(n)` time, where `n` is
    /// the degree of `from`.
    pub fn get_edge_from_neighbors(
        &self,
        from: FixedVertexHandle,
        to: FixedVertexHandle,
    ) -> Option<EdgeHandle<V>> {
        self.s.get_edge_from_neighbors(from, to)
    }

    /// Locates the nearest neighbor for a given point.
    ///
    /// Returns `None` if the triangulation is empty.
    pub fn nearest_neighbor(&self, point: &V::Point) -> Option<VertexHandle<V>> {
        if self.num_vertices() == 0 {
            return None;
        }
        let start = self.get_default_hint(point);
        let mut cur = self.vertex(start);
        let mut min_dist = cur.position().distance2(point);
        'outer: loop {
            for edge in cur.ccw_out_edges() {
                let neighbor = edge.to();
                let n_distance = neighbor.position().distance2(point);
                if n_distance < min_dist {
                    cur = neighbor;
                    min_dist = n_distance;
                    continue 'outer;
                }
            }
            break;
        }
        self.locate_structure.new_query_result(cur.fix());
        Some(cur)
    }

    /// Returns information about the location of a point in a triangulation.
    pub fn locate(
        &self,
        point: &V::Point,
    ) -> PositionInTriangulation<VertexHandle<V>, FaceHandle<V>, EdgeHandle<V>> {
        self.locate_with_hint_option(point, None)
    }

    /// Locates a vertex at a given position.
    ///
    /// Returns `None` if the point could not be found.
    pub fn locate_vertex(&self, point: &V::Point) -> Option<VertexHandle<V>> {
        if let Some(nn) = self.nearest_neighbor(point) {
            if &nn.position() == point {
                return Some(nn);
            }
        }
        None
    }

    /// Returns information about the location of a point in a triangulation.
    ///
    /// Additionally, a hint can be given to speed up computation. The hint should be a vertex close
    /// to the position that is being looked up.
    pub fn locate_with_hint(
        &self,
        point: &V::Point,
        hint: FixedVertexHandle,
    ) -> PositionInTriangulation<VertexHandle<V>, FaceHandle<V>, EdgeHandle<V>> {
        self.locate_with_hint_option(point, Some(hint))
    }

    /// Inserts a new vertex into the triangulation.
    ///
    /// This operation runs in `O(log(n))` on average when using a tree lookup to back up the
    /// triangulation, or in `O(sqrt(n))` when using a walk lookup. `n` denotes the number of vertices,
    /// the given running times assume that input data is given uniformly randomly distributed.
    /// If the point has already been contained in the triangulation, the old vertex is overwritten.
    ///
    /// Returns a handle to the new vertex. Use this handle with
    /// `DelaunayTriangulation::vertex(..)` to refer to it.
    pub fn insert(&mut self, t: V) -> FixedVertexHandle {
        self.insert_with_hint_option(t, None)
    }

    /// Inserts a new vertex into the triangulation.
    ///
    /// A hint can be given to speed up the process. The hint should be a handle of a vertex
    /// close to the new vertex. This method is recommended in combination with
    /// `DelaunayWalkLocate`, in this case the insertion time can be reduced
    /// to O(1) on average if the hint is close. If the hint is randomized, running time will be O(sqrt(n))
    /// on average with an O(n) worst case.
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
        BasicDelaunaySubdivision::remove(self, vertex)
    }
}

impl<V, K> DelaunayTriangulation<V, K, DelaunayTreeLocate<V::Point>>
where
    V: HasPosition2D,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    V::Point: TwoDimensional,
{
    /// Checks if the triangulation contains an object with a given coordinate.
    #[deprecated(since = "1.3.0", note = "Use locate_vertex instead")]
    pub fn lookup(&self, point: &V::Point) -> Option<VertexHandle<V>> {
        let handle = self.locate_structure.lookup(point);
        handle.map(|h| self.s.vertex(h.handle))
    }

    /// Removes a vertex at a given position from the triangulation.
    ///
    /// This operation runs in O(n²), where n is the degree of the
    /// removed vertex.
    ///
    /// # Handle invalidation
    /// This method will invalidate all vertex, edge and face handles upon
    /// successful removal.
    #[deprecated(since = "1.3.0", note = "Use locate_and_remove instead")]
    pub fn lookup_and_remove(&mut self, point: &V::Point) -> Option<V> {
        let handle = self.locate_vertex(point).map(|h| h.fix());
        handle.map(|h| self.remove(h))
    }
}

const INTPL_SMALLVEC_CAPACITY: usize = 8;

impl<V, K, L> DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    <V::Point as PointN>::Scalar: SpadeFloat,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    L: DelaunayLocateStructure<V::Point>,
    V::Point: TwoDimensional,
{
    /// Performs a barycentric interpolation.
    /// Returns `None` if the triangulation has no triangles yet.
    /// Points outside of the convex hull will be interpolated as well.
    /// The other interpolation methods are used very similarly, check their
    /// documentation for an example.
    pub fn barycentric_interpolation<F>(
        &self,
        point: &V::Point,
        f: F,
    ) -> Option<<V::Point as PointN>::Scalar>
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
    {
        let vertices: SmallVec<[_; 3]> = match self.locate(point) {
            PositionInTriangulation::NoTriangulationPresent => return None,
            PositionInTriangulation::OnPoint(v) => smallvec![v],
            PositionInTriangulation::OnEdge(e) => smallvec![e.from(), e.to()],
            PositionInTriangulation::InTriangle(f) => {
                let vs = f.as_triangle();
                smallvec![vs[0], vs[1], vs[2]]
            }
            PositionInTriangulation::OutsideConvexHull(e) => smallvec![e.from(), e.to()],
        };
        if vertices.len() == 1 {
            Some(f(&*vertices[0]))
        } else if vertices.len() == 2 {
            let p0 = vertices[0].position();
            let p1 = vertices[1].position();
            let one = <<V::Point as PointN>::Scalar>::one();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(zero(), edge.project_point(point), one);
            let w0 = one - w1;
            Some(w1 * f(&*vertices[1]) + w0 * f(&*vertices[0]))
        } else {
            let triangle = crate::primitives::SimpleTriangle::new(
                vertices[0].position(),
                vertices[1].position(),
                vertices[2].position(),
            );
            let b_coords = triangle.barycentric_interpolation(point);
            let w0 = f(&*vertices[0]);
            let w1 = f(&*vertices[1]);
            let w2 = f(&*vertices[2]);
            Some(w0 * b_coords[0] + w1 * b_coords[1] + w2 * b_coords[2])
        }
    }

    /// Performs a natural neighbor interpolation for a given position.
    ///
    /// Returns `None` if the triangulation has no triangles yet.
    /// Points outside of the convex hull will be interpolated as well.
    ///
    /// # Example
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    ///
    /// # use nalgebra::{Point2};
    /// use spade::HasPosition;
    /// use spade::delaunay::{FloatDelaunayTriangulation};
    ///
    /// struct PointWithHeight {
    ///   point: Point2<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Point = Point2<f64>;
    ///     fn position(&self) -> Point2<f64> {
    ///       self.point
    ///     }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
    ///   delaunay.insert(PointWithHeight { point: Point2::new(0.0, 0.0), height: 5. });
    ///   delaunay.insert(PointWithHeight { point: Point2::new(1.0, 0.0), height: 0. });
    ///   delaunay.insert(PointWithHeight { point: Point2::new(0.0, 1.0), height: 0. });
    ///   let lookup = Point2::new(0.2, 0.2);
    ///   // Interpolate the points height
    ///   let interpolated = delaunay.nn_interpolation(&lookup, |p| p.height).unwrap();
    ///   // and insert it afterwards.
    ///   delaunay.insert(PointWithHeight { point: lookup, height: interpolated });
    ///   // Data points themselves will always yield their own height
    ///   assert_eq!(delaunay.nn_interpolation(&Point2::new(0.0, 0.0), |p| p.height),
    ///              Some(5.0));
    /// }
    /// ```
    pub fn nn_interpolation<F>(
        &self,
        point: &V::Point,
        f: F,
    ) -> Option<<V::Point as PointN>::Scalar>
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
    {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);

        let mut sum = None;
        for (index, fixed_handle) in nns.iter().enumerate() {
            let new = f(&*self.s.vertex(*fixed_handle)) * ws[index];
            sum = match sum {
                Some(val) => Some(val + new),
                None => Some(new),
            }
        }
        sum
    }

    fn get_weights(
        &self,
        nns: &SmallVec<[FixedVertexHandle; INTPL_SMALLVEC_CAPACITY]>,
        point: &V::Point,
    ) -> SmallVec<[<V::Point as PointN>::Scalar; INTPL_SMALLVEC_CAPACITY]> {
        let mut result = SmallVec::new();
        if nns.len() == 1 {
            result.push(one());
            return result;
        }

        if nns.len() == 2 {
            let p0 = (*self.s.vertex(nns[0])).position();
            let p1 = (*self.s.vertex(nns[1])).position();
            let one = <<V::Point as PointN>::Scalar>::one();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(zero(), edge.project_point(point), one);
            let w0 = one - w1;
            result.push(w0);
            result.push(w1);
            return result;
        }
        let len = nns.len();

        // Get voronoi vertices of adjacent faces
        let mut point_cell: SmallVec<[_; 16]> = SmallVec::new();
        for (index, cur) in nns.iter().enumerate() {
            let cur_pos = (*self.s.vertex(*cur)).position();
            let next = (*self.s.vertex(nns[(index + 1) % len])).position();
            let triangle = SimpleTriangle::new(next, cur_pos, point.clone());
            point_cell.push(triangle.circumcenter());
        }

        let mut total_area = zero();
        let mut ccw_edge = self
            .s
            .get_edge_from_neighbors(*nns.first().unwrap(), *nns.last().unwrap())
            .unwrap();
        for (index, cur) in nns.iter().enumerate() {
            // Calculate area of voronois cells
            let cur_pos = (*self.s.vertex(*cur)).position();
            let next = nns[(index + 1) % len];

            let mut polygon_area = zero();
            let mut last = point_cell[((index + len) - 1) % len].clone();
            let first = point_cell[index].clone();
            loop {
                if ccw_edge.to().fix() == next {
                    break;
                }
                let cw_edge = ccw_edge.cw();
                let triangle = SimpleTriangle::new(
                    (*ccw_edge.to()).position().clone(),
                    (*cw_edge.to()).position().clone(),
                    cur_pos.clone(),
                );
                let cur = triangle.circumcenter();

                let tri = SimpleTriangle::new(first.clone(), cur.clone(), last);
                last = cur;
                polygon_area += tri.double_area();
                ccw_edge = cw_edge;
            }
            ccw_edge = ccw_edge.sym();

            total_area += polygon_area;
            result.push(polygon_area);
        }
        for area in &mut result {
            *area /= total_area;
        }
        result
    }

    fn get_natural_neighbors(
        &self,
        position: &V::Point,
    ) -> SmallVec<[FixedVertexHandle; INTPL_SMALLVEC_CAPACITY]> {
        match self.locate_with_hint_option_fixed(position, None) {
            PositionInTriangulation::InTriangle(face) => {
                let mut edges: SmallVec<_> = self
                    .face(face)
                    .adjacent_edges()
                    .rev()
                    .map(|e| e.sym().fix())
                    .collect();
                self.inspect_flips(&mut edges, position)
            }
            PositionInTriangulation::OnEdge(edge) => {
                let edge = self.edge(edge);
                if self.is_ch_edge(edge.fix()) || self.is_ch_edge(edge.sym().fix()) {
                    let mut vec = SmallVec::new();
                    vec.push(edge.from().fix());
                    vec.push(edge.to().fix());
                    vec
                } else {
                    let mut edges = SmallVec::new();
                    edges.push(edge.o_prev().sym().fix());
                    edges.push(edge.o_next().sym().fix());
                    edges.push(edge.sym().o_prev().sym().fix());
                    edges.push(edge.sym().o_next().sym().fix());
                    self.inspect_flips(&mut edges, position)
                }
            }
            PositionInTriangulation::OutsideConvexHull(edge) => {
                // Get the closest edge on the convex hull
                let edges = self.get_convex_hull_edges_for_point(edge, position);
                let mut min_dist = <V::Point as PointN>::Scalar::max_value();
                let mut min_edge = 0;
                for cur_edge in edges {
                    let edge = self.s.edge(cur_edge);
                    let simple = Self::to_simple_edge(edge);
                    let new_dist = simple.distance2(position);
                    if new_dist < min_dist {
                        min_edge = cur_edge;
                        min_dist = new_dist;
                    }
                }
                let min_edge = self.s.edge(min_edge);
                let from_pos = (*min_edge.from()).position();
                let to_pos = (*min_edge.to()).position();
                let simple = Self::to_simple_edge(min_edge);
                let mut result = SmallVec::new();
                if simple.is_projection_on_edge(position) {
                    result.push(min_edge.from().fix());
                    result.push(min_edge.to().fix());
                } else {
                    // Return the closer point
                    if from_pos.distance2(position) < to_pos.distance2(position) {
                        result.push(min_edge.from().fix());
                    } else {
                        result.push(min_edge.to().fix());
                    }
                }
                result
            }
            PositionInTriangulation::OnPoint(fixed_handle) => {
                let mut result = SmallVec::new();
                result.push(fixed_handle);
                result
            }
            _ => SmallVec::new(),
        }
    }

    fn inspect_flips(
        &self,
        edges: &mut SmallVec<[FixedEdgeHandle; INTPL_SMALLVEC_CAPACITY]>,
        position: &V::Point,
    ) -> SmallVec<[FixedVertexHandle; INTPL_SMALLVEC_CAPACITY]> {
        let mut result = SmallVec::new();

        while let Some(e) = edges.pop() {
            if self.is_ch_edge(e) {
                result.push(self.edge(e).from().fix());
            } else {
                let (v0, v1, v2, e1, e2);
                {
                    let edge = self.s.edge(e);
                    v0 = (*edge.from()).position();
                    v1 = (*edge.to()).position();
                    v2 = (*edge.ccw().to()).position();
                    e1 = edge.o_prev().sym().fix();
                    e2 = edge.o_next().sym().fix();
                }
                debug_assert!(K::is_ordered_ccw(&v1, &v2, &v0));
                debug_assert!(K::is_ordered_ccw(position, &v1, &v0));
                if K::contained_in_circumference(&v2, &v1, &v0, &position) {
                    // The edge is illegal
                    edges.push(e1);
                    edges.push(e2);
                } else {
                    result.push(self.edge(e).from().fix());
                }
            }
        }
        result
    }

    /// Estimates a normal value for a given vertex.
    ///
    /// This assumes that the triangulation models some kind of height field, given by the
    /// function `f`.
    /// The normal is the weighted and normalized average of the normals of all triangles
    /// adjacent to the given vertex.
    pub fn estimate_normal<F, RV>(&self, v: FixedVertexHandle, f: &F) -> RV
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
        RV: ThreeDimensional<Scalar = <V::Point as PointN>::Scalar>,
    {
        let mut v_pos = RV::new();
        let neighbor_positions: Vec<_> = {
            let handle = self.vertex(v);
            let v_2d = (*handle).position();
            *v_pos.nth_mut(0) = *v_2d.nth(0);
            *v_pos.nth_mut(1) = *v_2d.nth(1);
            *v_pos.nth_mut(2) = f(&*handle);

            handle
                .ccw_out_edges()
                .map(|e| {
                    let pos = (*e.to()).position();
                    let mut result = RV::new();
                    *result.nth_mut(0) = *pos.nth(0);
                    *result.nth_mut(1) = *pos.nth(1);
                    *result.nth_mut(2) = f(&*e.to());
                    result
                })
                .collect()
        };
        let mut final_normal = RV::new();
        for index in 0..neighbor_positions.len() {
            let p0 = neighbor_positions[index].clone();
            let p1 = neighbor_positions[(index + 1) % neighbor_positions.len()].clone();
            let d0 = v_pos.sub(&p0);
            let d1 = v_pos.sub(&p1);
            let normal = d0.cross(&d1);
            if *normal.nth(2) > zero() {
                final_normal = final_normal.add(&normal);
            }
        }
        // Normalize
        final_normal.div(final_normal.length2().sqrt())
    }

    /// Estimates and returns the gradient for a single vertex in this triangulation.
    pub fn estimate_gradient<F>(&self, v: FixedVertexHandle, f: &F) -> V::Point
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
    {
        use cgmath::Point3;
        let normal = self.estimate_normal::<_, Point3<_>>(v, f);
        // Calculate gradient from normal
        let mut gradient = V::Point::new();
        *gradient.nth_mut(0) = normal.x;
        *gradient.nth_mut(1) = normal.y;
        let g2 = gradient.length2();
        if g2 != zero() {
            let one = <V::Point as PointN>::Scalar::one();
            let d2 = one - normal.z * normal.z;
            let a2 = d2 / (one - d2);
            gradient = gradient.mul((a2 / g2).sqrt());
        }

        gradient
    }

    /// Interpolates a data point on this triangulation according to Sibson's c1 interpolant.
    ///
    /// The interpolation given by `nn_interpolation` is not differentiable at the triangulation's
    /// data points. Sibson introduced another interpolation scheme that takes the gradient of each
    /// data point into account and offers an interpolation that is differentiable (c1) at the data
    /// points.
    /// The interpolation needs to know the gradients of the points natural neighbors, though.
    /// Spade can estimate them automatically, see `estimate_gradient` and `estimate_gradients`.
    /// The value that should be interpolated is given by `f`, the gradient of a vertex must
    /// be given by `g`.
    ///
    /// # flatness
    /// An additional flatness factor determines how flat the triangulation will be around
    /// the datapoints. A flatness factor of 0.5 is the factor used in sibson's original interpolant.
    /// A flatness of 0.0 or lower is nearly identical to sibson's original interpolant
    /// (`nn_interpolation(..)`). A factor of (exactly) 1.0 should yield best performance since
    /// an exponentiation can be omitted.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    ///
    /// # use nalgebra::{Point2, Point3};
    /// use spade::delaunay::FloatDelaunayTriangulation;
    /// use spade::HasPosition;
    ///
    /// struct PointWithHeight {
    ///   point: Point2<f64>,
    ///   gradient: Point2<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Point = Point2<f64>;
    ///     fn position(&self) -> Point2<f64> {
    ///       self.point
    ///     }
    /// }
    ///
    /// impl PointWithHeight {
    ///   fn new(point: Point2<f64>, height: f64) -> PointWithHeight {
    ///     // Initialize the gradient to any value since it will be overwritten
    ///     PointWithHeight { point: point, height: height, gradient: Point2::new(0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
    ///   // Insert some points here... (skipped)
    ///   # delaunay.insert(PointWithHeight::new(Point2::new(0.0, 0.0), 5.));
    ///   # delaunay.insert(PointWithHeight::new(Point2::new(1.0, 0.0), 0.));
    ///   # delaunay.insert(PointWithHeight::new(Point2::new(0.0, 1.0), 2.));
    ///   // Estimate all gradients and store them:
    ///   delaunay.estimate_gradients(&(|v: &PointWithHeight| v.height),
    ///                               &(|v: &mut PointWithHeight, g| v.gradient = g));
    ///   
    ///   // Now we can use the gradients for interpolation, flatness is set to 2.0:
    ///   let interpolated = delaunay.nn_interpolation_c1_sibson(
    ///      &Point2::new(0.5, 0.2), 2.0, |v| v.height, |_, v| v.gradient);
    ///   println!("interpolated: {}", interpolated.unwrap());
    /// }
    /// ```
    pub fn nn_interpolation_c1_sibson<F, G>(
        &self,
        point: &V::Point,
        flatness: <V::Point as PointN>::Scalar,
        f: F,
        g: G,
    ) -> Option<<V::Point as PointN>::Scalar>
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
        G: Fn(&Self, &VertexHandle<V>) -> V::Point,
    {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        if ws.len() == 1 {
            return Some(f(&*self.vertex(*nns.first().unwrap())));
        }
        let mut sum_c0 = zero();
        let mut sum_c1 = zero();
        let mut sum_c1_weights = zero();
        let mut alpha = <V::Point as PointN>::Scalar::zero();
        let mut beta = <V::Point as PointN>::Scalar::zero();
        for (index, fixed_handle) in nns.iter().enumerate() {
            let handle = self.s.vertex(*fixed_handle);
            let pos_i = (*handle).position();
            let h_i = f(&*handle);
            let diff = pos_i.sub(point);
            let r_i2 = diff.length2();
            let r_i = r_i2.powf(flatness);
            let c1_weight_i = ws[index] / r_i;
            let grad_i = g(&self, &handle);
            let zeta_i = h_i + diff.dot(&grad_i);
            alpha += c1_weight_i * r_i;
            beta += c1_weight_i * r_i2;
            sum_c1_weights += c1_weight_i;
            sum_c1 += zeta_i * c1_weight_i;
            sum_c0 += h_i * ws[index];
        }
        alpha /= sum_c1_weights;
        sum_c1 /= sum_c1_weights;
        let result = (alpha * sum_c0 + beta * sum_c1) / (alpha + beta);
        Some(result)
    }

    /// Interpolates a data point on this triangulation using Farin's c1 interpolant.
    ///
    /// This method is used in the same way as `nn_interpolation_c1_sibson`. The resulting
    /// interpolant is very similar to sibson's c1 interpolant but uses a different algorithm.
    pub fn nn_interpolation_c1_farin<F, G>(
        &self,
        point: &V::Point,
        f: F,
        g: G,
    ) -> Option<<V::Point as PointN>::Scalar>
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
        G: Fn(&Self, &VertexHandle<V>) -> V::Point,
    {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        let handles: Vec<_> = nns.iter().map(|v| self.vertex(*v)).collect();

        let one = <V::Point as PointN>::Scalar::one();
        let two = one + one;
        let three = one + one + one;
        let four = two + two;
        let six = three * two;

        let get_edge_height = |i2: usize, i1: usize| {
            // Calculates an edge control point.
            // The edge is given by handles[i2] -> handles[i1]
            let p1 = (*handles[i1]).position();
            let p2 = (*handles[i2]).position();
            let diff = p1.sub(&p2);
            let grad = g(&self, &handles[i2]);
            f(&*handles[i2]) - grad.dot(&diff) / three
        };

        let mut result = <V::Point as PointN>::Scalar::zero();

        for first_index in 0..nns.len() {
            for second_index in first_index..nns.len() {
                for third_index in second_index..nns.len() {
                    // Determine control point "height"
                    let control_height;
                    let norm;
                    if first_index == second_index && second_index == third_index {
                        // Control point is one of the original data points
                        control_height = f(&*handles[first_index]);
                        norm = six;
                    } else if first_index == second_index
                        || first_index == third_index
                        || second_index == third_index
                    {
                        // Point lies on an edge of the bezier simplex
                        let (i2, i1) = if first_index == second_index {
                            (first_index, third_index)
                        } else if first_index == third_index {
                            (first_index, second_index)
                        } else {
                            (second_index, first_index)
                        };
                        control_height = get_edge_height(i2, i1);
                        norm = two;
                    } else {
                        // We have an inner control point, first != second != third
                        // Get all 6 edge control points of the triangle spanned by
                        // first_index, second_index and third_index
                        let cur_handles = [first_index, second_index, third_index];
                        const EDGE_POINTS: [(usize, usize); 6] =
                            [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
                        let mut edge_contrib = <V::Point as PointN>::Scalar::zero();
                        for &(ref i1, ref i2) in &EDGE_POINTS {
                            let i2 = cur_handles[*i2];
                            let i1 = cur_handles[*i1];
                            // TODO: We could cache the edge points in some way instead of
                            // calculating them anew for each inner point...
                            edge_contrib += get_edge_height(i2, i1);
                        }
                        edge_contrib /= four;

                        let inner_contrib = (f(&*handles[first_index])
                            + f(&*handles[second_index])
                            + f(&*handles[third_index]))
                            / six;
                        control_height = edge_contrib - inner_contrib;
                        norm = one;
                    }

                    // Add control height to result, weight it with the appropriate natural
                    // neighbor coordinates.
                    result +=
                        six * ws[first_index] * ws[second_index] * ws[third_index] * control_height
                            / norm;
                }
            }
        }
        Some(result)
    }
}

impl<V, K, L> DelaunayTriangulation<V, K, L>
where
    V: HasPosition2D,
    <V::Point as PointN>::Scalar: SpadeFloat,
    K: DelaunayKernel<<V::Point as PointN>::Scalar>,
    L: DelaunayLocateStructure<V::Point>,
    V::Point: TwoDimensional,
{
    /// Estimates a normal for each vertex in the triangulation.
    ///
    /// `f` must yield a 'height' value for each vertex,
    /// `g` is a callback function that can be used to store the calculated normals.
    ///
    /// ```
    /// # extern crate nalgebra;
    /// # extern crate spade;
    ///
    /// # use nalgebra::{Point2, Point3};
    /// use spade::delaunay::FloatDelaunayTriangulation;
    /// use spade::HasPosition;
    ///
    /// struct PointWithHeight {
    ///   point: Point2<f64>,
    ///   normal: Point3<f64>,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///   type Point = Point2<f64>;
    ///     fn position(&self) -> Point2<f64> {
    ///       self.point
    ///     }
    /// }
    /// impl PointWithHeight {
    ///   fn new(point: Point2<f64>, height: f64) -> PointWithHeight {
    ///     PointWithHeight { point: point, height: height, normal: Point3::new(0., 0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
    ///   // Insert some points here... (skipped)
    ///   # delaunay.insert(PointWithHeight { point: Point2::new(0.0, 0.0), height: 5., normal: Point3::new(0., 0., 0.)});
    ///   // Then, estimate all normals at once:
    ///   delaunay.estimate_normals(&(|v: &PointWithHeight| v.height), &(|v: &mut PointWithHeight, n| v.normal = n));
    ///   
    ///   // And print them
    ///   for vertex in delaunay.vertices() {
    ///      println!("vertex: {:?}, normal: {:?}", vertex.position(), vertex.normal);
    ///   }
    /// }
    /// ```
    pub fn estimate_normals<F, G, RV>(&mut self, f: &F, g: G)
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
        G: Fn(&mut V, RV),
        RV: ThreeDimensional<Scalar = <V::Point as PointN>::Scalar>,
    {
        for v in 0..self.num_vertices() {
            let normal = self.estimate_normal::<F, RV>(v, f);
            g(self.s.vertex_mut(v), normal);
        }
    }
    /// Estimates gradients for all vertices in this triangulation.
    ///
    /// `f` yields the value for which a gradient should be calculated,
    /// `g` can be used to store the gradient. See `estimate_normals` for a similar example.
    ///
    /// Internally, the normal for each vertex is calculated first. Then, an appropriate
    /// gradient is calculated.
    pub fn estimate_gradients<F, G>(&mut self, f: &F, g: &G)
    where
        F: Fn(&V) -> <V::Point as PointN>::Scalar,
        G: Fn(&mut V, V::Point),
    {
        for v in 0..self.num_vertices() {
            let gradient = self.estimate_gradient::<F>(v, f);
            g(self.s.vertex_mut(v), gradient);
        }
    }
}

#[cfg(test)]
mod test {
    use super::delaunay_basic::BasicDelaunaySubdivision;
    use super::{FloatDelaunayTriangulation, IntDelaunayTriangulation};
    use crate::testutils::*;
    use crate::traits::{HasPosition, SpatialObject};
    use cgmath::Point2;
    use rand::distributions::{Distribution, Range};
    use rand::{Rng, SeedableRng, XorShiftRng};

    #[test]
    fn test_insert_one_point() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        assert_eq!(d.num_vertices(), 0);
        d.insert(Point2::new(0f64, 0f64));
        assert_eq!(d.num_vertices(), 1);
        assert!(d.is_degenerate());
        d.sanity_check();
    }

    #[test]
    fn test_insert_three_points() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0f64, 0f64));
        d.insert(Point2::new(1f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));
        assert!(!d.is_degenerate());
        assert_eq!(d.num_vertices(), 3);
        d.sanity_check();
    }

    #[test]
    fn test_insert_three_points_cw() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));
        d.insert(Point2::new(1f64, 0f64));
        assert!(!d.is_degenerate());
        assert_eq!(d.num_vertices(), 3);
        d.sanity_check();
    }

    #[test]
    fn test_insert_four_points() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0f64, 0f64));
        d.insert(Point2::new(1f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));
        d.insert(Point2::new(1f64, 1f64));
        assert_eq!(d.num_vertices(), 4);
        d.sanity_check();
    }

    #[test]
    fn test_iterate_faces() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(-1f64, -1f64));
        d.insert(Point2::new(1f64, 1f64));
        d.insert(Point2::new(-1f64, 1f64));
        d.insert(Point2::new(1f64, -1f64));
        assert_eq!(d.triangles().count(), 2);
        const SIZE: usize = 1000;
        let points = random_points_with_seed::<f64>(SIZE, b"ilikechillychees");
        for p in points {
            d.insert(p);
        }
        assert_eq!(d.triangles().count(), SIZE * 2 + 2);
        for _ in 0..SIZE / 2 {
            d.remove(5);
        }
        assert_eq!(d.triangles().count(), SIZE + 2);
        for t in d.triangles() {
            let adj = t
                .adjacent_edge()
                .expect("Triangle must have an adjacent edge");
            let next = adj.o_next();
            let prev = adj.o_prev();
            // Check if the edge is really adjacent to a triangle
            assert_eq!(next.o_next(), prev);
            assert_eq!(prev.o_prev(), next);
            assert_eq!(next.to(), prev.from());
        }
    }

    #[test]
    fn test_insert_points() {
        // Just check if this won't crash
        const SIZE: usize = 10000;
        let mut points = random_points_with_seed::<f64>(SIZE, b"test_insert_poin");
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
        d.sanity_check();
    }

    #[test]
    fn test_insert_integral_points() {
        const SIZE: usize = 10000;
        let mut points = random_points_in_range(1000i64, SIZE, b"test insert inte");
        let mut d = IntDelaunayTriangulation::with_tree_locate();
        for p in points.drain(..) {
            d.insert(p);
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_outside_convex_hull() {
        const NUM: usize = 100;
        let mut rng = XorShiftRng::from_seed(*b"insert_outside_c");
        let range = Range::new(0., ::std::f64::consts::PI);
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for _ in 0..NUM {
            let ang = range.sample(&mut rng);
            let vec = Point2::new(ang.sin(), ang.cos()) * 100.;
            d.insert(vec);
        }
        assert_eq!(d.num_vertices(), NUM);
        d.sanity_check();
    }

    #[test]
    fn test_insert_same_point_small() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.2, 0.1),
        ];
        for p in &points {
            d.insert(*p);
        }
        for p in &points {
            d.insert(*p);
        }
        assert_eq!(d.num_vertices(), points.len());
        d.sanity_check();
    }

    #[test]
    fn test_insert_same_point() {
        const SIZE: usize = 300;
        let mut points = random_points_with_seed::<f64>(SIZE, b"insert same.poi!");
        let mut d = FloatDelaunayTriangulation::with_walk_locate();
        for p in &points {
            d.insert(*p);
        }
        for p in points.drain(..) {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
        d.sanity_check();
    }

    #[test]
    fn test_insert_point_on_ch_edge() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0., 0f64));
        d.insert(Point2::new(1., 0.));
        d.insert(Point2::new(0., 1.));

        d.insert(Point2::new(0., 0.4));
        d.sanity_check();
    }

    #[test]
    fn test_insert_on_edges() {
        // Just check if this won't crash
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0., 0f64));
        d.insert(Point2::new(1., 0.));
        d.insert(Point2::new(0., 1.));
        d.insert(Point2::new(1., 1.));
        d.insert(Point2::new(0.5, 0.5));
        d.insert(Point2::new(0.2, 0.2));
        d.insert(Point2::new(0., 0.4));
        d.insert(Point2::new(1., 0.5));
        d.insert(Point2::new(0.5, 1.));
        d.insert(Point2::new(0.7, 0.));
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0., 1f64));
        for i in -50..50 {
            d.insert(Point2::new(f64::from(i), 0.));
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line_2() {
        // This test inserts the line first
        let mut d = FloatDelaunayTriangulation::with_tree_locate();

        for i in -50..50 {
            d.insert(Point2::new(f64::from(i), 0.));
        }

        assert!(d.is_degenerate());

        for i in -10..10 {
            d.insert(Point2::new(f64::from(i), 0.5 * f64::from(i)));
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_grid() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0., 1f64));
        d.insert(Point2::new(0.0, 0.0));
        d.insert(Point2::new(1.0, 0.0));
        for y in 0..20 {
            for x in 0..7 {
                d.insert(Point2::new(f64::from(x), f64::from(y)));
            }
        }
        d.sanity_check();
    }

    #[test]
    fn test_create_edges_for_degenerate_triangulation() {
        let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
        delaunay.insert(Point2::new(0.0, 1.0));
        assert_eq!(delaunay.num_edges(), 0);
        let v0 = delaunay.insert(Point2::new(3.0, 0.0));
        assert_eq!(delaunay.num_edges(), 1);
        delaunay.sanity_check();
        delaunay.insert(Point2::new(-3.0, 2.0));
        delaunay.sanity_check();
        assert_eq!(delaunay.num_edges(), 2);
        delaunay.insert(Point2::new(9.0, -2.0));
        delaunay.sanity_check();
        assert_eq!(delaunay.num_edges(), 3);
        delaunay.insert(Point2::new(6.0, -1.0));
        assert_eq!(delaunay.insert(Point2::new(3.0, 0.0)), v0);
        assert_eq!(delaunay.num_edges(), 4);
        assert_eq!(delaunay.num_vertices(), 5);
        assert!(delaunay.is_degenerate());
        delaunay.sanity_check();
    }

    #[test]
    fn test_locate_for_degenerate_triangulation() {
        use super::PositionInTriangulation::*;
        let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
        let v0 = delaunay.insert(Point2::new(0.0, 1.0));
        let v1 = delaunay.insert(Point2::new(1.0, 0.0));
        delaunay.insert(Point2::new(2.0, -1.0));
        delaunay.insert(Point2::new(-1.0, 2.0));
        assert!(delaunay.is_degenerate());
        let v0 = delaunay.vertex(v0);
        let v1 = delaunay.vertex(v1);
        assert_eq!(delaunay.locate(&Point2::new(0.0, 1.0)), OnPoint(v0));
        assert_eq!(delaunay.locate(&Point2::new(1.0, 0.0)), OnPoint(v1));
        let edge = delaunay
            .get_edge_from_neighbors(v0.fix(), v1.fix())
            .unwrap();
        let locate = delaunay.locate(&Point2::new(0.5, 0.5));
        assert!(locate == OnEdge(edge) || locate == OnEdge(edge.sym()));
        let locate = delaunay.locate(&Point2::new(1.0, 1.0));
        assert!(locate == OutsideConvexHull(edge) || locate == OutsideConvexHull(edge.sym()));
        delaunay.sanity_check();
    }

    #[test]
    fn test_locate_for_simple_degenerate_triangulation() {
        use super::PositionInTriangulation::*;
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));
        assert_eq!(d.num_edges(), 1);
        let edge = d.edges().next().unwrap();
        let located = d.locate(&Point2::new(1.0, 0.0));
        assert!(located == OutsideConvexHull(edge) || located == OutsideConvexHull(edge.sym()));
        d.sanity_check();
    }

    struct PointWithHeight {
        point: Point2<f64>,
        height: f64,
    }

    impl HasPosition for PointWithHeight {
        type Point = Point2<f64>;
        fn position(&self) -> Point2<f64> {
            self.point
        }
    }

    impl PointWithHeight {
        fn new(x: f64, y: f64, height: f64) -> PointWithHeight {
            PointWithHeight {
                point: Point2::new(x, y),
                height,
            }
        }
    }

    #[test]
    fn test_natural_neighbor_interpolation() {
        let mut points = vec![
            PointWithHeight::new(0.0, 0.0, 0.0),
            PointWithHeight::new(1.0, 0.0, 0.0),
            PointWithHeight::new(0.0, 1.0, 0.0),
            PointWithHeight::new(1.0, 1.0, 0.0),
            PointWithHeight::new(2.0, 0.0, 0.0),
            PointWithHeight::new(2.0, 1.0, 0.0),
            PointWithHeight::new(3.0, 0.0, 1.0),
            PointWithHeight::new(3.0, 1.0, 1.0),
            PointWithHeight::new(4.0, 0.0, 1.0),
            PointWithHeight::new(4.0, 1.0, 1.0),
        ];
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in points.drain(..) {
            d.insert(point);
        }
        assert_eq!(
            d.nn_interpolation(&Point2::new(0.5, 0.5), |p| p.height),
            Some(0.0)
        );
        assert_eq!(
            d.nn_interpolation(&Point2::new(0.2, 0.8), |p| p.height),
            Some(0.0)
        );
        assert_eq!(
            d.nn_interpolation(&Point2::new(3.5, 1.), |p| p.height),
            Some(1.0)
        );
        assert_eq!(
            d.nn_interpolation(&Point2::new(-20., 0.2), |p| p.height),
            Some(0.0)
        );
        let height = d
            .nn_interpolation(&Point2::new(3.2, 0.9), |p| p.height)
            .unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        let height = d
            .nn_interpolation(&Point2::new(3.5, 0.5), |p| p.height)
            .unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        assert_eq!(
            d.nn_interpolation(&Point2::new(3.0, 0.0), |p| p.height),
            Some(1.0)
        );
    }

    #[test]
    fn test_insert_points_with_increasing_distance() {
        use cgmath::EuclideanSpace;
        let mut points = random_points_with_seed::<f64>(1000, b"with_increasing.");
        points.sort_by(|p1, p2| {
            p1.dot(p1.to_vec())
                .partial_cmp(&p2.dot(p2.to_vec()))
                .unwrap()
        });
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in &points {
            d.insert(*point);
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_grid_with_increasing_distance() {
        // This test inserts points on a grid with increasing distance from (0., 0.)
        use cgmath::EuclideanSpace;
        let mut points = Vec::new();
        const SIZE: i64 = 7;
        for x in -SIZE..SIZE {
            for y in -SIZE..SIZE {
                let point = Point2::new(x as f64, y as f64);
                points.push(point);
            }
        }
        points.sort_by(|p1, p2| {
            p1.dot(p1.to_vec())
                .partial_cmp(&p2.dot(p2.to_vec()))
                .unwrap()
        });
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in &points {
            d.insert(*point);
        }
        d.sanity_check();
    }

    #[test]
    fn test_remove_in_triangle() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(-1.0, 0.0f64));
        d.insert(Point2::new(1.0, 0.0f64));
        d.insert(Point2::new(0.0, 1.0f64));
        let to_remove = d.insert(Point2::new(0.0, 0.5));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 3);
        // Reinsert the last point, just to see if a crash occurs
        d.insert(Point2::new(0.0, 0.5));
        d.sanity_check();
    }

    #[test]
    fn test_remove_in_quad() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0.0, 0.0f64));
        d.insert(Point2::new(1.0, 0.0f64));
        d.insert(Point2::new(0.0, 1.0f64));
        d.insert(Point2::new(1.0, 1.0f64));
        let to_remove = d.insert(Point2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        let to_remove = d.insert(Point2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        d.insert(Point2::new(0.5, 0.6));
        d.sanity_check();
    }

    #[test]
    fn seven_point_test() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0.0, 0.0f64));
        d.insert(Point2::new(1.0, 0.0f64));
        d.insert(Point2::new(1.0, 1.0f64));
        d.insert(Point2::new(0.5, 1.5));
        d.insert(Point2::new(0.0, 1.0f64));
        d.insert(Point2::new(0.55, 0.5));
        let last = d.insert(Point2::new(0.4, 0.9));
        d.remove(last);
    }

    #[test]
    fn test_remove_inner() {
        use ::rand::{Rng, SeedableRng};

        let mut points = random_points_with_seed::<f64>(1000, b"remove inner!?&.");
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in &points {
            d.insert(*point);
        }
        // Insert an outer quad since we don't want to remove vertices from
        // the convex hull.
        d.insert(Point2::new(-2.0, -2.0));
        d.insert(Point2::new(-2.0, 2.0));
        d.insert(Point2::new(2.0, -2.0));
        d.insert(Point2::new(2.0, 2.0));
        // Now remove all inner points
        let mut rng = ::rand::XorShiftRng::from_seed(*b" next_seed/%&&2+");
        rng.shuffle(&mut points);
        assert_eq!(d.num_vertices(), 1004);
        for point in &points {
            let handle = d.locate_vertex(point).unwrap().fix();
            d.remove(handle);
        }
        assert_eq!(d.num_vertices(), 4);
        d.sanity_check();
    }

    #[test]
    fn test_remove_outer() {
        use cgmath::EuclideanSpace;
        let mut points = random_points_with_seed::<f64>(1000, b"-testremoveouter");
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in &points {
            d.insert(*point);
        }
        points.sort_by(|p1, p2| {
            p1.dot(p1.to_vec())
                .partial_cmp(&p2.dot(p2.to_vec()))
                .unwrap()
        });
        for point in points[3..].iter().rev() {
            let handle = d.locate_vertex(point).unwrap().fix();
            d.remove(handle);
        }
        d.sanity_check();
    }

    #[test]
    fn test_removal_and_insertion() {
        use cgmath::Point2;
        let points = random_points_with_seed::<f64>(1000, b"al and insertion");
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        for point in &points {
            d.insert(*point);
        }
        let mut rng = XorShiftRng::from_seed(*b"?i like fried&)%");
        for _ in 0..1000 {
            if rng.gen() {
                // Insert new random point
                let x = rng.gen();
                let y = rng.gen();
                d.insert(Point2::new(x, y));
            } else {
                // Remove random point
                let range = Range::new(0, d.num_vertices());
                let handle = range.sample(&mut rng);
                d.remove(handle);
            }
        }
        d.sanity_check();
    }

    #[test]
    fn test_remove_until_degenerate() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0., 0f64));
        d.insert(Point2::new(1., 0.));
        d.insert(Point2::new(0., 1.));
        d.insert(Point2::new(0., 0.5));
        d.insert(Point2::new(0., 0.25));
        d.insert(Point2::new(0., 0.75));
        assert_eq!(d.num_triangles(), 4);
        assert!(d.locate_and_remove(&Point2::new(1., 0.)).is_some());
        d.sanity_check();
        assert!(d.is_degenerate());
        while d.num_vertices() != 0 {
            d.remove(0);
            d.sanity_check();
        }
        assert!(d.is_degenerate());
        d.sanity_check();
        d.insert(Point2::new(0.5, 0.5));
        d.insert(Point2::new(0.2, 0.5));
        d.insert(Point2::new(1.5, 0.0));
        d.sanity_check();
    }

    #[test]
    fn test_remove_when_degenerate() {
        let mut d = FloatDelaunayTriangulation::with_tree_locate();
        d.insert(Point2::new(0.0, 0.0));
        let v1 = d.insert(Point2::new(1.0, 1.0));
        d.insert(Point2::new(2.0, 2.0));
        assert!(d.is_degenerate());
        assert_eq!(d.remove(v1), Point2::new(1.0, 1.0));
        d.sanity_check();
        let v0 = d.locate_vertex(&Point2::new(0.0, 0.0)).unwrap().fix();
        let v2 = d.locate_vertex(&Point2::new(2.0, 2.0)).unwrap().fix();
        assert!(d.get_edge_from_neighbors(v0, v2).is_some());
        assert_eq!(d.num_edges(), 1);
        assert_eq!(d.num_vertices(), 2);
    }

    #[test]
    fn test_nearest_neighbor_degenerate() {
        let mut d = FloatDelaunayTriangulation::with_walk_locate();
        assert!(d.nearest_neighbor(&Point2::new(3.0, 4.0)).is_none());
        let v0 = d.insert(Point2::new(0.0, 0.0));
        assert!(d.nearest_neighbor(&Point2::new(3.0, 4.0)).is_some());
        d.insert(Point2::new(0.0, 1.0));
        d.insert(Point2::new(0.0, 2.0));
        let v3 = d.insert(Point2::new(0.0, 3.0));
        let v4 = d.insert(Point2::new(0.0, 4.0));

        assert!(d.is_degenerate());

        let v0 = d.vertex(v0);
        let v3 = d.vertex(v3);
        let v4 = d.vertex(v4);
        assert_eq!(d.nearest_neighbor(&Point2::new(7.0, 3.2)), Some(v3));
        assert_eq!(d.nearest_neighbor(&Point2::new(4.0, 40.0)), Some(v4));
        assert_eq!(d.nearest_neighbor(&Point2::new(-4.0, -20.0)), Some(v0));
    }

    #[test]
    fn test_nearest_neighbor() {
        const SIZE: usize = 100;
        let points = random_points_with_seed::<f64>(SIZE, b"very nice test!_");
        let mut d = FloatDelaunayTriangulation::with_walk_locate();
        assert!(d.nearest_neighbor(&points[0]).is_none());
        for p in &points {
            d.insert(*p);
        }
        let sample_points = random_points_with_seed::<f64>(SIZE * 3, b"nearest_neighbor");
        for p in &sample_points {
            let nn_delaunay = d.nearest_neighbor(p);
            let nn_linear_search = points.iter().min_by(|l, r| {
                let d1 = l.distance2(p);
                let d2 = r.distance2(p);
                d1.partial_cmp(&d2).unwrap()
            });
            assert_eq!(nn_delaunay.map(|p| p.position()), nn_linear_search.cloned());
        }
    }

    #[test]
    #[cfg(feature = "serde_serialize")]
    fn test_serialization() {
        use crate::delaunay::{DelaunayTreeLocate, DelaunayWalkLocate};
        use serde_json;
        let mut t1 = IntDelaunayTriangulation::with_tree_locate();
        t1.insert([0i32, 12]);
        let mut t2 = IntDelaunayTriangulation::with_walk_locate();
        t2.insert([0i32, 22]);
        t2.insert([2, 10]);
        t2.insert([19, 29]);

        let json = serde_json::to_string(&t2).unwrap();
        let parsed: IntDelaunayTriangulation<[i32; 2], DelaunayWalkLocate> =
            serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.num_vertices(), 3);

        let json = serde_json::to_string(&t1).unwrap();
        let parsed: IntDelaunayTriangulation<[i32; 2], DelaunayTreeLocate<[i32; 2]>> =
            serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.num_vertices(), 1);
    }
}
