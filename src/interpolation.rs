// TODO: This needs to be overhauled properly and split off
// ParameterizedDelaunayTriangulation
const INTERPOLATION_SMALLVEC_CAPACITY: usize = 8;

/*
impl<V, E, T, K> ParameterizedDelaunayTriangulation<V, E, T, K>
where
    V: HasPosition,
    E: Default,
    T: Default,
    K: DelaunayKernel,
{
    /// Performs a barycentric interpolation.
    /// Returns `None` if the triangulation has no triangles yet.
    /// Points outside of the convex hull will be interpolated as well.
    /// The other interpolation methods are used very similarly, check their
    /// documentation for an example.
    pub fn barycentric_interpolation<W>(&self, point: Point, w: W) -> Option<f64>
    where
        W: Fn(&V) -> f64,
    {
        let vertices: SmallVec<[_; 3]> = match self.locate(point) {
            PositionInTriangulation::OnVertex(v) => smallvec![v],
            PositionInTriangulation::OnEdge(e) => smallvec![self.edge(e).from(), self.edge(e).to()],
            PositionInTriangulation::OnFace(f) => {
                let vs = self.face(f).vertices();
                smallvec![vs[0], vs[1], vs[2]]
            }
            PositionInTriangulation::OutsideOfConvexHull(_) => unimplemented!(),
        };
        if vertices.len() == 1 {
            Some(w(&*vertices[0]))
        } else if vertices.len() == 2 {
            let p0 = vertices[0].position();
            let p1 = vertices[1].position();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(0.0, edge.project_point(point), 1.0);
            let w0 = 1.0 - w1;
            Some(w1 * w(&*vertices[1]) + w0 * w(&*vertices[0]))
        } else {
            let triangle = crate::primitives::SimpleTriangle::new(
                vertices[0].position(),
                vertices[1].position(),
                vertices[2].position(),
            );
            let b_coords = triangle.barycentric_interpolation(point);
            let w0 = w(&*vertices[0]);
            let w1 = w(&*vertices[1]);
            let w2 = w(&*vertices[2]);
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
    /// use spade::{DelaunayTriangulation, Point, HasPosition};
    ///
    /// struct PointWithHeight {
    ///   point: Point,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///     fn position(&self) -> Point {
    ///       self.point
    ///     }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = DelaunayTriangulation::new();
    ///   delaunay.insert(PointWithHeight { point: Point::new(0.0, 0.0), height: 5. });
    ///   delaunay.insert(PointWithHeight { point: Point::new(1.0, 0.0), height: 0. });
    ///   delaunay.insert(PointWithHeight { point: Point::new(0.0, 1.0), height: 0. });
    ///   let lookup = Point::new(0.2, 0.2);
    ///   // Interpolate the points height
    ///   let interpolated = delaunay.nn_interpolation(lookup, |p| p.height).unwrap();
    ///   // and insert it afterwards.
    ///   delaunay.insert(PointWithHeight { point: lookup, height: interpolated });
    ///   // Data points themselves will always yield their own height
    ///   assert_eq!(delaunay.nn_interpolation(Point::new(0.0, 0.0), |p| p.height),
    ///              Some(5.0));
    /// }
    /// ```
    pub fn nn_interpolation<W>(&self, point: Point, w: W) -> Option<f64>
    where
        W: Fn(&V) -> f64,
    {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);

        let mut sum = None;
        for (index, fixed_handle) in nns.iter().enumerate() {
            let new = w(&*self.s.vertex(*fixed_handle)) * ws[index];
            sum = match sum {
                Some(val) => Some(val + new),
                None => Some(new),
            }
        }
        sum
    }

    fn get_weights(
        &self,
        nns: &SmallVec<[FixedVertexHandle; INTERPOLATION_SMALLVEC_CAPACITY]>,
        point: Point,
    ) -> SmallVec<[f64; INTERPOLATION_SMALLVEC_CAPACITY]> {
        let mut result = SmallVec::new();
        if nns.len() == 1 {
            result.push(1.0);
            return result;
        }

        if nns.len() == 2 {
            let p0 = (*self.s.vertex(nns[0])).position();
            let p1 = (*self.s.vertex(nns[1])).position();
            let edge = SimpleEdge::new(p0, p1);
            let w1 = ::clamp::clamp(0.0, edge.project_point(point), 1.0);
            let w0 = 1.0 - w1;
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
            let triangle = SimpleTriangle::new(next, cur_pos, point);
            point_cell.push(triangle.circumcenter());
        }

        let mut total_area = 0.0;
        let mut ccw_edge = self
            .s
            .get_edge_from_neighbors(*nns.first().unwrap(), *nns.last().unwrap())
            .unwrap();
        for (index, cur) in nns.iter().enumerate() {
            // Calculate area of voronoi cells
            let cur_pos = (*self.s.vertex(*cur)).position();
            let next = nns[(index + 1) % len];

            let mut polygon_area = 0.0;
            let mut last = point_cell[((index + len) - 1) % len];
            let first = point_cell[index];
            loop {
                if ccw_edge.to().fix() == next {
                    break;
                }
                let cw_edge = ccw_edge.cw();
                let triangle = SimpleTriangle::new(
                    (*ccw_edge.to()).position(),
                    (*cw_edge.to()).position(),
                    cur_pos,
                );
                let cur = triangle.circumcenter();

                let tri = SimpleTriangle::new(first, cur, last);
                last = cur;
                polygon_area += tri.double_area();
                ccw_edge = cw_edge;
            }
            ccw_edge = ccw_edge.rev();

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
        position: Point,
    ) -> SmallVec<[FixedVertexHandle; INTERPOLATION_SMALLVEC_CAPACITY]> {
        match self.locate_with_hint_option_fixed(position, None) {
            PositionInTriangulation::OnFace(face) => {
                let mut edges: SmallVec<_> = self
                    .face(face)
                    .adjacent_edges()
                    .iter()
                    .rev()
                    .map(|e| e.rev().fix())
                    .collect();
                self.inspect_flips(&mut edges, position)
            }
            PositionInTriangulation::OnEdge(edge) => {
                let edge = self.edge(edge);
                if self.is_ch_edge(edge.fix()) || self.is_ch_edge(edge.rev().fix()) {
                    let mut vec = SmallVec::new();
                    vec.push(edge.from().fix());
                    vec.push(edge.to().fix());
                    vec
                } else {
                    let mut edges = SmallVec::new();
                    edges.push(edge.prev().rev().fix());
                    edges.push(edge.next().rev().fix());
                    edges.push(edge.rev().prev().rev().fix());
                    edges.push(edge.rev().next().rev().fix());
                    self.inspect_flips(&mut edges, position)
                }
            }
            PositionInTriangulation::OnVertex(fixed_handle) => {
                let mut result = SmallVec::new();
                result.push(fixed_handle);
                result
            }
            PositionInTriangulation::OutsideOfConvexHull(_) => {
                unimplemented!();
            }
        }
    }

    fn inspect_flips(
        &self,
        edges: &mut SmallVec<[FixedDirectedEdgeHandle; INTERPOLATION_SMALLVEC_CAPACITY]>,
        position: Point,
    ) -> SmallVec<[FixedVertexHandle; INTERPOLATION_SMALLVEC_CAPACITY]> {
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
                    e1 = edge.prev().rev().fix();
                    e2 = edge.next().rev().fix();
                }
                debug_assert!(K::is_ordered_ccw(v1, v2, v0));
                debug_assert!(K::is_ordered_ccw(position, v1, v0));
                if K::contained_in_circumference(v2, v1, v0, position) {
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

    /// Estimates and returns the gradient for a single vertex in this triangulation.
    pub fn estimate_gradient<W>(&self, v: FixedVertexHandle, w: W) -> Gradient
    where
        W: Fn(&V) -> f64,
    {
        let normal = self.estimate_normal(v, w);
        // Calculate gradient from normal
        let mut flat = Point::new(normal.x, normal.y);
        let length = flat.length2();
        if length != 0.0 {
            let d2 = 1.0 - normal.z * normal.z;
            let a2 = d2 / (1.0 - d2);
            flat = flat.mul((a2 / length).sqrt());
        }
        Gradient {
            x: flat.x,
            y: flat.y,
        }
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
    /// the data points. A flatness factor of 0.5 is the factor used in sibson's original interpolant.
    /// A flatness of 0.0 or lower is nearly identical to sibson's original interpolant
    /// (`nn_interpolation(..)`). A factor of (exactly) 1.0 should yield best performance since
    /// an exponentiation can be omitted.
    ///
    /// # Example
    ///
    /// ```
    /// use spade::{DelaunayTriangulation, Point, Gradient, HasPosition};
    ///
    /// struct PointWithHeight {
    ///   point: Point,
    ///   gradient: Gradient,
    ///   height: f64,
    /// }
    ///
    /// impl HasPosition for PointWithHeight {
    ///     fn position(&self) -> Point {
    ///       self.point
    ///     }
    /// }
    ///
    /// impl PointWithHeight {
    ///   fn new(point: Point, height: f64) -> PointWithHeight {
    ///     // Initialize the gradient to any value since it will be overwritten
    ///     PointWithHeight { point: point, height: height, gradient: Gradient::new(0., 0.) }
    ///   }
    /// }
    ///
    /// fn main() {
    ///   let mut delaunay = DelaunayTriangulation::new();
    ///   // Insert some points here... (skipped)
    ///   # delaunay.insert(PointWithHeight::new(Point::new(0.0, 0.0), 5.));
    ///   # delaunay.insert(PointWithHeight::new(Point::new(1.0, 0.0), 0.));
    ///   # delaunay.insert(PointWithHeight::new(Point::new(0.0, 1.0), 2.));
    ///   // Estimate all gradients and store them:
    ///   for vertex in delaunay.vertices_fixed() {
    ///     let gradient = delaunay.estimate_gradient(vertex, |v| v.height);
    ///     delaunay.vertex_mut(vertex).gradient = gradient;
    ///   }
    ///
    ///   // Now we can use the gradients for interpolation, flatness is set to 2.0:
    ///   let interpolated = delaunay.nn_interpolation_c1_sibson(
    ///      Point::new(0.5, 0.2), 2.0, |v| v.height, |v| v.gradient);
    ///   println!("interpolated: {}", interpolated.unwrap());
    /// }
    /// ```
    pub fn nn_interpolation_c1_sibson<W, G>(
        &self,
        point: Point,
        flatness: f64,
        w: W,
        g: G,
    ) -> Option<f64>
    where
        W: Fn(&V) -> f64,
        G: Fn(&V) -> Gradient,
    {
        let nns = self.get_natural_neighbors(point);
        let ws = self.get_weights(&nns, point);
        if ws.is_empty() {
            return None;
        }
        if ws.len() == 1 {
            return Some(w(&*self.vertex(*nns.first().unwrap())));
        }
        let mut sum_c0 = 0.0;
        let mut sum_c1 = 0.0;
        let mut sum_c1_weights = 0.0;
        let mut alpha = 0.0;
        let mut beta = 0.0;
        for (index, fixed_handle) in nns.iter().enumerate() {
            let handle = self.s.vertex(*fixed_handle);
            let pos_i = (*handle).position();
            let h_i = w(&*handle);
            let diff = pos_i.sub(point);
            let r_i2 = diff.length2();
            let r_i = r_i2.powf(flatness);
            let c1_weight_i = ws[index] / r_i;
            let gradient = g(&*handle);
            let gradient = Point {
                x: gradient.x,
                y: gradient.y,
            };
            let zeta_i = h_i + diff.dot(gradient);
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
}

impl<V, E, T, K> DelaunayTriangulation<V, E, T, K>
where
    V: HasPosition,
    E: Default,
    T: Default,
    K: DelaunayKernel,
{
    /// Estimates a normal value for a given vertex.
    ///
    /// This assumes that the triangulation models some kind of height field, given by the
    /// function `f`.
    /// The normal is the weighted and normalized average of the normals of all triangles
    /// adjacent to the given vertex.
    pub fn estimate_normal<W>(&self, v: FixedVertexHandle, w: W) -> Normal
    where
        W: Fn(&V) -> f64,
    {
        let handle = self.vertex(v);
        let v_2d = (*handle).position();
        let v_pos = Normal::new(v_2d.x, v_2d.y, w(&*handle));

        let neighbor_positions: Vec<_> = {
            handle
                .out_edges()
                .map(|e| {
                    let pos = e.to().position();
                    Normal::new(pos.x, pos.y, w(&*e.to()))
                })
                .collect()
        };
        let mut final_normal = Normal::default();
        for index in 0..neighbor_positions.len() {
            let p0 = neighbor_positions[index];
            let p1 = neighbor_positions[(index + 1) % neighbor_positions.len()];
            let d0 = v_pos.sub(p0);
            let d1 = v_pos.sub(p1);
            let normal = d0.cross(d1);
            if normal.z > 0.0 {
                final_normal = final_normal.add(normal);
            }
        }
        // Normalize
        final_normal.mul(1. / final_normal.length2().sqrt())
    }
}

*/
