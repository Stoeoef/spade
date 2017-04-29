use delaunay::*;
use self::delaunay_basic::{BasicDelaunaySubdivision, HasSubdivision};
use self::dcel::*;
use traits::HasPosition2D;
use point_traits::{PointN, PointNExtensions, TwoDimensional};
use kernels::DelaunayKernel;
use primitives::SimpleEdge;

/// A two dimensional constrained delaunay triangulation.
///
/// A constrained delaunay triangulation is a triangulation that
/// can contain _constraint edges_. These edges will be present
/// in the resulting triangulation, the resulting triangulation
/// does not necessarily fulfill the delaunay property.
///
/// This implementation currently supports only _weakly intersecting_
/// constraints, thus, constraint edges are allowed to touch at
/// their start or end point, but are not allowed to intersect at
/// any interior point.
pub struct ConstrainedDelaunayTriangulation<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    s: DCEL<V>,
    all_points_on_line: bool,
    constraints: Vec<bool>,
    num_constraints: usize,
    __marker: ::std::marker::PhantomData<K>,
}

struct ConflictRegion {
    left_hull: Vec<(FixedVertexHandle, FixedVertexHandle)>,
    right_hull: Vec<(FixedVertexHandle, FixedVertexHandle)>,
    conflicts: Vec<FixedEdgeHandle>,
    constraint_end: FixedVertexHandle,
}

impl<V, K> BasicDelaunaySubdivision<V, K> for ConstrainedDelaunayTriangulation<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    fn is_defined_legal(&self, edge: FixedEdgeHandle) -> bool {
        self.constraints[edge]
    }
}

impl<V, K> HasSubdivision<V, K> for ConstrainedDelaunayTriangulation<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    fn s(&self) -> &DCEL<V> {
        &self.s
    }

    fn s_mut(&mut self) -> &mut DCEL<V> {
        &mut self.s
    }
}

impl<V, K> ConstrainedDelaunayTriangulation<V, K>
    where V: HasPosition2D,
          V::Point: TwoDimensional,
          K: DelaunayKernel<<V::Point as PointN>::Scalar>
{
    pub fn new() -> ConstrainedDelaunayTriangulation<V, K> {
        ConstrainedDelaunayTriangulation {
            s: DCEL::new(),
            all_points_on_line: true,
            __marker: Default::default(),
            constraints: Vec::new(),
            num_constraints: 0,
        }
    }

    pub fn insert(&mut self, vertex: V) -> FixedVertexHandle {
        let handle = self.insert_with_hint(vertex, 0);
        self.constraints.resize(self.s.num_edges() * 2, false);
        handle
    }

    pub fn num_constraints(&self) -> usize {
        self.num_constraints
    }

    /// Returns `true` if a given edge is a constraint edge.
    pub fn is_constraint_edge(&self, edge: FixedEdgeHandle) -> bool {
        self.constraints[edge]
    }

    /// Checks if two vertices are connected by a constraint edge.
    pub fn exists_constraint(&self, from: FixedVertexHandle, to: FixedVertexHandle) -> bool {
        self.get_edge_from_vertices(from, to)
            .map(|e| self.is_constraint_edge(e.fix()))
            .unwrap_or(false)
    }
    
    fn insert_with_hint(&mut self, t: V, hint: FixedVertexHandle) -> FixedVertexHandle {
        let pos = t.position();
        let position_in_triangulation = if self.all_points_on_line {
            PositionInTriangulation::NoTriangulationPresent
        } else {
            self.locate_with_hint_fixed(&pos, hint)
        };
        let insertion_result = match position_in_triangulation {
            PositionInTriangulation::OutsideConvexHull(edge) => {
                Result::Ok(self.insert_outside_convex_hull(edge, t))
            }
            PositionInTriangulation::InTriangle(face) => {
                Result::Ok(self.insert_into_triangle(face, t))
            }
            PositionInTriangulation::OnEdge(edge) => {
                // Check if the edge is a constraint edge and split the
                // constraint if necessary
                let is_constraint = self.is_constraint_edge(edge);
                let (from, to);
                {
                    let edge = self.s.edge(edge);
                    from = edge.from().fix();
                    to = edge.to().fix();
                }
                let new_handle = self.insert_on_edge(edge, t);
                self.constraints.resize(self.s.num_edges() * 2, false);
                if is_constraint {
                    let handles = {
                        let e1 = self.get_edge_from_vertices(from, new_handle).unwrap();
                        let e2 = self.get_edge_from_vertices(new_handle, to).unwrap();
                        assert!(self.is_constraint_edge(e1.fix()) ^ self.is_constraint_edge(e2.fix()));
                        [e1.fix(), e1.sym().fix(), e2.fix(), e2.sym().fix()]
                    };
                    for h in &handles {
                        self.constraints[*h] = true;
                    }
                    self.num_constraints += 1;
                }
                Result::Ok(new_handle)
            },
            PositionInTriangulation::OnPoint(vertex) => {
                self.s.update_vertex(vertex, t);
                Result::Err(vertex)
            }
            PositionInTriangulation::NoTriangulationPresent => self.initial_insertion(t),
        };
        match insertion_result {
            Result::Ok(new_handle) => new_handle,
            Result::Err(update_handle) => update_handle,
        }
    }

    fn initial_insertion(&mut self, t: V) -> Result<FixedVertexHandle, FixedVertexHandle> {
        assert!(self.all_points_on_line);
        // Inserts points if no points are present or if all points
        // lie on the same line
        let new_pos = t.position();
        for vertex in self.s.fixed_vertices() {
            let pos = (*self.s.vertex(vertex)).position();
            if pos == new_pos {
                self.s.update_vertex(vertex, t);
                return Result::Err(vertex);
            }
        }

        if self.s.num_vertices() <= 1 {
            return Result::Ok(self.s.insert_vertex(t));
        }

        // Check if the new point is on the same line as all points in the
        // triangulation
        let from = (*self.s.vertex(0)).position();
        let to = (*self.s.vertex(1)).position();
        let edge = SimpleEdge::new(from.clone(), to.clone());
        if K::side_query(&edge, &new_pos).is_on_line() {
            return Result::Ok(self.s.insert_vertex(t));
        }
        // The point does not lie on the same line as all other points.
        // Start creating a triangulation
        let dir = to.sub(&from);
        let mut vertices: Vec<_> = self.s.vertices()
            .map(|v| (v.fix(), dir.dot(&(*v).position())))
            .collect();
        // Sort vertices according to their position on the line
        vertices.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());

        // Create line
        let is_ccw = K::is_ordered_ccw(&new_pos, &from, &to);
        let mut last_edge = self.s.connect_two_isolated_vertices(vertices[0].0, vertices[1].0, 0);
        let mut edges = vec![last_edge];
        for v in vertices.iter().skip(2) {
            let edge = self.s.connect_edge_to_isolated_vertex(last_edge, v.0);
            edges.push(self.s.edge(edge).fix());
            last_edge = edge;
        }
        if is_ccw {
            edges.reverse();
        }
        let new_vertex = self.s.insert_vertex(t);
        // Connect all points on the line to the new vertex
        let mut last_edge = *edges.first().unwrap();
        if !is_ccw {
            last_edge = self.s.edge(last_edge).sym().fix();
        }
        last_edge = self.s.connect_edge_to_isolated_vertex(last_edge, new_vertex);
        for e in edges {
            let e = if !is_ccw {
                self.s.edge(e).sym().fix()
            } else {
                e
            };
            last_edge = self.s.create_face(last_edge, e);
            last_edge = self.s.edge(last_edge).sym().fix();
        }
        self.all_points_on_line = false;
        Result::Ok(new_vertex)
    }

    /// Insert two points and creates a constraint between them.
    ///
    /// Returns `true` if the constraint has not yet existed.
    /// This method will only return `false` if both given vertices
    /// have already been inserted and connected by a constraint edge.
    ///
    /// # Panics
    /// Panics if the new constraint edge intersects an existing 
    /// constraint edge.
    pub fn add_new_constraint_edge(&mut self, from: V, to: V) -> bool {
        let from_handle = self.insert(from);
        let to_handle = self.insert(to);
        self.add_constraint(from_handle, to_handle)
    }

    /// Adds a constraint edge between to vertices.
    ///
    /// Returns `true` if a new constraint was added.
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
        loop {
            let mut region = self.get_conflict_region(cur_from, to);
            let cur_to = region.constraint_end;
            if let Some((edge, sym)) = from_neighbors(&self.s, cur_from, cur_to)
                .map(|e| (e.fix(), e.sym().fix()))
            {
                // The edge already exists - probably mark it as a constraint edge
                assert_eq!(self.constraints[edge], self.constraints[sym]);
                if !self.constraints[edge] {
                    // The constraint edge does not already exists
                    self.num_constraints += 1;
                    self.constraints[edge] = true;
                    self.constraints[sym] = true;
                    result |= true;
                }
            } else {
                self.resolve_conflict_region(&mut region);
                result |= true;
            }
            if cur_to == to {
                break;
            }
            cur_from = cur_to;
        }
        result
    }

    fn resolve_conflict_region(&mut self, region: &mut ConflictRegion) {
        let &mut ConflictRegion {
            ref left_hull,
            ref right_hull,
            ref mut conflicts,
            ..
        } = region;
        assert!(!conflicts.is_empty());
        {
            // Assert that no conflict edge is a constraint edge
            assert!(conflicts.iter().all(|e| !self.constraints[*e]),
                    "Constraint intersects another constraint edge");
            // Remove edges from highest to lowest value to
            // prevent any index change of a conflicting edge
            conflicts.sort();
            // Delete edges
            for handle in conflicts.iter().rev().cloned() {
                let dual = self.s.edge(handle).sym().fix();
                self.s.remove_edge(handle, None);
                if handle < dual {
                    self.constraints.swap_remove(dual);
                    self.constraints.swap_remove(handle);
                } else {
                    self.constraints.swap_remove(handle);
                    self.constraints.swap_remove(dual);
                }
            }
        }

        self.constraints.truncate(self.s.num_edges() * 2);

        let prev_edge = right_hull[right_hull.len() - 1];
        let prev_edge = from_neighbors(&self.s, prev_edge.0, prev_edge.1).unwrap().fix();
        let next_edge = right_hull[0];
        let next_edge = from_neighbors(&self.s, next_edge.0, next_edge.1).unwrap().fix();

        let constraint_edge = self.s.create_face(prev_edge, next_edge);
        // Push constraint, both for constraint_edge and its dual
        self.constraints.push(true);
        self.constraints.push(true);

        // Retriangulate the areas
        let mut edges = Vec::new();
        for &(v0, v1) in right_hull {
            let edge = from_neighbors(&self.s, v0, v1).unwrap().fix();
            edges.push(edge);
        }
        edges.push(constraint_edge);
        self.fill_hole(edges);
        let mut edges = Vec::new();
        for &(v0, v1) in left_hull.iter().rev() {
            let edge = from_neighbors(&self.s, v0, v1).unwrap().fix();
            edges.push(edge);
        }
        edges.push(self.edge(constraint_edge).sym().fix());
        self.fill_hole(edges);
        self.constraints.resize(self.s.num_edges() * 2, false);
        self.num_constraints += 1;
    }
    
    fn get_conflict_region(&self,
                         v0: FixedVertexHandle,
                         v1: FixedVertexHandle)
                           -> ConflictRegion {
        let mut new_conflicts = Vec::new();
        let mut right_hull = Vec::new();
        let mut left_hull = Vec::new();
        let mut constraint_end = v1;
        let rotate_vertex = self.vertex(v0);
        let target = self.vertex(v1).position();
        let constraint_edge = SimpleEdge::new(rotate_vertex.position(), target.clone());

        let mut start_edge = None;
        for edge in rotate_vertex.ccw_out_edges() {
            if edge.face().fix() == 0 {
                // Make sure this edge is not part of the convex hull
                continue;
            }
            let simple = Self::to_simple_edge(edge.o_next());
            if simple.intersects_edge_non_colinear::<K>(&constraint_edge) {
                // Check if the edge just touches the constraint edge
                let from_query = constraint_edge.side_query::<K>(&simple.from);
                if from_query.is_on_line() {
                    return ConflictRegion {
                        conflicts: Vec::new(),
                        right_hull: Vec::new(),
                        left_hull: Vec::new(),
                        constraint_end: edge.to().fix(),
                    };
                }
                let to_query = constraint_edge.side_query::<K>(&simple.to);
                if to_query.is_on_line() {
                    return ConflictRegion {
                        conflicts: Vec::new(),
                        right_hull: Vec::new(),
                        left_hull: Vec::new(),
                        constraint_end: edge.o_next().to().fix(),
                    };
                }
                start_edge = Some(edge.o_next());
                right_hull.push((v0, edge.to().fix()));
                left_hull.push((edge.o_prev().from().fix(), v0));
                break;
            }
        }
        let mut cur_edge = start_edge.expect("Could not find start edge. This is a bug.").sym();
        loop {
            new_conflicts.push(cur_edge.fix());
            assert!(Self::to_simple_edge(cur_edge).side_query::<K>(&target).is_on_left_side());
            let e_prev = cur_edge.o_prev();
            let o_next = cur_edge.o_next();
            if e_prev.from().fix() == v1 {
                // We've reached the target point
                right_hull.push((cur_edge.to().fix(), v1));
                left_hull.push((v1, cur_edge.from().fix()));
                break;
            }
            // One of the edges of the left face must intersect the constraint edge,
            // find out which
            let e_prev_inter = Self::to_simple_edge(e_prev)
                .intersects_edge_non_colinear::<K>(&constraint_edge);
            let o_next_inter = Self::to_simple_edge(o_next)
                .intersects_edge_non_colinear::<K>(&constraint_edge);
            if e_prev_inter && o_next_inter {
                // Both edges intersect - this mean the intersection edge is cutting through
                // the common point of o_next and e_prev.
                // This splits the constraint into multiple parts
                constraint_end = e_prev.from().fix();
                println!("Found different inner constraint end");
                break;
            } else if e_prev_inter {
                // o_prev is intersecting, add o_next to right hull
                right_hull.push((cur_edge.to().fix(), e_prev.from().fix()));
                cur_edge = e_prev.sym();
            } else {
                // o_next is intersecting, add o_prev to left hull
                left_hull.push((cur_edge.o_prev().from().fix(), cur_edge.from().fix()));
                cur_edge = cur_edge.o_next().sym();
            }
        }
        ConflictRegion {
            left_hull: left_hull,
            right_hull: right_hull,
            conflicts: new_conflicts,
            constraint_end: constraint_end,
        }
    }

    #[cfg(test)]
    fn cdt_sanity_check(&self) {
        let count_constraints = self.constraints.iter().filter(|&c| *c).count();
        let num_constraints = self.num_constraints() * 2;
        assert_eq!(count_constraints, num_constraints);
        for edge in self.s.edges() {
            if self.is_constraint_edge(edge.fix()) {
                assert!(self.is_constraint_edge(edge.sym().fix()));
            }
        }
    }
}

#[cfg(test)]
mod test {
    use testutils::*;
    use super::ConstrainedDelaunayTriangulation;
    use super::{DelaunayTriangulation, DelaunayWalkLocate};
    use traits::HasPosition;
    use kernels::FloatKernel;
    use cgmath::Point2;

    type CDT = ConstrainedDelaunayTriangulation<Point2<f64>, FloatKernel>;
    type Delaunay = DelaunayTriangulation<Point2<f64>, FloatKernel, DelaunayWalkLocate>;

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
        // Use a delaunay triangulation to "generate" non overlapping constraint edges
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
                let did_insert = cdt.add_constraint(h0, h1);
                if did_insert {
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
        assert!(cdt.add_constraint(v0, v_last));
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
}
