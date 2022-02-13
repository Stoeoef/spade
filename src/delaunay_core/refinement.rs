use std::collections::{HashSet, VecDeque};

use num_traits::Float;

use crate::{
    delaunay_core::math, CdtEdge, ConstrainedDelaunayTriangulation, HasPosition, HintGenerator,
    Point2, PositionInTriangulation, SpadeNum, Triangulation,
};

use super::{
    FaceHandle, FixedFaceHandle, FixedUndirectedEdgeHandle, FixedVertexHandle, InnerTag,
    TriangulationExt, UndirectedEdgeHandle,
};

/// Contains details about the outcome of a refinement procedure.
///
/// *See [ConstrainedDelaunayTriangulation::refine]*
pub struct RefinementResult {
    /// A hash set containing all excluded faces at the end of the triangulation.
    ///
    /// This set will be empty unless [RefinementParameters::exclude_outer_faces] has been used during refinement.
    /// In this case, the set contains the outer faces at the end of the triangulation, including any outer faces
    /// that were created during the refinement.
    pub excluded_faces: HashSet<FixedFaceHandle<InnerTag>>,

    /// Set to `true` if the refinement could be completed regularly.
    ///
    /// This will be `false` if the refinement ran out of additional vertices
    /// (see [RefinementParameters::with_max_additional_vertices]). Consider adapting the refinement parameters in this case,
    /// either by using a higher additional vertex count or by e.g. lowering the [angle limit](RefinementParameters::with_angle_limit).
    pub refinement_complete: bool,
}

/// Specifies the minimum allowed angle that should be kept after a refinement procedure.
///
/// The refinement algorithm will attempt to keep the *minimum angle in the triangulation* greater than
/// an angle limit specified with this struct.
///
/// *See [ConstrainedDelaunayTriangulation::refine], [RefinementParameters::with_angle_limit]*
#[derive(Copy, Clone, PartialEq, PartialOrd)]
pub struct AngleLimit {
    radius_to_shortest_edge_limit: f64,
}

impl AngleLimit {
    /// Create a new angle limit from an angle given in degrees.
    ///
    /// Note that angles larger than 30 degrees will quickly lead to overrefinement as the algorithm
    /// cannot necessarily guarantee termination (other than limiting the number of additional inserted vertices).
    ///
    /// Defaults to 30°. An angle of 0 degrees will disable refining due to small angles.
    ///
    /// *See also [from_rad](crate::AngleLimit::from_rad)*
    pub fn from_deg(degree: f64) -> Self {
        Self::from_rad(degree.to_radians())
    }

    /// Create a new angle limit from an angle given in radians.
    ///
    /// Note angles larger than 30 degrees (≈0.52rad = PI / 6) will quickly lead to poor refinement quality.
    /// Passing in an angle of 0rad will disable refining due to small angles.
    ///
    /// *See also [from_deg](crate::AngleLimit::from_deg)*
    pub fn from_rad(rad: f64) -> Self {
        let sin = rad.sin();
        if sin == 0.0 {
            Self::from_radius_to_shortest_edge_ratio(f64::INFINITY)
        } else {
            Self::from_radius_to_shortest_edge_ratio(0.5 / sin)
        }
    }

    /// Returns the radius to shortest edge limit corresponding to this angle limit.
    ///
    /// See [from_radius_to_shortest_edge_ratio](crate::AngleLimit::from_radius_to_shortest_edge_ratio) for more
    /// information.
    pub fn radius_to_shortest_edge_limit(&self) -> f64 {
        self.radius_to_shortest_edge_limit
    }

    /// Creates a new angle limit by specifying the circumradius to shortest edge ratio that must be kept.
    ///
    /// For each face, this ratio is calculated by dividing the circumradius of the face by the length of its shortest
    /// edge.
    /// This ratio is related directly to the minimum allowed angle by the formula
    /// `ratio = 1 / (2 sin * (min_angle))`.
    /// The *larger* the allowed min angle is, the *smaller* will the ratio become.
    ///
    /// Larger ratio values will lead to a less refined triangulation. Passing in `f64::INFINITY` will disable
    /// refining due to small angles.
    ///
    /// Defaults to 1.0 (30 degrees).
    ///
    /// # Example values
    ///
    /// | ratio | Bound on smallest angle (deg) | Bound on smallest angle (rad) |
    /// |-------|-------------------------------|-------------------------------|
    /// | 0.58  |                        60.00° |                          1.05 |
    /// | 0.60  |                        56.44° |                          0.99 |
    /// | 0.70  |                        45.58° |                          0.80 |
    /// | 0.80  |                        38.68° |                          0.68 |
    /// | 0.90  |                        33.75° |                          0.59 |
    /// | 1.00  |                        30.00° |                          0.52 |
    /// | 1.10  |                        27.04° |                          0.47 |
    /// | 1.20  |                        24.62° |                          0.43 |
    /// | 1.30  |                        22.62° |                          0.39 |
    /// | +INF  |                            0° |                             0 |
    pub fn from_radius_to_shortest_edge_ratio(ratio: f64) -> Self {
        Self {
            radius_to_shortest_edge_limit: ratio,
        }
    }
}

impl std::fmt::Debug for AngleLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AngleLimit")
            .field(
                "angle limit (deg)",
                &(0.5 / self.radius_to_shortest_edge_limit)
                    .asin()
                    .to_degrees(),
            )
            .finish()
    }
}

impl Default for AngleLimit {
    fn default() -> Self {
        Self::from_radius_to_shortest_edge_ratio(1.0)
    }
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Hash)]
pub enum RefinementHint {
    Ignore,
    ShouldRefine,
    MustRefine,
}

/// Controls how a refinement is performed.
///
/// Refer to [ConstrainedDelaunayTriangulation::refine] and methods implemented by this type for more details
/// about which parameters are supported.
///
/// The following parameters will be used by `Self::default` and [Self::new]:
/// * `exclude_outer_faces`: false
/// * `keep_constraint_edges`: false
/// * `min_required_area`: disabled - no lower area limit is used
/// * `max_allowed_area`: disabled - no upper area limit is used
/// * `angle_limit`: 30 degrees by default.
/// * `num_additional_vertices`: 10 times the number of vertices in the triangulation
///
/// # Example
///
/// ```
/// use spade::{AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters};
///
/// fn refine_cdt(cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>) {
///     let params = RefinementParameters::<f64>::new()
///         .exclude_outer_faces(&cdt)
///         .keep_constraint_edges()
///         .with_min_required_area(0.0001)
///         .with_max_allowed_area(0.5)
///         .with_angle_limit(AngleLimit::from_deg(25.0));
///
///     cdt.refine(params);
/// }
/// ```
#[derive(Debug, PartialEq, Clone)]
pub struct RefinementParameters<S: SpadeNum + Float> {
    max_additional_vertices: Option<usize>,

    angle_limit: AngleLimit,
    min_area: Option<S>,
    max_area: Option<S>,
    keep_constraint_edges: bool,
    excluded_faces: HashSet<FixedFaceHandle<InnerTag>>,
}

impl<S: SpadeNum + Float> Default for RefinementParameters<S> {
    fn default() -> Self {
        Self {
            max_additional_vertices: None,
            angle_limit: AngleLimit::from_radius_to_shortest_edge_ratio(1.0),
            min_area: None,
            max_area: None,
            excluded_faces: HashSet::new(),
            keep_constraint_edges: false,
        }
    }
}

impl<S: SpadeNum + Float> RefinementParameters<S> {
    /// Creates a new set of `RefinementParameters`.
    ///
    /// Refer to the [struct definition](Self) for more information.
    pub fn new() -> Self {
        Self::default()
    }

    /// Specifies the smallest allowed inner angle in a refined triangulation.
    ///
    /// The refinement algorithm will attempt to insert additional points (called steiner points) until the
    /// minimum angle is larger than the angle bound specified by this function.
    ///
    /// Defaults to 30 degrees.
    ///
    /// Note that angle limits much larger than 30 degrees may not always terminate successfully - check
    /// [RefinementResult::refinement_complete] to make sure that the angle limit did not lead to an overly strong
    /// refined triangulation.
    ///
    /// # Examples of different angle limits
    /// <table>
    /// <tr><th>0° (no angle refinement)</th><th>20°</th><th>30°</th><th>34°</th></tr>
    /// <tr><td>
    #[doc = concat!(include_str!("../../images/angle_limit_00.svg"), "</td><td>")]
    #[doc = concat!(include_str!("../../images/angle_limit_20.svg"), "</td><td>")]
    #[doc = concat!(include_str!("../../images/angle_limit_30.svg"), "</td><td>")]
    #[doc = concat!(include_str!("../../images/angle_limit_34.svg"), "</td></tr></table>")]
    ///
    /// *See also [ConstrainedDelaunayTriangulation::refine]*
    pub fn with_angle_limit(mut self, angle_limit: AngleLimit) -> Self {
        self.angle_limit = angle_limit;
        self
    }

    /// Specifies a lower bound for a triangles area.
    ///
    /// The algorithm will attempt to ignore any triangle with an area below this limit. This can also prevent an
    /// exhaustion of additionally available vertices (see [Self::with_max_additional_vertices]).
    ///
    /// Note that there is no guarantee that no face below this area bound will be kept intact - in some cases, a split
    /// will still be required to restore the triangulation's Delaunay property. Also, this value does not specify a lower
    /// bound for the smallest possible triangle in the triangulation.
    ///
    /// Should be set to something lower than [Self::with_max_allowed_area]. If this method is not called, no lower
    /// bound check will be performed.
    pub fn with_min_required_area(mut self, min_area: S) -> Self {
        self.min_area = Some(min_area);
        self
    }

    /// Specifies an upper bound for triangle areas in the triangulation.
    ///
    /// By default, the refinement tries to be conservative in how many vertices it adds. This will lead to an uneven
    /// triangle size distribution - areas with larger feature will contain fewer, larger triangles whereas regions with
    /// small features will contain more densely packed triangles.
    /// By specifying an upper area bound for triangles, the resulting triangle sizes can be brought evened out a little
    /// as any large triangle above the bound will be split into smaller parts.
    ///
    /// Should be set to something larger than `with_min_required_area`. If this method is not called, no upper area bound
    /// check will be performed.
    pub fn with_max_allowed_area(mut self, max_area: S) -> Self {
        self.max_area = Some(max_area);
        self
    }

    /// Specifies how many additional vertices may be inserted during Delaunay refinement.
    ///
    /// Refinement may, in some cases, fail to terminate if the angle limit is set too high
    /// (see [with_angle_limit](Self::with_angle_limit)). Simply stopping the refinement after a certain number of vertices
    /// has been inserted is an easy way to enforce termination. However, the resulting mesh may exhibit very poor quality
    /// in this case - some areas may have become overly refined, others might be overlooked completely. Consider changing
    /// the parameters (most notably the angle limit) if the refinement runs out of vertices.
    ///
    /// Use [RefinementResult::refinement_complete] to check if the refinement has completed successfully.
    pub fn with_max_additional_vertices(mut self, max_additional_vertices: usize) -> Self {
        self.max_additional_vertices = Some(max_additional_vertices);
        self
    }

    /// Prevents constraint edges from being split during refinement.
    ///
    /// By default, constraint edges may be split in order to restore the triangulation's Delaunay property.
    /// The resulting two new edges will *become* new constraint edges, hence the original shape outlined by
    /// constraint edges remains the same - no "gaps" or deviations are introduced.
    ///
    /// Enabling this option will, in general, reduce the quality of the resulting mesh - it is not necessarily
    /// Delaunay anymore and faces adjacent to long constraint edges may violate the configured [AngleLimit].
    pub fn keep_constraint_edges(mut self) -> Self {
        self.keep_constraint_edges = true;
        self
    }

    /// Allows to exclude outer faces from the refinement process.
    ///
    /// This is useful if the constraint edges form a *closed shape* with a clearly defined inner and outer part.
    /// Spade will determine inner and outer faces by identifying which faces can be reached from the outer face
    /// without "crossing" a constraint edge, similar to a flood fill algorithm.
    ///
    /// Any holes in the triangulation will also be excluded. More specifically, any point with an odd winding number
    /// is considered to be inner (see e.g. [Wikipedia](https://en.wikipedia.org/wiki/Point_in_polygon#Winding_number_algorithm)).
    ///
    /// Note that excluded faces may still be subdivided if a neighboring edge needs to be split. However, they will never be the
    /// *cause* for a subdivision - their angle and area is ignored.
    ///
    /// The resulting outer faces of the triangulation are returned by the call to [refine](ConstrainedDelaunayTriangulation::refine),
    /// see [RefinementResult::excluded_faces].
    ///
    /// # Example
    /// <table>
    /// <tr><th>Unrefined</th><th>Refined</th></tr>
    /// <tr><td>
    #[doc = concat!(include_str!("../../images/exclude_outer_faces_unrefined.svg"), "</td><td>",include_str!("../../images/exclude_outer_faces_refined.svg"), " </td></tr></table>")]
    ///
    /// *A refinement operation configured to exclude outer faces. All colored faces are considered outer faces and are
    /// ignored during refinement. Note that the inner part of the "A" shape forms a hole and is also excluded.*
    pub fn exclude_outer_faces<V: HasPosition, DE: Default, UE: Default, F: Default>(
        mut self,
        triangulation: &ConstrainedDelaunayTriangulation<V, DE, UE, F>,
    ) -> Self {
        if triangulation.all_vertices_on_line() {
            return self;
        }

        // Determine excluded faces by "peeling of" outer layers and adding them to an outer layer set.
        // This needs to be done repeatedly to also get inner "holes" within the triangulation
        let mut inner_faces = HashSet::new();
        let mut outer_faces = HashSet::new();

        let mut current_todo_list: Vec<_> =
            triangulation.convex_hull().map(|edge| edge.rev()).collect();
        let mut next_todo_list = Vec::new();

        let mut return_outer_faces = true;

        loop {
            while let Some(next_edge) = current_todo_list.pop() {
                let (list, face_set) = if next_edge.is_constraint_edge() {
                    (&mut next_todo_list, &mut inner_faces)
                } else {
                    (&mut current_todo_list, &mut outer_faces)
                };

                if let Some(inner) = next_edge.face().as_inner() {
                    if face_set.insert(inner.fix()) {
                        list.push(next_edge.prev().rev());
                        list.push(next_edge.next().rev());
                    }
                }
            }

            if next_todo_list.is_empty() {
                break;
            }
            std::mem::swap(&mut inner_faces, &mut outer_faces);
            std::mem::swap(&mut next_todo_list, &mut current_todo_list);

            return_outer_faces = !return_outer_faces;
        }

        self.excluded_faces = if return_outer_faces {
            outer_faces
        } else {
            inner_faces
        };

        self
    }

    pub(crate) fn get_refinement_hint<V, DE, UE, F>(
        &self,
        face: FaceHandle<InnerTag, V, DE, UE, F>,
    ) -> RefinementHint
    where
        V: HasPosition<Scalar = S>,
    {
        if let Some(max_area) = self.max_area {
            if face.area() > max_area {
                return RefinementHint::MustRefine;
            }
        }

        if let Some(min_area) = self.min_area {
            if face.area() < min_area {
                return RefinementHint::Ignore;
            }
        }

        let (_, length2) = face.shortest_edge();
        let (_, radius2) = face.circumcircle();

        let ratio2 = radius2 / length2;
        let angle_limit = self.angle_limit.radius_to_shortest_edge_limit;
        if ratio2.into() > angle_limit * angle_limit {
            RefinementHint::ShouldRefine
        } else {
            RefinementHint::Ignore
        }
    }
}

impl<V, DE, UE, F, L> ConstrainedDelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition + From<Point2<<V as HasPosition>::Scalar>>,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
    <V as HasPosition>::Scalar: Float,
{
    /// Refines a triangulation by inserting additional points to improve the quality of its mesh.
    ///
    /// *Mesh quality*, when applied to constrained delaunay triangulations (CDT), usually refers to how skewed its
    /// triangles are. A skewed triangle is a triangle with very large or very small (acute) inner angles.
    /// Some applications (e.g. interpolation and finite element methods) perform poorly in the presence of skewed triangles.
    ///
    /// Refining by inserting additional points (called "steiner points") may increase the minimum angle. The given
    /// [RefinementParameters] should be used to fine-tune the refinement behavior.
    ///
    /// # General usage
    ///
    /// The vertex type must implement `From<Point2<...>>` - otherwise, Spade cannot construct new steiner points at a
    /// certain location. The refinement itself happens *in place* and will result in a valid CDT.
    ///
    /// ```
    /// use spade::{ConstrainedDelaunayTriangulation, RefinementParameters, Point2, InsertionError, Triangulation};
    ///
    /// fn get_refined_triangulation(vertices: Vec<Point2<f64>>) ->
    ///     Result<ConstrainedDelaunayTriangulation<Point2<f64>>,  InsertionError>
    /// {
    ///     let mut cdt = ConstrainedDelaunayTriangulation::bulk_load(vertices)?;
    ///     let result = cdt.refine(RefinementParameters::default());
    ///     if !result.refinement_complete {
    ///         panic!("Refinement failed - I should consider using different parameters.")
    ///     }
    ///
    ///     return Ok(cdt)
    /// }
    /// ```
    ///
    /// # Example image
    ///
    /// <table>
    /// <tr><th>Unrefined</th><th>Refined</th></tr>
    /// <tr><td>
    #[doc = concat!(include_str!("../../images/unrefined.svg"), "</td><td>",include_str!("../../images/refined.svg"), " </td></tr></table>")]
    ///
    /// *A refinement example. The CDT on the left has some acute angles and skewed triangles.
    /// The refined CDT on the right contains several additional points that prevents such triangles from appearing while keeping
    /// all input vertices and constraint edges.*
    ///
    /// # Quality guarantees
    ///
    /// Refinement will ensure that the resulting triangulation fulfills a few properties:
    ///  - The triangulation's minimum angle will be larger than the angle specified by
    ///    [with_angle_limit](crate::RefinementParameters::with_angle_limit).<br>
    ///    *Exception*: Acute input angles (small angles between initial constraint edges) cannot be resolved as the constraint edges
    ///    must be kept intact. The algorithm will, for the most part, leave those places unchanged.<br>
    ///    *Exception*: The refinement will often not be able to increase the minimal angle much above 30 degrees as every newly
    ///    inserted steiner point may create additional skewed triangles.
    ///  - The refinement will fullfil the *Delaunay Property*: Every triangle's circumcenter will not contain any other vertex.<br>
    ///    *Exception*: Refining with [keep_constraint_edges](crate::RefinementParameters::keep_constraint_edges) cannot restore
    ///    the Delaunay property if doing so would require splitting a constraint edge.<br>
    ///    *Exception*: Refining with [exclude_outer_faces](crate::RefinementParameters::exclude_outer_faces) will not
    ///    ensure the Delaunay property of any outer face.
    ///  - Spade allows to specify a [maximum allowed triangle area](crate::RefinementParameters::with_max_allowed_area).
    ///    The algorithm will attempt to subdivide any triangle with an area larger than this, independent of its smallest angle.
    ///  - Spade allows to specify a [minimum required triangle area](crate::RefinementParameters::with_min_required_area).
    ///    The refinement will attempt to ignore any triangle with an area smaller than this parameter. This can prevent the
    ///    refinement algorithm from over-refining in some cases.
    ///
    /// # General limitations
    ///
    /// The algorithm may fail to terminate in some cases for a minimum angle limit larger than 30 degrees. Such a limit can
    /// result in an endless loop: Every additionally inserted point creates more triangles that need to be refined.
    ///
    /// To prevent this, spade limits the number of additionally inserted steiner points
    /// (see [RefinementParameters::with_max_additional_vertices]). However, this may leave the refinement in an incomplete state -
    /// some areas of the input mesh may not have been triangulated at all, some will be overly refined.
    /// Use [RefinementResult::refinement_complete] to identify if a refinement operation has succeeded without running out of
    /// vertices.
    ///
    /// To prevent this from happening, consider either lowering the minimum angle limit
    /// (see [RefinementParameters::with_angle_limit]) or introduce a
    /// [minimum required area](RefinementParameters::with_min_required_area).
    ///
    /// # References
    ///
    /// This is an adaption of the classical refinement algorithms introduced by Jim Ruppert and Paul Chew.
    ///
    /// For a good introduction to the topic, refer to the slides from a short course at the thirteenth and fourteenth
    /// International Meshing Roundtables (2005) by Jonathan Richard Shewchuk:
    /// <https://people.eecs.berkeley.edu/~jrs/papers/imrtalk.pdf>
    ///
    /// Wikipedia: <https://en.wikipedia.org/wiki/Delaunay_refinement>
    ///
    ///
    #[doc(alias = "Refinement")]
    #[doc(alias = "Delaunay Refinement")]
    pub fn refine(&mut self, mut parameters: RefinementParameters<V::Scalar>) -> RefinementResult {
        use PositionInTriangulation::*;

        let mut excluded_faces = std::mem::take(&mut parameters.excluded_faces);

        let mut legalize_edges_buffer = Vec::with_capacity(20);
        let mut forcibly_splitted_segments_buffer = Vec::with_capacity(5);

        let mut encroached_segment_candidates =
            VecDeque::with_capacity(self.num_constraints() + self.convex_hull_size());

        encroached_segment_candidates.extend(
            self.undirected_edges()
                .filter(|edge| {
                    if parameters.keep_constraint_edges {
                        edge.is_part_of_convex_hull()
                    } else {
                        Self::is_fixed_edge(*edge)
                    }
                })
                .map(|edge| edge.fix()),
        );

        let mut encroached_face_candidates: VecDeque<_> = self.fixed_inner_faces().collect();

        let num_initial_vertices: usize = self.num_vertices();
        let num_additional_vertices = parameters
            .max_additional_vertices
            .unwrap_or_else(|| num_initial_vertices * 10);
        let max_allowed_vertices = num_initial_vertices + num_additional_vertices;

        let is_original_vertex = |vertex: FixedVertexHandle| vertex.index() < num_initial_vertices;

        let mut refinement_complete = true;

        'outer: loop {
            if self.num_vertices() >= max_allowed_vertices {
                refinement_complete = false;
                break;
            }

            if let Some(forcibly_split_segment) = forcibly_splitted_segments_buffer.pop() {
                self.resolve_encroachment(
                    &mut encroached_segment_candidates,
                    &mut encroached_face_candidates,
                    forcibly_split_segment,
                    is_original_vertex,
                    &mut excluded_faces,
                );
                continue;
            }

            // Take from encroached segments
            if let Some(segment_candidate) = encroached_segment_candidates.pop_front() {
                for edge in segment_candidate.directed_edges() {
                    let edge = self.directed_edge(edge);

                    let is_protected = edge
                        .face()
                        .as_inner()
                        .map(|face| excluded_faces.contains(&face.fix()))
                        .unwrap_or(true);

                    if is_protected {
                        continue;
                    }

                    if let Some(opposite_position) = edge.opposite_position() {
                        if is_encroaching_edge(
                            edge.from().position(),
                            edge.to().position(),
                            opposite_position,
                        ) {
                            // The edge is encroaching
                            self.resolve_encroachment(
                                &mut encroached_segment_candidates,
                                &mut encroached_face_candidates,
                                segment_candidate,
                                is_original_vertex,
                                &mut excluded_faces,
                            );
                            continue 'outer;
                        }
                    }
                }

                continue;
            }

            // Take the next triangle candidate
            if let Some(face) = encroached_face_candidates.pop_front() {
                if excluded_faces.contains(&face) {
                    continue;
                }

                let face = self.face(face);

                let (shortest_edge, _) = face.shortest_edge();
                let is_fixed_angle = [shortest_edge.next(), shortest_edge.prev()]
                    .into_iter()
                    .all(|edge| Self::is_fixed_edge(edge.as_undirected()));

                let refinement_hint = parameters.get_refinement_hint(face);

                if refinement_hint == RefinementHint::Ignore {
                    continue;
                }

                if refinement_hint == RefinementHint::ShouldRefine
                    && is_fixed_angle
                    && (!is_original_vertex(shortest_edge.from().fix())
                        || !is_original_vertex(shortest_edge.to().fix()))
                {
                    // The fixed segment has already been split. Don't attempt to subdivide it any further, this
                    // is as good as we can get.
                    continue;
                }

                let locate_hint = face.vertices()[0].fix();

                forcibly_splitted_segments_buffer.clear();
                legalize_edges_buffer.clear();

                let circumcenter = face.circumcenter();

                match self.locate_with_hint(circumcenter, locate_hint) {
                    OnEdge(edge) => {
                        let edge = self.directed_edge(edge);
                        if edge.is_part_of_convex_hull() {
                            if parameters.keep_constraint_edges && edge.is_constraint_edge() {
                                continue;
                            }

                            forcibly_splitted_segments_buffer.push(edge.fix().as_undirected());
                        } else {
                            legalize_edges_buffer.extend([
                                edge.next().fix(),
                                edge.prev().fix(),
                                edge.rev().next().fix(),
                                edge.rev().prev().fix(),
                            ]);
                        }
                    }
                    OnFace(face_under_circumcenter) => {
                        legalize_edges_buffer.extend(
                            self.face(face_under_circumcenter)
                                .adjacent_edges()
                                .map(|edge| edge.fix()),
                        );
                    }
                    OutsideOfConvexHull(_) => continue,
                    OnVertex(_) => continue,
                    NoTriangulation => unreachable!(),
                };

                let mut is_encroaching = false;

                while let Some(edge) = legalize_edges_buffer.pop() {
                    let edge = self.directed_edge(edge);
                    let [from, to] = edge.as_undirected().positions();
                    if Self::is_fixed_edge(edge.as_undirected()) {
                        if is_encroaching_edge(from, to, circumcenter) {
                            is_encroaching = true;

                            if !parameters.keep_constraint_edges || !edge.is_constraint_edge() {
                                // New circumcenter would encroach a constraint edge. Don't insert the circumcenter
                                // but force splitting the segment
                                forcibly_splitted_segments_buffer.push(edge.as_undirected().fix());
                            }
                        }
                        continue;
                    }

                    // Edge is not a constraint edge. Check if it needs to be legalized.
                    // We've already checked that this edge is not part of the convex hull - unwrap is safe
                    let opposite = edge.rev().opposite_position().unwrap();
                    let from = edge.from().position();
                    let to = edge.to().position();
                    let should_flip =
                        math::contained_in_circumference(opposite, to, from, circumcenter);

                    if should_flip {
                        let e1 = edge.rev().next().fix();
                        let e2 = edge.rev().prev().fix();

                        legalize_edges_buffer.push(e1);
                        legalize_edges_buffer.push(e2);
                    }
                }

                if !is_encroaching {
                    // The circumcenter doesn't encroach any segment. Continue really inserting it.
                    let new_vertex = self
                        .insert_with_hint(circumcenter.into(), locate_hint)
                        .expect("Failed to insert circumcenter, likely due to loss of precision. Consider refining with fewer additional vertices.");
                    encroached_face_candidates.extend(
                        self.vertex(new_vertex)
                            .out_edges()
                            .flat_map(|edge| edge.face().fix().as_inner()),
                    );
                } else if !forcibly_splitted_segments_buffer.is_empty() {
                    // Revisit this face later
                    encroached_face_candidates.push_back(face.fix());
                }
            } else {
                break;
            }
        }
        RefinementResult {
            excluded_faces,
            refinement_complete,
        }
    }

    fn is_fixed_edge(edge: UndirectedEdgeHandle<V, DE, CdtEdge<UE>, F>) -> bool {
        edge.is_constraint_edge() || edge.is_part_of_convex_hull()
    }

    fn resolve_encroachment(
        &mut self,
        encroached_segments_buffer: &mut VecDeque<FixedUndirectedEdgeHandle>,
        encroached_faces_buffer: &mut VecDeque<FixedFaceHandle<InnerTag>>,
        encroached_edge: FixedUndirectedEdgeHandle,
        is_original_vertex: impl Fn(FixedVertexHandle) -> bool,
        excluded_faces: &mut HashSet<FixedFaceHandle<InnerTag>>,
    ) {
        let segment = self.directed_edge(encroached_edge.as_directed());

        let [v0, v1] = segment.vertices();

        let half = Into::<V::Scalar>::into(0.5f32);

        let (w0, w1) = match (is_original_vertex(v0.fix()), is_original_vertex(v1.fix())) {
            (true, true) => {
                // Split the segment exactly in the middle if it has not been split before.
                (half, half)
            }
            (false, false) => {
                // If both ends are steiner points. We'll also split in the middle - this seems to
                // slightly improve the average angle
                (half, half)
            }
            (v0_is_original, _) => {
                let half_length = segment.length_2().sqrt() * half;

                // Round down to the nearest power of two to prevent runaway encroachment
                let nearest_power_of_two = nearest_power_of_two(half_length);
                let other_vertex_weight = half * nearest_power_of_two / half_length;
                let original_vertex_weight = Into::<V::Scalar>::into(1.0) - other_vertex_weight;

                if v0_is_original {
                    (original_vertex_weight, other_vertex_weight)
                } else {
                    (other_vertex_weight, original_vertex_weight)
                }
            }
        };

        let final_position = v0.position().mul(w0).add(v1.position().mul(w1));

        let is_left_side_protected = segment
            .face()
            .as_inner()
            .map(|face| excluded_faces.contains(&face.fix()));

        let is_right_side_protected = segment
            .rev()
            .face()
            .as_inner()
            .map(|face| excluded_faces.contains(&face.fix()));

        let is_constraint_edge = segment.is_constraint_edge();

        // Split the edge at its center
        let segment = segment.fix();
        let (new_vertex, [e1, e2]) = self.insert_on_edge(segment, final_position.into());

        if is_constraint_edge {
            self.handle_legal_edge_split([e1, e2].map(|edge| edge.as_undirected()));
        }

        let (h1, h2) = (self.directed_edge(e1), self.directed_edge(e2));

        if is_left_side_protected == Some(true) {
            excluded_faces.insert(h1.face().fix().as_inner().unwrap());
            excluded_faces.insert(h2.face().fix().as_inner().unwrap());
        }

        if is_right_side_protected == Some(true) {
            excluded_faces.insert(h1.rev().face().fix().as_inner().unwrap());
            excluded_faces.insert(h2.rev().face().fix().as_inner().unwrap());
        }

        self.legalize_vertex(new_vertex);

        encroached_faces_buffer.extend(
            self.vertex(new_vertex)
                .out_edges()
                .flat_map(|edge| edge.face().fix().as_inner()),
        );

        // Update encroachment candidates
        encroached_segments_buffer.push_back(e1.as_undirected());
        encroached_segments_buffer.push_back(e2.as_undirected());

        // Neighboring edges may have become encroached. Check if they need to be added to the encroached segment buffer.
        encroached_segments_buffer.extend(
            [e1.rev(), e2]
                .into_iter()
                .map(|edge| self.directed_edge(edge))
                .flat_map(|edge| [edge.next(), edge.rev().prev()].into_iter())
                .filter(|edge| Self::is_fixed_edge(edge.as_undirected()))
                .map(|edge| edge.as_undirected().fix()),
        );
    }
}

fn is_encroaching_edge<S: SpadeNum + Float>(
    edge_from: Point2<S>,
    edge_to: Point2<S>,
    query_point: Point2<S>,
) -> bool {
    let edge_center = edge_from.add(edge_to).mul(0.5f32.into());
    let radius_2 = edge_from.distance_2(edge_to) * 0.25.into();

    query_point.distance_2(edge_center) < radius_2
}

fn nearest_power_of_two<S: Float + SpadeNum>(input: S) -> S {
    input.log2().round().exp2()
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use crate::{
        test_utilities::{random_points_with_seed, SEED},
        AngleLimit, ConstrainedDelaunayTriangulation, InsertionError, Point2, RefinementParameters,
        Triangulation as _,
    };

    pub type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;

    #[test]
    fn test_zero_angle_limit_dbg() {
        let limit = AngleLimit::from_deg(0.0);
        let debug_string = format!("{:?}", limit);
        assert_eq!(debug_string, "AngleLimit { angle limit (deg): 0.0 }");
    }

    #[test]
    fn test_zero_angle_limit() -> Result<(), InsertionError> {
        let limit = AngleLimit::from_deg(0.0);

        assert_eq!(limit.radius_to_shortest_edge_limit(), f64::INFINITY);

        let mut vertices = random_points_with_seed(20, SEED);

        // Insert an artificial outer boundary that will prevent the convex hull from being encroached.
        // This should prevent any refinement.
        vertices.push(Point2::new(100.0, 100.0));
        vertices.push(Point2::new(100.0, 0.0));
        vertices.push(Point2::new(100.0, -100.0));
        vertices.push(Point2::new(0.0, -100.0));
        vertices.push(Point2::new(-100.0, -100.0));
        vertices.push(Point2::new(-100.0, 0.0));
        vertices.push(Point2::new(-100.0, 100.0));
        vertices.push(Point2::new(0.0, 100.0));

        let mut cdt = Cdt::bulk_load(vertices)?;

        let initial_num_vertices = cdt.num_vertices();
        cdt.refine(RefinementParameters::new().with_angle_limit(limit));

        assert_eq!(initial_num_vertices, cdt.num_vertices());

        cdt.refine(RefinementParameters::new());
        assert!(initial_num_vertices < cdt.num_vertices());

        Ok(())
    }

    #[test]
    fn test_nearest_power_of_two() {
        use super::nearest_power_of_two;

        for i in 1..400u32 {
            let log = (i as f64).log2() as u32;
            let nearest = nearest_power_of_two(i as f64);

            assert!((1 << log) as f64 == nearest || (1 << (log + 1)) as f64 == nearest);
        }

        assert_eq!(0.25, nearest_power_of_two(0.25));
        assert_eq!(0.25, nearest_power_of_two(0.27));
        assert_eq!(0.5, nearest_power_of_two(0.5));
        assert_eq!(1.0, nearest_power_of_two(0.75));
        assert_eq!(2.0, nearest_power_of_two(1.5));
        assert_eq!(2.0, nearest_power_of_two(2.5));
        assert_eq!(4.0, nearest_power_of_two(3.231));
        assert_eq!(4.0, nearest_power_of_two(4.0));
    }

    #[test]
    fn test_small_refinement() -> Result<(), InsertionError> {
        let vertices = random_points_with_seed(20, SEED);
        let mut cdt = Cdt::bulk_load(vertices)?;

        let mut peekable = cdt.fixed_vertices().peekable();
        while let (Some(p0), Some(p1)) = (peekable.next(), peekable.peek()) {
            if !cdt.can_add_constraint(p0, *p1) {
                cdt.add_constraint(p0, *p1);
            }
        }

        cdt.refine(Default::default());

        cdt.cdt_sanity_check();

        Ok(())
    }

    #[test]
    fn test_sharp_angle_refinement() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();

        // This sharp angle should only be subdivided once rather than infinitely often.
        cdt.add_constraint_edge(Point2::new(-40.0, 80.0), Point2::new(0.0, 0.0))?;
        cdt.add_constraint_edge(Point2::new(0.0, 0.0), Point2::new(4.0, 90.0))?;

        let refinement_parameters = RefinementParameters::new()
            .with_max_additional_vertices(10)
            .with_angle_limit(AngleLimit::from_radius_to_shortest_edge_ratio(1.0));

        let result = cdt.refine(refinement_parameters);

        assert!(result.refinement_complete);
        cdt.cdt_sanity_check();
        Ok(())
    }

    #[test]
    fn test_simple_exclude_outer_faces() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();
        let vertices = [
            Point2::new(-1.0, 1.0),
            Point2::new(0.0, 0.5),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, -1.0),
            Point2::new(-1.0, -1.0),
            Point2::new(-1.0, 1.0),
        ];

        for points in vertices.windows(2) {
            let p1 = points[0];
            let p2 = points[1];
            cdt.add_constraint_edge(p1, p2)?;
        }

        let excluded_faces = RefinementParameters::<f64>::new()
            .exclude_outer_faces(&cdt)
            .excluded_faces;

        let v2 = cdt.locate_vertex(Point2::new(1.0, 1.0)).unwrap().fix();
        let v0 = cdt.locate_vertex(Point2::new(-1.0, 1.0)).unwrap().fix();

        let excluded_face = cdt
            .get_edge_from_neighbors(v2, v0)
            .and_then(|edge| edge.face().as_inner())
            .unwrap()
            .fix();

        assert_eq!(
            excluded_faces,
            HashSet::from_iter(std::iter::once(excluded_face))
        );

        Ok(())
    }

    #[test]
    fn test_exclude_outer_faces() -> Result<(), InsertionError> {
        let mut cdt = Cdt::new();

        let shape = [
            Point2::new(-1.0, 1.0),
            Point2::new(0.0, 0.7),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
            Point2::new(0.5, -1.0),
            Point2::new(-0.5, -2.0),
            Point2::new(-1.0, 1.0),
        ];

        let scale_factors = [1.0, 2.0, 3.0, 4.0];

        let mut current_excluded_faces = RefinementParameters::<f64>::new()
            .exclude_outer_faces(&cdt)
            .excluded_faces;

        assert!(current_excluded_faces.is_empty());

        for factor in scale_factors {
            for points in shape.windows(2) {
                let p1 = points[0].mul(factor);
                let p2 = points[1].mul(factor);

                cdt.add_constraint_edge(p1, p2)?;
            }

            let next_excluded_faces = RefinementParameters::<f64>::new()
                .exclude_outer_faces(&cdt)
                .excluded_faces;

            current_excluded_faces = next_excluded_faces;
            assert!(current_excluded_faces.len() < cdt.num_inner_faces());
        }

        Ok(())
    }
}
