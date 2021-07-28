use super::delaunay_core::DCEL;
use crate::{
    handles::VertexHandle, triangulation::TriangulationExt, HasPosition, HintGenerator,
    LastUsedVertexHintGenerator, Point2, Triangulation,
};
use doc_comment::doc_comment;

#[cfg(feature = "serde")]
use serde_crate::{Deserialize, Serialize};

doc_comment! {
concat!(
"A two dimensional [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation).

A Delaunay triangulation  a triangulation that fulfills the *Delaunay Property*: No
vertex of the triangulation is contained in the
[circumcircle](https://en.wikipedia.org/wiki/Circumscribed_circle) of any triangle.
As a consequence, Delaunay triangulations are well suited to support interpolation
algorithms and nearest neighbor searches. It is often constructed in order to retrieve its dual
graph, the [Voronoi diagram](#voronoi-diagram).\n

",
include_str!("../images/circumcircle.svg"),
"\n
 *An example triangulation with a few circumcircles drawn. No circumcircle contains more than three vertices.*

 Most methods on this type require the [Triangulation] trait. Refer to its documentation
 for more details on how to use `DelaunayTriangulation`.

 # Basic Usage
 Vertices need to implement the [HasPosition] trait. Spade bundles
 the [Point2](crate::Point2) struct for basic use cases.

 ## Basic example
 ```
 use spade::{DelaunayTriangulation, Triangulation, Point2};

 let mut triangulation: DelaunayTriangulation<_> = DelaunayTriangulation::new();

 // Insert three vertices that span one triangle (face)
 triangulation.insert(Point2::new(0.0, 1.0));
 triangulation.insert(Point2::new(1.0, 1.0));
 triangulation.insert(Point2::new(0.5, -1.0));

 assert_eq!(triangulation.num_vertices(), 3);
 assert_eq!(triangulation.num_inner_faces(), 1);
 assert_eq!(triangulation.num_undirected_edges(), 3);
 ```
 ## Right handed and left handed coordinate systems
 For simplicity, all method names and their documentation assume that the underlying coordinate system
 is right handed (e.g. x axis points to the right, y axis points upwards). If a left handed system
 (lhs) is used, any term related to orientation needs to be reversed:
 - \"left\" becomes \"right\" (example: the face of a directed edge is on the right side for a lhs
 - \"counter clock wise\" becomes \"clockwise\" (example: the vertices of a face are returned in clock wise order for a lhs)
 
 <table>
 <tr><th>left handed system</th><th>right handed system</th></tr>
 <tr><td>",
 include_str!("../images/lhs.svg"), "</td><td>",include_str!("../images/rhs.svg"), " </td></tr></table>",
 "\n
 # Extracting geometry information
 Spade uses [handles](crate::handles) to extract the triangulation's geometry.
 Handles are usually retrieved by inserting a vertex or by iterating.
 
 ## Example
 ```
 use crate::spade::{DelaunayTriangulation, Triangulation, Point2};
 
 let mut triangulation: DelaunayTriangulation<Point2<f64>> = DelaunayTriangulation::new();
 
 triangulation.insert(Point2::new(0.0, 1.0));
 triangulation.insert(Point2::new(1.0, 1.0));
 triangulation.insert(Point2::new(0.5, -1.0));
 
 for face in triangulation.inner_faces() {
   // face is a FaceHandle
   // edges is an array containing 3 directed edge handles
   let edges = face.adjacent_edges();
   for edge in &edges {
     let from = edge.from();
     let to = edge.to();
     // from and to are vertex handles
     println!(\"found an edge: {:?} -> {:?}\", from, to);
   }
 
   // vertices is an array containing 3 vertex handles
   let vertices = face.vertices();
   for vertex in &vertices {
     println!(\"Found vertex with position {:?}\", vertex.position());
   }
 }
 ```
 
 # Type parameters
 The triangulation's vertices, edges and faces can contain custom data.
 By default, the edge and face types are set to `()`. The vertex type must
 be specified.
 
  * `V: HasPosition` The vertex type
  * `DE: Default` The directed edge type.
  * `UE: Default` The undirected edge type.
  * `F: Default` The face type.
 
 Only vertices can be inserted directly. Faces and edges are create via `Default`
 necessary. Usually, edge and face data will need to be modified in a separate pass.
 
 ## Example
 ```
 use crate::spade::{DelaunayTriangulation, Triangulation, Point2};
 
 // A custom undirected edge type used to cache the length of an edge
 #[derive(Default)]
 struct EdgeWithLength { length: f64 }
 
 // Creates a new triangulation with a custom undirected edge type
 let mut triangulation: DelaunayTriangulation<Point2<f64>, (), EdgeWithLength>
                          = DelaunayTriangulation::new();
 
 triangulation.insert(Point2::new(0.0, 1.0));
 triangulation.insert(Point2::new(1.0, 1.0));
 triangulation.insert(Point2::new(0.5, -1.0));
 
 for edge in triangulation.fixed_undirected_edges() {
   let positions = triangulation.undirected_edge(edge).positions();
   let length = positions[0].point_distance_2(positions[1]).sqrt();
   // Write length into the edge data
   triangulation.undirected_edge_data_mut(edge).length = length;
 }
 
 for edge in triangulation.undirected_edges() {
    let length = edge.data().length;
    assert!(length > 0.0);
 }
 ```
 
 # Outer face
 Every triangulation contains an *outer face* that is adjacent to all edges of the
 triangulation's convex hull. The outer face is even present for a triangulation without
 vertices. It is either referenced by calling [Triangulation::outer_face()] or
 [handles::OUTER_FACE](crate::handles::OUTER_FACE)
 ", include_str!("../images/outer_faces.svg"),
"\n
# Voronoi Diagram

the dual graph of the Delaunay triangulation is the *Voronoi diagram*. The Voronoi diagram
subdivides the plane into several *Voronoi cells* (called `VoronoiFace` by Spade for consistency)
which contain all points in the plane that share a common nearest neighbor.

Keep in mind that \"faces\" and \"vertices\" are swapped - an (inner) Voronoi *vertex*
corresponds to a single Delaunay *face*.
The position of an inner voronoi vertex is the *circumcenter* of its corresponding Delaunay
face.
",
include_str!("../images/basic_voronoi.svg"),
"\n
*A Delaunay triangulation (black lines) and its dual graph, the Voronoi diagram (blue lines)*

## Extracting the Voronoi Diagram
Spade defines various functions to extract information about the Voronoi diagram:

**Types**
 * [DirectedVoronoiEdge](./handles/type.DirectedVoronoiEdge.html)
 * [UndirectedVoronoiEdge](./handles/type.UndirectedVoronoiEdge.html)
 * [VoronoiVertex](./handles/type.VoronoiVertex.html)
 * [VoronoiFace](./handles/type.VoronoiFace.html)

**Iterators**
 * [Triangulation::directed_voronoi_edges()]
 * [Triangulation::undirected_voronoi_edges()]

 **Conversion**
 * [DirectedVoronoiEdge::as_undirected()](./handles/type.DirectedVoronoiEdge.html#method.as_undirected)
 * [UndirectedVoronoiEdge::as_directed()](./handles/type.UndirectedVoronoiEdge.html#method.as_directed)
 * [DirectedEdgeHandle::as_voronoi_edge()](./handles/type.DirectedEdgeHandle.html#method.as_voronoi_edge)
 * [DirectedVoronoiEdge::as_delaunay_edge()](./handles/type.DirectedVoronoiEdge.html#method.as_delaunay_edge)
 * [UndirectedEdgeHandle::as_voronoi_edge()](./handles/type.UndirectedEdgeHandle.html#method.as_voronoi_edge)
 * [UndirectedVoronoiEdge::as_delaunay_edge()](./handles/type.UndirectedVoronoiEdge.html#method.as_delaunay_edge)

## Extracting the Voronoi Diagram (Example)
Extracting the geometry of the voronoi diagram can be slightly tricky as some of the voronoi
extend into infinity (see the dashed lines in the example above).

```
use spade::{DelaunayTriangulation, Point2, Triangulation};
use spade::handles::{VoronoiVertex::Inner, VoronoiVertex::Outer};

// Prints out the location of all voronoi edges in a triangulation
fn log_voronoi_diagram(triangulation: &DelaunayTriangulation<Point2<f64>>) {
    for edge in triangulation.undirected_voronoi_edges() {
        match edge.vertices() {
            [Inner(from), Inner(to)] => {
                // \"from\" and \"to\" are inner faces of the Delaunay triangulation
                println!(\"Found voronoi edge between {:?} and {:?}\",
                    from.circumcenter(),
                    to.circumcenter());
            }
            [Inner(from), Outer(edge)] | [Outer(edge), Inner(from)] => {
                // Some lines don't have a finite end and extend into infinity.
                println!(\"Found infinite voronoi edge going out of {:?} into the direction {:?}\", 
                    from.circumcenter(), 
                    edge.direction_vector());
            }
            [Outer(_), Outer(_)] => {
                // This case only happens if all vertices of the triangulation lie on the
                // same line and can probably be ignored.
            }
        }
    }
}
```
# Performance considerations

The insertion and query speed of Delaunay triangulations heavily depends on the triangulations specific parameters and
the nature of the underlying point set. For further illustration, this section uses insertion of new vertices as 
an example. Lookup and nearest neighbor queries behave similarly.

In general, inserting a new vertex consists of two stages:
 - Finding the _site_ that a new vertex occupies. This is usually the face that contains the vertex.
 - Modifying the site to contain the new vertex.

The second step runs, on average and amortized, in constant time. For larger triangulations (a few hundred vertices 
suffice), the first step will quickly take the majority of time.
The naive approach for finding the site works by traversing the triangulation in the direction of the target position.
For a uniformly distributed set of points, this approach takes `O(sqrt(n))` time for each lookup.

The first step can be sped up significantly by providing a _hint_ - this is an arbitrary vertex of the triangulation
that is used as the starting vertex when traversing the insertion site. This can be done either explicitly by 
calling a method suffixed with `with_hint` (e.g. [Triangulation.insert_with_hint]) or by simply relying on the
triangulations [HintGenerator].

# See also
Delaunay triangulations are closely related to [constrained Delaunay triangulations](crate::ConstrainedDelaunayTriangulation)
"
),
#[doc(alias = "Voronoi")]
#[doc(alias = "Voronoi diagram")]
#[doc(alias = "Delaunay")]
#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(crate="serde_crate"))]
pub struct DelaunayTriangulation<V, DE = (), UE = (), F = (), L = LastUsedVertexHintGenerator>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    dcel: DCEL<V, DE, UE, F>,
    hint_generator: L,
}
}

impl<V, DE, UE, F, L> DelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    /// Returns the nearest neighbor of a given input vertex.
    ///
    /// Returns `None` if the triangulation is empty.
    ///
    /// # Runtime
    /// This method take O(sqrt(n)) on average where n is the number of vertices.
    pub fn nearest_neighbor(
        &self,
        position: Point2<<V as HasPosition>::Scalar>,
    ) -> Option<VertexHandle<V, DE, UE, F>> {
        if self.num_vertices() == 0 {
            return None;
        }

        let hint = self.hint_generator().get_hint(position);
        let hint = self.validate_vertex_handle(hint);

        Some(self.walk_to_nearest_neighbor(hint, position))
    }
}

impl<V, DE, UE, F, L> Default for DelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    fn default() -> Self {
        Self {
            dcel: Default::default(),
            hint_generator: Default::default(),
        }
    }
}

impl<V, DE, UE, F, L> Triangulation for DelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    type Vertex = V;
    type DirectedEdge = DE;
    type UndirectedEdge = UE;
    type Face = F;
    type HintGenerator = L;

    fn s(&self) -> &DCEL<V, DE, UE, F> {
        &self.dcel
    }

    fn s_mut(&mut self) -> &mut DCEL<V, DE, UE, F> {
        &mut self.dcel
    }

    fn hint_generator(&self) -> &Self::HintGenerator {
        &self.hint_generator
    }

    fn hint_generator_mut(&mut self) -> &mut Self::HintGenerator {
        &mut self.hint_generator
    }
}

impl<V, DE, UE, F, L> std::iter::FromIterator<V> for DelaunayTriangulation<V, DE, UE, F, L>
where
    V: HasPosition,
    DE: Default,
    UE: Default,
    F: Default,
    L: HintGenerator<<V as HasPosition>::Scalar>,
{
    fn from_iter<I: IntoIterator<Item = V>>(iterator: I) -> Self {
        let mut triangulation = Self::new();
        for vertex in iterator {
            triangulation.insert(vertex);
        }
        triangulation
    }
}

#[cfg(test)]
mod test {
    use std::iter::FromIterator;

    use crate::test_utilities::{random_points_with_seed, SEED};
    use crate::{DelaunayTriangulation, Point2, Triangulation};

    #[allow(unused)]
    #[cfg(all(feature = "serde_crate"))]
    // Just needs to compile
    fn check_serde<'a>() {
        use serde_crate::{Deserialize, Serialize};

        use crate::{HierarchyHintGenerator, LastUsedVertexHintGenerator, Point2};

        fn requires_serde<'de, T: Serialize + Deserialize<'de>>() {}

        type DT<L> = super::DelaunayTriangulation<Point2<f64>, (), (), (), L>;

        requires_serde::<'a, DT<LastUsedVertexHintGenerator>>();
        requires_serde::<'a, DT<HierarchyHintGenerator<f64>>>();
    }

    #[test]
    fn test_nearest_neighbor() {
        const SIZE: usize = 100;
        let points = random_points_with_seed(SIZE, SEED);
        let d = DelaunayTriangulation::<_>::from_iter(points.iter().copied());

        let sample_points = random_points_with_seed(SIZE * 3, SEED);
        for p in sample_points {
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
    fn test_nearest_neighbor_small() {
        let mut d = DelaunayTriangulation::<_>::new();
        assert_eq!(None, d.nearest_neighbor(Point2::new(0.0, 0.0)));

        d.insert(Point2::new(0.0, 0.0));
        assert!(d.nearest_neighbor(Point2::new(0.0, 1.0)).is_some());
    }
}
