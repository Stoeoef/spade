This document intends to give a very high level overview of all major components of Spade.

Please feel free to open an issue if anything appears to be incorrect, missing or outdated.

# Handles

All handles are implemented by just two structs: `FixedHandleImpl<Type, InnerOuterMarker>` and `DynamicHandleImpl<'a, V, DE, UE, F, Type, InnerOuter: InnerOuterMarker>`.

Public users should not use these directly but use their various type shorthands instead (e.g. `type FixedVertexHandle = FixedHandleImpl<VertexTag, InnerTag>`)

The type parameters control which API is implemented for them:
 - `Type`: Specifies to which element of the triangulation the handle refers to, e.g. `DirectedEdgeTag`, `VoronoiFaceTag`
 - `InnerOuterMarker`: For faces, this can be either `InnerTag` or `PossiblyOuterTag`. For all other types, `InnerTag` is used.
 - The ubiquitous type parameters `V, DE, UE, F` refer to the vertex, directed edge, undirected edge and face handles of a DCEL.

Internally, fixed handles use `u32` to index the appropriate element in a vector. Using `u32` reduces the required space and performs much better in practice.
Dynamic handles store both a fixed handle and a reference to the DCEL they belong to.

# Directed edge and undirected edge handles

Directed edge handles that just differ in the smallest bit of their `u32` index refer to the same *undirected* edge and encode both of its possible directions.
E.g. indices 0 and 1 will refer to a directed edge and its `rev()`, indices 2 and 3 refer to another edge and its `rev()`.
Flipping the smallest bit of an edge index will also flip its direction.

Indices for undirected are shifted one bit to the right and "forget" the directional information.
I.e. the directed edge indices 26 and 27 both refer to the undirected edge with index 13.

This system has several advantages:
 - Any edge and it's `rev()` are always located close to each other in memory, leading to a good cache utilization.
 - Many operations (e.g. `rev()`, `as_directed()` and `as_undirected()`) can also be performed on fixed vertex handles - no intermediate lookup of the handle is required.

# Triangulation and TriangulationExt

`Triangulation` contains *publicly visible* functions that are common for both Cdts and unconstrained triangulations.

`TriangulationExt` is defined as `pub(crate)` and contains helper methods that are used by both Cdts and DTs to fulfill the Delaunay property.

# DCEL

`DCEL` (in `dcel.rs`) implements the underlying doubly connected edge list and all of its major accessor methods (e.g., for accessing vertices, edges and faces).
`dcel_operations.rs` implements several common operations on a DCEL, e.g. flipping an edge.

`DCEL`, and `dcel_operations` should, in contrast to `Triangulation`, not require a `V: HasPosition` bound anywhere. In other words: DCEL's only care about the *topology* of their elements, not about their actual position. The position is only required once the Delaunay property is evaluated (as part of `TriangulationExt`)

# Module dependency graph

This graph omits transitive dependencies and may not exactly coincide with the actual type dependencies in code. Smaller, less important modules are also left out. An arrow "A-->B" should be interpreted as "A calls methods / uses types from B"
<pre>
                        dcel&handles
               point     &iterators
                 ▲          ▲
                 │          │
               math   dcel_operations
                 ▲          ▲
                 │          │
                 └────┬─────┘
                      │
               ┌──────┴───────┐
               │              │
           bulk_load    triangulation_ext
               ▲              ▲
               │              │
               └──────┬───────┘
                      │
intersection_    triangulation
  iterator            ▲
     ▲                │
     └──────┐ ┌───────┴──────┐
            │ │              │
            cdt     delaunay_triangulation
</pre>