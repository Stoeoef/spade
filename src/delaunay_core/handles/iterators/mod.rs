use doc_comment::doc_comment;

mod circular_iterator;
mod fixed_iterators;
mod hull_iterator;

use super::{handle_defs::*, InnerTag, PossiblyOuterTag};
use fixed_iterators::{DynamicHandleIterator, FixedHandleIterator};

pub use circular_iterator::{CircularIterator, NextBackFn};
pub use hull_iterator::HullIterator;

macro_rules! type_handle_iterator {
    ($doc_name:expr, $item_type:expr, $def:item) => {
        doc_comment! {
            concat!("An iterator over ", $doc_name, ".\n",
            "\n",
            "This iterator is created by [Triangulation::", $item_type, "()]",
            "(../trait.Triangulation.html#method.", $item_type, ")."
        ),
            $def
        }
    };
}

// Vertex iterators
type_handle_iterator!(
    "vertices",
    "vertices",
    pub type VertexIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, VertexTag, InnerTag>;
);

type_handle_iterator!(
    "fixed vertices",
    "fixed_vertices",
    pub type FixedVertexIterator = FixedHandleIterator<VertexTag, InnerTag>;
);

type_handle_iterator!(
    "directed edges",
    "directed_edges",
    pub type DirectedEdgeIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, DirectedEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "fixed directed edges",
    "fixed_directed_edges",
    pub type FixedDirectedEdgeIterator = FixedHandleIterator<DirectedEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "undirected edges",
    "undirected_edges",
    pub type UndirectedEdgeIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, UndirectedEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "voronoi faces",
    "voronoi_faces",
    pub type VoronoiFaceIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, VoronoiFaceTag, PossiblyOuterTag>;
);

type_handle_iterator!(
    "directed voronoi edges",
    "directed_voronoi_edges",
    pub type DirectedVoronoiEdgeIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, DirectedVoronoiEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "undirected voronoi edges",
    "undirected_voronoi_edges",
    pub type UndirectedVoronoiEdgeIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, UndirectedVoronoiEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "fixed undirected edges",
    "fixed_undirected_edges",
    pub type FixedUndirectedEdgeIterator = FixedHandleIterator<UndirectedEdgeTag, InnerTag>;
);

type_handle_iterator!(
    "faces",
    "faces",
    pub type FaceIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, FaceTag, PossiblyOuterTag>;
);

type_handle_iterator!(
    "inner faces",
    "inner_faces",
    pub type InnerFaceIterator<'a, V, DE, UE, F> =
        DynamicHandleIterator<'a, V, DE, UE, F, FaceTag, InnerTag>;
);

type_handle_iterator!(
    "fixed faces",
    "fixed_faces",
    pub type FixedFaceIterator = FixedHandleIterator<FaceTag, PossiblyOuterTag>;
);

type_handle_iterator!(
    "fixed inner faces",
    "fixed_inner_faces",
    pub type FixedInnerFacesIterator = FixedHandleIterator<FaceTag, InnerTag>;
);
