mod bulk_load;

#[cfg(test)]
mod bulk_load_fuzz_tests;

mod dcel;
pub mod dcel_operations;
mod delaunay_core_impl;
mod handles;
mod hint_generator;
mod line_side_info;
mod triangulation_ext;

pub mod math;

pub use bulk_load::bulk_load;

pub use triangulation_ext::{RemovalResult, TriangulationExt};

pub use dcel::Dcel;
pub use hint_generator::{
    HierarchyHintGenerator, HierarchyHintGeneratorWithBranchFactor, HintGenerator,
    LastUsedVertexHintGenerator,
};
pub use line_side_info::LineSideInfo;

pub use handles::iterators;
pub use handles::*;

pub use delaunay_core_impl::PositionInTriangulation;
pub use math::{
    validate_coordinate, validate_vertex, InsertionError, MAX_ALLOWED_VALUE, MIN_ALLOWED_VALUE,
};
