mod bulk_load;

#[cfg(test)]
mod bulk_load_fuzz_tests;

#[cfg(test)]
mod refinement_fuzz_tests;

mod dcel;
pub mod dcel_operations;
mod handles;
mod hint_generator;
mod line_side_info;
mod triangulation_ext;

pub mod refinement;

pub mod math;

pub use bulk_load::bulk_load;

pub use triangulation_ext::{RemovalResult, TriangulationExt};

pub use dcel::Dcel;
pub use hint_generator::{
    HierarchyHintGenerator, HierarchyHintGeneratorWithBranchFactor, HintGenerator,
    LastUsedVertexHintGenerator,
};

pub use refinement::{AngleLimit, RefinementParameters, RefinementResult};

pub use line_side_info::LineSideInfo;

pub use handles::iterators;
pub use handles::*;

pub use math::{
    validate_coordinate, validate_vertex, InsertionError, MAX_ALLOWED_VALUE, MIN_ALLOWED_VALUE,
};
