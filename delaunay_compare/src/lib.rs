pub mod cdt_crate;
pub mod delaunator_crate;
pub mod spade_crate;

/// Abstraction over different crates with Delaunay triangulations
pub trait DelaunayCrate: Default {
    type ResultType;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>);
    fn run_creation(&self) -> Self::ResultType;
}
