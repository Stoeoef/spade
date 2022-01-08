#[derive(Default)]
pub struct CdtCrate {
    vertices: Vec<(f64, f64)>,
}

impl crate::DelaunayCrate for CdtCrate {
    type ResultType = Vec<(usize, usize, usize)>;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>) {
        self.vertices = vertices.map(|array| (array[0], array[1])).collect()
    }

    fn run_creation(&self) -> Self::ResultType {
        cdt::triangulate_points(&self.vertices).unwrap()
    }
}
