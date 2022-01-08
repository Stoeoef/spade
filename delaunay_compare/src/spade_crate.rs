use spade::Triangulation;

type SpadePoint = spade::Point2<f64>;

#[derive(Default)]
pub struct SpadeCrateWithHintGenerator<HintGeneratorType> {
    vertices: Vec<SpadePoint>,
    _hint_generator_type: std::marker::PhantomData<HintGeneratorType>,
}

pub type SpadeCrate = SpadeCrateWithHintGenerator<spade::LastUsedVertexHintGenerator>;
pub type SpadeCrateWithHierarchy = SpadeCrateWithHintGenerator<spade::HierarchyHintGenerator<f64>>;

pub trait HintGeneratorWithMetadata {
    fn is_uniform_insertion_expensive() -> bool;
}

impl HintGeneratorWithMetadata for spade::HierarchyHintGenerator<f64> {
    fn is_uniform_insertion_expensive() -> bool {
        false
    }
}

impl HintGeneratorWithMetadata for spade::LastUsedVertexHintGenerator {
    fn is_uniform_insertion_expensive() -> bool {
        true
    }
}

impl<HintGeneratorType> crate::DelaunayCrate for SpadeCrateWithHintGenerator<HintGeneratorType>
where
    HintGeneratorType: spade::HintGenerator<f64> + HintGeneratorWithMetadata,
{
    type ResultType = spade::DelaunayTriangulation<SpadePoint, (), (), (), HintGeneratorType>;

    fn init(&mut self, vertices: impl IntoIterator<Item = [f64; 2]>) {
        self.vertices = vertices.into_iter().map(|vertex| vertex.into()).collect()
    }

    fn run_creation(&self) -> Self::ResultType {
        Self::ResultType::bulk_load(self.vertices.clone()).unwrap()
    }
}
