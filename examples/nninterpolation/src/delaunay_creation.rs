use rand::distributions::{IndependentSample, Range};
use rand::Rng;

use spade::delaunay::{DelaunayTriangulation, FloatDelaunayTriangulation, TriangulationWalkLookup};
use nalgebra::{Vector2, Vector3};
use spade::{HasPosition,};

use constants::*;

pub type Delaunay = FloatDelaunayTriangulation<
        VectorWithHeight, TriangulationWalkLookup<Vector2<f64>>>;

pub struct VectorWithHeight {
    point: Vector2<f64>,
    pub height: f64,
    pub gradient: Vector2<f64>,
    // We don't need the normal for interpolation purposes. We store it only for
    // visualization.
    pub normal: Vector3<f64>,
}

impl HasPosition for VectorWithHeight {
    type Vector = Vector2<f64>;
    fn position(&self) -> Vector2<f64> {
        self.point
    }
}

impl VectorWithHeight {
    pub fn position_3d(&self) -> Vector3<f64> {
        Vector3::new(self.point.x, self.point.y, self.height)
    }

    pub fn new(point: Vector2<f64>, height: f64) -> VectorWithHeight {
        VectorWithHeight {
            point: point,
            height: height,
            gradient: Vector2::new(0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, 0.0),
        }
    }
}

// Triangulation creation and normal estimation
pub fn generate_random_triangulation() -> Delaunay {

    let mut rng = ::rand::thread_rng();
    let seed = ::noise::Seed::new(rng.gen());
    let mut delaunay = DelaunayTriangulation::with_walk_lookup();

    let range = Range::new(-SAMPLE_REGION, SAMPLE_REGION);
    for _ in 0 .. NUM_POINTS {
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        let height = ::noise::open_simplex2(&seed, &[x * FREQUENCY, y * FREQUENCY]) * MAX_HEIGHT;
        // Try out some other height functions, like those:
        // let height = (x * x + y * y) * 0.3;
        // let height = (x * 3.).sin() + (y - 2.).exp();
        delaunay.insert(VectorWithHeight::new(Vector2::new(x, y), height));
    }

    // Note that, for interpolation, we only need the gradients. For visualization
    // purposes, the normals are also generated and stored within the vertices
    delaunay.estimate_gradients(&(|v| v.height), &(|v, g| v.gradient = g));
    delaunay.estimate_normals(&(|v| v.height), &(|v: &mut VectorWithHeight, n| v.normal = n));
    
    delaunay
}
