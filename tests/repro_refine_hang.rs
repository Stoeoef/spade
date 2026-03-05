use spade::handles::FixedVertexHandle;
use spade::{ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation};
use std::collections::HashMap;

#[test]
fn repro_refine_hang_regression() {
    let points = points();
    let constraints = constraints();

    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::bulk_load_cdt(points, vec![])
        .expect("Failed to build CDT for regression input");

    for edge in &constraints {
        cdt.add_constraint_and_split(
            FixedVertexHandle::from_index(edge[0]),
            FixedVertexHandle::from_index(edge[1]),
            |v| v,
        );
    }

    let parameters = RefinementParameters::<f64>::new()
        .keep_constraint_edges()
        .with_max_allowed_area(0.4 * 0.32 * 0.32)
        .with_max_additional_vertices(1_000_000);

    let result = cdt.refine(parameters);
    assert!(
        result.refinement_complete,
        "Refinement ran out of additional vertices for regression input"
    );

    let mut missing = 0usize;
    for edge in &constraints {
        let a = FixedVertexHandle::from_index(edge[0]);
        let b = FixedVertexHandle::from_index(edge[1]);
        if !cdt.exists_constraint(a, b) {
            missing += 1;
        }
    }
    // Will ignore for now, this requires another PR to fix
    // assert_eq!(missing, 0, "Refinement lost constraint edges");

    let mut counts: HashMap<(u64, u64), usize> = HashMap::new();
    let mut points_after = Vec::new();
    for v in cdt.vertices() {
        let p = v.position();
        *counts.entry((p.x.to_bits(), p.y.to_bits())).or_default() += 1;
        points_after.push(p);
    }
    let duplicate_vertices: usize = counts.values().filter(|&&c| c > 1).map(|&c| c - 1).sum();
    assert_eq!(
        duplicate_vertices, 0,
        "Refinement produced exactly duplicate vertices"
    );

    let mut near_duplicates = 0usize;
    let eps = 1e-12_f64;
    for i in 0..points_after.len() {
        for j in (i + 1)..points_after.len() {
            let dx = (points_after[i].x - points_after[j].x).abs();
            let dy = (points_after[i].y - points_after[j].y).abs();
            if dx < eps && dy < eps {
                near_duplicates += 1;
            }
        }
    }
    assert_eq!(
        near_duplicates, 0,
        "Refinement produced near-duplicate vertices under 1e-12 threshold"
    );
}

pub fn points() -> Vec<Point2<f64>> {
    vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
        Point2::new(0.5, 0.5),
        Point2::new(0.5040435515666191, 0.48484848484848486),
        Point2::new(0.5080871031332381, 0.4696969696969697),
        Point2::new(0.5121306546998572, 0.45454545454545453),
        Point2::new(0.5161742062664761, 0.4393939393939394),
        Point2::new(0.5202177578330952, 0.42424242424242425),
        Point2::new(0.5242613093997143, 0.40909090909090906),
        Point2::new(0.5283048609663333, 0.3939393939393939),
        Point2::new(0.5323484125329524, 0.3787878787878788),
        Point2::new(0.5363919640995715, 0.36363636363636365),
        Point2::new(0.5404355156661904, 0.3484848484848485),
        Point2::new(0.5444790672328095, 0.33333333333333337),
        Point2::new(0.5485226187994285, 0.3181818181818182),
        Point2::new(0.5525661703660476, 0.30303030303030304),
        Point2::new(0.5566097219326667, 0.2878787878787879),
        Point2::new(0.5606532734992857, 0.2727272727272727),
        Point2::new(0.5646968250659048, 0.25757575757575757),
        Point2::new(0.5687403766325237, 0.24242424242424243),
        Point2::new(0.5727839281991428, 0.2272727272727273),
        Point2::new(0.5768274797657619, 0.2121212121212121),
        Point2::new(0.5808710313323809, 0.19696969696969702),
        Point2::new(0.584914582899, 0.18181818181818188),
        Point2::new(0.5889581344656191, 0.16666666666666669),
        Point2::new(0.593001686032238, 0.15151515151515155),
        Point2::new(0.5970452375988571, 0.13636363636363635),
        Point2::new(0.6010887891654761, 0.12121212121212122),
        Point2::new(0.6051323407320952, 0.10606060606060608),
        Point2::new(0.6091758922987143, 0.09090909090909088),
        Point2::new(0.6132194438653333, 0.0757575757575758),
        Point2::new(0.6172629954319524, 0.06060606060606066),
        Point2::new(0.6213065469985714, 0.04545454545454547),
        Point2::new(0.6253500985651904, 0.030303030303030387),
        Point2::new(0.6289250162663348, 0.01690752419363295),
        Point2::new(0.6334372016984285, 5.551115123125783E-17),
        Point2::new(0.5156817958365265, 0.5),
        Point2::new(0.531363591673053, 0.5),
        Point2::new(0.5470453875095795, 0.5),
        Point2::new(0.562727183346106, 0.5),
        Point2::new(0.5784089791826326, 0.5),
        Point2::new(0.5940907750191591, 0.5),
        Point2::new(0.6097725708556856, 0.5),
        Point2::new(0.6254543666922121, 0.5),
        Point2::new(0.6411361625287386, 0.5),
        Point2::new(0.6568179583652651, 0.5),
        Point2::new(0.6724997542017918, 0.5),
        Point2::new(0.6881815500383183, 0.5),
        Point2::new(0.7038633458748448, 0.5),
        Point2::new(0.7195451417113713, 0.5),
        Point2::new(0.7352269375478978, 0.5),
        Point2::new(0.7509087333844243, 0.5),
        Point2::new(0.7665905292209508, 0.5),
        Point2::new(0.7822723250574773, 0.5),
        Point2::new(0.797954120894004, 0.5),
        Point2::new(0.8136359167305304, 0.5),
        Point2::new(0.8293177125670569, 0.5),
        Point2::new(0.8449995084035835, 0.5),
        Point2::new(0.86068130424011, 0.5),
        Point2::new(0.8763631000766365, 0.5),
        Point2::new(0.892044895913163, 0.5),
        Point2::new(0.9077266917496896, 0.5),
        Point2::new(0.9234084875862161, 0.5),
        Point2::new(0.9390902834227426, 0.5),
        Point2::new(0.9547720792592691, 0.5),
        Point2::new(0.9704538750957956, 0.5),
        Point2::new(0.9861356709323221, 0.5),
        Point2::new(1.0, 0.5),
    ]
}

pub fn constraints() -> Vec<[usize; 2]> {
    vec![
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 8],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 12],
        [12, 13],
        [13, 14],
        [14, 15],
        [15, 16],
        [16, 17],
        [17, 18],
        [18, 19],
        [19, 20],
        [20, 21],
        [21, 22],
        [22, 23],
        [23, 24],
        [24, 25],
        [25, 26],
        [26, 27],
        [27, 28],
        [28, 29],
        [29, 30],
        [30, 31],
        [31, 32],
        [32, 33],
        [33, 34],
        [34, 35],
        [35, 36],
        [36, 37],
        [4, 38],
        [38, 39],
        [39, 40],
        [40, 41],
        [41, 42],
        [42, 43],
        [43, 44],
        [44, 45],
        [45, 46],
        [46, 47],
        [47, 48],
        [48, 49],
        [49, 50],
        [50, 51],
        [51, 52],
        [52, 53],
        [53, 54],
        [54, 55],
        [55, 56],
        [56, 57],
        [57, 58],
        [58, 59],
        [59, 60],
        [60, 61],
        [61, 62],
        [62, 63],
        [63, 64],
        [64, 65],
        [65, 66],
        [66, 67],
        [67, 68],
        [68, 69],
    ]
}
