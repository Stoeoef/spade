use super::handles::*;

/// Describes a position in a triangulation.
///
/// The position is set in relation to the triangulation's vertices, edges and faces.
/// This type is usually the result of calling [trait.Triangulation.html#method.locate]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum PositionInTriangulation {
    /// A position lies exactly on an existing vertex. The verticis handle is given.
    OnVertex(FixedVertexHandle),

    /// A position lies exactly on an edge. The edge's handle is given.
    OnEdge(FixedDirectedEdgeHandle),

    /// A position lies in the interior of a face. The face's handle is given.
    OnFace(FixedFaceHandle<InnerTag>),

    /// A position lies outside the convex hull. The given edge handle refers to an edge
    /// of the convex hull which has both the point and an outer face on its left side.
    ///
    /// *Note*: The given edge is *not* necessarily the *closest* edge to a position.
    OutsideOfConvexHull(FixedDirectedEdgeHandle),

    /// The triangulation contains either no vertices or exactly one vertex which has a
    /// different position than the query point.
    NoTriangulation,
}

#[cfg(test)]
mod test {
    use std::iter::FromIterator;

    use super::PositionInTriangulation;
    use crate::test_utilities::SEED;
    use crate::test_utilities::*;
    use crate::triangulation::TriangulationExt;
    use crate::DelaunayTriangulation;
    use crate::Triangulation;
    use crate::{handles::FixedVertexHandle, Point2};
    use rand::distributions::{Distribution, Uniform};
    use rand::{seq::SliceRandom, Rng, SeedableRng};
    use rand_hc::Hc128Rng;

    #[test]
    fn test_empty() {
        let d = DelaunayTriangulation::<Point2<f32>>::default();
        assert_eq!(d.num_vertices(), 0);
        assert_eq!(d.num_all_faces(), 1);
        assert_eq!(d.num_undirected_edges(), 0);
    }

    #[test]
    fn test_insert_first() {
        let mut d = DelaunayTriangulation::<Point2<f32>>::default();
        d.insert(Point2::default());
        assert_eq!(d.num_vertices(), 1);
        assert_eq!(d.num_all_faces(), 1);
        assert_eq!(d.num_undirected_edges(), 0);
    }

    #[test]
    fn test_insert_second() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::default());
        d.insert(Point2::new(0.123, 1.234));
        assert_eq!(d.num_vertices(), 2);
        assert_eq!(d.num_all_faces(), 1);
        assert_eq!(d.num_undirected_edges(), 1);
        d.sanity_check();
    }

    #[test]
    fn test_insert_third_point() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(1f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));
        d.insert(Point2::new(1f64, 1f64));

        assert_eq!(d.num_vertices(), 3);
        assert_eq!(d.num_all_faces(), 2);
        d.sanity_check();
    }

    #[test]
    fn test_insert_five_points() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(1f64, 0f64));
        d.insert(Point2::new(0f64, 1f64));

        let v3 = Point2::new(0.433_833_144_214_401f64, 0.900_993_231_373_602_9f64);
        let v4 = Point2::new(2.0, 2.0);
        let v5 = Point2::new(0.5, 0.25);
        d.insert(v3);
        d.sanity_check();
        d.insert(v4);
        d.s().sanity_check();
        d.insert(v5);
        d.sanity_check();
    }

    #[test]
    fn test_small_triangulation_iterators() {
        let mut d = DelaunayTriangulation::<_>::default();
        assert_eq!(d.all_faces().count(), 1);
        assert_eq!(d.inner_faces().count(), 0);

        d.insert(Point2::new(1f64, 1f64));
        assert_eq!(d.all_faces().count(), 1);
        assert_eq!(d.inner_faces().count(), 0);

        d.insert(Point2::new(-1f64, 1f64));
        assert_eq!(d.all_faces().count(), 1);
        assert_eq!(d.inner_faces().count(), 0);
    }

    #[test]
    fn test_iterate_faces() {
        const SIZE: usize = 1000;
        let points = random_points_with_seed(SIZE, SEED);
        let mut d = DelaunayTriangulation::<Point2<f64>>::from_iter(points);
        d.sanity_check();

        assert_eq!(d.all_faces().count(), d.num_all_faces());
        assert_eq!(d.inner_faces().count(), d.num_inner_faces());

        for _ in 0..SIZE / 2 {
            d.remove(FixedVertexHandle::new(5));
        }

        assert_eq!(d.all_faces().count(), d.num_all_faces());
        assert_eq!(d.inner_faces().count(), d.num_inner_faces());

        d.sanity_check();
    }

    #[test]
    fn test_insert_many_points() {
        const SIZE: usize = 10000;
        let points = random_points_with_seed(SIZE, SEED);
        let d = DelaunayTriangulation::<_>::from_iter(points);
        d.sanity_check();
    }

    #[test]
    fn test_insert_outside_convex_hull() {
        const NUM: usize = 100;
        let mut rng = Hc128Rng::from_seed(*SEED);
        let range = Uniform::new(0., 2.0 * ::std::f64::consts::PI);

        let mut d = DelaunayTriangulation::<_>::default();

        // Insert points on a circle. Every new point lies outside the convex hull.
        for _ in 0..NUM {
            let ang = range.sample(&mut rng);
            let vec = Point2::new(ang.sin(), ang.cos());
            d.insert(vec);
        }
        assert_eq!(d.num_vertices(), NUM);
        d.sanity_check();
    }

    #[test]
    fn test_insert_same_point_small() {
        let points = vec![
            Point2::new(0.2, 0.1),
            Point2::new(1.3, 2.2),
            Point2::new(0.0, 0.0),
        ];
        let mut d = DelaunayTriangulation::<_>::from_iter(points.clone());

        for p in &points {
            d.insert(*p);
            d.sanity_check();
        }
        assert_eq!(d.num_vertices(), points.len());
        d.sanity_check();
    }

    #[test]
    fn test_insert_same_point() {
        const SIZE: usize = 300;
        let points = random_points_with_seed(SIZE, SEED);
        let mut d = DelaunayTriangulation::<_>::from_iter(points.clone());
        for p in points {
            d.insert(p);
        }
        assert_eq!(d.num_vertices(), SIZE);
        d.sanity_check();
    }

    #[test]
    fn test_insert_point_on_ch_edge() {
        let points = vec![
            Point2::new(0., 0f64),
            Point2::new(1., 0.),
            Point2::new(0., 1.),
            Point2::new(0., 0.4),
        ];
        let d = DelaunayTriangulation::<_>::from_iter(points);
        d.sanity_check();
    }

    #[test]
    fn test_insert_on_edges() {
        let points = vec![Point2::new(0., 0f64), Point2::new(1., 0.)];
        let mut d = DelaunayTriangulation::<_>::from_iter(points);

        d.insert(Point2::new(1., 1.));
        d.sanity_check();
        d.insert(Point2::new(0.5, 0.5));
        d.sanity_check();
        d.insert(Point2::new(0., 0.4));
        d.sanity_check();
        d.insert(Point2::new(1., 0.5));
        d.sanity_check();
        d.insert(Point2::new(0.5, 1.));
        d.sanity_check();
        d.insert(Point2::new(0.7, 0.));
        d.sanity_check();
    }

    #[test]
    fn test_degenerate_triangulation() {
        let mut d = DelaunayTriangulation::<_>::default();
        for i in -50..50 {
            d.insert(Point2::new(f64::from(i), 0.));
        }

        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(0.0, 1.0));
        for i in -50..50 {
            d.insert(Point2::new(f64::from(i), 0.));
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_line_2() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(0.0, 1.0));

        for i in -50..50 {
            d.insert(Point2::new(f64::from(i), 0.));
            d.sanity_check();
        }

        for i in -10..10 {
            d.insert(Point2::new(f64::from(i), 0.5 * f64::from(i)));
            d.sanity_check();
        }
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_grid2() {
        let mut d = DelaunayTriangulation::<_>::default();

        for y in 0..20 {
            for x in 0..7 {
                d.insert(Point2::new(f64::from(x), f64::from(y)));
                d.sanity_check();
            }
        }
        d.sanity_check();
    }

    /*
        TODO
        #[derive(Default)]
        struct PointWithHeight {
            point: Point2<f64>,
        }

        impl HasPosition for PointWithHeight {
            type Scalar = f64;

            fn position(&self) -> Point2<f64> {
                self.point
            }
        }

        impl PointWithHeight {
            fn new(x: f64, y: f64, height: f64) -> PointWithHeight {
                PointWithHeight {
                    point: Point2::new(x, y),
                }
            }
        }
    */
    /*
    #[test]
    fn test_natural_neighbor_interpolation() {
        let points = vec![
            PointWithHeight::new(0.0, 0.0, 0.0),
            PointWithHeight::new(1.0, 0.0, 0.0),
            PointWithHeight::new(0.0, 1.0, 0.0),
            PointWithHeight::new(1.0, 1.0, 0.0),
            PointWithHeight::new(2.0, 0.0, 0.0),
            PointWithHeight::new(2.0, 1.0, 0.0),
            PointWithHeight::new(3.0, 0.0, 1.0),
            PointWithHeight::new(3.0, 1.0, 1.0),
            PointWithHeight::new(4.0, 0.0, 1.0),
            PointWithHeight::new(4.0, 1.0, 1.0),
        ];
        let d = DelaunayTriangulation::<_>::from_iter(Default::default(), points).unwrap();
        assert_eq!(
            d.nn_interpolation(Point::new(0.5, 0.5), |p| p.height),
            Some(0.0)
        );
        assert_eq!(
            d.nn_interpolation(Point::new(0.2, 0.8), |p| p.height),
            Some(0.0)
        );
        assert_eq!(
            d.nn_interpolation(Point::new(3.5, 1.), |p| p.height),
            Some(1.0)
        );
        assert_eq!(
            d.nn_interpolation(Point::new(-20., 0.2), |p| p.height),
            Some(0.0)
        );
        let height = d
            .nn_interpolation(Point::new(3.2, 0.9), |p| p.height)
            .unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        let height = d
            .nn_interpolation(Point::new(3.5, 0.5), |p| p.height)
            .unwrap();
        assert!((height - 1.0).abs() < 0.00001);
        assert_eq!(
            d.nn_interpolation(Point::new(3.0, 0.0), |p| p.height),
            Some(1.0)
        );
    }
    */

    #[test]
    fn test_insert_points_with_increasing_distance() {
        let mut points = random_points_with_seed(1000, SEED);
        points.sort_by(|p1, p2| p1.length2().partial_cmp(&p2.length2()).unwrap());
        let d = DelaunayTriangulation::<_>::from_iter(points);
        d.sanity_check();
    }

    #[test]
    fn test_insert_points_on_grid_with_increasing_distance() {
        // This test inserts points on a grid with increasing distance from (0., 0.)
        let mut points = Vec::new();
        const SIZE: i64 = 7;
        for x in -SIZE..SIZE {
            for y in -SIZE..SIZE {
                let point = Point2::new(x as f64, y as f64);
                points.push(point);
            }
        }
        points.sort_by(|p1, p2| p1.length2().partial_cmp(&p2.length2()).unwrap());
        let d = DelaunayTriangulation::<_>::from_iter(points);
        d.sanity_check();
    }

    #[test]
    fn test_remove_in_triangle() {
        let points = vec![
            Point2::new(-1.0, 0.0f64),
            Point2::new(1.0, 0.0f64),
            Point2::new(0.0, 1.0f64),
        ];
        let mut d = DelaunayTriangulation::<_>::from_iter(points);
        let to_remove = d.insert(Point2::new(0.0, 0.5));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 3);
        // Reinsert the last point, just to see if a crash occurs
        d.insert(Point2::new(0.0, 0.5));
        d.sanity_check();
    }

    #[test]
    fn test_remove_complex_single_outer_vertex() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(0.0, 0.0));
        d.insert(Point2::new(0.0, 1.0));
        d.insert(Point2::new(4.0, 0.0));
        d.insert(Point2::new(4.0, 1.0));
        d.insert(Point2::new(2.0, 0.5));
        d.insert(Point2::new(1.0, 0.5));
        d.insert(Point2::new(3.0, 0.5));

        let v4_position = Point2::new(2.5, 2.0);
        let v4 = d.insert(v4_position);

        let removed = d.remove(v4);
        d.sanity_check();
        assert_eq!(removed, v4_position);
        assert_eq!(d.num_vertices(), 7);
    }

    #[test]
    fn test_remove_single_outer() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(0.0, 0.0));
        d.insert(Point2::new(0.0, 1.0));
        d.insert(Point2::new(1.0, 0.0));
        let v4_position = Point2::new(1.5, 1.5);
        let v4 = d.insert(v4_position);

        let removed = d.remove(v4);
        d.sanity_check();
        assert_eq!(removed, v4_position);
        assert_eq!(d.num_vertices(), 3);
    }

    #[test]
    fn test_remove_in_quad() {
        let points = vec![
            Point2::new(0.0, 0.0f64),
            Point2::new(1.0, 0.0f64),
            Point2::new(0.0, 1.0f64),
            Point2::new(1.0, 1.0f64),
        ];

        let mut d = DelaunayTriangulation::<_>::from_iter(points);

        let to_remove = d.insert(Point2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        let to_remove = d.insert(Point2::new(0.5, 0.6));
        d.remove(to_remove);
        assert_eq!(d.num_vertices(), 4);
        d.insert(Point2::new(0.5, 0.6));
        d.sanity_check();
    }

    #[test]
    fn test_remove_star_shaped() {
        let mut rng = Hc128Rng::from_seed(*SEED);
        let mut points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, -1.0),
            Point2::new(-1.0, 1.0),
            Point2::new(-1.0, -1.0),
            Point2::new(0.0, 3.0),
            Point2::new(0.0, -3.0),
            Point2::new(3.0, 0.0),
            Point2::new(-3.0, 0.0),
        ];

        points.shuffle(&mut rng);
        points.shuffle(&mut rng);
        points.shuffle(&mut rng);
        //for i in 0..20
        {
            //points.shuffle(&mut rng);
            let mut d = DelaunayTriangulation::<_>::from_iter(points.iter().copied());
            d.locate_and_remove(Point2::new(0.0, 0.0));
            d.sanity_check();
        }
    }

    #[test]
    fn test_remove_inner() {
        use ::rand::SeedableRng;
        use PositionInTriangulation::OnVertex;

        let mut points = random_points_with_seed(1000, SEED);
        let mut d = DelaunayTriangulation::<_>::from_iter(points.clone());

        // Insert an outer quad since we don't want to remove vertices from
        // the convex hull.
        d.insert(Point2::new(-2.0, -2.0));
        d.insert(Point2::new(-2.0, 2.0));
        d.insert(Point2::new(2.0, -2.0));
        d.insert(Point2::new(2.0, 2.0));
        // Now remove all inner points
        let mut rng = Hc128Rng::from_seed(*SEED);
        points.shuffle(&mut rng);
        assert_eq!(d.num_vertices(), 1004);
        for point in points {
            match d.locate(point) {
                OnVertex(handle) => {
                    d.remove(handle);
                }
                _ => panic!("Point lookup failed: {:?}", point),
            }
        }
        assert_eq!(d.num_vertices(), 4);
        d.sanity_check();
    }

    #[test]
    fn test_remove_outer() {
        use PositionInTriangulation::OnVertex;

        let mut points = random_points_with_seed(1000, SEED);
        let mut d = DelaunayTriangulation::<_>::from_iter(points.clone());

        points.sort_by(|p1, p2| p1.length2().partial_cmp(&p2.length2()).unwrap());
        for (index, point) in points[3..].iter().rev().enumerate() {
            match d.locate(*point) {
                OnVertex(handle) => {
                    d.remove(handle);
                    if index % 50 == 0 {
                        // Check only every 50 iterations to reduce test runtime
                        d.sanity_check();
                    }
                }
                _ => panic!("Point lookup failed: {:?}", point),
            }
        }
        d.sanity_check();
    }

    #[test]
    fn test_removal_and_insertion() {
        let points = random_points_with_seed(1000, SEED);
        let mut d = DelaunayTriangulation::<_>::from_iter(points);

        let mut rng = Hc128Rng::from_seed(*SEED);
        for _ in 0..1000 {
            if rng.gen() {
                // Insert new random point
                let x = rng.gen();
                let y = rng.gen();
                d.insert(Point2::new(x, y));
            } else {
                // Remove random point
                let range = Uniform::new(1, d.num_vertices());
                let handle = range.sample(&mut rng);
                d.remove(FixedVertexHandle::new(handle));
            }
        }
        d.sanity_check();
    }

    #[test]
    fn test_remove_until_empty() {
        let mut d = DelaunayTriangulation::<_>::from_iter(vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 2.0),
        ]);

        while let Some(to_remove) = d.vertices().next() {
            let to_remove = to_remove.fix();
            d.remove(to_remove);
            d.sanity_check();
        }

        assert_eq!(d.num_vertices(), 0);

        d.insert(Point2::new(1.0, 0.0));
        d.insert(Point2::new(1.0, 1.0));
        d.insert(Point2::new(1.0, 2.0));
        d.insert(Point2::new(2.3, 1.4));
        d.sanity_check();
    }

    #[test]
    fn test_remove_until_degenerate() {
        let points = vec![
            Point2::new(0., 0f64),
            Point2::new(1., 0.),
            Point2::new(0., 1.),
            Point2::new(0., 0.5),
            Point2::new(0., 0.25),
            Point2::new(0., 0.75),
        ];
        let mut d = DelaunayTriangulation::<_>::from_iter(points);
        assert_eq!(d.num_all_faces(), 5);
        d.remove(FixedVertexHandle::new(1));
        d.sanity_check();
        assert!(d.all_vertices_on_line());

        while let Some(to_remove) = d.vertices().next() {
            let to_remove = to_remove.fix();
            d.remove(to_remove);
            d.sanity_check();
        }

        d.sanity_check();
        d.insert(Point2::new(0.5, 0.5));
        d.insert(Point2::new(0.2, 0.5));
        d.insert(Point2::new(1.5, 0.0));
        d.sanity_check();
    }

    #[test]
    fn test_remove_few_points() {
        let mut triangulation = DelaunayTriangulation::<_>::from_iter(vec![
            Point2::new(0.0, 1.0),
            Point2::new(100.0, 1.0),
            Point2::new(0.0, 110.0),
            Point2::new(110.110, 110.0),
            Point2::new(50.0, 50.0),
            Point2::new(50.0, 80.0),
            Point2::new(75.0, 80.0),
        ]);

        triangulation.remove(FixedVertexHandle::new(5));
        triangulation.sanity_check();
        triangulation.remove(FixedVertexHandle::new(5));
        triangulation.sanity_check();
    }

    #[test]
    fn test_remove_on_line_small() {
        let mut triangulation = DelaunayTriangulation::<_>::from_iter(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0), // This point will be removed
            Point2::new(2.0, 0.0),
        ]);
        triangulation.remove(FixedVertexHandle::new(2));
        triangulation.sanity_check();
    }

    #[test]
    fn test_remove_on_line_big() {
        let mut triangulation = DelaunayTriangulation::<_>::default();
        for x in 2..100 {
            triangulation.insert(Point2::new(f64::from(x), 0.0));
        }
        let mut rng = Hc128Rng::from_seed(*SEED);
        while triangulation.num_vertices() > 3 {
            if rng.gen() {
                triangulation.remove(FixedVertexHandle::new(1));
            } else {
                triangulation.remove(FixedVertexHandle::new(2));
            }

            triangulation.sanity_check();
        }
    }

    #[test]
    fn test_small_insert_on_line() {
        let mut d = DelaunayTriangulation::<_>::default();
        d.insert(Point2::new(0.0, 0.0));
        d.insert(Point2::new(2.0, 0.0));
        d.insert(Point2::new(1.0, 0.0));
        d.sanity_check();
    }

    #[test]
    fn test_locate_when_empty() {
        let triangulation = DelaunayTriangulation::<Point2<_>>::default();
        assert_eq!(
            triangulation.locate(Point2::new(0.0, 0.0)),
            PositionInTriangulation::NoTriangulation
        )
    }

    #[test]
    fn test_locate_with_single_vertex() {
        let mut triangulation = DelaunayTriangulation::<_>::default();
        let point = Point2::new(0.0, 0.0);
        triangulation.insert(point);
        assert_eq!(
            triangulation.locate(point),
            PositionInTriangulation::OnVertex(FixedVertexHandle::new(0))
        );
        assert_eq!(
            triangulation.locate(Point2::new(1.0, 1.0)),
            PositionInTriangulation::NoTriangulation
        )
    }

    #[test]
    fn test_remove_on_line_end() {
        let mut triangulation = DelaunayTriangulation::<_>::from_iter(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
        ]);
        triangulation.remove(FixedVertexHandle::new(2));
        triangulation.sanity_check();
    }
}
