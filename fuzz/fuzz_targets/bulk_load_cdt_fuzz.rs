#![no_main]
mod fuzz_shared;
use fuzz_shared::FuzzPoint;
use libfuzzer_sys::fuzz_target;

use spade::{ConstrainedDelaunayTriangulation, HasPosition, Triangulation};

fuzz_target!(|input: (Vec<FuzzPoint>, Vec<[usize; 2]>)| {
    let (data, mut edges) = input;
    for p in &data {
        if spade::validate_coordinate(p.x).is_err() || spade::validate_coordinate(p.y).is_err() {
            return;
        }
        if p.x.abs() > 20.0 || p.y.abs() > 20.0 {
            return;
        }
    }

    for &[from, to] in &edges {
        if from >= data.len() || to >= data.len() || from == to {
            return;
        }
    }

    let mut reference_cdt =
        ConstrainedDelaunayTriangulation::<FuzzPoint>::bulk_load(data.clone()).unwrap();

    let mut last_index = 0;
    for (index, [from, to]) in edges.iter().copied().enumerate() {
        let from = reference_cdt
            .locate_vertex(data[from].position())
            .unwrap()
            .fix();
        let to = reference_cdt
            .locate_vertex(data[to].position())
            .unwrap()
            .fix();

        if reference_cdt.can_add_constraint(from, to) {
            reference_cdt.add_constraint(from, to);
        } else {
            last_index = index;
            break;
        }
    }

    edges.truncate(last_index);

    let bulk_loaded =
        ConstrainedDelaunayTriangulation::<FuzzPoint>::bulk_load_cdt(data, edges).unwrap();

    bulk_loaded.cdt_sanity_check();
});
