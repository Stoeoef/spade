use std::time::Instant;

use anyhow::{Context, Ok};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use tiny_skia::{Color, Paint, PathBuilder, Pixmap, Stroke, Transform};

fn main() -> anyhow::Result<()> {
    let shape_file_path = "./examples/Europe_coastline.shp";
    println!("Loading {} ...", shape_file_path);
    let mut reader = shapefile::Reader::from_path(shape_file_path)?;
    println!("Loaded!");

    println!("Extracting data...");
    let mut vertices = Vec::new();
    let mut edges = Vec::new();

    for shape_record in reader.iter_shapes_and_records() {
        let (shape, _) = shape_record?;

        match shape {
            shapefile::Shape::Polyline(line) => {
                for part in line.parts() {
                    let first_vertex = vertices.len();
                    vertices.extend(part.iter().map(|p| Point2::new(p.x, p.y)));
                    let last_vertex = vertices.len() - 1;
                    edges.extend((first_vertex..last_vertex).map(|i| [i, i + 1]));
                }
            }
            _ => unimplemented!(),
        }
    }

    println!("Done!");
    println!("{} vertices", vertices.len());
    println!("{} constraint edges", edges.len());

    vertices.shrink_to_fit();
    edges.shrink_to_fit();

    load_with_spade(&vertices, &edges)?;
    println!();

    load_with_cdt_crate(&vertices, &edges)?;

    Ok(())
}

fn load_with_spade(vertices: &Vec<Point2<f64>>, edges: &Vec<[usize; 2]>) -> anyhow::Result<()> {
    let vertices_clone = vertices.clone();
    let edges_clone = edges.clone();

    println!("Loading triangulation (spade)...");
    let now = Instant::now();
    let cdt =
        spade::ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt(vertices_clone, edges_clone)?;
    println!("Done!");
    println!("{} vertices (without duplicates)", cdt.num_vertices());
    println!("{} undirected edges", cdt.num_undirected_edges());
    println!("{} constraint edges", cdt.num_constraints());
    println!("{} triangles", cdt.num_inner_faces());
    println!("{} convex hull edges", cdt.convex_hull_size());
    println!();
    println!(
        "loading time (spade cdt bulk load): {}ms",
        now.elapsed().as_millis()
    );

    let vertices_clone = vertices.clone();
    let edges_clone = edges.clone();

    let now = Instant::now();
    spade::ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt_stable(
        vertices_clone,
        edges_clone,
    )?;

    println!(
        "loading time (spade cdt bulk load stable): {}ms",
        now.elapsed().as_millis()
    );

    let vertices_clone = vertices.clone();

    let now = Instant::now();
    spade::ConstrainedDelaunayTriangulation::<_>::bulk_load(vertices_clone)?;
    println!(
        "loading time (spade without constraints): {}ms",
        now.elapsed().as_millis()
    );

    println!();
    // Let's turn the triangulation into an image, just to verify that everything works!
    println!("Creating and saving output image...");
    draw_to_pixmap(cdt)?.save_png("examples/europe.png")?;
    println!("Done!");

    Ok(())
}

fn load_with_cdt_crate(vertices: &[Point2<f64>], edges: &[[usize; 2]]) -> anyhow::Result<()> {
    let vertices_clone = vertices.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>();

    println!("Loading cdt (cdt crate)");
    let edges = edges
        .iter()
        .map(|[from, to]| (*from, *to))
        .collect::<Vec<_>>();
    let now = Instant::now();
    cdt::triangulate_with_edges(&vertices_clone, &edges)?;
    println!("Done!");
    println!(
        "loading time (cdt crate with constraints): {}ms",
        now.elapsed().as_millis()
    );

    let now = Instant::now();
    cdt::triangulate_points(&vertices_clone)?;
    println!("Done!");
    println!(
        "loading time (cdt crate without constraints): {}ms",
        now.elapsed().as_millis()
    );

    Ok(())
}

fn draw_to_pixmap(cdt: ConstrainedDelaunayTriangulation<Point2<f64>>) -> anyhow::Result<Pixmap> {
    let mut min_x = f64::MAX;
    let mut min_y = f64::MAX;
    let mut max_x = f64::MIN;
    let mut max_y = f64::MIN;

    for position in cdt.convex_hull().map(|edge| edge.from().position()) {
        if position.x < min_x {
            min_x = position.x;
        }
        if position.y < min_y {
            min_y = position.y;
        }
        if position.x > max_x {
            max_x = position.x;
        }
        if position.y > max_y {
            max_y = position.y;
        }
    }

    let mut constraints = PathBuilder::new();
    let mut edges = PathBuilder::new();
    for ([from, to], is_constraint_edge) in cdt
        .undirected_edges()
        .map(|edge| (edge.positions(), edge.is_constraint_edge()))
    {
        if is_constraint_edge {
            constraints.move_to(from.x as f32, from.y as f32);
            constraints.line_to(to.x as f32, to.y as f32);
        } else {
            edges.move_to(from.x as f32, from.y as f32);
            edges.line_to(to.x as f32, to.y as f32);
        }
    }
    let constraints = constraints
        .finish()
        .context("Failed to finish constraint path")?;
    let edges = edges.finish().context("Failed to finish edge path")?;

    let res_x = 1024.0;
    let res_y = 1024.0;

    let mut pixmap = Pixmap::new(res_x as u32, res_y as u32).unwrap();

    let scale_x = res_x / (max_x - min_x);
    let scale_y = res_y / (max_y - min_y);

    let transform = Transform::from_translate(-min_x as f32, -min_y as f32)
        .post_scale(scale_x as f32, -scale_y as f32)
        .post_translate(0.0, res_y as f32);

    let mut stroke = Stroke::default();
    stroke.width = 0.1 / scale_x as f32;

    let mut paint = Paint::default();

    paint.set_color_rgba8(50, 127, 150, 255);

    pixmap.stroke_path(&edges, &paint, &stroke, transform, None);
    paint.set_color_rgba8(200, 126, 150, 255);
    stroke.width = 2.0 / scale_x as f32;
    pixmap.stroke_path(&constraints, &paint, &stroke, transform, None);

    Ok(pixmap)
}
