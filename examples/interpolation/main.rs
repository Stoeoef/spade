use anyhow::*;
use image::{ImageBuffer, Rgba};
use std::io::{Cursor, Write};

use base64::Engine;
use tiny_skia::*;

use spade::{DelaunayTriangulation, FloatTriangulation, HasPosition, Point2, Triangulation};

#[derive(Debug, Copy, Clone)]
struct PointWithHeight {
    position: Point2<f64>,
    height: f64,
}

impl PointWithHeight {
    const fn new(x: f64, y: f64, height: f64) -> Self {
        Self {
            position: Point2::new(x, y),
            height,
        }
    }
}

impl HasPosition for PointWithHeight {
    type Scalar = f64;

    fn position(&self) -> Point2<Self::Scalar> {
        self.position
    }
}

type TriangulationType = DelaunayTriangulation<PointWithHeight>;

const VERTICES: &[PointWithHeight] = &[
    PointWithHeight::new(20.0, 20.0, 0.9),
    PointWithHeight::new(20.0, -20.0, 0.9),
    PointWithHeight::new(-20.0, 20.0, 0.9),
    PointWithHeight::new(-20.0, -20.0, 0.9),
    PointWithHeight::new(20.0, 10.0, 0.0),
    PointWithHeight::new(10.0, 10.0, 0.5),
    PointWithHeight::new(0.0, 23.0, 0.1),
    PointWithHeight::new(0.0, -23.0, 0.1),
    PointWithHeight::new(-10.0, 10.0, 0.0),
    PointWithHeight::new(-10.0, -10.0, 0.1),
    PointWithHeight::new(-15.0, 20.0, 0.5),
    PointWithHeight::new(-20.0, 20.0, 0.0),
    PointWithHeight::new(-5.0, 0.0, 0.25),
    PointWithHeight::new(5.0, 7.0, 0.75),
    PointWithHeight::new(12.0, -10.0, 0.4),
    PointWithHeight::new(5.0, 10.0, 0.3),
    PointWithHeight::new(5.0, 1.0, 0.2),
];

pub fn main() -> Result<()> {
    let mut t = TriangulationType::default();
    for vertex in VERTICES {
        t.insert(*vertex)?;
    }

    const OFFSET: f32 = 25.0;
    const SCALAR: f32 = 50.0 / 512.0;

    fn px_to_coords(x: u32, y: u32) -> Point2<f64> {
        Point2::new(
            x as f64 * SCALAR as f64 - OFFSET as f64,
            y as f64 * SCALAR as f64 - OFFSET as f64,
        )
    }

    let dimensions = 512;
    let mut nn_c0_pixmap =
        Pixmap::new(dimensions, dimensions).context("Failed to allocate image")?;
    let mut nn_c1_pixmap =
        Pixmap::new(dimensions, dimensions).context("Failed to allocate image")?;
    let mut barycentric_pixmap =
        Pixmap::new(dimensions, dimensions).context("Failed to allocate image")?;
    let mut nearest_neighbor_pixmap =
        Pixmap::new(dimensions, dimensions).context("Failed to allocate image")?;

    fn set_pixel(pixmap: &mut Pixmap, x: u32, y: u32, value: Option<f64>) {
        let background_color = [255, 255, 255];
        let [r, g, b] = value.map(float_to_color).unwrap_or(background_color);
        let base = (y * pixmap.height() + x) as usize * 4;
        pixmap.data_mut()[base] = r;
        pixmap.data_mut()[base + 1] = g;
        pixmap.data_mut()[base + 2] = b;
        // Alpha
        pixmap.data_mut()[base + 3] = 255;
    }

    let nn = t.natural_neighbor();
    let barycentric = t.barycentric();

    for y in 0..dimensions {
        for x in 0..dimensions {
            let coords = px_to_coords(x, y);
            let value_nn_c0 = nn.interpolate(|v| v.data().height, coords);
            set_pixel(&mut nn_c0_pixmap, x, y, value_nn_c0);

            let value_nn_c1 =
                nn.interpolate_gradient(|v| v.data().height, |_| [0.0, 0.0], 1.0, coords);
            set_pixel(&mut nn_c1_pixmap, x, y, value_nn_c1);

            let value_barycentric = barycentric.interpolate(|v| v.data().height, coords);
            set_pixel(&mut barycentric_pixmap, x, y, value_barycentric);

            let value_nearest_neighbor = t.nearest_neighbor(coords).map(|v| v.data().height);
            set_pixel(&mut nearest_neighbor_pixmap, x, y, value_nearest_neighbor);
        }
    }

    let path = {
        let mut pb = PathBuilder::new();

        for edge in t.undirected_edges() {
            let [from, to] = edge.positions();
            pb.move_to(from.x as f32, from.y as f32);
            pb.line_to(to.x as f32, to.y as f32);
        }

        pb.finish().context("Could not build path")?
    };

    let stroke = Stroke {
        width: 0.15,
        ..Default::default()
    };

    let mut paint = Paint::default();
    paint.set_color_rgba8(0, 0, 0, 70);
    paint.anti_alias = true;

    for pixmap in [
        &mut nn_c0_pixmap,
        &mut nn_c1_pixmap,
        &mut barycentric_pixmap,
        &mut nearest_neighbor_pixmap,
    ] {
        pixmap.stroke_path(
            &path,
            &paint,
            &stroke,
            Transform::from_translate(OFFSET, OFFSET).post_scale(1.0 / SCALAR, 1.0 / SCALAR),
            None,
        );
    }

    fn save_pixmap(pixmap: Pixmap, name: &str) -> Result<()> {
        // tiny_skia doesn't support jpg encoding which is required for small file size when embedding this into
        // the documentation. We'll have to convert the data into ImageBuffer from the image crate and then do
        // the jpeg encoding.
        let (width, height) = (pixmap.width(), pixmap.height());
        let data = pixmap.take();
        let buffer = ImageBuffer::<Rgba<u8>, _>::from_vec(width, height, data)
            .context("Failed to convert to ImageBuffer")?;

        let mut data_jpeg: Cursor<Vec<u8>> = Cursor::new(Vec::new());
        buffer.write_to(&mut data_jpeg, image::ImageFormat::Jpeg)?;

        std::fs::write(format!("images/{}.jpg", name), data_jpeg.get_ref())?;

        // Encode image as <img> tag for inclusion into the documentation
        let encoded = base64::engine::general_purpose::STANDARD_NO_PAD.encode(data_jpeg.get_ref());
        let mut file = std::fs::File::create(format!("images/{}.img", name))?;
        write!(file, r#"<img src="data:image/jpg;base64,{}" />"#, encoded)?;
        Ok(())
    }

    save_pixmap(nn_c0_pixmap, "interpolation_nn_c0")?;
    save_pixmap(nn_c1_pixmap, "interpolation_nn_c1")?;
    save_pixmap(barycentric_pixmap, "interpolation_barycentric")?;
    save_pixmap(nearest_neighbor_pixmap, "interpolation_nearest_neighbor")?;

    Ok(())
}

fn float_to_color(value: f64) -> [u8; 3] {
    // mostly AI generated...
    // Converts a hue value in the range 0.0 ..= 1.0 from HLS to RGB
    let value = value.clamp(0.0, 1.0);

    const LIGHTNESS: f64 = 0.45;
    const SATURATION: f64 = 0.55;

    let c = (1.0 - (2.0 * LIGHTNESS - 1.0).abs()) * SATURATION;

    let hue = value * 360.0;
    let hs = hue / 60.0;

    let x = c * (1.0 - (hs % 2.0 - 1.0).abs());
    let m = LIGHTNESS - c * 0.5;

    let (r, g, b) = if (0.0..60.0).contains(&hue) {
        (c, x, 0.0)
    } else if (60.0..120.0).contains(&hue) {
        (x, c, 0.0)
    } else if (120.0..180.0).contains(&hue) {
        (0.0, c, x)
    } else if (180.0..240.0).contains(&hue) {
        (0.0, x, c)
    } else if (240.0..300.0).contains(&hue) {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };

    // Convert RGB to 8-bit values
    let r = ((r + m) * 255.0) as u8;
    let g = ((g + m) * 255.0) as u8;
    let b = ((b + m) * 255.0) as u8;

    [r, g, b]
}
