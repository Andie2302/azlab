//! Tests der Anzeige: Abbildungen und das Zeichnen in den Speicher.

use crate::raster::RasterCanvas;
use crate::{Batch, Camera, Canvas, Color, Plane};

/// Farbe des Bildpunkts `(x, y)`.
fn pixel(canvas: &RasterCanvas, x: u32, y: u32) -> Option<(u8, u8, u8)> {
    canvas.pixmap().pixel(x, y).map(|color| (color.red(), color.green(), color.blue()))
}

#[test]
fn colors() {
    let (a, b) = (Color::rgb(0, 100, 200), Color::rgb(200, 100, 0));
    assert_eq!(a.mix(b, 0.0), a);
    assert_eq!(a.mix(b, 1.0), b);
    assert_eq!(a.mix(b, 0.5), Color::rgb(100, 100, 100));
    // Außerhalb von 0 bis 1 wird begrenzt.
    assert_eq!(a.mix(b, 7.0), b);
    assert_eq!(Color::WHITE.shade(0.0), Color::rgb(0, 0, 0));
    assert_eq!(Color::RED.with_alpha(10).a, 10);
}

#[test]
fn raster() {
    let canvas = RasterCanvas::new(40, 30);
    assert!(RasterCanvas::new(0, 10).is_none());
    assert!(canvas.is_some_and(|mut canvas| {
        canvas.clear(Color::rgb(1, 2, 3));
        canvas.circle([10.0, 10.0], 4.0, Color::rgb(250, 0, 0));
        canvas.rect([20.0, 5.0], [10.0, 10.0], Color::rgb(0, 250, 0));
        canvas.triangle([2.0, 28.0], [12.0, 28.0], [7.0, 20.0], Color::rgb(0, 0, 250));
        canvas.line([30.0, 25.0], [38.0, 25.0], 3.0, Color::rgb(250, 250, 0));
        canvas.text([0.0, 0.0], "wird in Dateien ausgelassen", 10.0, Color::WHITE);
        (canvas.width(), canvas.height()) == (40.0, 30.0)
            && pixel(&canvas, 0, 0) == Some((1, 2, 3))
            && pixel(&canvas, 10, 10) == Some((250, 0, 0))
            && pixel(&canvas, 25, 10) == Some((0, 250, 0))
            && pixel(&canvas, 7, 25) == Some((0, 0, 250))
            && pixel(&canvas, 34, 25) == Some((250, 250, 0))
    }));
}

#[test]
fn plane() {
    let canvas = RasterCanvas::new(200, 100);
    assert!(canvas.is_some_and(|canvas| {
        // Ausschnitt 2 × 2 in ein Bild 200 × 100 mit 10 Bildpunkten Rand:
        // Die Höhe begrenzt, also 40 Bildpunkte je Einheit, waagrecht
        // mittig.
        let plane = Plane::fit(&canvas, [1.0, 1.0], [3.0, 3.0], 10.0);
        plane.length(1.0) == 40.0
            && plane.point([1.0, 1.0]) == [60.0, 90.0]
            && plane.point([3.0, 3.0]) == [140.0, 10.0]
    }));
}

#[test]
fn camera() {
    let canvas = RasterCanvas::new(200, 100);
    assert!(canvas.is_some_and(|mut canvas| {
        let camera = Camera::look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0);
        // Der Zielpunkt liegt in der Bildmitte, in der richtigen Tiefe.
        let centre = camera.project([0.0, 0.0, 0.0], &canvas);
        // Oben im Raum ist oben im Bild, rechts ist rechts.
        let above = camera.project([0.0, 1.0, 0.0], &canvas);
        let right = camera.project([1.0, 0.0, 0.0], &canvas);
        // Hinter der Kamera gibt es kein Bild.
        let behind = camera.project([0.0, 0.0, 6.0], &canvas);

        // Ein nahes Dreieck verdeckt ein fernes, gleich in welcher
        // Reihenfolge sie hinzugefügt wurden.
        let (near, far) = (Color::rgb(250, 0, 0), Color::rgb(0, 0, 250));
        let sheet = |z: f32| [[-1.0, -1.0, z], [1.0, -1.0, z], [0.0, 1.5, z]];
        let mut batch = Batch::new();
        batch.triangle(sheet(1.0), near, near);
        batch.triangle(sheet(-1.0), far, far);
        batch.quad([[9.0, 9.0, 0.0], [9.5, 9.0, 0.0], [9.5, 9.5, 0.0], [9.0, 9.5, 0.0]], 2, far, far);
        batch.sphere([0.0, -20.0, 0.0], 0.5, far);
        batch.floor(-30.0, 1.0, 2, far, far);
        canvas.clear(Color::rgb(0, 0, 0));
        batch.draw(&mut canvas, &camera, false);
        let front = pixel(&canvas, 100, 50);

        centre.is_some_and(|(point, depth)| point == [100.0, 50.0] && (depth - 5.0).abs() < 1.0e-5)
            && above.is_some_and(|(point, _)| point[0] == 100.0 && point[1] < 50.0)
            && right.is_some_and(|(point, _)| point[0] > 100.0 && point[1] == 50.0)
            && behind.is_none()
            && front.is_some_and(|(red, _, blue)| red > 100 && blue == 0)
    }));
}
