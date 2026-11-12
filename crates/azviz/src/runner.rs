//! Ablauf: Fenster oder Datei.

use crate::canvas::Canvas;
use crate::color::Color;
use crate::raster::RasterCanvas;
use crate::window::WindowCanvas;
use alloc::string::String;
use core::fmt;
use macroquad::prelude as mq;
use std::fs::File;
use std::io::BufWriter;

/// Eine Simulation, die sich zeigen lässt.
pub trait Scene {
    /// Titel des Fensters.
    fn title(&self) -> String;

    /// Rechnet um `dt` Sekunden weiter.
    fn update(&mut self, dt: f32);

    /// Zeichnet den aktuellen Zustand.
    fn draw(&self, canvas: &mut dyn Canvas);
}

/// Einstellungen der Anzeige.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Options {
    /// Breite in Bildpunkten.
    pub width: u32,
    /// Höhe in Bildpunkten.
    pub height: u32,
    /// Bilder je Sekunde; der Zeitschritt der Simulation ist ihr Kehrwert.
    pub fps: u32,
    /// Zahl der Bilder bei Ausgabe in eine Datei, wenn `--frames` fehlt.
    pub frames: u32,
}

impl Default for Options {
    fn default() -> Self {
        Self { width: 800, height: 600, fps: 60, frames: 180 }
    }
}

/// Fehler der Anzeige.
#[derive(Debug)]
pub enum VizError {
    /// Die Kommandozeile war nicht zu verstehen.
    Arguments(String),
    /// Die Datei ließ sich nicht schreiben.
    Io(std::io::Error),
    /// Das Bild ließ sich nicht kodieren.
    Encoding(String),
}

impl fmt::Display for VizError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Arguments(message) => write!(f, "Kommandozeile: {message}"),
            Self::Io(error) => write!(f, "Datei: {error}"),
            Self::Encoding(message) => write!(f, "Bild: {message}"),
        }
    }
}

impl core::error::Error for VizError {}

impl From<std::io::Error> for VizError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// Ziel der Ausgabe laut Kommandozeile.
enum Target {
    Window,
    Gif(String),
    Png(String),
}

fn parse(mut options: Options) -> Result<(Target, Options), VizError> {
    let mut target = Target::Window;
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        let mut value = |name: &str| {
            arguments.next().ok_or_else(|| VizError::Arguments(alloc::format!("{name} braucht einen Wert")))
        };
        match argument.as_str() {
            "--gif" => target = Target::Gif(value("--gif")?),
            "--png" => target = Target::Png(value("--png")?),
            "--frames" => {
                options.frames = value("--frames")?
                    .parse()
                    .map_err(|_| VizError::Arguments(String::from("--frames braucht eine Zahl")))?;
            }
            "--size" => {
                let text = value("--size")?;
                let size = text.split_once('x').and_then(|(w, h)| Some((w.parse().ok()?, h.parse().ok()?)));
                let Some((width, height)) = size else {
                    return Err(VizError::Arguments(String::from("--size braucht BREITExHÖHE, etwa 480x360")));
                };
                (options.width, options.height) = (width, height);
            }
            other => {
                return Err(VizError::Arguments(alloc::format!(
                    "unbekanntes Argument {other}; erlaubt sind --gif DATEI, --png DATEI, --frames N, --size BxH"
                )));
            }
        }
    }
    Ok((target, options))
}

/// Zeigt die Szene: ohne Argumente im Fenster, mit `--gif` oder `--png` als
/// Datei (siehe die Übersicht in der Crate-Dokumentation).
pub fn run<S: Scene + 'static>(scene: S, options: Options) -> Result<(), VizError> {
    let (target, options) = parse(options)?;
    let frames = options.frames;
    match target {
        Target::Window => {
            window(scene, options);
            Ok(())
        }
        Target::Gif(path) => gif(scene, options, frames, &path),
        Target::Png(path) => png(scene, options, frames, &path),
    }
}

fn window<S: Scene + 'static>(mut scene: S, options: Options) {
    let conf = mq::Conf {
        window_title: scene.title(),
        window_width: options.width as i32,
        window_height: options.height as i32,
        high_dpi: true,
        sample_count: 4,
        ..mq::Conf::default()
    };
    let dt = 1.0 / options.fps.max(1) as f32;
    macroquad::Window::from_config(conf, async move {
        let mut paused = false;
        let mut canvas = WindowCanvas;
        loop {
            if mq::is_key_pressed(mq::KeyCode::Escape) {
                break;
            }
            if mq::is_key_pressed(mq::KeyCode::Space) {
                paused = !paused;
            }
            if !paused {
                scene.update(dt);
            }
            scene.draw(&mut canvas);
            let hint = if paused { "angehalten – Leertaste: weiter, Esc: Ende" } else { "Leertaste: anhalten, Esc: Ende" };
            canvas.text([12.0, 22.0], hint, 20.0, Color::GRAY);
            mq::next_frame().await;
        }
    });
}

fn canvas(options: Options) -> Result<RasterCanvas, VizError> {
    RasterCanvas::new(options.width, options.height)
        .ok_or_else(|| VizError::Encoding(String::from("Bildgröße null")))
}

fn png<S: Scene>(mut scene: S, options: Options, frames: u32, path: &str) -> Result<(), VizError> {
    let dt = 1.0 / options.fps.max(1) as f32;
    for _ in 0..frames {
        scene.update(dt);
    }
    let mut canvas = canvas(options)?;
    scene.draw(&mut canvas);
    canvas.pixmap().save_png(path).map_err(|error| VizError::Encoding(alloc::format!("{error}")))
}

fn gif<S: Scene>(mut scene: S, options: Options, frames: u32, path: &str) -> Result<(), VizError> {
    let dt = 1.0 / options.fps.max(1) as f32;
    let encode = |error: gif::EncodingError| VizError::Encoding(alloc::format!("{error}"));
    let (width, height) = (options.width.min(0xFFFF) as u16, options.height.min(0xFFFF) as u16);
    let mut encoder = gif::Encoder::new(BufWriter::new(File::create(path)?), width, height, &[]).map_err(encode)?;
    encoder.set_repeat(gif::Repeat::Infinite).map_err(encode)?;
    let mut canvas = canvas(options)?;
    // GIF kennt Verzögerungen nur in Hundertstelsekunden: jedes zweite Bild
    // schreiben, damit die Datei klein bleibt und die Zeit stimmt.
    let skip = (options.fps / 30).max(1);
    let delay = (100 * skip / options.fps.max(1)).max(2) as u16;
    for frame in 0..frames {
        scene.update(dt);
        if frame % skip != 0 {
            continue;
        }
        scene.draw(&mut canvas);
        let mut pixels = canvas.pixmap().data().to_vec();
        let mut image = gif::Frame::from_rgba_speed(width, height, &mut pixels, 10);
        image.delay = delay;
        encoder.write_frame(&image).map_err(encode)?;
    }
    Ok(())
}
