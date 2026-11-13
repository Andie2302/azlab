//! Ablauf: Fenster oder Datei.

use crate::canvas::Canvas;
use crate::Scene;
use crate::color::Color;
use crate::raster::RasterCanvas;
use crate::window::WindowCanvas;
use alloc::string::String;
use core::fmt;
use macroquad::prelude as mq;
use std::fs::File;
use std::io::BufWriter;

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
    /// Fenster pausiert starten; Exporte rechnen unabhängig davon weiter.
    pub paused: bool,
    /// CPU-Zeiten für Physik und Zeichenaufbau im Fenster protokollieren.
    pub profile: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self { width: 800, height: 600, fps: 60, frames: 180, paused: false, profile: false }
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
            "--paused" => options.paused = true,
            "--profile" => options.profile = true,
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
                    "unbekanntes Argument {other}; erlaubt sind --gif DATEI, --png DATEI, --frames N, --size BxH, --paused, --profile"
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
        scene.set_profiling(options.profile);
        let mut paused = options.paused;
        let (mut update_ms, mut draw_ms, mut steps, mut frames) = (0.0, 0.0, 0u32, 0u32);
        let mut canvas = WindowCanvas::default();
        let mut state_message = String::new();
        loop {
            if mq::is_key_pressed(mq::KeyCode::Escape) {
                break;
            }
            if mq::is_key_pressed(mq::KeyCode::Space) {
                paused = !paused;
            }
            // N führt genau einen Schritt aus und lässt das Fenster pausiert.
            let mut single = mq::is_key_pressed(mq::KeyCode::N);
            if single { paused = true; }
            if mq::is_key_pressed(mq::KeyCode::G) { scene.toggle_ground(); }
            for (key, variant) in [(mq::KeyCode::Key1, 0), (mq::KeyCode::Key2, 1),
                (mq::KeyCode::Key3, 2), (mq::KeyCode::Key4, 3), (mq::KeyCode::Key5, 4), (mq::KeyCode::Key6, 5)] {
                if mq::is_key_pressed(key) { scene.select_variant(variant); state_message.clear(); paused = true; single = false; }
            }
            // Vor dem Update speichern/laden: L stellt exakt das gespeicherte
            // Bild wieder her und lässt die Simulation angehalten.
            if mq::is_key_pressed(mq::KeyCode::S) || mq::is_key_pressed(mq::KeyCode::L) {
                if let Some(path) = scene.state_path().map(std::path::PathBuf::from) {
                    let loading = mq::is_key_pressed(mq::KeyCode::L);
                    single = false;
                    paused = true;
                    let result = if loading { scene.load_state(&path) } else { scene.save_state(&path) };
                    state_message = match result {
                        Ok(()) => format!("{}: {}", if loading { "Geladen" } else { "Gespeichert" }, path.display()),
                        Err(error) => format!("Zustand: {error}"),
                    };
                    eprintln!("{state_message}");
                } else {
                    state_message = "Diese Szene unterstützt S/L nicht".into();
                }
            }
            let yaw = (i32::from(mq::is_key_down(mq::KeyCode::Right))
                - i32::from(mq::is_key_down(mq::KeyCode::Left))) as f32 * 0.025;
            let pitch = (i32::from(mq::is_key_down(mq::KeyCode::Up))
                - i32::from(mq::is_key_down(mq::KeyCode::Down))) as f32 * 0.025;
            let zoom = (i32::from(mq::is_key_down(mq::KeyCode::PageDown))
                - i32::from(mq::is_key_down(mq::KeyCode::PageUp))) as f32 * 0.025;
            if yaw != 0.0 || pitch != 0.0 || zoom != 0.0 { scene.orbit(yaw, pitch, zoom); }
            let start = std::time::Instant::now();
            if should_update(paused, single) {
                scene.update(dt);
                steps += 1;
            }
            let draw_start = std::time::Instant::now();
            if should_update(paused, single) {
                update_ms += (draw_start - start).as_secs_f64() * 1000.0;
            }
            scene.draw(&mut canvas);
            draw_ms += draw_start.elapsed().as_secs_f64() * 1000.0;
            frames += 1;
            if options.profile && frames == 120 {
                eprintln!("azviz CPU: Physik {:.3} ms/Schritt ({}), Zeichnen {:.3} ms/Bild; ohne GPU/Present",
                    update_ms / f64::from(steps.max(1)), steps, draw_ms / f64::from(frames));
                if let Some(details) = scene.profile_details() { eprintln!("{details}"); }
                (update_ms, draw_ms, steps, frames) = (0.0, 0.0, 0, 0);
            }
            let hint = if paused { "Pause: Leertaste weiter, N Schritt, Pfeile Kamera, Bild auf/ab Zoom, Esc Ende" } else { "Leertaste Pause, N Schritt, Pfeile Kamera, Esc Ende" };
            canvas.text([12.0, 22.0], hint, 20.0, Color::GRAY);
            if !state_message.is_empty() {
                let message_y = canvas.height() - 12.0;
                canvas.text([12.0, message_y], &state_message, 18.0, Color::WHITE);
            }
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

fn should_update(paused: bool, single: bool) -> bool {
    !paused || single
}

#[cfg(test)]
mod control_tests {
    use super::*;
    #[test]
    fn pause_and_single_step() {
        assert!(!should_update(true, false));
        assert!(should_update(true, true));
        assert!(should_update(false, false));
        assert!(!Options::default().paused);
    }
}
