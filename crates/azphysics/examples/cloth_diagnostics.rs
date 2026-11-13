//! Tuchdiagnose im Fenster oder als eigenständiger HTML-/CSV-Bericht.
//! `cargo demo cloth_diagnostics -- --help` zeigt die Aufrufe.

extern crate alloc;

use alloc::collections::VecDeque;
use azmath::Vector;
use azphysics::cloth::{Cloth, ClothDiagnostics, Collider, StepStats};
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};
use core::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

const HISTORY: usize = 1200;
const NAMES: [&str; 4] = [
    "Freier Flug",
    "Rotation",
    "Hängendes Tuch",
    "Fall auf Kugel",
];
const NOTES: [&str; 4] = [
    "Keine äußeren Kräfte, keine Dämpfung: Energie, Impuls und Drehimpuls sollten konstant bleiben.",
    "Freie Rotation ohne Dämpfung: Impuls und Drehimpuls sollten konstant bleiben; Energieänderungen zeigen numerische Dissipation oder Energiezufuhr.",
    "Haltepunkte und Dämpfung wirken von außen. Energie und Impuls müssen hier nicht erhalten bleiben.",
    "Schwerkraft, Reibung und unelastische Kontakte: Impuls ist nicht erhalten, mechanische Energie darf abnehmen.",
];

#[derive(Clone, Copy)]
struct Sample {
    time: f64,
    metrics: ClothDiagnostics,
    stats: StepStats,
    step_ms: f64,
}

struct Diagnostics {
    cloth: Cloth<f32>,
    colliders: Vec<Collider<f32>>,
    camera: Camera,
    variant: usize,
    time: f64,
    baseline: ClothDiagnostics,
    history: VecDeque<Sample>,
    limit: usize,
}

impl Diagnostics {
    fn new(variant: usize, limit: usize) -> Self {
        let variant = variant.min(3);
        let size = 14;
        let spacing = 1.2 / (size - 1) as f32;
        let (origin, down) = if variant == 2 {
            (Vector([-0.6, 1.5, 0.0]), Vector([0.0, -spacing, 0.0]))
        } else {
            (Vector([-0.6, 1.4, -0.6]), Vector([0.0, 0.0, spacing]))
        };
        let mut cloth = Cloth::grid(size, size, origin, Vector([spacing, 0.0, 0.0]), down, 0.01);
        cloth.params.gravity = Vector([0.0; 3]);
        cloth.params.damping = 0.0;
        cloth.params.friction = 0.0;
        cloth.params.substeps = 8;
        match variant {
            0 => cloth.set_all_velocities(Vector([0.7, 0.2, -0.3])),
            1 => {
                for index in 0..cloth.positions().len() {
                    let p = cloth.positions()[index];
                    cloth.set_velocity(index, Vector([2.0 * p[2], 0.0, -2.0 * p[0]]));
                }
            }
            2 => {
                cloth.pin(cloth.index(0, 0));
                cloth.pin(cloth.index(size - 1, 0));
                cloth.params.gravity = Vector([0.0, -9.81, 0.0]);
                cloth.params.damping = 0.3;
            }
            _ => {
                cloth.params.gravity = Vector([0.0, -9.81, 0.0]);
                cloth.params.friction = 0.3;
            }
        }
        let colliders = if variant == 3 {
            vec![
                Collider::Sphere {
                    center: Vector([0.0, 0.45, 0.0]),
                    radius: 0.45,
                },
                Collider::Plane {
                    point: Vector([0.0; 3]),
                    normal: Vector([0.0, 1.0, 0.0]),
                },
            ]
        } else {
            Vec::new()
        };
        let baseline = cloth.diagnostics();
        let mut scene = Self {
            cloth,
            colliders,
            variant,
            time: 0.0,
            baseline,
            camera: Camera::look_at([2.6, 2.4, 3.2], [0.0, 0.7, 0.0], 40.0),
            history: VecDeque::new(),
            limit: limit.max(2),
        };
        scene.record(StepStats::default(), 0.0);
        scene
    }

    fn record(&mut self, stats: StepStats, step_ms: f64) {
        if self.history.len() == self.limit {
            self.history.pop_front();
        }
        self.history.push_back(Sample {
            time: self.time,
            metrics: self.cloth.diagnostics(),
            stats,
            step_ms,
        });
    }

    fn reset_baseline(&mut self) {
        self.baseline = self.cloth.diagnostics();
        self.history.clear();
        self.record(StepStats::default(), 0.0);
    }

    fn csv(&self) -> String {
        let mut out = String::from(
            "time_s,finite,mass_kg,kinetic_J,gravity_J,elastic_J,mechanical_J,px_kg_m_s,py_kg_m_s,pz_kg_m_s,lx_com_kg_m2_s,ly_com_kg_m2_s,lz_com_kg_m2_s,max_strain,rms_strain,max_constraint_m,ccd_hits,ccd_undecided,full_rollbacks,repulsions,proximity_rebuilds,step_ms\n",
        );
        for s in &self.history {
            let m = s.metrics;
            let _ = writeln!(
                out,
                "{:.9},{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{},{},{},{},{},{:.6}",
                s.time,
                m.finite,
                m.mass,
                m.kinetic_energy,
                m.gravitational_energy,
                m.elastic_energy,
                m.mechanical_energy,
                m.momentum[0],
                m.momentum[1],
                m.momentum[2],
                m.angular_momentum[0],
                m.angular_momentum[1],
                m.angular_momentum[2],
                m.max_strain,
                m.rms_strain,
                m.max_constraint_error,
                s.stats.ccd_hits,
                s.stats.ccd_undecided,
                s.stats.full_rollbacks,
                s.stats.repulsions,
                s.stats.proximity_rebuilds,
                s.step_ms
            );
        }
        out
    }

    fn report(&self, path: &Path) -> Result<(), VizError> {
        let mut data = String::from("[");
        for (index, s) in self.history.iter().enumerate() {
            if index > 0 {
                data.push(',');
            }
            let m = s.metrics;
            data.push('[');
            for (column, value) in [
                s.time,
                m.kinetic_energy,
                m.gravitational_energy,
                m.elastic_energy,
                m.mechanical_energy,
                m.momentum[0],
                m.momentum[1],
                m.momentum[2],
                m.angular_momentum[0],
                m.angular_momentum[1],
                m.angular_momentum[2],
                m.max_strain * 100.0,
                m.rms_strain * 100.0,
                m.max_constraint_error,
                s.step_ms,
                f64::from(s.stats.ccd_hits),
                f64::from(s.stats.full_rollbacks),
                if m.finite { 1.0 } else { 0.0 },
            ]
            .iter()
            .enumerate()
            {
                if column > 0 {
                    data.push(',');
                }
                if value.is_finite() {
                    let _ = write!(data, "{value:.12e}");
                } else {
                    data.push_str("null");
                }
            }
            data.push(']');
        }
        data.push(']');
        let params = format!(
            "f32 · {} Teilchen · {} Teilschritte · {} Durchläufe · Dämpfung {:.4}/s · Selbstkollision {} · CCD {} · g={:?} m/s² · Nachgiebigkeit Dehnung/Biegung={:.3e}/{:.3e} · Reibung={} · Dicke={} m · dt={:.9} s",
            self.cloth.positions().len(),
            self.cloth.params.substeps,
            self.cloth.params.iterations,
            self.cloth.params.damping,
            self.cloth.params.self_collision,
            self.cloth.params.ccd,
            self.cloth.params.gravity.0,
            self.cloth.params.stretch_compliance,
            self.cloth.params.bend_compliance,
            self.cloth.params.friction,
            self.cloth.params.thickness,
            self.history
                .get(1)
                .map_or(0.0, |s| s.time - self.history[0].time)
        );
        let html = include_str!("common/cloth_diagnostics.html")
            .replace("__TITLE__", NAMES[self.variant])
            .replace("__NOTE__", NOTES[self.variant])
            .replace("__PARAMS__", &params)
            .replace("__DATA__", &data)
            .replace("__CSV__", &format!("{:?}", self.csv()));
        // create_new schützt vorhandene Messergebnisse vor versehentlichem Überschreiben.
        let csv_path = path.with_extension("csv");
        if path == csv_path {
            return Err(VizError::Arguments(
                "Bericht braucht eine andere Endung als .csv".into(),
            ));
        }
        let mut html_file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        let mut csv_file = match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&csv_path)
        {
            Ok(file) => file,
            Err(error) => {
                drop(html_file);
                let _ = std::fs::remove_file(path);
                return Err(error.into());
            }
        };
        html_file.write_all(html.as_bytes())?;
        csv_file.write_all(self.csv().as_bytes())?;
        Ok(())
    }
}

type Track = (&'static str, Color, fn(&ClothDiagnostics) -> f64);
const ENERGY: [Track; 4] = [
    ("Gesamt", Color::WHITE, |d| d.mechanical_energy),
    ("Kinetisch", Color::ORANGE, |d| d.kinetic_energy),
    ("Gravitation", Color::BLUE, |d| d.gravitational_energy),
    ("Elastisch", Color::GREEN, |d| d.elastic_energy),
];
const MOMENTUM: [Track; 3] = [
    ("x", Color::RED, |d| d.momentum[0]),
    ("y", Color::GREEN, |d| d.momentum[1]),
    ("z", Color::BLUE, |d| d.momentum[2]),
];
const ANGULAR: [Track; 3] = [
    ("x", Color::RED, |d| d.angular_momentum[0]),
    ("y", Color::GREEN, |d| d.angular_momentum[1]),
    ("z", Color::BLUE, |d| d.angular_momentum[2]),
];
const STRAIN: [Track; 2] = [
    ("Max %", Color::ORANGE, |d| d.max_strain * 100.0),
    ("RMS %", Color::BLUE, |d| d.rms_strain * 100.0),
];

fn range(values: impl Iterator<Item = f64>) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for x in values.filter(|x| x.is_finite()) {
        lo = lo.min(x);
        hi = hi.max(x);
    }
    if !lo.is_finite() {
        return (-1.0, 1.0);
    }
    let padding = ((hi - lo) * 0.08)
        .max(hi.abs().max(lo.abs()) * 1e-6)
        .max(1e-12);
    (lo - padding, hi + padding)
}

fn plot(
    canvas: &mut dyn Canvas,
    history: &VecDeque<Sample>,
    bounds: [f32; 4],
    title: &str,
    tracks: &[Track],
) {
    let [x, y, w, h] = bounds;
    canvas.rect([x, y], [w, h], Color::rgb(25, 32, 44));
    canvas.text([x + 10.0, y + 22.0], title, 18.0, Color::WHITE);
    let (lo, hi) = range(
        history
            .iter()
            .filter(|s| s.metrics.finite)
            .flat_map(|s| tracks.iter().map(move |t| t.2(&s.metrics))),
    );
    let Some(first) = history.front() else {
        return;
    };
    let Some(last) = history.back() else {
        return;
    };
    let span = (last.time - first.time).max(1e-9);
    let point = |s: &Sample, value: f64| {
        [
            x + 78.0 + ((s.time - first.time) / span) as f32 * (w - 94.0),
            y + 38.0 + ((hi - value) / (hi - lo)) as f32 * (h - 76.0),
        ]
    };
    for tick in 0..3 {
        let v = lo + (hi - lo) * f64::from(tick) / 2.0;
        let p = point(first, v);
        canvas.line(p, [x + w - 16.0, p[1]], 1.0, Color::rgb(45, 53, 68));
        canvas.text(
            [x + 5.0, p[1] + 4.0],
            &format!("{v:.2e}"),
            12.0,
            Color::GRAY,
        );
    }
    for (i, (name, color, read)) in tracks.iter().enumerate() {
        canvas.text(
            [x + 10.0 + i as f32 * 100.0, y + h - 8.0],
            name,
            13.0,
            *color,
        );
        let mut previous = None;
        for s in history {
            let value = read(&s.metrics);
            if !s.metrics.finite || !value.is_finite() {
                previous = None;
                continue;
            }
            let next = point(s, value);
            if let Some(p) = previous {
                canvas.line(p, next, 1.5, *color);
            }
            previous = Some(next);
        }
    }
    canvas.text(
        [x + 78.0, y + h - 27.0],
        &format!("{:.2} – {:.2} s", first.time, last.time),
        12.0,
        Color::GRAY,
    );
}

impl Scene for Diagnostics {
    fn title(&self) -> String {
        format!("azphysics – Tuchdiagnose: {}", NAMES[self.variant])
    }

    fn update(&mut self, dt: f32) {
        if self.history.back().is_some_and(|s| !s.metrics.finite) {
            return;
        }
        let start = std::time::Instant::now();
        let stats = self.cloth.step(dt, &self.colliders);
        let step_ms = start.elapsed().as_secs_f64() * 1000.0;
        self.time += f64::from(dt);
        self.record(stats, step_ms);
    }

    fn select_variant(&mut self, variant: u32) {
        if variant < 4 {
            *self = Self::new(variant as usize, self.limit);
        }
    }
    fn orbit(&mut self, yaw: f32, pitch: f32, zoom: f32) {
        self.camera.orbit(yaw, pitch, zoom);
    }
    fn state_path(&self) -> Option<&str> {
        Some("FreierFall-Daten.csv")
    }
    fn save_state(&self, path: &Path) -> Result<(), VizError> {
        std::fs::write(path, self.csv())?;
        Ok(())
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let Some(last) = self.history.back() else {
            return;
        };
        // Freie Translation in der Ansicht herausrechnen, nicht in der Messung.
        let offset = if self.variant == 0 {
            Vector([0.7, 0.2, -0.3]) * self.time as f32
        } else {
            Vector([0.0; 3])
        };
        if last.metrics.finite {
            let positions: Vec<_> = self
                .cloth
                .positions()
                .iter()
                .map(|p| (*p - offset).0)
                .collect();
            let mut batch = Batch::new();
            batch.mesh(
                &positions,
                self.cloth.triangles(),
                Color::ORANGE,
                Color::RED,
            );
            if self.variant == 3 {
                batch.sphere([0.0, 0.45, 0.0], 0.44, Color::GRAY);
            }
            batch.draw(canvas, &self.camera, true);
        }
        canvas.text(
            [12.0, 50.0],
            "1 Flug | 2 Rotation | 3 Hängend | 4 Kugel | S CSV (letzte 1200 Werte)",
            18.0,
            Color::WHITE,
        );
        canvas.text([12.0, 75.0], NAMES[self.variant], 22.0, Color::ORANGE);
        let d = last.metrics;
        let de = d.mechanical_energy - self.baseline.mechanical_energy;
        let dp: [f64; 3] = core::array::from_fn(|i| d.momentum[i] - self.baseline.momentum[i]);
        let dl: [f64; 3] =
            core::array::from_fn(|i| d.angular_momentum[i] - self.baseline.angular_momentum[i]);
        canvas.text([12.0, 102.0], &format!("t={:.3} s | ΔE={:+.3e} J | Δp=({:+.2e}, {:+.2e}, {:+.2e}) | ΔL=({:+.2e}, {:+.2e}, {:+.2e})",
            self.time, de, dp[0], dp[1], dp[2], dl[0], dl[1], dl[2]), 16.0, Color::WHITE);
        canvas.text(
            [12.0, 126.0],
            &format!(
                "Constraint max={:.3e} m | CCD={} | Rücksetzungen={} | Physik={:.2} ms | gültig={}",
                d.max_constraint_error,
                last.stats.ccd_hits,
                last.stats.full_rollbacks,
                last.step_ms,
                d.finite
            ),
            16.0,
            if d.finite { Color::WHITE } else { Color::RED },
        );
        let (w, h) = (canvas.width(), canvas.height());
        let pw = ((w - 36.0) / 2.0).max(100.0);
        let ph = ((h * 0.54 - 36.0) / 2.0).max(90.0);
        let y = h - 2.0 * ph - 40.0;
        plot(
            canvas,
            &self.history,
            [12.0, y, pw, ph],
            "Energie [J]",
            &ENERGY,
        );
        plot(
            canvas,
            &self.history,
            [24.0 + pw, y, pw, ph],
            "Impuls [kg m/s]",
            &MOMENTUM,
        );
        plot(
            canvas,
            &self.history,
            [12.0, y + ph + 12.0, pw, ph],
            "Drehimpuls um Schwerpunkt [kg m²/s]",
            &ANGULAR,
        );
        plot(
            canvas,
            &self.history,
            [24.0 + pw, y + ph + 12.0, pw, ph],
            "Relative Constraint-Abweichung [%]",
            &STRAIN,
        );
        canvas.text([12.0, y - 12.0], NOTES[self.variant], 13.0, Color::WHITE);
    }
}

fn args_error(message: &str) -> VizError {
    VizError::Arguments(message.into())
}

fn export(args: &[String]) -> Result<(), VizError> {
    let (mut path, mut frames, mut dt, mut variant) = (None, 600usize, 1.0f32 / 60.0, 1usize);
    let (mut substeps, mut iterations, mut damping) = (None, None, None);
    let (mut self_collision, mut ccd) = (true, true);
    let mut args = args.iter();
    while let Some(arg) = args.next() {
        let mut value = || {
            args.next()
                .ok_or_else(|| args_error("Argument braucht einen Wert"))
        };
        match arg.as_str() {
            "--report" => path = Some(value()?.clone()),
            "--frames" => {
                frames = value()?
                    .parse()
                    .map_err(|_| args_error("--frames braucht eine ganze Zahl"))?
            }
            "--dt" => {
                dt = value()?
                    .parse()
                    .map_err(|_| args_error("--dt braucht Sekunden"))?
            }
            "--scenario" => {
                variant = match value()?.as_str() {
                    "flight" => 0,
                    "rotation" => 1,
                    "hanging" => 2,
                    "sphere" => 3,
                    _ => {
                        return Err(args_error(
                            "--scenario: flight, rotation, hanging oder sphere",
                        ));
                    }
                }
            }
            "--substeps" => {
                substeps = Some(
                    value()?
                        .parse::<u32>()
                        .map_err(|_| args_error("--substeps braucht eine ganze Zahl"))?,
                )
            }
            "--iterations" => {
                iterations = Some(
                    value()?
                        .parse::<u32>()
                        .map_err(|_| args_error("--iterations braucht eine ganze Zahl"))?,
                )
            }
            "--damping" => {
                damping = Some(
                    value()?
                        .parse::<f32>()
                        .map_err(|_| args_error("--damping braucht eine Zahl"))?,
                )
            }
            "--no-self-collision" => self_collision = false,
            "--no-ccd" => ccd = false,
            _ => return Err(args_error("Unbekanntes Berichtsargument; siehe --help")),
        }
    }
    if frames == 0
        || frames > 1_000_000
        || !dt.is_finite()
        || dt <= 0.0
        || substeps == Some(0)
        || damping.is_some_and(|x| !x.is_finite() || x < 0.0)
    {
        return Err(args_error(
            "Ungültige Parameter: frames 1..1000000, dt/substeps > 0, damping endlich und >= 0",
        ));
    }
    let path = path.ok_or_else(|| args_error("--report braucht einen Dateipfad"))?;
    let mut scene = Diagnostics::new(variant, frames + 1);
    if let Some(n) = substeps {
        scene.cloth.params.substeps = n;
    }
    if let Some(n) = iterations {
        scene.cloth.params.iterations = n;
    }
    if let Some(d) = damping {
        scene.cloth.params.damping = d;
    }
    scene.cloth.params.self_collision = self_collision;
    scene.cloth.params.ccd = ccd;
    scene.reset_baseline();
    for _ in 0..frames {
        scene.update(dt);
        if scene.history.back().is_some_and(|s| !s.metrics.finite) {
            break;
        }
    }
    scene.report(Path::new(&path))?;
    eprintln!(
        "Bericht: {path}; Messwerte: {} ({} Stichproben einschließlich t=0)",
        Path::new(&path).with_extension("csv").display(),
        scene.history.len()
    );
    if scene.history.back().is_some_and(|s| !s.metrics.finite) {
        return Err(args_error(
            "Simulation wegen ungültiger Messwerte angehalten; Bericht wurde gespeichert",
        ));
    }
    Ok(())
}

fn main() -> Result<(), VizError> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help") {
        println!(
            "Tuchdiagnose\nFenster: cargo demo cloth_diagnostics\nBericht: cargo demo cloth_diagnostics -- --report diagnose.html --scenario rotation --frames 600\nBerichtsoptionen: --dt SEKUNDEN --substeps N --iterations N --damping WERT --no-self-collision --no-ccd\nSzenarien: flight, rotation (Standard), hanging, sphere\nFenster: 1–4 Versuch zurücksetzen; Leertaste Pause; N Einzelschritt; S CSV; Pfeile Kamera.\nBerichte überschreiben keine bestehenden Dateien. CSV enthält SI-Einheiten, Dehnungen als Verhältnis.\nFensteroptionen: --paused --png DATEI --gif DATEI --frames N --size BxH"
        );
        return Ok(());
    }
    if args.iter().any(|a| a == "--report") {
        return export(&args);
    }
    azviz::run(
        Diagnostics::new(1, HISTORY),
        Options {
            width: 1280,
            height: 900,
            paused: true,
            ..Options::default()
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn history_has_baseline_and_only_advances_with_steps() {
        let mut scene = Diagnostics::new(0, 3);
        let initial = scene.baseline;
        for _ in 0..4 {
            scene.update(1.0 / 60.0);
        }
        assert_eq!(scene.history.len(), 3);
        assert!(scene.history.front().is_some_and(|s| s.time > 0.0));
        assert!(scene.history.back().is_some_and(|last| {
            last.metrics.finite
                && (0..3)
                    .all(|axis| (last.metrics.momentum[axis] - initial.momentum[axis]).abs() < 1e-3)
        }));
        assert_eq!(scene.csv().lines().count(), 4);
        scene.select_variant(2);
        assert_eq!(scene.time, 0.0);
        assert_eq!(scene.history.len(), 1);
    }

    #[test]
    fn invalid_values_stop_recording() {
        let mut scene = Diagnostics::new(0, 3);
        scene.cloth.set_position(0, Vector([f32::NAN; 3]));
        scene.record(StepStats::default(), 0.0);
        scene.update(1.0 / 60.0);
        assert_eq!(scene.time, 0.0);
        assert!(scene.history.back().is_some_and(|s| !s.metrics.finite));
    }

    #[test]
    fn ranges_include_negative_constant_and_invalid_data() {
        let (lo, hi) = range([f64::NAN, -2.0, -2.0].into_iter());
        assert!(lo < -2.0 && hi > -2.0);
        assert_eq!(range([f64::NAN].into_iter()), (-1.0, 1.0));
    }
}
