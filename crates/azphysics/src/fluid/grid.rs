//! Navier-Stokes auf dem versetzten Gitter.

use super::Vec2;
use super::advect::{Flow, Layout, Scratch};
use super::pressure::{Preconditioner, PressureOperator, Stencil};
use alloc::vec;
use alloc::vec::Vec;
use azmath::Vector;
use azmath::scalar::RealElementary;
use azmath::solve::{SolveError, Tolerance, preconditioned_conjugate_gradient};

/// Was an einem Rand des Rechengebiets gilt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Edge<T> {
    /// Undurchlässige Wand.
    Wall,
    /// Einlass: Die Flüssigkeit strömt mit dieser Geschwindigkeit senkrecht
    /// zum Rand herein (negativ: hinaus).
    Inflow(T),
    /// Offen: Die Flüssigkeit kann frei ein- und ausströmen; der Druck
    /// draußen ist null.
    Open,
}

/// Die vier Ränder des Rechengebiets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edges<T> {
    /// Linker Rand.
    pub left: Edge<T>,
    /// Rechter Rand.
    pub right: Edge<T>,
    /// Unterer Rand.
    pub bottom: Edge<T>,
    /// Oberer Rand.
    pub top: Edge<T>,
}

impl<T> Edges<T> {
    /// Rundum geschlossener Kasten.
    pub const fn closed() -> Self {
        Self { left: Edge::Wall, right: Edge::Wall, bottom: Edge::Wall, top: Edge::Wall }
    }

    /// Kanal: Einlass links mit der Geschwindigkeit `speed`, rechts offen,
    /// oben und unten Wände.
    pub const fn channel(speed: T) -> Self {
        Self { left: Edge::Inflow(speed), right: Edge::Open, bottom: Edge::Wall, top: Edge::Wall }
    }
}

/// Einstellungen der Strömung.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridParams<T> {
    /// Kinematische Zähigkeit; null ist reibungsfrei.
    pub viscosity: T,
    /// Auftrieb: Beschleunigung nach oben je Einheit Rauchdichte.
    pub buoyancy: T,
    /// Anteil des Rauchs, der je Sekunde verschwindet.
    pub dissipation: T,
    /// Transport mit Fehlerkorrektur: schärfere Wirbel, etwas mehr Aufwand.
    pub maccormack: bool,
    /// Höchstzahl der Iterationen für die Druckgleichung.
    pub pressure_iterations: u32,
    /// Verlangte Genauigkeit der Druckgleichung, bezogen auf die rechte
    /// Seite.
    pub pressure_tolerance: T,
}

impl<T: RealElementary> Default for GridParams<T> {
    fn default() -> Self {
        Self {
            viscosity: T::ZERO,
            buoyancy: T::ZERO,
            dissipation: T::ZERO,
            maccormack: true,
            pressure_iterations: 200,
            pressure_tolerance: T::from_f64(1.0e-4),
        }
    }
}

/// Was die Druckkorrektur eines Zeitschritts ergab.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridStats<T> {
    /// Iterationen des Lösers.
    pub iterations: u32,
    /// Restfehler der Druckgleichung.
    pub residual: T,
    /// `false`, wenn die verlangte Genauigkeit nicht erreicht wurde; das
    /// Feld ist dann nur näherungsweise quellenfrei.
    pub converged: bool,
}

/// Geschwindigkeitskomponente auf dem versetzten Gitter.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Component {
    U,
    V,
}

/// Inkompressible Strömung in der Ebene.
///
/// Das Gebiet besteht aus `nx × ny` quadratischen Zellen; Zelle `(i, j)`
/// liegt in Spalte `i` von links und Zeile `j` von unten. Die
/// Geschwindigkeit sitzt auf den Zellkanten, Druck und Rauch in der
/// Zellmitte.
#[derive(Clone, Debug)]
pub struct GridFluid<T> {
    nx: usize,
    ny: usize,
    spacing: T,
    edges: Edges<T>,
    u: Vec<T>,
    v: Vec<T>,
    smoke: Vec<T>,
    solid: Vec<bool>,
    pressure: Vec<T>,
    // Arbeitsspeicher.
    next_u: Vec<T>,
    next_v: Vec<T>,
    next_smoke: Vec<T>,
    scratch: Scratch<T>,
    initial: Vec<T>,
    free: Vec<bool>,
    rhs: Vec<T>,
    stencil: Stencil<T>,
    /// Die Hindernisse haben sich geändert; die Druckgleichung muss neu
    /// aufgebaut werden.
    stale: bool,
    work: Vec<T>,
    /// Einstellungen; dürfen zwischen den Schritten geändert werden.
    pub params: GridParams<T>,
}

impl<T: RealElementary> GridFluid<T> {
    /// Ruhende Flüssigkeit auf `nx × ny` Zellen der Kantenlänge `spacing`.
    pub fn new(nx: usize, ny: usize, spacing: T, edges: Edges<T>) -> Self {
        let (nx, ny) = (nx.max(1), ny.max(1));
        let spacing = if spacing > T::ZERO && spacing.is_finite() { spacing } else { T::ONE };
        let cells = nx * ny;
        let mut fluid = Self {
            nx,
            ny,
            spacing,
            edges,
            u: vec![T::ZERO; (nx + 1) * ny],
            v: vec![T::ZERO; nx * (ny + 1)],
            smoke: vec![T::ZERO; cells],
            solid: vec![false; cells],
            pressure: vec![T::ZERO; cells],
            next_u: vec![T::ZERO; (nx + 1) * ny],
            next_v: vec![T::ZERO; nx * (ny + 1)],
            next_smoke: vec![T::ZERO; cells],
            scratch: Scratch::new((nx + 1) * (ny + 1)),
            initial: vec![T::ZERO; (nx + 1) * (ny + 1)],
            free: vec![false; (nx + 1) * (ny + 1)],
            rhs: vec![T::ZERO; cells],
            stencil: Stencil::new(nx, ny),
            stale: true,
            work: vec![T::ZERO; 4 * cells],
            params: GridParams::default(),
        };
        fluid.enforce_boundaries();
        fluid
    }

    /// Zahl der Zellen in der Breite.
    #[inline]
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Zahl der Zellen in der Höhe.
    #[inline]
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Kantenlänge einer Zelle.
    #[inline]
    pub fn spacing(&self) -> T {
        self.spacing
    }

    #[inline]
    fn cell(&self, i: usize, j: usize) -> Option<usize> {
        (i < self.nx && j < self.ny).then(|| j * self.nx + i)
    }

    #[inline]
    fn is_fluid(&self, i: usize, j: usize) -> bool {
        self.cell(i, j).is_some_and(|index| !self.solid[index])
    }

    /// Rauchdichte aller Zellen, zeilenweise von unten.
    #[inline]
    pub fn smoke(&self) -> &[T] {
        &self.smoke
    }

    /// Rauchdichte der Zelle; null außerhalb des Gebiets.
    pub fn smoke_at(&self, i: usize, j: usize) -> T {
        self.cell(i, j).map_or(T::ZERO, |index| self.smoke[index])
    }

    /// Gibt Rauch in eine Zelle.
    pub fn add_smoke(&mut self, i: usize, j: usize, amount: T) {
        if let Some(index) = self.cell(i, j) {
            if !self.solid[index] {
                self.smoke[index] += amount;
            }
        }
    }

    /// `true`, wenn die Zelle ein Hindernis ist.
    pub fn is_solid(&self, i: usize, j: usize) -> bool {
        self.cell(i, j).is_some_and(|index| self.solid[index])
    }

    /// Macht die Zelle zum Hindernis oder wieder frei.
    pub fn set_solid(&mut self, i: usize, j: usize, solid: bool) {
        if let Some(index) = self.cell(i, j) {
            self.solid[index] = solid;
            if solid {
                self.smoke[index] = T::ZERO;
                self.pressure[index] = T::ZERO;
            }
            self.stale = true;
            self.enforce_boundaries();
        }
    }

    /// Macht alle Zellen zum Hindernis, deren Mitte in der Kreisscheibe
    /// liegt.
    pub fn add_solid_disc(&mut self, centre: Vec2<T>, radius: T) {
        let half = T::from_f64(0.5);
        for index in 0..self.solid.len() {
            let x = (T::from_f64((index % self.nx) as f64) + half) * self.spacing;
            let y = (T::from_f64((index / self.nx) as f64) + half) * self.spacing;
            if (Vector([x, y]) - centre).norm() <= radius {
                self.solid[index] = true;
                self.smoke[index] = T::ZERO;
                self.pressure[index] = T::ZERO;
            }
        }
        self.stale = true;
        self.enforce_boundaries();
    }

    /// Geschwindigkeit in der Mitte der Zelle; null außerhalb des Gebiets.
    pub fn velocity(&self, i: usize, j: usize) -> Vec2<T> {
        if i >= self.nx || j >= self.ny {
            return Vector::ZERO;
        }
        let half = T::from_f64(0.5);
        Vector([
            (self.u[j * (self.nx + 1) + i] + self.u[j * (self.nx + 1) + i + 1]) * half,
            (self.v[j * self.nx + i] + self.v[(j + 1) * self.nx + i]) * half,
        ])
    }

    /// Ändert die Geschwindigkeit der Zelle um `delta`, etwa für ein Gebläse.
    pub fn add_velocity(&mut self, i: usize, j: usize, delta: Vec2<T>) {
        if !self.is_fluid(i, j) {
            return;
        }
        for offset in 0..2 {
            if !self.fixed(Component::U, i + offset, j) {
                self.u[j * (self.nx + 1) + i + offset] += delta[0];
            }
            if !self.fixed(Component::V, i, j + offset) {
                self.v[(j + offset) * self.nx + i] += delta[1];
            }
        }
    }

    /// Quellstärke der Zelle: was netto aus ihr herausströmt, bezogen auf
    /// ihre Fläche. Nach jedem Schritt nahe null.
    pub fn divergence(&self, i: usize, j: usize) -> T {
        if !self.is_fluid(i, j) {
            return T::ZERO;
        }
        (self.u[j * (self.nx + 1) + i + 1] - self.u[j * (self.nx + 1) + i] + self.v[(j + 1) * self.nx + i]
            - self.v[j * self.nx + i])
            / self.spacing
    }

    /// Betragsmäßig größte Quellstärke aller Zellen.
    pub fn max_divergence(&self) -> T {
        (0..self.nx * self.ny)
            .fold(T::ZERO, |worst, index| worst.max(self.divergence(index % self.nx, index / self.nx).abs()))
    }

    /// Wirbelstärke in der Mitte der Zelle: positiv dreht gegen den
    /// Uhrzeigersinn.
    pub fn vorticity(&self, i: usize, j: usize) -> T {
        if !self.is_fluid(i, j) {
            return T::ZERO;
        }
        let (left, right) = (i.saturating_sub(1), (i + 1).min(self.nx - 1));
        let (below, above) = (j.saturating_sub(1), (j + 1).min(self.ny - 1));
        let width = T::from_f64((right - left).max(1) as f64) * self.spacing;
        let height = T::from_f64((above - below).max(1) as f64) * self.spacing;
        (self.velocity(right, j)[1] - self.velocity(left, j)[1]) / width
            - (self.velocity(i, above)[0] - self.velocity(i, below)[0]) / height
    }

    /// `true`, wenn der Wert auf dieser Kante vorgegeben ist: Wand, Einlass
    /// oder Kante eines Hindernisses.
    fn fixed(&self, component: Component, i: usize, j: usize) -> bool {
        let open = |edge: Edge<T>| edge == Edge::Open;
        match component {
            Component::U => {
                if i == 0 {
                    !open(self.edges.left) || !self.is_fluid(0, j)
                } else if i >= self.nx {
                    !open(self.edges.right) || !self.is_fluid(self.nx - 1, j)
                } else {
                    !self.is_fluid(i - 1, j) || !self.is_fluid(i, j)
                }
            }
            Component::V => {
                if j == 0 {
                    !open(self.edges.bottom) || !self.is_fluid(i, 0)
                } else if j >= self.ny {
                    !open(self.edges.top) || !self.is_fluid(i, self.ny - 1)
                } else {
                    !self.is_fluid(i, j - 1) || !self.is_fluid(i, j)
                }
            }
        }
    }

    /// Setzt die vorgegebenen Geschwindigkeiten: null an Wänden und
    /// Hindernissen, den Einlasswert am Einlass.
    fn enforce_boundaries(&mut self) {
        let inflow = |edge: Edge<T>, sign: T| match edge {
            Edge::Inflow(speed) => speed * sign,
            _ => T::ZERO,
        };
        for j in 0..self.ny {
            for i in 0..=self.nx {
                if !self.fixed(Component::U, i, j) {
                    continue;
                }
                let value = if i == 0 && self.is_fluid(0, j) {
                    inflow(self.edges.left, T::ONE)
                } else if i == self.nx && self.is_fluid(self.nx - 1, j) {
                    inflow(self.edges.right, -T::ONE)
                } else {
                    T::ZERO
                };
                self.u[j * (self.nx + 1) + i] = value;
            }
        }
        for j in 0..=self.ny {
            for i in 0..self.nx {
                if !self.fixed(Component::V, i, j) {
                    continue;
                }
                let value = if j == 0 && self.is_fluid(i, 0) {
                    inflow(self.edges.bottom, T::ONE)
                } else if j == self.ny && self.is_fluid(i, self.ny - 1) {
                    inflow(self.edges.top, -T::ONE)
                } else {
                    T::ZERO
                };
                self.v[j * self.nx + i] = value;
            }
        }
    }

    /// Rechnet um `dt` weiter.
    ///
    /// Reihenfolge: Auftrieb, Transport der Geschwindigkeit, Zähigkeit,
    /// Druckkorrektur, Transport des Rauchs. Ein nicht positiver oder nicht
    /// endlicher Zeitschritt ändert nichts.
    pub fn step(&mut self, dt: T) -> GridStats<T> {
        if dt <= T::ZERO || !dt.is_finite() {
            return GridStats { iterations: 0, residual: T::ZERO, converged: true };
        }
        self.apply_buoyancy(dt);
        self.enforce_boundaries();
        self.advect_velocity(dt);
        self.enforce_boundaries();
        if self.params.viscosity > T::ZERO {
            self.diffuse(Component::U, dt);
            self.diffuse(Component::V, dt);
        }
        let stats = self.project();
        self.advect_smoke(dt);
        stats
    }

    fn apply_buoyancy(&mut self, dt: T) {
        if self.params.buoyancy == T::ZERO {
            return;
        }
        let factor = dt * self.params.buoyancy * T::from_f64(0.5);
        for j in 1..self.ny {
            for i in 0..self.nx {
                if !self.fixed(Component::V, i, j) {
                    let lift = self.smoke[(j - 1) * self.nx + i] + self.smoke[j * self.nx + i];
                    self.v[j * self.nx + i] += factor * lift;
                }
            }
        }
    }

    fn advect_velocity(&mut self, dt: T) {
        let flow = Flow { u: &self.u, v: &self.v, nx: self.nx, ny: self.ny, step: dt / self.spacing };
        let maccormack = self.params.maccormack;
        flow.advect(Layout::u(self.nx, self.ny), &self.u, &mut self.next_u, &mut self.scratch, maccormack);
        flow.advect(Layout::v(self.nx, self.ny), &self.v, &mut self.next_v, &mut self.scratch, maccormack);
        core::mem::swap(&mut self.u, &mut self.next_u);
        core::mem::swap(&mut self.v, &mut self.next_v);
    }

    fn advect_smoke(&mut self, dt: T) {
        let flow = Flow { u: &self.u, v: &self.v, nx: self.nx, ny: self.ny, step: dt / self.spacing };
        flow.advect(
            Layout::centre(self.nx, self.ny),
            &self.smoke,
            &mut self.next_smoke,
            &mut self.scratch,
            self.params.maccormack,
        );
        core::mem::swap(&mut self.smoke, &mut self.next_smoke);
        let keep = (-(self.params.dissipation.max(T::ZERO) * dt)).exp();
        for (value, &solid) in self.smoke.iter_mut().zip(&self.solid) {
            *value = if solid { T::ZERO } else { *value * keep };
        }
    }

    /// Zähigkeit, implizit: löst `(1 − ν·dt·Δ) q = q₀` näherungsweise mit
    /// Gauß-Seidel-Durchläufen. An Wänden haftet die Flüssigkeit.
    fn diffuse(&mut self, component: Component, dt: T) {
        let (width, height) = match component {
            Component::U => (self.nx + 1, self.ny),
            Component::V => (self.nx, self.ny + 1),
        };
        let amount = self.params.viscosity * dt / (self.spacing * self.spacing);
        let scale = (T::ONE + T::from_f64(4.0) * amount).recip();
        // Jeder Durchlauf verkleinert den Fehler mindestens um den Faktor
        // `rate`; so viele Durchläufe, bis er auf ein Zehntausendstel
        // gefallen ist.
        let rate = T::from_f64(4.0) * amount * scale;
        let sweeps = if rate > T::ZERO && rate < T::ONE {
            (T::from_f64(1.0e-4).ln() / rate.ln()).ceil().to_f64().clamp(2.0, 40.0)
        } else {
            40.0
        };
        let sweeps = sweeps as usize;
        let mut field = core::mem::take(match component {
            Component::U => &mut self.u,
            Component::V => &mut self.v,
        });
        self.initial[..field.len()].copy_from_slice(&field);
        for index in 0..field.len() {
            self.free[index] = !self.fixed(component, index % width, index / width);
        }
        // Quer zur Komponente liegende Ränder: Dort haftet die Strömung,
        // außer der Rand ist offen.
        let (low_open, high_open) = match component {
            Component::U => (self.edges.bottom == Edge::Open, self.edges.top == Edge::Open),
            Component::V => (self.edges.left == Edge::Open, self.edges.right == Edge::Open),
        };
        let ghost = |value: T, open: bool| if open { value } else { -value };
        for _ in 0..sweeps {
            for j in 0..height {
                for i in 0..width {
                    let index = j * width + i;
                    if !self.free[index] {
                        continue;
                    }
                    let centre = field[index];
                    let left = if i > 0 { field[index - 1] } else { centre };
                    let right = if i + 1 < width { field[index + 1] } else { centre };
                    let below = if j > 0 { field[index - width] } else { centre };
                    let above = if j + 1 < height { field[index + width] } else { centre };
                    let (left, right, below, above) = match component {
                        Component::U => (
                            left,
                            right,
                            if j > 0 { below } else { ghost(centre, low_open) },
                            if j + 1 < height { above } else { ghost(centre, high_open) },
                        ),
                        Component::V => (
                            if i > 0 { left } else { ghost(centre, low_open) },
                            if i + 1 < width { right } else { ghost(centre, high_open) },
                            below,
                            above,
                        ),
                    };
                    field[index] = (self.initial[index] + amount * (left + right + below + above)) * scale;
                }
            }
        }
        match component {
            Component::U => self.u = field,
            Component::V => self.v = field,
        }
    }

    /// Druckkorrektur: macht das Geschwindigkeitsfeld quellenfrei.
    fn project(&mut self) -> GridStats<T> {
        let (nx, ny) = (self.nx, self.ny);
        let open = [
            self.edges.left == Edge::Open,
            self.edges.right == Edge::Open,
            self.edges.bottom == Edge::Open,
            self.edges.top == Edge::Open,
        ];
        let mut sum = T::ZERO;
        let mut count = 0_usize;
        for index in 0..nx * ny {
            let (i, j) = (index % nx, index / nx);
            self.rhs[index] = if self.solid[index] {
                T::ZERO
            } else {
                count += 1;
                -(self.u[j * (nx + 1) + i + 1] - self.u[j * (nx + 1) + i] + self.v[(j + 1) * nx + i]
                    - self.v[j * nx + i])
            };
            sum += self.rhs[index];
        }
        // Ohne offenen Rand ist der Druck nur bis auf eine Konstante
        // bestimmt und die Gleichung nur lösbar, wenn sich die rechte Seite
        // zu null summiert. Das stellt der Abzug des Mittelwerts sicher.
        if !open.contains(&true) && count > 0 {
            let mean = sum / T::from_f64(count as f64);
            for (value, &solid) in self.rhs.iter_mut().zip(&self.solid) {
                if !solid {
                    *value -= mean;
                }
            }
        }

        if self.stale {
            self.stencil.rebuild(&self.solid, open);
            self.stale = false;
        }
        let tolerance = Tolerance::default()
            .with_rel_x(self.params.pressure_tolerance)
            .with_max_iterations(self.params.pressure_iterations);
        let result = preconditioned_conjugate_gradient(
            &PressureOperator(&self.stencil),
            &Preconditioner(&self.stencil),
            &self.rhs,
            &mut self.pressure,
            &mut self.work,
            &tolerance,
        );
        let stats = match result {
            Ok(state) => GridStats { iterations: state.iterations, residual: state.residual, converged: true },
            Err(SolveError::NotConverged(state)) => {
                GridStats { iterations: state.iterations, residual: state.residual, converged: false }
            }
            Err(_) => {
                // Kein brauchbarer Druck: Feld unverändert lassen und beim
                // nächsten Schritt neu beginnen.
                self.pressure.fill(T::ZERO);
                return GridStats { iterations: 0, residual: T::INFINITY, converged: false };
            }
        };

        for j in 0..ny {
            for i in 0..=nx {
                if self.fixed(Component::U, i, j) {
                    continue;
                }
                let left = if i > 0 { self.pressure[j * nx + i - 1] } else { T::ZERO };
                let right = if i < nx { self.pressure[j * nx + i] } else { T::ZERO };
                self.u[j * (nx + 1) + i] -= right - left;
            }
        }
        for j in 0..=ny {
            for i in 0..nx {
                if self.fixed(Component::V, i, j) {
                    continue;
                }
                let below = if j > 0 { self.pressure[(j - 1) * nx + i] } else { T::ZERO };
                let above = if j < ny { self.pressure[j * nx + i] } else { T::ZERO };
                self.v[j * nx + i] -= above - below;
            }
        }
        stats
    }
}
