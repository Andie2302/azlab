//! Farben.

/// Farbe aus Rot, Grün, Blau und Deckkraft, je 0 bis 255.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color {
    /// Rot.
    pub r: u8,
    /// Grün.
    pub g: u8,
    /// Blau.
    pub b: u8,
    /// Deckkraft; 255 ist undurchsichtig.
    pub a: u8,
}

impl Color {
    /// Dunkler Hintergrund.
    pub const BACKGROUND: Self = Self::rgb(18, 22, 30);
    /// Weiß.
    pub const WHITE: Self = Self::rgb(240, 240, 240);
    /// Grau.
    pub const GRAY: Self = Self::rgb(110, 118, 130);
    /// Rot.
    pub const RED: Self = Self::rgb(226, 86, 78);
    /// Orange.
    pub const ORANGE: Self = Self::rgb(240, 160, 60);
    /// Grün.
    pub const GREEN: Self = Self::rgb(96, 190, 120);
    /// Blau.
    pub const BLUE: Self = Self::rgb(80, 150, 235);

    /// Undurchsichtige Farbe.
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: 255 }
    }

    /// Dieselbe Farbe mit anderer Deckkraft.
    pub const fn with_alpha(self, a: u8) -> Self {
        Self { a, ..self }
    }

    /// Mischung zweier Farben; `t = 0` ergibt `self`, `t = 1` ergibt `other`.
    pub fn mix(self, other: Self, t: f32) -> Self {
        let t = t.clamp(0.0, 1.0);
        let blend = |from: u8, to: u8| (f32::from(from) + (f32::from(to) - f32::from(from)) * t) as u8;
        Self {
            r: blend(self.r, other.r),
            g: blend(self.g, other.g),
            b: blend(self.b, other.b),
            a: blend(self.a, other.a),
        }
    }

    /// Dieselbe Farbe, um den Faktor `factor` (0 bis 1) abgedunkelt.
    pub fn shade(self, factor: f32) -> Self {
        let factor = factor.clamp(0.0, 1.0);
        let scale = |channel: u8| (f32::from(channel) * factor) as u8;
        Self { r: scale(self.r), g: scale(self.g), b: scale(self.b), a: self.a }
    }
}
