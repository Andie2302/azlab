//! Versionierte Bewegungsdateien der beiden Diagnosebeispiele.
//! Geometrie und Parameter werden über die feste Beispielvariante rekonstruiert.
use azviz::VizError;
use core::fmt::Write as _;
use std::io::Read as _;
use std::path::Path;

pub struct Snapshot {
    pub kind: String,
    pub variant: u32,
    pub steps: u64,
    pub time: f64,
    pub positions: Vec<[f32; 3]>,
    pub velocities: Vec<[f32; 3]>,
}

fn invalid() -> VizError {
    VizError::Encoding("ungültige oder inkompatible Zustandsdatei".into())
}

impl Snapshot {
    pub fn encode(&self) -> Result<String, VizError> {
        if self.positions.len() != self.velocities.len() || self.positions.len() > 100_000
            || !self.time.is_finite() || self.time < 0.0
            || self.kind.split_whitespace().count() != 1
            || self.positions.iter().chain(&self.velocities).flatten().any(|v| !v.is_finite()) {
            return Err(invalid());
        }
        let mut text = format!("AZLAB_MOTION_V1 {} {} {} {:016x} {}\n",
            self.kind, self.variant, self.steps, self.time.to_bits(), self.positions.len());
        for (p, v) in self.positions.iter().zip(&self.velocities) {
            let _ = writeln!(text, "{:08x} {:08x} {:08x} {:08x} {:08x} {:08x}",
                p[0].to_bits(), p[1].to_bits(), p[2].to_bits(),
                v[0].to_bits(), v[1].to_bits(), v[2].to_bits());
        }
        Ok(text)
    }

    pub fn decode(text: &str) -> Result<Self, VizError> {
        let mut words = text.split_whitespace();
        if words.next() != Some("AZLAB_MOTION_V1") { return Err(invalid()); }
        let kind = words.next().ok_or_else(invalid)?.to_owned();
        let variant = words.next().ok_or_else(invalid)?.parse().map_err(|_| invalid())?;
        let steps = words.next().ok_or_else(invalid)?.parse().map_err(|_| invalid())?;
        let time = f64::from_bits(u64::from_str_radix(words.next().ok_or_else(invalid)?, 16).map_err(|_| invalid())?);
        let count: usize = words.next().ok_or_else(invalid)?.parse().map_err(|_| invalid())?;
        if count > 100_000 || !time.is_finite() || time < 0.0 { return Err(invalid()); }
        let mut positions = Vec::with_capacity(count);
        let mut velocities = Vec::with_capacity(count);
        for _ in 0..count {
            let mut data = [0.0f32; 6];
            for value in &mut data {
                *value = f32::from_bits(u32::from_str_radix(words.next().ok_or_else(invalid)?, 16).map_err(|_| invalid())?);
                if !value.is_finite() { return Err(invalid()); }
            }
            positions.push([data[0], data[1], data[2]]);
            velocities.push([data[3], data[4], data[5]]);
        }
        if words.next().is_some() { return Err(invalid()); }
        Ok(Self { kind, variant, steps, time, positions, velocities })
    }

    pub fn save(&self, path: &Path) -> Result<(), VizError> {
        let text = self.encode()?;
        std::fs::write(path, text)?;
        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self, VizError> {
        let mut text = String::new();
        std::fs::File::open(path)?.take(16 * 1024 * 1024 + 1).read_to_string(&mut text)?;
        if text.len() > 16 * 1024 * 1024 { return Err(invalid()); }
        Self::decode(&text)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_float_roundtrip_and_invalid_input() {
        let state = Snapshot { kind: "test".into(), variant: 2, steps: 42, time: 0.7,
            positions: vec![[f32::from_bits(1), -0.0, 1.25]],
            velocities: vec![[0.1, -2.0, 3.0]] };
        let encoded = state.encode().expect("kodiert");
        let decoded = Snapshot::decode(&encoded).expect("dekodiert");
        assert_eq!(decoded.steps, 42);
        assert_eq!(decoded.time.to_bits(), state.time.to_bits());
        for (a, b) in decoded.positions.iter().flatten().zip(state.positions.iter().flatten()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(decoded.velocities, state.velocities);
        assert!(Snapshot::decode(&format!("{encoded} extra")).is_err());
        assert!(Snapshot::decode("AZLAB_MOTION_V1 test 0 0 0 1").is_err());
        assert!(Snapshot::decode("AZLAB_MOTION_V1 test 0 0 0 1 7fc00000 0 0 0 0 0").is_err());
    }
}
