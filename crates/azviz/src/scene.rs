//! Öffentliche Schnittstelle zwischen Simulation und Anzeige.

use alloc::string::String;
use crate::{Canvas, VizError};

/// Eine Simulation, die sich zeigen lässt.
pub trait Scene {
    /// Titel des Fensters.
    fn title(&self) -> String;

    /// Rechnet um `dt` Sekunden weiter.
    fn update(&mut self, dt: f32);

    /// Aktiviert zusätzliche Messungen, wenn die Szene sie unterstützt.
    fn set_profiling(&mut self, _enabled: bool) {}

    /// Liefert und leert zusätzliche Messwerte nach einem Messfenster.
    fn profile_details(&mut self) -> Option<String> { None }

    /// Ändert die Ansicht, ohne die Simulation fortzuschreiben.
    /// Standardmäßig ohne Wirkung (etwa bei ebenen Flüssigkeitsbeispielen).
    fn orbit(&mut self, _yaw: f32, _pitch: f32, _zoom: f32) {}

    /// Blendet in räumlichen Beispielen den Boden ein oder aus.
    fn toggle_ground(&mut self) {}

    /// Zustandsdatei für S/L; nur Szenen mit expliziter Unterstützung.
    fn state_path(&self) -> Option<&str> { None }

    /// Speichert einen Zustand. Standardmäßig nicht unterstützt.
    fn save_state(&self, _path: &std::path::Path) -> Result<(), VizError> {
        Err(VizError::Arguments("Szene unterstützt keine Zustandsdateien".into()))
    }

    /// Lädt einen Zustand; Implementierungen müssen Fehler atomar behandeln.
    fn load_state(&mut self, _path: &std::path::Path) -> Result<(), VizError> {
        Err(VizError::Arguments("Szene unterstützt keine Zustandsdateien".into()))
    }

    /// Wählt eine Diagnosevariante (Tasten 1 bis 6).
    fn select_variant(&mut self, _variant: u32) {}

    /// Zeichnet den aktuellen Zustand.
    fn draw(&self, canvas: &mut dyn Canvas);
}

