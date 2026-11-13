//! Prüft die öffentliche API so, wie sie von azphysics-Beispielen benutzt wird.
//! Fehlende Trait-Methoden müssen hier als Compilerfehler auffallen.
use azviz::{Camera, Canvas, Scene, VizError};
use std::path::Path;

struct ExampleScene {
    camera: Camera,
    profiling: bool,
    ground: bool,
    variant: u32,
}

impl Scene for ExampleScene {
    fn title(&self) -> String { "API-Test".into() }
    fn update(&mut self, _dt: f32) {}
    fn draw(&self, _canvas: &mut dyn Canvas) {}
    fn set_profiling(&mut self, enabled: bool) { self.profiling = enabled; }
    fn profile_details(&mut self) -> Option<String> { self.profiling.then(|| "Profil".into()) }
    fn orbit(&mut self, yaw: f32, pitch: f32, zoom: f32) { self.camera.orbit(yaw, pitch, zoom); }
    fn toggle_ground(&mut self) { self.ground = !self.ground; }
    fn state_path(&self) -> Option<&str> { Some("api-test.azstate") }
    fn save_state(&self, _path: &Path) -> Result<(), VizError> { Ok(()) }
    fn load_state(&mut self, _path: &Path) -> Result<(), VizError> { Ok(()) }
    fn select_variant(&mut self, variant: u32) { self.variant = variant; }
}

#[test]
fn public_scene_api_matches_examples() {
    let mut example = ExampleScene {
        camera: Camera::look_at([0.0, 1.0, 3.0], [0.0, 0.0, 0.0], 40.0),
        profiling: false, ground: true, variant: 0,
    };
    let scene: &mut dyn Scene = &mut example;
    scene.set_profiling(true);
    assert_eq!(scene.profile_details().as_deref(), Some("Profil"));
    scene.orbit(0.1, 0.1, 0.1);
    scene.toggle_ground();
    scene.select_variant(2);
    assert_eq!(scene.state_path(), Some("api-test.azstate"));
    assert!(scene.save_state(Path::new("unused")).is_ok());
    assert!(scene.load_state(Path::new("unused")).is_ok());
    assert!(example.profiling && !example.ground && example.variant == 2);
}
