mod constants;
mod ground_track;
mod helpers;
mod renderer;
mod sgp4;
mod test_constants;
mod types;

fn main() -> eframe::Result {
    let renderer = renderer::Renderer::new(helpers::compute_points);

    renderer.run()
}
