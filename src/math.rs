const INV_SQRT_PI: f64 = 0.3989422804014327;

pub fn normal_pdf(x: f32, mu: f32, sigma: f32) -> f32 {
    let a = (x - mu) / sigma;

    INV_SQRT_PI as f32 / sigma * (-0.5 * a * a).exp()
}
