pub trait Estimate {
    fn generate_estimates(&self) -> Vec<(&[u8], f32)>;
}
