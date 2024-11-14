//! A trait for generating genome size estimates, and calculating the median of those estimates.

/// This trait provides methods to generate estimates and calculate the median
/// of those estimates, both with and without considering infinite values.
pub trait Estimate {
    /// Generate a list of genome size estimates.
    ///
    /// # Returns
    ///
    /// A `Vec<f32>` containing the generated estimates. These estimates may be finite or infinite.
    fn generate_estimates(&self) -> Vec<f32>;

    /// Generate an estimate of the genome size, taking the median of the finite estimates.
    ///
    /// Note that this method will return `None` if there are no finite estimates.
    fn estimate(&self) -> Option<f32> {
        let estimates = self.generate_estimates();
        let iter = estimates.iter().filter(|&x| x.is_finite()).copied();
        median(iter)
    }

    /// Generate an estimate of the genome size, taking the median of all estimates - infinity
    /// included.
    ///
    /// Note that this method will return `None` if there are no estimates.
    fn estimate_with_infinity(&self) -> Option<f32> {
        let estimates = self.generate_estimates();
        let iter = estimates.iter().copied();
        median(iter)
    }
}

fn median(iter: impl Iterator<Item = f32>) -> Option<f32> {
    let mut values: Vec<f32> = iter.collect();
    let len = values.len();

    if len == 0 {
        return None;
    }

    // we unwrap here because we don't expect NaN values
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if len % 2 == 1 {
        Some(values[len / 2])
    } else {
        let mid1 = values[len / 2 - 1];
        let mid2 = values[len / 2];
        Some((mid1 + mid2) / 2.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd_length() {
        let data = vec![1.0f32, 3.0, 5.0, 7.0, 9.0];
        assert_eq!(median(data.into_iter()), Some(5.0));
    }

    #[test]
    fn test_median_unsorted() {
        let data = vec![3.0f32, 1.0, 7.0, 5.0, 9.0];
        assert_eq!(median(data.into_iter()), Some(5.0));
    }

    #[test]
    fn test_median_even_length() {
        let data = vec![1.0f32, 3.0, 5.0, 7.0];
        assert_eq!(median(data.into_iter()), Some(4.0));
    }

    #[test]
    fn test_median_single_element() {
        let data = vec![10.0f32];
        assert_eq!(median(data.into_iter()), Some(10.0));
    }

    #[test]
    fn test_median_empty() {
        let data: Vec<f32> = vec![];
        assert_eq!(median(data.into_iter()), None);
    }

    #[test]
    fn test_median_with_negative_numbers() {
        let data = vec![-3.0f32, 1.0, 0.0, 3.0, -1.0];
        assert_eq!(median(data.into_iter()), Some(0.0));
    }

    #[test]
    fn test_median_with_positive_infinity() {
        let data = vec![1.0f32, 2.0, 3.0, f32::INFINITY];
        assert_eq!(median(data.into_iter()), Some(2.5));
    }

    #[test]
    fn test_median_with_negative_infinity() {
        let data = vec![f32::NEG_INFINITY, 1.0, 2.0, 3.0];
        assert_eq!(median(data.into_iter()), Some(1.5));
    }

    #[test]
    fn test_median_with_both_infinities() {
        let data = vec![f32::NEG_INFINITY, 1.0, 2.0, f32::INFINITY];
        assert_eq!(median(data.into_iter()), Some(1.5));
    }

    #[test]
    fn test_median_with_only_infinity() {
        let data = vec![f32::INFINITY, f32::INFINITY];
        assert_eq!(median(data.into_iter()), Some(f32::INFINITY));
    }

    #[test]
    fn test_median_with_only_negative_infinity() {
        let data = vec![f32::NEG_INFINITY, f32::NEG_INFINITY];
        assert_eq!(median(data.into_iter()), Some(f32::NEG_INFINITY));
    }

    #[test]
    fn test_median_with_inf_and_regular_values() {
        let data = vec![-1.0, f32::NEG_INFINITY, 0.0, 1.0, f32::INFINITY];
        assert_eq!(median(data.into_iter()), Some(0.0));
    }
}
