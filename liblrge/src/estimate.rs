//! A trait for generating genome size estimates, and calculating the median of those estimates.

/// The lower quantile we found to give the highest confidence in our analysis.
pub const LOWER_QUANTILE: f32 = 0.15;
/// The upper quantile we found to give the highest confidence in our analysis.
pub const UPPER_QUANTILE: f32 = 0.65;

pub struct EstimateResult {
    /// The lower quantile of the estimates
    pub lower: Option<f32>,
    /// The median of the estimates - this is the genome size estimate
    pub estimate: Option<f32>,
    /// The upper quantile of the estimates
    pub upper: Option<f32>,
    /// The number of reads that did not have an overlap
    pub no_mapping_count: u32,
}

/// This trait provides methods to generate estimates and calculate the median
/// of those estimates, both with and without considering infinite values.
pub trait Estimate {
    /// Generate a list of genome size estimates.
    ///
    /// # Returns
    ///
    /// A `Vec<f32>` containing the generated estimates. These estimates may be finite or infinite.
    fn generate_estimates(&mut self) -> crate::Result<(Vec<f32>, u32)>;

    /// Generate an estimate of the genome size, taking the median of the per-read estimates.
    ///
    /// # Arguments
    ///
    /// * `finite`: Whether to consider only finite estimates. We found setting this to `true` gave
    ///   more accurate results (see [the paper][doi]).
    /// * `lower_quant`: The lower percentile to calculate. If `None`, this will not be calculated.
    ///   This value should be between 0 and 0.5. So, for the 25th percentile, you would pass `0.25`.
    /// * `upper_quant`: The upper percentile to calculate. If `None`, this will not be calculated.
    ///   This value should be between 0.5 and 1.0. So, for the 75th percentile, you would pass `0.75`.
    ///
    /// In [our analysis][doi], we found that the 15th and 65th percentiles gave the highest confidence (~92%).
    /// If you want to use our most current recommended values, you can use the constants [`LOWER_QUANTILE`]
    /// and [`UPPER_QUANTILE`]. You can of course use any values you like.
    ///
    /// # Returns
    ///
    /// An [`EstimateResult`] containing the lower, median, and upper estimates, as well as the number
    /// of reads that did not have an overlap. This number can be important for quality control. For
    /// examples, if you have a high number of reads that did not overlap, you may want to investigate
    /// why that is (e.g., contamination, poor quality reads, etc.).
    ///
    /// The estimate will be `None` if there are no finite estimates when `finite` is `true`, or if
    /// there are no estimates at all.
    ///
    /// [doi]: https://doi.org/10.1101/2024.11.27.625777
    fn estimate(
        &mut self,
        finite: bool,
        lower_quant: Option<f32>,
        upper_quant: Option<f32>,
    ) -> crate::Result<EstimateResult> {
        let (estimates, no_mapping_count) = self.generate_estimates()?;

        let iter: Box<dyn Iterator<Item = f32>> = if finite {
            Box::new(estimates.iter().filter(|&x| x.is_finite()).copied())
        } else {
            Box::new(estimates.iter().copied())
        };

        let (lower, median, upper) = median(iter, lower_quant, upper_quant);

        Ok(EstimateResult {
            lower,
            estimate: median,
            upper,
            no_mapping_count,
        })
    }
}

fn median(
    iter: impl Iterator<Item = f32>,
    lower_quant: Option<f32>,
    upper_quant: Option<f32>,
) -> (Option<f32>, Option<f32>, Option<f32>) {
    let mut values: Vec<f32> = iter.collect();
    let len = values.len();

    if len == 0 {
        return (None, None, None);
    }

    // we unwrap here because we don't expect NaN values
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut quantiles = vec![0.5];
    if let Some(lower) = lower_quant {
        quantiles.push(lower);
    }
    if let Some(upper) = upper_quant {
        quantiles.push(upper);
    }
    let quantiles: Vec<_> = quantiles
        .iter()
        .map(|&q| calculate_quantile(&values, q))
        .collect();
    match (lower_quant, upper_quant) {
        (Some(_), Some(_)) => (quantiles[1], quantiles[0], quantiles[2]),
        (Some(_), None) => (quantiles[1], quantiles[0], None),
        (None, Some(_)) => (quantiles[0], quantiles[1], quantiles[2]),
        (None, None) => (None, quantiles[0], None),
    }
}

fn calculate_quantile(data: &[f32], quantile: f32) -> Option<f32> {
    if data.is_empty() {
        return None;
    }
    if !(0.0..=1.0).contains(&quantile) {
        panic!("Quantile must be between 0.0 and 1.0");
    }

    let n = data.len();
    let pos = quantile * (n - 1) as f32;
    let idx = pos.floor() as usize;
    let frac = pos - idx as f32;

    if idx + 1 < n {
        Some(data[idx] * (1.0 - frac) + data[idx + 1] * frac)
    } else {
        Some(data[idx])
    }
}

/// Estimate genome size using the formula from Equation 3 in [the paper][doi].
///
/// # Returns
///
/// A floating point number representing the estimated genome size. If the number of overlaps is 0,
/// this function will return [`f32::INFINITY`].
///
/// [doi]: https://doi.org/10.1101/2024.11.27.625777
pub(crate) fn per_read_estimate(
    read_len: usize,
    avg_target_len: f32,
    n_target_reads: usize,
    n_ovlaps: usize,
    ovlap_thresh: u32,
) -> f32 {
    if n_ovlaps == 0 {
        return f32::INFINITY;
    }

    let ovlap_ratio: f32 = n_target_reads as f32 / n_ovlaps as f32;

    read_len as f32 + ovlap_ratio * (read_len as f32 + avg_target_len - 2.0 * ovlap_thresh as f32 + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd_length() {
        let data = vec![1.0f32, 3.0, 5.0, 7.0, 9.0];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(5.0), None)
        );
    }

    #[test]
    fn test_median_unsorted() {
        let data = vec![3.0f32, 1.0, 7.0, 5.0, 9.0];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(5.0), None)
        );
    }

    #[test]
    fn test_median_even_length() {
        let data = vec![1.0f32, 3.0, 5.0, 7.0];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(4.0), None)
        );
    }

    #[test]
    fn test_median_single_element() {
        let data = vec![10.0f32];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(10.0), None)
        );
    }

    #[test]
    fn test_median_empty() {
        let data: Vec<f32> = vec![];
        assert_eq!(median(data.into_iter(), None, None), (None, None, None));
    }

    #[test]
    fn test_median_with_negative_numbers() {
        let data = vec![-3.0f32, 1.0, 0.0, 3.0, -1.0];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(0.0), None)
        );
    }

    #[test]
    fn test_median_with_positive_infinity() {
        let data = vec![1.0f32, 2.0, 3.0, f32::INFINITY];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(2.5), None)
        );
    }

    #[test]
    fn test_median_with_negative_infinity() {
        let data = vec![f32::NEG_INFINITY, 1.0, 2.0, 3.0];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(1.5), None)
        );
    }

    #[test]
    fn test_median_with_both_infinities() {
        let data = vec![f32::NEG_INFINITY, 1.0, 2.0, f32::INFINITY];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(1.5), None)
        );
    }

    #[test]
    fn test_median_with_only_infinity() {
        let data = vec![f32::INFINITY, f32::INFINITY];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(f32::INFINITY), None)
        );
    }

    #[test]
    fn test_median_with_only_negative_infinity() {
        let data = vec![f32::NEG_INFINITY, f32::NEG_INFINITY];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(f32::NEG_INFINITY), None)
        );
    }

    #[test]
    fn test_median_with_inf_and_regular_values() {
        let data = vec![-1.0, f32::NEG_INFINITY, 0.0, 1.0, f32::INFINITY];
        assert_eq!(
            median(data.into_iter(), None, None),
            (None, Some(0.0), None)
        );
    }

    #[test]
    fn test_median_with_quantiles() {
        let data = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        assert_eq!(
            median(data.into_iter(), Some(0.15), Some(0.65)),
            (Some(2.35), Some(5.5), Some(6.85))
        );
    }

    #[test]
    fn test_calculate_quantile_with_infinity_in_quantile() {
        let data = vec![
            1.0f32,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            f32::INFINITY,
            f32::INFINITY,
            f32::INFINITY,
            f32::INFINITY,
        ];
        assert_eq!(
            median(data.into_iter(), Some(0.15), Some(0.65)),
            (Some(2.35), Some(5.5), Some(f32::INFINITY))
        );
    }

    #[test]
    #[should_panic(expected = "Quantile must be between 0.0 and 1.0")]
    fn test_calculate_quantile_panics() {
        let data = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
        calculate_quantile(&data, 1.1);
    }

    #[test]
    fn test_per_read_estimate() {
        let read_len = 100;
        let avg_target_len = 200.0;
        let n_target_reads = 1000;
        let n_ovlaps = 100;
        let ovlap_thresh = 10;
        let expected = 2800.0;
        assert_eq!(
            per_read_estimate(
                read_len,
                avg_target_len,
                n_target_reads,
                n_ovlaps,
                ovlap_thresh
            ),
            expected
        );
    }

    #[test]
    fn test_per_read_estimate_zero_ovlaps() {
        let read_len = 100;
        let avg_target_len = 200.0;
        let n_target_reads = 1000;
        let n_ovlaps = 0;
        let ovlap_thresh = 10;
        let expected = f32::INFINITY;
        assert_eq!(
            per_read_estimate(
                read_len,
                avg_target_len,
                n_target_reads,
                n_ovlaps,
                ovlap_thresh
            ),
            expected
        );
    }
}
