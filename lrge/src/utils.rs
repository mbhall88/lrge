use anyhow::Result;
use std::path::PathBuf;

pub(crate) fn create_temp_dir(temp_dir: Option<&PathBuf>, keep: bool) -> Result<tempfile::TempDir> {
    let mut binding = tempfile::Builder::new();
    let builder = binding.keep(keep).prefix("lrge-");
    let tmpdir = match temp_dir {
        Some(path) => {
            if !path.exists() {
                std::fs::create_dir_all(path)?;
            }
            builder.tempdir_in(path)?
        }
        None => builder.tempdir()?,
    };
    Ok(tmpdir)
}

pub(crate) fn format_estimate(estimate: f32) -> String {
    if estimate.is_infinite() {
        return String::from("∞ bp");
    }

    // Define the metric suffixes and their corresponding powers of 10
    let units = [
        ("bp", 0),
        ("kbp", 1),
        ("Mbp", 2),
        ("Gbp", 3),
        ("Tbp", 4),
        ("Pbp", 5),
    ];

    // Determine the appropriate unit
    let mut value = estimate;
    let mut suffix = "bp";
    for (unit, power) in units {
        let threshold = 10f32.powi(power * 3); // 10^(power * 3) for 10^0, 10^3, etc.
        if estimate >= threshold {
            value = estimate / threshold;
            suffix = unit;
        } else {
            break;
        }
    }

    // Format the value with the determined suffix
    format!("{:.2} {}", value, suffix)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_create_temp_dir_with_none() {
        let result = create_temp_dir(None, false);
        assert!(result.is_ok());
        let temp_dir = result.unwrap();
        assert!(temp_dir.path().exists());
    }

    #[test]
    fn test_create_temp_dir_with_some() {
        let base_dir = TempDir::new().unwrap();
        let result = create_temp_dir(Some(&base_dir.path().to_path_buf()), false);
        assert!(result.is_ok());
        let temp_dir = result.unwrap();
        assert!(temp_dir.path().exists());
        assert!(temp_dir.path().starts_with(base_dir.path()));
    }

    #[test]
    fn test_create_temp_dir_with_keep() {
        let result = create_temp_dir(None, true);
        assert!(result.is_ok());
        let temp_dir = result.unwrap();
        let temp_path = temp_dir.path().to_path_buf();
        drop(temp_dir); // Dropping the TempDir to check if it is kept
        assert!(temp_path.exists());
        fs::remove_dir_all(temp_path).unwrap(); // Clean up
    }

    #[test]
    fn test_create_temp_dur_with_some_that_does_not_exist() {
        let base_dir = TempDir::new().unwrap();
        let non_existent_dir = base_dir.path().join("non_existent");
        let result = create_temp_dir(Some(&non_existent_dir), false);
        assert!(result.is_ok());
        let temp_dir = result.unwrap();
        assert!(temp_dir.path().exists());
        assert!(temp_dir.path().starts_with(non_existent_dir));
    }

    #[test]
    fn test_bp_range() {
        assert_eq!(format_estimate(0.0), "0.00 bp");
        assert_eq!(format_estimate(999.99), "999.99 bp");
    }

    #[test]
    fn test_kbp_range() {
        assert_eq!(format_estimate(1_000.0), "1.00 kbp");
        assert_eq!(format_estimate(1_234.56), "1.23 kbp");
        assert_eq!(format_estimate(999_999.99), "1.00 Mbp");
    }

    #[test]
    fn test_mbp_range() {
        assert_eq!(format_estimate(1_000_000.0), "1.00 Mbp");
        assert_eq!(format_estimate(1_500_000.0), "1.50 Mbp");
        assert_eq!(format_estimate(999_999_999.99), "1.00 Gbp");
    }

    #[test]
    fn test_gbp_range() {
        assert_eq!(format_estimate(1_000_000_000.0), "1.00 Gbp");
        assert_eq!(format_estimate(1_500_000_000.0), "1.50 Gbp");
        assert_eq!(format_estimate(999_999_999_999.99), "1.00 Tbp");
    }

    #[test]
    fn test_tbp_range() {
        assert_eq!(format_estimate(1_000_000_000_000.0), "1.00 Tbp");
        assert_eq!(format_estimate(1_500_000_000_000.0), "1.50 Tbp");
        assert_eq!(format_estimate(999_999_999_999_999.99), "1.00 Pbp");
    }

    #[test]
    fn test_pbp_range() {
        assert_eq!(format_estimate(1_000_000_000_000_000.0), "1.00 Pbp");
        assert_eq!(format_estimate(4_500_000_000_000_000.0), "4.50 Pbp");
    }

    #[test]
    fn test_infinity() {
        assert_eq!(format_estimate(f32::INFINITY), "∞ bp");
    }

    #[test]
    fn test_small_values() {
        assert_eq!(format_estimate(0.1), "0.10 bp");
        assert_eq!(format_estimate(10.0), "10.00 bp");
        assert_eq!(format_estimate(999.99), "999.99 bp");
    }
}
