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
}
