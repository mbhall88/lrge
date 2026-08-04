"""This script is use to replace existing estimates for lrge 2set to save having to 
rerun the entire pipeline. It is not intended to be used for any other purpose.
We discovered that using the median of the finite estimates was a better approach"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import math


# function to apply to each group. it produces the finite median for each group
def f(group):
    estimates = group["estimate"]
    median = estimates.median()
    finite = estimates[~np.isinf(estimates)]
    finite_median = finite.median()
    return pd.Series([median, finite_median])


def main():
    searchdir = Path(sys.argv[2])
    df = pd.read_csv(sys.argv[1], low_memory=False)
    # we only want to replace the 2set estimates
    df.query("method == '2set'", inplace=True)

    groups = df.groupby("accession")
    medians_df = groups.apply(f)
    medians_df.columns = ["median", "finite_median"]

    for acc, row in medians_df.iterrows():
        # replace the 2set estimates with the finite median
        median = row["median"]
        finite_median = row["finite_median"]
        dir1 = acc[:3]
        dir2 = acc[3:6]
        dir3 = acc[6:9]
        path_to_estimate_file = searchdir / f"{dir1}/{dir2}/{dir3}/{acc}/{acc}.size"
        assert path_to_estimate_file.exists(), f"{path_to_estimate_file} does not exist"
        estimate_in_file = float(path_to_estimate_file.read_text().strip())
        if math.isclose(estimate_in_file, finite_median, rel_tol=1e-5):
            print(f"Skipping {acc} as estimate is already {finite_median}")
            continue
        elif not math.isclose(estimate_in_file, median, rel_tol=1e-5):
            raise ValueError(
                f"Estimate in file {path_to_estimate_file} ({estimate_in_file}) is not the same as the median {median}"
            )
        path_to_estimate_file.write_text(str(finite_median))


if __name__ == "__main__":
    main()
