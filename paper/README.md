This directory contains all the code and data used in [the paper][doi]. The code is organised in the following way:

- [`workflow/`](./workflow): Contains the code to reproduce the analyses of the paper. It requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) to run.
- [`config/`](./config): Contains the configuration files for the workflow as well as the metadata for the samples. The final dataset used in the 
    paper is available is [`config/bacteria_lr_runs.filtered.tsv`](./config/bacteria_lr_runs.filtered.tsv). Various intermediate metadata files prior to filtering can be found 
    in [`config/`](./config). [`config/get_datasets.md`](./config/get_datasets.md) contains details how where the data were 
    downloaded from and how the metadata was generated and filtered.
- [`scripts/`](./scripts): Miscellaneous scripts used for the paper, but directly part of the workflow.
- [`notebooks/`](./notebooks): Jupyter notebooks used for the paper. These are not part of the workflow, but were used to generate figures and tables.
- [`results/`](./results): Contains the results of the workflow. The final estimates used in the paper are available in [`results/estimates/estimates.tsv`](./results/estimates/estimates.tsv). 
    The figures and tables for the paper are available in [`results/figures/`](./results/figures) and [`results/tables/`](./results/tables), respectively.

[doi]: https://doi.org/10.1101/2024.11.27.625777