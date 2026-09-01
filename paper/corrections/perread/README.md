# Per-read estimates

One file per (accession, variant). Each line is one query read's genome-size estimate as emitted at
`-vv`; `inf` means the read overlapped no target read and was excluded from the median. The
infinite-estimate counts in [`../NOTES.md`](../NOTES.md) are `grep -c inf` on these files.

All runs used `-T 10000 -Q 5000 -s 4556` unless the tag says otherwise. Flags below were recovered
from the `Args { ... }` line of each run's log (the logs themselves are not preserved).

## Tag vocabulary

| Tag | Binary | Flags |
|---|---|---|
| `default`, `F`, `F-strict`, `overhang-strict` | **pre-fix** (`lrge-0.3.0`, before #31) | see below |
| `clean-*` | post-fix, built from a worktree carrying only the #31 change | as below |
| `fixed-*` | post-fix, built from merged `main` | as below |
| `-default` | — | no `-F`, ratio 0.2 (default) |
| `-F` | — | `-F`, ratio 0.2 |
| `-F-strict` | — | `-F --max-overhang-ratio 0.05` |
| `-F-vloose` / `-F-loose` | — | `-F --max-overhang-ratio 0.5` / `1.0` |
| `overhang-strict` | pre-fix | `--max-overhang-ratio 0.05` **without** `-F` |
| `n<N>` | post-fix | default flags, input subsampled to `<N>` reads |

`clean-*` and `fixed-*` are both post-fix and agree wherever both exist (e.g. SRR8618952
`-default` = 5,890,324 and `-F` = 7,201,990 under both). The two prefixes exist only because the
first batch was run from an isolated worktree while another agent was editing the main tree.

## Notable files

- `SRR16767125.{default,F,F-strict,overhang-strict}` — the pre-fix baseline for #31. `F` = 197,672
  with 1,883 infinite estimates is the "`-F` makes it worse" evidence; `overhang-strict` is
  byte-identical to `default`, which is the no-op proof for defect 3.
- `SRR16767125.fixed-F-{vloose,loose}` — the post-fix threshold sweep (0.5, 1.0).
- `SRR8618952.n*` — the read-count control series (see NOTES.md §3.5).
- `SRR26465560.fixed-F*`, `SRR26465526.fixed-F*` — the ONT probe of §4.3, run with `-P ont`. All
  other files are `-P pb`.
- The 13 ONT accessions of §4.4/§4.5 were run with `-P ont` by `../ont_31.sh`, which also recorded
  each run's internal-match fraction; those figures live in `../rerun_estimates.tsv`, not here.
