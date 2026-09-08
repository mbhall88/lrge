# Changelog

## [2.0.0](https://github.com/mbhall88/lrge/compare/lrge-1.0.0...lrge-2.0.0) (2026-09-08)


### ⚠ BREAKING CHANGES

* PacBio estimates change, because -P pb now selects the ava-pb overlap preset instead of being ignored. The reported interval changes on both platforms, because its quantiles are refitted and now depend on -P. liblrge::estimate::LOWER_QUANTILE and UPPER_QUANTILE keep their values but are deprecated in favour of Platform::interval_quantiles. The estimate line prints 95% CI or a named quantile pair where it used to print IQR. Anyone reproducing the published results should pin the 0.3.x series.
* **liblrge:** -F/--filter-contained no longer excludes every internal match. It now excludes them only when they exceed the share given to it, which defaults to 0.8, and that share must be written after an equals sign, as -F=0.5. Pass -F=0 for the previous behaviour. In liblrge, Builder::remove_internal is deprecated in favour of Builder::internal_filter, which takes an optional threshold rather than a flag; remove_internal still compiles, and asking it to filter is exactly a threshold of zero. Separately, an input too small to fill both read sets is now divided between them in the requested ratio rather than taking the whole shortfall from the target set, per #62, so every input holding fewer reads than -T and -Q ask for gets a different estimate from 0.3.0 on default flags; use --shortfall target for the previous behaviour.
* `lrge --max-overhang-ratio <FLOAT>` without `-F/--filter-contained` now exits with an error instead of silently ignoring the value, and `-F` filters the opposite set of alignments to before.

### Features

* a more accurate calculation of per read estimation ([5d2e211](https://github.com/mbhall88/lrge/commit/5d2e211e976ab7df9d3db05243175c9fda210b90))
* add option '--use-min-ref' for using smaller reference ([4791ec0](https://github.com/mbhall88/lrge/commit/4791ec066dd204c97cd91884b473c3bcf50e2931))
* add option '-F' to remove internal matches ([b23a8c7](https://github.com/mbhall88/lrge/commit/b23a8c7fea82aff1344c4f108628cf6044245245))
* add quantile logging and options ([7f7ad53](https://github.com/mbhall88/lrge/commit/7f7ad53856fa67cafaad03136e92ff5cc675c364))
* **bin:** add long threads options ([9ac48ec](https://github.com/mbhall88/lrge/commit/9ac48eca5e1ced2746c3697d4b68faa9928a43f3))
* **bin:** add output option ([1f11f7b](https://github.com/mbhall88/lrge/commit/1f11f7b0ebc964aaab5ed5c8f41259b346ff1b53))
* change default '--max-overhang-ratio' to 0.2 ([eb941d4](https://github.com/mbhall88/lrge/commit/eb941d4a4f054eedc1681b126ee465948ace3c83))
* completed ava strategy ([2a2fc35](https://github.com/mbhall88/lrge/commit/2a2fc35da9acd6120f30d2179c06b3c21fcaad57))
* completed twoset strategy ([190a4fa](https://github.com/mbhall88/lrge/commit/190a4fa60aad586b76213fb6750b9d08ce274f30))
* fit the reported interval quantiles per platform ([27cf272](https://github.com/mbhall88/lrge/commit/27cf272b9b45b7d7221e03353e02fcc54e757478))
* **liblrge:** cap the memory depth normalization spends on reads ([#49](https://github.com/mbhall88/lrge/issues/49)) ([796ae02](https://github.com/mbhall88/lrge/commit/796ae025984e885aa74e320a4330fd4a0700dc8f))
* **liblrge:** filter internal matches above a share of a run's overlaps ([437ec27](https://github.com/mbhall88/lrge/commit/437ec27f86ef20469b03f308f57c6cef78e7780d))
* **liblrge:** floor the number of reads skew detection samples ([#55](https://github.com/mbhall88/lrge/issues/55)) ([77efe95](https://github.com/mbhall88/lrge/commit/77efe952f0465757c103fedeeda1f4350315e6dd))
* **liblrge:** normalize depth-skewed read selection ([#45](https://github.com/mbhall88/lrge/issues/45)) ([c440d19](https://github.com/mbhall88/lrge/commit/c440d196ed01a61098aeb8142bbfcd6c3dac93bd))
* **lib:** output overlap PAF ([1e905fb](https://github.com/mbhall88/lrge/commit/1e905fb70e0e50c0862e5c637489e3118fb39317))
* **lib:** return number of no mapping reads ([917b786](https://github.com/mbhall88/lrge/commit/917b7869081d627edb0041fe343a5171af366640))
* remove option --max-overhang-size ([8753864](https://github.com/mbhall88/lrge/commit/875386408f7d718816a08bda6eb6398b22230bf0))
* update installation, minimap2, and add BAM/CRAM/SAM support ([#24](https://github.com/mbhall88/lrge/issues/24)) ([3ae5767](https://github.com/mbhall88/lrge/commit/3ae57671c47c639ae50ab907313045deb02039fd))


### Bug Fixes

* correct inverted internal-match filtering ([fc09157](https://github.com/mbhall88/lrge/commit/fc09157743d2df7b56af9942850f21f1a3d2d90d)), closes [#31](https://github.com/mbhall88/lrge/issues/31)
* **lib:** error on duplicate read IDs in ava ([86481e0](https://github.com/mbhall88/lrge/commit/86481e0f871253cd36c8667dffb83dbf61209b17))
* **liblrge:** keep the requested target:query ratio when the input is short ([#62](https://github.com/mbhall88/lrge/issues/62)) ([dd1bfd3](https://github.com/mbhall88/lrge/commit/dd1bfd37c1074485647da8f216a5a4ed6699d15e))
* make the documented install command work ([98ca62f](https://github.com/mbhall88/lrge/commit/98ca62f52dcb2bbc8b49024d7f27f59b0bafc3d6))


### Performance Improvements

* **liblrge:** score reads for retention in parallel ([#57](https://github.com/mbhall88/lrge/issues/57)) ([b96e482](https://github.com/mbhall88/lrge/commit/b96e4826d2b9e26358cbfa3d36f49fc8a1bbbcf1))
* **liblrge:** split depth detection from depth profiling ([#50](https://github.com/mbhall88/lrge/issues/50)) ([cb1b17d](https://github.com/mbhall88/lrge/commit/cb1b17d2c3205837c3b04c8131ee505dd091557c))

## [1.0.0](https://github.com/mbhall88/lrge/compare/lrge-0.3.0...lrge-1.0.0) (2026-09-08)


### ⚠ BREAKING CHANGES

* PacBio estimates change, because -P pb now selects the ava-pb overlap preset instead of being ignored. The reported interval changes on both platforms, because its quantiles are refitted and now depend on -P. liblrge::estimate::LOWER_QUANTILE and UPPER_QUANTILE keep their values but are deprecated in favour of Platform::interval_quantiles. The estimate line prints 95% CI or a named quantile pair where it used to print IQR. Anyone reproducing the published results should pin the 0.3.x series.
* **liblrge:** -F/--filter-contained no longer excludes every internal match. It now excludes them only when they exceed the share given to it, which defaults to 0.8, and that share must be written after an equals sign, as -F=0.5. Pass -F=0 for the previous behaviour. In liblrge, Builder::remove_internal is deprecated in favour of Builder::internal_filter, which takes an optional threshold rather than a flag; remove_internal still compiles, and asking it to filter is exactly a threshold of zero. Separately, an input too small to fill both read sets is now divided between them in the requested ratio rather than taking the whole shortfall from the target set, per #62, so every input holding fewer reads than -T and -Q ask for gets a different estimate from 0.3.0 on default flags; use --shortfall target for the previous behaviour.
* `lrge --max-overhang-ratio <FLOAT>` without `-F/--filter-contained` now exits with an error instead of silently ignoring the value, and `-F` filters the opposite set of alignments to before.

### Features

* fit the reported interval quantiles per platform ([27cf272](https://github.com/mbhall88/lrge/commit/27cf272b9b45b7d7221e03353e02fcc54e757478))
* **liblrge:** cap the memory depth normalization spends on reads ([#49](https://github.com/mbhall88/lrge/issues/49)) ([796ae02](https://github.com/mbhall88/lrge/commit/796ae025984e885aa74e320a4330fd4a0700dc8f))
* **liblrge:** filter internal matches above a share of a run's overlaps ([437ec27](https://github.com/mbhall88/lrge/commit/437ec27f86ef20469b03f308f57c6cef78e7780d))
* **liblrge:** floor the number of reads skew detection samples ([#55](https://github.com/mbhall88/lrge/issues/55)) ([77efe95](https://github.com/mbhall88/lrge/commit/77efe952f0465757c103fedeeda1f4350315e6dd))
* **liblrge:** normalize depth-skewed read selection ([#45](https://github.com/mbhall88/lrge/issues/45)) ([c440d19](https://github.com/mbhall88/lrge/commit/c440d196ed01a61098aeb8142bbfcd6c3dac93bd))


### Bug Fixes

* correct inverted internal-match filtering ([fc09157](https://github.com/mbhall88/lrge/commit/fc09157743d2df7b56af9942850f21f1a3d2d90d)), closes [#31](https://github.com/mbhall88/lrge/issues/31)
* **liblrge:** keep the requested target:query ratio when the input is short ([#62](https://github.com/mbhall88/lrge/issues/62)) ([dd1bfd3](https://github.com/mbhall88/lrge/commit/dd1bfd37c1074485647da8f216a5a4ed6699d15e))
* make the documented install command work ([98ca62f](https://github.com/mbhall88/lrge/commit/98ca62f52dcb2bbc8b49024d7f27f59b0bafc3d6))


### Performance Improvements

* **liblrge:** score reads for retention in parallel ([#57](https://github.com/mbhall88/lrge/issues/57)) ([b96e482](https://github.com/mbhall88/lrge/commit/b96e4826d2b9e26358cbfa3d36f49fc8a1bbbcf1))
* **liblrge:** split depth detection from depth profiling ([#50](https://github.com/mbhall88/lrge/issues/50)) ([cb1b17d](https://github.com/mbhall88/lrge/commit/cb1b17d2c3205837c3b04c8131ee505dd091557c))


### Dependencies

* The following workspace dependencies were updated
  * dependencies
    * liblrge bumped from 0.3.0 to 1.0.0

## [0.3.0](https://github.com/mbhall88/lrge/compare/lrge-0.2.1...lrge-0.3.0) (2026-05-01)


### Features

* update installation, minimap2, and add BAM/CRAM/SAM support ([#24](https://github.com/mbhall88/lrge/issues/24)) ([3ae5767](https://github.com/mbhall88/lrge/commit/3ae57671c47c639ae50ab907313045deb02039fd))


### Dependencies

* The following workspace dependencies were updated
  * dependencies
    * liblrge bumped from 0.2.1 to 0.3.0
