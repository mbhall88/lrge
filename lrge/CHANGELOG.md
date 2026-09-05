# Changelog

## [1.0.0](https://github.com/mbhall88/lrge/compare/lrge-0.3.0...lrge-1.0.0) (2026-09-05)


### ⚠ BREAKING CHANGES

* `lrge --max-overhang-ratio <FLOAT>` without `-F/--filter-contained` now exits with an error instead of silently ignoring the value, and `-F` filters the opposite set of alignments to before.

### Features

* **liblrge:** cap the memory depth normalization spends on reads ([#49](https://github.com/mbhall88/lrge/issues/49)) ([796ae02](https://github.com/mbhall88/lrge/commit/796ae025984e885aa74e320a4330fd4a0700dc8f))
* **liblrge:** floor the number of reads skew detection samples ([#55](https://github.com/mbhall88/lrge/issues/55)) ([77efe95](https://github.com/mbhall88/lrge/commit/77efe952f0465757c103fedeeda1f4350315e6dd))
* **liblrge:** normalize depth-skewed read selection ([#45](https://github.com/mbhall88/lrge/issues/45)) ([c440d19](https://github.com/mbhall88/lrge/commit/c440d196ed01a61098aeb8142bbfcd6c3dac93bd))


### Bug Fixes

* correct inverted internal-match filtering ([fc09157](https://github.com/mbhall88/lrge/commit/fc09157743d2df7b56af9942850f21f1a3d2d90d)), closes [#31](https://github.com/mbhall88/lrge/issues/31)


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
