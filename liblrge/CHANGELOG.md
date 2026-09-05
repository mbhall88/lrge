# Changelog

## [1.0.0](https://github.com/mbhall88/lrge/compare/liblrge-0.3.0...liblrge-1.0.0) (2026-09-05)


### ⚠ BREAKING CHANGES

* `lrge --max-overhang-ratio <FLOAT>` without `-F/--filter-contained` now exits with an error instead of silently ignoring the value, and `-F` filters the opposite set of alignments to before.

### Features

* **liblrge:** cap the memory depth normalization spends on reads ([#49](https://github.com/mbhall88/lrge/issues/49)) ([796ae02](https://github.com/mbhall88/lrge/commit/796ae025984e885aa74e320a4330fd4a0700dc8f))
* **liblrge:** detect depth skew from minimizer counts ([#44](https://github.com/mbhall88/lrge/issues/44)) ([b49f345](https://github.com/mbhall88/lrge/commit/b49f345eeca2e2f96f8a4ae9958574a78ccec374))
* **liblrge:** floor the number of reads skew detection samples ([#55](https://github.com/mbhall88/lrge/issues/55)) ([77efe95](https://github.com/mbhall88/lrge/commit/77efe952f0465757c103fedeeda1f4350315e6dd))
* **liblrge:** normalize depth-skewed read selection ([#45](https://github.com/mbhall88/lrge/issues/45)) ([c440d19](https://github.com/mbhall88/lrge/commit/c440d196ed01a61098aeb8142bbfcd6c3dac93bd))
* **liblrge:** support weighted read selection ([28dad47](https://github.com/mbhall88/lrge/commit/28dad4748fc5eeb3eee5330a76dd585d4a556bad))


### Bug Fixes

* correct inverted internal-match filtering ([fc09157](https://github.com/mbhall88/lrge/commit/fc09157743d2df7b56af9942850f21f1a3d2d90d)), closes [#31](https://github.com/mbhall88/lrge/issues/31)
* **liblrge:** draw the depth sample from the reads detection sees ([#59](https://github.com/mbhall88/lrge/issues/59)) ([f1e4189](https://github.com/mbhall88/lrge/commit/f1e4189d843c74dc91be3394bf81e3a62010795a))


### Performance Improvements

* **liblrge:** cut the per-minimizer cost of the depth sketch ([#52](https://github.com/mbhall88/lrge/issues/52)) ([67f74f9](https://github.com/mbhall88/lrge/commit/67f74f93faaa6be9fd208bf654364654e6a077a5))
* **liblrge:** read minimizer values from a packed copy of the read ([#58](https://github.com/mbhall88/lrge/issues/58)) ([b607171](https://github.com/mbhall88/lrge/commit/b6071717ab733ff4036e66cb39dcd8c6d2bd1b77))
* **liblrge:** score reads for retention in parallel ([#57](https://github.com/mbhall88/lrge/issues/57)) ([b96e482](https://github.com/mbhall88/lrge/commit/b96e4826d2b9e26358cbfa3d36f49fc8a1bbbcf1))
* **liblrge:** split depth detection from depth profiling ([#50](https://github.com/mbhall88/lrge/issues/50)) ([cb1b17d](https://github.com/mbhall88/lrge/commit/cb1b17d2c3205837c3b04c8131ee505dd091557c))

## [0.3.0](https://github.com/mbhall88/lrge/compare/liblrge-0.2.1...liblrge-0.3.0) (2026-05-01)


### Features

* update installation, minimap2, and add BAM/CRAM/SAM support ([#24](https://github.com/mbhall88/lrge/issues/24)) ([3ae5767](https://github.com/mbhall88/lrge/commit/3ae57671c47c639ae50ab907313045deb02039fd))
