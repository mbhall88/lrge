## 0.2.1 ([12443db..14061b7](https://github.com/mbhall88/lrge/compare/lrge-0.2.0...lrge-0.2.1))
#### Bug Fixes
- **(liblrge)** change warn logging to info for no overlaps - ([2a11150](https://github.com/mbhall88/lrge/commit/2a11150fe39bec395d9780570d9df84641e03b00)) - [@mbhall88](https://github.com/mbhall88)
#### Build system
- **(liblrge)** update crossbeam - ([bb24d78](https://github.com/mbhall88/lrge/commit/bb24d7829a63aa8cd7a9adaaadd5c36d50f5e9ce)) - [@mbhall88](https://github.com/mbhall88)
#### CI
- **(deps)** bump actions/attest-build-provenance from 2 to 3 (#11) - ([53023f8](https://github.com/mbhall88/lrge/commit/53023f8f16f96eb842596131dd321f372a2b9322)) - dependabot[bot]
- **(deps)** bump actions/checkout from 4 to 5 (#9) - ([1faca39](https://github.com/mbhall88/lrge/commit/1faca39cf3494e2108e2f69897142aa20f35b784)) - dependabot[bot]
#### Documentation
- add guidance for selecting readset sizes - ([bc30ba0](https://github.com/mbhall88/lrge/commit/bc30ba0c78bee5f559aaf7f3b5b843a128370e43)) - [@mbhall88](https://github.com/mbhall88)

## 0.2.0 ([[30c8e5d..12443d](https://github.com/mbhall88/lrge/compare/lrge-0.1.3...lrge-0.2.0))
#### Bug Fixes
- aviod minimap2 multi-part indexing - ([e5e9a77](https://github.com/mbhall88/lrge/commit/e5e9a773336638a8f57384d76913e34cae483a66)) - Chenxi Zhou
#### Build system
- **(deps)** update proc-macro2 for docs - ([cf5b55b](https://github.com/mbhall88/lrge/commit/cf5b55bd3bc9cdc4ed3b80007ed45076dd034002)) - [@mbhall88](https://github.com/mbhall88)
- **(deps)** update deps - ([8b78ed1](https://github.com/mbhall88/lrge/commit/8b78ed110e967f09ed4481ec51dbdb25e4bf860f)) - [@mbhall88](https://github.com/mbhall88)
#### Features
- change default '--max-overhang-ratio' to 0.2 - ([eb941d4](https://github.com/mbhall88/lrge/commit/eb941d4a4f054eedc1681b126ee465948ace3c83)) - Chenxi Zhou
- add option '--use-min-ref' for using smaller reference - ([4791ec0](https://github.com/mbhall88/lrge/commit/4791ec066dd204c97cd91884b473c3bcf50e2931)) - Chenxi Zhou
- a more accurate calculation of per read estimation - ([5d2e211](https://github.com/mbhall88/lrge/commit/5d2e211e976ab7df9d3db05243175c9fda210b90)) - Chenxi Zhou
- remove option --max-overhang-size - ([8753864](https://github.com/mbhall88/lrge/commit/875386408f7d718816a08bda6eb6398b22230bf0)) - Chenxi Zhou
- use smaller dataset Q/T as reference - ([f4a4b73](https://github.com/mbhall88/lrge/commit/f4a4b73283062435c8ca8958e4fe71fd796371d7)) - Chenxi Zhou
- add option '-F' to remove internal matches - ([b23a8c7](https://github.com/mbhall88/lrge/commit/b23a8c7fea82aff1344c4f108628cf6044245245)) - Chenxi Zhou
#### Refactoring
- testing for internal alignments is done on PafRecord - ([e5dbc3f](https://github.com/mbhall88/lrge/commit/e5dbc3f2454cbca1753a76aed0a884c703f2dea5)) - [@mbhall88](https://github.com/mbhall88)

## 0.1.3 ([2467dab..01de783](https://github.com/mbhall88/lrge/compare/2467dab..01de783))
#### Bug Fixes
- handle multi-member gzip compressed files - ([2412f4a](https://github.com/mbhall88/lrge/commit/2412f4ab52c9fd8baf3e310fd637f389228ac7f6)) - [@mbhall88](https://github.com/mbhall88) - fixes <https://github.com/gbouras13/hybracter/issues/110>

## 0.1.2 ([495d160..dbe9bf3](https://github.com/mbhall88/lrge/compare/495d160..dbe9bf3))
#### Bug Fixes
- **(liblrge)** FASTQ read ID was not splitting on tabs - ([9fd66b7](https://github.com/mbhall88/lrge/commit/9fd66b77fa971fc17086686faf3c1da8fd8111d4)) - [@mbhall88](https://github.com/mbhall88) with thanks to bug hunter 🐛 [@gbouras13](https://github.com/gbouras13) 🐛 
#### Build system
- **(deps)** update minimap2-sys to 0.1.20 - ([d3ed478](https://github.com/mbhall88/lrge/commit/d3ed478b86ba813a22017b5d8f4546f5cb35338a)) - [@mbhall88](https://github.com/mbhall88)
#### Documentation
- add link to preprint - ([17285e8](https://github.com/mbhall88/lrge/commit/17285e86a1f72bcb52e3f260ff1a0f6af438b529)) - [@mbhall88](https://github.com/mbhall88)
#### Refactoring
- **(liblrge)** use CString for read ID to ensure null termination - ([28d343b](https://github.com/mbhall88/lrge/commit/28d343b4fb0b4974dad8f9e1af6d4061e7d2077c)) - [@mbhall88](https://github.com/mbhall88)

## 0.1.1 ([33bfd4b..6236d64](https://github.com/mbhall88/lrge/compare/33bfd4b..6236d64))
#### Bug Fixes
- **(liblrge)** dont assume 8-bit intege size - use c_char - ([5e630e7](https://github.com/mbhall88/lrge/commit/5e630e76def0d01592896334816a972019c71a9f)) - [@mbhall88](https://github.com/mbhall88)
#### Build system
- **(deps)** bump docker/build-push-action from 5 to 6 (#3) - ([0058b4b](https://github.com/mbhall88/lrge/commit/0058b4b818de05d76c4826876b4b4f58a96052b7)) - dependabot[bot]