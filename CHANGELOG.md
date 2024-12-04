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