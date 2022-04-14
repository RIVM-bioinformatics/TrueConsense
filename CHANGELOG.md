# Changelog

## [0.4.0](https://www.github.com/RIVM-bioinformatics/TrueConsense/compare/v0.3.0...v0.4.0) (2022-04-14)


### ⚠ BREAKING CHANGES

* Use only a single minimum coverage level

### Features

* Use only a single minimum coverage level ([cdd29b7](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/cdd29b7a41f78680c56b932ccbc9efe40a705f60))


### Bug Fixes

* Fix WriteOutputs number of args ([8bf72f8](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/8bf72f875c850520a16ad1484c7d4019d86f83a9))


### Performance Improvements

* Use pysam directly in stead of pysamstats ([afd8c0a](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/afd8c0a9769d344ab32d90cb81b9c9368742e92d))


### Dependencies

* remove parmap as a dependency ([0aaed4d](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/0aaed4d2d803227ee1e000cea2369bf91fae0b93))


### Documentation

* add docstring to `TrueConsense.Outputs` ([4855ee9](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/4855ee94e2cc3876bdc3d89af3b51243580ecbe6))
* add docstrings to functions in `TrueConsense.Ambig` ([4e8a138](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/4e8a138f24ca253a4b73f5c9a6671c77883f34e8))
* add docstrings to functions in `TrueConsense.Coverage` ([47d8c89](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/47d8c89a2a220842559a426b193d49d562530051))
* add docstrings to functions in `TrueConsense.Events` ([6d754e3](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/6d754e3e82419e93192ae2a79440f49310062cb8))
* add docstrings to functions in `TrueConsense.indexing` ([650a5d1](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/650a5d1e41f9acca588fa7c8eed71b41b77bebd0))
* add docstrings to functions in `TrueConsense.Sequences` ([5abdb81](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/5abdb81cb5bf44846fb912b2fd4f8c2344335b89))
* add docstrings to functions of `TrueConsense.ORFs` ([b52a750](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/b52a750ed9094004e1bc90af89922d98c911b847))
* update readme and docs-index to support a single coverage level ([f9fd73e](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/f9fd73e99090ec3a70842b3df321778b29e1880c))
* update TrueConsense docs to only feature a single coverage level ([53cecc3](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/53cecc3e4fdabd50a7daa642f4afa2d91b974e86))

## [0.3.0](https://www.github.com/RIVM-bioinformatics/TrueConsense/compare/v0.2.0...v0.3.0) (2021-11-25)


### Features

* add arg to override index with pre-made index on certain positions ([9b3227c](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/9b3227c74a420dd3861cacba4058b6ab34e315a3))


### Documentation

* add changelog page ([2f29130](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/2f2913028bdd13bf067a7b27be1746241867ac3a))
* add explanation of `--index-override` flag to docs ([ac4ac4b](https://www.github.com/RIVM-bioinformatics/TrueConsense/commit/ac4ac4b9c5ce49336602fa364ccc04203741d522))
