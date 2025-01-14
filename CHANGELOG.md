# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2021-07-04

### Added

- `-v,--verbose` flag to control logging
- Ability to specify the MAFFT binary [[#1][1]]

### Changed

- Logging switched to `loguru` [[#7][7]]

### Fixed

- Use the parent path of the update data structure to locate pickle objects

## [0.2.0] - 2021-05-07

### Changed

- Version bump from `0.2.0_prototype` to `0.2.0`

## [0.2.0_prototype] - 2021-03-22

### Added

- `update` command to `make_prg`, which is able to update PRGs based on denovo sequences;
- A script to build portable binary executables;
- A sample example;

### Changed

- `make_prg from_msa` now receives as input a directory with fasta files to build PRGs from;
- `make_prg from_msa` is now multithreaded;
- `make_prg from msa` output format changed: it now creates a single `.prg.fa` file and a `.update_DS` file,
  which contains data structures that make the PRG updateable;

### Removed
- Docker recipe;
- Nextflow pipeline;

## [0.1.1] - 2021-01-22
### Added
- Dockerfile
- `-V` option to get version

### Changed
- A test that was clustering all unique 5-mers was reduced to all 4-mers as the memory
  usage of all 5-mers was causing a segfault when trying to run the tests during the
  docker image build.

### Removed
- Singularity file as it is redundant with the new Dockerfile (that will be hosted on
  quay.io)
- `scipy` dependency. We never actually explicitly use `scipy`.

## [0.1.0] - 2021-01-20
### Added
- This CHANGELOG file to hopefully serve as an evolving example of a standardized open
  source project CHANGELOG.


[Unreleased]: https://github.com/rmcolq/make_prg/compare/v0.1.1...HEAD
[0.3.0]: https://github.com/leoisl/make_prg/releases/tag/v0.3.0
[0.2.0]: https://github.com/leoisl/make_prg/releases/tag/v0.2.0
[0.1.1]: https://github.com/rmcolq/make_prg/releases/v0.1.1
[0.1.0]: https://github.com/rmcolq/make_prg/releases/v0.1.0

[1]: https://github.com/leoisl/make_prg/issues/1
[7]: https://github.com/leoisl/make_prg/issues/7
