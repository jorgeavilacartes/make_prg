# make_prg

A tool to create and update PRGs for input to [Pandora][pandora] and [Gramtools][gramtools] from a set of 
Multiple Sequence Alignments.

[TOC]: #

## Table of Contents
- [Install](#install)
  - [No installation needed - precompiled portable binary](#conda)
  - [pip](#pip)
- [Running on a sample example](#usage)
- [Usage](#usage)

## Install

### No installation needed - precompiled portable binary

You can use `make_prg` with no installation at all by simply downloading the precompiled binary, and running it.
In this binary, all libraries are linked statically.

* **Requirements**: None

* **Download**:
  ```
  wget https://github.com/leoisl/make_prg/releases/download/v0.2.0_prototype/make_prg_0.2.0_prototype
  ```
* **Running**:
```
chmod +x make_prg_0.2.0_prototype
./make_prg_0.2.0_prototype -h
```

* **Credits**:
  * Compilation is done using [PyInstaller](https://github.com/pyinstaller/pyinstaller).

* **Notes**:
  * We provide precompiled binaries for Linux OS only;


### pip

* **Requirements**: `python>=3`

* **Installing**:
```sh
pip install git+https://github.com/leoisl/make_prg
```

* **Running**:
```
./make_prg -h
```

## Running on a sample example

See [sample example](sample_example).

## Usage

```
$ make_prg --help
usage: make_prg <subcommand> <options>

Subcommand entrypoint

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  show program's version number and exit

Available subcommands:
  
    from_msa     Make PRG from multiple sequence alignment dir
    update       Update PRGs given new sequences output by pandora.

```

#### `from_msa`

```
$ make_prg from_msa --help
usage: make_prg from_msa

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input dir: all files in this will try to be read as the supported alignment_format. If not aligned in fasta alignment_format, use -f to input the alignment_format type
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output prefix: prefix for the output files
  -t THREADS, --threads THREADS
                        Number of threads
  -f ALIGNMENT_FORMAT, --alignment_format ALIGNMENT_FORMAT
                        Alignment format of MSA, must be a biopython AlignIO input alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta
  --max_nesting MAX_NESTING
                        Maximum number of levels to use for nesting. Default: 5
  --min_match_length MIN_MATCH_LENGTH
                        Minimum number of consecutive characters which must be identical for a match. Default: 7
```

#### `update`

```
$ make_prg update --help
usage: make_prg update_prg

optional arguments:
  -h, --help            show this help message and exit
  -u UPDATE_DS, --update_DS UPDATE_DS
                        Filepath to the update data structures. Should point to a file *.update_DS.
  -d DENOVO_PATHS, --denovo_paths DENOVO_PATHS
                        Filepath containing denovo sequences output by pandora. Should point to a denovo_paths.txt file.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output prefix: prefix for the output files
  -t THREADS, --threads THREADS
                        Number of threads
  --keep_temp           Keep temp files.
```

[gramtools]: https://github.com/iqbal-lab-org/gramtools
[pandora]: https://github.com/rmcolq/pandora
