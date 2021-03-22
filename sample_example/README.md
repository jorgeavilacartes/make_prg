# Toy example

Here we present a walkthrough of running `make_prg` on a toy example.
We run:
1) `make_prg from_msa` to create PRGs from MSAs;
2) `make_prg update` to update PRGs with denovo paths;

## Input data description

```
msas/ : contains the MSAs of 4 genes we are using as toy example here;
denovo_paths/denovo_paths.txt : contains some denovo paths on 2 of these 7 genes;
```

## Dependencies

* `MAFFT` has to be in your `PATH`. It can be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using conda: `conda install -c bioconda mafft`;
* `md5sum`;
* `wget`;

## Running

```
./run_make_prg_on_sample_example.sh
```

### Quick look at the output

1. `msas_output/sample.prg.fa`: output from `make_prg from_msa`. It contains the PRGs built from the MSAs;

2. `msas_output/sample.update_DS`: output from `make_prg from_msa`. It contains data structures that make the PRGs updateable;

3. `msas_updated/updated_sample.prg.fa`: output from `make_prg update`. It contains the PRGs updated with denovo paths.

