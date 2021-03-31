# Toy example

Here we present a walkthrough of running `make_prg` on a toy example.
We run:
1) `make_prg from_msa` to create PRGs from MSAs;
2) `make_prg update` to update PRGs with denovo paths;

## Input data description

```
msas/ : contains the MSAs of 4 genes we are using as toy example here;
denovo_paths/denovo_paths.txt : contains some denovo paths on 3 of these 4 genes. This file is produced by pandora;
```

## Dependencies

* **There is no need to have `make_prg` installed. The running script will automatically download
  and run the precompiled binary**;
* `MAFFT` has to be in your `PATH`. It can be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using conda: `conda install -c bioconda mafft`;
* `wget`;

## Running

```
./run_make_prg_on_sample_example.sh
```

### Quick look at the output

1. `msas_output/sample.prg.fa`: output from `make_prg from_msa`. It contains the PRGs built from the MSAs;

2. `msas_output/sample.update_DS`: output from `make_prg from_msa`. It contains data structures that make the PRGs updateable;

3. `msas_updated/updated_sample.prg.fa`: output from `make_prg update`. It contains the PRGs updated with denovo paths.

Diffing `msas_output/sample.prg.fa` and `msas_updated/updated_sample.prg.fa`, we can see that the updated PRGs have more sites/alleles:

```
$ diff msas_output/sample.prg.fa msas_updated/updated_sample.prg.fa 
4c4
<  5 AATAGGCCG 7  9 GATGCAGTTCAA 10 GATGCGGCGTA 9 AACGCCTTATCCGGCATACGA 11 ATTTATT 12 TTTTATT 11  8  13 G 14 A 13 ATGCGGCGTACGAATTTAT 15 T 16 C 15  7 CGGCCTGGCTCCCCGTAGGCCG 17 A 18 G 17 ATAAGATGCGCCAGCATCGCATCCGGCTATAATGC 19 G 20 A 19  6  21 TTCATTGG 22 TTCAATG 23 G 24 A 23  21 TTTATAATGCCTGATAAACGCACGGTCGATCCCCTCGCCCCTTCGGGGAGAGGATTAGGGTGAGGGGGTACAAGCCAGCCAGAGACCAGGCAA 25 TGACATG 26 CGACATG 25  5 CACATAACC 27 TCT 28 ACC 27 TGAAACT 29  31 CTTT 32 CGTC 32 CTTC 32 CATC 31 CCCAGAGCCTCTT 33 CAGC 34 TAGC 34 CAGT 33 CATCTATT 35 CA 36 TG 35  30  37 ACATCTCTTCA 38 ACATTTCTTCA 37  29 GGAGCAAACAATTTCAT 39 GCCAACTC 40 TCCAACTC 40 TCCAACTT 40 ACCAACTC 39 ATAACCCCAGCATATAAATCCAG 41 T 42 A 41 TGGTAACTTTT 43 A 44 C 43 TTTAACCT 45 G 46 A 45 AAACCAGTTT 47 TATCCAC 49 T 50 C 49  48 AATCCACC 47 ATTTATAAAATTATGTGAAGCATTTCATAGAAGAAAAATCACTGGC 51 C 52 T 51 TAAACATTAT 53 C 54 T 53 CCCTTTTTGC 55 CTGG 56 CTGA 56 ATGA 56 CTAG 56 CTTA 56 CTGT 55 TTTTTGACCATTTCCG 57 C 58 T 57 GATTTGTTACACATTGAAATATCACTTTTGCTGTGCGTAATATGGCTATTCGTTAGC 59 C 60 A 59 AAAAAATAAGAAAAGAT 61 T 62 A 61 
---
>  5 AATAGGCCG 7  9 GATGCAGTTCAA 10 GATGCGGCGTA 9 AACGCCTTATCCGGCATACGA 11 ATTTATT 12 TTTTATT 11  8  13 G 14 A 13 ATGCGGCGTACGAATTTAT 15 T 16 C 15  7 CGGCCTGGCTCCCCGTAGGCCG 17 A 18 G 17 ATAAGATGCGCCAGCATCGCATCCGGCTATAATGC 19 G 20 A 19  6  21 TTCATTGG 22 TTCAATG 23 G 24 A 23  21 TTTATAATGCCTGATAAACGCACGGTCGATCCCCTCGCCCCTTCGGGGAGAGGATTAGGGTGAGGGGGTACAAGCCAGCCAGAGACCAGGCAA 25 TGACATG 26 CGACATG 25  5 CACATAACC 27 TCT 28 ACC 27 TGAAACT 29  31 CTTT 32 CGTC 32 CTTC 32 CATC 31 CCCAGAGCCTCTT 33 CAGC 34 TAGC 34 CAGT 33 CATCTATT 35 CA 36 TG 35  30  37 ACATCTCTTCA 38 ACATTTCTTCA 37  29 GGAGCAAACAATTTCAT 39 GCCAACTC 40 TCCAACTC 40 TCCAACTT 40 ACCAACTC 39 ATAACCCCAGCATATAAATCCAG 41 T 42 A 41 TGGTAACTTTT 43 A 44 C 43 TTTAACCT 45 G 46 A 45 AAACCAGTTT 47 TATCCAC 49 T 50 C 49  48 AATCCACC 47 ATTTATAAAATTATGTGAAGCATTTCATAGAAGAAAAATCACTGGC 51 C 52 T 51 TAAACATTAT 53 C 54 T 53 CCCTTTTTGC 55 CTGG 56 CTGA 56 ATGA 56 CTAG 56 CTTA 56 CTGT 55 TTTTTGACCATTTCC 57 G 58 A 57  59 C 60  61 T 62 C 61  59 GATTTGTTACACATTGAAATATCACTTTTGCTGTGCGTAATATGGCTATTCGTTAGC 63 C 64 A 63 AAAAAATAAGAAAAGAT 65 T 66 A 65 
6c6
< TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTAATTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG 5 T 6 C 5 CTCCCTGGCTAAT 7 C 8 A 7 ACCACATTGGCATTTATGGAGCACATCACAATATTTCAATACCATTAAAGCACTGCA 9 C 10 T 9 CAAAATGAAACACTGCGA 11 C 12 T 11 ATTAAAATT 13 A 14 C 13 TTTCAATT
---
> TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTA 5 A 6 G 5 TTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG 7 T 8 C 7 CTCCCTGGCTAAT 9 C 10 A 9 ACCACATTGGCATTTATGGAGCACATCACAATATTTCAATACCATTAAAGCACTGCA 11 C 12 T 11 CAAAATGAAACACTGCGA 13 C 14 T 13 ATTAAAATT 15 A 16 C 15 TTTCAATT
8c8
< ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC 5 C 6 G 5 GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC 7 A 8 G 7 CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA 9 GCTA 10 ACTG 9 CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA 11 A 12 G 11 GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGGAGGGGGCACCAAAGGGGAAATGA
---
> ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCA 5 C 6 T 5 CGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC 7 C 8 G 7 GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC 9 A 10 G 9 CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA 11 GCTA 12 ACTG 11 CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA 13 A 14 G 13 GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGG 15 A 16 T 15 GGGGGCACCAAAGGGGAAATGA

```