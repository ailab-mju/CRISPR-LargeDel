# CRISPR-LargeDel 
Short description

## Installation
In order to run a large deletion program, you must proceed with complie.

```bash
git clone https://github.com/ailab-mju/CRISPR-LargeDel.git
cd CRISPR-LargeDel
make
```

## Run
To run a program, you must run it with parameters in the shell script.
bash run_kmer_algin.sh {Genename} {Cellname} {Target start position} {Target end position} {Clevage site from start position} {Control sample} {Cas9 Mutation sample} {K}
```bash
bash run_kmer_align.sh HPRT1 TEST 4652 4671 3 test_con test_mut 11
