# CRISPR-LargeDel 
Short description

## Installation
In order to run a large deletion program, you must proceed with complie.

```bash
git clone https://github.com/ailab-mju/CRISPR-LargeDel.git
cd CRISPR-LargeDel
make
```

## Usage
You need to provide reference fasta file, control and cas9 mutation paired end fastq files. <br>
The following sample files are included in the test_data directory: <br>
test_data/HPRT1.fa <br>
test_data/test_con_1.fastq <br>
test_data/test_con_2.fastq <br>
test_data/test_mut_1.fastq <br>
test_data/test_mut_2.fastq <br>
To run a program, you must run it with parameters in the shell script. <br>

bash run_kmer_algin.sh {Genename} {Cellname} {Target start position} {Target end position} {Clevage site from start position} {Control sample} {Cas9 Mutation sample} {K} <br>
```bash
bash run_kmer_align.sh HPRT1 TEST 4652 4671 3 test_con test_mut 11
