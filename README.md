# CRISPR-LargeDel 
Short description

## Installation
To compile and run the large deletion program, follow these steps:

```bash
git clone https://github.com/ailab-mju/CRISPR-LargeDel.git
cd CRISPR-LargeDel
make
```

## Usage

To utilize the program, you'll need a reference fasta file and paired-end fastq files for control and Cas9 mutation samples. Sample files are available in the `test_data` directory:

- `test_data/HPRT1.fa`
- `test_data/test_con_1.fastq`
- `test_data/test_con_2.fastq`
- `test_data/test_mut_1.fastq`
- `test_data/test_mut_2.fastq`

Run the program using the provided shell script with the following parameters:
```bash
bash run_kmer_align.sh {Genename} {Cellname} {Target start position} {Target end position} {Cleavage site from start position} {Control sample} {Cas9 Mutation sample} {K}
```
For example:
```bash
bash run_kmer_align.sh HPRT1 TEST 4652 4671 3 test_con test_mut 11
```


