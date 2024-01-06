# CRISPR-LargeDel 
Short description

## Installation
To compile and run the CRISPR-LargeDel, follow these steps:

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

After running the CRISPR-LargeDel program, the following result files will be generated:
```bash
cd result
cd HPRT1_TEST_11
```

- `HPRT1_TEST_ctr.sam`: SAM file containing alignment results for the control sample.
- `HPRT1_TEST_ctr.sam`: SAM file containing alignment results for the Cas9 mutatoin sample.
- `HPRT1_TEST_ctr_only.sam`: SAM file containing large deletion unique alignments from the control sample.
- `HPRT1_TEST_ctr_only.sam`: SAM file containing large deletion unique alignments from the Cas9 mutation sample. 
- `HPRT1_TEST_ctr_common.sam`: SAM file containing large deletion common alignments between control and Cas9 mutation samples.
- `HPRT1_TEST_mut_common.sam`: SAM file containing large deletion common alignments between control and Cas9 mutation samples.
- `HPRT1_TEST_site_ctr.sam`: SAM file containing alignments around the target site for the control sample.
- `HPRT1_TEST_site_mut.sam`: SAM file containing alignments around the target site for the Cas9 mutation sample.
- `HPRT1_TEST_site_ctr_sm.sam`: SAM file with additional small indels around the target site for the control sample.
- `HPRT1_TEST_site_mut_sm.sam`: SAM file with additional small indels around the target site for the Cas9 mutation sample.

These files provide detailed information about the sequence alignments and can be useful for further analysis or visualization.

## License

