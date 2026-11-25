# CRISPR-LargeDel 
Short description

## Installation

Requirements:

- a C++ compiler with C++14 support (g++ 7 or newer, or clang++) and `make`
- optional: htslib, for BAM output (`sudo apt install libhts-dev`, or
  `conda install -c bioconda htslib`). Without it the build writes SAM and says so
  at run time; the calls are identical.

```bash
make                 # uses htslib automatically when it is found
make HTS=0           # force a SAM-only build
make CPP=clang++     # another compiler
```

Tested 2026-10-23 with g++ 11.4 (no htslib), g++ 13 (with and without htslib) and
clang++ 18: all four builds give identical output on the HPRT1 test data shipped then. The language level
is pinned to `-std=gnu++14` because `struct data` collides with `std::data` under
C++17. g++-7 is no longer required.

Files in this folder:

| file | what it is |
|---|---|
| `*.cpp`, `*.h`, `dkm.hpp` | the caller (`kmer_align`) |
| `Makefile` | build |
| `run_kmer_align.sh` | wrapper that reads `test_data/` and writes `result/` |
| `test_data/` | simulated TRAC example, the same as the web page's **Load example** (5,899 bp amplicon, 53,363 read pairs per sample, plain FASTQ) |
| `result/` | written by the example run |
| `ka_report.py` | optional: combine several result folders into one HTML report (Python 3, no extra packages; `./kmer_align report` does the same natively) |
| `docs/` | background notes on the fixes (`FIXES.md`, `INSERTION_FIX.md`, `LARGEDEL_CLUSTER_INSTABILITY.md`) |

## Usage

To utilize the program, you'll need a reference fasta file and paired-end fastq files for control and Cas9 mutation samples. Sample files are available in the `test_data` directory:

- `test_data/TRAC.fa`
- `test_data/TRAC_SIM_ctr_1.fastq`
- `test_data/TRAC_SIM_ctr_2.fastq`
- `test_data/TRAC_SIM_mut_1.fastq`
- `test_data/TRAC_SIM_mut_2.fastq`

The reads are simulated on a TRAC amplicon (SpCas9 target AGAGCAACAGTGCTGTGGCC on the plus
strand), so the numbers show how the output looks, not a real editing outcome.

Run the program using the provided shell script with the following parameters:
```bash
bash run_kmer_align.sh {Genename} {Cellname} {Target start position} {Target end position} {Cleavage site from start position} {Control sample} {Cas9 Mutation sample} {K}
```
For example:
```bash
bash run_kmer_align.sh TRAC SIM 2952 2971 17 TRAC_SIM_ctr TRAC_SIM_mut 11
```
The wrapper reads `<prefix>_1.fastq` when it exists and `<prefix>_1.fastq.gz` otherwise.
The target start and end are 0-based; the cleavage offset is 17 for an SpCas9 target on the
plus strand (3 on the minus strand). It takes a few seconds.

FASTQ inputs may be gzip, bzip2 or xz compressed - pass `reads.fastq.gz` directly,
with no decompression step (it is in fact slightly faster, since the decompressor
runs in its own process). The reference FASTA may be wrapped at any line width;
earlier versions silently used only its first sequence line, which produced a 0%
alignment rate with no error message. The reference length actually read is now
printed as `[ref] <file>: <N> bp` at the start of every run.

Every run also writes `<gene>_<sample>_report.html` next to its other outputs:
a single self-contained HTML report (no external scripts) with the indel size
spectrum, the pattern composition, the large-deletion size distribution, and a
table of insertions of 8 bp and over, each checked against the control's raw
reads so that co-amplified background is flagged rather than reported as an
edit. Argument 21 (`report`, default 1) turns it off. To combine several runs,
or to compare against the same samples run with `ins_fix=0`:

```bash
./kmer_align report -o PD1.html --compare result/*_OFF -- result/*_ON
```

Arguments 19 and 20 control insertion handling. `ins_fix` defaults to 1, which
records a query-side gap between two seed blocks as an `I` in the CIGAR instead of
sliding the next block over the inserted bases; this is what lets insertions of
9 bp and up align at all, and what gives them a size. Setting it to 0 restores the
pre-2026 behaviour. `ins_max` caps the insertion length that will be accepted
(0 = no cap). See `docs/INSERTION_FIX.md` for the measurements behind this.

### Optional arguments by name

Every optional argument can be given by name, in any order, as `--name value`
or `--name=value`; anything not named keeps its default, so there is no need
to fill in the positions in front of it:

```bash
./kmer_align test_data/TRAC.fa test_data/TRAC_SIM_ctr_1.fastq test_data/TRAC_SIM_ctr_2.fastq \
    test_data/TRAC_SIM_mut_1.fastq test_data/TRAC_SIM_mut_2.fastq \
    2952 2971 17 TRAC S1 11 out --out_fmt sam --clip_rescue 1
```

| option | position | default |
|---|---|---|
| `--dist_cut` | 13 | 8000 |
| `--full_sam` | 14 | 1 |
| `--out_fmt` | 15 | bam |
| `--min_gap` | 16 | 400 |
| `--anchor_min` | 17 | 13 |
| `--cut_between` | 18 | 1 |
| `--ins_fix` | 19 | 1 |
| `--ins_max` | 20 | 0 |
| `--report` | 21 | 1 |
| `--large_del_len` | 22 | 50 |
| `--clip_rescue` | 23 | 1 |
| `--common_mode` | 24 | local |

`--dist_cut` is the largest distance between the cut site and the nearest read
end of a large-deletion call (split read or read pair). The bigger a deletion,
the farther its reads sit from the cut, so this value is effectively a ceiling
on the deletion size. On the noPCR simulations (HPRT1, TRAC, BCL11A, 1-7 kb
deletions) 4000, 8000 and 20000 give identical calls, while 1000 keeps 0-62 %
of the 5-7 kb deletions instead of 78-91 %; small indels are unaffected. The
default is therefore 8000 (it was 1000 before 2026-10-23 - pass
`--dist_cut 1000` to reproduce earlier results).

| HPRT1 noPCR, deletions kept | 1000 | 2000 | 4000 | 8000 |
|---|---|---|---|---|
| 3 kb | 58.6 % | 85.8 % | 86.2 % | 86.2 % |
| 5 kb | 13.2 % | 56.4 % | 77.8 % | 77.8 % |
| 7 kb | 0.1 % | 0.1 % | 82.0 % | 82.0 % |

`--large_del_len` is the deletion length inside one read that makes it a
split-read large deletion (internally `MIN_skip`). The default is 50 bp, the
original tool's compile-time value; on the stress simulation 50-99 bp deletions
were found 69 % of the time at 50 and 40 % at 100. Pass `--large_del_len 100`
for the paper's large-deletion definition (the default during 2026-10).
The old positional form (arguments 13-24 in this order) still works and gives
identical results; an option given both ways must have the same value, and an
unknown `--option` is an error.

Both the script and the binary print this argument list themselves when they
are run with no arguments, with `-h`, or with `--help`:
```bash
bash run_kmer_align.sh          # wrapper usage
./kmer_align                    # full argument list, incl. every --option
```
Every argument is validated before the run starts: missing FASTA/FASTQ files,
non-numeric positions, `target_end <= target_start`, a `K` outside 5..31 and a
missing or unwritable output directory are each reported by argument number,
followed by the usage block. Exit status is 0 for a successful run or a plain
usage request, and 1 for a bad argument.

After running the CRISPR-LargeDel program, the following result files will be generated:
```bash
cd result
cd TRAC_SIM_11
```

- `TRAC_SIM_ctr.sam`: SAM file containing alignment results for the control sample.
- `TRAC_SIM_mut.sam`: SAM file containing alignment results for the Cas9 mutation sample.
- `TRAC_SIM_ctr_only.sam`: SAM file containing large deletion unique alignments from the control sample.
- `TRAC_SIM_mut_only.sam`: SAM file containing large deletion unique alignments from the Cas9 mutation sample.
- `TRAC_SIM_ctr_common.sam`: SAM file containing large deletion common alignments between control and Cas9 mutation samples.
- `TRAC_SIM_mut_common.sam`: SAM file containing large deletion common alignments between control and Cas9 mutation samples.
- `TRAC_SIM_site_ctr.sam`: SAM file containing alignments around the target site for the control sample.
- `TRAC_SIM_site_mut.sam`: SAM file containing alignments around the target site for the Cas9 mutation sample.
- `TRAC_SIM_site_ctr_sm.sam`: SAM file with additional small indels around the target site for the control sample.
- `TRAC_SIM_site_mut_sm.sam`: SAM file with additional small indels around the target site for the Cas9 mutation sample.
- `TRAC_SIM_site_ctr_clip.sam` / `TRAC_SIM_site_mut_clip.sam`: SAM files with reads carrying a terminal soft clip at the cut site, i.e. an insertion at the edge of the read.
- `TRAC_SIM_mutation_summary.tsv`: per mutation pattern read pair counts and percentages for both samples.
- `TRAC_SIM_report.html` / `TRAC_SIM_report.tsv`: the standalone HTML report and the run's parameters and counts.

On the example, the Cas9-treated sample has 2,423 read pairs used: 16.9% WT, 21.8% small
insertions, 23.6% small deletions, 31.4% large deletions, 5.2% terminal insertions and 1.0%
substitutions; the control is 94.8% WT.
All of these are written as `.bam` + `.bam.bai` by default; `OUT_FMT=sam` gives the `.sam` files named above instead.

These files provide detailed information about the sequence alignments and can be useful for further analysis or visualization.

## Mutation pattern summary

At the end of a run the same table is printed to stdout for each sample and
written to `{Genename}_{Cellname}_mutation_summary.tsv`:

```
###### Cas9 Mutation sample mutation pattern summary
  Pattern                          Read pairs    Percent
  ---------------------------------------------------------
  WT (no mutation)                         10     55.56%
  Small insertion                           1      5.56%
  Small deletion                            1      5.56%
  Terminal insertion (clipped)              1      5.56%
  Substitution                              0      0.00%
  Large deletion                            5     27.78%
  ---------------------------------------------------------
  Total analyzed                           18    100.00%
  Not covered (excluded)                    0          -
  Large del (common, excluded)              0          -
  Aligned pairs (all positions)          1411          -
```

Counts are **read pairs**, and each row corresponds one-to-one to a SAM file,
so the numbers can be reproduced with `wc -l file.sam` divided by two:

| Row | Source file |
| --- | --- |
| WT (no mutation) | `_site_mut.sam` / `_site_ctr.sam` |
| Small insertion / deletion / substitution | `_site_mut_sm.sam` / `_site_ctr_sm.sam`, split by CIGAR |
| Terminal insertion (clipped) | `_site_mut_clip.sam` / `_site_ctr_clip.sam` |
| Large deletion | `_mut_only.sam` / `_ctr_only.sam` |
| Large del (common, excluded) | `_mut_common.sam` / `_ctr_common.sam` |
| Aligned pairs (all positions) | `_mut.sam` / `_ctr.sam` |

`Total analyzed` is the read pairs that actually interrogate the cut site,
plus the sample specific large deletions. Every row shares this one
denominator, so they add up to 100%.

### What counts as interrogating the cut site

A pair is only counted when **at least one of its two mates covers the cut
window itself** - either its alignment spans cut +/-13 bp, or it ends in a
qualifying terminal clip right at the cut. Two mates sitting either side of the
cut with the variant hidden in the unsequenced gap between them say nothing
about that molecule, so they are reported on the `Not covered (excluded)` line
and stay out of every rate.

This matters more than it sounds. The test used to be applied to the pair's
outer coordinates, so those bracketing pairs passed, aligned perfectly, and were
counted as WT. On a homozygous edited sample that produced a spurious ~47% WT.
The fraction of bracketing pairs is set by library insert size rather than by
biology and differs a lot between samples, so leaving them in made the control
and treated numbers incomparable as well.

Large deletions are unaffected by this test - they are called from the pair
geometry before it is reached. Clusters shared with the other sample are
excluded from the denominator as PCR artefacts and are reported on their own
line.

### BAM output — one binary, no samtools

`kmer_align` links htslib, so it writes **coordinate-sorted BAM with a `.bai`
index by itself**. There is no samtools call anywhere in the pipeline and no
multi-GB intermediate SAM: records are bucketed by position into small temp
files as the run proceeds and merged back in order at close, so memory stays
bounded whatever the depth.

```bash
make                       # prints whether htslib was found
bash run_kmer_align.sh TRAC SIM 2952 2971 17 TRAC_SIM_ctr TRAC_SIM_mut 11
# -> TRAC_SIM_site_mut.bam + TRAC_SIM_site_mut.bam.bai, straight into IGV

OUT_FMT=sam  bash run_kmer_align.sh ...   # plain SAM instead
OUT_FMT=both bash run_kmer_align.sh ...   # both, identical content
```

The raw binary takes the format as `--out_fmt sam|bam|both` (or as its 15th
positional argument), and `./kmer_align` with no arguments prints the whole list.

## Large deletion filters

Two filters run by default on the large-deletion path. Both were measured
against the stress simulation (ground truth in the read names) and against PRR5,
a line carrying one +16 bp insertion and **no** real large deletions, so every
large deletion it reports is an artefact.

```bash
MIN_GAP=400 ANCHOR_MIN=13 CUT_BETWEEN=1  bash run_kmer_align.sh ...   # defaults
ANCHOR_MIN=0 CUT_BETWEEN=0               bash run_kmer_align.sh ...   # old behaviour
ANCHOR_MIN=16                            bash run_kmer_align.sh ...   # stricter
```

`CUT_BETWEEN=1` requires the cleavage site to fall strictly between the two
mates before a non-split pair counts as a large deletion — read-pair category
(2) of the paper, which the code previously approximated with "the gap overlaps
the protospacer".

`ANCHOR_MIN=13` requires the worst aligned block of a split read to score
`length - 4 x mismatches >= 13`, so a 13 bp anchor must be perfect and a 17 bp
anchor may carry one mismatch. The score is always written as an `XA:i:` tag.

Together they cost 0.31 pp of small-deletion sensitivity on the simulation
(large deletions gain 0.96 pp) and remove 40 % of PRR5's false large deletions
while leaving all 5432 of its real terminal-insertion pairs untouched.

`MIN_GAP` is an **insert size** floor, not the paper's 100 bp large-deletion
definition — a non-split pair only reveals the distance between its mates, and a
~300 bp library produces several hundred bp of that on ordinary WT fragments.
Setting it to 100 quadruples the false-call rate on the simulation and takes
PRR5 from 8.8 % to 40.1 % large deletions. Leave it at roughly
(max insert size) + (smallest deletion worth calling).

`FIXES.md` has the full comparison tables, including the two filters that were
measured and rejected (RF-orientation removal and contradictory-pair removal:
both strip artefacts from the control as well, which is exactly what the
k-means common-cluster step needs the control to keep).

**Building.** htslib is detected automatically; `make` says which way it went.

```bash
sudo apt-get install libhts-dev     # Debian/Ubuntu
make                                # BAM enabled
make HTS=0                          # force SAM-only
make HTS_CFLAGS=-I$HOME/tools/include \
     HTS_LDFLAGS='-L$HOME/tools/lib -Wl,-rpath,$HOME/tools/lib'   # custom prefix
```

Without htslib the same code path still sorts, it just emits SAM — which is
also what the browser build does, since `public/viewer.js` reads SAM text.

Every file carries a real header, so samtools, pysam and IGV read it directly
and `Issue/fix_sam.sh` is only needed for results produced before this change:

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:TRAC	LN:5899
@PG	ID:kmer_align	PN:CRISPR-LargeDel	CL:./kmer_align ...
```

The optional fields are spec-valid `TAG:TYPE:VALUE`. The old four-letter
spellings (`MISS:`, `Dist:`, `D_size:`, `CLIP:`) and the stray empty column
after QUAL are gone:

| tag | meaning |
| --- | --- |
| `XM:i:` | mismatches inside the aligned blocks |
| `XT:i:` | of those, how many fall within the cut site +/-5 bp |
| `XD:i:` | distance from the cut site to the nearest breakpoint |
| `XL:i:` | reference span of the alignment (deletion size for a split read) |
| `XC:i:` | terminal soft clip length, when one qualifies |
| `XS:A:` | which side that clip is on, `L` or `R` in reference orientation |
| `XV:Z:` | the pattern this pair was classified as — same value as the QNAME suffix |
| `YC:Z:` | k-means cluster colour, for viewers that honour it |

### Read names carry the call

Every SAM this tool writes appends the pattern to QNAME, so any record can be
traced back without re-deriving the call:

```
INS_L1_2966_00024|m0_P0015923|f2761|L236|205:AM:J12MnR0/1|INS  99  TRAC  2762  60  84M1I65M
DEL_L5000_401_5401|m0_P0003881|f271|L412|130:AM:JdMnR0/1|LDEL_ONLY  99  TRAC  272  60  132M5000D18M
```

`WT`, `INS`, `DEL`, `SUB`, `CLIP`, `LDEL_ONLY`, `LDEL_COMMON`, plus `NOCOV`
(straddles the cut, neither mate covers it) and `OUT` (nowhere near the cut) for
the pairs that appear only in the full alignment SAM. Both mates always carry
the same suffix, so pairing is unaffected.

### A reference gap needs real anchors on both sides

`side_align()` re-places the first or last aligned block up to 30 bp away to
shave off mismatches, and used to push the result with no length check. A 2-3 bp
coincidental match downstream of the cut was therefore enough to fabricate a
large `D`. Blocks shorter than `MIN_match` (13 bp) at either end are now dropped
before the CIGAR is built, so what is left becomes a soft clip and is picked up
as a terminal insertion instead.

### Terminal insertion (clipped)

An insertion only becomes an `I` in the CIGAR when there is aligned sequence on
both sides of it, because `make_cigar()` emits `I` from the overlap between two
aligned blocks. An insertion at the very edge of a read has aligned sequence on
one side only, so it ends up in a soft clip and the pair used to be reported as
WT. A read whose terminal soft clip is at least **5 bp** long and whose clip
boundary lies within **+/-13 bp** of the cut site is now counted as
`Terminal insertion (clipped)` and written to `_site_{ctr,mut}_clip.sam`, with a
`CLIP:<length>:<L|R>` tag in the optional SAM field (`L` and `R` are the side in
reference orientation). Thresholds are `MIN_CLIP` and `CLIP_WIN` in `align.cpp`.

This category is deliberately kept out of the legacy
`Cas9 Mutation sample small_indel rate` line, so numbers reported before this
change do not move. Two limits are worth knowing:

- The alignment extension allows one mismatch, so it absorbs one or two of the
  inserted bases; an N bp terminal insertion surfaces as roughly an
  `(N-1)S`..`(N-2)S` clip. Insertions shorter than ~6 bp are also frequently
  absorbed as a short `D` by the chunk merging and land in `Small deletion`.
- A terminal clip longer than 20% of the read length fails the
  `r_match < 0.8*query_len` check in `make_cigar()` and the pair is discarded
  entirely (for 150 bp reads that is an insertion of 32 bp or more).

### Small indels at the read edge (clip rescue)

The 13 bp anchor rule above has a cost: a real 1-10 bp deletion that sits within
13 bp of a read end keeps only a 5-12 bp anchor on that side, the anchor is
dropped, and the read becomes `7S143M` - counted as a *terminal insertion*. On
the noPCR simulation (HPRT1, TRAC, BCL11A) that was about 10 % of all 1-10 bp
deletions, which put small-deletion detection *below* the original code
(87.5 % -> 81.9 % for 1 bp).

Whether a short anchor can be trusted is decided from the sample, not from the
read. Every small indel seen with full anchors on both sides is collected while
the sample is aligned; terminal-clip pairs are held back to the end, and a
clipped mate is re-assigned to a known allele when

- that allele was seen in at least 3 pairs with full anchors,
- laid onto the allele, the read has no more mismatches than its aligned part
  already had against the reference (the clipped bases are explained exactly),
- at least 5 read bases sit on each side of the event, and
- no allele of the other kind (insertion vs deletion) explains the read as well
  - that is the PRR5 case, inserted bases that repeat downstream of the cut, and
  such a read stays a clip.

The read is rewritten with a real CIGAR (`7S143M` -> `7M10D143M`) and an
`XR:Z:clip_rescue:D10@4646/n71` tag (allele, 1-based position, supporting
pairs). When two placements of the same insertion (or deletion) explain the
read with the same support, the one with the smaller allele key wins, so the
native and the WebAssembly build pick the same one. The run prints how many
clip pairs were re-assigned. Only the label of
an already mutated pair changes - a WT pair can never become a mutation, and
large deletions are untouched (byte-identical `_mut_only`).

| 3 genes, noPCR | original | before rescue | with rescue |
|---|---|---|---|
| deletion 1-10 bp, right type | 86.6-88.0 % | 81.4-83.4 % | **90.0-91.4 %** |
| insertion 1-3 bp, exact size | 80.2-82.9 % | 80.1-82.2 % | **85.8-87.7 %** |
| stress_edge, deletion 1-10 bp | - | 44.0-46.6 % | **68.9-70.5 %** |
| stress_homology, false deletions on HOM inserts | - | 2.3 / 8.5 % | 2.3 / 8.5 % (no change) |

`--clip_rescue 0` (or `CLIP_RESCUE=0`) restores the previous behaviour.

### Shared with the control? A per-pair breakpoint test (`--common_mode`)

The paper removes large deletions that are also present in the control (PCR
and alignment artefacts) by clustering both samples with k-means on
(start, end) and matching cluster centroids. For deletions of 2 kb and up that
split was not reproducible: k is picked from the grid 1,3,6,9,... by a penalised
mean error, so removing 1 % of the reads, or only shuffling their order, picks
another k or other seeds, and one centroid can swallow real deletions together
with an artefact cloud next to them. On the HPRT1 simulation the share of 5 kb
deletions kept as CRISPR-derived read 58.6 / 46.9 / 43.0 % for the full input,
a 99 % subsample and the same reads shuffled.

The default is now a per-pair test with no clustering:

    tol  = max(20 bp, 2 % of the pair's span)
    n_c  = control pairs with both breakpoints within tol
    n_m  = treated pairs within tol (this one included)
    common  <=>  n_c > 3  and  n_c / N_c >= 0.5 * n_m / N_m

`N_c`, `N_m` are the aligned pairs of each sample: an artefact turns up at a
similar rate in both, a CRISPR deletion is enriched in the treated sample even
with a few artefacts close by. The result depends only on the set of
breakpoints - not on read order, seeds or k.

Why 2 % and not the 5 % of the centroid test: a single pair's breakpoints are
far more precise than a cluster centroid, and 5 % is +/-250 bp at 5 kb. The WT
HPRT1 control carries alignment artefacts (mates mis-placed in a repeat near
both amplicon ends) that land within that window of real 5-7 kb deletions, so
at 5 % the test removed 29 % of real 5 kb deletions; at 2 % it removes 10 %.
The 20 bp floor keeps short deletions from being matched to unrelated
artefacts a few tens of bp away.

| | k-means | per-pair test |
|---|---|---|
| HPRT1, 2 / 5 / 7 kb kept, full vs 99 % vs shuffled | 69.8-88.7 / 43.0-58.6 / 67.8-79.0 % | 87.9-88.0 / 77.6-77.8 / 81.9-82.0 % |
| stress_cluster, shared artefacts wrongly kept (250 / 300 bp) | 66.9 / 43.1 % | 0.5 / 0.1 % |
| stress_cluster, chimeras wrongly kept | 20.0 % | 9.4 % |
| stress_cluster, real 500-2000 bp deletions kept | 45.7-74.4 % | 80.3-83.1 % |
| post-processing time, HPRT1 | ~200 s | 0.2 s |

Viewer colours (`YC:Z:`) now come from a deterministic grouping of nearby
breakpoints. `--common_mode kmeans` (or `COMMON_MODE=kmeans`) restores the
paper's method.

## License

