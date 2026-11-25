#!/bin/bash
# Wrapper around ./kmer_align: fills in the test_data/result paths and sorts
# the resulting SAM files by position.

usage() {
cat <<'HELPTXT'
CRISPR-LargeDel : large deletion / small indel caller for Cas9 amplicon data

USAGE
  bash run_kmer_align.sh <gene> <sample> <target_start> <target_end> \
                         <cleavage_offset> <ctr_prefix> <mut_prefix> <K>

ARGUMENTS
  gene             reference name; reads test_data/<gene>.fa
  sample           sample label; results land in result/<gene>_<sample>_<K>/
  target_start     protospacer start, 0-based offset into the reference
  target_end       protospacer end, must be greater than target_start
  cleavage_offset  cut site measured from target_start (SpCas9, 20 bp
                   protospacer: 17 for a target on the plus strand, 3 on the
                   minus strand)
  ctr_prefix       control FASTQ prefix; reads test_data/<ctr_prefix>_1.fastq
                   and test_data/<ctr_prefix>_2.fastq (or .fastq.gz)
  mut_prefix       Cas9 treated FASTQ prefix, same _1/_2 convention
  K                k-mer length for seeding, 5..31 (11 works for 150 bp reads)

EXAMPLE
  bash run_kmer_align.sh TRAC SIM 2952 2971 17 TRAC_SIM_ctr TRAC_SIM_mut 11
  (the simulated TRAC amplicon in test_data/, the web page's Load example)

OUTPUT
  Coordinate-sorted BAM + .bai for every output, written by kmer_align itself -
  no samtools needed.  Drop them straight into IGV.
    OUT_FMT=sam    plain SAM instead
    OUT_FMT=both   write each output as both

LARGE DELETION FILTERS (environment variables, defaults shown)
  DIST_CUT=8000    max distance (bp) between the cut site and the nearest read
                   end of a large-deletion call. 1000 (the default before
                   2026-10-23) drops most deletions of 3 kb and more.
  MIN_GAP=400      inner gap between mates before a non-split pair counts as a
                   large deletion.  An insert-size floor, NOT the paper's 100 bp
                   deletion-length definition - do not set it to 100.
  ANCHOR_MIN=13    split-read anchor quality floor, (length - 4*mismatches) of
                   the worst aligned block.  0 disables, 16 is stricter.
  CUT_BETWEEN=1    require the cut site to fall between the two mates.
  LARGE_DEL_LEN=50 deletion length inside one read that counts as a large
                   (split-read) deletion. 50 = the original tool; the paper's
                   definition is 100; see --help.
  CLIP_RESCUE=1    re-assign a small deletion/insertion at the very end of a
                   read (seen only as a soft clip) to the same allele seen in
                   full reads of the sample. 0 = behaviour before 2026-10-23.
  COMMON_MODE=local  how a large deletion is judged shared with the control:
                   local = per-pair breakpoint comparison (deterministic),
                   kmeans = the paper's cluster comparison (before 2026-10-23).

For the raw argument list (custom paths, full_sam) run ./kmer_align
with no arguments.
HELPTXT
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  usage
  exit 0
fi
if [ "$#" -lt 8 ]; then
  echo "[argument error] $# argument(s) given, 8 required." >&2
  echo >&2
  usage >&2
  exit 1
fi

genename=$1
cellname=$2
tg_st=$3
tg_ed=$4
cle_pos=$5
ctr_sample=$6
mut_sample=$7
K=$8

base="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
res=$base'/result/'
data_dir=$base'/test_data/'


# kmer_align sorts and indexes its own output now (htslib is linked in), so
# there is no samtools step here any more.  OUT_FMT is passed straight through:
#   bam   (default) coordinate-sorted BAM + .bai only
#   sam             plain SAM
#   both            each output written twice
OUT_FMT=${OUT_FMT:-bam}

# Large deletion calling knobs, see ./kmer_align --help (arguments 16-18).
#   MIN_GAP      inner-gap floor for a read-pair geometry call.  This is an
#                INSERT SIZE floor, not the paper's 100 bp deletion length -
#                setting it to 100 makes ordinary WT pairs look deleted.
#   ANCHOR_MIN   split-read anchor quality floor, (block length - 4*mismatches).
#                13 = a 13 bp anchor must be perfect.  0 turns the test off.
#                16 is worth trying on artefact heavy amplicons.
#   CUT_BETWEEN  1 = the cut site must fall between the two mates (paper
#                category 2).  0 restores the pre-2026 behaviour.
DIST_CUT=${DIST_CUT:-8000}
MIN_GAP=${MIN_GAP:-400}
ANCHOR_MIN=${ANCHOR_MIN:-13}
CUT_BETWEEN=${CUT_BETWEEN:-1}
# Arguments 19-24.  LARGE_DEL_LEN = deletion length (bp) inside one read that makes
# it a split-read large deletion (default 50 = the original tool; 100 = the
# paper's definition).
INS_FIX=${INS_FIX:-1}
INS_MAX=${INS_MAX:-0}
REPORT=${REPORT:-1}
LARGE_DEL_LEN=${LARGE_DEL_LEN:-50}
CLIP_RESCUE=${CLIP_RESCUE:-1}
COMMON_MODE=${COMMON_MODE:-local}

run_align(){
genename=$1
cellname=$2
tg_st=$3
tg_ed=$4
k=$5
dist2cs=$6
ctr=$7
mut=$8
cle_pos=${9}
resdir=$res/$genename''_$cellname''_$k
# <prefix>_1.fastq if it exists, else <prefix>_1.fastq.gz (kmer_align reads both)
fq(){ if [ -f "$data_dir/$1.fastq" ]; then echo "$data_dir/$1.fastq"; else echo "$data_dir/$1.fastq.gz"; fi; }
dataname=$genename''_$cellname
echo $dataname
echo $resdir
mkdir -p $resdir
./kmer_align $data_dir/$genename.fa "$(fq $ctr''_1)" "$(fq $ctr''_2)" "$(fq $mut''_1)" "$(fq $mut''_2)" $tg_st $tg_ed $cle_pos $genename $cellname $k $resdir $dist2cs 1 $OUT_FMT $MIN_GAP $ANCHOR_MIN $CUT_BETWEEN $INS_FIX $INS_MAX $REPORT $LARGE_DEL_LEN $CLIP_RESCUE $COMMON_MODE
# an empty SAM makes the guard above fail; do not let that become the
# script's exit status
return 0
}
run_align $genename $cellname $tg_st $tg_ed $K $DIST_CUT $ctr_sample $mut_sample $cle_pos
