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
dataname=$genename''_$cellname
echo $dataname
echo $resdir
mkdir $resdir
./kmer_align $data_dir/$genename.fa $data_dir/$ctr''_1.fastq $data_dir/$ctr''_2.fastq $data_dir/$mut''_1.fastq $data_dir/$mut''_2.fastq $tg_st $tg_ed $cle_pos $genename $cellname $k $resdir $dist2cs #> $resdir/$dataname''_log
sort -k4n $resdir/$dataname''_ctr_only.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_ctr_only.sam
sort -k4n $resdir/$dataname''_ctr_common.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_ctr_common.sam
sort -k4n $resdir/$dataname''_mut_only.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_mut_only.sam
sort -k4n $resdir/$dataname''_mut_common.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_mut_common.sam
sort -k4n $resdir/$dataname''_mut.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_mut.sam
sort -k4n $resdir/$dataname''_ctr.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_ctr.sam
sort -k4n $resdir/$dataname''_site_ctr.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_site_ctr.sam
sort -k4n $resdir/$dataname''_site_mut.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_site_mut.sam
sort -k4n $resdir/$dataname''_site_ctr_sm.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_site_ctr_sm.sam
sort -k4n $resdir/$dataname''_site_mut_sm.sam > $resdir/t
mv $resdir/t $resdir/$dataname''_site_mut_sm.sam
}
run_align $genename $cellname $tg_st $tg_ed $K 1000 $ctr_sample $mut_sample $cle_pos
