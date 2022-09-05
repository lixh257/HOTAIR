###批量计算 FRiP, fraction of reads in called peak regions
###sample_OE-puro-3-2_peaks.narrowPeak
#####sample_OE-puro-3-2.pe.q10.sort.rmdup.bed###这个bed文件是mapping里面的bam之后的bed文件

cd /md01/lixh/HOTAIR/ATAC/output/Peak_calling/
ls *narrowPeak|while  read id;
do
echo $id
bed=$(basename $id "_peaks.narrowPeak").pe.q10.sort.rmdup.bed
summ=$(basename $id "_peaks.narrowPeak")_summits.bed
rud=$(basename $id "_peaks.narrowPeak").filterBL.bed
bam=$(basename $id "_peaks.narrowPeak").bam
#ls  -lh $bed 
totalreads=$(samtools flagstat /md01/lixh/HOTAIR/ATAC/output/Mapping/$bam|awk '{print $1}')
peakReads=$(bedtools intersect -a /md01/lixh/HOTAIR/ATAC/output/Mapping/$bed -b /md01/lixh/HOTAIR/ATAC/output/Peak_calling/$id |wc -l|awk '{print $1}')
rmdupReads=$(wc -l /md01/lixh/HOTAIR/ATAC/output/Mapping/$bed|awk '{print $1}')
allpeaks=$(wc -l /md01/lixh/HOTAIR/ATAC/output/Peak_calling/$summ|awk '{print $1}')
cleanpeaks=$(wc -l /md01/lixh/HOTAIR/ATAC/output/Peak_calling/$srud|awk '{print $1}')
echo $totalreads  $rmdupReads $peakReads $allpeaks $cleanpeaks
echo '==> FRiP value:' $(bc <<< "scale=2;100*$cleanpeaks/$allpeaks")'%'
done

####结果
sample_C-1_peaks.narrowPeak
25356378 49266282
==> FRiP value: 51.46%
sample_C-2_peaks.narrowPeak
31885110 53390668
==> FRiP value: 59.72%
sample_C-DF3-2-1_peaks.narrowPeak
19145428 40169114
==> FRiP value: 47.66%
sample_C-DF3-2-2_peaks.narrowPeak
23113596 39440636
==> FRiP value: 58.60%
sample_C-DF3-4-1_peaks.narrowPeak
9199972 19677786
==> FRiP value: 46.75%
sample_C-DF3-4-2_peaks.narrowPeak
15968607 29953744
==> FRiP value: 53.31%
sample_C-NT-1-1-1_peaks.narrowPeak
20631031 40005126
==> FRiP value: 51.57%
sample_C-NT-1-1-2_peaks.narrowPeak
31357548 51311390
==> FRiP value: 61.11%
sample_C-puro-8-2_peaks.narrowPeak
33891674 55159110
==> FRiP value: 61.44%
sample_OE-1_peaks.narrowPeak
24239654 47096722
==> FRiP value: 51.46%
sample_OE-2_peaks.narrowPeak
43906729 70458014
==> FRiP value: 62.31%
sample_OE-DF3-1-1_peaks.narrowPeak
25097552 46547060
==> FRiP value: 53.91%
sample_OE-DF3-1-2_peaks.narrowPeak
57595110 91627456
==> FRiP value: 62.85%
sample_OE-DF3-2-1_peaks.narrowPeak
42610364 73316006
==> FRiP value: 58.11%
sample_OE-DF3-2-2_peaks.narrowPeak
47722188 72730274
==> FRiP value: 65.61%
sample_OE-DF3-4-1_peaks.narrowPeak
28908616 50259784
==> FRiP value: 57.51%
sample_OE-DF3-4-2_peaks.narrowPeak
32061108 51495898
==> FRiP value: 62.25%
sample_OE-DF3-5-1_peaks.narrowPeak
32448542 54802138
==> FRiP value: 59.21%
sample_OE-DF3-5-2_peaks.narrowPeak
30742024 52970780
==> FRiP value: 58.03%
sample_OE-puro-3-1_peaks.narrowPeak
25283918 45325342
==> FRiP value: 55.78%
sample_OE-puro-3-2_peaks.narrowPeak
31953614 52934640
==> FRiP value: 60.36%
###peak的数量结果不大对，有点奇怪
wc -l sample_OE-puro-3-2_summits.bed
130852 sample_OE-puro-3-2_summits.bed
wc -l sample_OE-puro-3-2.filterBL.bed
130708 sample_OE-puro-3-2.filterBL.bed
samtools view -c -L /md01/lixh/HOTAIR/ATAC/output/Peak_calling/sample_OE-puro-3-2_peaks.narrowPeak -F 0x0100 sample_OE-puro-3-2.pe.q10.sort.rmdup.bam
samtools view -c -L /md01/lixh/HOTAIR/ATAC/output/Peak_calling/sample_OE-puro-3-2_peaks.narrowPeak -F 0x0100 sample_OE-puro-3-2.pe.q10.sort.rmdup.bam
32833368
samtools flagstat sample_OE-puro-3-2.pe.q10.sort.rmdup.bam
52934640
FRip=32833368/52934640=0.6202
 

bedtools intersect -a /md01/lixh/HOTAIR/ATAC/output/Mapping/sample_OE-puro-3-2.pe.q10.sort.rmdup.bed -b /md01/lixh/HOTAIR/ATAC/output/Peak_calling/sample_OE-puro-3-2_peaks.narrowPeak |wc -l
31953614
wc -l /md01/lixh/HOTAIR/ATAC/output/Mapping/sample_OE-puro-3-2.pe.q10.sort.rmdup.bed
52934640
FRip=31953614/52934640=0.603
