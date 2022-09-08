#!/usr/bin/perl
#use strict;
#use warnings;

############################################
#ATACseq_pipleline_Jin
#Contact : xujin937@gmail.com
#update: 2015-07-08
#update: 2016-04-18: Add TSS enrichment profile 
#upadte: 2018-02-15: Correction TSS enrichment plot to strand specific and pdf format
#update: 2018-04-04: Clean up the comments
#update:2021-05-26:Lixh
############################################

#####	Parameters	###

if(@ARGV<1)
{
	print "Usage perl $0 <configure.txt>\n";
	print "congfigure.txt example: \n";
	print "\tref_index=path	\n";
	print "\tref_size=path	\n";
	print "\tmax_thread=n\n";
	print "\tgsize=hs/mm/ce/dm\n";
	print "\tRead1\tRead2\t ouput dir1\t outdir/$output file prefix1 \n";
	print "\tRead1\tRead2\t ouput dir2\t outdir/$output file prefix2 \n";
	print "\t.\n";
	print "\t.\n";
	print "\t.\n";

	exit;
}
open CON, "$ARGV[0]" or die "can not open $ARGV[0]\n";

my $total_thread=1;
my $ref;
my $ref_size;
my $gsize;
my @samples;

while(<CON>)
{
chomp;
#print $_,"\n";
if(/^ref_index=(\S+)/)
{
	$ref=$1;
	#print $ref,"\n";
}
elsif(/^ref_size=(\S+)/)
{
	$ref_size=$1;
	#print $ref_size,"\n";
}
elsif(/^max_thread=(\d+)/)
	{
		$total_thread=$1;
		#print $total_thread,"\n";
			
	}
elsif(/^gsize=(\w+)/)
	{
		$gsize=$1;
	}
else
{
#	print $_,"\n";
	push @sample, $_;
}

}


foreach my $item(@sample)
{
chomp($item);
my @a=split(/\s+/,$item);
#print join("\t",@a),"\n";
my $file1=$a[0];
$file1 =~ /(.*)\./;
my $file1Base = $1;
$file1Base =~ s/.gz//g;
$file1Base =~ s/.fastq//g;
my $file2=$a[1];
$file2 =~ /(.*)\./;
my $file2Base = $1;
$file2Base =~ s/.gz//g;
$file2Base =~ s/.fastq//g;

my $outdir=$a[2];
#print $outdir,"\n";
my $output=$a[3];
my $thread=int($total_thread/($#sample+1));
if($thread<=0)
{
	print "Error message : Too much jobs,no enough resource \n";
	exit;
}

#print $#sample,"\n";
if(! -d $outdir)
{
mkdir $outdir
}


my $outdir_map="$outdir/Mapping";
my $outdir_qc="$outdir/QC";
my $outdir_peak="$outdir/Peak_calling";

if(! -d $outdir_map)
{
mkdir $outdir_map
}

if(! -d $outdir_qc)
{
mkdir $outdir_qc
}

if(! -d $outdir_peak)
{
mkdir $outdir_peak
}

my $script=$outdir."/".$output.".sh";
#print $script,"\n";
open OUT, ">$script" or die "can not open $script \n";
my $peakFile = $output . "_peaks.narrowPeak";


print OUT   qq(
#####	Mapping	#####

#PBS -N $output
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=12
#PBS -l mem=1000m

source /etc/profile.d/modules.sh
module load bowtie2/2.1.0 
module load samtools


# 1. Adapter Trimming
echo "trimming adaptor"
time=`date`
echo \$time
/md01/lixh/HOTAIR/ATAC/Code/adapterTrimmingModified $file1 $file2
echo "trimming adaptor finished "
time=`date`
echo \$time

# 2. Mapping
echo "Mapping by bowtie2"
time=`date`
echo \$time
/public/software/bowtie2-2.1.0/bowtie2 -p $thread --very-sensitive -x $ref -1 $file1Base.trim.fastq -2 $file2Base.trim.fastq  -S  $outdir_map/$output.sam 
awk '\$3!="chrM"' $outdir_map/$output.sam |samtools view -S -b -f 0x2 -q 10 - |samtools sort - -o $outdir_map/$output.pe.q10.sort.bam
samtools view -Sb $outdir_map/$output.sam  > $outdir_map/$output.bam
#rm $file1Base.trim.fastq
#rm $file2Base.trim.fastq
time=`date`
echo \$time

# 3. ChrM 
echo "ChrM"
time=`date`
echo \$time
awk '\$3=="chrM" || NF<10 ' $outdir_map/$output.sam |samtools view -S -b -  > $outdir_map/$output.chrM.bam
#rm $outdir_map/$output.sam
echo "ChrM Mapping by bowtie2 finished"
time=`date`
echo \$time

# 4. Picard (duplicate removal) 
echo "Duplicate removal by Picard"
time=`date`
echo \$time
java -Xmx4G -jar /public/home/jinxu/software/picard-tools-2.4.1/picard.jar MarkDuplicates  INPUT=$outdir_map/$output.pe.q10.sort.bam OUTPUT=$outdir_map/$output.pe.q10.sort.rmdup.bam METRICS_FILE=$outdir_qc/$output.Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &> $outdir_qc/$output.Picard.log
echo "Duplicate removal finished"
time=`date`
echo \$time


# 5. bedgraph, normalized bedgraph, normalized bw 
echo "convert bam 2 bed file , shift to cleavage sites, bedGraph"
/md01/lixh/HOTAIR/ATAC/Code/bam2bed_shift.pl  $outdir_map/$output.pe.q10.sort.rmdup.bam 


echo "make bedGraph"
bedtools  genomecov  -bg -split -i $outdir_map/$output.pe.q10.sort.rmdup.bed  -g $ref_size > $outdir_map/$output.bedGraph
echo "make normalized bedGraph"
/md01/lixh/HOTAIR/ATAC/Code/norm_bedGraph.pl $outdir_map/$output.bedGraph $outdir_map/$output.norm.bedGraph &>  $outdir_map/$output.norm.bedGraph.log
#bedGraphToBigWig $outdir_map/$output.norm.bedGraph $ref_size  $outdir_map/$output.norm.bw 
#rm $outdir_map/$output.norm.bedGraph.log
echo "finish bam2bed , bedGraph, normalized bedGraph and normalized bigwig"
time=`date`
echo \$time


echo "count of total reads after QC filter"
bedtools bamtobed -i $outdir_map/$output.pe.q10.sort.bam | wc -l

echo "count of chrM reads"
bedtools bamtobed -i $outdir_map/$output.chrM.bam | wc -l

echo "final mapped reads"
wc -l $outdir_map/$output.pe.q10.sort.rmdup.bed

echo "count of total reads before QC filter"
bedtools bamtobed -i $outdir_map/$output.bam | wc -l


#####	QC	#####

# 6. Preseq (library complexity) 
echo "Library complexity estimation by Preseq"
time=`date`
echo \$time
#/seq/preseq/preseq lc_extrap -B -P -o $outdir_qc/Preseq1.1.2_$output.pe.q10.sort.bam.txt $outdir_map/$output.pe.q10.sort.bam 
echo "Library complexity estimation finished"
time=`date`
echo \$time

# 7. Count reads in black list region
echo "Count the total number of reads in black list region"
time=`date`
echo \$time
echo "Total number of reads in the black list"
bedtools intersect -a $outdir_map/$output.pe.q10.sort.rmdup.bed -b /md01/lixh/HOTAIR/ATAC/Reference/hg19_kundaje_blacklist.bed  -u | wc -l
echo "Finished counting the total number of reads in black list region"
time=`date`
echo \$time

# 8. Enrichment score 
if [[ $gsize == "hs" ]]
then
echo "Enrichment score calculation"
time=`date`
echo \$time
samtools index $outdir_map/$output.pe.q10.sort.rmdup.bam
python /md01/lixh/HOTAIR/ATAC/Code/pyMakeVplot.py -a $outdir_map/$output.pe.q10.sort.rmdup.bam -b  /md01/lixh/HOTAIR/ATAC/Reference/hg19_refseq_genes_TSS.bed.filter  -p ends -e 2000 -u -v -o $outdir_qc/$output.TSSenrich
samtools idxstats $outdir_map/$output.pe.q10.sort.rmdup.bam > $outdir_qc/$output.samtools.idxstats.txt
samtools flagstat $outdir_map/$output.pe.q10.sort.rmdup.bam > $outdir_qc/$output.samtools.flagstat.txt
echo "Number of reads in +/- 2kb window around TSS" >> $outdir_qc/$output.4kb.hg19TSS.txt 
bedtools intersect -a $outdir_map/$output.pe.q10.sort.rmdup.bed -b /md01/lixh/HOTAIR/ATAC/Reference/4kbwindow_hg19TSS.txt -wa -u |wc -l >> $outdir_qc/$output.4kb.hg19TSS.txt
echo "Enrichment score calculation finished"
time=`date`
echo \$time
fi

if [[ $gsize == "mm" ]]
then
echo "Enrichment score calculation"
time=`date`
echo \$time
samtools index $outdir_map/$output.pe.q10.sort.rmdup.bam
python /md01/lixh/HOTAIR/ATAC/Code/pyMakeVplot.py -a $outdir_map/$output.pe.q10.sort.rmdup.bam -b /md01/lixh/HOTAIR/ATAC/Reference/mm9_refseq_genes_TSS.txt -p ends -e 2000 -u -v  -o  $outdir_qc/$output.TSSenrich
samtools idxstats $outdir_map/$output.pe.q10.sort.rmdup.bam > $outdir_qc/$output.samtools.idxstats.txt
samtools flagstat $outdir_map/$output.pe.q10.sort.rmdup.bam > $outdir_qc/$output.samtools.flagstat.txt
echo "Number of reads in +/- 2kb window around TSS" >> $outdir_qc/$output.4kb.mm9TSS.txt
bedtools intersect -a $outdir_map/$output.pe.q10.sort.rmdup.bed -b /md01/lixh/HOTAIR/ATAC/Reference/4kbwindow_mm9TSS.txt -wa -u |wc -l >> $outdir_qc/$output.4kb.mm9TSS.txt
echo "Enrichment score calculation finished"
time=`date`
echo \$time
fi

# 9. Read length distribution 
echo "Read length distribution"
time=`date`
echo \$time
perl /md01/lixh/HOTAIR/ATAC/Code/fragment_length_dist.pl $outdir_map/$output.pe.q10.sort.rmdup.bam $outdir_qc/$output.fragL.txt
sort -n $outdir_qc/$output.fragL.txt | uniq -c > $outdir_qc/$output.frag.sort.txt
Rscript /md01/lixh/HOTAIR/ATAC/Code/fragment_length_dist.R $outdir_qc/$output.fragL.txt $outdir_qc/$output.frag.sort.txt $outdir_qc/$output.fragment_length_distribution.pdf $outdir_qc/$output.fragment_length_distribution.txt
#rm $outdir_qc/$output.fragL.txt
#rm $outdir_qc/$output.frag.sort.txt
echo "Read length distribution finished"
time=`date`
echo \$time


#####	Peak Calling	#####

# 10. peak calling
echo "peak calling "
time=`date`
echo \$time
macs2  callpeak -t  $outdir_map/$output.pe.q10.sort.rmdup.bed -f BED  -g $gsize --outdir $outdir_peak  -q 0.01 -n $outdir_peak/$output --nomodel  --shift 0  
echo "finish peak calling by Macs2"
time=`date`
echo \$time

# 11. Remove peak in the black list region
if [[ $gsize == "hs" ]]
then
echo "Remove peaks in the black list region"
time=`date`
echo \$time
bedtools intersect -a $outdir_peak/$peakFile -b /md01/lixh/HOTAIR/ATAC/Reference/hg19_kundaje_blacklist.bed  -v > $outdir_peak/$output$filename.filterBL.bed
echo "finish filtering peaks in black list region"
time=`date`
echo \$time
fi

if [[ $gsize == "mm" ]]
then
echo "Remove peaks in the black list region"
time=`date`
echo \$time
bedtools intersect -a $outdir_peak/$peakFile -b /md01/lixh/HOTAIR/ATAC/Reference/mm9_kundaje_blacklist.bed  -v > $outdir_peak/$output$filename.filterBL.bed
echo "finish filtering peaks in black list region"
time=`date`
echo \$time
fi
);


close OUT;

#`chmod +x $script`;

my $cmd ="sh $script 1>$script.log 2>$script.err ";
print $cmd,"\n";
#system($cmd);

}
