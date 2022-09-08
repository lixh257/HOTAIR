#!/usr/bin/perl
#use strict;
#use warnings;

############################################
#RNAseq_pipleline_hui
#Contact : 18846044714@163.com
#update:2021-06-24:Lixh
############################################/public/home/jinxu/DB/hg19/

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


my $outdir_Trim="$outdir/Trim";
my $outdir_Map="$outdir/MAp";
my $outdir_Count="$outdir/Count";

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



#1.过滤（Trimmomatics去掉接头和低质量数据）
echo "Trimmomatics"
time=`date`
echo \$time
java -jar /md01/lixh/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33
/md01/lixh/HOTAIR/RNA-seq/OE-YTHDF3-2-1/OE-YTHDF3-2-1_1.fq.gz /md01/lixh/HOTAIR/RNA-seq/OE-YTHDF3-2-1/OE-YTHDF3-2-1_2.fq.gz
/md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_1_paired.fq.gz /md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_1_unpaired.fq.gz
/md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_2_paired.fq.gz /md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_2_unpaired.fq.gz
ILLUMINACLIP:/md01/lixh/software/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
time=`date`
echo \$time

#2.HISAT2进行序列比对
#a.构建索引
###①官网直接下载https://daehwankimlab.github.io/hisat2/download/
###下载对应的基因组文件索引，以GRCh38为例：https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz

##切换到该目录
#$ cd /md01/lixh/HOTAIR/RNA-seq/Ref/
##下载索引文件
#$ wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
##解压得到hg19目录
#$ tar -zxvf *.tar.gz /hg19

##查看
#ls hg19
###得到的文件
 #  (base) [zhouj@alpha new]$ ls grch38
  # genome.1.ht2  genome.2.ht2  genome.3.ht2  genome.4.ht2  genome.5.ht2  genome.6.ht2  genome.7.ht2  genome.8.ht2  make_grch38.sh

#b.序列比对
echo "Mapping"
time=`date`
echo \$time
/md01/lixh/software/hisat2-2.2.1/hisat2 -t -p 8 -x /md01/lixh/HOTAIR/RNA-seq/Ref/hg19/genome \
-1 /md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_1_paired.fq.gz \
-2 /md01/lixh/HOTAIR/RNA-seq/Trim/OE-YTHDF3-2-1_1_paired.fq.gz \
-S /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1.sam

#c.SAM文件转换为BAM文件并sorting
 samtools view -S /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1.sam -b > /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1.bam
# 将所有的bam文件按默认的染色体位置进行排序,并且会自动将多个bam文件merge到一起。merge后的名字为B1-1_L2_164A64_sorted.bam
 samtools sort /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1.bam -o /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1_sorted.bam
time=`date`
echo \$time
 #3.featureCount进行reads计数
echo "Count"
time=`date`
echo \$time
/md01/lixh/software/subread-2.0.2-Linux-x86_64/bin/featureCounts -T 6 -p -t exon -g gene_id -a /md01/lixh/HOTAIR/RNA-seq/Ref/Homo_sapiens.GRCh37.75.gtf.gz -o /md01/lixh/HOTAIR/RNA-seq/Count/OE-YTHDF3-2-1_sorted_counts.txt /md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1_sorted.bam
time=`date`
echo \$time

#4.计算TPM(R)
TPM(/md01/lixh/HOTAIR/RNA-seq/Map/OE-YTHDF3-2-1_sorted_counts.txt)

);


close OUT;

#`chmod +x $script`;

my $cmd ="sh $script 1>$script.log 2>$script.err ";
print $cmd,"\n";
#system($cmd);

}
