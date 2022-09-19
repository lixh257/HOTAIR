# HOTAIR
YTHDF3 &amp; HOTAIR
######RNA-seq


Major steps to run RNA-seq pipeline 

1. Trimming 

#use trimmomatic-0.38.jar

#see an example of this step in RNA-Trim.pbs

qsub RNA-Trim.pbs

2. Mapping #hisat2-2.2.1

#download the software in https://daehwankimlab.github.io/hisat2/download/

#download the reference

wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz

tar -zxvf *.tar.gz /hg19

##see an example of this step in RNA-hisat.pbs

qsub RNA-hisat.pbs

3. sorting

#converted the SAM files to BAM files and sorting them

##see an example of this step in RNA-bam-sort.pbs

qsub RNA-bam-sort.pbs

4. Count reads 

#subread-2.0.2-Linux-x86_64/featureCount

##see an example of this step in RNA-featureCounts.pbs

qsub RNA-featureCounts.pbs

5. The downstream analysis

#see the code of  this step in RNA-downstream analysis.R




######ATAC-seq


Major steps to run ATAC-seq pipeline 

1. Set up a configuration.txt file

#See an example of this step in ATAC-seq/configure.txt

** No empty row allowed in configure


ref_index=/path/to/bowite2 index/ 
ref_size=/path/to/genome_size
max_thread=xx (the maximun thread per sample, default is 8)
gsize=hs/mm/ce/dm (-g parameter in macs2, and species code to select blacklist)
/path/to/sample1_R1.fastq	/path/to/sample1_R2.fastq	/path/to/output	sample_prefix


ref_index=/path/to/bowite2 index/ 
ref_size=/path/to/genome_size
max_thread=xx (the maximun thread for all samples. Can be calculated as 8 * number of samples, default is 8.)
gsize=hs/mm/ce/dm (-g parameter in macs2, and species code to select blacklist)
/ke sure that you ONLY submit ONE job at a time. DO NOT simultaneously
submit multiple jobs!!path/to/sample1_R1.fastq	/path/to/sample1_R2.fastq	/path/to/output	sample1
/path/to/sample2_R1.fastq	/path/to/sample2_R2.fastq	/path/to/output	sample2
/path/to/sample3_R1.fastq	/path/to/sample3_R2.fastq	/path/to/output	sample3


2. Run ATAC-seq pipeline

(1) Create a shell script for each sample 
(The shell script once created, will be located at the same directory as the configure.txt)
** Max_thread is calculated as 10 * number of sample
 
perl /md01/lixh/HOTAIR/ATAC/Code/ATACseq_pipeline_hui.pl /md01/lixh/HOTAIR/ATAC/Code/configure.txt

(2) submit shell scripts

** Make sure that you ONLY submit ONE job at a time. DO NOT simultaneously
submit multiple jobs!!

You can run shell script using full path 
sh /path/to/sample_prefix.sh 1> /path/to/sample_prefix.sh.log 2>/path/to/sample_prefix.sh.err 
#sh /md01/lixh/HOTAIR/ATAC/output/sample_OE-DF3-5-2.sh 1>/md01/lixh/HOTAIR/ATAC/output/sample_OE-DF3-5-2.sh.log 2>/md01/lixh/HOTAIR/ATAC/output/sample_OE-DF3-5-2.sh.err

Output files
Three folders, Mapping, QC, and Peak_calling, will be created for each sample.

3 Generate QC table

   (1) First make sure that you finished running all the samples listed in configure file successfully. 
   (2) Run the following script which gives default quality control output file named qcTable.txt.
 python /md01/lixh/HOTAIR/ATAC/Code/AgetQCTable.py /md01/lixh/HOTAIR/ATAC/Code/configure.txt
 
 4 The downstream analysis
 
#see the code of  this step in ATAC-downstream analysis.R
   
