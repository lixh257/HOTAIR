# HOTAIR
YTHDF3 &amp; HOTAIR
############################################
#ATACseq_pipleline_v4.0.2


############################################
Major steps to run ATAC-seq pipeline v4.0.2

1. Set up a configuration.txt file
(See an example of this step in /seq/ATAC-seq/Example/configuration.txt.)
** No empty row allowed in configure

# for single sample
ref_index=/path/to/bowite2 index/ 
ref_size=/path/to/genome_size
max_thread=xx (the maximun thread per sample, default is 8)
gsize=hs/mm/ce/dm (-g parameter in macs2, and species code to select blacklist)
/path/to/sample1_R1.fastq	/path/to/sample1_R2.fastq	/path/to/output	sample_prefix

# for multiple sample
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
(The shell script once created, will be located at the same directory as the configuration.txt
See an example of this step in /seq/ATAC-seq/Example/command.sh.)
** Max_thread is calculated as 10 * number of sample
 
perl /seq/ATAC-seq/Code/ATACseq_pipeline_v4.0.2.pl /path/to/configuration.txt
##2021-05-26#perl /md01/lixh/HOTAIR/ATAC/Code/ATACseq_pipeline_hui.pl /md01/lixh/HOTAIR/ATAC/Code/configure.txt###

(2) submit shell scripts
** Make sure that you ONLY submit ONE job at a time. DO NOT simultaneously
submit multiple jobs!!

You can run shell script using full path 
sh /path/to/sample_prefix.sh 1> /path/to/sample_prefix.sh.log 2>/path/to/sample_prefix.sh.err 

OR you can use relative path
sh ./sample_prefix.sh 1> ./sample_prefix.log 2> ./sample_prefix.err 


Output files
Three folders, Mapping, QC, and Peak_calling, will be created for each sample.

Mapping/
- sample_1M.bam (chr1-22,X,Y bam)
- sample_1M_Rep2.pe.q10.sort.bam (MAPQ>10)
- sample_1M_Rep2.pe.q10.sort.rmdup.bam (MAPQ>10 and PCR duplicates removed)
- sample_1M_Rep2.pe.q10.sort.rmdup.bam.bai (bam index)
- sample_1M_Rep2.pe.q10.sort.rmdup.bed (shift 4bp on + strand; 5bp on -strand; expand 25bp l/r from the shifted sites)
- sample_1M_Rep2.chrM.bam (chrM bam)
- sample_1M_Rep2.bedGraph (bedGraph from sample_1M_Rep2.pe.q10.sort.rmdup.bed)
- sample_1M_Rep2.norm.bedGraph (normalized bedGraph)
- sample_1M_Rep2.norm.bw (normalized bw) 

QC/
- sample_1M_Rep2.Picard_Metrics_unfiltered_bam.txt (Picard output)
- Preseq1.1.2_sample_1M_Rep2.pe.q10.sort.bam.txt (Preseq output)
- sample_1M_Rep2.pe.q10.sort.rmdup.bam.hg19_refseq_genes_TSS.txt.vect (ATAC-seq signal around TSS)
- sample_1M_Rep2.pe.q10.sort.rmdup.bam.hg19_refseq_genes_TSS.txt.vect.pdf (ATAC-seq signal around TSS)
- sample_1M_Rep2.4kb.hg19TSS.txt (ATAC-seq signal around +/- 2kb of TSS)
- sample_1M_Rep2.fragment_length_distribution.pdf (read length distribution plot)
- sample_1M_Rep2.fragment_length_distribution.txt (read length distribution statistics)
- sample_1M_Rep2.samtools.flagstat.txt (samtools flagstat output)
- sample_1M_Rep2.samtools.idxstats.txt (samtools idxstat output)

Peak_calling/
- sample_1M_Rep2_peaks.narrowPeak (macs2 output)
- sample_1M_Rep2_peaks.xls (macs2 output)
- sample_1M_Rep2_summits.bed (macs2 output)
- sample_1M_Rep2.filterBL.bed (peak list after filtering peaks in 1.ENCODE blacklist regions made by kundaje lab, and 2.ucsc NumtS regions)

3 Generate QC table
   (1) First make sure that you finished running all the samples listed in configure file successfully. 
   (2) Run the following script which gives default quality control output file named qcTable.txt.
           python /seq/ATAC-seq/Code/getQCTable.py -i /path/to/configuration.txt
       If you want to give your output file a specific name, you can add -o to add your preferred file name
           python /seq/ATAC-seq/Code/getQCTable.py -i /path/to/configuration.txt -o your_preferred_name
##2021-05-26# python /md01/lixh/HOTAIR/ATAC/Code/AgetQCTable.py /md01/lixh/HOTAIR/ATAC/Code/configure.txt###
