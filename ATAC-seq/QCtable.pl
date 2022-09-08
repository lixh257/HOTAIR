#!/usr/bin/perl
# get QC for mutiple samples 
# The last step for the ATACseq pipeline, which based on the previous output from the pipeline
# Author: xujin937@gmail.com
# Date : 2015-12-31
#
use strict;
use warnings;
if(@ARGV<1)
{
	print "Usage perl $0 <configure file> \n";
	exit;
}

# print table head 
print join("\t",("sample","Total number of reads","number of trimmed reads","mapping rate","Q10 reads","chrM reads","final reads","median fragment length","TSS enrichment")),"\n";

open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";
while(<IN>)
{
my @a=split; 
if($#a>1)
{
my $dir=$a[2];
my $file=$a[3];
my $log_file=$a[2]."/".$a[3].".sh.err";
#print $log_file,"\n";
open LOG,"$log_file" or die "can not open $log_file\n";
my $total_reads=0;
my $total_trimmed=0;
my $mapping_rate=0;

while(<LOG>)
{

	if(/TOTAL number of reads = (\d+)/)
	{
		$total_reads=$1*2;	
	}
	if(/Total number of trimmed reads = (\d+)/)
	{
		$total_trimmed=$1;
	}
	if(/(\d+\.\d+)% overall alignment rate/)
	{
		$mapping_rate=$1;
	}
}	
close LOG;

print $a[2],"\t",$total_reads,"\t",$total_trimmed,"\t",$mapping_rate,"\t";

 $log_file=$a[2]."/".$a[3].".sh.log";
#print $log_file,"\n";
open LOG,"$log_file" or die "can not open $log_file\n";
my $total_q10=0;
my $total_chrM=0;
my $total_rmdup=0;

while(<LOG>)
{
if(/count of total reads after QC filter/)
{
	$total_q10=<LOG>;
	chomp($total_q10);
}	
if(/count of chrM reads/)
{
	$total_chrM=<LOG>;
	chomp($total_chrM);
}
if(/final mapped reads/)
{
	$_=<LOG>;
	my @a=split;
	$total_rmdup=$a[0];
}
}
close LOG;

print $total_q10,"\t",$total_q10/$total_reads,"\t",$total_chrM,"\t",$total_chrM/$total_reads,"\t",$total_rmdup,"\t",$total_rmdup/$total_reads,"\t";

$log_file=$a[2]."/"."QC/".$a[3].".fragment_length_distribution.txt";
open LOG,"$log_file" or die "can not open $log_file\n";
while(<LOG>)
{
	my $_=<LOG>;
	my @b=split;
	print $b[2],"\t";
}

close LOG;


$log_file=$a[2]."/"."QC/".$a[3].".TSSenrich";
open OUT, ">Score.R" or die "can not open \n";
print OUT qq(
a<-read.table(\"$log_file\")
p=max(a\$V1)/mean(a[1:100,1])
print(p)
);
system("Rscript Score.R >tmp");
open LOG, "tmp" or die "can not open tmp\n";
while(<LOG>)
{
	my @a=split;	
	if($#a>=1)
	{
		print  $a[1],"\n";
	}
}



}
}
