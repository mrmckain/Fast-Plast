#!/usr/bin/perl -w
use strict;

####Script to take coverage and identify stretches that need to be fixed###

####Usage: 1-Coverage file 2-Kmer Size 3-Coverage Min

my $kmer=$ARGV[1];
my $mincov=$ARGV[2];
my $exp_id;
$ARGV[0] =~ /(.*?).coverage_25kmer.txt/;
$exp_id=$1;

my $current_kmer=0;
my $current_kmer_start;
my $current_kmer_end;


open my $out, ">", $exp_id  . "_problem_regions_plastid_assembly.txt";
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	
	if($tarray[2] <= $mincov){
		if($current_kmer_start){
			$current_kmer_end=$tarray[1];
			$current_kmer++;
		}
		else{
			$current_kmer_start=$tarray[1];
			$current_kmer_end=$tarray[1];
			$current_kmer++;
		}
	}
	
	else{
		if($current_kmer >= $kmer){	
			print $out "$current_kmer_start\t$current_kmer_end\n";
			$current_kmer=0;
			$current_kmer_end=();
			$current_kmer_start=();
		}
		else{
			$current_kmer=0;
			$current_kmer_end=();
			$current_kmer_start=();
		}
	}
}


