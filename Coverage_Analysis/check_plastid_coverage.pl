#!/usr/bin/perl -w
use strict;

####Script to take coverage and identify stretches that need to be fixed###

####Usage: 1-Coverage file 2-Kmer Size 3-Coverage Min

my $kmer=$ARGV[1];
my $mincov;
if($ARGV[2]){
	$mincov=$ARGV[2];
}
my $exp_id;
$ARGV[0] =~ /(.*?).coverage_25kmer.txt/;
$exp_id=$1;

open my $file, "<", $ARGV[0];
my $total_cov;
my $total_windows;
while(<$file>){
        chomp;
        my @tarray=split/\s+/;
	$total_cov +=$tarray[2];
	$total_windows++;
}

my $avg_cov = $total_cov/$total_windows;
unless($mincov){
	$mincov = $avg_cov*0.25;
}
print "$mincov\n";
my $current_kmer=0;
my $current_kmer_start;
my $current_kmer_end;


open my $out, ">", $exp_id  . "_problem_regions_plastid_assembly.txt";
open $file, "<", $ARGV[0];
my $run_cov;
while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	
	if($tarray[2] <= $mincov && $tarray[0] !~ /N/i){
		if($current_kmer_start){
			$current_kmer_end=$tarray[1];
			$current_kmer++;
			$run_cov+=$tarray[2];
		}
		else{
			$current_kmer_start=$tarray[1];
			$current_kmer_end=$tarray[1];
			$current_kmer++;
			$run_cov+=$tarray[2];
		}
	}
	
	else{
		if($current_kmer >= $kmer){
			my $run_avg = $run_cov/($current_kmer_end-$current_kmer_start);	
			print $out "$current_kmer_start\t$current_kmer_end\t$run_avg\n";
			$current_kmer=0;
			$current_kmer_end=();
			$current_kmer_start=();
			$run_cov=0;
		}
		else{
			$current_kmer=0;
			$current_kmer_end=();
			$current_kmer_start=();
			$run_cov=0;
		}
	}
}


