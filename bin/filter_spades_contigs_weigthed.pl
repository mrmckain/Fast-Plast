#!/usr/bin/env perl -w
use strict;

my $total_length;
my $accumulated_coverage;
my %coverage_lengths;
my %cov_seqid;
my %sequences;

open my $file, "<", $ARGV[0];
my $sid;
while(<$file>){
	chomp;
	if(/>/){
		$_ =~ /length_(\d+)_cov_(.*?)(_|$)/;
		$sid=$_;
		my $cov = $2;
		my $len_true=$1;
		$total_length+=$1;
		$coverage_lengths{$cov}=$len_true;
		$accumulated_coverage+=$cov*$len_true;
		$cov_seqid{$cov}=$sid;
	}
	else{
			$sequences{$sid}.=$_;
	}
}

my $weighted_average = $accumulated_coverage/$total_length;

my $variance;
my $std_count;

for my $tcov (keys %coverage_lengths){
	for(my $i=0; $i < $coverage_lengths{$tcov}; $i++){
		my $tempstd = ($tcov-$weighted_average)**2;
		$variance+=$tempstd;
		$std_count++;
	}
}

my $stdev = sqrt($variance/$std_count);

my $min_cov = $weighted_average-$stdev;
my $max_cov	= $weighted_average+2.5*$stdev;

open my $out, ">", "filtered_spades_contigs.fsa";
for my $tcov (sort keys %cov_seqid){
	if($tcov <= $max_cov && $tcov >= $min_cov){
		print $out "$cov_seqid{$tcov}\n$sequences{$cov_seqid{$tcov}}\n";
	}
}
