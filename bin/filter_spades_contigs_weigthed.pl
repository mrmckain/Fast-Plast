#!/usr/bin/env perl
use strict;
use warnings;

my $total_length;
my $accumulated_coverage;
my %coverage_lengths;
my %cov_seqid;
my %sequences;
my $forced_minimum;
if($ARGV[1]){	
	my $forced_minimum=$ARGV[1]; #Optional
}
open my $file, "<", $ARGV[0];
my $sid;
while(<$file>){
	chomp;
	if(/>/){
		$_ =~ /length_(\d+)_cov_(.*?)(_|$)/;
		$sid=$_;
		my $cov = $2;
		my $len_true=$1;
		if ($1 < 1000){
			next;
		}
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

my $min_cov;
if($forced_minimum){
	$min_cov = $forced_minimum;
 }
else{
	$min_cov = $weighted_average-$stdev;
}
my $max_cov	= $weighted_average+2.5*$stdev; #Can add a factor here to change the amount of STDEV allowed for coverage.

open my $out, ">", "filtered_spades_contigs.fsa";
for my $tcov (sort keys %cov_seqid){
    if(scalar keys %sequences >=4){
        if($tcov <= $max_cov && $tcov >= $min_cov){
            print $out "$cov_seqid{$tcov}\n$sequences{$cov_seqid{$tcov}}\n";
        }
    }
    else{
        print $out "$cov_seqid{$tcov}\n$sequences{$cov_seqid{$tcov}}\n";
    }
}
