#!/usr/bin/perl -w

use strict;

my %coverage_seqs;
open my $file, "<", $ARGV[0];
my $sid;
my $cov;
my $max_cov=0;
my $min_cov=100000000;
while(<$file>){
	chomp;
	if(/>/){
		$_ =~ /cov_(.*?)[_|$]/;
		$cov = $1;
		if($max_cov <= $cov){
			$max_cov = $cov;
		}
		if($min_cov >= $cov){
			$min_cov = $cov;
		}
		$sid = $_;
	}
	else{
		$coverage_seqs{$cov}{$sid}.=$_;
	}
}

my $average = ($max_cov+$min_cov)/4;

open my $out, ">", "filtered_spades_contigs.fsa";
for my $tcov (keys %coverage_seqs){
	if($tcov >= $average){
		for my $id (keys %{$coverage_seqs{$tcov}}){
			print $out "$id\n$coverage_seqs{$tcov}{$id}\n";
		}
	}
}
