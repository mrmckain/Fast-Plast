#!/usr/bin/env perl
use strict;
use warnings;

# Filter SPAdes contigs by k-mer coverage.
#
# Usage: filter_spades_contigs_weigthed.pl <contigs.fasta> [min_coverage]
#
# Contigs shorter than 1000 bp are ignored. Among the rest, the length-weighted
# median coverage is taken as single-copy depth, and contigs outside
# [0.25 x median, 5 x median] are dropped. If min_coverage is given it
# replaces the lower bound. Fewer than four contigs are always passed through.
# Output: filtered_spades_contigs.fsa in the current directory.

my ($infile, $forced_minimum) = @ARGV;
defined $infile or die "Usage: $0 <contigs.fasta> [min_coverage]\n";

my $total_length = 0;
my $accumulated_coverage = 0;
my %cov_of;        # seqid -> coverage   (keyed by id: two contigs with the same
my %len_of;        # seqid -> length      coverage value must not collapse)
my %sequences;
my @order;

open my $file, "<", $infile or die "ERROR: cannot open $infile: $!\n";
my $sid;
while(<$file>){
	chomp;
	if(/^>/){
		if(/length_(\d+)_cov_([0-9.eE+-]+)/){
			my ($len, $cov) = ($1, $2);
			if($len < 1000){
				$sid = undef;
				next;
			}
			$sid = $_;
			push @order, $sid;
			$len_of{$sid} = $len;
			$cov_of{$sid} = $cov;
			$total_length += $len;
			$accumulated_coverage += $cov * $len;
		}
		else{
			warn "WARNING: skipping header without SPAdes length/cov fields: $_\n";
			$sid = undef;
		}
	}
	elsif(defined $sid){
		$sequences{$sid} .= $_;
	}
}
close $file;

open my $out, ">", "filtered_spades_contigs.fsa" or die "ERROR: cannot write filtered_spades_contigs.fsa: $!\n";

if(!@order){
	warn "WARNING: no SPAdes contigs >= 1000 bp; nothing to filter.\n";
	close $out;
	exit 0;
}

# Baseline = length-weighted MEDIAN coverage. The single-copy regions are about
# three quarters of a plastome, so the median sits at single-copy depth, with
# the IR near 2x and its rRNA stretch often at 3-4x. The old bounds, mean - 1 sd to mean + 2.5 sd,
# were pulled around by that spread: they dropped the rRNA part of the IR for
# being too high and single-copy contigs for being too low, which left holes
# afin could not close. Keep 0.25x to 5x the median instead; mitochondrial and
# nuclear plastid-like contigs sit well below 0.25x.
my $LOW_FACTOR  = 0.25;
my $HIGH_FACTOR = 5;
my $median = $cov_of{$order[0]};
{
	my $acc = 0;
	for my $id (sort { $cov_of{$a} <=> $cov_of{$b} } @order){
		$acc += $len_of{$id};
		if($acc >= $total_length / 2){ $median = $cov_of{$id}; last; }
	}
}
my $min_cov = defined $forced_minimum ? $forced_minimum : $LOW_FACTOR * $median;
my $max_cov = $HIGH_FACTOR * $median;

for my $id (@order){
	if(@order >= 4){
		next unless $cov_of{$id} <= $max_cov && $cov_of{$id} >= $min_cov;
	}
	print $out "$id\n$sequences{$id}\n";
}
close $out;
