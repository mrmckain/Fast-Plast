#!/usr/bin/env perl
use strict;
use warnings;

# Filter SPAdes contigs by k-mer coverage.
#
# Usage: filter_spades_contigs_weigthed.pl <contigs.fasta> [min_coverage]
#
# Contigs shorter than 1000 bp are ignored. Among the rest, a length-weighted
# mean and standard deviation of coverage are computed, and contigs outside
# [mean - 1 sd, mean + 2.5 sd] are dropped. If min_coverage is given it
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

my $weighted_average = $accumulated_coverage / $total_length;

# Length-weighted variance (each contig contributes once per base).
my $variance = 0;
for my $id (@order){
	$variance += $len_of{$id} * ($cov_of{$id} - $weighted_average) ** 2;
}
my $stdev = sqrt($variance / $total_length);

my $min_cov = defined $forced_minimum ? $forced_minimum : $weighted_average - $stdev;
my $max_cov = $weighted_average + 2.5 * $stdev;  # factor here controls the upper tolerance

for my $id (@order){
	if(@order >= 4){
		next unless $cov_of{$id} <= $max_cov && $cov_of{$id} >= $min_cov;
	}
	print $out "$id\n$sequences{$id}\n";
}
close $out;
