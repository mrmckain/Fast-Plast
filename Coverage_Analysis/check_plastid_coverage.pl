#!/usr/bin/perl -w
use strict;

# Identify stretches of the assembly that the reads do not support.
#
# Usage: check_plastid_coverage.pl <coverage_file> <kmer_size> [min_coverage]
#
# The coverage file (from new_window_coverage.pl) has one line per assembly
# position: k-mer, position, count of that k-mer in the reads. A problem
# region is a run of at least <kmer_size> consecutive windows at or below the
# minimum coverage.
#
# Minimum coverage: the third argument if given; otherwise 0.15 x the mean
# window coverage, capped at MAX_MIN_COVERAGE. Without the cap a deep library
# (mean 5,900x) treated everything under 886x as unsupported, which is absurd:
# a window seen a few hundred times is supported no matter what the mean is.
# A real assembly error leaves its k-mers at roughly the sequencing error
# rate, i.e. under 1% of the mean, so the cap costs no sensitivity.
#
# Low-complexity dropouts: Illumina coverage collapses in extremely AT-rich
# stretches and long homopolymers, to single-digit k-mer counts even at
# thousands-fold mean depth. Those regions are real sequence that SPAdes
# assembled from the few reads that exist; flagging them as problems sends
# the run into needless reassembly. A candidate region whose sequence at the
# point of lowest coverage (the 25 k-mers centred on the minimum) is at least
# LOW_COMPLEXITY_AT AT- or GC-skewed, or contains a homopolymer of at least
# LOW_COMPLEXITY_HOMOPOLYMER bases, is therefore written to
# <id>_low_complexity_dropouts.txt (for the record) instead of the problem
# file. The first line of STDOUT is the minimum coverage used (the driver
# reads it); a second line summarises the counts.

my $MAX_MIN_COVERAGE          = 100;
my $LOW_COMPLEXITY_AT         = 0.85;   # fraction of A/T (or of G/C) over the region
my $LOW_COMPLEXITY_HOMOPOLYMER = 10;    # bases

my ($covfile, $kmer, $mincov) = @ARGV;
defined $kmer or die "Usage: $0 <coverage_file> <kmer_size> [min_coverage]\n";

my $exp_id;
$covfile =~ /(.*?)\.coverage_25kmer\.txt/;
$exp_id = $1;

# ---- pass 1: mean window coverage
open my $file, "<", $covfile or die "ERROR: cannot open $covfile: $!\n";
my $total_cov = 0;
my $total_windows = 0;
while(<$file>){
	chomp;
	my @tarray = split /\s+/;
	next unless defined $tarray[2] && $tarray[2] =~ /^\d/;
	$total_cov += $tarray[2];
	$total_windows++;
}
close $file;
$total_windows or die "ERROR: $covfile is empty.\n";

my $avg_cov = $total_cov / $total_windows;
unless(defined $mincov){
	$mincov = $avg_cov * 0.15;
	$mincov = $MAX_MIN_COVERAGE if $mincov > $MAX_MIN_COVERAGE;
}
print "$mincov\n";

# ---- pass 2: runs of low windows, classified by sequence complexity
my $current_kmer = 0;
my $current_kmer_start;   # positions start at 0, so test with defined(), not truthiness
my $current_kmer_end;
my $run_cov = 0;
my @run_kmers;            # k-mers of the current run, for the complexity test
my @run_covs;             # their coverages, to find where the run bottoms out

open my $out, ">", $exp_id . "_problem_regions_plastid_assembly.txt" or die "ERROR: cannot write problem-region file: $!\n";
open my $lc,  ">", $exp_id . "_low_complexity_dropouts.txt"       or die "ERROR: cannot write dropout file: $!\n";
my ($n_problem, $n_lowc) = (0, 0);

# Sequence of a run from its k-mers: first k-mer, then the last base of each.
sub run_sequence {
	my @k = @_;
	return "" unless @k;
	my $s = $k[0];
	$s .= substr($_, -1) for @k[1..$#k];
	return $s;
}
sub max_homopolymer {
	my ($s) = @_;
	my $max = 0;
	while($s =~ /((.)\2*)/g){ $max = length($1) if length($1) > $max; }
	return $max;
}

# Emit the current low-coverage run if it spans at least $kmer windows.
sub flush_run {
	if(defined $current_kmer_start && $current_kmer >= $kmer){
		my $span = $current_kmer_end - $current_kmer_start;
		my $run_avg = $span > 0 ? $run_cov / $span : $run_cov;
		# judge the sequence where coverage bottoms out (the 25 k-mers centred on
		# the minimum), not the whole run, whose flanks are ordinary sequence
		my $imin = 0;
		for my $i (1..$#run_covs){ $imin = $i if $run_covs[$i] < $run_covs[$imin]; }
		my $lo = $imin - 12; $lo = 0 if $lo < 0;
		my $hi = $imin + 12; $hi = $#run_kmers if $hi > $#run_kmers;
		my $seq = run_sequence(@run_kmers[$lo..$hi]);
		my $len = length($seq) || 1;
		my $at  = ($seq =~ tr/ATat//) / $len;
		my $gc  = ($seq =~ tr/GCgc//) / $len;
		my $hp  = max_homopolymer($seq);
		if($at >= $LOW_COMPLEXITY_AT || $gc >= $LOW_COMPLEXITY_AT || $hp >= $LOW_COMPLEXITY_HOMOPOLYMER){
			printf $lc "%d\t%d\t%.2f\tAT %.0f%%\tGC %.0f%%\tlongest homopolymer %d\n", $current_kmer_start, $current_kmer_end, $run_avg, 100*$at, 100*$gc, $hp;
			$n_lowc++;
		}
		else{
			print $out "$current_kmer_start\t$current_kmer_end\t$run_avg\n";
			$n_problem++;
		}
	}
	$current_kmer = 0;
	$current_kmer_start = undef;
	$current_kmer_end = undef;
	$run_cov = 0;
	@run_kmers = ();
	@run_covs = ();
}

open $file, "<", $covfile or die "ERROR: cannot reopen $covfile: $!\n";
while(<$file>){
	chomp;
	my @tarray = split /\s+/;
	next unless defined $tarray[2] && $tarray[2] =~ /^\d/;

	if($tarray[2] <= $mincov && $tarray[0] !~ /N/i){
		$current_kmer_start = $tarray[1] unless defined $current_kmer_start;
		$current_kmer_end = $tarray[1];
		$current_kmer++;
		$run_cov += $tarray[2];
		push @run_kmers, $tarray[0];
		push @run_covs, $tarray[2];
	}
	else{
		flush_run();
	}
}
close $file;
# a run that reaches the last window must be reported too
flush_run();
close $out;
close $lc;

printf STDOUT "problem regions: %d; low-complexity dropouts (not counted): %d; mean window coverage %.0f\n", $n_problem, $n_lowc, $avg_cov;
