#!/usr/bin/perl -w
use strict;

####Script to take coverage and identify stretches that need to be fixed###

####Usage: 1-Coverage file 2-Kmer Size 3-Coverage Min (optional; default 0.15 * mean window coverage)

my ($covfile, $kmer, $mincov) = @ARGV;
defined $kmer or die "Usage: $0 <coverage_file> <kmer_size> [min_coverage]\n";

my $exp_id;
$covfile =~ /(.*?)\.coverage_25kmer\.txt/;
$exp_id = $1;

open my $file, "<", $covfile or die "ERROR: cannot open $covfile: $!\n";
my $total_cov = 0;
my $total_windows = 0;
while(<$file>){
	chomp;
	my @tarray = split /\s+/;
	$total_cov += $tarray[2];
	$total_windows++;
}
close $file;
$total_windows or die "ERROR: $covfile is empty.\n";

my $avg_cov = $total_cov / $total_windows;
$mincov = $avg_cov * 0.15 unless defined $mincov;
print "$mincov\n";

my $current_kmer = 0;
my $current_kmer_start;   # positions start at 0, so test with defined(), not truthiness
my $current_kmer_end;
my $run_cov = 0;

open my $out, ">", $exp_id . "_problem_regions_plastid_assembly.txt" or die "ERROR: cannot write problem-region file: $!\n";

# Emit the current low-coverage run if it spans at least $kmer windows.
sub flush_run {
	if(defined $current_kmer_start && $current_kmer >= $kmer){
		my $span = $current_kmer_end - $current_kmer_start;
		my $run_avg = $span > 0 ? $run_cov / $span : $run_cov;
		print $out "$current_kmer_start\t$current_kmer_end\t$run_avg\n";
	}
	$current_kmer = 0;
	$current_kmer_start = undef;
	$current_kmer_end = undef;
	$run_cov = 0;
}

open $file, "<", $covfile or die "ERROR: cannot reopen $covfile: $!\n";
while(<$file>){
	chomp;
	my @tarray = split /\s+/;

	if($tarray[2] <= $mincov && $tarray[0] !~ /N/i){
		$current_kmer_start = $tarray[1] unless defined $current_kmer_start;
		$current_kmer_end = $tarray[1];
		$current_kmer++;
		$run_cov += $tarray[2];
	}
	else{
		flush_run();
	}
}
close $file;
# A low-coverage run that reaches the last window was previously never
# reported, because it was only flushed when a good window followed it.
flush_run();
close $out;
