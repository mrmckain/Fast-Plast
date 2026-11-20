#!/usr/bin/env perl
use strict;
use warnings;

# coverage_ir_expand.pl -- rebuild a plastome from one contig that carries a
# single IR copy.
#
# Usage: coverage_ir_expand.pl <contig.fa> <out.fa> <report.txt> <reads.fq> [reads.fq ...]
# Exit 0 = rebuilt (out.fa written), 2 = not applicable / not supported (no output).
#
# afin can return a whole plastome as one contig LSC-IR-SSC: the two IR copies
# collapse into one in the assembly graph, so the walk passes through the IR
# once. Such a contig has no inverted repeat, and sequence_based_ir_id.pl
# cannot find an IR in it. The reads still show the IR: its k-mers are covered
# by reads from both copies, about twice the single-copy depth. This script
#   1. finds the IR as the longest stretch at >= $IR_RATIO x the median k-mer
#      coverage, bridging dips of up to $MAX_DIP windows (length $MIN_IR..$MAX_IR);
#   2. stops (exit 2) if the contig already holds a second, inverted copy of it;
#   3. trims overhangs: a contig end that continued into the other IR copy is
#      the reverse complement of part of the IR (>= $MIN_OVH bp, found inside
#      the IR stretch);
#   4. places the IR ends exactly where reads join a contig end to the reverse
#      complement of the IR (the two junctions the rebuilt circle adds: end of
#      the second single-copy region -> second IR copy -> start of the first);
#      each needs >= $MIN_JUNCTION supporting reads and a clear best position;
#   5. writes SC1 + IR + SC2 + revcomp(IR).
# The output is one linear sequence holding both IR copies, for the usual
# orientation step.

my $K            = 31;
my $WIN          = 500;
my $IR_RATIO     = 1.6;
my $MIN_IR       = 8000;
my $MAX_IR       = 50000;
my $MIN_OVH      = 25;
my $SEARCH       = 2500;   # +/- bases around the coverage step searched for the exact IR end
my $MIN_JUNCTION = 3;

my ($in, $out, $report, @reads) = @ARGV;
@reads or die "Usage: $0 <contig.fa> <out.fa> <report.txt> <reads.fq> [reads.fq ...]\n";

open my $rep, ">", $report or die "ERROR: cannot write $report: $!\n";
sub stop { print $rep "NOT_APPLIED\t$_[0]\n"; close $rep; print STDERR "coverage_ir_expand.pl: not applied: $_[0]\n"; exit 2; }

my (@ids, %seq);
{
	open my $f, "<", $in or die "ERROR: cannot read $in: $!\n";
	my $id;
	while(<$f>){
		chomp;
		if(/^>(\S+)/){ $id = $1; push @ids, $id; $seq{$id} = ""; }
		elsif(defined $id){ $seq{$id} .= uc $_; }
	}
	close $f;
}
stop("expects one contig, found " . scalar(@ids)) unless @ids == 1;
my $C = $seq{$ids[0]};
stop("contig contains N") if $C =~ /N/;
my $L = length $C;
stop("contig too short ($L bp)") if $L < 3 * $MIN_IR;

sub rc { my $s = reverse $_[0]; $s =~ tr/ACGTN/TGCAN/; return $s; }

# all read k-mers, both strands
my %cnt;
for my $fq (@reads){
	open my $f, "<", $fq or die "ERROR: cannot read $fq: $!\n";
	while(<$f>){
		my $s = <$f>; <$f>; <$f>;
		last unless defined $s;
		chomp $s; $s = uc $s;
		for my $r ($s, rc($s)){
			for(my $i = 0; $i + $K <= length $r; $i++){
				$cnt{substr($r, $i, $K)}++;
			}
		}
	}
	close $f;
}
sub sup { return $cnt{$_[0]} || 0; }

# ---- 1. coverage step -------------------------------------------------------
my @cov = map { sup(substr($C, $_, $K)) } 0 .. $L - $K;
sub median { my @s = sort { $a <=> $b } @_; return @s ? $s[int(@s / 2)] : 0; }
my @w;
for(my $i = 0; $i < @cov; $i += $WIN){
	my $j = $i + $WIN - 1; $j = $#cov if $j > $#cov;
	push @w, median(@cov[$i .. $j]);
}
my $M = median(@w);
stop("no read coverage on the contig") unless $M > 0;
my @hi = map { $_ >= $IR_RATIO * $M ? 1 : 0 } @w;
# runs of high windows, merged across dips of up to $MAX_DIP windows (AT-rich or
# otherwise poorly covered stretches inside the IR); the longest merged run is the IR
my $MAX_DIP = 4;
my @runs;
for my $i (0 .. $#hi){
	next unless $hi[$i];
	if(@runs && $i - $runs[-1][1] - 1 <= $MAX_DIP){ $runs[-1][1] = $i; }
	else{ push @runs, [$i, $i]; }
}
my ($bs, $be) = (-1, -1);
for my $r (@runs){ ($bs, $be) = @$r if $bs < 0 || $r->[1] - $r->[0] > $be - $bs; }
stop(sprintf("no stretch at >= %.1fx the median k-mer coverage (%d)", $IR_RATIO, $M)) if $bs < 0;
my ($s0, $e0) = ($bs * $WIN, ($be + 1) * $WIN);
$e0 = $L if $e0 > $L;
my $irlen0 = $e0 - $s0;
stop("high-coverage stretch is $irlen0 bp, not IR-sized ($MIN_IR-$MAX_IR)") if $irlen0 < $MIN_IR || $irlen0 > $MAX_IR;
my $ratio = median(@w[$bs .. $be]) / $M;

# ---- 2. already two copies? ------------------------------------------------
{
	my $mid = substr($C, $s0 + int($irlen0 / 2) - 500, 1000);
	my $rcm = rc($mid);
	stop("the contig already holds an inverted copy of the IR") if index($C, $rcm) >= 0;
}

# ---- 3. overhangs into the other IR copy ----------------------------------
my $irzone = substr($C, ($s0 > $SEARCH ? $s0 - $SEARCH : 0), $irlen0 + 2 * $SEARCH);
my ($ov1, $ov2) = (0, 0);
# start: longest prefix whose reverse complement lies in the IR stretch
if(index($irzone, rc(substr($C, 0, $MIN_OVH))) >= 0){
	my ($lo, $hi) = ($MIN_OVH, ($e0 - $s0 < $s0 ? $e0 - $s0 : $s0));
	while($lo < $hi){ my $m = int(($lo + $hi + 1) / 2); if(index($irzone, rc(substr($C, 0, $m))) >= 0){ $lo = $m; } else{ $hi = $m - 1; } }
	$ov1 = $lo;
}
if(index($irzone, rc(substr($C, $L - $MIN_OVH))) >= 0){
	my ($lo, $hi) = ($MIN_OVH, ($e0 - $s0 < $L - $e0 ? $e0 - $s0 : $L - $e0));
	while($lo < $hi){ my $m = int(($lo + $hi + 1) / 2); if(index($irzone, rc(substr($C, $L - $m))) >= 0){ $lo = $m; } else{ $hi = $m - 1; } }
	$ov2 = $lo;
}
my $Ct = substr($C, $ov1, $L - $ov1 - $ov2);
my ($s1, $e1) = ($s0 - $ov1, $e0 - $ov1);
my $Lt = length $Ct;

# ---- 4. exact IR ends from the junction reads ------------------------------
# The circle adds two junctions:
#   end of the contig -> revcomp(IR end):   Ct[.. t] + rc(Ct[e-16 .. e-1])
#   revcomp(IR start) -> start of contig:   rc(Ct[s .. s+15]) + Ct[h ..]
# afin stops a contig where reads thin out or turn noisy, so the last/first
# bases can be wrong or run a few bases into the other IR copy. Each end is
# therefore searched jointly over a trim of the contig tip (0..$MAX_TRIM) and
# the IR boundary position; the smallest trim whose junction support is at
# least half the best is taken. Support must reach $MIN_JUNCTION reads and
# $MIN_JUNCTION_FRAC of the median single-copy k-mer count, and the best
# boundary must be clear of the runner-up.
my $MAX_TRIM = 300;
my $MIN_JUNCTION_FRAC = 0.2;
sub best_junction {
	my ($kmer_at, $from, $to) = @_;   # $kmer_at->(trim, pos)
	my %best_by_trim;                 # trim -> [pos, support, second-best support]
	for my $t (0 .. $MAX_TRIM){
		my ($bp, $bs, $b2) = (undef, 0, 0);
		for my $p ($from .. $to){
			my $k = $kmer_at->($t, $p); next unless defined $k;
			my $c = sup($k);
			if($c > $bs){ ($b2, $bs, $bp) = ($bs, $c, $p); } elsif($c > $b2){ $b2 = $c; }
		}
		$best_by_trim{$t} = [$bp, $bs, $b2] if defined $bp;
	}
	my $max = 0; $max < $_->[1] and $max = $_->[1] for values %best_by_trim;
	return () unless $max;
	for my $t (sort { $a <=> $b } keys %best_by_trim){
		my ($bp, $bs, $b2) = @{$best_by_trim{$t}};
		return ($t, $bp, $bs, $b2) if $bs >= $max / 2;
	}
	return ();
}
my $need = $MIN_JUNCTION; $need = int($MIN_JUNCTION_FRAC * $M) if $MIN_JUNCTION_FRAC * $M > $need;
my $lo_e = $e1 - $SEARCH; $lo_e = 16 if $lo_e < 16;
my $hi_e = $e1 + $SEARCH; $hi_e = $Lt - $MAX_TRIM - 15 if $hi_e > $Lt - $MAX_TRIM - 15;
my ($te, $e, $se, $se2) = best_junction(sub { my ($t, $p) = @_; return substr($Ct, $Lt - $t - 15, 15) . rc(substr($Ct, $p - 16, 16)); }, $lo_e, $hi_e);
my $lo_s = $s1 - $SEARCH; $lo_s = $MAX_TRIM + 15 if $lo_s < $MAX_TRIM + 15;
my $hi_s = $s1 + $SEARCH; $hi_s = $Lt - 16 if $hi_s > $Lt - 16;
my ($ts, $s, $ss, $ss2) = best_junction(sub { my ($t, $p) = @_; return rc(substr($Ct, $p, 16)) . substr($Ct, $t, 15); }, $lo_s, $hi_s);
stop("no read support for the junction second single-copy region -> second IR copy") unless defined $e;
stop("no read support for the junction second IR copy -> first single-copy region") unless defined $s;
stop("junction second single-copy region -> second IR copy has $se reads (need $need)") if $se < $need;
stop("junction second IR copy -> first single-copy region has $ss reads (need $need)") if $ss < $need;
stop("IR end ambiguous (junction support $se vs $se2)") if $se2 >= $se / 2;
stop("IR start ambiguous (junction support $ss vs $ss2)") if $ss2 >= $ss / 2;
stop("IR ends out of order (start $s, end $e)") unless $s < $e;
# drop the trimmed tips (they are wrong or belong to the other IR copy)
$Ct = substr($Ct, $ts, $Lt - $ts - $te);
($s, $e) = ($s - $ts, $e - $ts);
$Lt = length $Ct;
my $IR = substr($Ct, $s, $e - $s);
my $irlen = length $IR;
stop("IR from the junctions is $irlen bp, not IR-sized") if $irlen < $MIN_IR || $irlen > $MAX_IR;

# ---- 5. rebuild ------------------------------------------------------------
my $new = $Ct . rc($IR);
open my $o, ">", $out or die "ERROR: cannot write $out: $!\n";
print $o ">$ids[0]_IR_rebuilt\n";
for(my $i = 0; $i < length $new; $i += 80){ print $o substr($new, $i, 80), "\n"; }
close $o;
printf $rep "APPLIED\tcontig %d bp; IR %d-%d of the trimmed contig (%d bp, %.2fx median k-mer coverage %d); overhangs removed %d + %d bp, tips trimmed %d (start) + %d (end) bp; junction reads %d (end) %d (start); SC1 %d bp, SC2 %d bp; rebuilt %d bp\n",
	$L, $s + 1, $e, $irlen, $ratio, $M, $ov1, $ov2, $ts, $te, $se, $ss, $s, $Lt - $e, length $new;
close $rep;
printf STDERR "coverage_ir_expand.pl: IR %d bp (%.2fx coverage), rebuilt %d bp\n", $irlen, $ratio, length $new;
exit 0;
