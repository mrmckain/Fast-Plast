#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);

# bin_sam_windows.pl -- partition mapped reads into overlapping reference
# windows for per-window local assembly (experimental).
#
# Reads a SAM produced by Fast-Plast's bowtie2 step (2_Bowtie_Mapping/*.sam),
# picks one reference sequence (the one with the most mapped reads unless --ref
# is given), tiles it with overlapping windows, and writes a FASTQ set per
# window. A pair is assigned to every window its fragment overlaps, so the two
# mates always travel together; an unmapped mate is carried along with its
# mapped partner. Reads are written in their original orientation (the SAM
# stores reverse-strand reads reverse-complemented).
#
# A final "wrap" window joins the two ends of the reference so a circular
# plastome's arbitrary origin does not become an assembly break.
#
# Usage:
#   bin_sam_windows.pl --sam mapping.sam --out windows/ [--window 10000] [--step 5000]
#                      [--ref NAME] [--max-tlen 2000] [--min-mapq 0]
#
# Output:
#   <out>/reference_counts.tsv        mapped reads per reference sequence
#   <out>/windows.tsv                 one row per window: reads, pairs, est. depth
#   <out>/w<NNN>_<start>-<end>/       R1.fq R2.fq (pairs) and U.fq (singles)
#
# No samtools or BAM support: it streams the SAM once (plus one pass to count
# per-reference reads), which is fast enough for the read counts Fast-Plast
# produces after mapping.

my ($sam, $out, $ref_pick);
my $window   = 10000;
my $step     = 5000;
my $max_tlen = 2000;
my $min_mapq = 0;
GetOptions(
    'sam=s'      => \$sam,
    'out=s'      => \$out,
    'window=i'   => \$window,
    'step=i'     => \$step,
    'ref=s'      => \$ref_pick,
    'max-tlen=i' => \$max_tlen,
    'min-mapq=i' => \$min_mapq,
) or die "bad options\n";
$sam && $out or die "Usage: $0 --sam <file.sam> --out <dir> [--window N --step N --ref NAME]\n";
$step > 0 && $step <= $window or die "--step must be > 0 and <= --window\n";

# ---------------------------------------------------------------- pass 1
# Reference lengths from @SQ and mapped-read tally per reference.
my %reflen;
my %refcount;
open my $in, "<", $sam or die "cannot open $sam: $!\n";
while (<$in>) {
    if (/^\@/) {
        if (/^\@SQ\tSN:(\S+)\tLN:(\d+)/) { $reflen{$1} = $2; }
        next;
    }
    my ($qname, $flag, $rname) = split /\t/, $_, 4;
    next if $flag & 4;          # unmapped
    next if $flag & 0x900;      # secondary / supplementary
    $refcount{$rname}++;
}
close $in;
%refcount or die "no mapped reads in $sam\n";

make_path($out);
open my $rc, ">", "$out/reference_counts.tsv" or die $!;
print $rc "reference\tlength\tmapped_reads\n";
for my $r (sort { $refcount{$b} <=> $refcount{$a} } keys %refcount) {
    print $rc join("\t", $r, ($reflen{$r} // "NA"), $refcount{$r}), "\n";
}
close $rc;

my $ref = $ref_pick // (sort { $refcount{$b} <=> $refcount{$a} } keys %refcount)[0];
exists $reflen{$ref} or die "reference '$ref' not in SAM header\n";
my $L = $reflen{$ref};
printf STDERR "reference: %s (%d bp, %d mapped reads, %.1f%% of all mapped)\n",
    $ref, $L, $refcount{$ref} // 0,
    100 * ($refcount{$ref} // 0) / (eval { my $t = 0; $t += $_ for values %refcount; $t } || 1);

# ---------------------------------------------------------------- windows
# [start,end) half-open, 0-based. The last regular window is extended to L.
# The wrap window covers the final half-window and the first half-window.
my @win;   # [id, start, end, is_wrap]
for (my $s = 0; $s < $L; $s += $step) {
    my $e = $s + $window;
    if ($e >= $L) { push @win, [scalar(@win), $s, $L, 0]; last; }
    push @win, [scalar(@win), $s, $e, 0];
}
my $half = int($window / 2);
push @win, [scalar(@win), $L - $half, $half, 1] if $L > $window;

sub windows_overlapping {
    my ($a, $b) = @_;    # fragment span [a,b), 0-based
    my @hit;
    for my $w (@win) {
        my (undef, $s, $e, $wrap) = @$w;
        if ($wrap) {
            push @hit, $w if $b > $s || $a < $e;
        }
        else {
            push @hit, $w if $a < $e && $b > $s;
        }
    }
    return @hit;
}

# ---------------------------------------------------------------- helpers
sub aln_end {                    # 0-based exclusive end from POS + CIGAR
    my ($pos, $cigar) = @_;
    my $len = 0;
    while ($cigar =~ /(\d+)([MDN=X])/g) { $len += $1; }
    return $pos - 1 + $len;
}
sub revcomp { my $s = reverse $_[0]; $s =~ tr/ACGTNacgtn/TGCANtgcan/; return $s; }
sub fastq_record {               # original orientation
    my ($qname, $flag, $seq, $qual, $mate) = @_;
    if ($flag & 16) { $seq = revcomp($seq); $qual = reverse $qual; }
    my $suffix = $mate ? "/$mate" : "";
    return "\@$qname$suffix\n$seq\n+\n$qual\n";
}

# open output handles lazily per window
my %fh;
sub handle {
    my ($w, $kind) = @_;
    my $dir = sprintf("%s/w%03d_%d-%d%s", $out, $w->[0], $w->[1], $w->[2], $w->[3] ? "_wrap" : "");
    my $key = "$dir/$kind";
    unless ($fh{$key}) {
        make_path($dir);
        open $fh{$key}, ">", "$key.fq" or die "cannot write $key.fq: $!\n";
    }
    return $fh{$key};
}
my %stat;   # window id -> {pairs, singles, bases}

sub emit_pair {
    my ($a, $b, @spans) = @_;    # two SAM records (arrayrefs), one or more [s,e] spans
    # union of windows over all spans, so a pair is written once per window
    my %seen;
    my @targets = grep { !$seen{$_->[0]}++ } map { windows_overlapping(@$_) } @spans;
    for my $w (@targets) {
        my ($r1, $r2) = ($a->[1] & 64) ? ($a, $b) : ($b, $a);
        print { handle($w, "R1") } fastq_record(@{$r1}[0,1,9,10], 1);
        print { handle($w, "R2") } fastq_record(@{$r2}[0,1,9,10], 2);
        $stat{$w->[0]}{pairs}++;
        $stat{$w->[0]}{bases} += length($r1->[9]) + length($r2->[9]);
    }
}
sub emit_single {
    my ($a, $span_s, $span_e) = @_;
    for my $w (windows_overlapping($span_s, $span_e)) {
        print { handle($w, "U") } fastq_record(@{$a}[0,1,9,10], 0);
        $stat{$w->[0]}{singles}++;
        $stat{$w->[0]}{bases} += length($a->[9]);
    }
}

# ---------------------------------------------------------------- pass 2
# Pair records by QNAME. bowtie2 writes mates consecutively, but buffering by
# name is cheap for the mapped subset and does not depend on that.
my %pending;
my ($n_pairs, $n_singles, $n_other_ref) = (0, 0, 0);
open $in, "<", $sam or die $!;
while (<$in>) {
    next if /^\@/;
    chomp;
    my @f = split /\t/;
    my ($qname, $flag, $rname, $pos, $mapq) = @f[0..4];
    next if $flag & 0x900;
    my $mapped = !($flag & 4);
    my $on_ref = $mapped && $rname eq $ref;
    if ($mapped && !$on_ref) { $n_other_ref++; }

    unless ($flag & 1) {                              # unpaired read
        next unless $on_ref && $mapq >= $min_mapq;
        emit_single(\@f, $pos - 1, aln_end($pos, $f[5]));
        $n_singles++;
        next;
    }

    # paired: skip pairs where neither mate is mapped
    next if (($flag & 4) && ($flag & 8));

    if (my $m = delete $pending{$qname}) {
        my ($a, $b) = ($m, \@f);
        my @spans;
        for my $r ($a, $b) {
            next if $r->[1] & 4;
            next unless $r->[2] eq $ref && $r->[4] >= $min_mapq;
            push @spans, [$r->[3] - 1, aln_end($r->[3], $r->[5])];
        }
        next unless @spans;                           # neither mate on chosen ref
        my $tlen = abs($a->[8] || 0);
        my ($s, $e);
        if (@spans == 2 && $tlen > 0 && $tlen <= $max_tlen) {
            $s = ($spans[0][0] < $spans[1][0]) ? $spans[0][0] : $spans[1][0];
            $e = ($spans[0][1] > $spans[1][1]) ? $spans[0][1] : $spans[1][1];
            emit_pair($a, $b, [$s, $e]);
            $n_pairs++;
        }
        else {
            # discordant / distant / one mate unmapped: place the pair by each
            # mapped mate separately (a mate may sit in another window or be
            # unplaced; keeping both reads together still helps the assembler).
            emit_pair($a, $b, @spans);
            $n_pairs++;
        }
    }
    else {
        $pending{$qname} = \@f;
    }
}
close $in;
close $_ for values %fh;

my $orphans = scalar keys %pending;
printf STDERR "pairs placed: %d   singles placed: %d   mapped to other references: %d   unpaired leftovers: %d\n",
    $n_pairs, $n_singles, $n_other_ref, $orphans;

open my $ws, ">", "$out/windows.tsv" or die $!;
print $ws "window\tstart\tend\twrap\tpairs\tsingles\tbases\test_depth\n";
for my $w (@win) {
    my $st  = $stat{$w->[0]} || {};
    my $len = $w->[3] ? ($L - $w->[1]) + $w->[2] : $w->[2] - $w->[1];
    printf $ws "w%03d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\n",
        $w->[0], $w->[1], $w->[2], $w->[3], $st->{pairs} // 0, $st->{singles} // 0,
        $st->{bases} // 0, ($st->{bases} // 0) / $len;
}
close $ws;
print STDERR "wrote ", scalar(@win), " windows under $out/\n";
