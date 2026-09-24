#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# tile_scaffold.pl -- reference-ordered tiling of contigs into one linear
# scaffold (experimental; a manual alternative to afin fusion when afin has
# joined contigs through a repeat).
#
# Given a FASTA of candidate contigs, a BLAST database of a reference plastome
# and a plan listing the tiles in reference order, the script orients each
# tile, places it on the reference with blastn, and splices each tile onto
# the previous one at an exact 40-mer anchor inside their shared reference
# interval. Where the two tiles do not share an exact anchor it looks for an
# exact end overlap of the raw sequences, and failing that bridges them with N
# (as many as the reference implies, at least 100), trimming any unaligned
# tail bases on both sides.
#
# Repeating the first tile at the end (as the Pogonachne example does with its
# wrap-around window) closes the circle and leaves an afin-style overhang so
# Fast-Plast's IR identification and orientation scripts find both IR copies.
#
# Usage:
#   tile_scaffold.pl --tiles tiles.fa --db ref_blastdb --plan plan.tsv --out scaffold.fa [--name ID]
#
# Plan file: one tile per line, tab-separated:
#   prefix   orientation   left_window   right_window
#   prefix       the start of the tile's FASTA id, unique among the tiles
#   orientation  + or - (reverse-complement before use)
#   left_window  reference interval (lo-hi) of the placement used to splice
#                this tile onto the previous one, or - for the first tile
#   right_window reference interval used to splice the next tile on, or -
#                for the last tile
#
# Example plan (Pogonachne racemosa on Sorghum bicolor NC_008602; NODE_* are
# error-corrected SPAdes contigs, w* are per-window SPAdes contigs from
# bin_sam_windows.pl that bridge junctions SPAdes broke):
#   w028_   -   -               1-5207
#   NODE_4_ -   404-13551       404-13551
#   NODE_1_ +   13538-59538     13538-59538
#   w011_   +   54834-65159     54834-65159
#   NODE_2_ -   59856-85451     59856-85451
#   w016_   +   79874-90157     79874-90157
#   NODE_3_ +   85447-105900    85447-105900
#   w020_   -   99841-110157    99841-110157
#   NODE_5_ -   105995-118495   105995-118495
#   w023_   +   114841-125183   114841-125183
#   NODE_3_ -   118588-139041   118588-139041
#   w028_   -   135594-140754   -
#
# Requires blastn on PATH.

my ($tiles_fa, $blastdb, $plan_file, $out, $name);
$name = "tiled_scaffold";
GetOptions('tiles=s' => \$tiles_fa, 'db=s' => \$blastdb, 'plan=s' => \$plan_file, 'out=s' => \$out, 'name=s' => \$name)
    or die "bad options\n";
$tiles_fa && $blastdb && $plan_file && $out
    or die "Usage: $0 --tiles tiles.fa --db ref_blastdb --plan plan.tsv --out scaffold.fa [--name ID]\n";

sub rc { my $s = reverse $_[0]; $s =~ tr/ACGTNacgtn/TGCANtgcan/; $s }

# ---- tiles
my (%seq, $id);
open my $f, "<", $tiles_fa or die "cannot open $tiles_fa: $!\n";
while (<$f>) { chomp; if (/^>(\S+)/) { $id = $1 } else { $seq{$id} .= $_ } }
close $f;
sub full {
    my $k = shift;
    my @m = grep { index($_, $k) == 0 } keys %seq;
    die "no tile with prefix $k\n" unless @m;
    die "prefix $k matches several tiles: @m\n" if @m > 1;
    $m[0];
}

# ---- plan
my @tiles;
open my $p, "<", $plan_file or die "cannot open $plan_file: $!\n";
while (<$p>) {
    chomp; next if /^\s*(#|$)/;
    my ($k, $o, $lw, $rw) = split /\s+/;
    die "bad plan line: $_\n" unless defined $rw && $o =~ /^[+-]$/;
    my $win = sub { my $w = shift; return undef if $w eq "-"; $w =~ /^(\d+)-(\d+)$/ or die "bad window '$w'\n"; [$1, $2] };
    push @tiles, [$k, $o, $win->($lw), $win->($rw)];
}
close $p;
@tiles >= 2 or die "the plan needs at least two tiles\n";

# ---- orient tiles and place each oriented tile once on the reference
my (%ori, %hsp);
for my $t (@tiles) {
    my ($k, $o) = @$t; my $key = "$k$o"; next if $ori{$key};
    my $s = $seq{ full($k) }; $s = rc($s) if $o eq "-"; $ori{$key} = $s;
    open my $q, ">", "tmp_tile.fa" or die; print $q ">$key\n$s\n"; close $q;
    my @h = `blastn -query tmp_tile.fa -db $blastdb -outfmt "6 qstart qend sstart send length" -evalue 1e-50 2>/dev/null`;
    for (@h) { chomp; my @c = split /\t/; next if $c[4] < 500 || $c[2] > $c[3]; push @{$hsp{$key}}, [@c] }   # + strand only
    die "no plus-strand placement for $key (wrong orientation?)\n" unless $hsp{$key};
}
unlink "tmp_tile.fa";
sub pick {
    my ($key, $win) = @_;
    my @c = grep { $_->[2] >= $win->[0] - 500 && $_->[3] <= $win->[1] + 500 } @{$hsp{$key}};
    @c or die "no placement of $key inside reference window $win->[0]-$win->[1]\n";
    sort { $a->[2] <=> $b->[2] } @c;
}

# ---- tile
my $scaf = ""; my $off = 0; my ($ckey, $cwin); my @log;
for my $i (0..$#tiles) {
    my ($k, $o, $lw, $rw) = @{$tiles[$i]}; my $key = "$k$o"; my $B = $ori{$key};
    if ($i == 0) { $scaf = $B; $off = 0; ($ckey, $cwin) = ($key, $rw); push @log, "start with $key (" . length($B) . " bp)"; next }
    die "tile $key has no left window\n" unless $lw;
    die "tile $ckey has no right window\n" unless $cwin;
    my @A = pick($ckey, $cwin); my $Alast = $A[-1];
    my @Bh = pick($key, $lw);    my $Bfirst = $Bh[0];
    my ($eA, $sB) = ($Alast->[3], $Bfirst->[2]);
    my $done = 0;
    if ($sB <= $eA) {                                           # reference overlap: exact 40-mer anchor
        for (my $P = int(($sB + $eA) / 2); $P <= $eA - 40 && $P >= $sB; $P += 50) {
            my $qa = $Alast->[0] + ($P - $Alast->[2]);
            my $anchor = substr($ori{$ckey}, $qa, 40);
            my $qb_est = $Bfirst->[0] + ($P - $Bfirst->[2]);
            my ($best, $bd); my $pos = -1;
            while (($pos = index($B, $anchor, $pos + 1)) >= 0) { my $d = abs($pos - $qb_est); ($best, $bd) = ($pos, $d) if !defined $bd || $d < $bd }
            next unless defined $best && $bd < 2000;
            $scaf = substr($scaf, 0, $off + $qa) . substr($B, $best);
            push @log, sprintf("%s -> %s: spliced at reference ~%d (exact 40-mer anchor)", $ckey, $key, $P);
            $off = $off + $qa - $best; $done = 1; last;
        }
    }
    unless ($done) {                                            # exact end overlap of the raw sequences
        my $found = 0;
        for (my $L = 5000; $L >= 30; $L--) {
            if (length($scaf) >= $L && length($B) >= $L && substr($scaf, -$L) eq substr($B, 0, $L)) {
                $scaf .= substr($B, $L); $off = length($scaf) - length($B);
                push @log, "$ckey -> $key: joined by exact $L bp end overlap"; $found = 1; last;
            }
        }
        unless ($found) {                                       # N gap
            my $trimA = length($ori{$ckey}) - $Alast->[1];
            $scaf = substr($scaf, 0, length($scaf) - $trimA) if $trimA > 0;
            my $gap = $sB - $eA - 1; $gap = 100 if $gap < 100;
            my $head = $Bfirst->[0] - 1;
            $scaf .= ("N" x $gap) . substr($B, $head);
            $off = length($scaf) - (length($B) - $head);
            push @log, sprintf("%s -> %s: no anchor; %d bp N gap (trimmed %d + %d unaligned tail bases)", $ckey, $key, $gap, $trimA, $head);
            $_->[0] -= $head, $_->[1] -= $head for @{$hsp{$key}}; $ori{$key} = substr($B, $head);
        }
    }
    ($ckey, $cwin) = ($key, $rw);
}

open my $o, ">", $out or die "cannot write $out: $!\n";
print $o ">$name\n$scaf\n"; close $o;
print "$_\n" for @log;
my $n = () = $scaf =~ /N+/g;
printf "scaffold: %d bp, %d N gap(s), %d N total -> %s\n", length $scaf, $n, ($scaf =~ tr/N//), $out;
