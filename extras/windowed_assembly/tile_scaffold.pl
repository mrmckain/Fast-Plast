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
# interval. Where the two tiles do not share an exact anchor it accepts a
# reference overlap of 10 bp or more that is identical in both tiles, then an
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
    # The join is made at the end of the previous tile and the start of this
    # one, so use placements that lie there (within $END_SLACK bp of the tile's
    # end / start): of those, the previous tile's placement that reaches
    # furthest along the reference (a short second placement of a repeat
    # inside the window must not be taken for the tile's end) and this tile's
    # that starts earliest. A tile with an IR flank at its far end has a
    # placement in the window there too; joining on it would throw the body
    # of the tile away. Without a placement near the end, take the nearest.
    my $END_SLACK = 3000;
    my @A = pick($ckey, $cwin); my $lenA = length($ori{$ckey});
    my @Aend = grep { $lenA - $_->[1] <= $END_SLACK } @A;
    my $Alast = @Aend ? (sort { $a->[3] <=> $b->[3] } @Aend)[-1] : (sort { $a->[1] <=> $b->[1] } @A)[-1];
    my @Bh = pick($key, $lw);
    my @Bstart = grep { $_->[0] <= $END_SLACK } @Bh;
    my $Bfirst = @Bstart ? $Bstart[0] : (sort { $a->[0] <=> $b->[0] } @Bh)[0];
    my ($eA, $sB) = ($Alast->[3], $Bfirst->[2]);
    my $done = 0;
    if ($sB <= $eA) {                                           # reference overlap: exact 40-mer anchor
        # Try anchors from the middle of the shared interval to its end, then
        # from its start to the middle: when the previous tile runs past the
        # true overlap on a diverged copy of a repeat, only the first part of
        # the interval is shared sequence.
        my $mid = int(($sB + $eA) / 2);
        my @try;
        for (my $P = $mid; $P <= $eA - 40; $P += 50) { push @try, $P }
        for (my $P = $sB;  $P <  $mid && $P <= $eA - 40; $P += 50) { push @try, $P }
        # Every usable anchor in the shared interval, in the order tried.
        my @anc;                                                # [P, position in tile A, position in tile B]
        for my $P (@try) {
            next if $P < $sB;
            my $qa = $Alast->[0] + ($P - $Alast->[2]);
            next if $qa < 0 || $qa + 40 > length($ori{$ckey});
            my $anchor = substr($ori{$ckey}, $qa, 40);
            next if length($anchor) < 40 || $anchor =~ /N/i;       # an empty anchor would match everywhere, for ever
            my $qb_est = $Bfirst->[0] + ($P - $Bfirst->[2]);
            my ($best, $bd); my $pos = -1;
            while (($pos = index($B, $anchor, $pos + 1)) >= 0) { my $d = abs($pos - $qb_est); ($best, $bd) = ($pos, $d) if !defined $bd || $d < $bd }
            next unless defined $best && $bd < 2000;
            push @anc, [$P, $qa, $best];
        }
        if (@anc) {
            # The first anchor tried is the default splice point: tile A before
            # it, tile B after it. Where the two tiles differ in length between
            # neighbouring anchors (an indel in the overlap), keep the longer
            # version whichever tile holds it: assemblers collapse sequence at
            # hairpins and short direct repeats far more often than they invent
            # it (the 65 bp IR deletion in Johnsongrass is such a collapse).
            my $default = $anc[0];
            my @chain;                                          # anchors in order, colinear in both tiles
            for my $c (sort { $a->[1] <=> $b->[1] } @anc) {
                push @chain, $c if !@chain || $c->[2] > $chain[-1][2];
            }
            @chain = ($default) unless grep { $_ == $default } @chain;
            my $merged = substr($scaf, 0, $off + $chain[0][1]);
            my @kept;
            for my $j (0 .. $#chain - 1) {
                my ($c, $d) = @chain[$j, $j + 1];
                my $la = $d->[1] - $c->[1]; my $lb = $d->[2] - $c->[2];
                my $useA = $la > $lb ? 1 : $lb > $la ? 0 : ($d->[1] <= $default->[1] ? 1 : 0);
                $merged .= $useA ? substr($ori{$ckey}, $c->[1], $la) : substr($B, $c->[2], $lb);
                push @kept, sprintf("kept the longer version (%d bp against %d bp, from %s) at reference ~%d-%d", ($la > $lb ? ($la, $lb, $ckey) : ($lb, $la, $key)), $c->[0], $d->[0]) if $la != $lb;
            }
            my $last = $chain[-1];
            $off = length($merged) - $last->[2];
            $scaf = $merged . substr($B, $last->[2]);
            push @log, sprintf("%s -> %s: spliced at reference ~%d (exact 40-mer anchor)", $ckey, $key, $default->[0]);
            push @log, "    $_" for @kept;
            $done = 1;
        }
    }
    if (!$done && $sB <= $eA && $eA - $sB + 1 >= 10) {          # short reference overlap, identical in both tiles
        my $ov = $eA - $sB + 1;
        my $a_end   = $Alast->[1];                              # tile A bases aligned up to here (1-based)
        my $b_start = $Bfirst->[0] - 1;                         # tile B aligned from here (0-based)
        if ($a_end >= $ov && substr($ori{$ckey}, $a_end - $ov, $ov) eq substr($B, $b_start, $ov)) {
            $scaf = substr($scaf, 0, $off + $a_end) . substr($B, $b_start + $ov);
            push @log, sprintf("%s -> %s: joined at an identical %d bp reference overlap", $ckey, $key, $ov);
            $off = $off + $a_end - ($b_start + $ov); $done = 1;
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
            push @log, sprintf("%s -> %s: no anchor; %d bp N gap (trimmed %d + %d unaligned tail bases); reference %d to %d", $ckey, $key, $gap, $trimA, $head, $eA, $sB);
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
