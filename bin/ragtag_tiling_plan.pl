#!/usr/bin/env perl
use strict;
use warnings;

# ragtag_tiling_plan.pl -- turn a RagTag scaffold into a plan for
# extras/windowed_assembly/tile_scaffold.pl.
#
# Usage: ragtag_tiling_plan.pl <ragtag.scaffold.agp> <ragtag.scaffold.asm.paf> <plan.tsv> <unplaced.txt>
#
# RagTag orders and orients the contigs on the reference and joins them with
# N. Neighbouring contigs often overlap on the reference (afin stopped at a
# repeat, or the two contigs differ inside the overlap), and tile_scaffold.pl
# can splice those at an exact shared 40-mer instead. This script writes the
# plan it needs: one line per placed contig, in scaffold order, with the
# contig's orientation and the reference interval RagTag placed it on (used
# as both the left and the right splice window). Ids of contigs RagTag left
# unplaced go to <unplaced.txt>.
#
# The orientation of a contig that holds a whole single-copy region between
# two IR ends is arbitrary (the plastome's two isomers), and RagTag's choice
# can leave the IR sequence that joins it to its neighbour on the far side of
# the contig. So for the first and the last contig the script also writes
# plans with that contig reverse-complemented: <plan>.flipfirst,
# <plan>.fliplast and <plan>.flipboth, and where the scaffold has a wide
# uncovered stretch a <plan>.rotated (see below). The caller tiles each plan
# and keeps the best.
#
# Exits non-zero without writing a plan when there is nothing to tile: fewer
# than two placed contigs, more than one scaffold, or a placed contig with no
# usable alignment in its scaffold orientation.

my $WINDOW_MARGIN = 1500;
my ($agp, $paf, $plan, $unplaced) = @ARGV;
defined $unplaced or die "Usage: $0 <agp> <paf> <plan.tsv> <unplaced.txt>\n";

my (@placed, @loose, %objects);
open my $a, "<", $agp or die "cannot open $agp: $!\n";
while(<$a>){
	next if /^#/;
	chomp;
	my @f = split /\t/;
	next unless @f >= 9 && $f[4] eq "W";
	if($f[0] =~ /_RagTag$/){
		$objects{$f[0]} = 1;
		push @placed, [$f[5], $f[8]];
	}
	else{
		push @loose, $f[5];
	}
}
close $a;
exit 2 unless @placed >= 2 && keys %objects == 1;

# Reference extent of each contig on each strand, from the alignments RagTag
# itself uses (at least 1000 bp, mapping quality at least 10).
my %extent;
open my $p, "<", $paf or die "cannot open $paf: $!\n";
while(<$p>){
	my @f = split /\t/;
	next unless @f >= 12 && $f[10] >= 1000 && $f[11] >= 10;
	my $e = $extent{$f[0]}{$f[4]} ||= [$f[7] + 1, $f[8]];
	$e->[0] = $f[7] + 1 if $f[7] + 1 < $e->[0];
	$e->[1] = $f[8]     if $f[8] > $e->[1];
}
close $p;

my @lines;
for my $c (@placed){
	my ($id, $strand) = @$c;
	my $e = $extent{$id}{$strand} or exit 3;
	# tile_scaffold.pl places tiles with blastn, which can carry an alignment
	# further than minimap2 did (through a diverged repeat copy at a contig
	# end), so give the window some room.
	my $lo = $e->[0] - $WINDOW_MARGIN; $lo = 1 if $lo < 1;
	my $hi = $e->[1] + $WINDOW_MARGIN;
	push @lines, join("\t", $id, $strand, "$lo-$hi", "$lo-$hi");
}
my @full = @lines;   # every tile with both windows, for the rotated plan below
# tile_scaffold.pl wants no left window on the first tile and no right window on the last.
{ my @f = split /\t/, $lines[0];  $f[2] = "-"; $lines[0]  = join("\t", @f); }
{ my @f = split /\t/, $lines[-1]; $f[3] = "-"; $lines[-1] = join("\t", @f); }

sub write_plan {
	my ($file, @l) = @_;
	open my $o, ">", $file or die "cannot write $file: $!\n";
	print $o "$_\n" for @l;
	close $o;
}
write_plan($plan, @lines);

# Alternative plans with an end contig reverse-complemented. minimap2 reports
# no second alignment for the short IR end that would join it, so the window
# is simply the stretch of reference around the neighbouring contig's near
# edge; tile_scaffold.pl places the flipped contig there with blastn, and the
# plan fails harmlessly if it has no placement in that window. The placement
# has to fit inside the window, and the flipped contig can carry a whole IR
# copy (up to ~28 kb in most plastomes), so reach that far either side.
my $FLIP_REACH = 30000;
sub flipped_line {
	my ($idx, $neighbour) = @_;
	my ($id, $strand) = @{$placed[$idx]};
	my $other = $strand eq "+" ? "-" : "+";
	my $ne = $extent{$placed[$neighbour][0]}{$placed[$neighbour][1]};
	my $target = $idx < $neighbour ? $ne->[0] : $ne->[1];   # the neighbour's near edge
	my $lo = $target - $FLIP_REACH; $lo = 1 if $lo < 1;
	my $hi = $target + $FLIP_REACH;
	return $idx < $neighbour ? join("\t", $id, $other, "-", "$lo-$hi") : join("\t", $id, $other, "$lo-$hi", "-");
}
my $ff = flipped_line(0, 1);
my $fl = flipped_line($#placed, $#placed - 1);
if(@placed == 2){
	# With two contigs each line carries only one window, so both flips are simple swaps.
	write_plan("$plan.flipfirst", $ff, $lines[1]) if $ff;
	write_plan("$plan.fliplast",  $lines[0], $fl) if $fl;
	write_plan("$plan.flipboth",  $ff, $fl) if $ff && $fl;
}
else{
	write_plan("$plan.flipfirst", $ff, @lines[1..$#lines]) if $ff;
	write_plan("$plan.fliplast",  @lines[0..$#lines-1], $fl) if $fl;
	write_plan("$plan.flipboth",  $ff, @lines[1..$#lines-1], $fl) if $ff && $fl;
}

# Rotated plan. A contig that runs across the reference's first base (the end
# of one IR copy plus the start of the LSC) is placed by its LSC part, at the
# start of the scaffold, and its IR end hangs off before the reference. The
# contig at the far end of the scaffold can then only join it around the
# circle, and the scaffold shows a stretch of N the size of the IR copy the
# assembly collapsed. When neighbouring contigs leave at least $ROTATE_GAP bp
# of reference uncovered, also write <plan>.rotated: the same tiles started
# after that stretch, so the old last tile is spliced onto the old first one.
my $ROTATE_GAP = 5000;
{
	my ($cut, $widest) = (undef, $ROTATE_GAP - 1);
	for my $i (0 .. $#placed - 1){
		my $e1 = $extent{$placed[$i][0]}{$placed[$i][1]};
		my $e2 = $extent{$placed[$i+1][0]}{$placed[$i+1][1]};
		my $g = $e2->[0] - $e1->[1] - 1;
		($cut, $widest) = ($i, $g) if $g > $widest;
	}
	if(defined $cut){
		my @r = map { [split /\t/] } @full[$cut+1 .. $#full, 0 .. $cut];
		my $wrap = $#full - $cut;                               # where the old first tile now sits
		my $le = $extent{$placed[-1][0]}{$placed[-1][1]};       # the old last tile's reference extent
		my $lo = $le->[1] - $FLIP_REACH; $lo = 1 if $lo < 1;
		$r[$wrap][2] = $lo . "-" . ($le->[1] + $FLIP_REACH);
		$r[0][2] = "-"; $r[-1][3] = "-";
		write_plan("$plan.rotated", map { join("\t", @$_) } @r);
	}
}
open my $u, ">", $unplaced or die "cannot write $unplaced: $!\n";
print $u "$_\n" for @loose;
close $u;
