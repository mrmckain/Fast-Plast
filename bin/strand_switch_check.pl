#!/usr/bin/env perl
use strict;
use warnings;

# Place contigs on a reference plastome and report strand switches inside a
# single-copy region, which afin produces when it joins contigs through a short
# inverted repeat (the segment between the two repeat copies comes out
# inverted). Optionally write a corrected copy of the contigs:
#   - a reversed segment that starts and ends in the two copies of an inverted
#     repeat, or that is bounded by N gaps or contig ends, is turned round in
#     place. Turning it round leaves the sequence at both junctions unchanged
#     (the repeat reads the same from either side), so this is exact;
#   - any other flagged contig is cut at the switch points, so scaffolding can
#     put the pieces back in order.
#
# Usage: strand_switch_check.pl <contigs.fa> <reference.fa> <report.txt> [fixed_out.fa]
#
# blastn and makeblastdb must be on PATH. Work files are written next to the
# contigs: <contigs>.strand_check.blastn, strand_check_ref_db.*,
# strand_check_ref_self.blastn.
#
# Report lines:
#   Reference: <path>
#   Contigs placed: <n>
#   FLAG  <text>     strand switch inside a single-copy region
#   NOTE  <text>     single-copy regions in different orientations (IR isomer)
#   FLIP  <contig> <start> <end>            segment turned round in fixed_out.fa
#   SPLIT <contig> <start> <end> <new id>   piece written to fixed_out.fa
#   OK    no strand switches
#
# fixed_out.fa is only written when at least one contig was changed. It holds
# every contig of the input. Neighbouring pieces of a cut contig keep the
# repeat copy they were joined through, so they overlap and can be spliced.

my $MIN_HSP     = 2000;   # alignment length for a contig placement
my $MIN_IR_HIT  = 1000;   # reference-to-itself, minus strand: part of the IR
my $MIN_PIECE   = 1000;   # placement left inside one single-copy region
my $NEIGHBOUR   = 2000;   # a neighbouring placement may overlap a switch point by this much
my $MIN_END_SEQ = 1000;   # an end piece shorter than this is not written
my $MIRROR_K    = 30;     # bases that must mirror at the two ends of a segment turned round in place
my $SNAP        = 300;    # a segment end this close to an N gap or contig end is moved onto it

my ($contigs, $ref, $report, $fixed_out) = @ARGV;
defined $report or die "Usage: $0 <contigs.fa> <reference.fa> <report.txt> [fixed_out.fa]\n";
-s $contigs or die "ERROR: $contigs is missing or empty\n";
-s $ref     or die "ERROR: $ref is missing or empty\n";

my $db = "strand_check_ref_db";
system("makeblastdb -in \"$ref\" -dbtype nucl -out $db > /dev/null 2>&1") == 0
	or die "ERROR: makeblastdb failed on $ref\n";
my $bl = "$contigs.strand_check.blastn";
system("blastn -query \"$contigs\" -db $db -outfmt \"6 qseqid qlen qstart qend sstart send length slen\" -evalue 1e-50 > \"$bl\" 2>/dev/null") == 0
	or die "ERROR: blastn of $contigs against the reference failed\n";
my $self = "strand_check_ref_self.blastn";
system("blastn -query \"$ref\" -db $db -strand minus -outfmt \"6 qstart qend sstart send length slen\" -evalue 1e-50 > $self 2>/dev/null") == 0
	or die "ERROR: blastn of the reference against itself failed\n";

# ---- reference regions --------------------------------------------------
# IR = everything the reference matches on its own minus strand over at least
# MIN_IR_HIT bp. The single-copy regions are what is left. The reference is a
# circle cut at an arbitrary point, so a single-copy stretch at its start and
# one at its end are the same region.
my $reflen = 0;
my (@ir, @short);
open my $sf, "<", $self or die "ERROR: cannot read $self: $!\n";
while(<$sf>){
	chomp; my @c = split /\t/;
	$reflen = $c[5];
	my @pair = ([$c[0], $c[1]], ($c[2] <= $c[3]) ? [$c[2], $c[3]] : [$c[3], $c[2]]);
	if($c[4] >= $MIN_IR_HIT){ push @ir, @pair; }
	elsif($c[4] >= 100){ push @short, \@pair; }
}
close $sf;
# A shorter hit is IR too when both of its ends continue an IR stretch already
# found. This is the piece of IR cut off where the reference sequence starts.
{
	my $touches = sub {
		my ($iv) = @_;
		for my $k (@ir){
			return 1 if abs($iv->[0] - $k->[1]) <= 50 || abs($k->[0] - $iv->[1]) <= 50;
			return 1 if $iv->[0] <= 50 && $k->[1] >= $reflen - 50;      # across the start point
			return 1 if $iv->[1] >= $reflen - 50 && $k->[0] <= 50;
		}
		return 0;
	};
	my $added = 1;
	while($added){
		$added = 0;
		for my $p (@short){
			next if !$p || !($touches->($p->[0]) && $touches->($p->[1]));
			push @ir, @$p; $p = undef; $added = 1;
		}
	}
}
my @merged;
for my $i (sort { $a->[0] <=> $b->[0] } @ir){
	if(@merged && $i->[0] <= $merged[-1][1] + 50){ $merged[-1][1] = $i->[1] if $i->[1] > $merged[-1][1]; }
	else{ push @merged, [@$i]; }
}
my @sc;      # [lo, hi, region number]
{
	my $pos = 1;
	for my $i (@merged, [$reflen + 1, $reflen + 1]){
		push @sc, [$pos, $i->[0] - 1] if $i->[0] - $pos >= $MIN_PIECE;
		$pos = $i->[1] + 1;
	}
	my $n = 0;
	$_->[2] = $n++ for @sc;
	if(@sc > 1 && $sc[0][0] == 1 && $sc[-1][1] == $reflen){ $sc[-1][2] = $sc[0][2]; }
}
my %region_len;
$region_len{$_->[2]} += $_->[1] - $_->[0] + 1 for @sc;

# ---- contig placements --------------------------------------------------
my (%hits, %qlen);
open my $in, "<", $bl or die "ERROR: cannot read $bl: $!\n";
while(<$in>){
	chomp; my @c = split /\t/;
	$reflen ||= $c[7];
	next if $c[6] < $MIN_HSP;
	my ($qs, $qe, $ss, $se) = @c[2..5];
	my $str = ($ss <= $se) ? "+" : "-";
	($ss, $se) = ($se, $ss) if $ss > $se;
	push @{$hits{$c[0]}}, { qs => $qs, qe => $qe, ss => $ss, se => $se, str => $str };
	$qlen{$c[0]} = $c[1];
}
close $in;
@sc = ([1, $reflen, 0]) and $region_len{0} = $reflen if !@sc && $reflen;   # no IR found in the reference

sub covered {
	# bases of the query covered by a list of pieces (union of intervals)
	my ($cov, $lo, $hi) = (0);
	for my $p (sort { $a->{qs} <=> $b->{qs} } @_){
		if(defined $hi && $p->{qs} <= $hi){ $hi = $p->{qe} if $p->{qe} > $hi; }
		else{ $cov += $hi - $lo + 1 if defined $hi; ($lo, $hi) = ($p->{qs}, $p->{qe}); }
	}
	$cov += $hi - $lo + 1 if defined $hi;
	return $cov;
}

my (%seq, @order);
{
	my $id;
	open my $cf, "<", $contigs or die "ERROR: cannot read $contigs: $!\n";
	while(<$cf>){
		chomp;
		if(/^>(\S+)/){ $id = $1; push @order, $id; $seq{$id} = ""; }
		elsif(defined $id){ $seq{$id} .= $_; }
	}
	close $cf;
}

my (@flags, @notes, %cuts, @flip_lines, @split_lines, $changed);
for my $q (sort keys %hits){
	my @hs = @{$hits{$q}};
	# Cut each placement down to the part inside each single-copy region. What
	# falls in the IR says nothing about strand: it matches both copies.
	my %pieces;     # region -> strand -> pieces
	for my $h (@hs){
		for my $r (@sc){
			my $lo = $h->{ss} > $r->[0] ? $h->{ss} : $r->[0];
			my $hi = $h->{se} < $r->[1] ? $h->{se} : $r->[1];
			next if $hi - $lo + 1 < $MIN_PIECE;
			my ($qs, $qe) = $h->{str} eq "+"
				? ($h->{qs} + ($lo - $h->{ss}), $h->{qe} - ($h->{se} - $hi))
				: ($h->{qs} + ($h->{se} - $hi), $h->{qe} - ($lo - $h->{ss}));
			$qs = $h->{qs} if $qs < $h->{qs};
			$qe = $h->{qe} if $qe > $h->{qe};
			next if $qe <= $qs;
			push @{$pieces{$r->[2]}{$h->{str}}}, { qs => $qs, qe => $qe, ss => $lo, se => $hi, hsp => $h };
		}
	}
	my %major;
	my @minor_hsps;
	for my $reg (sort { $a <=> $b } keys %pieces){
		my $plus  = covered(@{$pieces{$reg}{"+"} || []});
		my $minus = covered(@{$pieces{$reg}{"-"} || []});
		$major{$reg} = ($plus >= $minus) ? "+" : "-";
		next unless $plus && $minus;
		my $minor = ($plus >= $minus) ? "-" : "+";
		my @seg = sort { $a->{qs} <=> $b->{qs} } @{$pieces{$reg}{$minor}};
		my $desc = join("; ", map { sprintf("contig %d-%d = reference %d-%d", $_->{hsp}{qs}, $_->{hsp}{qe}, $_->{hsp}{ss}, $_->{hsp}{se}) } @seg);
		push @flags, "$q: reversed segment(s) $desc, while the rest of the same single-copy region maps to the other strand. Likely an afin misjoin at a repeat; check this contig before trusting the assembly.";
		for my $s (@seg){ $s->{hsp}{minor} = 1; push @minor_hsps, $s->{hsp}; }
	}
	# Single-copy regions facing different ways: the largest sets the direction.
	my @regs = sort { $region_len{$b} <=> $region_len{$a} } keys %major;
	for my $reg (@regs[1..$#regs]){
		next if $major{$reg} eq $major{$regs[0]};
		push @notes, "$q: a single-copy region (" . join(", ", map { "reference $_->[0]-$_->[1]" } grep { $_->[2] == $reg } @sc)
			. ") is in the opposite orientation to the largest one. This is a normal IR isomer, not an error.";
	}

	next unless @minor_hsps;
	my %seen;
	my @segs = grep { !$seen{"$_->{qs}-$_->{qe}"}++ } sort { $a->{qs} <=> $b->{qs} } @minor_hsps;
	my @rest = grep { !$_->{minor} } @hs;
	next unless $fixed_out && exists $seq{$q};
	my $S = $seq{$q};
	my $L = length $S;

	# First choice: turn each reversed segment round in place.
	my @flip;
	for my $s (@segs){
		my ($a, $b) = ($s->{qs}, $s->{qe});
		# move an end onto a nearby N gap or contig end
		my $lw = substr($S, ($a > $SNAP ? $a - 1 - $SNAP : 0), ($a > $SNAP ? $SNAP : $a - 1));
		if($a <= $SNAP && $lw !~ /N/i){ $a = 1; }
		elsif($lw =~ /N[^Nn]*$/i){ $a -= length($lw) - $+[0] + length($&) - 1; }
		my $rw = substr($S, $b, $SNAP);
		if($b + $SNAP >= $L && $rw !~ /N/i){ $b = $L; }
		elsif($rw =~ /^([^Nn]*)N/i){ $b += length($1); }
		my $free_l = ($a == 1  || substr($S, $a - 2, 1) =~ /N/i);
		my $free_r = ($b == $L || substr($S, $b, 1)     =~ /N/i);
		if($free_l && $free_r){ push @flip, [$a, $b]; next; }
		last if $free_l || $free_r;          # one loose end only: cut instead
		# both ends inside the contig: they must sit in the two copies of an
		# inverted repeat, at mirrored positions
		my $found;
		SEARCH: for my $d (0..10){
			for my $da (-$d..$d){
				for my $db (-($d - abs($da)), ($d - abs($da))){
					my ($x, $y) = ($a + $da, $b + $db);
					next if $x < 1 || $y > $L || $y - $x + 1 < 2 * $MIRROR_K;
					my $head = substr($S, $x - 1, $MIRROR_K);
					my $tail = reverse substr($S, $y - $MIRROR_K, $MIRROR_K);
					$tail =~ tr/ACGTacgt/TGCAtgca/;
					if(uc $head eq uc $tail && $head !~ /N/i){ $found = [$x, $y]; last SEARCH; }
				}
			}
		}
		last unless $found;
		push @flip, $found;
	}
	if(@flip == @segs){
		for my $f (@flip){
			my ($a, $b) = @$f;
			my $rc = reverse substr($S, $a - 1, $b - $a + 1);
			$rc =~ tr/ACGTacgtRYKMrykmBVDHbvdh/TGCAtgcaYRMKyrmkVBHDvbdh/;
			substr($S, $a - 1, $b - $a + 1) = $rc;
			push @flip_lines, "$q\t$a\t$b";
		}
		$seq{$q} = $S;
		$changed = 1;
		next;
	}

	# Otherwise cut at the switch points. Each reversed placement becomes its
	# own piece; the pieces around it run up to the end of the neighbouring
	# placement, which reaches through the repeat copy the two were joined in.
	my @out;
	my $start = 1;
	for my $s (@segs){
		my $left_end = $s->{qs} - 1;
		for my $h (@rest){
			next unless $h->{qs} < $s->{qs} && $h->{qe} >= $s->{qs} - 10 && $h->{qe} <= $s->{qs} + $NEIGHBOUR;
			$left_end = $h->{qe} if $h->{qe} > $left_end;
		}
		push @out, [$start, $left_end] if $left_end >= $start && $s->{qs} > $start;
		push @out, [$s->{qs}, $s->{qe}];
		$start = $s->{qe} + 1;
		for my $h (@rest){
			next unless $h->{qe} > $s->{qe} && $h->{qs} <= $s->{qe} + 10 && $h->{qs} >= $s->{qe} - $NEIGHBOUR;
			$start = $h->{qs} if $h->{qs} < $start;
		}
	}
	push @out, [$start, $L] if $start <= $L;
	# an end piece too short to place is left out
	shift @out if @out > 1 && $out[0][1]  - $out[0][0]  + 1 < $MIN_END_SEQ;
	pop   @out if @out > 1 && $out[-1][1] - $out[-1][0] + 1 < $MIN_END_SEQ;
	if(@out > 1){ $cuts{$q} = \@out; $changed = 1; }
}

# ---- corrected FASTA ----------------------------------------------------
if($fixed_out && $changed){
	open my $so, ">", $fixed_out or die "ERROR: cannot write $fixed_out: $!\n";
	for my $c (@order){
		if(!$cuts{$c}){ print $so ">$c\n$seq{$c}\n"; next; }
		my $k = 0;
		for my $p (@{$cuts{$c}}){
			$k++;
			my ($s, $e) = @$p;
			my $len = $e - $s + 1;
			my $new = ($c =~ /^(.*?)_len_\d+(.*)$/) ? "$1s${k}_len_$len$2" : "${c}_s${k}_len_$len";
			print $so ">$new\n" . substr($seq{$c}, $s - 1, $len) . "\n";
			push @split_lines, "$c\t$s\t$e\t$new";
		}
	}
	close $so;
}

open my $rep, ">", $report or die "ERROR: cannot write $report: $!\n";
print $rep "Reference: $ref\n";
print $rep "Contigs placed: " . scalar(keys %hits) . "\n";
print $rep "FLAG\t$_\n"  for @flags;
print $rep "NOTE\t$_\n"  for @notes;
print $rep "FLIP\t$_\n"  for @flip_lines;
print $rep "SPLIT\t$_\n" for @split_lines;
print $rep "OK\tno strand switches\n" unless @flags || @notes;
close $rep;
