#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);

use lib (".");
use File::Basename;
use lib dirname(__FILE__);
require "kabsch.pm";
require "matrix.pm";

my $dssp                 = "/work/dimaio/bin/dssp";
my $RMS_ALIGN = 2.5;

###############################
###############################

my %one_to_three = ( 
	'G' => 'GLY', 'A' => 'ALA', 'V' => 'VAL', 'L' => 'LEU',
	'I' => 'ILE', 'P' => 'PRO', 'C' => 'CYS', 'M' => 'MET',
	'H' => 'HIS', 'F' => 'PHE', 'Y' => 'TYR', 'W' => 'TRP',
	'N' => 'ASN', 'Q' => 'GLN', 'S' => 'SER', 'T' => 'THR',
	'K' => 'LYS', 'R' => 'ARG', 'D' => 'ASP', 'E' => 'GLU' );

###############################
###############################

if ($#ARGV < 1) {
	print STDERR "usage: $0 <fastafile> <pdb>*\n";
	exit -1;
}

# read sequence
my $fastafile = shift @ARGV;
open (FASTA, $fastafile) || print STDERR "Cannot open $_";
my @fastalines = <FASTA>;
close FASTA;
my $seq="";
foreach my $line (@fastalines) {
	chomp $line;
	if ($line !~ /^>/) {
		$seq = $seq.$line;
	}
}
$seq =~ s/\s//g; #remove whitespace
my $nres = length( $seq );
print STDERR "Read sequence: $seq\n";

my $currfrag_lines;
my $currfrag_ids;
my %allfrags;
my %bbatoms;
my %resids;
my %fraglens;
my %ss;

# read in all PDBs
my $minres = 9999;
my $maxres = 0;
foreach my $pdb (@ARGV) {
	open (PDB, $pdb) || print STDERR "Cannot open $_";

	my $start_new_frag = 0;
	$bbatoms{$pdb} = {};
	$resids{$pdb} = {};
	$currfrag_lines = [];
	$currfrag_ids = [];
	$allfrags{$pdb} = [];

	# get SS parse
	%ss = getSecStruct( $pdb );

	my $curr_frag_start = -999;
	my $last_res_read = -999;
	while (my $line = <PDB>) {
		if ($line =~ /^TER/) {
			$start_new_frag = 1;
		}
		next if (! $line =~ /^ATOM/);

		my $atom = substr ($line, 12, 4);
		my $chain = substr($line, 21, 1);
		my $resid = substr($line, 17, 3);
		my $res = int( substr($line, 22, 4) );
		$minres = min($minres,$res);
		$maxres = max($maxres,$res);

		if ($curr_frag_start == -999) {
			$curr_frag_start = $minres;
		}

		# only care about bb heavy atoms + CB
		next if ($atom ne " CA " && $atom ne " C  " && $atom ne " O  " && $atom ne " N  " && $atom ne " CB ");

		my $x = substr ($line, 30, 8);
		my $y = substr ($line, 38, 8);
		my $z = substr ($line, 46, 8);

		my $id = $atom.$chain.$res;
		$bbatoms{$pdb}->{ $id } = [$x,$y,$z];
		$resids{$pdb}->{ $id } = $resid;

		# check for chainbreaks
		if ($atom eq " N  ") {
			if ($start_new_frag != 1) { # no TER just read
				my $previd = " C  ".$chain.($res-1);
				if (!defined $bbatoms{$pdb}->{ $previd }) { # no previous residue
					$start_new_frag = 1;
				} else {
					my $dist = dist($bbatoms{$pdb}->{ $id },$bbatoms{$pdb}->{ $previd });
					if ($dist > 2.5) { # ????
						$start_new_frag = 1;
					} else {
						$start_new_frag = 0;
					}
				}
			}
		} else {
			$start_new_frag = 0;
		}

		if ($start_new_frag && scalar @{ $currfrag_lines } > 0) {
			push @{ $allfrags{$pdb} }, $currfrag_lines;
			foreach my $id_i (@{ $currfrag_ids }) {
				$fraglens{$pdb}->{ $id_i } = scalar @{ $currfrag_ids };
			}

			# make subfrags
	 		my $in_helix=0; my $in_strand=0; my $ss_start=1;
			my @ss_elts = ();
			for (my $ii=$curr_frag_start; $ii<=$last_res_read+1; ++$ii) {
				if (!defined $ss{$ii} || $ss{$ii} ne 'E') {
					if ($in_strand==1) { # strand end
						$in_strand=0;
						if ($ii-$ss_start >= 3) {
							push @ss_elts, [$ss_start, $ii-1];
						}
					}
				}
				if (!defined $ss{$ii} || $ss{$ii} ne 'H') {
					if ($in_helix == 1) { # helix end
						$in_helix = 0;
						if ($ii-$ss_start >= 5) {
							push @ss_elts, [$ss_start, $ii-1];
						}
					}
				}
				last if (!defined $ss{$ii});
				if ($ss{$ii} eq 'E' && $in_strand==0)  { # strand start
					$in_strand = 1;
					$ss_start = $ii;
				}
				if ($ss{$ii} eq 'H' && $in_helix==1)  { # helix start
					$in_helix = 1;
					$ss_start = $ii;
				}
	 		}
			if (scalar(@ss_elts) > 1) {
				print STDERR "SUBDIVIDING $pdb fragment [".$curr_frag_start.",".$last_res_read."] into ".scalar(@ss_elts)." pieces\n";
				my $prev_cut = $curr_frag_start;
				foreach my $ii (0..$#ss_elts) {
					my $cut;
					if ($ii == $#ss_elts) {
						$cut = $last_res_read;
					} else {
						$cut = floor( 0.5 * ($ss_elts[$ii]->[1] + $ss_elts[$ii+1]->[0]) + 0.5 );
					}
					print STDERR "   fragment $ii [".$prev_cut.",".$cut."]\n";
					my $subfrag_lines = [];
					foreach my $miniline (@{ $currfrag_lines }) {
						my $minires = int( substr($miniline, 22, 4) );
						if ($minires >=$prev_cut && $minires <= $cut) {
							push @{ $subfrag_lines }, $miniline;
						}
					}
					push @{ $allfrags{$pdb} }, $subfrag_lines;
					$prev_cut=$cut+1;
				}
			}
			$curr_frag_start = $res;
			$currfrag_lines = [];
			$currfrag_ids = [];
		}

		$last_res_read = $res;

		push @{ $currfrag_lines }, $line;
		push @{ $currfrag_ids }, $id;
	}

	if (scalar(@{ $currfrag_lines }) > 0) {
		push @{ $allfrags{$pdb} }, $currfrag_lines;
		foreach my $id_i (@{ $currfrag_ids }) {
			$fraglens{$pdb}->{ $id_i } = scalar @{ $currfrag_ids };
		}
		$currfrag_lines = [];
		$currfrag_ids = [];
	}

	print STDERR "Built ".scalar( @{ $allfrags{$pdb} } )." fragments from ".$pdb."\n";

	close (PDB);
}

# align all->all
my $npdbs = @ARGV;
my $rmsds = [];
my $nalignlen = [];
my $Rs = [];
my $comis = [];
my $comjs = [];

foreach my $i (0..$#ARGV) {
foreach my $j ($i..$#ARGV) {
	# intersection of $bbatoms{i} and {j}
	my @common_atoms = ();
	foreach ( keys %{ $bbatoms{$ARGV[$i]} } ) {
		push @common_atoms, $_ if exists $bbatoms{$ARGV[$j]}->{$_};
	}

	# atom lists
	my $atoms_i = [];
	my $atoms_j = [];
	foreach ( @common_atoms ) {
		push @{ $atoms_i }, deep_copy( $bbatoms{$ARGV[$i]}->{$_} );
		push @{ $atoms_j }, deep_copy( $bbatoms{$ARGV[$j]}->{$_} );
	}

	# realign
	my $done = 0;
	while ($done == 0 && scalar(@common_atoms) > 4) {
		($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
		$nalignlen->[$i][$j] = scalar(@common_atoms);

		$done = 1;
		my @new_common_atoms = ();
		$atoms_i = [];
		$atoms_j = [];
		foreach ( @common_atoms ) {
			my $x = $bbatoms{$ARGV[$i]}->{$_};
			my $y = vadd( mapply( $Rs->[$i][$j] , vsub( $bbatoms{$ARGV[$j]}->{$_} , $comis->[$i][$j] ) ), 
			              vadd($comis->[$i][$j] , $comjs->[$i][$j]) );
			my $dist = dist( $x,$y );
			if ($dist > $RMS_ALIGN) {
				$done = 0;
			} else {
				push @new_common_atoms, $_;
				push @{ $atoms_i }, deep_copy( $bbatoms{$ARGV[$i]}->{$_} );
				push @{ $atoms_j }, deep_copy( $bbatoms{$ARGV[$j]}->{$_} );
			}
		}
		@common_atoms = @new_common_atoms;
		if ($done == 1) {print STDERR "RMS(".$ARGV[$i].",".$ARGV[$j].") over ".@common_atoms." atoms is ".$rmsds->[$i][$j]."\n";}
	}
}
}


foreach my $i (0..$#ARGV) {
foreach my $j (0..$i-1) {
	 $Rs->[$i][$j] =  minv( $Rs->[$j][$i] );
	 $rmsds->[$i][$j] =  deep_copy( $rmsds->[$j][$i] );
	 $comis->[$i][$j] =  vadd( $comis->[$j][$i], $comjs->[$j][$i] );
	 $comjs->[$i][$j] =  vscale( -1 , $comjs->[$j][$i] );
	 $nalignlen->[$i][$j] = $nalignlen->[$j][$i];
}
}

# greedy align

# find min RMS in table
my $maxlen = -999;
my $minI = -1;
my $minJ = -1;
foreach my $i (0..$#ARGV) {
foreach my $j ($i+1..$#ARGV) {
	if ( $nalignlen->[$i][$j] > $maxlen ) {
		$maxlen = $nalignlen->[$i][$j];
		$minI = $i; $minJ = $j;
	}
}
}

# use i's coordinate frame
# seed j
my $globalRs = [];
$globalRs->[$minI] = [[1,0,0],[0,1,0],[0,0,1]];
$globalRs->[$minJ] = deep_copy( $Rs->[$minI][$minJ] );
my $preTs = [];
$preTs->[$minI] = deep_copy( $comis->[$minJ][$minI] );
$preTs->[$minJ] = deep_copy( $comis->[$minI][$minJ] );
my $postTs = [];
$postTs->[$minI] = $preTs->[$minI];
$postTs->[$minJ] = vadd(  $comjs->[$minI][$minJ], $preTs->[$minJ] );

# for i = 2 to $#ARGV find best alignment to already-aligned model
my @aligned = (0) x scalar(@ARGV);
$aligned[$minI] = 1; $aligned[$minJ] = 1;

print STDERR "align ".$ARGV[$minI].",".$ARGV[$minJ]."\n";

foreach my $cycle (2..$#ARGV) {
	$maxlen = -999; $minI = -1; $minJ = -1;
	foreach my $i (0..$#ARGV) {
	foreach my $j (0..$#ARGV) {
		next if ($i==$j);
		next if ($aligned[$i] == 0 || $aligned[$j] == 1);
		# i is aligned, j is not
		if ( $nalignlen->[$i][$j] > $maxlen ) {
			$maxlen = $nalignlen->[$i][$j];
			$minI = $i; $minJ = $j;
		}
	}
	}

	# add j to alignment
	print STDERR "align ".$ARGV[$minJ]." (from".$ARGV[$minI].")\n";
	$globalRs->[$minJ] = mmult( $globalRs->[$minI], $Rs->[$minI][$minJ]  );
	$preTs->[$minJ] = deep_copy( $comis->[$minI][$minJ] );

	# post-rotation translation is a bit tricky
	# we need to transform to j's coordinate frame, 
	#   then apply j's transformation to this result
	my $post_to_j = vadd( $comjs->[$minI][$minJ], $preTs->[$minJ] );
	my $post_j_to_global = vadd( mapply( $globalRs->[$minI] , vsub( $post_to_j , $preTs->[$minI]) ) , 
	                             $postTs->[$minI] );
	$postTs->[$minJ] = $post_j_to_global;

	$aligned[$minJ] = 1;
}


# apply rotation; dump primary fragments to DB
mkdir "fragments";
mkdir "aligned_templates";
foreach my $i (0..$#ARGV) {
	my $pdb = $ARGV[$i];
	my $nfrags = scalar( @{ $allfrags{$pdb} } );

	my $outpdb = $pdb;
	$outpdb =~ s/.*\///;
	$outpdb =~ s/\.pdb$/_aln.$i.pdb/;
	print STDERR "writing aligned_templates/$outpdb\n";
	open (PDBF, ">aligned_templates/$outpdb") || print STDERR "Cannot open $_";

	foreach my $frag (0..$nfrags-1) {
		my $outfile = $pdb;
		$outfile =~ s/.*\///;
		$outfile =~ s/\.pdb$/_$frag.$i.pdb/;
		open (PDB, ">fragments/$outfile") || print STDERR "Cannot open $_";
		foreach my $line (@{ $allfrags{$pdb}->[$frag] }) {
			my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
			$newX = vsub( $newX, $preTs->[$i]);
			$newX = mapply( $globalRs->[$i], $newX );
			$newX = vadd( $newX, $postTs->[$i]);

			substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
			print PDB $line;
			print PDBF $line;
		}
		close PDB;
	}
	close PDBF;
}


# prepare pseudo-backbone trace
my @pseudotrace;
my @gaps;
my $last_ungapped = 0;
foreach my $i (1..$nres) {
	my $id = " CA "."A".$i;

	# (opt 1) find longest fragment at each position
#	my $maxfraglen = 0;
#	my $maxid = -1;
#	foreach my $j (0..$#ARGV) {
#		my $pdb = $ARGV[$j];
#		if ( defined $fraglens{$pdb}->{$id} && $fraglens{$pdb}->{$id} > $maxfraglen) {
#			$maxfraglen = $fraglens{$pdb}->{$id};
#			$maxid = $j;
#		}
#	}

	# (opt 2) find residue closest to mean
 	my @cas_i = ();
 	my @ids_i = ();
	my $maxid = -1;
	foreach my $j (0..$#ARGV) {
 		my $pdb = $ARGV[$j];
 		if (defined $bbatoms{$pdb}->{$id} ) {
 			push @ids_i, $j;
			my $newX = deep_copy($bbatoms{$pdb}->{$id});
			$newX = vsub( $newX, $preTs->[$j]);
			$newX = mapply( $globalRs->[$j], $newX );
			$newX = vadd( $newX, $postTs->[$j]);
 			push @cas_i, $newX;
 		}
 	}
	my $com = [0,0,0];
	if (scalar(@cas_i) > 0) {
		foreach (@cas_i) { $com = vadd( $com, $_ ); }
		$com = vscale( 1/scalar(@cas_i), $com);
		my $mindist = 999;
		my $minId;
		foreach my $j (0..$#cas_i) {
			my $thisdist = vdist( $cas_i[$j], $com );
			if ($mindist > $thisdist) {
				$mindist = $thisdist;
				$maxid = $ids_i[$j];
			}
		}
	}
 
	# add it to pseudotrace
	if ($maxid >= 0) {
		$gaps[$i] = 0;
		$last_ungapped = $i;
		foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
			my $x = $bbatoms{ $ARGV[$maxid] }->{$atm."A".$i};
			$x = vsub( $x, $preTs->[$maxid]);
			$x = mapply( $globalRs->[$maxid], $x );
			$x = vadd( $x, $postTs->[$maxid]);

			$pseudotrace[$i]->{$atm} = $x;
		}
	} else {
		$gaps[$i] = 999;
		if ($i>1) {
			$gaps[$i] = $gaps[$i-1]+1; # right extension distance
		}
	}
}

# fill in gaps
my $maxgap = $gaps[$nres];
for (my $i=$last_ungapped-1; $i >= 1; $i--) {
	$gaps[$i] = min( $gaps[$i+1]+1, $gaps[$i] );
	$maxgap = max($maxgap, $gaps[$i]);
}
foreach my $i (1..$nres) { print STDERR $gaps[$i]; } print STDERR "\n";
for my $cycle (1..$maxgap) {
	# extend
	foreach my $i (1..$nres) {
		if ($gaps[$i] == $cycle) {
			print STDERR "Rebuild $i...\n";
			my $dir=-1;
			$dir=1 if ($i > $nres-2 || $gaps[$i+1] < $cycle);

			my $Xp1=[];
			my $Xp2=[];
			foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
				push @{$Xp1}, deep_copy($pseudotrace[$i-1*$dir]->{$atm});
				push @{$Xp2}, deep_copy($pseudotrace[$i-2*$dir]->{$atm});
			}

			# align xp2->xp1
			my ($R21, $rmsd21, $com2, $com1) = rms_align( $Xp2 , $Xp1 );

			# apply to xp1
			foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
				my $newX = vsub( $pseudotrace[$i-1*$dir]->{$atm}, $com2);
				$newX = mapply( $R21, $newX );
				$newX = vadd( $newX, $com2);
				$newX = vadd( $newX, $com1);
				$pseudotrace[$i]->{$atm} = $newX;
			}
		}
	}
}
	

# write starting model
my $atmidx = 1;
foreach my $i (1..$nres) {
	my $restype = $one_to_three{ substr($seq, $i-1, 1) };
	foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
		my $x = $pseudotrace[$i]->{$atm};
		printf "ATOM   %4d %s %s %s %3d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
		       $atmidx++, $atm, $restype, 'A', $i, $x->[0], $x->[1], $x->[2];
	}
}




exit 0;


###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}

#sub max ($$) { $_[$_[0] < $_[1]] }
#sub min ($$) { $_[$_[0] > $_[1]] }

sub mapSScode {
    my $incode = shift;
    $incode = uc $incode;
    my $newcode = undef;

    my %mapping = ( 'H' => 'H', 'G' => 'H', 'E' => 'E', 'B' => 'E' );

    return (defined $mapping{$incode}) ? $mapping{$incode} : 'L';
}

sub getSecStruct {
	my $pdbfile = $_[0];
	my $dsspfile = "$pdbfile.dssp";
	`$dssp -na $pdbfile $dsspfile &> /dev/null`;
	
	# read dssp
	open (FILE, $dsspfile) || die "$0: unable to open file $dsspfile for reading";
	my @dssp_buf = <FILE>;
    close (FILE);	

	my $reading = undef;
	my %ss;
	for (my $i=0; $i <= $#dssp_buf; ++$i) {
		if ($dssp_buf[$i] =~ /^  \#  RESIDUE AA STRUCTURE BP1 BP2  ACC/) {
			$reading = 'true';
			next;
		}
		next if (! $reading);
	
		my $res_n   = substr ($dssp_buf[$i],  6, 4);
		my $aa      = substr ($dssp_buf[$i], 13, 1);
		my $dssp_ss = substr ($dssp_buf[$i], 16, 1);

		next if ($res_n =~ /^\s*$/ || $aa eq '!' || $aa =~ /[a-z]/);

		$ss{int($res_n)} = &mapSScode ($dssp_ss);
	}

	return %ss;
}
