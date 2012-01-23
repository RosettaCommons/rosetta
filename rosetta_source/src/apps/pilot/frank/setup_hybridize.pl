#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);

use lib (".");
use File::Basename;
use Cwd;
use lib dirname(__FILE__);
require "kabsch.pm";
require "matrix.pm";

###############################
###############################
## CST files
my $CSTFILENAMES = "../align_new/cluster%d.filt.dist_csts";
#my $CSTFILENAMES = "NONE";

###############################
###############################
## alignment parameters
my @RMS_CUTOFFS = (10,5,4,3,2,1.5,1);
my $RMSCUTOFF = 999;

my $CLUSTERCUTOFF = 0.40;
my $ALIGNCUTOFF   = 0.2;  # to get better superpositions, trade coverage for alignment


###############################
###############################
## probability correct
my %template_probs =  (
	401 => 0.16318, 201 => 0.15481, 301 => 0.13389,
	302 => 0.09623, 303 => 0.05858, 304 => 0.05439,
	403 => 0.05021, 402 => 0.03766, 202 => 0.03347,
	305 => 0.02929, 203 => 0.02510, 404 => 0.02092,
	405 => 0.01674, 309 => 0.01255, 308 => 0.01255,
	409 => 0.01255, 408 => 0.01255, 407 => 0.01255,
	410 => 0.01255, 204 => 0.01255, 306 => 0.00837,
	310 => 0.00837, 205 => 0.00837, 307 => 0.00418,
	406 => 0.00418, 206 => 0.00418,
	0 => 0.00209    # default
);

###############################
###############################
## amino acids
my %one_to_three = ( 
	'G' => 'GLY', 'A' => 'ALA', 'V' => 'VAL', 'L' => 'LEU',
	'I' => 'ILE', 'P' => 'PRO', 'C' => 'CYS', 'M' => 'MET',
	'H' => 'HIS', 'F' => 'PHE', 'Y' => 'TYR', 'W' => 'TRP',
	'N' => 'ASN', 'Q' => 'GLN', 'S' => 'SER', 'T' => 'THR',
	'K' => 'LYS', 'R' => 'ARG', 'D' => 'ASP', 'E' => 'GLU' );

###############################
###############################
## main()
if ($#ARGV < 2) {
	print STDERR "usage: $0 <fasta-file> <cluster-file> <pdb>*\n";
	exit -1;
}

my $fastafile = shift @ARGV;
my $clusterfile = shift @ARGV;

## read sequence
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


## read cluster file to get seeds
## remember which constraint file each refers to
my %james_clustermap;
my $firstcluster;
open (CL, $clusterfile) || print STDERR "Cannot open $_";
my @cllines = <CL>;
close CL;
foreach my $line (@cllines) {
	chomp $line;
	my @fields = split / /, $line;
	if ($#fields > 0) {
		if (scalar keys %james_clustermap == 0) {
			$firstcluster = $fields[0];
		}
		foreach my $i (1..$#fields) {
			$james_clustermap{ $fields[$i] } = $fields[0];
		}
	}
}

print STDERR "First cluster : $firstcluster\n";

# initialize "james cluster" assignments
my @james_clusterid = (-1) x ($#ARGV+1);
foreach my $member (keys %james_clustermap) {
	my $membertag = $member;
	$membertag =~ s/\.pdb//;
	foreach my $i (0..$#ARGV) {
		if ($ARGV[$i] =~ /$membertag/) {
			$james_clusterid[$i] = $james_clustermap{ $member };
		}
	}
}
print STDERR "Cluster assignments: ";
print STDERR join ' ',@james_clusterid;
print STDERR "\n";


##
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
	print STDERR $pdb."\n";
	open (PDB, $pdb) || print STDERR "Cannot open $pdb";

	my $start_new_frag = 0;
	$bbatoms{$pdb} = {};
	$resids{$pdb} = {};
	$currfrag_lines = [];
	$currfrag_ids = [];
	$allfrags{$pdb} = [];

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

############
# align all->all
############

my $npdbs = @ARGV;
my $rmsds = [];
my $nalignlen = [];
my $Rs = [];
my $comis = [];
my $comjs = [];

my $overlapscore = []; # for clustering

foreach my $i (0..$#ARGV) {
foreach my $j ($i..$#ARGV) {
	# intersection of $bbatoms{i} and {j}
	my @common_atoms = ();
	foreach ( keys %{ $bbatoms{$ARGV[$i]} } ) {
		push @common_atoms, $_ if exists $bbatoms{$ARGV[$j]}->{$_};
	}
	my $ncommon_atoms = scalar( @common_atoms );

	if ($ncommon_atoms < 6) {
		$rmsds->[$i][$j] = 9999;
		$Rs->[$i][$j] = [[1,0,0],[0,1,0],[0,0,1]];
		$nalignlen->[$i][$j] = scalar(@common_atoms);
		next;
	}

	# atom lists
	my $atoms_i = [];
	my $atoms_j = [];
	foreach ( @common_atoms ) {
		push @{ $atoms_i }, deep_copy( $bbatoms{$ARGV[$i]}->{$_} );
		push @{ $atoms_j }, deep_copy( $bbatoms{$ARGV[$j]}->{$_} );
	}

	# initial alignment
	($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
	$nalignlen->[$i][$j] = scalar(@common_atoms);

	next if ($i==$j);

	# realign after trimming outliers
	foreach my $RMS_ALIGN (@RMS_CUTOFFS) {
		my @new_common_atoms = ();
		$atoms_i = [];
		$atoms_j = [];
		foreach ( @common_atoms ) {
			my $x = $bbatoms{$ARGV[$i]}->{$_};
			my $y = vadd( mapply( $Rs->[$i][$j] , vsub( $bbatoms{$ARGV[$j]}->{$_} , $comis->[$i][$j] ) ), 
			              vadd($comis->[$i][$j] , $comjs->[$i][$j]) );
			my $dist = dist( $x,$y );
			if ($dist <= $RMS_ALIGN) {
				push @new_common_atoms, $_;
				push @{ $atoms_i }, deep_copy( $bbatoms{$ARGV[$i]}->{$_} );
				push @{ $atoms_j }, deep_copy( $bbatoms{$ARGV[$j]}->{$_} );
			}
		}

		my $natoms_tot = min( scalar(keys %{ $bbatoms{$ARGV[$i]} }),  scalar(keys %{ $bbatoms{$ARGV[$j]} }) );
		if ($RMS_ALIGN == $RMS_CUTOFFS[1]) {
			#$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $ncommon_atoms;
			$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $natoms_tot;
			print STDERR "overlap(".$ARGV[$i].",".$ARGV[$j].") over ".@common_atoms."/".$natoms_tot." atoms is ".$overlapscore->[$i][$j]."\n";
		}

		if ( scalar( @new_common_atoms ) > $natoms_tot*$ALIGNCUTOFF ) {
			@common_atoms = @new_common_atoms;
			($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
			$nalignlen->[$i][$j] = scalar(@common_atoms);
		} else {
			#print STDERR "RMS(".$ARGV[$i].",".$ARGV[$j].") over ".@common_atoms." atoms is ".$rmsds->[$i][$j]."\n";
			last;
		}
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
	 $overlapscore->[$i][$j] = $overlapscore->[$j][$i];
}
}

# clustering
my $clusterid = [];


# start with highest probability member of first cluster as the "seed"
my $maxprob = 0;
my $minI = -1;
foreach my $i (0..$#ARGV) {
	if ($james_clusterid[$i] == $firstcluster) {
		my $tag = $ARGV[$i]; $tag =~ s/.*_(\d\d\d).*/$1/;
		my $prob = $template_probs{ $tag };
		$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
		if ($prob > $maxprob) {
			$maxprob = $prob;
			$minI = $i;
		}
	}
}
print STDERR "SEED $minI (".$ARGV[$minI].") with probability $maxprob\n";


my $nextclusterid = 1;
$clusterid->[$minI] = $nextclusterid;
$nextclusterid++;

# use i's coordinate frame
my $globalRs = [];
$globalRs->[$minI] = [[1,0,0],[0,1,0],[0,0,1]];
my $preTs = [];
$preTs->[$minI] = deep_copy( $comis->[$minI][$minI] );
my $postTs = [];
$postTs->[$minI] = $preTs->[$minI];

# for remaining structs find best alignment to already-aligned model
my @aligned = (0) x scalar(@ARGV);
$aligned[$minI] = 1;

my $seedPDB = $minI;

my $minJ;
my $minRMS;
my $maxOverlap;
foreach my $cycle (1..$#ARGV) {
	$minRMS = 999;
	$maxOverlap = 0;
	$minI = -1; $minJ = -1;
	foreach my $i (0..$#ARGV) {
	foreach my $j (0..$#ARGV) {
		next if ($i==$j);
		next if ($aligned[$i] == 0 || $aligned[$j] == 1);
		# i is aligned, j is not
		#if ( $rmsds->[$i][$j] < $minRMS ) {
		#	$minRMS = $rmsds->[$i][$j];
		#	$minI = $i; $minJ = $j;
		#}
		if ( $overlapscore->[$i][$j] > $maxOverlap ) {
			$maxOverlap = $overlapscore->[$i][$j];
			$minI = $i; $minJ = $j;
		}
	}
	}

	if ($overlapscore->[$minI][$minJ] < $CLUSTERCUTOFF) {
		$maxprob = 0;
		$minJ = -1;
		# find the highest prob model from 'james_cluster' still not placed
		foreach my $i (0..$#ARGV) {
			if ($james_clusterid[$i] != -1 && $aligned[$i] == 0) {
				my $tag = $ARGV[$i]; $tag =~ s/.*_(\d\d\d).*/$1/;
				my $prob = $template_probs{ $tag };
				$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
				if ($prob > $maxprob) {
					$maxprob = $prob;
					$minJ = $i;
				}
			}
		}

		print STDERR "Adding new cluster $minJ with prob $maxprob\n";
		last if ($minJ == -1); # we're done!

		# find the closest aligned thing to 'minJ'
		$minRMS = 999; $minI = -1;
		$maxOverlap = 0;
		foreach my $i (0..$#ARGV) {
			next if ($i==$minJ || $aligned[$i] == 0);
			# i is aligned, j is not
			#if ( $rmsds->[$i][$minJ] < $minRMS ) {
			#	$minRMS = $rmsds->[$i][$minJ];
			#	$minI = $i;
			#}
			if ( $overlapscore->[$i][$minJ] > $maxOverlap ) {
				$maxOverlap = $overlapscore->[$i][$minJ];
				$minI = $i;
			}
		}
	}

	# add j to alignment
	print STDERR "ALIGN ".$ARGV[$minJ]." (from".$ARGV[$minI].")\n";
	$globalRs->[$minJ] = mmult( $globalRs->[$minI], $Rs->[$minI][$minJ]  );
	$preTs->[$minJ] = deep_copy( $comis->[$minI][$minJ] );

	# james cluster id
	# since we only seed with things that have a james_cluster id, we will never have minI => -1 and minJ => -1
	if ($james_clusterid[$minJ] == -1 && $james_clusterid[$minI] >= 0) {
		$james_clusterid[$minJ] = $james_clusterid[$minI];
	} elsif ($james_clusterid[$minI] == -1 && $james_clusterid[$minJ] >= 0) {
		$james_clusterid[$minI] = $james_clusterid[$minJ];
	}

	# check to see if this belongs in a new cluster
	if ($overlapscore->[$minI][$minJ] >= $CLUSTERCUTOFF) {
		$clusterid->[$minJ] = $clusterid->[$minI];
	} else {
		$clusterid->[$minJ] = $nextclusterid;
		$nextclusterid++;
	}

	# post-rotation translation is a bit tricky
	# we need to transform to j's coordinate frame,
	#   then apply j's transformation to this result
	my $post_to_j = vadd( $comjs->[$minI][$minJ], $preTs->[$minJ] );
	my $post_j_to_global = vadd( mapply( $globalRs->[$minI] , vsub( $post_to_j , $preTs->[$minI]) ) , 
	                             $postTs->[$minI] );
	$postTs->[$minJ] = $post_j_to_global;

	$aligned[$minJ] = 1;
}

foreach my $i (0..$#ARGV) {
	if ($aligned[$i] == 0) {
		print STDERR "ELIMINATED ".$ARGV[$i]."\n";
	}
}

# apply transformation; dump aligned templates
mkdir "aligned_templates";
foreach my $i (0..$#ARGV) {
	next if ($aligned[$i] == 0);

	my $pdb = $ARGV[$i];
	my $nfrags = scalar( @{ $allfrags{$pdb} } );

	my $outpdb = $pdb;
	$outpdb =~ s/.*\///;
	$outpdb =~ s/\.pdb$/_aln.$i.pdb/;
	#$outpdb = "cl".$clusterid->[$i].".".$outpdb;
	print STDERR "writing aligned_templates/$outpdb\n";
	open (PDBF, ">aligned_templates/$outpdb") || print STDERR "Cannot open $_";

	foreach my $frag (0..$nfrags-1) {
		my $outfile = $pdb;
		next if scalar( @{ $allfrags{$pdb}->[$frag] } <=5 );
		foreach my $line (@{ $allfrags{$pdb}->[$frag] }) {
			my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
			$newX = vsub( $newX, $preTs->[$i]);
			$newX = mapply( $globalRs->[$i], $newX );
			$newX = vadd( $newX, $postTs->[$i]);

			substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
			print PDBF $line;
		}
		close PDB;
	}
	close PDBF;
}


# make full-length input file
# use seed
my $pdb = $ARGV[$seedPDB];
my @pseudotrace;
my @gaps;
my $last_ungapped = 1;
foreach my $i (1..$nres) {
	my $id = " CA "."A".$i;
	
	# add it to pseudotrace
	if (defined $bbatoms{$pdb}->{$id} && $fraglens{$pdb}->{ $id } > 5 ) {
		$gaps[$i] = 0;
		$last_ungapped = $i;
		foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
			my $x;
			$x = $bbatoms{$pdb}->{$atm."A".$i};
			$x = vsub( $x, $preTs->[$seedPDB]);
			$x = mapply( $globalRs->[$seedPDB], $x );
			$x = vadd( $x, $postTs->[$seedPDB]);
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
			my $dir=-1; # extend backwards
			if ( ($i != $nres && $gaps[$i+1] < $cycle)) {
				$dir=1; #extend forwards
			}
			print STDERR "Rebuild $i ($dir)...\n";

			my $Xp1=[];
			my $Xp2=[];
			foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
				push @{$Xp1}, deep_copy($pseudotrace[$i+1*$dir]->{$atm});
				push @{$Xp2}, deep_copy($pseudotrace[$i+2*$dir]->{$atm});
			}

			# align xp2->xp1
			my ($R21, $rmsd21, $com2, $com1) = rms_align( $Xp1 , $Xp2 );

			# apply to xp1
			foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
				my $newX = vsub( $pseudotrace[$i+1*$dir]->{$atm}, $com2);
				$newX = mapply( $R21, $newX );
				$newX = vadd( $newX, $com2);
				$newX = vadd( $newX, $com1);
				$pseudotrace[$i]->{$atm} = $newX;
			}
		}
	}
}
		
# write model
print STDERR "writing input.pdb\n";
open (PDB, ">input.pdb") || print STDERR "Cannot open $_";
my $atmidx = 1;
foreach my $i (1..$nres) {
	my $restype = $one_to_three{ substr($seq, $i-1, 1) };
	foreach my $atm (" N  ", " CA ", " C  ", " O  ") {
		my $x = $pseudotrace[$i]->{$atm};
		printf PDB "ATOM   %4d %s %s %s %3d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
			   $atmidx++, $atm, $restype, 'A', $i, $x->[0], $x->[1], $x->[2];
	}
}
close(PDB);



# write config file
#open (HYB, ">hybrid.config") || print STDERR "Cannot open $_";
foreach my $i (0..$#ARGV) {
	next if ($aligned[$i] == 0);
	my $dir = getcwd;
	my $id = $ARGV[$i];
	$id =~ s/.*\///;
	$id =~ s/\.pdb$/_aln.$i.pdb/;
	$id = "$dir/aligned_templates/$id";
	my $cstfile = sprintf $CSTFILENAMES, $james_clusterid[$i];
	if ($cstfile ne "NULL") {
		$cstfile = $dir."/$cstfile";
	}
	my $clusterid = $clusterid->[$i];
	my $tag = $ARGV[$i]; $tag =~ s/.*_(\d\d\d).*/$1/;
	my $prob = $template_probs{ $tag };
	$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };

	my $outline = sprintf "%20s %30s %4d %.5f\n", $id, $cstfile, $clusterid, $prob;
	#print HYB $outline;
	print $outline;
}

exit 0;

###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

