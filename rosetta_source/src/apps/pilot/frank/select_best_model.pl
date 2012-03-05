#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);

use lib (".");
use File::Basename;
use File::Path qw(rmtree);
use Cwd;
use Storable qw(dclone);
use lib dirname(__FILE__);
require "kabsch.pm";
require "matrix.pm";

## programs
my $EXTRACTAPP = "/work/dimaio/rosetta/rosetta_source/bin/extract_pdbs.default.linuxgccrelease".
                 " -database /work/dimaio/rosetta/rosetta_database ".
				 "-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer  VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm";

#####
#####

## structure alignments
my $NMODELSOUT=5;

my @RMS_CUTOFFS = (8,5,2.5);
my $CLUSTERCUTOFF = 0.0;
my $BOLTZMANTEMP = 20;

## remove neighbors
my $REMOVENEIGHBORS = 10;

## main()
if ($#ARGV < 1) {
	print STDERR "usage: $0 <topN> <silent_file>\n";
	exit -1;
}

my $N = shift @ARGV;
my $silentfile = shift @ARGV;

## STAGE 0 -- make output directory
my $WORKDIR = $silentfile;
$WORKDIR =~ s/\.out//;
$WORKDIR =~ s/\.silent//;
rmtree($WORKDIR);
mkdir $WORKDIR;

# scorecut
chdir $WORKDIR;
open( INSILENT , "../$silentfile" ) || die "Cannot open $silentfile\n";
my @filelines = <INSILENT>;
chomp (@filelines);

my @scores;
my %scoremap;
my %gdtmap;
my $counter = 0;
foreach my $line (@filelines) {
	if ($line =~ /^SCORE:/ && $line !~ /description/) {
		my @fields = split ' ', $line;
		my $score = $fields[1] - $fields[16] - $fields[22];
		$scores[$counter] = $score;
		$scoremap{ $fields[$#fields] } = $score;
		$gdtmap{ $fields[$#fields] } = $fields[25];
		$counter++;
	}
}
my $num_read = $counter;
$counter = 0;

## get cutoff
my $cutat = $N-1;
if ($cutat > $#scores) { $cutat = $#scores; }
my @scores_sort = sort { $a <=> $b } @scores;
my $scorecut = $scores_sort[$cutat];

## make pruned silent file
my $outfile = $silentfile;
$outfile =~ s/\.[\w]*$/.top$N\.silent/;
`head -2 ../$silentfile > $outfile`;
open( OUTSILENT , ">>$outfile" );

my $printstate = 0;
$counter = 0;
my $accept = 0;
foreach my $line (@filelines) {
	if ($line =~ /^SCORE:/ && $line !~ /description/) {
		if ($scores[$counter] <= $scorecut) { $printstate = 1; $accept++; } else { $printstate = 0; }
		$counter++;
	}
	if ($printstate == 1) {
		print OUTSILENT $line."\n";
	}
}

# weights
my %weightmap;
foreach my $tag (keys %scoremap) {
	if ($BOLTZMANTEMP < 0 ) { 
		$weightmap{ $tag } = 1;
	} else {
		$weightmap{ $tag } = exp( ( $scorecut - $scoremap{$tag} ) / $BOLTZMANTEMP );
	}
}

# extract
my $cmd = "$EXTRACTAPP -in:file:silent $outfile -in:file:silent_struct_type binary &> /dev/null";
#print STDERR $cmd."\n";
system($cmd);

## STAGE 2 -- clustering + alignment
##
my @THREADED_MDLS=<*.pdb>;
my %allfrags;
my %bbatoms;
my %caatoms;
my %resids;
my %fraglens;
my %ss;

# read in all PDBs
my $minres = 9999;
my $maxres = 0;
foreach my $pdb (@THREADED_MDLS) {
	open (PDB, $pdb) || die "Cannot open $pdb";

	$bbatoms{$pdb} = {};
	$caatoms{$pdb} = {};
	$resids{$pdb} = {};
	$allfrags{$pdb} = [];

	my $curr_frag_start = -999;
	my $last_res_read = -999;
	while (my $line = <PDB>) {
		next if ($line !~ /^ATOM/);

		my $atom = substr ($line, 12, 4);
		my $chain = substr($line, 21, 1);
		my $resid = substr($line, 17, 3);
		my $res = int( substr($line, 22, 4) );
		$minres = min($minres,$res);
		$maxres = max($maxres,$res);

		my $x = substr ($line, 30, 8);
		my $y = substr ($line, 38, 8);
		my $z = substr ($line, 46, 8);

		my $id = $atom.$chain.$res;
		$bbatoms{$pdb}->{ $id } = [$x,$y,$z];
		if ($atom eq " CA ") {
			$caatoms{$pdb}->{ $id } = [$x,$y,$z];
		}
		$resids{$pdb}->{ $id } = $resid;
		push @{$allfrags{$pdb}}, $line;
	}

	close (PDB);
}

# align all->all
my $npdbs = @THREADED_MDLS;
my $rmsds = []; my $nalignlen = []; my $Rs = [];
my $comis = []; my $comjs = [];
my $overlapscore = []; # for clustering

foreach my $i (0..$#THREADED_MDLS) {
foreach my $j ($i..$#THREADED_MDLS) {
	# intersection of $bbatoms{i} and {j}
	my @common_atoms = ();
	foreach ( keys %{ $caatoms{$THREADED_MDLS[$i]} } ) {
		push @common_atoms, $_ if exists $caatoms{$THREADED_MDLS[$j]}->{$_};
	}
	my $ncommon_atoms = scalar( @common_atoms );

	if ($ncommon_atoms < 4) {
		$rmsds->[$i][$j] = 9999;
		$Rs->[$i][$j] = [[1,0,0],[0,1,0],[0,0,1]];
		$nalignlen->[$i][$j] = scalar(@common_atoms);
	 	$overlapscore->[$i][$j] = 0;
		next;
	}

	# atom lists
	my $atoms_i = [];
	my $atoms_j = [];
	foreach ( @common_atoms ) {
		push @{ $atoms_i }, deep_copy( $caatoms{$THREADED_MDLS[$i]}->{$_} );
		push @{ $atoms_j }, deep_copy( $caatoms{$THREADED_MDLS[$j]}->{$_} );
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
			my $x = $caatoms{$THREADED_MDLS[$i]}->{$_};
			my $y = vadd( mapply( $Rs->[$i][$j] , vsub( $caatoms{$THREADED_MDLS[$j]}->{$_} , $comis->[$i][$j] ) ), 
			              vadd($comis->[$i][$j] , $comjs->[$i][$j]) );
			my $dist = dist( $x,$y );
			if ($dist <= $RMS_ALIGN) {
				push @new_common_atoms, $_;
				push @{ $atoms_i }, deep_copy( $caatoms{$THREADED_MDLS[$i]}->{$_} );
				push @{ $atoms_j }, deep_copy( $caatoms{$THREADED_MDLS[$j]}->{$_} );
			}
		}

		$ncommon_atoms = scalar( @new_common_atoms );

		my $natoms_tot = min( scalar(keys %{ $caatoms{$THREADED_MDLS[$i]} }),  scalar(keys %{ $caatoms{$THREADED_MDLS[$j]} }) );

		if ($RMS_ALIGN == $RMS_CUTOFFS[$#RMS_CUTOFFS]) {
			$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $natoms_tot;
			last;
		}

		if ($ncommon_atoms < 4) {
		 	$overlapscore->[$i][$j] = 0;
			next;
		}
		@common_atoms = @new_common_atoms;

		($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
		$nalignlen->[$i][$j] = scalar(@common_atoms);
	}
}
}

foreach my $i (0..$#THREADED_MDLS) {
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

# start with highest probability member as the seed
my $maxprob = 1;
my $minI = 0;

my $nextclusterid = 1;
$clusterid->[$minI] = $nextclusterid;
$nextclusterid++;

# use i's coordinate frame
my $globalRs = [];
$globalRs->[$minI] = [[1,0,0],[0,1,0],[0,0,1]];
my $preTs = [];
$preTs->[$minI] = deep_copy( $comis->[$minI][$minI] );
my $postTs = [];
$postTs->[$minI] = [0,0,0];  # center at origin

# for remaining structs find best alignment to already-aligned model
my @aligned = (0) x scalar(@THREADED_MDLS);
$aligned[$minI] = 1;

my $seedPDB = $minI;

my $minJ;
my $minRMS;
my $maxOverlap;
foreach my $cycle (1..$#THREADED_MDLS) {
	$minRMS = 999;
	$maxOverlap = 0;
	$minI = -1; $minJ = -1;
	foreach my $i (0..$#THREADED_MDLS) {
	foreach my $j (0..$#THREADED_MDLS) {
		next if ($i==$j);
		next if ($aligned[$i] == 0 || $aligned[$j] == 1);
		# i is aligned, j is not
		if ( $overlapscore->[$i][$j] > $maxOverlap ) {
			$maxOverlap = $overlapscore->[$i][$j];
			$minI = $i; $minJ = $j;
		}
	}
	}

	# add j to alignment
	$globalRs->[$minJ] = mmult( $globalRs->[$minI], $Rs->[$minI][$minJ]  );
	$preTs->[$minJ] = deep_copy( $comis->[$minI][$minJ] );

	# check to see if this belongs in a new cluster
	if ($overlapscore->[$minI][$minJ] >= $CLUSTERCUTOFF) {
		$clusterid->[$minJ] = $clusterid->[$minI];
	} else {
		$clusterid->[$minJ] = $nextclusterid;
		$nextclusterid++;
	}

	my $post_to_j = vadd( $comjs->[$minI][$minJ], $preTs->[$minJ] );
	my $post_j_to_global = vadd( mapply( $globalRs->[$minI] , vsub( $post_to_j , $preTs->[$minI]) ) , 
	                             $postTs->[$minI] );
	$postTs->[$minJ] = $post_j_to_global;

	$aligned[$minJ] = 1;
}

# sort clusters by p(correct)
my %cluster_members;
my @clusterprobs;
foreach my $i (1..$nextclusterid-1) {
	my @currclust = ();
	$clusterprobs[$i-1] = 0;
	foreach my $j (0..$#THREADED_MDLS) {
		if ($clusterid->[$j] == $i) {
			push @currclust, $j;

			my $pdbidstem = $THREADED_MDLS[$j];
			$pdbidstem =~ s/\.pdb//g;
			my $prob = $weightmap{$pdbidstem};
			$clusterprobs[$i-1] += $prob;
		}
	}
	$cluster_members{ $clusterprobs[$i-1] } = \@currclust;
}


@clusterprobs = sort {$b <=> $a} @clusterprobs;
$counter = 0;
my $nclust = $#clusterprobs;

# store new cluster assignments
foreach my $j (0..$#THREADED_MDLS) {
	$clusterid->[$j] = -1;
}
foreach my $i (0..$nclust) {
	my $prob = $clusterprobs[$i];
	foreach my $member (@{ $cluster_members{ $prob } }) {
		$clusterid->[ $member ] = $i;
	}
}

## get low energy & median of each cluster
my @towrite;
foreach my $i (0..$#THREADED_MDLS) {
	$towrite[$i] = 0;
}

my @chooseGDT = (0) x $NMODELSOUT;
my @chooseModel = (-1) x $NMODELSOUT;
my $modelswritten = 0;

while ($modelswritten < $NMODELSOUT) {
	foreach my $clid (0,0,1,2,3,4) {
		last if ($modelswritten == $NMODELSOUT);
	
		# once we sample every cluster go back to cluster 1
		my $effclid = $clid;
		if ($effclid >= $nclust) {
			$effclid = 0;
		}
	
		$chooseGDT[$modelswritten] = -1;
		my $lowmaxsub = -999999;
		outer: foreach my $i (0..$#THREADED_MDLS) {
			next if ($towrite[$i] != 0);
			next if ($clusterid->[$i] != $effclid);
	
			# remove identical
			#foreach my $j (0..$modelswritten-1) {
			#	next outer if ($overlapscore->[$i][$chooseModel[$j]] == 1);
			#}
	
			my $pdbidstem = $THREADED_MDLS[$i];
			$pdbidstem =~ s/\.pdb//g;
			my $score = $weightmap{$pdbidstem};
			foreach my $j (0..$#THREADED_MDLS) {
				next if ($clusterid->[$j] == -1);  # use the other clusters to decide?
				next if ($i==$j);

				my $pdbidstemJ = $THREADED_MDLS[$j];
				$pdbidstemJ =~ s/\.pdb//g;
				my $prob = $weightmap{$pdbidstemJ};
				$score += $prob*$overlapscore->[$i][$j];
			}
	
			if ($score > $lowmaxsub) {
				$chooseModel[$modelswritten] = $i;
				$chooseGDT[$modelswritten] = $gdtmap { $pdbidstem };
				$lowmaxsub = $score;
			}
		}
		if ($chooseModel[$modelswritten] != -1) {
			$towrite[$chooseModel[$modelswritten]] = $modelswritten+1;

			# remove this model and its neighbors
			$clusterid->[$chooseModel[$modelswritten]] = -1;

			removeloop: foreach my $k (1..$REMOVENEIGHBORS) {
				$clusterid->[$chooseModel[$modelswritten]] = -1;
				my ($maxScore,$maxJ) = (0,-1);
				foreach my $j (0..$#THREADED_MDLS) {
					next if ($clusterid->[$j] != $effclid);
					next if ($chooseModel[$modelswritten]==$j);
					my $score = $overlapscore->[$chooseModel[$modelswritten]][$j];
					if ($score > $maxScore) { $maxScore = $score; $maxJ = $j; }
				}
				if ($maxJ != -1) { $clusterid->[$maxJ] = -1; }
			}
			$modelswritten++;
		}
	}
}


# model 1 - median of largest cluster
# my $lowmaxsub = -999999;
# foreach my $i (0..$#THREADED_MDLS) {
# 	next if ($clusterid->[$i] != 0);
# 
# 	my $pdbidstem = $THREADED_MDLS[$i];
# 	$pdbidstem =~ s/\.pdb//g;
# 	my $score = $weightmap{$pdbidstem};
# 
# 	foreach my $j (0..$#THREADED_MDLS) {
# 		next if ($clusterid->[$j] == -1);
# 		next if ($i==$j);
# 
# 		my $pdbidstemJ = $THREADED_MDLS[$j];
# 		$pdbidstemJ =~ s/\.pdb//g;
# 		my $prob = $weightmap{$pdbidstemJ};
# 
# 		$score += $prob*$overlapscore->[$i][$j];
# 	}
# 
# 	if ($score > $lowmaxsub) {
# 		$chooseModel[0] = $i;
# 		$chooseGDT[0] = $gdtmap { $pdbidstem };
# 		$lowmaxsub = $score;
# 	}
# }
# $towrite[$chooseModel[0]] = 1;
# $modelswritten++;

# models 2-5 - low energy cl 1-4
# foreach my $clid (0..4) {
# 	last if ($modelswritten == 5);
# 
# 	# once we sample every cluster go back to cluster 1
# 	my $effclid = $clid;
# 	if ($effclid >= $nclust) {
# 		$effclid = 0;
# 	}
# 
# 	$chooseGDT[$modelswritten] = -1;
# 	my $lowenergy = 999999;
# 	outer: foreach my $i (0..$#THREADED_MDLS) {
# 		next if ($towrite[$i] != 0);
# 		next if ($clusterid->[$i] != $effclid);
# 
# 		# remove identical
# 		foreach my $j (0..$modelswritten-1) {
# 			#print $i." ".$chooseModel[$j]." ".$overlapscore->[$i][$chooseModel[$j]]."\n";
# 			next outer if ($overlapscore->[$i][$chooseModel[$j]] == 1);
# 		}
# 
# 		my $pdbidstem = $THREADED_MDLS[$i];
# 		$pdbidstem =~ s/\.pdb//g;
# 		my $score = $scoremap { $pdbidstem };
# 		if ($score < $lowenergy) {
# 			$chooseModel[$modelswritten] = $i;
# 			$chooseGDT[$modelswritten] = $gdtmap { $pdbidstem };
# 			$lowenergy = $score;
# 		}
# 	}
# 	if ($chooseModel[$modelswritten] != -1) {
# 		$towrite[$chooseModel[$modelswritten++]] = $modelswritten;
# 	}
# }

print "$WORKDIR ".(join ' ',@chooseGDT)."\n";

## STAGE 3 -- output
foreach my $i (0..$#THREADED_MDLS) {
	unlink ($THREADED_MDLS[$i]); # clean up

	next if ($towrite[$i] == 0);
	my $pdb = $THREADED_MDLS[$i];

	my $outpdb = "$WORKDIR.model".$towrite[$i].".pdb";
	#print STDERR "writing $outpdb\n";
	open (PDBF, ">$outpdb") || die "Cannot open $_";

	foreach my $line (@{ $allfrags{$pdb} }) {
		my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
		$newX = vsub( $newX, $preTs->[$i]);
		$newX = mapply( $globalRs->[$i], $newX );
		$newX = vadd( $newX, $postTs->[$i]);
		substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
		substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
		substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
		print PDBF $line;
	}
	close PDBF;
}



exit 0;

###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}
