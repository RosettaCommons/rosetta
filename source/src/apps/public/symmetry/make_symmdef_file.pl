#!/usr/bin/perl
##
##
###############################################################################

use strict;
use Math::Trig;	  # inv trig ops
use POSIX qw(ceil floor fmod fabs);
use List::Util qw(max min);
use File::Basename;
#use Getopt::Long qw(permute);
use constant PI	   => 4 * atan2(1, 1);

#use lib (".");
#use lib dirname(__FILE__);

###############################################################################

if ($#ARGV < 0) {
	print STDERR "usage: $0 [options]\n";
	print STDERR "\n";
	print STDERR "example:	 $0 -m NCS -a A -i B C -r 12.0 -p mystructure.pdb\n";
	print STDERR "example:	 $0 -m HELIX -a A -b B -i G -r 12.0 -p mystructure.pdb\n";
	print STDERR "example:	 $0 -m CRYST -r 12.0 -c 42.4 41.2 88.6 90.0 90.0 90.0 -s P 1 21 1 -p mystructure.pdb\n";
	print STDERR "\n";
	print STDERR "common options: \n";
	print STDERR "	  -m (NCS|CRYST|HELIX|PSEUDO) : [default NCS] which symmetric mode to run\n";
	print STDERR "			  NCS: generate noncrystallographic (point) symmetries from multiple chains in a PDB file\n";
	print STDERR "			  CRYST: generate crystallographic symmetry (fixed unit cell) from the CRYST1 line in a PDB file\n";
	print STDERR "			  HELIX: generate helical/fiber symmetry from multiple chains in a PDB file\n";
	print STDERR "			  PSEUDO: (EXPERIMENTAL) generate pseudo-symmetric system\n";
	print STDERR "	  -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
	print STDERR "	  -r <real>	  : [default 1e6 for NCS+helix, 12.0 for cryst] the max CA-CA distance between two interacting chains\n";
	print STDERR "	  -q		  : [default false] quiet mode (no files are output)\n";
	print STDERR "\n";
	print STDERR "NCS-specific options: \n";
	print STDERR "	  -a <char>	  : [default A] the chain ID of the main chain\n";
	print STDERR "	  -d <char>*  : the chain ID of other chains to keep in the output file\n";
	print STDERR "	  -i <char>*  : [default B] the chain IDs of one chain in each symmetric subcomplex\n";
	print STDERR "	  -e		  : [default false] allow rigid body minimization of complete system\n";
	print STDERR "\n";
	print STDERR "CRYST-specific options: \n";
	print STDERR "	  -c <real>x6 : override the unit cell parameters in the PDB with these values\n";
	print STDERR "	  -g <real>x3 : perturb the unit cell parameters (A,B,C only) in the PDB with these values\n";
	print STDERR "	  -s <string> : override the spacegroup in the PDB with these values\n";
	print STDERR "	  -h		  : [default false] dont restrict translational motion along axes that do not change the system\n";
	print STDERR "\n";
	print STDERR "HELIX-specific options: \n";
	print STDERR "	  -a <char>	  : [default A] the chain ID of the main chain\n";
	print STDERR "	  -b <char>[*]: [default B] the chain ID of the next chain along the fiber/helix\n";
	print STDERR "				  : Note that three optional arguments may be provided, chain:symm_order:rise:axis\n";
	print STDERR "				  : This will force a symm order (#subunits per rotation), a rise per subunit,\n";
	print STDERR "					and a helical axis ('x','y', or 'z'), respectively.\n";
	print STDERR "	  -i <char>	  : the chain ID of a chain in -a's point symmetry group\n";
	print STDERR "	  -t <real>	  : [default 4] the number of subunits to generate along the -b direction\n";
	print STDERR "	  -e		  : [default false] allow rigid body minimization of complete system\n";
	print STDERR "\n";
	print STDERR "PSEUDO-specific options: \n";
	print STDERR "	  -a <char>	  : [default A] the chain ID of the main chain\n";
	exit -1;
}


##	set default options
##
my $pdbfile = '';
my $interact_dist = 0.0;
my $primary_chain = 'A';
my %keep_chains = ();
my @secondary_chains = ();
my $helical_chain = 'B';
my $rbminAll = 0;
my @cell_new;
my @cell_offset;
my $spacegp_new;
my $modestring = "NCS";
my $nturns = 4;
my $quietMode = 0;
my $restrictCrystTrans = 0;


## parse options (do this by hand since Getopt does not handle this well)
##
my $inlinefull = (join ' ',@ARGV)." ";
my @suboptions = split /(-[a-z|A-Z] )/, $inlinefull;
for ( my $i=0; $i<=$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-m " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$modestring = $suboptions[++$i];
		$modestring =~ s/\s*(\S+)\s*/$1/;
	}
	elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
		$pdbfile =~ s/\s*(\S+)\s*/$1/;
	}
	elsif ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$interact_dist = $suboptions[++$i];
	}
	elsif ($suboptions[$i] =~ /^-q/ ) {
		$quietMode = 1;
	}
	elsif ($suboptions[$i] =~ /^-e/ ) {
		$rbminAll = 1;
	}
	elsif ($suboptions[$i] =~ /^-h/ ) {
		$restrictCrystTrans = 1;
	}
	# cryst
	elsif ($suboptions[$i] eq "-c " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@cell_new = split /[, ]/,$suboptions[++$i];
	}
	elsif ($suboptions[$i] eq "-g " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@cell_offset = split /[, ]/,$suboptions[++$i];
	}
	elsif ($suboptions[$i] eq "-s " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$spacegp_new = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-t " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$nturns = ( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-a " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$primary_chain = $suboptions[++$i];
		$primary_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-b " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$helical_chain = $suboptions[++$i];
		$helical_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-i " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@secondary_chains = split /[, ]/,$suboptions[++$i];
	} elsif ($suboptions[$i] eq "-d " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		my @keep_chain_list = split /[, ]/,$suboptions[++$i];
		foreach my $keepchain ( @keep_chain_list ) {
			$keep_chains{ $keepchain } = 1;
		}
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

##
if ($pdbfile eq '') {
	print STDERR "Must provide an input PDB file with -p\n";
	exit -1;
}

## parse mode string
##	  + set mode-specific defaults
##
my ($cryst_mode, $ncs_mode, $helix_mode, $pseudo_mode) = (0,0,0,0);
if ($modestring eq "CRYST" || $modestring eq "cryst" || $modestring eq "Cryst") {
	$cryst_mode = 1;
} elsif ($modestring eq "NCS" || $modestring eq "ncs" || $modestring eq "Ncs") {
	$ncs_mode = 1;
	if ( scalar(@secondary_chains) == 0) {
		@secondary_chains = ('B');
	}
} elsif ($modestring eq "HELIX" || $modestring eq "helix" || $modestring eq "Helix") {
	$helix_mode = 1;
} elsif ($modestring eq "PSEUDO" || $modestring eq "pseudo" || $modestring eq "Pseudo") {
	$pseudo_mode = 1;
} else {
	print STDERR "Unrecognized mode string '$modestring'\n";
	exit -1;
}
if ($quietMode!= 1) { print STDERR "Running in mode $modestring.\n"; }

## set default interaction radius (mode-specific)
if ($interact_dist == 0) {
	if ($cryst_mode) {
		$interact_dist = 12.0;
	} else {
		$interact_dist = 1e6;
	}
}

## substitute'_' -> ' '
if ($primary_chain eq '_') {
	$primary_chain = ' ';
}
if ($helical_chain eq '_') {
	$primary_chain = ' ';
}
foreach my $i (0..$#secondary_chains) {
	if ($secondary_chains[$i] eq '_') {
		$secondary_chains[$i] = ' ';
	}
}


###
### Read input PDB file
my %chains;		 # input as separate chains
my @chaintrace;	 # input as a single monomer
my @filebuf;
my @altfilebuf;
my $minRes = 1;
my $monomerRadius = 0;

# crystinfo
my $spacegp = "P 1";
my ($gpid, $nsymm, $Rs, $Ts, $Cs, $cheshire);
my ($A, $B, $C, $alpha, $beta, $gamma) = (0,0,0,90,90,90);
my ($f2c,$c2f);

my $CoM = [0,0,0];
my $anchor_ca = [0,0,0];
my $CoM;

if ($cryst_mode == 1) {
	if ( length ($pdbfile) > 0 ) {
		###
		### Input PDB file
		if ($pdbfile =~ "\.gz\$") {
			open (PDB, "gunzip -c $pdbfile | ") || die "Cannot open $pdbfile.";
		} else {
			open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
		}
		while (<PDB>) {
			chomp;
			if (/^CRYST1/) {
				$A = substr ($_,  6, 9);
				$B = substr ($_, 15, 9);
				$C = substr ($_, 24, 9);
				$alpha = substr ($_, 33, 7);
				$beta  = substr ($_, 40, 7);
				$gamma = substr ($_, 47, 7);
				$spacegp =	substr ($_, 55, 11);
			} elsif (/^ATOM/ || /^HETATM/) {
				my $atom = substr ($_, 12, 4);
				if ($atom eq " CA " || $atom eq " P  ") {
					push @chaintrace, [substr ($_, 30, 8),substr ($_, 38, 8),substr ($_, 46, 8)];
					$CoM = vadd( $chaintrace[ $#chaintrace ], $CoM );
				}

				push @filebuf, $_;
			}
		}
		close (PDB);

		# override crystal params from cmd line
		if ($#cell_new >= 5) {
			($A,$B,$C,$alpha,$beta,$gamma) = @cell_new;
		}
		if ($#cell_offset >= 2) {
			$A += $cell_offset[0];
			$B += $cell_offset[1];
			$C += $cell_offset[2];
		}
		if (length($spacegp_new) > 0) {
			$spacegp = $spacegp_new;
		}

		# initialize symmops
		($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp_lookup( $spacegp );
		($f2c,$c2f) = crystparams($A,$B,$C,$alpha,$beta,$gamma);
	}

	# find residue closest to CoM of the system
	## find the radius of the molecule (for -f)
	my $maxDist2 = 0;
	my $minDist2 = 9999999;
	my $nAtms = $#chaintrace+1;
	foreach my $i (0..$nAtms) {
		my $dist2 = vnorm2( vsub( $CoM, $chaintrace[ $i ] ) );
		if ($dist2 < $minDist2) {
			$minDist2 = $dist2;
			$minRes = $i+1;
			$anchor_ca = deep_copy($chaintrace[ $i ]);
		}
		if ($dist2 > $maxDist2) {
			$maxDist2 = $dist2;
		}
	}
	$monomerRadius = sqrt( $maxDist2 );

	# get CoM
	$CoM = [ $CoM->[0]/$nAtms , $CoM->[1]/$nAtms , $CoM->[2]/$nAtms ];
} else {
    # not cryst mode
	if ($pdbfile =~ "\.gz\$") {
		open (PDB, "gunzip -c $pdbfile | ") || die "Cannot open $pdbfile.";
	} else {
		open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
	}
	while (<PDB>) {
		chomp;
		if (/^ATOM/ || /^HETATM/) {
			my $chnid = substr ($_, 21, 1);
			my $atom  = substr ($_, 12, 4);
			if ($atom eq " CA ") {
				if (!defined $chains{ $chnid } ) {
					$chains{ $chnid } = [];
				}
				my $CA_i = [substr ($_, 30, 8),substr ($_, 38, 8),substr ($_, 46, 8)];
				push @{ $chains{ $chnid } }, $CA_i;
			}
			if ($primary_chain eq $chnid) {
				push @filebuf, $_;
			} elsif ( defined $keep_chains{ $chnid } ) {
				push @altfilebuf, $_;
			}
		}
	}
	close (PDB);


	## recenter the primary chain
	##
	if ( ! defined $chains{ $primary_chain } ) {
		die "Chain '$primary_chain' not in input!\n";
	}
	## find residue closest to CoM of the system
	## find the radius of the molecule (for -f)
	my $maxDist2 = 0;
	my $minDist2 = 9999999;
	foreach my $i ( 0..scalar( @{ $chains{ $primary_chain } })-1 ) {
		my $dist2 = vnorm2( $chains{ $primary_chain }->[$i] );
		if ($dist2 < $minDist2) {
			$minDist2 = $dist2;
			$minRes = $i+1;
			$anchor_ca = deep_copy( $chains{ $primary_chain }->[$i] );
		}
		if ($dist2 > $maxDist2) {
			$maxDist2 = $dist2;
		}
	}
	$monomerRadius = sqrt( $maxDist2 );
}


####################################################################################
####################################################################################
###	  mode-specific stuff
####################################################################################
####################################################################################

if ($ncs_mode == 1) {
	my $NCS_ops = {};
	my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
	my $COM_0 = recenter( $chains{ $primary_chain } );
	$NCS_ops->{R} = $R_0;
	$NCS_ops->{T} = $COM_0;
	$NCS_ops->{PATH} = "";
	$NCS_ops->{CHILDREN} = [];
	$NCS_ops->{AXIS} = [0,0,1];


	my @allQs;
	my @allCOMs;
	my @sym_orders;

	foreach my $sec_chain (@secondary_chains) {
		my @sec_chain_ids = split( ':', $sec_chain );

		if ( ! defined $chains{ $sec_chain_ids[0] } ) {
			die "Chain $sec_chain not in input!\n";
		}

		## CA check
		if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $sec_chain_ids[0] } } ) ) {
			print STDERR "ERROR! chains '$primary_chain' and '$sec_chain' have different residue counts! (".
						 scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $sec_chain_ids[0] } } ).")\n";
			die "Chain length mismatch!\n";
		}

		# get superposition
		my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $sec_chain_ids[0] } );
		if ($quietMode != 1) {
			print STDERR "Aligning $primary_chain and $sec_chain wth RMS=$rmsd.\n";
		}

		if ( is_identity( $R ) ) {
			print STDERR "Chains $primary_chain and $sec_chain related by transformation only! Aborting.\n";
			exit 1;
		}

		my $del_COM = vsub ($COM_i, $COM_0);
		push @allCOMs, $del_COM;

		my ($X,$Y,$Z,$W)=R2quat($R);

		my $Worig = $W;
		my $Wmult = 1;
		if ($W < 0) { $W = -$W; $Wmult = -1; }
		my $omega = acos($W);
		my $sym_order = int(PI/$omega + 0.5);

		# optionally allow input to 'force' a symmetric order
		if ($#sec_chain_ids > 0) {
			$sym_order = $sec_chain_ids[1];
		}
		push @sym_orders, $sym_order;

		my $rotaxis = [$X,$Y,$Z]; normalize( $rotaxis );
		if ($quietMode != 1) {
    		print STDERR "Found ".$sym_order."-fold (".(PI/$omega).") axis to ".$sec_chain_ids[0]." ";
	    	print STDERR ": ".$rotaxis->[0]." ".$rotaxis->[1]." ".$rotaxis->[2]."\n";
        }

		# now make perfectly symmetrical version of superposition
		my $newW = -$Wmult *cos( PI/$sym_order );
		my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
		my $newQ = [$X*$S , $Y*$S, $Z*$S, $newW];
		push @allQs, $newQ;
	}

	# enforce special conditions
	# 1) Dn symmetries
	if ($#allQs == 1 && ($sym_orders[ 0 ] == 2 || $sym_orders[ 1 ] == 2) ) {
		# the two transformations must be about perpendicular axes
		my $X = [ $allQs[0]->[0],  $allQs[0]->[1],	$allQs[0]->[2] ];
		my $Y = [ $allQs[1]->[0],  $allQs[1]->[1],	$allQs[1]->[2] ];

		normalize($X); normalize($Y);
		my $Xtgt = vsub( $X , vscale(dot($X,$Y),$Y) );
		my $Ytgt = vsub( $Y , vscale(dot($X,$Y),$X) );
		normalize( $Xtgt ); normalize( $Ytgt );

		my $X0 = [ ($X->[0]+$Xtgt->[0])/2 , ($X->[1]+$Xtgt->[1])/2 , ($X->[2]+$Xtgt->[2])/2 ];
		my $Y0 = [ ($Y->[0]+$Ytgt->[0])/2 , ($Y->[1]+$Ytgt->[1])/2 , ($Y->[2]+$Ytgt->[2])/2 ];
		my $W_x = $allQs[0]->[3];
		my $W_y = $allQs[1]->[3];
		my $S_x = sqrt ( (1-$W_x*$W_x)/vnorm2($X0) );
		my $S_y = sqrt ( (1-$W_y*$W_y)/vnorm2($Y0) );

		$allQs[0] = [ $X0->[0]*$S_x , $X0->[1]*$S_x, $X0->[2]*$S_x, $W_x];
		$allQs[1] = [ $Y0->[0]*$S_y , $Y0->[1]*$S_y, $Y0->[2]*$S_y, $W_y];

		# if we're Dn symmetry, always expand the Cn before the C2
		if ($sym_orders[ 0 ] == 2 && $sym_orders[ 1 ] != 2) {
			my $temp;
			$temp = $allQs[1]; $allQs[1] = $allQs[0]; $allQs[0] = $temp;
			$temp = $allCOMs[1]; $allCOMs[1] = $allCOMs[0]; $allCOMs[0] = $temp;
			$temp = $sym_orders[1]; $sym_orders[1] = $sym_orders[0]; $sym_orders[0] = $temp;
			$temp = $secondary_chains[1]; $secondary_chains[1] = $secondary_chains[0]; $secondary_chains[0] = $temp;
		}
	}

	# Now create the symmetry tree
	foreach my $i (0..$#allQs) {
		my $newQ = $allQs[ $i ];
		my $del_COM = $allCOMs[ $i ];
		my $sym_order = $sym_orders[ $i ];

		my $newR = quat2R( $newQ->[0], $newQ->[1], $newQ->[2], $newQ->[3] );

		my $com_complex = $NCS_ops->{T};
		my $com_secondary;
		if ($i == 0 ) {
			$com_secondary = vadd( $del_COM, $COM_0 );
		} else {
			$com_secondary = vadd( $del_COM, vadd( $COM_0, mapply( $newR , vsub ( $com_complex , $COM_0 ) ) ) );
		}
		my $newDelCOM = vsub ( $com_secondary , $com_complex );

		# symmetrize newT
		my $err_pos = [0,0,0];
		my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
		foreach my $j (1..$sym_order) {
			$err_pos = vadd( $err_pos, mapply( $R_i,$newDelCOM ) );
			$R_i = mmult($newR, $R_i);
		}

		# special case for Dn symmetry
		if ( $i == 1 && $#allQs == 1 ) {
			my $axis_proj_i =  [ $allQs[1-$i]->[0],	 $allQs[1-$i]->[1],	 $allQs[1-$i]->[2] ];
			# project $delCOM along this axis
			normalize( $axis_proj_i );
			my $del_COM_inplane	   = vsub( $newDelCOM , vscale(dot($newDelCOM,$axis_proj_i),$axis_proj_i) );
			$err_pos = vscale( $sym_orders[ $i ] , $del_COM_inplane );
		}
		print STDERR "translation error = ".vnorm( $err_pos )."\n";

		# get the center of the symmgp
		my $adj_newDelCOM = vsub( $newDelCOM , [ $err_pos->[0]/$sym_order , $err_pos->[1]/$sym_order , $err_pos->[2]/$sym_order ] );

		my $axis_i =  [ $newQ->[0], $newQ->[1], $newQ->[2] ];
		normalize($axis_i);

		expand_symmops_by_split( $NCS_ops, $newR, $adj_newDelCOM, $sym_order, $axis_i);
	}
	print STDERR "system center	 ".$NCS_ops->{T}->[0]." ".$NCS_ops->{T}->[1]." ".$NCS_ops->{T}->[2]."\n";

	my ($nnodes,$nleaves) = tree_size( $NCS_ops );
	if ($quietMode != 1) {
		print STDERR "Found a total of $nleaves monomers in the symmetric complex.\n";
		print STDERR "Placing $nnodes virtual residues.\n";
	}

	## dist checks
	my $symops = tree_traverse( $NCS_ops );
	my %symminterface = ();
	my $counter = 0;
	foreach my $X_i ( @{ $chains{ $primary_chain } } ) {
		foreach my $Y_i (  @{ $chains{ $primary_chain } } ) {
			foreach my $symop (@{ $symops }) {
				my $id = $symop->{PATH};
				next if (defined $symminterface{ $id });

				#	x_i = R_i * (x_0 - COM_0) + COM_i
				#	The rms function already ofsets x_0 by -COM_0
				my $rX_i = vadd( $X_i , $COM_0 );
				my $rY_j = vadd( mapply($symop->{R}, $Y_i) , $symop->{T} );
				my $delXY = vsub( $rY_j,$rX_i );
				my $dist2XY = vnorm2( $delXY );

				if ($dist2XY < $interact_dist*$interact_dist) {
					# we have a hit! tag NCS copy $j_symm as a non-symmetic interface
					print STDERR " Adding interface '".$symop->{PATH}."'\n";
					$symminterface{ $id } = $counter++;
				}
			}
		}
	}

	# find the equation for the energy of the complex
	# we need to identify duplicates
	my @syminterfaces = sort { $symminterface{$a} <=> $symminterface{$b} } keys %symminterface;
	my %energy_counter;
	foreach my $interface (@syminterfaces) {
		$energy_counter{ $interface } = 1;
	}
	delete ($energy_counter{ $syminterfaces[0] });
	# delete self-interface(???)
	delete ($energy_counter{ $syminterfaces[0] });

	##
	OUTER1: foreach my $i (1..$#syminterfaces) {
		next if (!defined  $energy_counter{ $syminterfaces[$i] } );

		my $subtree_i = get_subtree ( $NCS_ops, $syminterfaces[ $i ] );
		my $R_i = $subtree_i->{R};
		my $T_i = vsub( $subtree_i->{T} , $COM_0 );

		foreach my $j ($i+1..$#syminterfaces) {
			my $subtree_j = get_subtree ( $NCS_ops, $syminterfaces[ $j ] );
			my $R_j = $subtree_j->{R};
			my $T_j = vsub( $subtree_j->{T} , $COM_0 );

			# if transform i is the inverse of transform j we have our (i,j) pair
			if ( is_inverse( $R_i,$T_i, $R_j,$T_j ) ) {
				$energy_counter{ $syminterfaces[$i] } = 2;
				delete ($energy_counter{ $syminterfaces[$j] });
				next OUTER1;
			}
		}


	}

	# now, for each unique 1->n interface
	# sum up the number of equivalent i->j interfaces
	foreach my $i (1..scalar(@{ $symops })-1) {  # no need to test 0
		my $symop_i = $symops->[$i];
		my $R_i = $symop_i->{R};
		my $T_i = vsub( $symop_i->{T} , $COM_0 );

		foreach my $j ($i+1..scalar(@{ $symops })-1) {
			my $symop_j = $symops->[$j];
			my $R_j = $symop_j->{R};
			my $T_j = vsub( $symop_j->{T} , $COM_0 );

			foreach my $inter_n (keys %energy_counter) {
				my $subtree_n = get_subtree ( $NCS_ops, $inter_n );
				my $R_n = $subtree_n->{R};
				my $T_n = vsub( $subtree_n->{T} , $COM_0 );

				# is the transform (Rn,Tn) equivalent to the transform (Ri,Ti)->(Rj,Tj)
				if ( is_equivalent( $R_n,$T_n, $R_i,$T_i, $R_j,$T_j ) ) {
					$energy_counter{ $inter_n }++;
					next;
				}
				# how about (Rj,Tj)->(Ri,Ti)
				if ( is_equivalent( $R_n,$T_n, $R_j,$T_j, $R_i,$T_i ) ) {
					$energy_counter{ $inter_n }++;
					next;
				}
			}
		}
	}

	## symm file gen
	## write output symm file
	my $symmname = $pdbfile;
	$symmname =~ s/\.pdb$//;
	$symmname = $symmname."_".get_topology( $NCS_ops );
	print "symmetry_name $symmname\n";
	print "E = ".($nleaves)."*VRT".$syminterfaces[0]."_base";
	foreach my $complex (sort { $symminterface{$a} <=> $symminterface{$b} } keys %energy_counter) {
		print " + ".$energy_counter{$complex}."*(VRT".$syminterfaces[0]."_base".":VRT".$complex."_base".")";
	}
	print "\n";
	print "anchor_residue COM\n";

	# virtual_coordinates_start
	# xyz VRT1 -1,0,0 0,1,0 0,0,0
	# xyz VRT2 0,-1,0 -1,0,0 0,0,0
	# xyz VRT3 1,0,0 0,-1,0 0,0,0
	# xyz VRT4 0,1,0 1,0,0 0,0,0
	# virtual_coordinates_stop
	my $vrts_by_depth = tree_traverse_by_depth( $NCS_ops );
	my ($vrt_lines, $connect_lines, $dof_lines) = fold_tree_from_ncs( $NCS_ops , $vrts_by_depth, \%symminterface );
	print "virtual_coordinates_start\n";
	foreach my $vrt_line (@{ $vrt_lines }) {
		print $vrt_line."\n";
	}
	# connect_virtual JUMP1 VRT1 VRT2
	# connect_virtual JUMP2 VRT1 VRT3
	# connect_virtual JUMP3 VRT1 VRT4
	print "virtual_coordinates_stop\n";
	foreach my $connect_line (@{ $connect_lines }) {
		print $connect_line."\n";
	}
	if ($rbminAll == 1) {
		print "set_dof JUMP0 x y z angle_x angle_y angle_z\n";
	}
	foreach my $dof_line (@{ $dof_lines }) {
		print $dof_line."\n";
	}
	my $counter = 0;
	foreach my $vrt_level ( @{ $vrts_by_depth } ) {
		if ($counter > 1) {
			print "set_jump_group JUMPGROUP$counter";
			foreach my $tag (@{ $vrt_level }) {
				# jump to _first_ child in each group
				if ($tag =~ /_0$/) {
					print " JUMP$tag";
				}
			}
			print "\n";
		}
		$counter++;
	}
	print "set_jump_group JUMPGROUP$counter";
	foreach my $leaf ( @{ $symops } ) {
		my $tag = $leaf->{PATH};
		print " JUMP$tag"."_to_com";
	}
	print "\n";
	$counter++;
	print "set_jump_group JUMPGROUP$counter";
	foreach my $tag ( sort (keys %symminterface) ) {
		print " JUMP$tag"."_to_subunit";
	}
	print "\n";

	if ($quietMode == 0) {

		# output complete symm complex
		########################################
		## write output pdb
		my $outpdb = $pdbfile;
		my $outmon = $pdbfile;
		my $outmdl = $pdbfile;

		my $suffix = "_model_$primary_chain";
		foreach my $chn (@secondary_chains) {
			$suffix = $suffix.$chn;
		}
		$suffix =~ s/://g;

		if ($outpdb =~ /\.pdb$/) {
			$outpdb =~ s/\.pdb$/_symm.pdb/;
			$outmdl =~ s/\.pdb$/$suffix.pdb/;
			$outmon =~ s/\.pdb$/_INPUT.pdb/;
		} else {
			$outpdb = $outpdb."_symm.pdb";
			$outmdl = $outpdb."_model.pdb";
			$outmon = $outmon."_INPUT.pdb";
		}
		open (OUTPDB, ">$outpdb");
		open (OUTMON, ">$outmon");
		open (OUTMDL, ">$outmdl");

		my $chnidx = 0;
		my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()-=_+;:,.<>";
		foreach my $symop (@{ $symops }) {
			foreach my $line (@filebuf) {
				my $linecopy = $line;

				my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
				my $X_0 = vsub($X,$COM_0);
				my $rX = vadd( mapply($symop->{R}, $X_0) , $symop->{T} );

				substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
				substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
				substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
				substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);

				print OUTPDB $linecopy."\n";

				if (defined $symminterface{ $symop->{PATH} }) {
					print OUTMDL $linecopy."\n";
				}
			}
			print OUTPDB "TER	\n";
			print STDERR "Writing interface ".$symop->{PATH}." as chain ".substr ($chains, $chnidx, 1)."\n";
			$chnidx++;
		}

		######################################
		foreach my $line (@filebuf) {
			my $linecopy = $line;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $X->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $X->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $X->[2]);

			print OUTMON $linecopy."\n";
		}
		foreach my $line (@altfilebuf) {
			my $linecopy = $line;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $X->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $X->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $X->[2]);

			print OUTMON $linecopy."\n";
		}
	}
}


if ($cryst_mode == 1) {
	###############
	### CRYST mode!
	###############
	print STDERR "Read: $A x $B x $C x ($alpha , $beta , $gamma)   $spacegp\n";

	my $nAtms = $#chaintrace+1;

	# frac coords of chaintrace
	my @fchaintrace;
	foreach my $X_i (@chaintrace) {
		my $fX_i = mapply( $c2f , $X_i );
		push @fchaintrace , $fX_i;
	}

	# frac CoM
	my $fCoM = mapply( $c2f , $CoM );

	# do the symmetric expansion
	my $Ts_expand = [
		[-1,-1,-1],[-1,-1,0],[-1,-1,1],	 [-1,0,-1],[-1,0,0],[-1,0,1],  [-1,1,-1],[-1,1,0],[-1,1,1],
		[0,-1,-1] ,[0,-1,0] ,[0,-1,1] ,	 [0,0,-1] ,		   ,[0,0,1] ,  [0,1,-1] ,[0,1,0] ,[0,1,1],
		[1,-1,-1] ,[1,-1,0] ,[1,-1,1] ,	 [1,0,-1] ,[1,0,0] ,[1,0,1] ,  [1,1,-1] ,[1,1,0] ,[1,1,1]
	];

	# nmr check
	my %symminterface = ();
	my $nmr = 0;
	if ($A == 1 && $B == 1 && $C == 1 && $alpha == 90 && $beta == 90 && $gamma == 90) {
		$symminterface{ "0_0_0_0" } = 0;
	} else {
		foreach my $j_symm (0..($nsymm-1)) {
			foreach my $fY_i (@fchaintrace) {
				my $fY_j = vadd( mapply($Rs->[$j_symm],$fY_i) , $Ts->[$j_symm] );
				foreach my $fX_i (@fchaintrace) {
					my $delfXY = vsub( $fY_j,$fX_i );
					my $delfXY_min = vminmod( $delfXY , [1.0,1.0,1.0] );
					#my $delfXY_min = $delfXY;

					# distance check ...
					my $delXY = mapply( $f2c, $delfXY_min );
					my $dist2XY = vnorm2( $delXY );

					if ($dist2XY < $interact_dist*$interact_dist) {
						my $shiftXY = vsub( $delfXY_min , $delfXY );  # should be all integers

						# we have a hit! tag symmop $jsymm and offset $shiftXY as a symmetic interface!
						my $id = $j_symm."_".floor($shiftXY->[0]+0.5)."_".floor($shiftXY->[1]+0.5)."_".floor($shiftXY->[2]+0.5);
						if (!defined $symminterface{ $id }) {
							$symminterface{ $id } = $j_symm+vnorm2($shiftXY);
						}

						# ok ... for really small unit cells (e.g. 1YJP) it may be there is
						#	 another offset within the interaction radius
						#	 but not closer than this offset for every (i,j) pair
						# explicitly search these offsets
						foreach my $T_offset (@{ $Ts_expand }) {
							my $delfXY_min_shift = vadd( $delfXY_min , $T_offset );

							# distance check ...
							$delXY = mapply( $f2c, $delfXY_min_shift );
							$dist2XY = vnorm2( $delXY );
							if ($dist2XY < $interact_dist*$interact_dist) {
								my $newShiftXY = vadd( $T_offset, $shiftXY );
								$id = $j_symm."_".floor($newShiftXY->[0]+0.5)."_".
												  floor($newShiftXY->[1]+0.5)."_".
												  floor($newShiftXY->[2]+0.5);
								if (!defined $symminterface{ $id }) {
									$symminterface{ $id } = $j_symm+vnorm2($newShiftXY);
								}
							}
						}
					}
				}
			}
		}
	}

	# symmetric interfaces may come in pairs
	# find these pairs and place in syminterfaces_paired like [0, 1, 1', 2, 2', ...]
	my @syminterfaces = sort { $symminterface{$a} <=> $symminterface{$b} } keys %symminterface;
	my @syminterfaces_paired;
	my @syminterfaces_unpaired;
	push @syminterfaces_unpaired, $syminterfaces[0];   # first the origin
	my %paired = (0=>1);

	OUTER: foreach my $i (1..$#syminterfaces) {
		next if (defined $paired{$i});

		my ($symm_i,$shiftX_i,$shiftY_i,$shiftZ_i) = split '_',$syminterfaces[ $i ];
		my $R_i = $Rs->[$symm_i];
		my $T_i =vadd( $Ts->[$symm_i], [$shiftX_i,$shiftY_i,$shiftZ_i] );

		foreach my $j ($i+1..$#syminterfaces) {
			my ($symm_j,$shiftX_j,$shiftY_j,$shiftZ_j) = split '_',$syminterfaces[ $j ];
			my $R_j = $Rs->[$symm_j];
			my $T_j =vadd( $Ts->[$symm_j], [$shiftX_j,$shiftY_j,$shiftZ_j] );

			# if transform i is the inverse of transform j we have our pair
			if ( is_inverse( $R_i,$T_i, $R_j,$T_j ) ) {
				$paired{ $j } = 1;
				push @syminterfaces_paired, $syminterfaces[$i];
				push @syminterfaces_paired, $syminterfaces[$j];
				next OUTER;
			}
		}
		push @syminterfaces_unpaired, $syminterfaces[$i];	# unpaired
	}

	#######################################
	##
	## write output symm file
	# symmetry_name c4
	# subunits 4
	# number_of_interfaces 2
	# E = 3*E2
	# anchor_residue 17
	my $symmname = $pdbfile;
	$symmname =~ s/\.pdb$//;
	$symmname = $symmname." $spacegp";
	$symmname =~ s/\s*$//;
	$symmname =~ s/\s/_/g;
	print "symmetry_name $symmname\n";

	#print "E = ".(scalar(@syminterfaces)-1)."*VRT_0_0_0_0_base";
	print "E = 2*VRT_0_0_0_0_base";
	foreach my $i (1..($#syminterfaces_unpaired)) {
		my $estring = " + 1*(VRT_0_0_0_0_base:VRT_".$syminterfaces_unpaired[ $i ]."_base)";
		$estring =~ s/_-(\d)/_n\1/g;
		print $estring;
	}
	foreach my $i (0..($#syminterfaces_paired-1)/2) {
		my $estring = " + 2*(VRT_0_0_0_0_base:VRT_".$syminterfaces_paired[ 2*$i ]."_base)";
		$estring =~ s/_-(\d)/_n\1/g;
		print $estring;
	}
	print "\n";
	#print "anchor_residue COM\n";
	print "anchor_residue 1\n";

	# virtual_coordinates_start
	# xyz VRT1 -1,0,0 0,1,0 0,0,0
	# xyz VRT2 0,-1,0 -1,0,0 0,0,0
	# xyz VRT3 1,0,0 0,-1,0 0,0,0
	# xyz VRT4 0,1,0 1,0,0 0,0,0
	# virtual_coordinates_stop
	print "virtual_coordinates_start\n";
	my @syminterfaces_all;
	push @syminterfaces_all, @syminterfaces_unpaired;
	push @syminterfaces_all, @syminterfaces_paired;
	foreach my $symmkey (@syminterfaces_all) {
		# crystal lattice
		my ($j_symm,$shiftX,$shiftY,$shiftZ) = split '_',$symmkey;

		my $realrot = mmult( $f2c, mmult( $Rs->[$j_symm], $c2f ) );

		my $d = mdet( $realrot );
		if ($d>-0.999 && $d<0.999) {
			die "Unexpected determinant $d (should be +/-1)\n";
		}

		my $xyzline = "xyz VRT_".$symmkey;
		if ($d<-0.999) { $xyzline = "xyzM VRT_".$symmkey; }
		$xyzline =~ s/_-(\d)/_n\1/g;

		# X
		my $fX	= mapply( $c2f , [1,0,0] );
		my $sfX	 = mapply($Rs->[$j_symm],$fX);
		my $sX = mapply( $f2c , $sfX );
		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;

		# Y
		$fX	 = mapply( $c2f , [0,1,0] );
		$sfX  = mapply($Rs->[$j_symm],$fX);
		$sX = mapply( $f2c , $sfX );
		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;

		# orig
		$fX	 = mapply( $c2f , vadd( $CoM, [1,1,1] ) ) ;
		$sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
		$sfX = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
		$sX = mapply( $f2c , $sfX );

		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;
		print "$xyzline\n";

		# centers of mass
		# X
		$xyzline = "xyz VRT_".$symmkey."_base";
		if ($d<-0.999) { $xyzline = "xyzM VRT_".$symmkey."_base"; }
		$xyzline =~ s/_-(\d)/_n\1/g;
		my $fX	= mapply( $c2f , [1,0,0] );
		my $sfX	 = mapply($Rs->[$j_symm],$fX);
		my $sX = mapply( $f2c , $sfX );
		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;

		# Y
		$fX	 = mapply( $c2f , [0,1,0] );
		$sfX  = mapply($Rs->[$j_symm],$fX);
		$sX = mapply( $f2c , $sfX );
		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;

		# orig
		$fX	 = mapply( $c2f , $CoM );
		$sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
		$sfX = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
		$sX = mapply( $f2c , $sfX );

		my $string = sprintf("%.6f,%.6f,%.6f", $sX->[0], $sX->[1], $sX->[2]);
		$xyzline = $xyzline." ".$string;
		print "$xyzline\n";
	}
	print "virtual_coordinates_stop\n";

	# connect_virtual BASEJUMP VRT1 SUBUNIT
	# connect_virtual JUMP2 VRT2 SUBUNIT
	# connect_virtual JUMP10 VRT1 VRT2
	#print "connect_virtual BASEJUMP VRT_0_0_0_0 SUBUNIT\n";  # t_0
	#jump from com to subunit
	foreach my $i (0..$#syminterfaces_all ) {
		my $id = $syminterfaces_all[$i];
		$id =~ s/_-(\d)/_n\1/g;
		print "connect_virtual JUMP_".$id."_to_subunit VRT_".$id."_base SUBUNIT\n";
	}
	#jump to com
	foreach my $i (0..$#syminterfaces_all) {
		my $id = $syminterfaces_all[$i];
		$id =~ s/_-(\d)/_n\1/g;
		print "connect_virtual JUMP_".$id."_to_com VRT_".$id." VRT_".$id."_base\n";
	}
	#jump from base unit
	foreach my $i (1..$#syminterfaces_all ) {
		my $id = $syminterfaces_all[$i];
		$id =~ s/_-(\d)/_n\1/g;
		print "connect_virtual JUMP_".$id." VRT_0_0_0_0 VRT_".$id."\n";
	}
	if ($restrictCrystTrans == 1) {
		print "set_dof JUMP_$syminterfaces_all[0]"."_to_com x y z\n";
	} else {
		if ($cheshire->[0][0] != $cheshire->[0][1] || $cheshire->[1][0] != $cheshire->[1][1] || $cheshire->[2][0] != $cheshire->[2][1] ) {
			print "set_dof JUMP_$syminterfaces_all[0]"."_to_com";
			if ($cheshire->[0][0] != $cheshire->[0][1]) { print " x"; }
			if ($cheshire->[1][0] != $cheshire->[1][1]) { print " y"; }
			if ($cheshire->[2][0] != $cheshire->[2][1]) { print " z"; }
			print "\n";
		}
	}
	print "set_dof JUMP_$syminterfaces_all[0]"."_to_subunit angle_x angle_y angle_z\n";

	#define jumpgroups
	print "set_jump_group JUMPGROUP1 ";
	foreach my $i (0..$#syminterfaces_all ) {
		my $id = $syminterfaces_all[$i];
		$id =~ s/_-(\d)/_n\1/g;
		print " JUMP_".$id."_to_subunit";
	}
	print "\n";
	if ($cheshire->[0][0] != $cheshire->[0][1] || $cheshire->[1][0] != $cheshire->[1][1] || $cheshire->[2][0] != $cheshire->[2][1] ) {
		print "set_jump_group JUMPGROUP2 ";
		foreach my $i (0..$#syminterfaces_all ) {
			my $id = $syminterfaces_all[$i];
			$id =~ s/_-(\d)/_n\1/g;
			print " JUMP_".$id."_to_com";
		}
		print "\n";
	}
	# set_dof 1 x angle_x angle_y angle_z

	########################################
	##
	## write output pdb
	my $outpdb = $pdbfile;
	if ($outpdb =~ /\.pdb$/) {
		$outpdb =~ s/\.pdb$/_symm.pdb/;
	} else {
		$outpdb = $outpdb."_symm.pdb";
	}
	open (OUTPDB, ">$outpdb");


	my $mdlidx = 1;
	my $chnidx = 0;
	print OUTMDL "MODEL $mdlidx\n";
	my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()-=_+;:,.<>";
	my $first_ASU = 1; #True
	my $max_chain_index = -1;
	foreach my $symmkey (@syminterfaces_all) {
		my ($j_symm,$shiftX,$shiftY,$shiftZ) = split '_',$symmkey;
		my $prev_chain_ID = $primary_chain;
		my $current_chain_ID = $primary_chain;
		foreach my $line (@filebuf) {
			# Update current chain ID in original line
			$current_chain_ID = substr ($line, 21, 1);
			if ($current_chain_ID ne $prev_chain_ID) {
				$chnidx++;
			}
			$prev_chain_ID = $current_chain_ID;
			my $linecopy = $line;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
			my $fX = mapply( $c2f , $X );
			my $sfX	 = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
			my $sfX1 = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
			my $sX = mapply( $f2c , $sfX1 );

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $sX->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $sX->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $sX->[2]);
			# keep chain IDs of first ASU
			if (!$first_ASU) {
				substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);
			} else { #is first ASU
				# get chain's chnidx
				my $chain_index = index($chains, $current_chain_ID);
				# store the greater index
				$max_chain_index = $chain_index > $max_chain_index ? $chain_index : $max_chain_index;
			}
			print OUTPDB $linecopy."\n";
		}
		print OUTPDB "TER	\n";
		if ($first_ASU) {
			$chnidx = $max_chain_index;
		}
		$chnidx++;
		if ($chnidx >= length($chains)) {
			$chnidx=0;
			$mdlidx++;
			print OUTPDB "ENDMDL\n";
			print OUTPDB "MODEL $mdlidx\n";
		}
		$first_ASU = 0;
	}
	close(OUTPDB);
}


if ($helix_mode == 1) {
	###############
	### HELIX mode!
	###############

	my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
	my $COM_0 = recenter( $chains{ $primary_chain } );

	## first expand helical symm
	my @helical_chain_split = split( ':', $helical_chain );
	$helical_chain = $helical_chain_split[0];
	my $force_symm_order = 0;
	my $force_rise = 0;
	my $force_axis = '-';

	# optionally allow input to 'force' a symmetric order
	if ($#helical_chain_split > 0) {
		$force_symm_order = $helical_chain_split[1];
	}
	if ($#helical_chain_split > 1) {
		$force_rise = $helical_chain_split[2];
	}
	if ($#helical_chain_split > 2) {
		$force_axis = $helical_chain_split[3];
	}

	## make sure chains exist & are the same length
	if ( ! defined $chains{ $helical_chain } ) {
		die "Chain $helical_chain not in input!\n";
	}
	if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $helical_chain } } ) ) {
		print STDERR "ERROR! chains '$primary_chain' and '$helical_chain' have different residue counts! (".
					 scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $helical_chain } } ).")\n";
		die "Chain length mismatch!\n";
	}

	# superimpose
	my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $helical_chain } );
	my $del_COM = vsub ($COM_i, $COM_0);
	my $subunits_per_turn = 1;

	# special case for identity (fiber symmetry)
	# set helical axis == translation axis
	my ($X,$Y,$Z,$W) = (0,0,0,1);
	my $Wmult = 1;
	if ($force_symm_order == 1) {
		$R = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	}

	if ( is_identity( $R ) && $force_symm_order <= 1  ) {
		($X,$Y,$Z,$W) = ($del_COM->[0],$del_COM->[1],$del_COM->[2],1);
	} else {
		($X,$Y,$Z,$W)=R2quat($R);
		if ($W < 0) { $W = -$W; $Wmult = -1; }

		# force # of subunits
		if ( $force_symm_order > 1 ) {
			$W = cos( PI/$force_symm_order );
		}

		# force axis
		if ( $force_axis eq 'X' ) {
			$X=1; $Y=0; $Z=0;
		} elsif ( $force_axis eq 'Y' ) {
			$X=0; $Y=1; $Z=0;
		} elsif ( $force_axis eq 'Z' ) {
			$X=0; $Y=0; $Z=1;
		}

		my $omega = acos($W);
		$subunits_per_turn = (PI/$omega);
		my $newW = -$Wmult * $W;
		my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
		$R = quat2R( $X*$S , $Y*$S, $Z*$S, $newW );
	}

	# project $del_COM to helical axis
	my $helical_axis = [$X,$Y,$Z];
	normalize( $helical_axis );
	my $del_COM_inplane	   = vsub( $del_COM , vscale(dot($del_COM,$helical_axis),$helical_axis) );
	my $del_COM_alonghelix = vsub( $del_COM , $del_COM_inplane );

	if ($force_rise > 0) {
		my $rise = vnorm($del_COM_alonghelix);
		$del_COM_alonghelix = vscale( $force_rise/$rise , $del_COM_alonghelix );
	}

	# helix handedness
	my $right_handed = 1;
	if ( dot( $del_COM_alonghelix, $helical_axis ) < 0 ) { $right_handed = -1 };

	print STDERR "Found helical symmetry at chain ".$helical_chain."\n";
	print STDERR "   right_handed? ".$right_handed."\n";
	print STDERR "	 subunits per turn = ".$subunits_per_turn."\n";
	print STDERR "	 rise  = ".vnorm($del_COM_alonghelix)."\n";

	my $helix_center;
	my $omega;
	if ( $force_symm_order == 0 ) {
		$omega = acos( $W );
	} else {
		$omega = PI/$force_symm_order;
	}
	# get the center of the helix from angle
	my $helix_R = 0;
	if (sin($omega) > 1e-6) {
		$helix_R = vnorm($del_COM_inplane) / (2*sin( $omega ));
	}

	my $COM_0_5 = vadd( $COM_0 , vscale( 0.5, $del_COM_inplane ) ); # 1/2way between subunits
	my $d_0_5_center = $helix_R * cos( $omega );   # distance from here to helical center

	if ($d_0_5_center == 0) {
		$helix_center = $COM_0_5;
	} else {
		# direction from here to center
		my $center_dir = vscale ( $Wmult*$right_handed , cross( $del_COM_inplane , $del_COM_alonghelix ) );
		normalize( $center_dir );
		$helix_center = vadd( $COM_0_5 , vscale( $d_0_5_center, $center_dir ) );
	}
	my $helix_R = vnorm( vsub( $COM_0, $helix_center ) );

	print STDERR "	 omega = ".(2*$omega * 180.0/PI)." degrees\n";
	print STDERR "	 helical center: ".$helix_center->[0]." , ".$helix_center->[1]." , ".$helix_center->[2]."\n";
	print STDERR "	 helical axis: ".$helical_axis->[0]." , ".$helical_axis->[1]." , ".$helical_axis->[2]."\n";

	# helical NCS
	my $NCS_ops = {};
	$NCS_ops->{R} = [ [1,0,0], [0,1,0], [0,0,1] ];
	$NCS_ops->{T} = $COM_0;
	$NCS_ops->{PATH} = "0";
	$NCS_ops->{CHILDREN} = [];
	$NCS_ops->{AXIS} = [0,0,1];
	my $Toffset = [ 0,0,0 ]; # accounts for deviation from point symm COM and helix axis

	if ($#secondary_chains >= 0) {
		# expand point symmetry (if provided).	First detect symmetries.
		my (@allQs_ncs, @allCOMs_ncs, @sym_orders_ncs, @is_parallel_ncs);

		foreach my $sec_chain (@secondary_chains) {
			my @sec_chain_ids = split( ':', $sec_chain );

			if ( ! defined $chains{ $sec_chain_ids[0] } ) {
				die "Chain $sec_chain not in input!\n";
			}

			## CA check
			if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $sec_chain_ids[0] } } ) ) {
				print STDERR "ERROR! chains '$primary_chain' and '$sec_chain' have different residue counts! (".
							 scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $sec_chain_ids[0] } } ).")\n";
				die "Chain length mismatch!\n";
			}

			# get superposition
			my ($R_ncs,$rmsd_ncs, $COM_i_ncs, $COM_ij_ncs) = rms_align( $chains{ $primary_chain } , $chains{ $sec_chain_ids[0] } );
			if ($quietMode != 1) {
				print STDERR "Aligning $primary_chain and $sec_chain wth RMS=$rmsd.\n";
			}
			if ( is_identity( $R_ncs ) ) {
				print STDERR "Chains $primary_chain and $sec_chain related by transformation only! Aborting.\n";
				exit 1;
			}

			my $del_COM_ncs = vsub ($COM_i_ncs, $COM_0);
			push @allCOMs_ncs, $del_COM_ncs;

			my ($X_ncs,$Y_ncs,$Z_ncs,$W_ncs)=R2quat($R_ncs);

			my $Worig_ncs = $W_ncs;
			my $Wmult_ncs = 1;
			if ($W_ncs < 0) { $W_ncs = -$W_ncs; $Wmult_ncs = -1; }
			my $omega_ncs = acos($W_ncs);
			my $sym_order_ncs = int(PI/$omega_ncs + 0.5);

			# optionally allow input to 'force' a symmetric order
			if ($#sec_chain_ids > 0) {
				$sym_order_ncs = $sec_chain_ids[1];
			}
			push @sym_orders_ncs, $sym_order_ncs;

			# determine if parallel or perpendicular to helix axis
			my $ncs_axis = [$X_ncs,$Y_ncs,$Z_ncs]; normalize( $ncs_axis );
			my $angle = acos( dot($ncs_axis,$helical_axis) );
			my $is_parallel = 0;
			if ($angle > PI/2) {
				$angle = PI - $angle;
			}
			if ($angle < 0.3) {
				$is_parallel = 1;
			} elsif ($angle < PI/2-0.3) {
				print STDERR "Point symmetry between $primary_chain and $sec_chain incompatible with helical symmetry! Aborting.\n";
				print STDERR "Angle = ".($angle*180.0/PI)."\n";
				exit 1;
			}
			push @is_parallel_ncs, $is_parallel;

            if ($is_parallel==0 && $sym_order_ncs!=2) {
				print STDERR "Only a 2-fold axis can be perpendicular to helical axis!\n";
                exit -1;
            }

			my $newW_ncs = -$Wmult_ncs *cos( PI/$sym_order_ncs );
			my $S_ncs = sqrt ( (1-$newW_ncs*$newW_ncs)/($X_ncs*$X_ncs+$Y_ncs*$Y_ncs+$Z_ncs*$Z_ncs) );
			my $newQ_ncs = [$X_ncs*$S_ncs , $Y_ncs*$S_ncs, $Z_ncs*$S_ncs, $newW_ncs];

			push @allQs_ncs, $newQ_ncs;

			if ($quietMode != 1) {
				print STDERR "Found ".$sym_order_ncs."-fold (".(PI/$omega_ncs).") axis to ".$sec_chain_ids[0]." ";
				if ($is_parallel) {
					print STDERR "parallel to helical axis\n";
				} else {
					print STDERR "perpendicular to helical axis\n";
				}
			}
		}

		# Next resolve the point symmetry and helical symmetry.	 There are three cases to handle:
		#	1) Cn symmetry parallel the helical axis
		#	2) C2 symmetry perpindicular to the helical axis
		#	3) Dn symmetry with the n-fold parallel to the helical axis
		if ($#allQs_ncs == 1) {
			# put parallel group (the n-fold) first
			if ($is_parallel_ncs[1]) {
				my $temp;
				$temp = $allQs_ncs[1]; $allQs_ncs[1] = $allQs_ncs[0]; $allQs_ncs[0] = $temp;
				$temp = $allCOMs_ncs[1]; $allCOMs_ncs[1] = $allCOMs_ncs[0]; $allCOMs_ncs[0] = $temp;
				$temp = $sym_orders_ncs[1]; $sym_orders_ncs[1] = $sym_orders_ncs[0]; $sym_orders_ncs[0] = $temp;
				$temp = $secondary_chains[1]; $secondary_chains[1] = $secondary_chains[0]; $secondary_chains[0] = $temp;
				$temp = $is_parallel_ncs[1]; $is_parallel_ncs[1] = $is_parallel_ncs[0]; $is_parallel_ncs[0] = $temp;
			}

			# validate
			if ($is_parallel_ncs[1] || $sym_orders_ncs[1]!=2) {
				print STDERR "Invalid D symmetry! Aborting.\n";
				exit -1;
			}
		}

		# if omega (of the helix) is close to zero, use point symmetry to set axis
		if ( ($omega < 0.05) && $is_parallel_ncs[0]) {
			$helical_axis = [ $allQs_ncs[0][0],$allQs_ncs[0][1],$allQs_ncs[0][2] ];
			normalize( $helical_axis );
			if ($quietMode != 1) {
				print STDERR "UPDATE helical_axis = ".$helical_axis->[0]." , ".$helical_axis->[1]." , ".$helical_axis->[2]."\n";
			}
			$del_COM_inplane	= vsub( $del_COM , vscale(dot($del_COM,$helical_axis),$helical_axis) );
			$del_COM_alonghelix = vsub( $del_COM , $del_COM_inplane );
		} 

		foreach my $i (0..$#allQs_ncs) {
			my ($X_ncs,$Y_ncs,$Z_ncs,$W_ncs)=@{$allQs_ncs[$i]};
			my $axis_orig = [ $X_ncs,$Y_ncs,$Z_ncs ];
			normalize( $axis_orig );

			if ($is_parallel_ncs[$i]) {
				my $newS = sqrt ( (1-$W_ncs*$W_ncs)/
					($helical_axis->[0]*$helical_axis->[0] + 
					 $helical_axis->[1]*$helical_axis->[1] + 
					 $helical_axis->[2]*$helical_axis->[2]) );
				my $newQ = [ $helical_axis->[0]*$newS , $helical_axis->[1]*$newS, $helical_axis->[2]*$newS, $W_ncs];
				$allQs_ncs[$i] = $newQ;
			} else {
				# parallel axis
				my $axis_new = vsub($axis_orig , vscale(dot( $axis_orig, $helical_axis ),$helical_axis));
				normalize( $axis_new );
				my $newS = sqrt ( 1-$W_ncs*$W_ncs );
				my $newQ = [ $axis_new->[0]*$newS , $axis_new->[1]*$newS, $axis_new->[2]*$newS, $W_ncs];
				$allQs_ncs[$i] = $newQ;
			}
		}

		# Finally create the symmetry tree, resolving translation errors
		foreach my $i (0..$#allQs_ncs) {
			my $newQ_ncs = $allQs_ncs[ $i ];
			my $del_COM_ncs = $allCOMs_ncs[ $i ];
			my $sym_order_ncs = $sym_orders_ncs[ $i ];

			my $newR_ncs = quat2R( $newQ_ncs->[0], $newQ_ncs->[1], $newQ_ncs->[2], $newQ_ncs->[3] );

			my $com_complex_ncs = $NCS_ops->{T};
			my $com_secondary_ncs;
			if ($i == 0 ) {
				$com_secondary_ncs = vadd( $del_COM_ncs, $COM_0 );
			} else {
				$com_secondary_ncs = vadd( $del_COM_ncs, vadd( $COM_0, mapply( $newR_ncs , vsub ( $com_complex_ncs , $COM_0 ) ) ) );
			}
			my $newDelCOM_ncs = vsub ( $com_secondary_ncs , $com_complex_ncs );

			# symmetrize newT
			my $err_pos = [0,0,0];
			my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
			foreach my $j (1..$sym_order_ncs) {
				$err_pos = vadd( $err_pos, mapply( $R_i,$newDelCOM_ncs ) );
				$R_i = mmult($newR_ncs, $R_i);
			}

			# special case for Dn symmetry
			if ( $i == 1 && $#allQs_ncs == 1 ) {
				my $axis_proj_i =  [ $allQs_ncs[1-$i]->[0],	 $allQs_ncs[1-$i]->[1],	 $allQs_ncs[1-$i]->[2] ];
				# project $delCOM along this axis
				normalize( $axis_proj_i );
				my $del_COM_inplane_ncs	   = vsub( $newDelCOM_ncs , vscale(dot($newDelCOM_ncs,$axis_proj_i),$axis_proj_i) );
				$err_pos = vscale( $sym_orders_ncs [ $i ] , $del_COM_inplane_ncs );
			}
			if ($quietMode != 1) {
		        print STDERR "        NCS Translation error = ".vnorm( $err_pos )."\n";
			}

			# get the center of the symmgp
			my $adj_newDelCOM_ncs = vsub( $newDelCOM_ncs , 
				[ $err_pos->[0]/$sym_order_ncs , $err_pos->[1]/$sym_order_ncs , $err_pos->[2]/$sym_order_ncs ] );
			my $axis_i =  [ $newQ_ncs->[0], $newQ_ncs->[1], $newQ_ncs->[2] ];
			normalize($axis_i);
			expand_symmops_by_split( $NCS_ops, $newR_ncs, $adj_newDelCOM_ncs, $sym_order_ncs, $axis_i);
		}
		if ($quietMode != 1) {
    		print STDERR "        NCS system center = ".$NCS_ops->{T}->[0]." ".$NCS_ops->{T}->[1]." ".$NCS_ops->{T}->[2]."\n";
        }

        # use point symmetry to reassign helix center:
        if ($is_parallel_ncs[0]) {
    		$helix_center = $NCS_ops->{T};
        } else {
            # a bit confusing: take a point on the helical axis that intersects rotation axis
            my $delR = vsub( $NCS_ops->{T}, $helix_center );
            $helix_center = vadd( $helix_center, vscale( dot($delR,$helical_axis), $helical_axis) )
        }
    	if ($quietMode != 1) {
   	    	print STDERR "NEW: helix_center = ".$helix_center->[0]." , ".$helix_center->[1]." , ".$helix_center->[2]."\n";
        }
	}
	print STDERR "Radius (to CoM) = ".vnorm( vsub( $helix_center, vadd( $Toffset, $COM_0) ) )."\n";

	##########################
	my $omega = acos($W);
	if ($omega < 1e-6) { $omega = pi; }
	my $nsubunits_to_gen = $nturns;

	# store transforms
	my %Rs;
	my %Ts;
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		$Rs{$subunit} = {};
		$Ts{$subunit} = {};
    }

	my $R_helix;
	my $T_helix;
	my $Rinv_helix;
	my $Tinv_helix;

	#transform to helical frame (+ direction)
	$R_helix = [[1,0,0],[0,1,0],[0,0,1]];
	$T_helix = [0,0,0];

	#transfrom to helical frame (- direction)
	$Rinv_helix = [[1,0,0],[0,1,0],[0,0,1]];
	$Tinv_helix = [0,0,0];

	###
	print STDERR "Generating [-".$nsubunits_to_gen." to +".$nsubunits_to_gen."]\n";
	my $symops = tree_traverse( $NCS_ops );
	my $Rinv = minv($R);

	foreach my $subunit (0 .. $nsubunits_to_gen) {
		# gen in + direction
        $Rs{ $subunit }{ 'base' } = $R_helix;
        $Ts{ $subunit }{ 'base' } = vscale( $subunit, $del_COM_alonghelix );
		foreach my $symop (@{ $symops }) {
			my $id = $symop->{PATH};

            # distance from helix axis t
            my $T_local = mapply( $R_helix, vsub($symop->{T}, $helix_center) );
            my $R_local = mmult( $R_helix, $symop->{R} );

            $Rs{ $subunit }{ $id } = $R_local;
			$Ts{ $subunit }{ $id } =
                vadd( $helix_center, vadd( $T_local, vscale( $subunit, $del_COM_alonghelix ) ) );
		}

		# gen in - direction
		if ($subunit != 0) {
            $Rs{ -$subunit }{ 'base' } = $Rinv_helix;
            $Ts{ -$subunit }{ 'base' } = vscale( -$subunit, $del_COM_alonghelix );
			foreach my $symop (@{ $symops }) {
				my $id = $symop->{PATH};

                # distance from helix axis t
                my $T_local = vadd( mapply( $Rinv_helix, vsub($symop->{T}, $helix_center) ) );
                my $R_local = mmult( $Rinv_helix, $symop->{R} );

				$Rs{ -$subunit }{ $id } = $R_local;
				$Ts{ -$subunit }{ $id } =
                    vadd( $helix_center, vadd( $T_local, vscale( -$subunit, $del_COM_alonghelix ) ) );
			}
		}

		# update the transform to the 1st subunit in next helical layer
		$R_helix = mmult( $R, $R_helix );
		$Rinv_helix = mmult( $Rinv, $Rinv_helix );
	}


	#############################
	#############################
	my %symminterface = ();
	my $counter = 0;
	foreach my $subunit (0,-$nsubunits_to_gen .. -1,1..$nsubunits_to_gen) {
		$symminterface{$subunit} = {};
    }
	foreach my $X_i ( @{ $chains{ $primary_chain } } ) {
		foreach my $Y_i (  @{ $chains{ $primary_chain } } ) {
			foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
				foreach my $symop (@{ $symops }) {
					my $id = $symop->{PATH};
					next if (defined $symminterface{ $subunit }{ $id });

					#	x_i = R_i * (x_0 - COM_0) + COM_i
					#	The rms function already ofsets x_0 by -COM_0
					my $rX_i = vadd( $X_i , $COM_0 );
					my $rY_j = vadd( mapply($Rs{ $subunit }{ $id }, $Y_i) , $Ts{ $subunit }{ $id } );
					my $delXY = vsub( $rY_j,$rX_i );
					my $dist2XY = vnorm2( $delXY );

					if ($dist2XY < $interact_dist*$interact_dist) {
						# we have a hit! tag NCS copy $j_symm as a non-symmetic interface
						print STDERR " Adding interface '".$subunit."_".$id."'\n";
						$symminterface{ $subunit }{ $id } = $counter++;
					}
				}
			}
		}
	}

	# find the equation for the energy of the complex
	# simple version: "2" for self, "1" for all others
	my %energy_counter;
    my %syminterfaces_sorted;
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
        $energy_counter{$subunit} = {};
        my @syminterface_s_sorted = 
            sort { $symminterface{$subunit}{$a} <=> $symminterface{$subunit}{$b} } 
            keys %{ $symminterface{$subunit} };
        $syminterfaces_sorted{$subunit} = \@syminterface_s_sorted;
    	foreach my $interface (@syminterface_s_sorted) {
    		$energy_counter{ $subunit }{ $interface } = 1;
        }
	}
   	delete ($energy_counter{0}{ $syminterfaces_sorted{0}[0] });

	## symm file gen
	## write output symm file
	my $symmname = $pdbfile;
	$symmname =~ s/\.pdb$//;
	$symmname = $symmname."_helix_".get_topology( $NCS_ops );
	print "symmetry_name $symmname\n";
	print "E = 2*VRT_0_".$syminterfaces_sorted{0}[0];
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
        my @subs_with_nonzeroE = 
            sort { $symminterface{$subunit}{$a} <=> $symminterface{$subunit}{$b} } 
            keys %{ $energy_counter{$subunit} };
    	foreach my $complex (@subs_with_nonzeroE) {
	    	print " + ".$energy_counter{$subunit}{$complex}."*".
            "(VRT_0_".$syminterfaces_sorted{$subunit}[0].":VRT_".$subunit."_".$complex.")";
    	}
    }
	print "\n";
	print "anchor_residue COM\n";
 
    # generate connections from NCS operators
	my $vrts_by_depth = tree_traverse_by_depth( $NCS_ops );
    my %all_ncs_subunits;
    $counter = 0;
    foreach my $symop (@{ $symops }) {
        $all_ncs_subunits{$symop->{PATH}} = $counter++;
    }
    my ($vrt_lines, $connect_lines, $dof_lines) = fold_tree_from_ncs(
        $NCS_ops, $vrts_by_depth, \%all_ncs_subunits);

	print "virtual_coordinates_start\n";
    print "xyz VRT 1.0,0.0,0.0   0.0,1.0,0.0   0.0,0.0,0.0\n";

    my $refZ = deep_copy($helical_axis);
    my $delR = vsub( $COM_0, $helix_center );
    my $refX = vsub( $delR, vscale( dot($delR,$helical_axis), $helical_axis) );
    normalize($refX);
    my $refY = cross($refX, $refZ);

	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		foreach my $symop (@{ $symops }) {
			my $id = $symop->{PATH};

            my $myOri = vadd( $helix_center, vscale( $subunit, $del_COM_alonghelix ) );
            my $myX = mapply( $Rs{$subunit}{$id}, $refX );
            my $myY = mapply( $Rs{$subunit}{$id}, $refY );

            print "xyz VRT_".$subunit."_".$id
                ." ".$myX->[0].",".$myX->[1].",".$myX->[2]."  "
                ." ".$myY->[0].",".$myY->[1].",".$myY->[2]."  "
                ." ".$myOri->[0].",".$myOri->[1].",".$myOri->[2]."\n";
		}
    }
	print "virtual_coordinates_stop\n";

    my $basetag = $symops->[0]->{PATH};
    print "connect_virtual JUMP0 VRT VRT_".(-$nsubunits_to_gen)."_".$basetag."\n";
    # 1) build array of virtuals running up the helix
	foreach my $subunit (-$nsubunits_to_gen+1 .. $nsubunits_to_gen) {
    	print "connect_virtual JUMP_".$subunit."_".$basetag." VRT_".($subunit-1)."_".$basetag." VRT_".$subunit."_".$basetag."\n";
	}


    # 2) build each layer.  Make sure center subunit is generated first
	foreach my $subunit (0 .. $nsubunits_to_gen) {
		# plus direction
		foreach my $symop (@{ $symops }) {
            if ($symop->{PATH} ne $basetag) {
                print "connect_virtual JUMP_".$subunit."_".$symop->{PATH}." VRT_".$subunit."_".$basetag." VRT_".$subunit."_".$symop->{PATH}."\n";
            }
            if (defined $symminterface{$subunit}{$symop->{PATH}}) {
                print "connect_virtual JUMP_".$subunit."_".$symop->{PATH}."_to_sub VRT_".$subunit."_".$symop->{PATH}." SUBUNIT\n";
            }
        }

        # minus direction
        if ($subunit != 0) {
			foreach my $symop (@{ $symops }) {
	            if ($symop->{PATH} ne $basetag) {
	                print "connect_virtual JUMP_".-$subunit."_".$symop->{PATH}." VRT_".-$subunit."_".$basetag." VRT_".-$subunit."_".$symop->{PATH}."\n";
	            }
	            if (defined $symminterface{$subunit}{$symop->{PATH}}) {
	                print "connect_virtual JUMP_".-$subunit."_".$symop->{PATH}."_to_sub VRT_".-$subunit."_".$symop->{PATH}." SUBUNIT\n";
	            }
	        }
        }
	}

	if ($rbminAll == 1) {
		print "set_dof JUMP0 x y z angle_x angle_y angle_z\n";
	}
    if ($nsubunits_to_gen > 0) {
    	print "set_dof JUMP_".(-$nsubunits_to_gen+1)."_".$basetag." z angle_z\n";
    }
	print "set_dof JUMP_0_".$basetag."_to_sub x y z angle_x angle_y angle_z\n";

    if ($nsubunits_to_gen > 0) {
        print "set_jump_group JUMPGROUP0";
        foreach my $subunit (-$nsubunits_to_gen+1 .. $nsubunits_to_gen) {
            print " JUMP_".$subunit."_".$basetag;
        }
        print "\n";
    }

    print "set_jump_group JUMPGROUP1";
    foreach my $subunit (0, -$nsubunits_to_gen .. -1, 1 .. $nsubunits_to_gen) {
    	foreach my $tag ( keys %{ $symminterface{$subunit} } ) {
        	print " JUMP_".$subunit."_".$tag."_to_sub";
    	}
    }
	print "\n";

	if ($quietMode == 0) {
    	########################################
	    # output complete symm complex
    	########################################
    	my $outpdb = $pdbfile;
	    my $outmon = $pdbfile;
    	my $outmdl = $pdbfile;
 
    	my $suffix = "_model_$primary_chain"."$helical_chain";
    	foreach my $chn (@secondary_chains) {
    		$suffix = $suffix.$chn;
    	}
    	$suffix =~ s/://g;
 
        if ($outpdb =~ /\.pdb$/) {
            $outpdb =~ s/\.pdb$/_symm.pdb/;
            $outmdl =~ s/\.pdb$/$suffix.pdb/;
            $outmon =~ s/\.pdb$/_INPUT.pdb/;
        } else {
            $outpdb = $outpdb."_symm.pdb";
            $outmdl = $outpdb."_model.pdb";
            $outmon = $outmon."_INPUT.pdb";
        }
    	open (OUTPDB, ">$outpdb");
    	open (OUTMON, ">$outmon");
    	open (OUTMDL, ">$outmdl");
 
	    my $chnidx = 0;
    	my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()-=_+;:,.<>";
 
    	foreach my $subunit (0 .. $nsubunits_to_gen) {
        	foreach my $dir (1,-1) {
                next if ($subunit==0 && $dir==-1);
                foreach my $symop (@{ $symops }) {
                    my $id = $symop->{PATH};
                    foreach my $line (@filebuf) {
                        my $linecopy = $line;

                        my $R_i = $Rs{ $dir*$subunit }{ $id };
                        my $T_i = $Ts{ $dir*$subunit }{ $id };

                        my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
                        my $X_0 = vsub($X,$COM_0);
                        my $rX = vadd( mapply($R_i, $X_0) , $T_i );

                        substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
                        substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
                        substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
                        substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);

                        print OUTPDB $linecopy."\n";

                        if (defined $symminterface{ $subunit }{ $id }) {
                            print OUTMDL $linecopy."\n";
                        }
                    }
                    print OUTPDB "TER	\n";
                    print STDERR "Writing interface ".
                        ($dir*$subunit)."_".$id." as chain ".substr ($chains, $chnidx, 1)."\n";
                    $chnidx++;
                }
            }
        }

        ######################################
		foreach my $line (@filebuf) {
			my $linecopy = $line;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $X->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $X->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $X->[2]);

			print OUTMON $linecopy."\n";
		}
		foreach my $line (@altfilebuf) {
			my $linecopy = $line;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $X->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $X->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $X->[2]);

			print OUTMON $linecopy."\n";
		}

        close (OUTPDB);
        close (OUTMON);
        close (OUTMDL);
    }
}


if ($pseudo_mode == 1) {
	###############
	### PSEUDO mode!
	###############
	my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
	my $COM_0 = recenter( $chains{ $primary_chain } );

	my @Rs;
	my @Ts;
	push @Rs, $R_0;
	push @Ts, $COM_0;


	foreach my $sec_chain (sort keys %chains) {
		next if ($sec_chain eq $primary_chain);
		print STDERR "Chain $sec_chain -> $primary_chain\n";

		# superpose, steal R and T
		my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $sec_chain } );
		my $del_COM = vsub ($COM_i, $COM_0);
		my $Rinv = minv( $R );
		push @Rs, $Rinv;
		push @Ts, $COM_i;
	}


	#######################################
	##
	## write output symm file
	# symmetry_name c4
	# subunits 4
	# number_of_interfaces 2
	# E = 3*E2
	# anchor_residue 17
	my $symmname = $pdbfile;
	$symmname =~ s/\.pdb$//;
	$symmname = $symmname."_pseudo".scalar(keys %chains)."fold";
	print "symmetry_name $symmname\n";

	print "E = 2*VRT_0_base";
	foreach my $i (1..($#Rs)) {
		my $estring = " + 1*(VRT_0_base:VRT_".($i)."_base)";
		$estring =~ s/_-(\d)/_n\1/g;
		print $estring;
	}
	print "\n";
	print "anchor_residue COM\n";

	# virtual_coordinates_start
	# xyz VRT1 -1,0,0 0,1,0 0,0,0
	# xyz VRT2 0,-1,0 -1,0,0 0,0,0
	# xyz VRT3 1,0,0 0,-1,0 0,0,0
	# xyz VRT4 0,1,0 1,0,0 0,0,0
	# virtual_coordinates_stop
	print "virtual_coordinates_start\n";
	foreach my $i (0..($#Rs)) {
		# crystal lattice
		my $xyzline = "xyz VRT_".$i;
		$xyzline =~ s/_-(\d)/_n\1/g;

		# X
		my $rX	= mapply( $Rs[$i] , [1,0,0] );
		my $string = sprintf("%.6f,%.6f,%.6f", $rX->[0], $rX->[1], $rX->[2]);
		$xyzline = $xyzline." ".$string;

		# Y
		my $rY	= mapply( $Rs[$i] , [0,1,0] );
		$string = sprintf("%.6f,%.6f,%.6f", $rY->[0], $rY->[1], $rY->[2]);
		$xyzline = $xyzline." ".$string;

		# orig
		my $ori = $Ts[$i];
		$string = sprintf("%.6f,%.6f,%.6f", $ori->[0], $ori->[1], $ori->[2]);
		$xyzline = $xyzline." ".$string;
		print "$xyzline\n";

		# centers of mass
		# X
		my $xyzline = "xyz VRT_".$i."_base";
		$xyzline =~ s/_-(\d)/_n\1/g;
		my $string = sprintf("%.6f,%.6f,%.6f", $rX->[0], $rX->[1], $rX->[2]);
		$xyzline = $xyzline." ".$string;

		# Y
		$string = sprintf("%.6f,%.6f,%.6f", $rY->[0], $rY->[1], $rY->[2]);
		$xyzline = $xyzline." ".$string;

		# orig
		my $ori = $Ts[$i];
		$string = sprintf("%.6f,%.6f,%.6f", $ori->[0], $ori->[1], $ori->[2]);
		$xyzline = $xyzline." ".$string;
		print "$xyzline\n";
	}
	print "virtual_coordinates_stop\n";

	# connect_virtual BASEJUMP VRT1 SUBUNIT
	# connect_virtual JUMP2 VRT2 SUBUNIT
	# connect_virtual JUMP10 VRT1 VRT2
	#print "connect_virtual BASEJUMP VRT_0_0_0_0 SUBUNIT\n";  # t_0
	#jump from com to subunit
	foreach my $i (0..($#Rs)) {
		print "connect_virtual JUMP_".$i."_to_subunit VRT_".$i."_base SUBUNIT\n";
	}
	#jump to com
	foreach my $i (0..($#Rs)) {
		print "connect_virtual JUMP_".$i."_to_com VRT_".$i." VRT_".$i."_base\n";
	}
	#jump from base unit
	foreach my $i (1..($#Rs)) {
		print "connect_virtual JUMP_".$i." VRT_0 VRT_".$i."\n";
	}
	print "set_dof JUMP_0_to_com x y z\n";
	print "set_dof JUMP_0_to_subunit angle_x angle_y angle_z\n";

	#define jumpgroups
	print "set_jump_group JUMPGROUP1 ";
	foreach my $i (0..($#Rs)) {
		print " JUMP_".$i."_to_subunit";
	}
	print "\n";
	print "set_jump_group JUMPGROUP2 ";
	foreach my $i (0..($#Rs)) {
		print " JUMP_".$i."_to_com";
	}
	print "\n";


	########################################
	##
	## write output pdb
	my $outpdb = $pdbfile;
	my $outmon = $pdbfile;
	if ($outpdb =~ /\.pdb$/) {
		$outpdb =~ s/\.pdb$/_symm.pdb/;
		$outmon =~ s/\.pdb$/_INPUT.pdb/;
	} else {
		$outpdb = $outpdb."_symm.pdb";
		$outmon = $outpdb."_INPUT.pdb";
	}
	open (OUTPDB, ">$outpdb");
	open (OUTMON, ">$outmon");
	my $chnidx = 0;
	my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()-=_+;:,.<>";
	foreach my $i (0..($#Rs)) {
		foreach my $line (@filebuf) {
			my $linecopy = $line;

			my $X = vsub([substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)] , $COM_0);
			my $rX = vadd( mapply($Rs[$i],$X) , $Ts[$i] );

			substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
			substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
			substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
			substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);

			print OUTPDB $linecopy."\n";
		}
		print OUTPDB "TER	\n";
		$chnidx++;
	}
	close(OUTPDB);

	######################################

	foreach my $line (@filebuf) {
		my $linecopy = $line;

		my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];

		substr ($linecopy, 30, 8) = sprintf ("%8.3f", $X->[0]);
		substr ($linecopy, 38, 8) = sprintf ("%8.3f", $X->[1]);
		substr ($linecopy, 46, 8) = sprintf ("%8.3f", $X->[2]);
		substr ($linecopy, 21, 1) = "A";

		print OUTMON $linecopy."\n";
	}

}


#############################################################################
#############################################################################
###
###	 NCS symmetry tree generation fns.
###
#############################################################################
#############################################################################

# expand_symmops_by_split( $NCS_ops, $newR, $adj_newDelCOM, $sym_order)
sub expand_symmops_by_split {
	my ( $tree, $newR, $newDelT, $sym_order, $axis_i ) = @_;

	my $newNCSops = {};
	my $COM_0 = [0 , 0 , 0];
	$newNCSops->{R} = [ [1,0,0], [0,1,0], [0,0,1] ];
	$newNCSops->{CHILDREN} = [];
	$newNCSops->{AXIS} = deep_copy($axis_i); #axis of children

	my $COM_i = [ 0,0,0 ];
	my $R_i	  = [ [1,0,0], [0,1,0], [0,0,1] ];
	my $newCOM0 = [0,0,0];

	foreach my $i (0..$sym_order-1) {
		#my $newNCSops_i = deep_copy( $NCS_ops );
		my $newNCSops_i = deep_copy( $tree );

		# rotate about the center of mass of the subtree
		#	 then translate to the new CoM
		#apply_transformation ( $newNCSops_i, $R_i, $NCS_ops->{T}, $COM_i, $i );
		apply_transformation ( $newNCSops_i, $R_i, $tree->{T}, $COM_i, $i );
		$newNCSops_i->{AXIS} = deep_copy($axis_i);
		push @{ $newNCSops->{CHILDREN} }, $newNCSops_i;

		$newCOM0 = vadd( $newCOM0, $newNCSops_i->{T} );

		$COM_i = vadd( $COM_i , mapply( $R_i, $newDelT ) );
		$R_i = mmult( $newR , $R_i );
	}
	$newNCSops->{T} = [ $newCOM0->[0]/$sym_order , $newCOM0->[1]/$sym_order , $newCOM0->[2]/$sym_order ];
	$newNCSops->{PATH} = "";
	$_[0] = $newNCSops;
}


#############
# apply a transform ... rotate about $T_about, applying a post_transform $T_post
sub apply_transformation {
	my ($tree, $R, $T_about, $T_post, $prefix) = @_;

	$tree->{R} = mmult( $R, $tree->{R} );
	$tree->{T} = vadd( vadd($T_about, $T_post) , mapply( $R, vsub( $tree->{T} , $T_about ) ) );
	my $newPath = $prefix;
	if (length($tree->{PATH}) > 0) { $newPath = $newPath.'_'.$tree->{PATH}; }
	$tree->{PATH} = $newPath;
	foreach my $child ( @{ $tree->{CHILDREN} } ) {
		apply_transformation( $child, $R, $T_about, $T_post, $prefix);
	}
}


##########
## ($nnodes,$nleaves) = tree_size( $tree )
##	  get the size of a tree
sub tree_size {
	my $tree = shift;
	my ($nnodes,$nleaves) = (0,0);
	tree_size_recursive( $tree, $nnodes,$nleaves);
	return ($nnodes,$nleaves);
}

sub tree_size_recursive {
	my ($tree, $nnodes, $nleaves) = @_;
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );

	$nnodes++;
	if ($nchildren == 0) {
		$nleaves++;
	} else {
		foreach my $child ( @{ $tree->{CHILDREN} } ) {
			tree_size_recursive( $child, $nnodes, $nleaves);
		}
	}

	# pass-by-ref
	@_[1] = $nnodes;
	@_[2] = $nleaves;
}


############
## my $leaves = tree_traverse( $tree)
##	 traverse the tree
sub tree_traverse {
	my $tree = shift;
	my $leaves = [];
	tree_traverse_recursive( $tree, $leaves );
	return $leaves;
}

sub tree_traverse_recursive {
	my ($tree,$leaves) = @_;
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );
	if ($nchildren == 0) {
		push @{ $leaves }, $tree;
	} else {
		foreach my $child ( @{ $tree->{CHILDREN} } ) {
			tree_traverse_recursive( $child, $leaves );
		}
	}
}

############
## my $vrts_by_depth = tree_traverse_by_depth( $NCS_ops );
##	   traverse every node in the tree returning them sorted by depth
sub tree_traverse_by_depth {
	my $tree = shift;
	my $depth = get_depth( $tree );
	my $nodes_by_depth = [	];
	tree_traverse_by_depth_recursive( $tree, $nodes_by_depth, 0 );
	return $nodes_by_depth;
}

sub tree_traverse_by_depth_recursive {
	my ( $tree, $nodes_by_depth , $curr_depth) = @_;
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );

	if (!defined $nodes_by_depth->[$curr_depth] ) {
		$nodes_by_depth->[$curr_depth] = [];
	}
	push @{ $nodes_by_depth->[$curr_depth] }, $tree->{PATH};

	foreach my $child ( @{ $tree->{CHILDREN} } ) {
		tree_traverse_by_depth_recursive( $child, $nodes_by_depth , $curr_depth+1 );
	}
}


##########
## get a string describing the topology of the tree
sub get_topology {
	my $tree = shift;
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );
	my $retval = "";
	if ($nchildren != 0) {
		# just look at first child
		$retval = $retval."_$nchildren".get_topology( $tree->{CHILDREN}->[0] );
	}
	return $retval;
}

##########
## get the depth
sub get_depth {
	my $tree = shift;
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );
	my $retval = 0;
	if ($nchildren != 0) {
		$retval = 1+get_depth( $tree->{CHILDREN}->[0] );
	}
	return $retval;
}

#########
## access a subtree
sub get_subtree {
	my ($tree,$accessor_string) = @_;
	my @accessor_list = split '_',$accessor_string;
	#shift @accessor_list; # pop initial '0' off accessor list
	my $subtree_j = get_subtree_recursive( $tree, \@accessor_list );
}

sub get_subtree_recursive {
	my $tree = shift;
	my $list = shift;
	if (scalar( @{ $list } ) == 0 || $list->[0] eq "" ) {
		return $tree;
	} else {
		my $idx = shift @{ $list };
		return get_subtree_recursive( $tree->{CHILDREN}->[$idx] , $list );
	}
	print "ERROR!  Subtree undefined\n";
	exit -1 ;
}


#########
## ($vrt_lines, $connect_lines, $dof_lines) = fold_tree_from_ncs( $NCS_ops );
##
## This function recreates an NCS coordinate system from an NCS tree
sub fold_tree_from_ncs {
	my ($tree, $nodes_by_depth, $connected_subunits ) = @_;

	# root doesnt have parents or siblings
	# use 1st child instead
	my $axis = $tree->{AXIS};

	# parent com (arbitrarily) perpendicular to rot axis
	my $paxis;
	if ($axis->[1] != 0 || $axis->[2] != 0) {
		$paxis = cross([1,0,0], $axis);
	} else {
		$paxis = cross([0,1,0], $axis);
	}
	my $parent_com = vadd( $paxis ,get_com( $tree ) );

	my $nsiblings = 1;
	my $vrt_lines = [];
	my $connect_lines = [];
	my $dof_lines = [];

	fold_tree_from_ncs_recursive( $tree , $parent_com, $nsiblings, $vrt_lines,
								  $connect_lines, $dof_lines, $nodes_by_depth, $connected_subunits );
	return ( $vrt_lines, $connect_lines, $dof_lines);
}


sub fold_tree_from_ncs_recursive {
	my ($tree,$parent_com,$nsiblings,
		$vrt_lines,$connect_lines,$dof_lines,$nodes_by_depth,$connected_subunits) = @_;

	my $origin = get_com( $tree );
	my $id_string = $tree->{PATH};

	# x points from origin to parent CoM
	my $myX = vsub( $parent_com , $origin );
	my $myZ = $tree->{AXIS}; #mapply( $tree->{R}, [0,0,1]);

	# y is whatever is left
	my $myY = cross( $myZ, $myX );

	# at top level, if symm axis is X, this can happen
	if (abs($myY->[0]) < 1e-6 && abs($myY->[1]) < 1e-6 && abs($myY->[2]) < 1e-6 ) {
		$myZ = mapply( $tree->{R}, [1,0,0]);
		$myY = cross( $myZ, $myX );
	}
	normalize( $myX );
	normalize( $myY );

	#
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );

	## recursive call
	foreach my $child_idx ( 0..$nchildren-1 ) {
		my $child_idstring = $tree->{CHILDREN}->[$child_idx]->{PATH};

		# now recursively call each child
		fold_tree_from_ncs_recursive( $tree->{CHILDREN}->[$child_idx], $origin, $nchildren ,
									  $vrt_lines, $connect_lines, $dof_lines, $nodes_by_depth,$connected_subunits);
	}

	# xyz VRT1 -1,0,0 0,1,0 0,0,0
	push @{ $vrt_lines }, "xyz VRT$id_string  ".
				sprintf("%.7f,%.7f,%.7f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.7f,%.7f,%.7f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.7f,%.7f,%.7f", $parent_com->[0], $parent_com->[1], $parent_com->[2]);

	## set up vrts + jumps
	if ($nchildren > 0) {
		## jump to first child
		my $child0_idstring = $tree->{CHILDREN}->[0]->{PATH};
		push @{ $connect_lines }, "connect_virtual JUMP$child0_idstring VRT$id_string VRT$child0_idstring";

		# is this jump the controlling jump?
		my $is_controlling = -1;
		my $depth = 0;
		foreach my $nodes_i (@{ $nodes_by_depth }) {
			if ($child0_idstring eq $nodes_i->[0]) {
				$is_controlling = $depth;
			}
			$depth++;
		}

		if ($is_controlling >= 1) {	 # level 1 jumps don't move, hence '>'
			my $x_dist = vnorm( vsub( $parent_com,$origin) ) ;

			# in icosohedral (and maybe octohedral/tetraherdal) symm we have
			#	  2-fold symm operators where x==0 and must ==0 to keep other symmetries intact
			# for these cases we cannot allow x/angle_x to move
			if ( $x_dist > 1e-3 ) {
				if ($nsiblings == 2) {
					push @{ $dof_lines },	  "set_dof JUMP$child0_idstring x($x_dist) angle_x";
					#push @{ $dof_lines },	   "set_dof JUMPGROUP$is_controlling x angle_x";
				} elsif ($nsiblings > 2) {
					push @{ $dof_lines },	  "set_dof JUMP$child0_idstring x($x_dist)";
					#push @{ $dof_lines },	   "set_dof JUMPGROUP$is_controlling x";
				}
			}
		}

		foreach my $child_idx ( 1..$nchildren-1 ) {
			my $child_idstring = $tree->{CHILDREN}->[$child_idx]->{PATH};
			push @{ $connect_lines }, "connect_virtual JUMP$child_idstring VRT$child0_idstring VRT$child_idstring";
		}
	} else { # $nchildren == 0
		##
		## now vrts + jumps for child nodes
		# 1 -- jump to the COM
		# is this jump the controlling jump?
		my $is_controlling = -1;
		my $depth = 0;
		foreach my $nodes_i (@{ $nodes_by_depth }) {
			if ($id_string eq $nodes_i->[0]) {
				$is_controlling = $depth;
			}
			$depth++;
		}
		push @{ $connect_lines }, "connect_virtual JUMP".$id_string."_to_com VRT$id_string VRT$id_string"."_base";

		if ($is_controlling >= 1) {
			my $x_dist = vnorm( vsub( $parent_com,$origin) ) ;
			# in icosohedral (and maybe octohedral/tetraherdal) symm we have
			#	  2-fold symm operators where x==0 and must ==0 to keep other symmetries intact
			# for these cases we cannot allow x/angle_x to move
			if ( $x_dist > 1e-3 ) {
				if ($nsiblings == 2) {
					push @{ $dof_lines },	  "set_dof JUMP$id_string"."_to_com x($x_dist) angle_x";
				} elsif ($nsiblings > 2) {
					push @{ $dof_lines },	  "set_dof JUMP$id_string"."_to_com x($x_dist)";
				}
			}
		}

		# 2 -- if an interface res, jump to the subunit
		if (defined $connected_subunits->{ $tree->{PATH} }) {
			# if this is "subunit 0"
			if ($connected_subunits->{ $tree->{PATH} } == 0) {
				push @{ $dof_lines },	  "set_dof JUMP".$id_string."_to_subunit angle_x angle_y angle_z";
			}

			push @{ $connect_lines }, "connect_virtual JUMP".$id_string."_to_subunit VRT$id_string"."_base SUBUNIT";

		}


		push @{ $vrt_lines }, "xyz VRT$id_string"."_base  ".
					sprintf("%.7f,%.7f,%.7f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.7f,%.7f,%.7f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.7f,%.7f,%.7f", $origin->[0], $origin->[1], $origin->[2]);

	}
}

##########
## my $x = get_com( $subtree )
##	  get the center of mass of a subtree
##	  do this by summing the CoM's of all subunits
sub get_com {
	my $tree = shift;
	return $tree->{T};
}


###################################################################################
###################################################################################
###
### Kabsch fast RMS alignment
###
###################################################################################
###################################################################################

# my ($R,$rmsd, $Ycom, $Ycom_to_Xcom) = rms_align( $x,$y );
sub rms_align {
	my ($X,$Y) = @_;

	my ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $X, $Y );
	my ($U,$residual) = calculate_rotation_matrix($R,$E0);

	my $rmsd = sqrt( fabs($residual*2.0/$nlist) );
	return ($U,$rmsd,$mov_com, $mov_to_ref);
}

# my $com = recenter( $x );
sub recenter {
	my ($X) = @_;
	my $Natms = scalar(@{ $X });
	my $com = [0,0,0];

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$com->[$i] += $X->[$n][$i];
		}
	}

	foreach my $i (0..2) {
		$com->[$i] /= $Natms;
	}

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$X->[$n][$i] -= $com->[$i];
		}
	}
	return $com;
}


# normalize($a)
sub normalize {
	my $a = shift;
	my $b = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
	if ($b > 1e-6) {
		$a->[0] /= $b; $a->[1] /= $b; $a->[2] /= $b;
	}
}


# my $a_dot_b = dot($a,$b)
sub dot {
	my ($a,$b) = @_;
	return ($a->[0]*$b->[0] + $a->[1]*$b->[1] + $a->[2]*$b->[2]);
}


# my $a = cross ( b , c )
sub cross {
	my ($b,$c) = @_;
	my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
			  $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
			  $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
	return $a;
}


# ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $ref_xlist, $mov_xlist )
sub setup_rotation {
	my ( $ref_xlist, $mov_xlist ) = @_;

	my $nlist = min( scalar(@{ $ref_xlist }) , scalar(@{ $mov_xlist }) );
	my $ref_com = [0,0,0];
	my $mov_com = [0,0,0];
	my $mov_to_ref = [0,0,0];

	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_com->[$i] += $mov_xlist->[$n][$i];
			$ref_com->[$i] += $ref_xlist->[$n][$i];
		}
	}
	foreach my $i (0..2) {
		$mov_com->[$i] /= $nlist;
		$ref_com->[$i] /= $nlist;
		$mov_to_ref->[$i] = $ref_com->[$i] - $mov_com->[$i];
	}


	# shift mov_xlist and ref_xlist to centre of mass */
	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_xlist->[$n][$i] -= $mov_com->[$i];
			$ref_xlist->[$n][$i] -= $ref_com->[$i];
		}
	}

	# initialize
	my $R = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $E0 = 0.0;

	foreach my $n (0..$nlist-1) {
		# E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n)
		foreach my $i (0..2) {
		  $E0 +=  $mov_xlist->[$n][$i] * $mov_xlist->[$n][$i]
				+ $ref_xlist->[$n][$i] * $ref_xlist->[$n][$i];
		}

		# R[i,j] = sum(over n): y(n,i) * x(n,j)
		foreach my $i (0..2) {
			foreach my $j (0..2) {
			   $R->[$i][$j] += $mov_xlist->[$n][$i] * $ref_xlist->[$n][$j];
			}
		}
	}
	$E0 *= 0.5;

	return ($nlist,$mov_com, $mov_to_ref, $R, $E0);
}

# helper funct
sub j_rotate {
	my ($a,$i,$j,$k,$l,$s,$tau) = @_;
	my $g = $a->[$i][$j];
	my $h = $a->[$k][$l];
	$a->[$i][$j] = $g-$s*($h+$g*$tau);
	$a->[$k][$l] = $h+$s*($g-$h*$tau);
}

# ($d,$v,$nrot) = jacobi3($a)
#	 computes eigenval and eigen_vec of a real 3x3
#	 symmetric matrix. On output, elements of a that are above
#	 the diagonal are destroyed. d[1..3] returns the
#	 eigenval of a. v[1..3][1..3] is a matrix whose
#	 columns contain, on output, the normalized eigen_vec of a.
#	 n_rot returns the number of Jacobi rotations that were required
sub jacobi3 {
	my $a = shift;

	my $v = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $b = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $d = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $z = [0,0,0];
	my $n_rot = 0;
	my $thresh = 0;

	# 50 tries!
	foreach my $count (0..49) {

		# sum off-diagonal elements
		my $sum = fabs($a->[0][1])+fabs($a->[0][2])+fabs($a->[1][2]);

		# if converged to machine underflow
		if ($sum == 0.0) {
		  return($d,$v,$n_rot);
		}

		# on 1st three sweeps..
		my $thresh = 0;
		if ($count < 3) {
			$thresh = $sum * 0.2 / 9.0;
		}

		foreach my $i (0,1) {
			foreach my $j ($i+1..2) {
				my $g = 100.0 * fabs($a->[$i][$j]);

				# after four sweeps, skip the rotation if
				# the off-diagonal element is small
				if ( $count > 3
					  && fabs($d->[$i])+$g == fabs($d->[$i])
					  && fabs($d->[$j])+$g == fabs($d->[$j]) ) {
					$a->[$i][$j] = 0.0;
				} elsif (fabs($a->[$i][$j]) > $thresh) {
					my $h = $d->[$j] - $d->[$i];
					my ($t,$s,$tau,$theta);

					if (fabs($h)+$g == fabs($h)) {
						$t = $a->[$i][$j] / $h;
					} else {
						$theta = 0.5 * $h / ($a->[$i][$j]);
						$t = 1.0 / ( fabs($theta) + sqrt(1.0 + $theta*$theta) );
						if ($theta < 0.0) { $t = -$t; }
					}

					my $c = 1.0 / sqrt(1 + $t*$t);
					$s = $t * $c;
					$tau = $s / (1.0 + $c);
					$h = $t * $a->[$i][$j];

					$z->[$i] -= $h;
					$z->[$j] += $h;
					$d->[$i] -= $h;
					$d->[$j] += $h;

					$a->[$i][$j] = 0.0;

					foreach my $k (0..$i-1) {
						j_rotate($a, $k, $i, $k, $j, $s, $tau);
					}
					foreach my $k ($i+1..$j-1) {
						j_rotate($a, $i, $k, $k, $j, $s, $tau);
					}
					foreach my $k ($j+1..2) {
						j_rotate($a, $i, $k, $j, $k, $s, $tau);
					}
					foreach my $k (0..2) {
						j_rotate($v, $k, $i, $k, $j, $s, $tau);
					}
					$n_rot++;
				}
			}
		}

		foreach my $i (0..2) {
			$b->[$i] += $z->[$i];
			$d->[$i] = $b->[$i];
			$z->[$i] = 0.0;
		}
	}

	print STDERR "WARNING: Too many iterations in jacobi3!	You're bad and you should feel bad.\n";
	exit -1;
}


# ($eigen_vec, $eigenval) = diagonalize_symmetric( $matrix )
sub diagonalize_symmetric {
	my $matrix = shift;
	my $n_rot = 0;

	my ($eigenval,$vec,$n_rot) = jacobi3($matrix);

	# sort solutions by eigenval
	foreach my $i (0..2) {
		my $k = $i;
		my $val = $eigenval->[$i];

		foreach my $j ($i+1..2) {
			if ($eigenval->[$j] >= $val) {
				$k = $j;
				$val = $eigenval->[$k];
			}
		}

		if ($k != $i) {
			$eigenval->[$k] = $eigenval->[$i];
			$eigenval->[$i] = $val;
			foreach my $j (0..2) {
				$val = $vec->[$j][$i];
				$vec->[$j][$i] = $vec->[$j][$k];
				$vec->[$j][$k] = $val;
			}
		}
	}

	# transpose
	my $eigen_vec = [ [$vec->[0][0],$vec->[1][0],$vec->[2][0]] ,
					  [$vec->[0][1],$vec->[1][1],$vec->[2][1]] ,
					  [$vec->[0][2],$vec->[1][2],$vec->[2][2]] ];
	return ($eigen_vec, $eigenval);
}


# ($U,$residual) = calculate_rotation_matrix($R,$E0)
sub calculate_rotation_matrix {
	my ($R,$E0) = @_;

	my $RtR = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $left_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $right_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $eigenval = 0;

	 # Rt <- transpose of R
	my $Rt = [ [$R->[0][0],$R->[1][0],$R->[2][0]] ,
			   [$R->[0][1],$R->[1][1],$R->[2][1]] ,
			   [$R->[0][2],$R->[1][2],$R->[2][2]] ];

	# make symmetric RtR = Rt X R
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$RtR->[$i][$j] = 0.0;
			foreach my $k (0..2) {
				$RtR->[$i][$j] += $Rt->[$k][$i] * $R->[$j][$k];
			}
		}
	}

	($right_eigenvec, $eigenval) = diagonalize_symmetric( $RtR );

	# right_eigenvec's should be an orthogonal system but could be left
	#	or right-handed. Let's force into right-handed system.
	$right_eigenvec->[2] = cross($right_eigenvec->[0], $right_eigenvec->[1]);

	# From the Kabsch algorithm, the eigenvec's of RtR
	#	are identical to the right_eigenvec's of R.
	#	This means that left_eigenvec = R x right_eigenvec
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$left_eigenvec->[$i][$j] = dot($right_eigenvec->[$i], $Rt->[$j]);
		}
	}

	foreach my $i (0..2) {
		normalize($left_eigenvec->[$i]);
	}

	# Force left_eigenvec[2] to be orthogonal to the other vectors.
	# First check if the rotational matrices generated from the
	#	orthogonal eigenvectors are in a right-handed or left-handed
	#	coordinate system - given by sigma. Sigma is needed to
	#	resolve this ambiguity in calculating the RMSD.
	my $sigma = 1.0;
	my $v = cross($left_eigenvec->[0], $left_eigenvec->[1]);
	if (dot($v, $left_eigenvec->[2]) < 0.0) {
		$sigma = -1.0;
	}
	foreach my $i (0..2) {
		$left_eigenvec->[2][$i] = $v->[$i];
	}

	# calc optimal rotation matrix U that minimises residual
	my $U = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			foreach my $k (0..2) {
				$U->[$i][$j] += $left_eigenvec->[$k][$i] * $right_eigenvec->[$k][$j];
			}
		}
	}

	my $residual = $E0 - sqrt(fabs($eigenval->[0]))
					   - sqrt(fabs($eigenval->[1]))
					   - $sigma * sqrt(fabs($eigenval->[2]));

	return ($U,$residual);
}


###################################################################################
###################################################################################
###
### Vector and matrix ops
###
###################################################################################
###################################################################################

sub deep_copy {
	my $this = shift;
	if (not ref $this) {
		$this;
	} elsif (ref $this eq "ARRAY") {
		[map deep_copy($_), @$this];
	} elsif (ref $this eq "HASH") {
		+{map { $_ => deep_copy($this->{$_}) } keys %$this};
	} else { die "what type is $_?" }
}


# rotation from euler angles
sub euler {
	my ($aa, $bb, $gg) = @_;
	my $MM;

	$MM->[0][0] = (-sin($aa)*cos($bb)*sin($gg) + cos($aa)*cos($gg));
	$MM->[0][1] = ( cos($aa)*cos($bb)*sin($gg) + sin($aa)*cos($gg));
	$MM->[0][2] = ( sin($bb)*sin($gg));
	$MM->[1][0] = (-sin($aa)*cos($bb)*cos($gg) - cos($aa)*sin($gg));
	$MM->[1][1] = ( cos($aa)*cos($bb)*cos($gg) - sin($aa)*sin($gg));
	$MM->[1][2] = ( sin($bb)*cos($gg));
	$MM->[2][0] = ( sin($aa)*sin($bb));
	$MM->[2][1] = (-cos($aa)*sin($bb));
	$MM->[2][2] = ( cos($bb));

	return $MM;
}

# my ($X,$Y,$Z,$W)=R2quat($M)
sub R2quat {
	my $R = shift;
	my ($S,$X,$Y,$Z,$W);
	if ( $R->[0][0] > $R->[1][1] && $R->[0][0] > $R->[2][2] )  {
		$S	= sqrt( 1.0 + $R->[0][0] - $R->[1][1] - $R->[2][2] ) * 2;
		$X = 0.25 * $S;
		$Y = ($R->[1][0] + $R->[0][1] ) / $S;
		$Z = ($R->[2][0] + $R->[0][2] ) / $S;
		$W = ($R->[2][1] - $R->[1][2] ) / $S;
	} elsif ( $R->[1][1] > $R->[2][2] ) {
		$S	= sqrt( 1.0 + $R->[1][1] - $R->[0][0] - $R->[2][2] ) * 2;
		$X = ($R->[1][0] + $R->[0][1] ) / $S;
		$Y = 0.25 * $S;
		$Z = ($R->[2][1] + $R->[1][2] ) / $S;
		$W = ($R->[0][2] - $R->[2][0] ) / $S;
	} else {
		$S	= sqrt( 1.0 + $R->[2][2] - $R->[0][0] - $R->[1][1] ) * 2;
		$X = ($R->[0][2] + $R->[2][0] ) / $S;
		$Y = ($R->[2][1] + $R->[1][2] ) / $S;
		$Z = 0.25 * $S;
		$W = ($R->[1][0] - $R->[0][1]) / $S;
	}
	return ($X,$Y,$Z,$W);
}

# my ($R)=R2quat($X,$Y,$Z,$W)
sub quat2R {
	my ($X,$Y,$Z,$W) = @_;
	my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
	my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
	my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
	my $R = [ [ 1 - 2 * ( $yy+$zz ) ,	  2 * ( $xy-$zw ) ,		2 * ( $xz+$yw ) ] ,
			  [		2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,		2 * ( $yz-$xw ) ] ,
			  [		2 * ( $xz-$yw ) ,	  2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
	return $R;
}

# my ($R)=R2quat($X,$Y,$Z,$W)
sub quatnorm {
	my ($X,$Y,$Z,$W) = @_;
	my $S = sqrt( $X*$X+$Y*$Y+$Z*$Z+$W*$W );
	return [ $X/$S , $Y/$S , $Z/$S , $W/$S ];
}

#####################################
#####################################

# vector addition
sub vstr {
	my ($x) = @_;
	return "[".($x->[0]).",".($x->[1]).",".($x->[2])."]";
}

# vector addition
sub vadd {
	my ($x, $y) = @_;
	return [ $x->[0]+$y->[0], $x->[1]+$y->[1], $x->[2]+$y->[2] ];
}

# vector subtraction
sub vsub {
	my ($x, $y) = @_;
	return [ $x->[0]-$y->[0], $x->[1]-$y->[1], $x->[2]-$y->[2] ];
}

# mult vector by scalar
sub vscale {
	my ($x, $y) = @_;
	return [ $x*$y->[0], $x*$y->[1], $x*$y->[2] ];
}

# "min mod"
sub minmod {
	my ($x,$y) = @_;
	my $r = fmod($x,$y);
	if ($r < -fabs( $y/2.0 ) ) { $r += fabs( $y ); }
	elsif ($r >	 fabs( $y/2.0 ) ) { $r -= fabs( $y ); }
	return $r;
}

# vector min-modulus
sub vminmod {
	my ($x,$y) = @_;
	return [ minmod($x->[0],$y->[0]), minmod($x->[1],$y->[1]), minmod($x->[2],$y->[2]) ];
}


#####################################
#####################################

# raise a matrix to a power
# dumb way of doing it
sub mpow {
	my ($mat, $pow) = @_;
	my $matpow = $mat;
	foreach my $i (2..$pow) {
		$matpow = mmult( $mat, $matpow );
	}
	return $matpow;
}

# matrix x vector mult
sub mapply {
	my ($rotmat, $cart) = @_;
	my $out = [0, 0, 0];
	my ($i, $j);
	for ($i=0; $i < 3; ++$i) {
		for ($j=0; $j < 3; ++$j) {
			$out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
		}
	}
	return $out;
}

# matrix x matrix mult
sub mmult {
	my ($m1, $m2) = @_;
	my $out = [ [0,0,0], [0,0,0], [0,0,0] ];
	my ($i, $j, $k);
	for ($i=0; $i<3; ++$i) {
		for ($j=0; $j<3; ++$j) {
			for ($k=0; $k<3; ++$k) {
				$out->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];
			}
		}
	}
	return $out;
}


# matrix inversion
sub minv {
	my $M = shift;
	my $Minv = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $D = $M->[0][0] * ( $M->[1][1]*$M->[2][2] - $M->[2][1]*$M->[1][2] ) -
			$M->[0][1] * ( $M->[1][0]*$M->[2][2] - $M->[1][2]*$M->[2][0] ) +
			$M->[0][2] * ( $M->[1][0]*$M->[2][1] - $M->[1][1]*$M->[2][0] );
	if ($D == 0)  {
		print STDERR "ERROR ... Inversion of singular matrix!\n";
		exit -1;
	}

	$Minv->[0][0] =	 ($M->[1][1]*$M->[2][2]-$M->[1][2]*$M->[2][1])/$D;
	$Minv->[0][1] = -($M->[0][1]*$M->[2][2]-$M->[0][2]*$M->[2][1])/$D;
	$Minv->[0][2] =	 ($M->[0][1]*$M->[1][2]-$M->[0][2]*$M->[1][1])/$D;
	$Minv->[1][0] = -($M->[1][0]*$M->[2][2]-$M->[2][0]*$M->[1][2])/$D;
	$Minv->[1][1] =	 ($M->[0][0]*$M->[2][2]-$M->[0][2]*$M->[2][0])/$D;
	$Minv->[1][2] = -($M->[0][0]*$M->[1][2]-$M->[0][2]*$M->[1][0])/$D;
	$Minv->[2][0] =	 ($M->[1][0]*$M->[2][1]-$M->[2][0]*$M->[1][1])/$D;
	$Minv->[2][1] = -($M->[0][0]*$M->[2][1]-$M->[0][1]*$M->[2][0])/$D;
	$Minv->[2][2] =	 ($M->[0][0]*$M->[1][1]-$M->[0][1]*$M->[1][0])/$D;

	return $Minv;
}

# matrix determinant
sub mdet {
	my $M = shift;
	my $D = $M->[0][0] * ( $M->[1][1]*$M->[2][2] - $M->[2][1]*$M->[1][2] ) -
			$M->[0][1] * ( $M->[1][0]*$M->[2][2] - $M->[1][2]*$M->[2][0] ) +
			$M->[0][2] * ( $M->[1][0]*$M->[2][1] - $M->[1][1]*$M->[2][0] );
	return $D;
}

# vector norm
sub vnorm {
	my $x = shift;
	return sqrt( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# vector norm^2
sub vnorm2 {
	my $x = shift;
	return ( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# cart distance
sub vdist {
	my ($x1, $x2) = @_;
	return sqrt (
			   ($x1->[0]-$x2->[0])*($x1->[0]-$x2->[0]) +
			   ($x1->[1]-$x2->[1])*($x1->[1]-$x2->[1]) +
			   ($x1->[2]-$x2->[2])*($x1->[2]-$x2->[2])
				);
}

# are two transformations the inverso of one another
sub is_inverse {
	my $tol = 1e-8;
	my ( $R_i,$T_i, $R_j,$T_j ) = @_;

	my $testR = mmult( $R_i , $R_j );
	my $testT = vadd( mapply( $R_j,$T_i ) , $T_j );

	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
			   square($testR->[0][1])	+ square($testR->[0][2])   + square($testR->[1][2]) +
			   square($testR->[1][0])	+ square($testR->[2][0])   + square($testR->[2][1]);
	my $errT = square($testT->[0])+square($testT->[1])+square($testT->[2]);

#print " (0) $errR	 $errT\n";
	if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

#
sub is_identity {
	my $testR = shift;
	my $tol = 1e-8;
	if (scalar( @_ ) >= 1) {
		$tol = shift;
	}
	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
			   square($testR->[0][1])	+ square($testR->[0][2])   + square($testR->[1][2]) +
			   square($testR->[1][0])	+ square($testR->[2][0])   + square($testR->[2][1]);

	if ($errR < $tol) { return 1; }
	return 0;
}

# is the transform (Rn,Tn) equivalent to the transform (Ri,Ti)->(Rj,Tj)
sub is_equivalent {
	my $tol = 1e-8;
	my ( $R_n,$T_n, $R_i,$T_i, $R_j,$T_j ) = @_;

	my $R_i_inv = minv( $R_i );
	my $T_i_inv = [ -$T_i->[0], -$T_i->[1], -$T_i->[2] ];

	my $R_test = mmult( $R_i_inv, $R_j );
	my $T_test = mapply( $R_i_inv, vsub( $T_j, $T_i ) );

	my $errR = square($R_test->[0][0]-$R_n->[0][0])+square($R_test->[0][1]-$R_n->[0][1])+square($R_test->[0][2]-$R_n->[0][2])
			 + square($R_test->[1][0]-$R_n->[1][0])+square($R_test->[1][1]-$R_n->[1][1])+square($R_test->[1][2]-$R_n->[1][2])
			 + square($R_test->[2][0]-$R_n->[2][0])+square($R_test->[2][1]-$R_n->[2][1])+square($R_test->[2][2]-$R_n->[2][2]);
	my $errT = square($T_test->[0]-$T_n->[0])+square($T_test->[1]-$T_n->[1])+square($T_test->[2]-$T_n->[2]);

	##
	## don't think we need to check T ...
	if ($errR < $tol) { return 1; }
	#if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

###########
# f2c/c2f #
###########

sub d2r { return (@_[0]*PI/180); }
sub square { return @_[0]*@_[0]; }

# my ($f2c,$c2f) = crystparams($a,$b,$c,$alpha,$beta,$gamma)
sub crystparams {
	my ($a,$b,$c,$alpha,$beta,$gamma) = @_;

	if ($a*$b*$c == 0) {
		print STDERR "Must provide valid crystal parameters!\n";
		exit -1;
	}

	my $f2c = [ [0,0,0] , [0,0,0] , [0,0,0] ];

	my $ca = cos(d2r($alpha)); my $cb = cos(d2r($beta)); my $cg = cos(d2r($gamma));
	my $sa = sin(d2r($alpha)); my $sb = sin(d2r($beta)); my $sg = sin(d2r($gamma));
	$f2c->[0][0] = $a;
	$f2c->[0][1] = $b * $cg;
	$f2c->[0][2] = $c * $cb;
	$f2c->[1][1] = $b * $sg;
	$f2c->[1][2] = $c * ($ca - $cb*$cg) / $sg;
	$f2c->[2][2] = $c * $sb * sqrt(1.0 - square(($cb*$cg - $ca)/($sb*$sg)));

	my $c2f = minv($f2c);

	return ($f2c,$c2f);
}


###################################################################################
###################################################################################
###
### Crystallographic symmetry groups
###
###################################################################################
###################################################################################

# my ($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp( $spacegpid );
sub spacegp_lookup {
	my $id = shift @_;
	my ($gpid, $nsymm, @Rs, @Ts, @Cs, @R_all, @T_all, $cheshire);

	# strip leading/trailing spaces
	$id =~ s/^\s+//;
	$id =~ s/\s+$//;

	if ($id eq "P 1") {
		$gpid = 1;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,0] , [0,0] ];
	}
	elsif ($id eq "P -1") {
		$gpid = 2;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2 1") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 2 1 1") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21 1") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 21 1 1") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "B 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "C 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 m 1") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 m") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 1 1") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 c 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 n 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 a 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 a") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 n") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 b") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 c 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 n 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 a 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 a 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 n 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 c 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 a") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 n") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 b") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 b") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 n") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 a") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C n 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B n 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 2/m 1") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/m") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/m 1 1") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/m 1") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/m") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/m 1 1") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/c 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/n 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/a 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/a") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/n") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/b") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/b 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/n 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/c 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/c 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/n 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/a 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/a") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/n") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/b") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/b 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/n 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/c 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/c 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/n 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/a 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/a 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/n 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/c 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/a") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/n") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/b") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/b") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/n") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/a") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/b 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/n 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/c 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/c 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/n 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/b 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 2 2") {
		$gpid = 16;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 2 21") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 2 2") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 21 2") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 21 2") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 21 21") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 2 21") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 21 21") {
		$gpid = 19;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2 2 21") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 21 2 2") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2 21 2") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F 2 2 2") {
		$gpid = 22;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2 2 2") {
		$gpid = 23;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 21 21 21") {
		$gpid = 24;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m 2") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 m m") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 2 m") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m c 21") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c m 21") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 m a") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 a m") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 21 m") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 21 b") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c c 2") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 a a") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 2 b") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m a 2") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b m 2") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 m b") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 c m") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 2 m") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 2 a") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c a 21") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b c 21") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 a b") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 c a") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 21 b") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 21 a") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n c 2") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c n 2") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 n a") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 a n") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 2 n") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 2 b") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m n 21") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n m 21") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 m n") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 n m") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 21 m") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 21 n") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b a 2") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 c b") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 2 a") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n a 21") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b n 21") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 n b") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 c n") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 21 n") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 21 a") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n n 2") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 n n") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 2 n") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m m 2") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 2 m m") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 2 m") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m c 21") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c m 21") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 21 m a") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 21 a m") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 21 m") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 21 b") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c c 2") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 2 a a") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 2 b") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m m 2") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m m 2") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 m m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 m m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 2 m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m 2 m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A b m 2") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m a 2") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 c m") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 m b") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 2 a") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A c 2 m") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m a 2") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b m 2") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 m b") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 c m") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 2 m") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m 2 a") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A b a 2") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b a 2") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 c b") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 c b") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 2 a") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A c 2 a") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F m m 2") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F 2 m m") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F m 2 m") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F d d 2") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F 2 d d") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F d 2 d") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m m 2") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 m m") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 2 m") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b a 2") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 c b") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 2 a") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m a 2") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b m 2") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 m b") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 c m") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 2 m") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 2 a") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m m m") {
		$gpid = 47;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n n :1") {
		$gpid = 48;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n n :2") {
		$gpid = 48;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c m") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a a") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m b") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a n :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a n :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c b :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c b :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n a :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n a :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m a") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m b") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n a") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n b") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n a") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m b") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m n") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n m") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c m") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a n") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a m") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c b") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m a") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c n") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a a") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n b") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c m") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a m") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c a") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a b") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m a") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m b") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n m") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n n") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m n") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m n :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m n :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m m :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m m :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n m :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n m :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c n") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a n") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c a") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a b") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n a") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n b") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c a") {
		$gpid = 61;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a b") {
		$gpid = 61;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m a") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n b") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n m") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m n") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c n") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a m") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m c m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c m m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m m a") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m a m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b m m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m m b") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m c a") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c m b") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b m a") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a m") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c m") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m a b") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c m") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m a a") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b m b") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m a") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m b") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b m m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c m m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m c m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m a m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b a a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b a a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b a b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b a b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F m m m") {
		$gpid = 69;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d d d :1") {
		$gpid = 70;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d d d :2") {
		$gpid = 70;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m m") {
		$gpid = 71;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b a m") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m c b") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c m a") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b c a") {
		$gpid = 73;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c a b") {
		$gpid = 73;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m a") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m b") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b m m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c m m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m c m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m a m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4") {
		$gpid = 75;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 41") {
		$gpid = 76;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42") {
		$gpid = 77;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 43") {
		$gpid = 78;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4") {
		$gpid = 79;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41") {
		$gpid = 80;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P -4") {
		$gpid = 81;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4") {
		$gpid = 82;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m") {
		$gpid = 83;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m") {
		$gpid = 84;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n :1") {
		$gpid = 85;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n :2") {
		$gpid = 85;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n :1") {
		$gpid = 86;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n :2") {
		$gpid = 86;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m") {
		$gpid = 87;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a :1") {
		$gpid = 88;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a :2") {
		$gpid = 88;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 2 2") {
		$gpid = 89;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 21 2") {
		$gpid = 90;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 41 2 2") {
		$gpid = 91;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 41 21 2") {
		$gpid = 92;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42 2 2") {
		$gpid = 93;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42 21 2") {
		$gpid = 94;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 43 2 2") {
		$gpid = 95;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 43 21 2") {
		$gpid = 96;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4 2 2") {
		$gpid = 97;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41 2 2") {
		$gpid = 98;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 m m") {
		$gpid = 99;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 b m") {
		$gpid = 100;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 c m") {
		$gpid = 101;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 n m") {
		$gpid = 102;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 c c") {
		$gpid = 103;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 n c") {
		$gpid = 104;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 m c") {
		$gpid = 105;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 b c") {
		$gpid = 106;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4 m m") {
		$gpid = 107;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4 c m") {
		$gpid = 108;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41 m d") {
		$gpid = 109;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41 c d") {
		$gpid = 110;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P -4 2 m") {
		$gpid = 111;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 2 c") {
		$gpid = 112;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 21 m") {
		$gpid = 113;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 21 c") {
		$gpid = 114;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 m 2") {
		$gpid = 115;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 c 2") {
		$gpid = 116;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 b 2") {
		$gpid = 117;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 n 2") {
		$gpid = 118;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 m 2") {
		$gpid = 119;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 c 2") {
		$gpid = 120;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 2 m") {
		$gpid = 121;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 2 d") {
		$gpid = 122;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m m m") {
		$gpid = 123;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m c c") {
		$gpid = 124;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n b m :1") {
		$gpid = 125;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n b m :2") {
		$gpid = 125;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n n c :1") {
		$gpid = 126;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n n c :2") {
		$gpid = 126;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m b m") {
		$gpid = 127;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m n c") {
		$gpid = 128;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n m m :1") {
		$gpid = 129;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n m m :2") {
		$gpid = 129;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n c c :1") {
		$gpid = 130;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n c c :2") {
		$gpid = 130;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m m c") {
		$gpid = 131;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m c m") {
		$gpid = 132;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n b c :1") {
		$gpid = 133;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n b c :2") {
		$gpid = 133;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n n m :1") {
		$gpid = 134;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n n m :2") {
		$gpid = 134;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m b c") {
		$gpid = 135;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m n m") {
		$gpid = 136;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n m c :1") {
		$gpid = 137;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n m c :2") {
		$gpid = 137;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n c m :1") {
		$gpid = 138;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n c m :2") {
		$gpid = 138;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m m m") {
		$gpid = 139;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m c m") {
		$gpid = 140;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a m d :1") {
		$gpid = 141;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a m d :2") {
		$gpid = 141;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a c d :1") {
		$gpid = 142;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.75];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a c d :2") {
		$gpid = 142;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 3") {
		$gpid = 143;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 31") {
		$gpid = 144;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 32") {
		$gpid = 145;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 :H"|| $id eq "H 3") {
		$gpid = 146;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 :R") {
		$gpid = 146;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P -3") {
		$gpid = 147;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 :H") {
		$gpid = 148;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 :R") {
		$gpid = 148;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 3 1 2") {
		$gpid = 149;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 3 2 1") {
		$gpid = 150;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 31 1 2") {
		$gpid = 151;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 31 2 1") {
		$gpid = 152;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 32 1 2") {
		$gpid = 153;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 32 2 1") {
		$gpid = 154;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R 3 2 :H" || $id eq "H 3 2") {
		$gpid = 155;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R 3 2 :R") {
		$gpid = 155;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 3 m 1") {
		$gpid = 156;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 3 1 m") {
		$gpid = 157;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 3 c 1") {
		$gpid = 158;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 3 1 c") {
		$gpid = 159;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "R 3 m :H") {
		$gpid = 160;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 m :R") {
		$gpid = 160;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 c :H") {
		$gpid = 161;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 c :R") {
		$gpid = 161;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P -3 1 m") {
		$gpid = 162;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 1 c") {
		$gpid = 163;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 m 1") {
		$gpid = 164;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 c 1") {
		$gpid = 165;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 m :H") {
		$gpid = 166;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 m :R") {
		$gpid = 166;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 c :H") {
		$gpid = 167;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 c :R") {
		$gpid = 167;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6") {
		$gpid = 168;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 61") {
		$gpid = 169;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 65") {
		$gpid = 170;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 62") {
		$gpid = 171;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 64") {
		$gpid = 172;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63") {
		$gpid = 173;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P -6") {
		$gpid = 174;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m") {
		$gpid = 175;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m") {
		$gpid = 176;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6 2 2") {
		$gpid = 177;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 61 2 2") {
		$gpid = 178;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 65 2 2") {
		$gpid = 179;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 62 2 2") {
		$gpid = 180;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 64 2 2") {
		$gpid = 181;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 63 2 2") {
		$gpid = 182;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6 m m") {
		$gpid = 183;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 6 c c") {
		$gpid = 184;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63 c m") {
		$gpid = 185;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63 m c") {
		$gpid = 186;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P -6 m 2") {
		$gpid = 187;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 c 2") {
		$gpid = 188;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 2 m") {
		$gpid = 189;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 2 c") {
		$gpid = 190;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m m m") {
		$gpid = 191;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m c c") {
		$gpid = 192;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m c m") {
		$gpid = 193;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m m c") {
		$gpid = 194;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 2 3") {
		$gpid = 195;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F 2 3") {
		$gpid = 196;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2 3") {
		$gpid = 197;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 21 3") {
		$gpid = 198;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 21 3") {
		$gpid = 199;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3") {
		$gpid = 200;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 :1") {
		$gpid = 201;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 :2") {
		$gpid = 201;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F m -3") {
		$gpid = 202;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 :1") {
		$gpid = 203;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 :2") {
		$gpid = 203;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m -3") {
		$gpid = 204;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P a -3") {
		$gpid = 205;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I a -3") {
		$gpid = 206;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 4 3 2") {
		$gpid = 207;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 42 3 2") {
		$gpid = 208;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F 4 3 2") {
		$gpid = 209;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F 41 3 2") {
		$gpid = 210;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4 3 2") {
		$gpid = 211;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 43 3 2") {
		$gpid = 212;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 41 3 2") {
		$gpid = 213;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 41 3 2") {
		$gpid = 214;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P -4 3 m") {
		$gpid = 215;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F -4 3 m") {
		$gpid = 216;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 3 m") {
		$gpid = 217;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P -4 3 n") {
		$gpid = 218;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F -4 3 c") {
		$gpid = 219;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 3 d") {
		$gpid = 220;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3 m") {
		$gpid = 221;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 n :1") {
		$gpid = 222;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 n :2") {
		$gpid = 222;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3 n") {
		$gpid = 223;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 m :1") {
		$gpid = 224;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 m :2") {
		$gpid = 224;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F m -3 m") {
		$gpid = 225;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F m -3 c") {
		$gpid = 226;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 m :1") {
		$gpid = 227;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 m :2") {
		$gpid = 227;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 c :1") {
		$gpid = 228;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 c :2") {
		$gpid = 228;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.5,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.5,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0,0.25];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m -3 m") {
		$gpid = 229;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I a -3 d") {
		$gpid = 230;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 1 21 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "C 1 21 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "B 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	else {
		print STDERR "SPACEGROUP $id undefined!\n";
		exit -1;
	}

	# expand cenops
	foreach my $C (@Cs) {
		foreach my $i (0..$#Ts) {
			push @R_all, $Rs[$i];
			push @T_all, vadd( $Ts[ $i ] , $C );
		}
	}
	my $nsymm = $#T_all + 1;

	return ($gpid, $nsymm, \@R_all, \@T_all, $cheshire);
}

#######
# end #
#######
