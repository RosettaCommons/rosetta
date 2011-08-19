#!/usr/bin/perl
##
##
###############################################################################

use strict;
use Math::Trig;   # inv trig ops
use POSIX qw(ceil floor fmod fabs);
#use Getopt::Long qw(permute);
use constant PI    => 4 * atan2(1, 1);

use lib (".");
use File::Basename;
use lib dirname(__FILE__);

## user modules
require "kabsch.pm";
require "matrix.pm";

###############################################################################

if ($#ARGV < 0) {
	print STDERR "usage: $0 [options]\n";
	print STDERR "example:   $0 -a A -b B -i G -r 12.0 -p mystructure.pdb\n";
	print STDERR "options: \n";
	print STDERR "    -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
	print STDERR "    -r <real>   : [default 8.0] the max CA-CA distance between two interacting chains\n";
	print STDERR "    -a <char>   : [default A] the chain ID of the main chain\n";
	print STDERR "    -b <char>   : [default B] the chain ID of the next chain along the fiber\n";
	print STDERR "    -c <char>   : the chain ID of the chain perpindicular to the -b chain (2D symmetry case)\n";
	print STDERR "    -i <char>   : the chain ID of a chain in the (non-helical) symmetric complex\n";
	print STDERR "    -t <real>   : [default 2] the number of helical turns to generate along the -b direction\n";
	print STDERR "    -u <real>    : [default 2] the number of repeats to generate along the -c direction\n";
	exit -1;
}

my $pdbfile;
my $interact_dist = 8.0;  # min interaction distance
my $primary_chain = 'A';
my $helical_chain = 'B';
my $ncs_chain = '';
my $nturns = 1;

my $perp_chain = '';
my $nperp_repeats = 1;


## parse options (do this by hand since Getopt does not handle this well)
my @suboptions = split /(-[a-z|A-Z] )/, (join ' ',@ARGV);
for ( my $i=0; $i<=$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$interact_dist = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-t " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$nturns = ( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-u " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$nperp_repeats = ( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-a " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$primary_chain = $suboptions[++$i];
		$primary_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-b " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$helical_chain = $suboptions[++$i];
		$helical_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-c " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$perp_chain = $suboptions[++$i];
		$perp_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-i " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$ncs_chain = $suboptions[++$i];
		$ncs_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
		$pdbfile =~ s/\s*(\S+)\s*/$1/;
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

###
### Input PDB file
my %chains;
my @filebuf;
open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
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
		}
	}
}
close (PDB);

# recenter the primary chain
if ( ! defined $chains{ $primary_chain } ) {
	die "Chain '$primary_chain' not in input!\n";
}

my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
my $COM_0 = recenter( $chains{ $primary_chain } );

# find residue closest to CoM of the system
my $minDist2 = 9999999;
my $minRes = 1;
foreach my $i ( 0..scalar( @{ $chains{ $primary_chain } })-1 ) {
	my $dist2 = vnorm2( $chains{ $primary_chain }->[$i] );
	if ($dist2 < $minDist2) {
		$minDist2 = $dist2;
		$minRes = $i+1;
	}
}


###############################
###############################
##
## first expand helical symm

my @helical_chain_split = split( ':', $helical_chain );
$helical_chain = $helical_chain_split[0];
my $force_symm_order = 0;
# optionally ... allow input to 'force' a symmetric order
# NOTE THAT THIS MAY RESULT IS A SYSTEM QUITE FAR FROM THE INPUT SYSTEM
if ($#helical_chain_split > 0) {
	$force_symm_order = $helical_chain_split[1];
}


if ( ! defined $chains{ $helical_chain } ) {
	die "Chain $helical_chain not in input!\n";
}

## CA check
if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $helical_chain } } ) ) {
	print STDERR "ERROR! chains '$primary_chain' and '$helical_chain' have different residue counts! (".
	             scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $helical_chain } } ).")\n";
	die "Chain length mismatch!\n";
}

# get superposition
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

	# if we have 2D lattice symmetry only allow 1 or 2 subunits per turn
	if ( $perp_chain ne '' ) {
		my $omega = acos($W);
		$subunits_per_turn = int(PI/$omega + 0.5);

		if ($subunits_per_turn == 1) {
			($X,$Y,$Z,$W) = ($del_COM->[0],$del_COM->[1],$del_COM->[2],1);
			$R = [ [1,0,0] , [0,1,0] , [0,0,1] ];
		} elsif ($subunits_per_turn == 2) {
			$W = cos( PI/$subunits_per_turn );
			my $newW = -$Wmult * $W;
			my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
			$R = quat2R( $X*$S , $Y*$S, $Z*$S, $newW );
		} else {
			print STDERR "If lattice symm op defined, helix must have 1 or 2 subunits per turn (found ".$subunits_per_turn.")\n";
		}
	} else {
		my $omega = acos($W);
		$subunits_per_turn = (PI/$omega);
		my $newW = -$Wmult * $W;
		my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
		$R = quat2R( $X*$S , $Y*$S, $Z*$S, $newW );
	}

}


# project $del_COM to helical axis
my $helical_axis = [$X,$Y,$Z];
normalize( $helical_axis );
my $del_COM_inplane    = vsub( $del_COM , vscale(dot($del_COM,$helical_axis),$helical_axis) );
my $del_COM_alonghelix = vsub( $del_COM , $del_COM_inplane );

#print STDERR "W_orig = $Worig\nW = $W\nomega = $omega\n";
print STDERR "Found helical symmetry at chain ".$helical_chain."\n";
print STDERR "   subunits per turn = ".$subunits_per_turn."\n";
print STDERR "   rise  = ".vnorm($del_COM_alonghelix)."\n";
#print STDERR "   inplane delta = ".vnorm($del_COM_inplane)."\n";

# get the center of the helix
my $omega = acos($W);
if ($omega < 1e-6) { $omega = pi; }
my $helix_R = vnorm($del_COM_inplane) / (2*sin( $omega ));
my $helix_center;

my $COM_0_5 = vadd( $COM_0 , vscale( 0.5, $del_COM_inplane ) ); # 1/2way between subunits
my $d_0_5_center = $helix_R * cos( $omega );   # distance from here to helical center

if ($d_0_5_center == 0) {
	$helix_center = $COM_0_5;
} else {
	# direction from here to center
	my $center_dir = vscale ( $Wmult , cross( $del_COM_inplane , $del_COM_alonghelix ) );
	normalize( $center_dir );
	$helix_center = vadd( $COM_0_5 , vscale( $d_0_5_center, $center_dir ) );
}

#print STDERR "   radius (to CoM) = ".$helix_R."\n";
#print STDERR "   d_0.5 (to CoM) = ".$d_0_5_center."\n";
#print STDERR "   helix_center = ".$helix_center->[0]." , ".$helix_center->[1]." , ".$helix_center->[2]."\n";


###############################
###############################
##
## expand perpendicular (2D lattice) symmetry
## CA check
my $del_COM_perp = [0,0,0];
if ( $perp_chain ne '' ) {
	if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $perp_chain } } ) ) {
		print STDERR "ERROR! chains '$primary_chain' and '$helical_chain' have different residue counts! (".
					 scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $perp_chain } } ).")\n";
		die "Chain length mismatch!\n";
	}

	# get superposition
	my ($R_perp,$rmsd_perp, $COM_i_perp, $COM_ij_perp) = rms_align( $chains{ $primary_chain } , $chains{ $perp_chain } );
	$del_COM_perp = vsub ($COM_i_perp, $COM_0);

	# this must always be a translation only
	if ( !is_identity( $R_perp , 0.01 ) ) {
		print STDERR "Trasformation to chain $perp_chain must be a translation only!\n";
		print STDERR "   >> R = [".$R->[0][0].",".$R->[0][1].",".$R->[0][2]." ; ".
		                           $R->[1][0].",".$R->[1][1].",".$R->[1][2]." ; ".
		                           $R->[2][0].",".$R->[2][1].",".$R->[2][2]." ]\n";
		exit(1);
	}

	# project $del_COM_perp normal to $helical_axis (already normalized)
	my $del_COM_perp_fix = vsub( $del_COM_perp , vscale(dot($del_COM_perp,$helical_axis),$helical_axis) );

	#print STDERR "W_orig = $Worig\nW = $W\nomega = $omega\n";
	print STDERR "Found perpendicular symmetry at chain ".$perp_chain."\n";
	print STDERR "   shift    = ".vnorm($del_COM_perp_fix)."\n";
	print STDERR "   [error1] = ".vnorm(vsub($del_COM_perp,$del_COM_perp_fix))."\n";

	$del_COM_perp = $del_COM_perp_fix;
} else {
	$nperp_repeats = 0;
}


##
##
# next expand point symmetry (Cn point group as each helical subunit)
# TO DO: nonpolar helical symmetry! (that is, Dn point group as each helical subunit)
#
my ($sym_order_ncs, $R_ncs, $T_ncs, $global_T) = (1,[[1,0,0],[0,1,0],[0,0,1]], [0,0,0], [0,0,0]);
if ($ncs_chain ne '' ) {
	if ( ! defined $chains{ $ncs_chain } ) {
		die "Chain $ncs_chain not in input!\n";
	}

	## CA check
	if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $ncs_chain } } ) ) {
		print STDERR "ERROR! chains '$primary_chain' and '$helical_chain' have different residue counts! (".
					 scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $ncs_chain } } ).")\n";
		die "Chain length mismatch!\n";
	}

	# get superposition
	my ($rmsd_ncs, $COM_i_ncs, $COM_ij_ncs);
	($R_ncs,$rmsd_ncs, $COM_i_ncs, $COM_ij_ncs) = rms_align( $chains{ $primary_chain } , $chains{ $ncs_chain } );
	my $del_COM_ncs = vsub ($COM_i_ncs, $COM_0);

	my ($X_ncs,$Y_ncs,$Z_ncs,$W_ncs)=R2quat($R_ncs);
	my $Worig_ncs = $W_ncs;
	my $Wmult_ncs = 1;
	if ($W_ncs < 0) { $W_ncs = -$W_ncs; $Wmult_ncs = -1; }
	my $omega_ncs = acos($W_ncs);
	$sym_order_ncs = int(PI/$omega_ncs + 0.5);
	print STDERR "Found ".$sym_order_ncs."-fold (".(PI/$omega_ncs).") symmetric complex at chain ".$ncs_chain."\n";

	# if we have a 2D lattice only allow C1/C2 symmetry 
	if ( $perp_chain ne '' && $sym_order_ncs > 2) {
		print STDERR "If lattice symm op defined, pt symm must have 1 or 2 subunits (found ".$sym_order_ncs.")\n";
		exit (1);
	}

	# make perfectly symmetrical version of superposition
	# adjust rotation axis to be parallel to helical axis
	my $newW = -$Wmult *cos( PI/$sym_order_ncs );
	my $newS = sqrt ( (1-$newW*$newW) );
	my $newQ = [ $helical_axis->[0]*$newS , $helical_axis->[1]*$newS, $helical_axis->[2]*$newS, $newW];
	$R_ncs = quat2R( $newQ->[0], $newQ->[1], $newQ->[2], $newQ->[3] );

	# get transform
	# symmetrize it
	my $err_pos = [0,0,0];
	my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $i (1..$sym_order_ncs) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$del_COM_ncs ) );
		$R_i = mmult($R_ncs, $R_i);
	}
	print STDERR "NCS Translation error = ".vnorm( $err_pos )."\n";

	# projection to symm plane
	$T_ncs = vsub( $del_COM_ncs , [ $err_pos->[0]/$sym_order_ncs , $err_pos->[1]/$sym_order_ncs , $err_pos->[2]/$sym_order_ncs ] );

	# make sure center of NCS is on helical axis???
	# find center of symm complex (in coord frame of subunit 1)
	my $CoM_cplx = [0,0,0];
	$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $i (1..$sym_order_ncs-1) {
		$CoM_cplx = vadd( $CoM_cplx, vadd( $CoM_cplx, mapply( $R_i,$T_ncs ) ) );
		$R_i = mmult($R_ncs, $R_i);
	}
	$CoM_cplx = vscale (1/$sym_order_ncs, $CoM_cplx) ;
	$CoM_cplx = vadd ($COM_0, $CoM_cplx) ;

	# transform to helical Center
	$global_T = vsub ( $helix_center, $CoM_cplx );
	#print "global_T = ".$global_T->[0].",".$global_T->[1].",".$global_T->[2]."\n";
}

print STDERR "Radius (to CoM) = ".vnorm( vsub( $helix_center, vadd( $global_T, $COM_0) ) )."\n";


##########################
##########################
my $nsubunits_to_gen = ceil( $nturns*PI/$omega );

# store transforms
my $Rs = {};
my $Ts = {};

my $T_sec;
my $R_helix;
my $T_helix;
my $Rinv_helix;
my $Tinv_helix;

foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	$T_sec = vscale( $sec_shift , $del_COM_perp );

	#transform to helical frame (+ direction)
	$R_helix = [[1,0,0],[0,1,0],[0,0,1]];
	$T_helix = [0,0,0];
	
	#transfrom to helical frame (- direction)
	$Rinv_helix = [[1,0,0],[0,1,0],[0,0,1]];
	$Tinv_helix = [0,0,0];

	###
	my $Rinv = minv($R);
	my $T = $del_COM_inplane;    # points from n->(n+1)
	my $Tinv = vscale( -1, mapply( $Rinv, $T ));  # points from n->(n-1)
		
	
	#
	print STDERR "Generating [-".$nsubunits_to_gen." to +".$nsubunits_to_gen."] at offset $sec_shift\n";
	foreach my $subunit (0 .. $nsubunits_to_gen) {
		#my $T_about    = vadd( $helix_center , vscale( $subunit, $del_COM_alonghelix ) );
		#my $Tinv_about = vadd( $helix_center , vscale(-$subunit, $del_COM_alonghelix ) );
	
		# gen in + direction
		my $R_i = [[1,0,0],[0,1,0],[0,0,1]];
		my $T_i = [0,0,0];
		#my $T_i = $T_sec;
		foreach my $i (0..$sym_order_ncs-1) {
			# rotate about helical axis then translate to the new CoM
			my $R_helix_i = mmult( $R_i, $R_helix  );
			my $T_helix_i = vadd( $T_helix, mapply( $R_helix, $T_i) );
	
			my $R_global_T = mapply( $R_helix, $global_T );
	
			#$Rs->{ $subunit."_".$i } = $R_helix_i;
			#$Ts->{ $subunit."_".$i } = vadd( $R_global_T, vadd( $COM_0, vadd( $T_helix_i, vscale( $subunit, $del_COM_alonghelix ) ) ) );
			$Rs->{ $sec_shift."_".$subunit."_".$i } = $R_helix_i;
			$Ts->{ $sec_shift."_".$subunit."_".$i } = 
				vadd( $T_sec, vadd( $R_global_T, vadd( $COM_0, vadd( $T_helix_i, vscale( $subunit, $del_COM_alonghelix ) ) ) ) );

			$T_i = vadd( $T_i , mapply( $R_i, $T_ncs ) );
			$R_i = mmult( $R_ncs , $R_i );
		}
	
		# gen in - direction
		if ($subunit != 0) {
			$R_i = [[1,0,0],[0,1,0],[0,0,1]];
			my $T_i = [0,0,0];
			#my $T_i = $T_sec;
			foreach my $i (0..$sym_order_ncs-1) {
				# rotate about helical axis then translate to the new CoM
				my $R_helix_i = mmult( $R_i, $Rinv_helix );
				my $T_helix_i = vadd( $Tinv_helix, mapply( $Rinv_helix, $T_i) );
	
				my $R_global_T = mapply( $Rinv_helix, $global_T );
	
				#$Rs->{ -$subunit."_".$i } = $R_helix_i;
				#$Ts->{ -$subunit."_".$i } = vadd( $R_global_T, vadd( $COM_0, vadd( $T_helix_i, vscale( -$subunit, $del_COM_alonghelix ) ) ) );
				$Rs->{ $sec_shift."_".-$subunit."_".$i } = $R_helix_i;
				$Ts->{ $sec_shift."_".-$subunit."_".$i } = 
					vadd( $T_sec, vadd( $R_global_T, vadd( $COM_0, vadd( $T_helix_i, vscale( -$subunit, $del_COM_alonghelix ) ) ) ) );
	
				$T_i = vadd( $T_i , mapply( $R_i, $T_ncs ) );
				$R_i = mmult( $R_ncs , $R_i );
			}
		}
	
		# update the transform to the 1st subunit in next helical layer
		$T_helix = vadd( $T_helix, mapply( $R_helix, $T ) );
		$Tinv_helix = vadd( $Tinv_helix, mapply( $Rinv_helix, $Tinv ) );
		$R_helix = mmult( $R, $R_helix );
		$Rinv_helix = mmult( $Rinv, $Rinv_helix );
	}

}


#############################
#############################
my %symminterface = ( '0_0_0'=>1 );
my $counter = 2;
foreach my $X_i ( @{ $chains{ $primary_chain } } ) {
  	foreach my $Y_i (  @{ $chains{ $primary_chain } } ) {
		foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
			foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
				foreach my $i (0..$sym_order_ncs-1) {
					my $id = $sec_shift."_".$subunit."_".$i;
					next if (defined $symminterface{ $id });
	
					#   x_i = R_i * (x_0 - COM_0) + COM_i
					#   The rms function already ofsets x_0 by -COM_0
					my $rX_i = vadd( $X_i , $COM_0 );
					#my $rX_i = $X_i;
					my $rY_j = vadd( mapply($Rs->{ $id }, $Y_i) , $Ts->{ $id } );
					my $delXY = vsub( $rY_j,$rX_i );
					my $dist2XY = vnorm2( $delXY );
	
					if ($dist2XY < $interact_dist*$interact_dist) {
						# we have a hit! tag NCS copy $j_symm as a non-symmetic interface
						print STDERR " Adding interface '".$id."'\n";
						$symminterface{ $id } = $counter++;
					}
				}
			}
		}
	}
}


# 
#######################################
## symm file gen
## write output symm file
# symmetry_name c4
# E = 2*VRT2
# anchor_residue 17
my $symmname = $pdbfile;
$symmname =~ s/\.pdb$//;
$symmname = $symmname."_helix_C".$sym_order_ncs;
print "symmetry_name $symmname\n";
print "E = ".($sym_order_ncs)."*VRT_0_0_0_base";
foreach my $complex (keys %symminterface) {
	my ($sec_shift,$subunit,$i) = split '_', $complex;

	# ($subunit , $i) symmetric with (-$subunit , -$i)
	next if ($sec_shift < 0);
	next if ($sec_shift == 0 && $subunit < 0);
	next if ($sec_shift == 0 && $subunit == 0 && $i > $sym_order_ncs/2);
	next if ($sec_shift == 0 && $subunit == 0 && $i == 0);

	my $cplx_string = $complex;
	$cplx_string =~ s/_-(\d)/_n\1/g;

	#if ($subunit == 0 && $i == $sym_order_ncs/2) {
		print " + 1*(VRT_0_0_0_base:VRT_".$cplx_string."_base)";
	#} else {
	#	print " + 2*(VRT_0_0_0_base:VRT_".$cplx_string."_base)";
	#}
}
print "\n";
print "anchor_residue $minRes\n";

# virtual_coordinates_start
# xyz VRT1 -1,0,0 0,1,0 0,0,0
# xyz VRT2 0,-1,0 -1,0,0 0,0,0
# xyz VRT3 1,0,0 0,-1,0 0,0,0
# xyz VRT4 0,1,0 1,0,0 0,0,0
# virtual_coordinates_stop
my @fakepdblines = ();
print "virtual_coordinates_start\n";

my $xyzline = sprintf("xyz VRT_0 %.3f,%.3f,%.3f %.3f,%.3f,%.3f %.3f,%.3f,%.3f",
                       1.0,0.0,0.0,  0.0,1.0,0.0,
                       $helix_center->[0],$helix_center->[1],$helix_center->[2] );
print "$xyzline\n";


foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {

	my $helix_shifted    = vadd( $helix_center , vscale( $sec_shift , $del_COM_perp ) );

	#####
	# controlling vrts along 2ary axis
	#####
	if ($nperp_repeats > 0) {
		my $subunit = -$nsubunits_to_gen;
		my $i = 0;

		my $T_about    = vadd( $helix_shifted , vscale( $subunit, $del_COM_alonghelix ) );
	
		my $id = $sec_shift."_".$subunit."_".$i;
	
		my $xyzline = "xyz VRT_intra_".$id;
		$xyzline =~ s/_-(\d)/_n\1/g;

		# X points to along 2ary axis
		my $myX = deep_copy($del_COM_perp);
		normalize( $myX );

		# Y points up helical axis
		my $myY = deep_copy($del_COM_alonghelix);
		normalize( $myY );
		my $string = sprintf("%.3f,%.3f,%.3f", $myX->[0], $myX->[1], $myX->[2]);
		$xyzline = $xyzline." ".$string;
		$string = sprintf("%.3f,%.3f,%.3f", $myY->[0], $myY->[1], $myY->[2]);
		$xyzline = $xyzline." ".$string;

		# orig
		my $origin  = $T_about;
		$string = sprintf("%.3f,%.3f,%.3f", $origin->[0], $origin->[1], $origin->[2]);
		$xyzline = $xyzline." ".$string;
		print "$xyzline\n";
		
		my $fakePDBline1 = sprintf "ATOM    %3d  C   ORI Z   1     %7.3f %7.3f %7.3f  1.00  0.00\n", 1,
							 $origin->[0], $origin->[1], $origin->[2];
		my $fakePDBline2 = sprintf "ATOM    %3d  O   X   Z   2     %7.3f %7.3f %7.3f  1.00  0.00\n", 2,
							 $origin->[0]+$myX->[0], $origin->[1]+$myX->[1], $origin->[2]+$myX->[2];
		my $fakePDBline3 = sprintf "ATOM    %3d  C   Y   Z   3     %7.3f %7.3f %7.3f  1.00  0.00\n", 3,
							 $origin->[0]+$myY->[0], $origin->[1]+$myY->[1], $origin->[2]+$myY->[2];
		push @fakepdblines, $fakePDBline1;
		push @fakepdblines, $fakePDBline2;
		push @fakepdblines, $fakePDBline3;
	}

	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		#####
		# controlling vrts on the helical axis
		#####
		my $T_about    = vadd( $helix_shifted , vscale( $subunit, $del_COM_alonghelix ) );
	
		foreach my $i (0..$sym_order_ncs-1) {
			my $id = $sec_shift."_".$subunit."_".$i;
	
			my $xyzline = "xyz VRT_".$id;
			$xyzline =~ s/_-(\d)/_n\1/g;


			# X --> points towards the subunit
			#my $myX = vsub( $Ts->{ $id } , $T_about );
			my $myX = vsub( $T_about,  $Ts->{ $id } );
			if ( vnorm($myX) < 1e-6 ) {
				# just pick any direction perpendicular to helical axis
				if ( $del_COM_alonghelix->[1] == 0 && $del_COM_alonghelix->[2] == 0) {
					$myX = [0,1,0];
				} else {
					$myX = [ 0, -$del_COM_alonghelix->[2], $del_COM_alonghelix->[1]];
				}
			}
			normalize( $myX );

			my $string = sprintf("%.3f,%.3f,%.3f", $myX->[0], $myX->[1], $myX->[2]);
			$xyzline = $xyzline." ".$string;
			# Y --> Z points up helical axis
			my $myZ = $del_COM_alonghelix;
			my $myY = cross( $myZ, $myX );
			normalize( $myY );
			$string = sprintf("%.3f,%.3f,%.3f", $myY->[0], $myY->[1], $myY->[2]);
			$xyzline = $xyzline." ".$string;
			# orig
			my $origin  = $T_about;
			$string = sprintf("%.3f,%.3f,%.3f", $origin->[0], $origin->[1], $origin->[2]);
			$xyzline = $xyzline." ".$string;
			print "$xyzline\n";
		

			my $fakePDBline1 = sprintf "ATOM    %3d  C   ORI Z   1     %7.3f %7.3f %7.3f  1.00  0.00\n", 1,
								 $origin->[0], $origin->[1], $origin->[2];
			my $fakePDBline2 = sprintf "ATOM    %3d  O   X   Z   2     %7.3f %7.3f %7.3f  1.00  0.00\n", 2,
								 $origin->[0]+$myX->[0], $origin->[1]+$myX->[1], $origin->[2]+$myX->[2];
			my $fakePDBline3 = sprintf "ATOM    %3d  C   Y   Z   3     %7.3f %7.3f %7.3f  1.00  0.00\n", 3,
								 $origin->[0]+$myY->[0], $origin->[1]+$myY->[1], $origin->[2]+$myY->[2];
			push @fakepdblines, $fakePDBline1;
			push @fakepdblines, $fakePDBline2;
			push @fakepdblines, $fakePDBline3;
		}
	
	
		#####
		# COM vrts
		#####
		foreach my $i (0..$sym_order_ncs-1) {
			my $id = $sec_shift."_".$subunit."_".$i;
	
			my $xyzline = "xyz VRT_".$id."_base";
			$xyzline =~ s/_-(\d)/_n\1/g;
	
			# X
			my $myX  = mapply( $Rs->{ $id } , [1,0,0] );
			my $string = sprintf("%.3f,%.3f,%.3f", $myX->[0], $myX->[1], $myX->[2]);
			$xyzline = $xyzline." ".$string;
			# Y
			my $myY  = mapply( $Rs->{ $id } , [0,1,0] );
			$string = sprintf("%.3f,%.3f,%.3f", $myY->[0], $myY->[1], $myY->[2]);
			$xyzline = $xyzline." ".$string;
			# orig
			my $origin  = $T_about;
			$string = sprintf("%.3f,%.3f,%.3f", $origin->[0], $origin->[1], $origin->[2]);
			$xyzline = $xyzline." ".$string;
			print "$xyzline\n";
	
			my $fakePDBline1 = sprintf "ATOM    %3d  C   ORI Z   1     %7.3f %7.3f %7.3f  1.00  0.00\n", 1,
								 $origin->[0], $origin->[1], $origin->[2];
			my $fakePDBline2 = sprintf "ATOM    %3d  O   X   Z   2     %7.3f %7.3f %7.3f  1.00  0.00\n", 2,
								 $origin->[0]+$myX->[0], $origin->[1]+$myX->[1], $origin->[2]+$myX->[2];
			my $fakePDBline3 = sprintf "ATOM    %3d  C   Y   Z   3     %7.3f %7.3f %7.3f  1.00  0.00\n", 3,
								 $origin->[0]+$myY->[0], $origin->[1]+$myY->[1], $origin->[2]+$myY->[2];
			push @fakepdblines, $fakePDBline1;
			push @fakepdblines, $fakePDBline2;
			push @fakepdblines, $fakePDBline3;
		}
	}
}
print "virtual_coordinates_stop\n";

## connect_virtual JUMP1 VRT1 VRT2
## connect_virtual JUMP2 VRT1 VRT3
## connect_virtual JUMP3 VRT1 VRT4
# (1) connect bottoms of helices
my $subunit = -$nsubunits_to_gen;
if ($nperp_repeats > 0) {
	foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
		my $id1 = ($sec_shift-1)."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
		my $id2 = $sec_shift."_".$subunit."_0"; $id2 =~ s/-(\d)/n\1/g;
		if ($sec_shift == -$nperp_repeats) {
			print "connect_virtual JUMP_0 VRT_0 VRT_intra_$id2\n";
		} else {
			print "connect_virtual JUMP_intra_$id1 VRT_intra_$id1 VRT_intra_$id2\n";
		}
		print "connect_virtual JUMP_to_helix_$id1 VRT_intra_$id2 VRT_$id2\n";
	}
} else {
	my $id2 = "0_".$subunit."_0"; $id2 =~ s/-(\d)/n\1/g;
	print "connect_virtual JUMP_0 VRT_0 VRT_$id2\n";
}

foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (-$nsubunits_to_gen+1 .. $nsubunits_to_gen) {
		my $id1 = $sec_shift."_".($subunit-1)."_0"; $id1 =~ s/-(\d)/n\1/g;
		my $id2 = $sec_shift."_".$subunit."_0"; $id2 =~ s/-(\d)/n\1/g;
		print "connect_virtual JUMP_$id1 VRT_$id1 VRT_$id2\n";
	}
}

# if point symm, jump to each of pt symm groups
foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		my $id1 = $sec_shift."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
		foreach my $i (1..$sym_order_ncs-1) {
			my $id2 = $sec_shift."_".$subunit."_".$i; $id2 =~ s/-(\d)/n\1/g;
			print "connect_virtual JUMP_pt_$id2 VRT_$id1 VRT_$id2\n";
		}
	}
}

#jump from helical axis to com
foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		foreach my $i (0..$sym_order_ncs-1) {
			my $id = $sec_shift."_".$subunit."_".$i; $id =~ s/-(\d)/n\1/g;
			print "connect_virtual JUMP_".$id."_to_com VRT_".$id." VRT_".$id."_base\n";
		}
	}
}

#jump from com to subunit
foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (-$nsubunits_to_gen .. $nsubunits_to_gen) {
		foreach my $i (0..$sym_order_ncs) {
			my $id = $sec_shift."_".$subunit."_".$i;
			if (defined $symminterface{ $id }) {
				$id =~ s/-(\d)/n\1/g;
				print "connect_virtual JUMP_".$id."_to_subunit VRT_".$id."_base SUBUNIT\n";
			}
		}
	}
}

##################
#  DOFs
##################
# jumps between helical bases
my $sec_shift = 0; #-$nperp_repeats;
my $subunit = -$nsubunits_to_gen;
my $id1 = $sec_shift."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
if ($nperp_repeats > 0) {
	print "set_dof JUMP_intra_$id1 x(".-vnorm($del_COM_perp).")\n";
}

# jumps up the helical axis
$sec_shift = 0;
$subunit = 0;
$id1 = $sec_shift."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
if ($nperp_repeats > 0) {
	print "set_dof JUMP_$id1 z(".-vnorm($del_COM_alonghelix).")\n";
} else {
	print "set_dof JUMP_$id1 z(".-vnorm($del_COM_alonghelix).") angle_z\n";
}

#jump from helical axis to com (not in the case of fiber symmetry)
my $distX = vnorm( vsub( $helix_center, vadd( $global_T, $COM_0) ) );
if (! (is_identity($R) && $sym_order_ncs==1) ) {
#	print "set_dof JUMP_$id1"."_to_com x(".-$distX.")\n";
	print "set_dof JUMP_$id1"."_to_com x(0)\n";
}

# jump about CoM
print "set_dof JUMP_0_0_0_to_subunit angle_x angle_y angle_z\n";

#define jumpgroups
if ($nperp_repeats > 0) {
	print "set_jump_group JUMPGROUP0";
	foreach my $sec_shift ( -$nperp_repeats+1..$nperp_repeats ) {
		my $id1 = ($sec_shift-1)."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
		print "  JUMP_intra_$id1";

		# weight
		if ($sec_shift > 0) { print ":".($sec_shift); }
	}
	print "\n";
}

# jumps up the helical axis
print "set_jump_group JUMPGROUP1";
foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (0 .. $nsubunits_to_gen) {
		my $id1 = $sec_shift."_".$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
		if ($subunit != $nsubunits_to_gen) {
			print "  JUMP_$id1";
			# weight
			if ($sec_shift >= 0) { print ":".($subunit+1); }
		}
	
		if ($subunit != 0) {
			$id1 = $sec_shift."_".-$subunit."_0"; $id1 =~ s/-(\d)/n\1/g;
			print "  JUMP_$id1";
		}
	}
}
print "\n";

#jump from helical axis to coms
print "set_jump_group JUMPGROUP2";
foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
	foreach my $subunit (0 .. $nsubunits_to_gen) {
		foreach my $i (0..$sym_order_ncs-1) {
			my $id = $sec_shift."_".$subunit."_".$i; $id =~ s/-(\d)/n\1/g;
			print " JUMP_".$id."_to_com";
	
			if ($subunit != 0) {
				$id = $sec_shift."_".-$subunit."_".$i; $id =~ s/-(\d)/n\1/g;
				print "  JUMP_".$id."_to_com";
			}
		}
	}
}
print "\n";

#jump from com to subunits
print "set_jump_group JUMPGROUP3";
#foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {
my @repeat_order = (0..$nperp_repeats,-$nperp_repeats..-1);
foreach my $sec_shift ( @repeat_order ) {
	foreach my $subunit (0 .. $nsubunits_to_gen) {
		foreach my $i (0..$sym_order_ncs) {
			my $id = $sec_shift."_".$subunit."_".$i;
			if (defined $symminterface{ $id }) {
				$id =~ s/-(\d)/n\1/g;
				print " JUMP_".$id."_to_subunit";
			}
	
			if ($subunit != 0) {
				$id = $sec_shift."_".-$subunit."_".$i;
				if (defined $symminterface{ $id }) {
					$id =~ s/-(\d)/n\1/g;
					print " JUMP_".$id."_to_subunit";
				}
			}
		}
	}
}
print "\n";



########################################
# output complete symm complex
########################################
########################################
## write output pdb
my $outpdb = $pdbfile;
my $outmon = $pdbfile;
my $outmdl = $pdbfile;

my $suffix = "_model_$primary_chain"."$helical_chain"."$perp_chain"."$ncs_chain";
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
my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

foreach my $sec_shift ( -$nperp_repeats..$nperp_repeats ) {

	foreach my $subunit (0 .. $nsubunits_to_gen) {
		foreach my $i (0..$sym_order_ncs-1) {
			my $id =  $sec_shift."_".$subunit."_".$i;
	
			foreach my $line (@filebuf) {
				my $linecopy = $line;
		
				my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
				my $X_0 = vsub($X,$COM_0);
				my $rX = vadd( mapply($Rs->{ $id }, $X_0) ,$Ts->{ $id } );
		
				substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
				substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
				substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
				substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);
		
				print OUTPDB $linecopy."\n";
		
				if (defined $symminterface{ $id }) {
					print OUTMDL $linecopy."\n";
				}
			}
			print OUTPDB "TER   \n";
			$chnidx++;
		}
		next if ($subunit == 0);
	
		foreach my $i (0..$sym_order_ncs-1) {
			my $id = $sec_shift."_".-$subunit."_".$i;
			foreach my $line (@filebuf) {
				my $linecopy = $line;
		
				my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
				my $X_0 = vsub($X,$COM_0);
				my $rX = vadd( mapply($Rs->{ $id }, $X_0) ,$Ts->{ $id } );
		
				substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
				substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
				substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
				substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);
		
				print OUTPDB $linecopy."\n";
		
				if (defined $symminterface{ $id }) {
					print OUTMDL $linecopy."\n";
				}
			}
			print OUTPDB "TER   \n";
			$chnidx++;
		}
	}

}
foreach my $debug_line (@fakepdblines) {
	print OUTPDB $debug_line;
}

######################################
foreach my $line (@filebuf) {
	my $linecopy = $line;
	my $id = "0_0_0";

	my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
	my $X_0 = vsub($X,$COM_0);
	my $rX = vadd( mapply($Rs->{ $id }, $X_0) ,$Ts->{ $id } );

	substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
	substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
	substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
	substr ($linecopy, 21, 1) = "A";

	print OUTMON $linecopy."\n";
}
#foreach my $debug_line (@{ $debug_lines }) {
#	print OUTPDB $debug_line;
# }
# 
close(OUTPDB);
close(OUTMON);
close (OUTMDL);




#############################################################################
#############################################################################

###
###   SUBROUTINES
###

#############################################################################
#############################################################################
