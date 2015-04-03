#!/usr/bin/perl -w

use strict;
use Storable qw(dclone);
use List::Util qw[min max];
use POSIX;
use Getopt::Long;

# dump stack on crash
$SIG{ __DIE__ } = sub { show_call_stack() };


if ($#ARGV < 1) {
	print STDERR "usage: $0 (I|O|T|I213|P4132|P23) [--eat N] [--cc] <pdbfile> <trimer-file> <dimer|trimer-file>\n";
	exit 1;
}

my ($EAT_INTO) = (0);
GetOptions(
	"eat=i" => \$EAT_INTO
);

my $symmtype = shift @ARGV;
my $pdbfile = shift @ARGV;
my $file1 = shift @ARGV;
my $file2 = shift @ARGV;

# from Will
my $HELIX_D_PER_RES = 1.53663757571;
my $HELIX_CA_RADIUS = 2.27024021944;
my $M_PI = 3.14159;

##
## PARAMETERS
##
my $NTHETA1STEPS = 6;         # bruteforce theta1 sampling (cryst symm only!)
my $NGRIDSEARCHSTEPS_th=120;  # grid seach for helix vector alignment, use this many ang. samples

my $MAX_ANGLE   = 9999;       # max angle offset in initial screen
my $MAX_TRANS   = 9999;       # max trans offset in initial screen
my $MAX_RMS     = 0.6;       # max rms of helical align
my $INIT_RMS    = 2.0;        # early RMS check (after 10 cycles opt)
my $NCLASH      = 0;          # clash check criteria (#Ca's w/i 3A)
my $NCYC        = 50;         # cycles of rms optimization
my $CLASHDIST   = 3.0;        # ca-ca clash distance  (<=0 means do not clash check)
my $TERMTRIM    = 8;          # trim this many terminal residues before clash checks
my $FORCE_Z     = 30;         # force Z to be at least this large (to favor large unit cell dimensions)
my $RSCALE      = 4;          # scale step size in refinement (4 seems to give better results)
my $VERBOSE     = 1;

# size limit! if two-comp, this is res limit per component
my $MAXRES = 99999;

# helix extension parameters
my $HELIX_ALIGN_WINDOW = 6;

## helix extension
## seperate parameters for dimer and trimer extensions match helical supercoil params from coiled coils
my ($HELIX2,$HELIX3);

## 11-res helix extension
## seperate parameters for dimer and trimer extensions match helical supercoil params from coiled coils
my $HELIX2cc = [
	[ -2.815 , -6.016 , 17.190 ],[ -3.736 , -5.792 , 16.093 ],[ -3.196 , -4.802 , 15.074 ],[ -3.351 , -4.974 , 13.866 ],
	[ -2.592 , -3.745 , 15.563 ],[ -1.968 , -2.759 , 14.702 ],[ -0.828 , -3.341 , 13.883 ],[ -0.656 , -3.027 , 12.707 ],
	[ -0.026 , -4.156 , 14.511 ],[ 1.064 , -4.826 , 13.828 ],[ 0.578 , -5.738 , 12.713 ],[ 1.171 , -5.812 , 11.639 ],
	[ -0.476 , -6.473 , 12.981 ],[ -1.085 , -7.329 , 11.981 ],[ -1.627 , -6.546 , 10.797 ],[ -1.506 , -6.957 , 9.645 ],
	[ -2.264 , -5.430 , 11.084 ],[ -2.762 , -4.546 , 10.048 ],[ -1.647 , -3.971 , 9.190 ],[ -1.770 , -3.849 , 7.973 ],
	[ -0.575 , -3.582 , 9.829 ],[ 0.586 , -3.074 , 9.125 ],[ 1.197 , -4.106 , 8.190 ],[ 1.617 , -3.796 , 7.077 ],
	[ 1.290 , -5.323 , 8.656 ],[ 1.794 , -6.413 , 7.843 ],[ 0.931 , -6.671 , 6.618 ],[ 1.429 , -6.938 , 5.526 ],
	[ -0.370 , -6.637 , 6.806 ],[ -1.307 , -6.785 , 5.708 ],[ -1.186 , -5.667 , 4.686 ],[ -1.264 , -5.888 , 3.479 ],
	[ -1.033 , -4.458 , 5.173 ],[ -0.833 , -3.310 , 4.309 ],[ 0.443 , -3.414 , 3.490 ],[ 0.481 , -3.061 , 2.313 ],
	[ 1.495 , -3.859 , 4.120 ],[ 2.759 , -4.062 , 3.438 ],[ 2.659 , -5.094 , 2.326 ],[ 3.236 , -4.937 , 1.251 ],
	[ 1.967 , -6.174 , 2.595 ],[ 1.732 , -7.202 , 1.599 ],[ 0.932 , -6.688 , 0.413 ],[ 1.202 , -7.025 , -0.738 ]
];
my $HELIX3cc = [
	[ 1.318 , 7.995 , -26.543 ],[ 0.500 , 8.338 , -25.396 ],[ 0.422 , 7.211 , -24.379 ],[ 0.463 , 7.432 , -23.171 ],
	[ 0.276 , 6.002 , -24.872 ],[ 0.261 , 4.829 , -24.019 ],[ 1.571 , 4.637 , -23.273 ],[ 1.594 , 4.271 , -22.099 ],
	[ 2.658 , 4.844 , -23.960 ],[ 3.971 , 4.749 , -23.350 ],[ 4.164 , 5.758 , -22.230 ],[ 4.749 , 5.460 , -21.190 ],
	[ 3.718 , 6.972 , -22.455 ],[ 3.772 , 8.008 , -21.442 ],[ 2.940 , 7.667 , -20.216 ],[ 3.338 , 7.915 , -19.080 ],
	[ 1.759 , 7.129 , -20.449 ],[ 0.897 , 6.679 , -19.373 ],[ 1.515 , 5.547 , -18.569 ],[ 1.411 , 5.497 , -17.345 ],
	[ 2.127 , 4.623 , -19.261 ],[ 2.812 , 3.520 , -18.614 ],[ 3.955 , 3.987 , -17.726 ],[ 4.173 , 3.467 , -16.634 ],
	[ 4.713 , 4.933 , -18.209 ],[ 5.797 , 5.510 , -17.438 ],[ 5.309 , 6.190 , -16.169 ],[ 5.925 , 6.093 , -15.110 ],
	[ 4.219 , 6.923 , -16.282 ],[ 3.598 , 7.558 , -15.135 ],[ 3.095 , 6.550 , -14.115 ],[ 3.219 , 6.742 , -12.907 ],
	[ 2.494 , 5.490 , -14.604 ],[ 2.030 , 4.415 , -13.747 ],[ 3.167 , 3.736 , -13.001 ],[ 3.050 , 3.392 , -11.826 ],
	[ 4.249 , 3.506 , -13.690 ],[ 5.425 , 2.915 , -13.082 ],[ 5.994 , 3.775 , -11.965 ],[ 6.421 , 3.277 , -10.925 ],
	[ 6.049 , 5.066 , -12.193 ],[ 6.500 , 6.004 , -11.184 ],[ 5.604 , 6.012 , -9.956 ],[ 6.070 , 6.090 , -8.821 ]
];
my $HELIXidl = [
	[ 0.370 , 3.598 , 2.351 ],[ 1.039 , 3.821 , 3.641 ],[ 2.491 , 3.342 , 3.566 ],[ 3.418 , 4.049 , 3.988 ],
	[ 2.634 , 2.146 , 3.027 ],[ 3.942 , 1.496 , 2.858 ],[ 4.862 , 2.382 , 2.015 ],[ 6.030 , 2.606 , 2.366 ],
	[ 4.295 , 2.858 , 0.922 ],[ 4.999 , 3.730 , -0.029 ],[ 5.506 , 4.984 , 0.686 ],[ 6.667 , 5.386 , 0.524 ],
	[ 4.606 , 5.561 , 1.460 ],[ 4.881 , 6.777 , 2.238 ],[ 6.063 , 6.539 , 3.181 ],[ 6.987 , 7.362 , 3.267 ],
	[ 5.989 , 5.409 , 3.859 ],[ 7.017 , 4.983 , 4.820 ],[ 8.378 , 4.895 , 4.125 ],[ 9.391 , 5.393 , 4.639 ],
	[ 8.350 , 4.257 , 2.970 ],[ 9.545 , 4.058 , 2.137 ],[ 10.169 , 5.410 , 1.784 ],[ 11.389 , 5.600 , 1.900 ],
	[ 9.299 , 6.308 , 1.361 ],[ 9.683 , 7.672 , 0.969 ],[ 10.384 , 8.373 , 2.135 ],[ 11.443 , 8.995 , 1.962 ],
	[ 9.761 , 8.246 , 3.291 ],[ 10.260 , 8.840 , 4.539 ],[ 11.671 , 8.328 , 4.836 ],[ 12.578 , 9.107 , 5.161 ],
	[ 11.804 , 7.020 , 4.711 ],[ 13.074 , 6.319 , 4.948 ],[ 14.161 , 6.878 , 4.028 ],[ 15.282 , 7.172 , 4.469 ],
	[ 13.785 , 7.007 , 2.770 ],[ 14.672 , 7.525 , 1.718 ],[ 15.161 , 8.926 , 2.090 ],[ 16.359 , 9.230 , 1.993 ],
	[ 14.206 , 9.736 , 2.506 ],[ 14.455 , 11.126 , 2.913 ],[ 15.476 , 11.162 , 4.052 ],[ 16.434 , 11.949 , 4.026 ]
];
my $NHELIX = 11;  # length of all helices

my ($TARGET_TH, $Rglobal, $pointsymm);
if ($symmtype eq "I213") {
	# "my" canonic to "crystal" canonic
	#    rotate [0,0,1] to [1,1,1]
	#    rotate [$sinth,0,$costh] to [1,0,0]
	$Rglobal = [
	  [ 8.16496580927726e-01 ,  5.55111512312578e-17 ,  5.77350269189626e-01],
	  [-4.08248290463863e-01 ,  7.07106781186547e-01 ,  5.77350269189626e-01],
	  [-4.08248290463863e-01 , -7.07106781186548e-01 ,  5.77350269189626e-01]
	];
	$TARGET_TH = 54.736; # acos(1/sqrt(3))
	$pointsymm = 0;
} elsif ($symmtype eq "P4132") {
	# "my" canonic to "crystal" canonic
	#    rotate [0,0,1] to [1,1,1]
	#    rotate [$sinth,0,$costh] to [1,-1,1]
	$Rglobal = [
		[4.08248290463863e-01 , -7.07106781186548e-01 ,  5.77350269189626e-01],
		[4.08248290463863e-01 ,  7.07106781186548e-01 ,  5.77350269189626e-01],
		[-8.16496580927726e-01 ,  2.78911008871896e-17 ,  5.77350269189626e-01]
	];
	$TARGET_TH = 35.264; # acos(sqrt(2/3))
	$pointsymm = 0;
}  elsif ($symmtype eq "P23") {
	# "my" canonic to "crystal" canonic
	#    rotate [0,0,1] to [1,1,1]
	#    rotate [$sinth,0,$costh] to [1,1,-1]
	$Rglobal = [
	  [ 4.08248290463863e-01 ,  7.07106781186548e-01 ,  5.77350269189626e-01],
	  [-8.16496580927726e-01 ,  1.44253099049196e-16 ,  5.77350269189626e-01],
	  [ 4.08248290463863e-01 , -7.07106781186548e-01 ,  5.77350269189626e-01]
	];
	$TARGET_TH = 70.529; # acos(1/3)
	$pointsymm = 0;
} elsif ($symmtype eq 'T') {
	$TARGET_TH = 54.736; # acos(1/sqrt(3))
	$Rglobal = [ [1,0,0],[0,1,0],[0,0,1] ];
	$pointsymm = 1;
} elsif ($symmtype eq 'O') {
	$TARGET_TH = 35.264; # acos(sqrt(2/3))
	$Rglobal = [ [1,0,0],[0,1,0],[0,0,1] ];
	$pointsymm = 1;
} elsif( $symmtype eq 'I') {
	$TARGET_TH = 20.905; # 90 - acos(-sqrt(5)/3) / 2
	$Rglobal = [ [1,0,0],[0,1,0],[0,0,1] ];
	$pointsymm = 1;
} else {
	die "Unrecognized symm type = $symmtype\n";
}


my $SINTH = sin( $TARGET_TH*3.141593/180.0 );
my $COSTH = cos( $TARGET_TH*3.141593/180.0 );


open (FILE1, $file1) || die "$0: unable to open file $file1 for reading";
my @file1 = <FILE1>;
close (FILE1);

open (FILE2, $file2) || die "$0: unable to open file $file2 for reading";
my @file2 = <FILE2>;
close (FILE2);


# read PDB; get angle and transform between termini
my ($pdbvec,$pdbmid,$inpdb_N,$inpdb_C, $all_cas3) = getTermHelixVectors($pdbfile, $HELIX_ALIGN_WINDOW);

# /work/dimaio/databases/symm_complexes/C3/1N1QA_C3.pdb  124  148 73.2680520552772 7.72399542581519 C
my $nhits = 0;
foreach my $line1 (@file1) {
	chomp $line1; #$line1 =~ s/-/ /;
	my @fields1 = split ' ', $line1;

	# truncate model
	my ($Ntrunc_at3,$Ctrunc_at3) = (-99999,99999);
	if    ($fields1[6] eq "N") { $Ntrunc_at3 = $fields1[1]; }
	elsif ($fields1[6] eq "C") { $Ctrunc_at3 = $fields1[2]; }

	# hacky
	if ($fields1[0] =~ /woolfson/) { # CC
		if ($symmtype eq "P23") {
			$HELIX3 = dclone($HELIX3cc);
		} else {
			$HELIX3 = dclone($HELIX3cc);
		}
	} else {
		$HELIX3 = dclone($HELIXidl);
	}	

	my ($helix_vector1,$helix_cas1,$helix_aas1,$all_cas1) = getHelixVector( $fields1[0],$fields1[3], [$fields1[1],$fields1[2]],$fields1[6] );

	# size filter
	next if ( scalar(@{$all_cas1}) > $MAXRES );

	my ($symmaxis1,$symmcenter1) = getSymmAxis( $fields1[0], "B" );

	# NOTE: c3_is_nterm means the n-term _extension_ comes off C3, actually making it the C-term of the fusion
	my $c3_is_nterm = ($fields1[6] eq "N");
	my $inv_fields1_6;
	if ($c3_is_nterm) {
		$inv_fields1_6 = "C";
	} else {
		$inv_fields1_6 = "N";
	}

	# transformation to canonic frame
	my $Rflip = [ [1,0,0],[0,-1,0],[0,0,-1] ];
	my $axisZ1 = $symmaxis1;
	my $axisX1scale =  1; if ($fields1[6] eq "N") { $axisX1scale = -1; }
	my $axisX1 = vscale( $axisX1scale, norm ( vsub( $helix_vector1 , vscale( vdot( $helix_vector1, $axisZ1), $axisZ1) ) ) );
	my $axisY1 = vcross( $axisZ1, $axisX1 );
	my $R1 = minv([[$axisX1->[0],$axisY1->[0],$axisZ1->[0]],[$axisX1->[1],$axisY1->[1],$axisZ1->[1]],[$axisX1->[2],$axisY1->[2],$axisZ1->[2]]]);
	my $angle1 = 180/3.141593*vangle( $helix_vector1, $symmaxis1 );

	my $compute_C3_flip = 0;
	if ($pointsymm == 0) {$compute_C3_flip = 1;}

	foreach my $C3flip (0..$compute_C3_flip) {
		if ($C3flip==1) { $R1 = mmult( $Rflip, $R1 ); }

		# compute transformed coordinates for 1st-pass clash check
		my $all_cas1_xform = [];
		foreach my $i (0..scalar(@{ $all_cas1 })-1) {
			push @{ $all_cas1_xform } , mapply( $R1, vsub( $all_cas1->[$i], $symmcenter1 ) );
		}

		# extend helix
		my ($c3_helix_end, $c3_helix_stub) = ([],[]);
		foreach my $i (0..$HELIX_ALIGN_WINDOW-1) {
			foreach my $j (0..3) {
				my $idx = 4*$i+$j;
				if ($c3_is_nterm) {
					push @{ $c3_helix_end }, mapply( $R1 , vsub( $helix_aas1->[$idx], $symmcenter1) );
					push @{ $c3_helix_stub }, dclone( $HELIX3->[-4*$HELIX_ALIGN_WINDOW+$idx] );
				} else {
					push @{ $c3_helix_end }, mapply( $R1 , vsub( $helix_aas1->[-4*$HELIX_ALIGN_WINDOW+$idx], $symmcenter1) );
					push @{ $c3_helix_stub }, dclone( $HELIX3->[$idx] );
				}
			}
		}
		my ($R_c3,$rmsd_c3, $COM_i_c3, $COM_ij_c3) = rms_align( $c3_helix_end, $c3_helix_stub  ); # potential could RMS filt here
		my ($c3_extend) = ([]);
		foreach my $i (0..$NHELIX-1) {
			foreach my $j (0..3) {
				my $idx = 4*$i+$j;
				push @{ $c3_extend }, vadd( mapply( $R_c3, vsub( $HELIX3->[$idx], $COM_i_c3 ) ), vadd( $COM_i_c3, $COM_ij_c3) );
			}
		}

		foreach my $C3i (0..$NHELIX-$HELIX_ALIGN_WINDOW) {
			my ($c3_extend_window) = ([]);
			foreach my $C3j ($C3i..$C3i+$HELIX_ALIGN_WINDOW-1) {
				foreach my $C3j_atm (0..3) {
					push @{$c3_extend_window}, dclone( $c3_extend->[4*$C3j+$C3j_atm]);
				}
			}
			my ($c3_helixvec, $c3_helixend) = getHelixVectorParsed( $c3_extend_window, $fields1[6] );
			my ($c3_helixvecI, $c3_helixendI) = getHelixVectorParsed( $c3_extend_window, $inv_fields1_6 );

			# align linker onto extension window
			my ($lnkr_c3_bundle, $lnkr_c2_bundle);
			if ($c3_is_nterm) {
				$lnkr_c3_bundle = dclone( $inpdb_C );
				$lnkr_c2_bundle = dclone( $inpdb_N );
			} else {
				$lnkr_c3_bundle = dclone( $inpdb_N );
				$lnkr_c2_bundle = dclone( $inpdb_C );
			}
			my $c3_extend_window_copy = dclone( $c3_extend_window );
			my ($R_lnkr,$rmsd_lnkr, $COM_i_lnkr, $COM_ij_lnkr) = rms_align( $c3_extend_window_copy, $lnkr_c3_bundle  );


			# transform the linker, and clash check linker+trimer
			# this will throw out buried helices right away
			my $all_cas3_xform = [];
			foreach my $i (0..scalar(@{ $all_cas3 })-1) {
				push @{ $all_cas3_xform } , vadd( mapply( $R_lnkr, vsub( $all_cas3->[$i], $COM_i_lnkr ) ), vadd( $COM_i_lnkr, $COM_ij_lnkr ) );
			}
			my $nclashAB = 0;
			if ($CLASHDIST > 0) {
				outerloop: foreach my $i ($TERMTRIM..scalar(@{ $all_cas1_xform })-$TERMTRIM-1) {  # ignore terminal clashes
					foreach my $j ($TERMTRIM..scalar(@{ $all_cas3_xform })-$TERMTRIM-1) {  # ignore terminal clashes
						my $d = len2( vsub( $all_cas1_xform->[$i], $all_cas3_xform->[$j] ) );
						if ($d<$CLASHDIST*$CLASHDIST) {
							$nclashAB++;
							last outerloop if ($nclashAB > $NCLASH);
						}
					}
				}
			}
			next if ($nclashAB > $NCLASH);


			# apply this transformation to the linkers
			foreach my $i (0..scalar(@{ $lnkr_c3_bundle })-1) {
				$lnkr_c3_bundle->[$i] = vadd( mapply( $R_lnkr, vsub( $lnkr_c3_bundle->[$i], $COM_i_lnkr ) ), vadd( $COM_i_lnkr, $COM_ij_lnkr ) );
				$lnkr_c2_bundle->[$i] = vadd( mapply( $R_lnkr, vsub( $lnkr_c2_bundle->[$i], $COM_i_lnkr ) ), vadd( $COM_i_lnkr, $COM_ij_lnkr ) );
			}

			# get statistics on this vector
			my ($c3_linker_helixvec, $c3_linker_helixend) = getHelixVectorParsed( $lnkr_c2_bundle, $fields1[6] );
			my ($c3_linker_helixvecI, $c3_linker_helixendI) = getHelixVectorParsed( $lnkr_c2_bundle, $inv_fields1_6 );

			my $thetaC3 = 180/3.141593*vangle( $c3_linker_helixvec, [0,0,1] );
			my $dC3 = vdot( $c3_linker_helixend , norm(vcross( $c3_linker_helixvec, [0,0,1] )) );

			# now scan through file2 to find c2's compatible with this model
			foreach my $line2 (@file2) {
				chomp $line2;
				my @fields2 = split ' ', $line2;

				next if ( $fields1[6] eq $fields2[6] ); # must be N<->C

				if ($fields2[0] =~ /woolfson/) { # CC
					if ($symmtype eq "P23") {
						$HELIX2 = dclone($HELIX3cc);
					} else {
						$HELIX2 = dclone($HELIX2cc);
					}
				} else {
					$HELIX2 = dclone($HELIXidl);
				}	

				# truncate model
				my ($Ntrunc_at2,$Ctrunc_at2) = (-99999,99999);
				if    ($fields2[6] eq "N") { $Ntrunc_at2 = $fields2[1]; }
				elsif ($fields2[6] eq "C") { $Ctrunc_at2 = $fields2[2]; }

				### TODO: replace with Will's fastcheck
				my $th_offset;
				if ($c3_is_nterm) { $th_offset = (-$fields2[4] + $thetaC3 + $TARGET_TH); }
				else              { $th_offset = ( $fields2[4] - $thetaC3 + $TARGET_TH); }
				my $d_offset = $fields2[5] - $dC3;
				next if (abs($th_offset) > $MAX_ANGLE || abs($d_offset) > $MAX_TRANS);

				# align model 2 to proper symm axis
				my ($symmaxis2,$symmcenter2) = getSymmAxis( $fields2[0], "B" );
				my ($helix_vector2,$helix_cas2,$helix_aas2,$all_cas2) = getHelixVector( $fields2[0],$fields2[3], [$fields2[1],$fields2[2]],$fields2[6] );
				next if ( scalar(@{$all_cas2}) > $MAXRES );

				# rotate symm axes to tgt orientation
				my $axisZ2 = $symmaxis2;
				my $axisX2scale = -1; if ($fields2[6] eq "N") { $axisX2scale =  1; }
				my $axisX2 = vscale( $axisX2scale, norm ( vsub( $helix_vector2 , vscale( vdot( $helix_vector2, $axisZ2), $axisZ2) ) ) );
				my $axisY2 = vcross( $axisZ2, $axisX2 );
				my $R2_onZ =  minv([ [$axisX2->[0],$axisY2->[0],$axisZ2->[0]] ,
									 [$axisX2->[1],$axisY2->[1],$axisZ2->[1]] ,
									 [$axisX2->[2],$axisY2->[2],$axisZ2->[2]] ] );
				my $angle2 = 180/3.141593*vangle( $helix_vector2, $symmaxis2 );

				foreach my $C2flip (0..1) {
					if ($C2flip==1) {  $R2_onZ = mmult( $Rflip, $R2_onZ ); }
					my $tgt_angle = 3.141593 * $TARGET_TH/(180.0);
					my $Rtheta = quat2R( 0 , sin($tgt_angle/2), 0 , cos($tgt_angle/2));
					my $Rinvtheta = minv( $Rtheta );
					my $R2 = mmult( $Rtheta, $R2_onZ );

					# extend helix ...
					my ($c2_helix_end, $c2_helix_stub) = ([],[]);
					foreach my $i (0..$HELIX_ALIGN_WINDOW-1) {
						foreach my $j (0..3) {
							my $idx = 4*$i+$j;
							if ($c3_is_nterm) {
								push @{ $c2_helix_end }, mapply( $R2 , vsub( $helix_aas2->[-4*$HELIX_ALIGN_WINDOW+$idx], $symmcenter2) );
								push @{ $c2_helix_stub }, dclone( $HELIX2->[$idx] );
							} else {
								push @{ $c2_helix_end }, mapply( $R2 , vsub( $helix_aas2->[$idx], $symmcenter2) );
								push @{ $c2_helix_stub }, dclone( $HELIX2->[-4*$HELIX_ALIGN_WINDOW+$idx] );
							}
						}
					}
					my ($R_c2,$rmsd_c2, $COM_i_c2, $COM_ij_c2) = rms_align( $c2_helix_end, $c2_helix_stub  );

					## We _might_ want to filter RMS here to avoid helix distortion

					my ($c2_extend) = ([]);
					foreach my $i (0..$NHELIX-1) {
						foreach my $j (0..3) {
							my $idx = 4*$i+$j;
							push @{ $c2_extend }, vadd( mapply( $R_c2, vsub( $HELIX2->[$idx], $COM_i_c2 ) ), vadd( $COM_i_c2, $COM_ij_c2) );
						}
					}

					# find the best c2 offset
					my ($A,$B,$Z)=(0,0,0);
					my ($R1f,$R2f,$R3, $bundle_center, $bundle_trans);

					foreach my $C2i (0..$NHELIX-$HELIX_ALIGN_WINDOW) {
						my ($c2_extend_window) = ([]);
						foreach my $C2j ($C2i..$C2i+$HELIX_ALIGN_WINDOW-1) {
							foreach my $C2j_atm (0..3) {
								push @{$c2_extend_window}, dclone( $c2_extend->[4*$C2j+$C2j_atm]);
							}
						}
						my ($c2_helixvec, $c2_helixend) = getHelixVectorParsed( $c2_extend_window, $fields2[6] );
						#invert
						my ($c2_helixvecI, $c2_helixendI) = getHelixVectorParsed( $c2_extend_window, $fields1[6] );

						# initial guess of _5_ params
						my ($Rth1 , $Rth2);
						my ($theta1,$theta2,$B_i,$A_i,$Z_i);

						# **** set 5 parameters initially using a grid search with a 2-point estimate
						foreach my $i_th1 (0..$NTHETA1STEPS-1) {
							my $bestR = 1e30;
							my @bestParams;

							# precompute rotated points
							my (@c3_end_0s, @c3_end_1s, @c2_end_0s, @c2_end_1s);
						#foreach my $i_th1 (0..$NGRIDSEARCHSTEPS_th-1) {
							$theta1 =  $i_th1 * 2*$M_PI/$NTHETA1STEPS;
							$Rth1 = quat2R( 0, 0, sin($theta1/2), cos($theta1/2) );
							my $c3_linker_helixend_rot0 = mapply( $Rth1, vadd( $c3_linker_helixend, vscale( -3*$HELIX_D_PER_RES, $c3_linker_helixvec ) ) );
							my $c3_linker_helixend_rot1 = mapply( $Rth1, vadd( $c3_linker_helixend, vscale(  3*$HELIX_D_PER_RES, $c3_linker_helixvec ) ) );
							#push @c3_end_0s, $c3_linker_helixend_rot0;
							#push @c3_end_1s, $c3_linker_helixend_rot1;
						#}
							foreach my $i_th2 (0..$NGRIDSEARCHSTEPS_th-1) {
								$theta2 =  $i_th2 * 2*$M_PI/$NGRIDSEARCHSTEPS_th;
								$Rth2 = mmult( mmult( $Rtheta, quat2R( 0, 0, sin($theta2/2), cos($theta2/2) ) ), $Rinvtheta);
								my $c2_helixendI_rot0 = mapply( $Rth2, vadd( $c2_helixendI, vscale( -3*$HELIX_D_PER_RES, $c2_helixvecI ) ) );
								my $c2_helixendI_rot1 = mapply( $Rth2, vadd( $c2_helixendI, vscale(  3*$HELIX_D_PER_RES, $c2_helixvecI ) ) );
								push @c2_end_0s, $c2_helixendI_rot0;
								push @c2_end_1s, $c2_helixendI_rot1;
							}

						#	foreach my $i_th1 (0..$NGRIDSEARCHSTEPS_th-1) {
						#		$theta1 =  $i_th1 * 2*$M_PI/$NGRIDSEARCHSTEPS_th;
								foreach my $i_th2 (0..$NGRIDSEARCHSTEPS_th-1) {
									$theta2 =  $i_th2 * 2*$M_PI/$NGRIDSEARCHSTEPS_th;
									# solve for A,B,Z
									#my $c3_0 = $c3_end_0s[$i_th1];
									#my $c3_1 = $c3_end_1s[$i_th1];
									my $c3_0 = $c3_linker_helixend_rot0;
									my $c3_1 = $c3_linker_helixend_rot1;
									my $c2_0 = $c2_end_0s[$i_th2];
									my $c2_1 = $c2_end_1s[$i_th2];
									my $c3M = vscale(0.5,vadd( $c3_0, $c3_1 ));
									my $c2M = vscale(0.5,vadd( $c2_0, $c2_1 ));

									$B_i = ($c3M->[0] - $c2M->[0]) / $SINTH;
									$A_i = $c2M->[2] - $c3M->[2] + $B_i * $COSTH;
									if ($pointsymm == 0) {
										$Z_i = $c3M->[1] - $c2M->[1];
									} else {
										$Z_i = 0;
									}
									#$Z_i = $Z_base;

									my $res0 = vsub( vadd ($c3_0, [0,0,$A_i] ) , vadd( $c2_0 , [$B_i*$SINTH,$Z_i,$B_i*$COSTH])  );
									my $res1 = vsub( vadd ($c3_1, [0,0,$A_i] ) , vadd( $c2_1 , [$B_i*$SINTH,$Z_i,$B_i*$COSTH])  );
									my $residual = len2($res0) + len2($res1);
									if ($residual < $bestR) {
										$bestR = $residual;
										@bestParams = ($theta1,$theta2, $A_i,$B_i,$Z_i);
									}
								}
						#	}
							($theta1,$theta2, $A_i,$B_i,$Z_i) = @bestParams;
							#print "($theta1,$theta2, $A_i,$B_i,$Z_i) --> $bestR\n";

							# refinement
							my $bestRMS = 99.9;
							foreach my $cyc (1..$NCYC) {
								$Rth1 = quat2R( 0, 0, sin($theta1/2), cos($theta1/2) );
								$Rth2 = mmult( mmult( $Rtheta, quat2R( 0, 0, sin($theta2/2), cos($theta2/2) ) ), $Rinvtheta);

								# [1] align c3 CA's to c2 CAs
								my ($symm_helices, $linker_helices) = ([],[]);
								foreach my $jj (0..4*$HELIX_ALIGN_WINDOW-1) {
									my $c3shift;
									$c3shift = vadd( $c3_extend_window->[$jj], [0,0,$A_i] );
									$c3shift = mapply( $Rth1, $c3shift );

									push @{ $symm_helices }, $c3shift;
									if ($c3_is_nterm) {
										push @{ $linker_helices }, dclone( $inpdb_C->[$jj] );
									} else {
										push @{ $linker_helices }, dclone( $inpdb_N->[$jj] );
									}
								}
								foreach my $jj (0..4*$HELIX_ALIGN_WINDOW-1) {
									my $c2shift;
									$c2shift = vadd( $c2_extend_window->[$jj], [$B_i*$SINTH,0,$B_i*$COSTH] );
									$c2shift = mapply( $Rth2, $c2shift );
									$c2shift = vadd( $c2shift, [0,$Z_i,0] );

									push @{ $symm_helices }, $c2shift;
									if ($c3_is_nterm) {
										push @{ $linker_helices }, dclone( $inpdb_N->[$jj] );
									} else {
										push @{ $linker_helices }, dclone( $inpdb_C->[$jj] );
									}
								}

								# save best RMS
								my ($R_f,$rms_f, $COM_f, $COM_ij_f) = rms_align( $symm_helices , $linker_helices );

								if (($cyc == $NCYC) && $VERBOSE != 0) {
									printf "cycle %d: rms=%.2f c3=%d/%d c2=%d/%d: (A,B,th1,th2,Z) = %7.3f %7.3f %7.3f %7.3f %7.3f \n", 
										$cyc, $rms_f, $C3i, $C3flip, $C2i, $C2flip, $A_i, $B_i, $theta1, $theta2, $Z_i;
								}
								#print "\n" if ($cyc == $NCYC);
								#pseudoPDB( $symm_helices, [0,0,0] , "cage_".$nhits."D.pdb" );
								#pseudoPDB( $linker_helices, [0,0,0], "cage_".$nhits."E.pdb" );

								if ($rms_f < $bestRMS) {
									$bestRMS = $rms_f;
									$A = $A_i; $B = $B_i; $Z=$Z_i;
									$R1f = mmult( $Rth1, $R1 ); $R2f = mmult( $Rth2, $R2 ); $R3 = dclone( $R_f );
									$bundle_center = dclone( $COM_f ); $bundle_trans = vadd( $COM_f, $COM_ij_f );
								}
								last if (($cyc==10 && $rms_f > $INIT_RMS) || $cyc==$NCYC);

								# [2] translate/rotatate to minimize rms
								my ($r_ai,$r_theta1,$theta_wtsum1,$r_bi,$r_zi,$r_theta2,$theta_wtsum2) = (0,0,0,0,0,0,0);
								foreach my $jj (0..4*$HELIX_ALIGN_WINDOW-1) {
									my ($targetA,$targetB,$templateA,$templateB, $del_th, $rdel_th);
									$targetA = vadd( $c3_extend_window->[$jj], [0,0,$A_i] );
									$targetA = mapply( $Rth1, $targetA );

									$targetB = vadd( $c2_extend_window->[$jj], [$B_i*$SINTH,0,$B_i*$COSTH] );
									$targetB = mapply( $Rth2, $targetB );
									$targetB = vadd( $targetB, [0,$Z_i,0] );

									if ($c3_is_nterm) {
										$templateA =  vadd( mapply( $R_f, vsub($inpdb_C->[$jj], $COM_f ) ), vadd( $COM_f, $COM_ij_f ) );
										$templateB =  vadd( mapply( $R_f, vsub($inpdb_N->[$jj], $COM_f ) ), vadd( $COM_f, $COM_ij_f ) );
									} else {
										$templateA =  vadd( mapply( $R_f, vsub($inpdb_N->[$jj], $COM_f ) ), vadd( $COM_f, $COM_ij_f ) );
										$templateB =  vadd( mapply( $R_f, vsub($inpdb_C->[$jj], $COM_f ) ), vadd( $COM_f, $COM_ij_f ) );
									}
									## A
									$r_ai += ($targetA->[2] - $templateA->[2]);

									## THETA_A
									$del_th = atan2( $targetA->[1], $targetA->[0]) - atan2( $templateA->[1], $templateA->[0]);
									$del_th += 2*$M_PI if ($del_th < -$M_PI); $del_th -= 2*$M_PI if ($del_th > $M_PI);
									my $thwt = sqrt( $targetA->[0]*$targetA->[0] + $targetA->[1]*$targetA->[1] );
									$theta_wtsum1 += $thwt;
									$r_theta1 += $thwt*($del_th);

									## B (rotate 2-fold to Z)
									my $r_target = mapply($Rinvtheta, $targetB);
									my $r_template = mapply($Rinvtheta, $templateB);
									$r_bi += ($r_target->[2] - $r_template->[2]);

									## THETA_B
									$rdel_th = atan2( $r_target->[1], $r_target->[0]) - atan2( $r_template->[1], $r_template->[0]);
									$rdel_th += 2*$M_PI if ($rdel_th < -$M_PI); $rdel_th -= 2*$M_PI if ($rdel_th > $M_PI);
									my $rthwt = sqrt( $r_target->[0]*$r_target->[0] + $r_target->[1]*$r_target->[1] );
									$theta_wtsum2 += $rthwt;
									$r_theta2 += $rthwt*($rdel_th);

									## Z
									$r_zi += ($r_target->[1] - $r_template->[1]);
								}
								$A_i -= $RSCALE * $r_ai/(4*$HELIX_ALIGN_WINDOW) if ($cyc % 2 == 0);
								$theta1 -= $RSCALE * $r_theta1/$theta_wtsum1    if ($cyc % 2 == 0);
								$B_i -= $RSCALE * $r_bi/(4*$HELIX_ALIGN_WINDOW) if ($cyc % 2 == 1);
								$theta2 -= $RSCALE * $r_theta2/$theta_wtsum2    if ($cyc % 2 == 1);
								$Z_i -= $r_zi/$HELIX_ALIGN_WINDOW               if ($cyc % 2 == 0 && $pointsymm == 0);
							} # for cyc in [1,NCYC]

							if ($bestRMS > $MAX_RMS) {
								#print "Fail! rms=$bestRMS > $MAX_RMS\n";
							} else {
								#print "Pass! rms=$bestRMS < $MAX_RMS\n";
								my $all_cas1c = transform_ca( $all_cas1, $R1f, $symmcenter1, [0,0,$A], $Z);
								my $all_cas2c = transform_ca( $all_cas2, $R2f, $symmcenter2, [$B_i*$SINTH,$Z,$B_i*$COSTH], $Z);
								my $all_cas3c = transform_ca( $all_cas3, $R3, $bundle_center, $bundle_trans, $Z);

								my $complex_ca = build_complex_ca(  $all_cas1c , $all_cas2c, $all_cas3c );
								my $nclash = 0;

								if ($CLASHDIST > 0) {
									$nclash = clash_check( $complex_ca, $Z );
									#my $clashchecked = clash_check_verbose( $complex_ca, $Z );
									#pseudoPDB( $clashchecked, [0,0,0] , "mdls".$nhits.".pdb" );
									#print "$nclash \n";
									#exit;
								}

								##
								if ($nclash > $NCLASH) {
									print "Fail! clash>=$nclash with rms=$bestRMS\n";
								} else {
									$nhits++;
									my ($Acryst,$offset) = get_A_offset_fromZ($Z);


									# transform model 1 by R1, model 2 by R2
									my $pdb1 = $fields1[0];# $pdb1 =~ s/(_C[2-3])/$1_0001/;
									my $pdb2 = $fields2[0];# $pdb2 =~ s/(_C[2-3])/$1_0001/;

									my $pdb1stem = $pdb1; $pdb1stem =~ s/.*\/([^\/]+)/$1/; $pdb1stem =~ s/([^\/]+)_0001\.pdb/$1/;
									my $pdb2stem = $pdb2; $pdb2stem =~ s/.*\/([^\/]+)/$1/; $pdb2stem =~ s/([^\/]+)_0001\.pdb/$1/;
									my $pdb3stem = $pdbfile; $pdb3stem =~ s/.*\/([^\/]+)/$1/; $pdb3stem =~ s/([^\/]+)\.pdb/$1/;
									my $outfile = $symmtype.'__'.$pdb1stem.'__'.$pdb3stem.'__'.$pdb2stem.'__e'.$EAT_INTO.'_'.$nhits;

									my ($overlapN, $overlapC);
									if ($c3_is_nterm) {
										$overlapN = 7-$C2i; #NOTE: this should probably be updated if ideal helix length changes!
										$overlapC = $C3i+1;
										transformPDBs(
											$pdb2, $R2f, $symmcenter2, [$B*$SINTH,$Z,$B*$COSTH],
											$Ctrunc_at2,
											$pdbfile, $R3, $bundle_center, $bundle_trans,
											$pdb1, $R1f, $symmcenter1, [0,0,$A],
											$Ntrunc_at3,
											$overlapN, $overlapC, $Z, $bestRMS, $nclash, "comp_$outfile.pdb");
									} else {
										$overlapN = 7-$C3i; #NOTE: this should probably be updated if ideal helix length changes!
										$overlapC = $C2i+1;
										transformPDBs(
											$pdb1, $R1f, $symmcenter1, [0,0,$A],
											$Ctrunc_at3,
											$pdbfile, $R3, $bundle_center, $bundle_trans,
											$pdb2, $R2f, $symmcenter2, [$B*$SINTH,$Z,$B*$COSTH],
											$Ntrunc_at2,
											$overlapN, $overlapC, $Z, $bestRMS, $nclash, "comp_$outfile.pdb");
									}

									print "Model $nhits [".$Acryst."] $c3_is_nterm $C3i $C2i".
										" ".$fields1[0].":".$fields1[3].":".$fields1[6].
										" ".$fields2[0].":".$fields2[3].":".$fields2[6]." clash=$nclash rms=$bestRMS\n";
								}
							}
						}  #ith1
					} # for C2i in helix ext
				} # for C2flip=(0,1)
			} # for line in file 2
		} # for C3i in helix ext
	} # for C3flip=(0,[0,1])
} # for line in file 1


##
## debugging purposes: output pseudoatomic model
##
sub pseudoPDB {
	my ($coords,$T12,$outpdb) = @_;
	open (OUTPDB, ">$outpdb");
	my $ctr = 1;
	foreach my $line (@{$coords}) {
		printf OUTPDB "ATOM   %4d  CA  ALA G%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", $ctr%10000, $ctr%10000, $line->[0]+$T12->[0],  $line->[1]+$T12->[1],  $line->[2]+$T12->[2];
		$ctr++;
	}
}


##
##
##
sub transformPDB {
	my ($inpdb1,$R1,$T11,$T12,$outpdb) = @_;
	open (OUTPDB, ">$outpdb");
	open (INPDB1, $inpdb1);
	while (<INPDB1>) {
		chomp;
		if (/^ATOM/ || /^HETATM/) {
			my $chnid = substr ($_, 21, 1);
			my $atom  = substr ($_, 12, 4);
			my $line = $_;

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
			my $X_0 = vsub($X,$T11);
			my $rX = vadd(  mapply($R1, $X_0) , $T12 );

			substr ($line, 30, 8) = sprintf ("%8.3f", $rX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $rX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $rX->[2]);

			print OUTPDB $line."\n";
		}
	}
}

sub transformPDBs {
	my ($inpdb1,$R1,$T11,$T12,$ctrunc,$inpdb2,$R2,$T21,$T22,$inpdb3,$R3,$T31,$T32,$ntrunc,$overlapN,$overlapC,$Z,$bestRMS,$clash,$outpdb) = @_;
	open (OUTPDB, ">$outpdb");

	# trim so no overlap
	my $cut1C = floor($overlapN/2);
	my $cut2N = $overlapN-$cut1C;

	my $cut3N = floor($overlapC/2);
	my $cut2C = $overlapC-$cut3N;

	my ($A,$offset) = get_A_offset_fromZ($Z);

	printf OUTPDB "REMARK  RMS=$bestRMS clash=$clash\n", $A, $A, $A;

	if ($symmtype eq "I213") {
		printf OUTPDB "CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 I 21 3\n", $A, $A, $A;
	} elsif ($symmtype eq "P23") {
		printf OUTPDB "CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 2 3\n", $A, $A, $A;
	} elsif ($symmtype eq "P4132") {
		if ($Z > 0) {
			printf OUTPDB "CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 41 3 2\n", $A, $A, $A;
		} else {
			printf OUTPDB "CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 43 3 2\n", $A, $A, $A;
		}
	}

	open (INPDB1, $inpdb1);
	my @lines1 = <INPDB1>;
	chomp @lines1;
	close (INPDB1);
	my 	($firstresid,$lastresid) = (99999,0);
	#foreach my $line (@lines1) {
	#	if ($line =~/^ATOM/ || $line =~/^HETATM/) {
	#		my $resid = int( substr ($line, 22, 4) ) ;
	#		$lastresid = max($lastresid,$resid);
	#	}
	#}
	foreach my $line (@lines1) {
		if ($line =~/^ATOM/ || $line =~/^HETATM/) {
			my $chnid = substr ($line, 21, 1);
			my $atom  = substr ($line, 12, 4);

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
			my $X_0 = vsub($X,$T11);
			my $rX = vadd( $offset , mapply( $Rglobal, vadd(  mapply($R1, $X_0) , $T12 ) ) );
			my $resid = int( substr ($line, 22, 4) ) ;

			#next if ($resid>$lastresid-$cut1C);
			next if ($resid>$ctrunc-$cut1C+1);

			#substr ($line, 21, 1) = "A";
			substr ($line, 30, 8) = sprintf ("%8.3f", $rX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $rX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $rX->[2]);
			substr ($line, 22, 4) = sprintf "%4d", $resid;

			print OUTPDB $line."\n";
		}
	}

	open (INPDB2, $inpdb2);
	my @lines2 = <INPDB2>;
	chomp @lines2;
	close (INPDB2);
	($firstresid,$lastresid) = (99999,0);
	foreach my $line (@lines2) {
		if ($line =~/^ATOM/ || $line =~/^HETATM/) {
			my $resid = int( substr ($line, 22, 4) ) ;
			$lastresid = max($lastresid,$resid);
			$firstresid = min($firstresid,$resid);
		}
	}
	foreach my $line (@lines2) {
		if ($line =~/^ATOM/ || $line =~/^HETATM/) {
			my $chnid = substr ($line, 21, 1);
			my $atom  = substr ($line, 12, 4);

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
			my $X_0 = vsub($X,$T21);
			my $rX = vadd( $offset , mapply( $Rglobal, vadd(  mapply($R2, $X_0) , $T22 ) ) );
			my $resid = int( substr ($line, 22, 4) ) ;

			next if ($resid>$lastresid-$cut2C-$EAT_INTO);
			next if ($resid<$firstresid+$cut2N+$EAT_INTO);

			substr ($line, 30, 8) = sprintf ("%8.3f", $rX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $rX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $rX->[2]);
			substr ($line, 22, 4) = sprintf "%4d", $resid+1000;

			print OUTPDB $line."\n";
		}
	}

	open (INPDB3, $inpdb3);
	my @lines3 = <INPDB3>;
	chomp @lines3;
	close (INPDB2);
	($firstresid,$lastresid) = (99999,0);
	#foreach my $line (@lines3) {
	#	if ($line =~/^ATOM/ || $line =~/^HETATM/) {
	#		my $resid = int( substr ($line, 22, 4) ) ;
	#		$firstresid = min($firstresid,$resid);
	#	}
	#}
	foreach my $line (@lines3) {
		if ($line =~/^ATOM/ || $line =~/^HETATM/) {
			my $chnid = substr ($line, 21, 1);
			my $atom  = substr ($line, 12, 4);

			my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
			my $X_0 = vsub($X,$T31);
			my $rX = vadd( $offset ,  mapply( $Rglobal, vadd(  mapply($R3, $X_0) , $T32 ) ) );
			my $resid = int( substr ($line, 22, 4) ) ;

			#next if ($resid<$firstresid+$cut3N);
			next if ($resid<$ntrunc+$cut3N);

			if (substr ($line, 21, 1) eq "B") {
				substr ($line, 21, 1) = "D";
			} elsif (substr ($line, 21, 1) eq "C") {
				substr ($line, 21, 1) = "E";
			}
			substr ($line, 30, 8) = sprintf ("%8.3f", $rX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $rX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $rX->[2]);
			substr ($line, 22, 4) = sprintf "%4d", $resid+2000;

			print OUTPDB $line."\n";
		}
	}
	close (OUTPDB);
}

##
##  converting for our representation to ideal
##
sub get_A_offset_fromZ {
	my ($Z) = @_;

	my ($A,$offset);
	if ($symmtype eq "I213") {
		$A = abs(4*sqrt(2)*$Z);
		$offset = [-$Z*sqrt(1/2),-$Z*sqrt(1/2),-$Z*sqrt(1/2)];
	} elsif ($symmtype eq "P23" && $Z>0) {
		$A = abs($Z*sqrt(2));
		$offset = [$A*(1/4),$A*(1/4),$A*(1/4)];
	} elsif ($symmtype eq "P4132" && $Z>0) {
		$A = abs(4*sqrt(2)*$Z);
		$offset = [$A*(1/8),$A*(1/8),$A*(1/8)];
	} elsif ($symmtype eq "P4132" && $Z<0) {
		$A = abs(4*sqrt(2)*$Z);
		$offset = [$A*(3/8),$A*(3/8),$A*(3/8)];
	} else {
		$A = 0;
		$offset = [0,0,0];
	}

	return ($A,$offset);
}

##
##
##
sub get_symmops {
	my ($symmtype, $Z) = @_;
	my ($symmRs,$symmTs) = ([],[]);

	if ($symmtype eq "I213") {
		$symmRs = [
			[ [1,0,0] , [0,1,0] , [0,0,1] ] ,
			[ [-1,0,0] , [0,-1,0] , [0,0,1] ] ,
			[ [1,0,0] , [0,-1,0] , [0,0,-1] ] ,
			[ [-1,0,0] , [0,1,0] , [0,0,-1] ] ,
			[ [0,0,1] , [1,0,0] , [0,1,0] ] ,
			[ [0,0,-1] , [-1,0,0] , [0,1,0] ] ,
			[ [0,0,1] , [-1,0,0] , [0,-1,0] ] ,
			[ [0,0,-1] , [1,0,0] , [0,-1,0] ] ,
			[ [0,1,0] , [0,0,1] , [1,0,0] ] ,
			[ [0,1,0] , [0,0,-1] , [-1,0,0] ] ,
			[ [0,-1,0] , [0,0,1] , [-1,0,0] ] ,
			[ [0,-1,0] , [0,0,-1] , [1,0,0] ] ,
			[ [1,0,0] , [0,1,0] , [0,0,1] ] ,
			[ [-1,0,0] , [0,-1,0] , [0,0,1] ] ,
			[ [1,0,0] , [0,-1,0] , [0,0,-1] ] ,
			[ [-1,0,0] , [0,1,0] , [0,0,-1] ] ,
			[ [0,0,1] , [1,0,0] , [0,1,0] ] ,
			[ [0,0,-1] , [-1,0,0] , [0,1,0] ] ,
			[ [0,0,1] , [-1,0,0] , [0,-1,0] ] ,
			[ [0,0,-1] , [1,0,0] , [0,-1,0] ] ,
			[ [0,1,0] , [0,0,1] , [1,0,0] ] ,
			[ [0,1,0] , [0,0,-1] , [-1,0,0] ] ,
			[ [0,-1,0] , [0,0,1] , [-1,0,0] ] ,
			[ [0,-1,0] , [0,0,-1] , [1,0,0] ]
		];
		$symmTs = [
			[0,0,0] ,
			[0,0.5,0] ,
			[0,0,0.5] ,
			[0,0.5,0.5] ,
			[0,0,0] ,
			[0,0.5,0] ,
			[0,0,0.5] ,
			[0,0.5,0.5] ,
			[0,0,0] ,
			[0,0,0.5] ,
			[0,0.5,0.5] ,
			[0.5,0,0.5],
			[0.5,0.5,0.5] ,
			[0.5,0,0.5] ,
			[0.5,0.5,0] ,
			[0.5,0,0] ,
			[0.5,0.5,0.5] ,
			[0.5,0,0.5] ,
			[0.5,0.5,0] ,
			[0.5,0,0] ,
			[0.5,0.5,0.5] ,
			[0.5,0.5,0] ,
			[0.5,0,0] ,
			[0,0.5,0]
		];
	} elsif ($symmtype eq "P23") { # P 2 3
		$symmRs = [
			[ [1,0,0] , [0,1,0] , [0,0,1] ],
			[ [-1,0,0] , [0,-1,0] , [0,0,1] ],
			[ [1,0,0] , [0,-1,0] , [0,0,-1] ],
			[ [-1,0,0] , [0,1,0] , [0,0,-1] ],
			[ [0,0,1] , [1,0,0] , [0,1,0] ],
			[ [0,0,-1] , [-1,0,0] , [0,1,0] ],
			[ [0,0,1] , [-1,0,0] , [0,-1,0] ],
			[ [0,0,-1] , [1,0,0] , [0,-1,0] ],
			[ [0,1,0] , [0,0,1] , [1,0,0] ],
			[ [0,1,0] , [0,0,-1] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,1] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,-1] , [1,0,0] ]
		];
		$symmTs = [
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0],
			[0,0,0]
		];
	} elsif ($symmtype eq "P4132" && $Z>0) {  # P 41 3 2
		$symmRs = [
			[ [1,0,0] , [0,1,0] , [0,0,1] ],
			[ [0,-1,0] , [1,0,0] , [0,0,1] ],
			[ [-1,0,0] , [0,-1,0] , [0,0,1] ],
			[ [0,1,0] , [-1,0,0] , [0,0,1] ],
			[ [1,0,0] , [0,-1,0] , [0,0,-1] ],
			[ [0,1,0] , [1,0,0] , [0,0,-1] ],
			[ [-1,0,0] , [0,1,0] , [0,0,-1] ],
			[ [0,-1,0] , [-1,0,0] , [0,0,-1] ],
			[ [0,0,1] , [1,0,0] , [0,1,0] ],
			[ [-1,0,0] , [0,0,1] , [0,1,0] ],
			[ [0,0,-1] , [-1,0,0] , [0,1,0] ],
			[ [1,0,0] , [0,0,-1] , [0,1,0] ],
			[ [0,0,1] , [-1,0,0] , [0,-1,0] ],
			[ [1,0,0] , [0,0,1] , [0,-1,0] ],
			[ [0,0,-1] , [1,0,0] , [0,-1,0] ],
			[ [-1,0,0] , [0,0,-1] , [0,-1,0] ],
			[ [0,1,0] , [0,0,1] , [1,0,0] ],
			[ [0,1,0] , [0,0,-1] , [-1,0,0] ],
			[ [0,0,1] , [0,1,0] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,1] , [-1,0,0] ],
			[ [0,0,-1] , [0,-1,0] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,-1] , [1,0,0] ],
			[ [0,0,1] , [0,-1,0] , [1,0,0] ],
			[ [0,0,-1] , [0,1,0] , [1,0,0] ]
		];

		$symmTs = [
			[0,0,0],
			[0.25,0.75,0.25],
			[0.5,0,0.5],
			[0.25,0.25,0.75],
			[0.5,0.5,0],
			[0.75,0.25,0.25],
			[0,0.5,0.5],
			[0.75,0.75,0.75],
			[0,0,0],
			[0.25,0.75,0.25],
			[0.5,0,0.5],
			[0.25,0.25,0.75],
			[0.5,0.5,0],
			[0.75,0.25,0.25],
			[0,0.5,0.5],
			[0.75,0.75,0.75],
			[0,0,0],
			[0.5,0.5,0],
			[0.75,0.25,0.25],
			[0,0.5,0.5],
			[0.75,0.75,0.75],
			[0.5,0,0.5],
			[0.25,0.25,0.75],
			[0.25,0.75,0.25]
		];
	} elsif ($symmtype eq "P4132" && $Z<0) { # P 43 3 2
		$symmRs = [
			[ [1,0,0] , [0,1,0] , [0,0,1] ],
			[ [0,-1,0] , [1,0,0] , [0,0,1] ],
			[ [-1,0,0] , [0,-1,0] , [0,0,1] ],
			[ [0,1,0] , [-1,0,0] , [0,0,1] ],
			[ [1,0,0] , [0,-1,0] , [0,0,-1] ],
			[ [0,1,0] , [1,0,0] , [0,0,-1] ],
			[ [-1,0,0] , [0,1,0] , [0,0,-1] ],
			[ [0,-1,0] , [-1,0,0] , [0,0,-1] ],
			[ [0,0,1] , [1,0,0] , [0,1,0] ],
			[ [-1,0,0] , [0,0,1] , [0,1,0] ],
			[ [0,0,-1] , [-1,0,0] , [0,1,0] ],
			[ [1,0,0] , [0,0,-1] , [0,1,0] ],
			[ [0,0,1] , [-1,0,0] , [0,-1,0] ],
			[ [1,0,0] , [0,0,1] , [0,-1,0] ],
			[ [0,0,-1] , [1,0,0] , [0,-1,0] ],
			[ [-1,0,0] , [0,0,-1] , [0,-1,0] ],
			[ [0,1,0] , [0,0,1] , [1,0,0] ],
			[ [0,1,0] , [0,0,-1] , [-1,0,0] ],
			[ [0,0,1] , [0,1,0] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,1] , [-1,0,0] ],
			[ [0,0,-1] , [0,-1,0] , [-1,0,0] ],
			[ [0,-1,0] , [0,0,-1] , [1,0,0] ],
			[ [0,0,1] , [0,-1,0] , [1,0,0] ],
			[ [0,0,-1] , [0,1,0] , [1,0,0] ]
		];
		$symmTs = [
			[0,0,0],
			[0.75,0.25,0.75],
			[0.5,0,0.5],
			[0.75,0.75,0.25],
			[0.5,0.5,0],
			[0.25,0.75,0.75],
			[0,0.5,0.5],
			[0.25,0.25,0.25],
			[0,0,0],
			[0.75,0.25,0.75],
			[0.5,0,0.5],
			[0.75,0.75,0.25],
			[0.5,0.5,0],
			[0.25,0.75,0.75],
			[0,0.5,0.5],
			[0.25,0.25,0.25],
			[0,0,0],
			[0.5,0.5,0],
			[0.25,0.75,0.75],
			[0,0.5,0.5],
			[0.25,0.25,0.25],
			[0.5,0,0.5],
			[0.75,0.75,0.25],
			[0.75,0.25,0.75]
		];
	} else {
		## point symm
		my $symmops_out;
		if ($symmtype eq 'T') {
			$symmops_out = [
				[[1,0,0],[0,1,0],[0,0,1]],
				[[-0.1667,0.2887,0.9428],[0.8660,0.5000,0],[-0.4714,0.8165,-0.3333]],
				[[-0.1667,0.8660,-0.4714],[0.2887,0.5000,0.8165],[0.9428,0.0000,-0.3333]],
				[[-0.1667 ,  -0.8660  , -0.4714],[ -0.2887 ,   0.5000 ,  -0.8165],[ 0.9428  , -0.0000  , -0.3333]]
			];
		} elsif ($symmtype eq 'O') {
			$symmops_out = [
				[[1,0,0],[0,1,0],[0,0,1]],
				[[0.16667,-0.28868,0.94281],[0.86603,0.50000,0.00000],[-0.47140,0.81650,0.33333]],
				[[ -6.6667e-01,5.7735e-01,4.7140e-01],[5.7735e-01,4.7184e-16,8.1650e-01],[4.7140e-01,8.1650e-01,-3.3333e-01]],
				[[ 1.6667e-01,8.6603e-01,-4.7140e-01],[-2.8868e-01,5.0000e-01,8.1650e-01],[9.4281e-01,-2.7756e-17,3.3333e-01]],
				[[ 0.50000,-0.86603,0.00000],[-0.86603 , -0.50000 , -0.00000],[0.00000  , 0.00000  ,-1.00000]],
				[[ 0.33333  , 0.00000 , -0.94281],[ 0.00000,  -1.00000  ,-0.00000],[-0.94281 ,  0.00000 , -0.33333]],
				[[ -0.83333 ,  0.28868  ,-0.47140],[ 0.28868 , -0.50000 , -0.81650],[-0.47140 , -0.81650  , 0.33333]],
				[[ -6.6667e-01 , -5.7735e-01 ,  4.7140e-01],[-5.7735e-01 , -4.1633e-16  ,-8.1650e-01],[ 4.7140e-01 , -8.1650e-01  ,-3.3333e-01]]
			];
		} elsif( $symmtype eq 'I') {
			$symmops_out = [
				[[1,0,0],[0,1,0],[0,0,1]],
				[[3.7268e-01,6.4550e-01,6.6667e-01],[-8.6603e-01,5.0000e-01,-4.3698e-17],[-3.3333e-01,-5.7735e-01,7.4536e-01]],
				[[-0.64235,0.17841,0.74536],[-0.75576,-0.30902,-0.57735],[0.12732,-0.93417,0.33333]],
				[[-0.64235,-0.75576,0.12732],[0.17841,-0.30902,-0.93417],[0.74536,-0.57735,0.33333]],
				[[3.7268e-01,-8.6603e-01,-3.3333e-01],[6.4550e-01,5.0000e-01,-5.7735e-01],[6.6667e-01,3.7817e-16,7.4536e-01]],
				[[0.56366,-0.75576,-0.33333],[0.75576,0.30902,0.57735],[-0.33333,-0.57735,0.74536]],
				[[0.47568,-0.46709,0.74536],[-0.11026,0.80902,0.57735],[-0.87268,-0.35682,0.33333]],
				[[-0.47568,0.11026,0.87268],[-0.46709,0.80902,-0.35682],[-0.74536,-0.57735,-0.33333]],
				[[-0.97568,0.17841,-0.12732],[ 0.17841,0.30902,-0.93417],[-0.12732,-0.93417,-0.33333]],
				[[-3.3333e-01,-3.5682e-01,-8.7268e-01],[ 9.3417e-01,2.7756e-16,-3.5682e-01],[1.2732e-01,-9.3417e-01,3.3333e-01]],
				[[-0.14235,-0.46709,-0.87268],[0.46709,-0.80902,0.35682],[-0.87268,-0.35682,0.33333]],
				[[-3.3333e-01,-9.3417e-01,1.2732e-01],[3.5682e-01,1.2768e-15,9.3417e-01],[-8.7268e-01,3.5682e-01,3.3333e-01]],
				[[-0.47568,-0.11026,0.87268],[0.46709,0.80902,0.35682],[-0.74536,0.57735,-0.33333]],
				[[-3.7268e-01,8.6603e-01,3.3333e-01],[6.4550e-01,5.0000e-01,-5.7735e-01],[-6.6667e-01,-1.8874e-15,-7.4536e-01]],
				[[-0.16667,0.64550,-0.74536],[0.64550,-0.50000,-0.57735],[-0.74536,-0.57735,-0.33333]],
				[[0.64235,-0.75576,-0.12732],[-0.17841,-0.30902,0.93417],[-0.74536,-0.57735,-0.33333]],
				[[-3.7268e-01,-8.6603e-01,3.3333e-01],[ -6.4550e-01,5.0000e-01,5.7735e-01],[-6.6667e-01,-9.9920e-16,-7.4536e-01]],
				[[-1.0000e+00,-6.9389e-16,-2.7756e-16],[-6.1062e-16,1.0000e+00,-3.3307e-16],[3.0531e-16,-7.7716e-16,-1.0000e+00]],
				[[-3.7268e-01,6.4550e-01,-6.6667e-01],[8.6603e-01,5.0000e-01,-5.5511e-17],[3.3333e-01,-5.7735e-01,-7.4536e-01]],
				[[0.64235,0.17841,-0.74536],[ 0.75576 , -0.30902,0.57735],[ -0.12732,-0.93417,-0.33333]]
			];
		} else {
			print STDERR "UNDEFINED! s=$symmtype z=$Z\n";
			die;
		}

		my $symmops_in = [
			[[1,0,0],[0,1,0],[0,0,1]],
			[[ -0.5000,-0.8660,0],[0.8660,-0.5000,0],[0,0,1.0000]],
			[[-0.5000,0.8660,0],[-0.8660,-0.5000,0],[0,0,1.0000]]
		];

		foreach my $symmop1 (@{$symmops_out}) {
			foreach my $symmop2 (@{$symmops_in}) {
				my $symmop = mmult( $symmop1,$symmop2 );
				push @{ $symmRs }, $symmop;
				push @{ $symmTs }, [0,0,0];
			}
		}
	}

	return ($symmRs,$symmTs);
}

###
###
###
sub build_complex_ca {
	my ($all_cas1 , $all_cas2, $all_cas3,$Z) = @_;
	my (@complex);

	# ignore termini!
	foreach my $i (0..scalar(@{$all_cas1})-1) {
		my $cashift = $all_cas1->[$i];
		push @complex, [$cashift->[0],$cashift->[1],$cashift->[2],1];
	}
	foreach my $i (0..scalar(@{$all_cas2})-1) {
		my $cashift = $all_cas2->[$i];
		push @complex, [$cashift->[0],$cashift->[1],$cashift->[2],2];
	}
	foreach my $i (0..scalar(@{$all_cas3})-1) {
		my $cashift = $all_cas3->[$i];
		push @complex, [$cashift->[0],$cashift->[1],$cashift->[2],3];
	}

	return (\@complex);
}

###
###
###
sub clash_check {
	my ($incas,$Z) = @_;

	my ($A,$offset) = get_A_offset_fromZ($Z);

	my $ncas = scalar( @{$incas} );
	my $nclash = 0;

	my ($symmRs,$symmTs) = get_symmops($symmtype,$Z);

	# self clashes
	foreach my $i( 0..$ncas-1 ) {
		my $fI = $incas->[$i];
		foreach my $j( $i+3..$ncas-1 ) {
			my ($fJ,$delfX,$d2);
			$fJ = $incas->[$j];
			next if ($fI->[3] != $fJ->[3] && ($j<$i+2*$TERMTRIM || $j<$TERMTRIM || $i<$TERMTRIM));
			$delfX = vsub($fI , $fJ);
			$d2 = len2( $delfX );
			if ($d2<$CLASHDIST*$CLASHDIST) {
				$nclash++;
				if ($nclash > $NCLASH) {
					return $nclash;
				}
			}
		}
	}

	# symm clashes
	foreach my $s (1..scalar(@{$symmRs})-1) {
		my $symmop = $symmRs->[$s];
		my $transop = $symmTs->[$s];
		foreach my $i( 0..$ncas-1 ) {
			my $fI;
			if ($A>0) { $fI = vscale(1/$A, $incas->[$i]); }
			else      { $fI = $incas->[$i]; }
			my $sfI = vadd( mapply($symmop, $fI), $transop);
			foreach my $j( 0..$ncas-1 ) {
				my ($fJ,$delfX,$d2);
				if ($A>0) { # cryst
					$fJ = vscale(1/$A, $incas->[$j]);
					$delfX = min_mod1(  vsub($sfI , $fJ) );
					$d2 = $A*$A*len2( $delfX );
				} else {    # cage
					$fJ = $incas->[$j];
					$delfX = vsub($sfI , $fJ);
					$d2 = len2( $delfX );
				}

				if ($d2<$CLASHDIST*$CLASHDIST) {
					$nclash++;
					if ($nclash > $NCLASH) {
						return $nclash;
					}
				}
			}
		}
	}


	return $nclash;
}

###
sub clash_check_verbose {
	my ($incas,$Z) = @_;

	my ($A,$offset) = get_A_offset_fromZ($Z);

	my $ncas = scalar( @{$incas} );
	my ($symmRs,$symmTs) = get_symmops($symmtype,$Z);

	my $retval = dclone( $incas );

	foreach my $s (1..scalar(@{$symmRs})-1) {
		my $symmop = $symmRs->[$s];
		my $transop = $symmTs->[$s];
		foreach my $i( 0..$ncas-1 ) {
			if ($A>0) {
				my $fI = vscale(1/$A, $incas->[$i]);
				my $sfI = vadd( mapply($symmop, $fI), $transop);

				my $fJ = vscale(1/$A, $incas->[int($ncas/2)]);
				my $delfX = min_mod1(  vsub($sfI , $fJ) );
				my $newX = vscale( $A, vadd($fJ,$delfX) );
				push @{ $retval }, $newX;
			} else {
				push @{ $retval }, mapply($symmop, $incas->[$i]);
			}
		}
	}

	return $retval;
}

sub min_mod1 {
	my ($X) = @_;
	my $Y = [fmod($X->[0],1.0) , fmod($X->[1],1.0) , fmod($X->[2],1.0) ];
	$Y->[0] += 1 if ($Y->[0]<-0.5); $Y->[1] += 1 if ($Y->[1]<-0.5); $Y->[2] += 1 if ($Y->[2]<-0.5);
	$Y->[0] -= 1 if ($Y->[0]>0.5); $Y->[1] -= 1 if ($Y->[1]>0.5); $Y->[2] -= 1 if ($Y->[2]>0.5);
	return $Y;
}


###
###
###
sub transform_ca {
	my ($incas,$R,$T1,$T2,$Z) = @_;

	my ($A,$offset) = get_A_offset_fromZ($Z);

	my @xform;
	foreach my $ca (@{ $incas }) {
		my $X0 = vadd( $offset , mapply( $Rglobal, vadd( mapply($R,vsub($ca,$T1)) , $T2 ) ) );
		push @xform, [$X0->[0],$X0->[1],$X0->[2],$ca->[3]];
	}
	return \@xform;
}


# my ($R)=R2quat($X,$Y,$Z,$W)
sub quat2R {
	my ($X,$Y,$Z,$W) = @_;
	my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
	my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
	my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
	my $R = [ [ 1 - 2 * ( $yy+$zz ) ,     2 * ( $xy-$zw ) ,     2 * ( $xz+$yw ) ] ,
	          [     2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,     2 * ( $yz-$xw ) ] ,
	          [     2 * ( $xz-$yw ) ,     2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
	return $R;
}

# my ($R)=R2quat($X,$Y,$Z,$W)
sub quatnorm {
	my ($X,$Y,$Z,$W) = @_;
	my $S = sqrt( $X*$X+$Y*$Y+$Z*$Z+$W*$W );
	return [ $X/$S , $Y/$S , $Z/$S , $W/$S ];
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
		exit 1;
	}

	$Minv->[0][0] =  ($M->[1][1]*$M->[2][2]-$M->[1][2]*$M->[2][1])/$D;
	$Minv->[0][1] = -($M->[0][1]*$M->[2][2]-$M->[0][2]*$M->[2][1])/$D;
	$Minv->[0][2] =  ($M->[0][1]*$M->[1][2]-$M->[0][2]*$M->[1][1])/$D;
	$Minv->[1][0] = -($M->[1][0]*$M->[2][2]-$M->[2][0]*$M->[1][2])/$D;
	$Minv->[1][1] =  ($M->[0][0]*$M->[2][2]-$M->[0][2]*$M->[2][0])/$D;
	$Minv->[1][2] = -($M->[0][0]*$M->[1][2]-$M->[0][2]*$M->[1][0])/$D;
	$Minv->[2][0] =  ($M->[1][0]*$M->[2][1]-$M->[2][0]*$M->[1][1])/$D;
	$Minv->[2][1] = -($M->[0][0]*$M->[2][1]-$M->[0][1]*$M->[2][0])/$D;
	$Minv->[2][2] =  ($M->[0][0]*$M->[1][1]-$M->[0][1]*$M->[1][0])/$D;

	return $Minv;
}

sub norm {
	my ($u) = @_;
	my $unorm = sqrt( $u->[0]*$u->[0] + $u->[1]*$u->[1] + $u->[2]*$u->[2] );
	return [$u->[0]/$unorm, $u->[1]/$unorm, $u->[2]/$unorm];
}

sub len {
	my ($u) = @_;
	my $unorm = sqrt( $u->[0]*$u->[0] + $u->[1]*$u->[1] + $u->[2]*$u->[2] );
	return $unorm;
}

sub len2 {
	my ($u) = @_;
	my $unorm = ( $u->[0]*$u->[0] + $u->[1]*$u->[1] + $u->[2]*$u->[2] );
	return $unorm;
}

sub vsub {
	my ($b,$c) = @_;
	my $a = [ $b->[0] - $c->[0] ,
			  $b->[1] - $c->[1] ,
			  $b->[2] - $c->[2] ];
	return $a;
}

sub vadd {
	my ($b,$c) = @_;
	my $a = [ $b->[0] + $c->[0] ,
			  $b->[1] + $c->[1] ,
			  $b->[2] + $c->[2] ];
	return $a;
}

sub vdot {
	my ($b,$c) = @_;
	my $a = $b->[0]*$c->[0] + $b->[1]*$c->[1] + $b->[2]*$c->[2];
	return $a;
}

sub vangle {
	my ($u,$v) = @_;

	my $unorm = sqrt( $u->[0]*$u->[0] + $u->[1]*$u->[1] + $u->[2]*$u->[2] );
	my $vnorm = sqrt( $v->[0]*$v->[0] + $v->[1]*$v->[1] + $v->[2]*$v->[2] );
	my $udotv = $u->[0]*$v->[0] + $u->[1]*$v->[1] + $u->[2]*$v->[2];

	my $angle = acos( $udotv / ($unorm*$vnorm) );

	return $angle;
}

sub vcross {
        my ($b,$c) = @_;
        my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
                  $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
                  $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
        return $a;
}
sub vscale {
        my ($b,$c) = @_;
        my $a = [ $b*$c->[0],  $b*$c->[1],  $b*$c->[2] ];
        return $a;
}

sub getSymmAxis {
	my $pdbfile = $_[0];
	my $alnchain = $_[1];
	$pdbfile =~ s/\.pdb/_0001.pdb/;

	# get calpha traces
	my %chains;
	my ($massSum, $natoms) = ([0,0,0],0);

	open (PDB, $pdbfile) || die "Cannot open $pdbfile";
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
				$massSum = vadd( $massSum, $CA_i ); $natoms++;
			}
		}
	}
	close (PDB);

	$massSum = vscale( 1.0/$natoms, $massSum);

	# align A to B, get axis of rotation
	my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{"A"} , $chains{$alnchain} );
	my ($X,$Y,$Z,$W)=R2quat($R);

	my $symmaxis = norm([$X,$Y,$Z]);

	return ($symmaxis,$massSum);
}


###
###
###
sub getHelixVector {
	my ($pdbfile, $tgtchn, $termhelix, $norc) = @_;

	# read pdb file
	open (FILE, $pdbfile) || die "$0: unable to open file $pdbfile for reading";
	my @pdb_buf = <FILE>;
    close (FILE);

	my $helix_vectors=[];
	my $helix_ids=[];
	my $helix_cas=[];
	my $helix_aas=[];
	my $all_cas=[];

	my %resmask;
	my $counter = 0;
	my $chainbreak = -1;

	# get all helix residues
	foreach my $line (@pdb_buf) {
		next if ($line !~ /^ATOM/);

		my $resid = int( substr ($line, 22, 4) );
		my $atom = substr ($line, 13, 3);
		next if ($atom ne "CA " && $atom ne "C  " && $atom ne "N  " && $atom ne "O  ");

		my $chn = substr ($line, 21, 1);
		next if ($chn ne 'A');

		my ($x,$y,$z) = (substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8));

		if ( ($norc eq 'N' && $resid >= $termhelix->[0]) || ($norc eq 'C' && $resid <= $termhelix->[1]) ) {
			push @{ $all_cas }, [$x,$y,$z,$chn] if ($atom eq "CA ");
		}

		next if ($chn ne $tgtchn);
		if ( $resid >= $termhelix->[0] && $resid <= $termhelix->[1]) {
			push @{ $helix_aas }, [$x,$y,$z];
			push @{ $helix_cas }, [$x,$y,$z] if ($atom eq "CA ")
		}
	}
	if (scalar(@{ $helix_aas }) < 24) {
		die "ERR at $pdbfile\n";
	}

	my ($helixvec,$helixend) = getHelixVectorParsed( $helix_aas, $norc );
	return ($helixvec,$helix_cas,$helix_aas,$all_cas,$chainbreak);
}

sub getHelixVectorParsed {
	my ($helix_cas, $norc) = @_;

	if (scalar( @{$helix_cas} ) < 24) { die "ERROR: ".scalar( @{$helix_cas} )." <=24 in getHelixVectorParsed!\n"; }

	my (@ca_0,@ca_1,@ca_2,@ca_3,@ca_4);
 	if ($norc eq 'N') {
 		@ca_0 = @{ $helix_cas->[5] }; # CA 2
 		@ca_1 = @{ $helix_cas->[9] };
 		@ca_2 = @{ $helix_cas->[13] };
 		@ca_3 = @{ $helix_cas->[17] };
 		@ca_4 = @{ $helix_cas->[21] };
 	} else {
 		@ca_0 = @{ $helix_cas->[-23] };
 		@ca_1 = @{ $helix_cas->[-19] };
 		@ca_2 = @{ $helix_cas->[-15] };
 		@ca_3 = @{ $helix_cas->[-11] };
 		@ca_4 = @{ $helix_cas->[-7] };
 	}

	# now do calc
	my @ca_02 = ( 0.5*($ca_2[0]+$ca_0[0]) , 0.5*($ca_2[1]+$ca_0[1]) , 0.5*($ca_2[2]+$ca_0[2]) );
	my @ca_1_tomid = ( $ca_02[0]-$ca_1[0] , $ca_02[1]-$ca_1[1] , $ca_02[2]-$ca_1[2] );
	my $ca_1_tomid_norm = sqrt ( $ca_1_tomid[0]*$ca_1_tomid[0] + $ca_1_tomid[1]*$ca_1_tomid[1] + $ca_1_tomid[2]*$ca_1_tomid[2] );
	$ca_1_tomid[0]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
	$ca_1_tomid[1]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
	$ca_1_tomid[2]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;

	my @ca_24 = ( 0.5*($ca_4[0]+$ca_2[0]) , 0.5*($ca_4[1]+$ca_2[1]) , 0.5*($ca_4[2]+$ca_2[2]) );
	my @ca_3_tomid = ( $ca_24[0]-$ca_3[0] , $ca_24[1]-$ca_3[1] , $ca_24[2]-$ca_3[2] );
	my $ca_3_tomid_norm = sqrt ( $ca_3_tomid[0]*$ca_3_tomid[0] + $ca_3_tomid[1]*$ca_3_tomid[1] + $ca_3_tomid[2]*$ca_3_tomid[2] );
	$ca_3_tomid[0]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
	$ca_3_tomid[1]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
	$ca_3_tomid[2]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;

	my $helixvec = norm( [ ($ca_3[0]+$ca_3_tomid[0]) - ($ca_1[0]+$ca_1_tomid[0]) ,
					       ($ca_3[1]+$ca_3_tomid[1]) - ($ca_1[1]+$ca_1_tomid[1]) ,
					       ($ca_3[2]+$ca_3_tomid[2]) - ($ca_1[2]+$ca_1_tomid[2]) ] );


	my $helixend = [ 0.5*(($ca_3[0]+$ca_3_tomid[0]) + ($ca_1[0]+$ca_1_tomid[0])) ,
					 0.5*(($ca_3[1]+$ca_3_tomid[1]) + ($ca_1[1]+$ca_1_tomid[1])) ,
					 0.5*(($ca_3[2]+$ca_3_tomid[2]) + ($ca_1[2]+$ca_1_tomid[2])) ];

 	if ($norc eq 'N') {
		$helixend = vadd( $helixend, vscale( -3*$HELIX_D_PER_RES, $helixvec ));  # ????
	} else {
		$helixend = vadd( $helixend, vscale(  3*$HELIX_D_PER_RES, $helixvec ));  # ????
	}

	#print "$norc:  ".$helixend->[0].",".$helixend->[1].",".$helixend->[2]
	#	."   ".$helix_cas->[0][0].",".$helix_cas->[0][1].",".$helix_cas->[0][2]
	#	."   ".$helix_cas->[-1][0].",".$helix_cas->[-1][1].",".$helix_cas->[-1][2]."\n";
	return ($helixvec,$helixend);
}


###
###
###
sub getTermHelixVectors {
	my ($pdbfile,$window) = @_;

	# read pdb file
	open (FILE, $pdbfile) || die "$0: unable to open file $pdbfile for reading";
	my @pdb_buf = <FILE>;
    close (FILE);

	my (@CAs_f,@AAs_f);
	foreach my $line (@pdb_buf) {
		next if ($line !~ /^ATOM/);
		my $resid = int( substr ($line, 22, 4) );
		my $atom = substr ($line, 13, 3);
		next if ($atom ne "CA " && $atom ne "N  " && $atom ne "C  " && $atom ne "O  ");

		my $x = substr ($line, 30, 8);
		my $y = substr ($line, 38, 8);
		my $z = substr ($line, 46, 8);
		my $chn = substr ($line, 21, 1);
		push @AAs_f, [$x,$y,$z,$chn];
		push @CAs_f, [$x,$y,$z,$chn] if ($atom eq "CA ");
	}
	my $nCAs_f = scalar( @CAs_f );
	my $nAAs_f = scalar( @AAs_f );

	# remove 'eat_into' at termini
	my @CAs=@CAs_f[  $EAT_INTO..$nCAs_f-  $EAT_INTO-1];
	my @AAs=@AAs_f[4*$EAT_INTO..$nAAs_f-4*$EAT_INTO-1];
	my $nCAs = scalar( @CAs );
	my $nAAs = scalar( @AAs );

	($nAAs == 4*$nCAs) || die "ERROR in input linker";

	my ($nterm,$cterm,$nterm_aa,$cterm_aa) = ([],[],[],[]);
	foreach my $i (0..$window-1) {
		push @{$nterm}, $CAs[$i];
		push @{$cterm}, $CAs[$nCAs-$window+$i];
		foreach my $j (0..3) {
			push @{$nterm_aa}, $AAs[4*($i)+$j];
			push @{$cterm_aa}, $AAs[$nAAs+4*(-$window+$i)+$j];
		}
	}

	# nterm vec
	my $helix_vectors=[];
	my $helix_midpts=[];
	{
		my @ca_0 = @{$CAs[1]};
		my @ca_1 = @{$CAs[2]};
		my @ca_2 = @{$CAs[3]};
		my @ca_3 = @{$CAs[4]};
		my @ca_4 = @{$CAs[5]};

		# now do calc
		my @ca_02 = ( 0.5*($ca_2[0]+$ca_0[0]) , 0.5*($ca_2[1]+$ca_0[1]) , 0.5*($ca_2[2]+$ca_0[2]) );
		my @ca_24 = ( 0.5*($ca_4[0]+$ca_2[0]) , 0.5*($ca_4[1]+$ca_2[1]) , 0.5*($ca_4[2]+$ca_2[2]) );

		my @ca_1_tomid = ( $ca_02[0]-$ca_1[0] , $ca_02[1]-$ca_1[1] , $ca_02[2]-$ca_1[2] );
		my $ca_1_tomid_norm = sqrt ( $ca_1_tomid[0]*$ca_1_tomid[0] + $ca_1_tomid[1]*$ca_1_tomid[1] + $ca_1_tomid[2]*$ca_1_tomid[2] );
		$ca_1_tomid[0]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
		$ca_1_tomid[1]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
		$ca_1_tomid[2]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;

		my @ca_3_tomid = ( $ca_24[0]-$ca_3[0] , $ca_24[1]-$ca_3[1] , $ca_24[2]-$ca_3[2] );
		my $ca_3_tomid_norm = sqrt ( $ca_3_tomid[0]*$ca_3_tomid[0] + $ca_3_tomid[1]*$ca_3_tomid[1] + $ca_3_tomid[2]*$ca_3_tomid[2] );
		$ca_3_tomid[0]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
		$ca_3_tomid[1]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
		$ca_3_tomid[2]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;

		my $helixvec = norm( [ ($ca_3[0]+$ca_3_tomid[0]) - ($ca_1[0]+$ca_1_tomid[0]) ,
						       ($ca_3[1]+$ca_3_tomid[1]) - ($ca_1[1]+$ca_1_tomid[1]) ,
						       ($ca_3[2]+$ca_3_tomid[2]) - ($ca_1[2]+$ca_1_tomid[2]) ] );
		push @{ $helix_vectors }, $helixvec;

		my $helixmid = [ 0.5*(($ca_3[0]+$ca_3_tomid[0]) + ($ca_1[0]+$ca_1_tomid[0])) ,
						 0.5*(($ca_3[1]+$ca_3_tomid[1]) + ($ca_1[1]+$ca_1_tomid[1])) ,
						 0.5*(($ca_3[2]+$ca_3_tomid[2]) + ($ca_1[2]+$ca_1_tomid[2])) ];
		push @{ $helix_midpts }, $helixmid;
	}
	# cterm vec
	{
		my @ca_0 = @{$CAs[$nCAs-6]};
		my @ca_1 = @{$CAs[$nCAs-5]};
		my @ca_2 = @{$CAs[$nCAs-4]};
		my @ca_3 = @{$CAs[$nCAs-3]};
		my @ca_4 = @{$CAs[$nCAs-2]};

		# now do calc
		my @ca_02 = ( 0.5*($ca_2[0]+$ca_0[0]) , 0.5*($ca_2[1]+$ca_0[1]) , 0.5*($ca_2[2]+$ca_0[2]) );
		my @ca_24 = ( 0.5*($ca_4[0]+$ca_2[0]) , 0.5*($ca_4[1]+$ca_2[1]) , 0.5*($ca_4[2]+$ca_2[2]) );

		my @ca_1_tomid = ( $ca_02[0]-$ca_1[0] , $ca_02[1]-$ca_1[1] , $ca_02[2]-$ca_1[2] );
		my $ca_1_tomid_norm = sqrt ( $ca_1_tomid[0]*$ca_1_tomid[0] + $ca_1_tomid[1]*$ca_1_tomid[1] + $ca_1_tomid[2]*$ca_1_tomid[2] );
		$ca_1_tomid[0]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
		$ca_1_tomid[1]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;
		$ca_1_tomid[2]  *= $HELIX_CA_RADIUS / $ca_1_tomid_norm;

		my @ca_3_tomid = ( $ca_24[0]-$ca_3[0] , $ca_24[1]-$ca_3[1] , $ca_24[2]-$ca_3[2] );
		my $ca_3_tomid_norm = sqrt ( $ca_3_tomid[0]*$ca_3_tomid[0] + $ca_3_tomid[1]*$ca_3_tomid[1] + $ca_3_tomid[2]*$ca_3_tomid[2] );
		$ca_3_tomid[0]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
		$ca_3_tomid[1]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;
		$ca_3_tomid[2]  *= $HELIX_CA_RADIUS / $ca_3_tomid_norm;

		my $helixvec = norm( [ ($ca_3[0]+$ca_3_tomid[0]) - ($ca_1[0]+$ca_1_tomid[0]) ,
						       ($ca_3[1]+$ca_3_tomid[1]) - ($ca_1[1]+$ca_1_tomid[1]) ,
						       ($ca_3[2]+$ca_3_tomid[2]) - ($ca_1[2]+$ca_1_tomid[2]) ] );
		push @{ $helix_vectors }, $helixvec;

		my $helixmid = [ 0.5*(($ca_3[0]+$ca_3_tomid[0]) + ($ca_1[0]+$ca_1_tomid[0])) ,
						 0.5*(($ca_3[1]+$ca_3_tomid[1]) + ($ca_1[1]+$ca_1_tomid[1])) ,
						 0.5*(($ca_3[2]+$ca_3_tomid[2]) + ($ca_1[2]+$ca_1_tomid[2])) ];
		push @{ $helix_midpts }, $helixmid;
	}

	return ($helix_vectors,$helix_midpts,$nterm_aa,$cterm_aa, \@CAs);
}

# my ($X,$Y,$Z,$W)=R2quat($M)
sub R2quat {
	my $R = shift;
	my ($S,$X,$Y,$Z,$W);
	if ( $R->[0][0] > $R->[1][1] && $R->[0][0] > $R->[2][2] )  {
		$S  = sqrt( 1.0 + $R->[0][0] - $R->[1][1] - $R->[2][2] ) * 2;
		$X = 0.25 * $S;
		$Y = ($R->[1][0] + $R->[0][1] ) / $S;
		$Z = ($R->[2][0] + $R->[0][2] ) / $S;
		$W = ($R->[2][1] - $R->[1][2] ) / $S;
	} elsif ( $R->[1][1] > $R->[2][2] ) {
		$S  = sqrt( 1.0 + $R->[1][1] - $R->[0][0] - $R->[2][2] ) * 2;
		$X = ($R->[1][0] + $R->[0][1] ) / $S;
		$Y = 0.25 * $S;
		$Z = ($R->[2][1] + $R->[1][2] ) / $S;
		$W = ($R->[0][2] - $R->[2][0] ) / $S;
	} else {
		$S  = sqrt( 1.0 + $R->[2][2] - $R->[0][0] - $R->[1][1] ) * 2;
		$X = ($R->[0][2] + $R->[2][0] ) / $S;
		$Y = ($R->[2][1] + $R->[1][2] ) / $S;
		$Z = 0.25 * $S;
		$W = ($R->[1][0] - $R->[0][1]) / $S;
	}
	return ($X,$Y,$Z,$W);
}


# my ($R,$rmsd, $Ycom, $Ycom_to_Xcom) = rms_align( $x,$y );
sub rms_align {
	my ($X,$Y) = @_;

	if (scalar( @{$X} ) != scalar( @{$Y} ) ) {die "rms_align mismatch!\n";}

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


###
###
###   kabsch alignment
###
###
# ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $ref_xlist, $mov_xlist )
sub setup_rotation {
	my ( $ref_xlist, $mov_xlist ) = @_;

	if (scalar(@{ $ref_xlist }) != scalar(@{ $mov_xlist }) ) {
		die "ERROR ".scalar(@{ $ref_xlist })." !=".scalar(@{ $mov_xlist })." in setup_rotation!\n";
	}

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
#    computes eigenval and eigen_vec of a real 3x3
#    symmetric matrix. On output, elements of a that are above
#    the diagonal are destroyed. d[1..3] returns the
#    eigenval of a. v[1..3][1..3] is a matrix whose
#    columns contain, on output, the normalized eigen_vec of a.
#    n_rot returns the number of Jacobi rotations that were required
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

	print STDERR "WARNING: Too many iterations in jacobi3!  You're bad and you should feel bad.\n";
	exit 1;
}


# ($eigen_vec, $eigenval) = diagonalize_symmetric( $matrix )
sub diagonalize_symmetric {
	my $matrix = shift;
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
	#   or right-handed. Let's force into right-handed system.
	$right_eigenvec->[2] = cross($right_eigenvec->[0], $right_eigenvec->[1]);

	# From the Kabsch algorithm, the eigenvec's of RtR
	#   are identical to the right_eigenvec's of R.
	#   This means that left_eigenvec = R x right_eigenvec
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
	#   orthogonal eigenvectors are in a right-handed or left-handed
	#   coordinate system - given by sigma. Sigma is needed to
	#   resolve this ambiguity in calculating the RMSD.
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


sub show_call_stack {
	my ( $path, $line, $subr );
	my $max_depth = 30;
	my $i = 1;
	print("--- Begin stack trace ---\n");
    while ( (my @call_details = (caller($i++))) && ($i<$max_depth) ) {
      print("$call_details[1] line $call_details[2] in function $call_details[3]\n");
    }
    print("--- End stack trace ---\n");
}
#######
# end #
#######

