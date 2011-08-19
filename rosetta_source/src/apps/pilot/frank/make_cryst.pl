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
require "symmops.pm";
require "matrix.pm";
require "kabsch.pm";

###############################################################################

if ($#ARGV < 0) {
	print STDERR "usage: $0 [options]\n";
	print STDERR "example:   $0 -r 12.0 -c 42.4 41.2 88.6 90.0 90.0 90.0 -s P 1 21 1 -p mystructure.pdb\n";
	print STDERR "options: \n";
	print STDERR "    -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
#	print STDERR "    -m (NCS|CRYST) : [default CRYST] generate NCS symm file or crystallographic symm file\n";
	print STDERR "    -r <real>   : [default 8.0] the max CA-CA distance between two interacting chains\n";
#	print STDERR "NCS-mode specific options:\n";
#	print STDERR "    -a <char>   : [default A] the chain ID of the main chain\n";
#	print STDERR "    -i <char>   : [default B] the chain IDs of one chain in each symmetric complex\n";
#	print STDERR "CRYST-mode specific options:\n";
	print STDERR "    -b <real>   : Approximate the protein as a ball of this radius (only if no '-p'!)\n";
	print STDERR "    -c <real>x6 : override the unit cell parameters in the PDB with these values\n";
	print STDERR "    -s <string> : override the spacegroup in the PDB with these values\n";
	exit -1;
}

my $pdbfile;
my $interact_dist = 8.0;  # min interaction distance
my $sphere_size = -1; # approx protein as a sphere this size
my @cell_new;
my $spacegp_new;
my $modestring = "CRYST";

## parse options (do this by hand since Getopt does not handle this well)
my @suboptions = split /(-[a-z|A-Z] )/, (join ' ',@ARGV);
for ( my $i=0; $i<$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		#$interact_dist = int( $suboptions[++$i] );
		$interact_dist = $suboptions[++$i];
		$interact_dist =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-b " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$sphere_size = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-c " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@cell_new = split /[, ]/,$suboptions[++$i];
	} elsif ($suboptions[$i] eq "-s " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$spacegp_new = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
		$pdbfile =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-m " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$modestring = $suboptions[++$i];
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

# parse mode string
my ($cryst_mode, $ncs_mode);
if ($modestring eq "CRYST" || $modestring eq "cryst") {
	$cryst_mode = 1; $ncs_mode = 0;
} elsif ($modestring eq "NCS" || $modestring eq "ncs") {
	$cryst_mode = 0; $ncs_mode = 1;
} else {
	print STDERR "Unrecognized mode string '$modestring'\n";
	exit -1;
}

print STDERR "Running in mode $modestring.\n";

###############################################################################

# crystinfo
my $spacegp = "P 1";
my ($gpid, $nsymm, $Rs, $Ts, $Cs, $cheshire);
my ($A, $B, $C, $alpha, $beta, $gamma) = (0,0,0,90,90,90);
my ($f2c,$c2f);

my @chaintrace;
my @filebuf;
my $CoM = [0,0,0];

if ( length ($pdbfile) > 0 ) {
	###
	### Input PDB file
	open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
	while (<PDB>) {
		chomp;
		if (/^CRYST1/) {
			$A = substr ($_,  6, 9);
			$B = substr ($_, 15, 9);
			$C = substr ($_, 24, 9);
			$alpha = substr ($_, 33, 7);
			$beta  = substr ($_, 40, 7);
			$gamma = substr ($_, 47, 7);
			$spacegp =  substr ($_, 55, 11);
		} elsif (/^ATOM/ || /^HETATM/) {
			my $atom = substr ($_, 12, 4);
			if ($atom eq " CA ") {
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
	if (length($spacegp_new) > 0) {
		$spacegp = $spacegp_new;
	}
	
	# initialize symmops
	($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp_lookup( $spacegp );
	($f2c,$c2f) = crystparams($A,$B,$C,$alpha,$beta,$gamma);
} else {
	###
	### Interaction radius only
	if ($sphere_size < 0) {
		print STDERR "Must provide an input pdb file with -p _or_ a nonnegative approximate radius with -b!\n";
		exit -1;
	} else {
		print STDERR "Warning! No input structure provided ... generating symmetry mates using interaction radius only!\n";
	}

	# override crystal params from cmd line
	if ($#cell_new >= 5) {
		($A,$B,$C,$alpha,$beta,$gamma) = @cell_new;
	}
	if (length($spacegp_new) > 0) {
		$spacegp = $spacegp_new;
	}
	
	# initialize symmops
	($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp_lookup( $spacegp );
	($f2c,$c2f) = crystparams($A,$B,$C,$alpha,$beta,$gamma);

	# fake atom
	my $fake_fX = [ ($cheshire->[0][0]+$cheshire->[0][1])/2 ,
	                ($cheshire->[1][0]+$cheshire->[1][1])/2 ,
	                ($cheshire->[2][0]+$cheshire->[2][1])/2 ];
	my $fake_cX = mapply( $f2c , $fake_fX );
	push @chaintrace, $fake_cX;
	my $fakePDBline = sprintf "ATOM      1  CA  ALA A   1     %7.3f %7.3f %7.3f  1.00  0.00           C  ", 
	                          $fake_cX->[0], $fake_cX->[1], $fake_cX->[2];
	push @filebuf, $fakePDBline;
	$CoM = $fake_cX;

	if ($sphere_size > 0) {
		my $plusZed = [0,0,$sphere_size];
		my $nBetaSteps = int( PI*$sphere_size / 8.0 )+1;
		my $counter = 2;
		foreach my $i (1..$nBetaSteps) {
			my $beta = PI*($i-0.5)/($nBetaSteps);
			my $nGammaSteps = int( PI*$sphere_size*sin( $beta )/ 4.0 )+1;

			foreach my $j (1..$nGammaSteps) {
				my $gamma = 2*PI*($j-0.5)/($nGammaSteps);

				my $rot_ij = euler( 0,$beta,$gamma );
				my $cX_ij = vadd( mapply( $rot_ij , $plusZed ) , $fake_cX );
				push @chaintrace, $cX_ij;
				$fakePDBline = sprintf "ATOM    %3d  H   ALA A   1     %7.3f %7.3f %7.3f  1.00  0.00           H  ", 
										  $counter++, $cX_ij->[0], $cX_ij->[1], $cX_ij->[2];
				push @filebuf, $fakePDBline;
			}
		}
	}
}

print STDERR "Read: $A x $B x $C x ($alpha , $beta , $gamma)   $spacegp\n";

my $nAtms = $#chaintrace+1;
$CoM = [ $CoM->[0]/$nAtms , $CoM->[1]/$nAtms , $CoM->[2]/$nAtms ];

# find residue closest to CoM of the system
my $minDist2 = 9999999;
my $minRes = 1;
foreach my $i (0..$nAtms) {
	my $dist2 = vnorm2( vsub( $CoM, $chaintrace[ $i ] ) );
	if ($dist2 < $minDist2) {
		$minDist2 = $dist2;
		$minRes = $i+1;
	}
}

# frac coords of chaintrace
my @fchaintrace;
foreach my $X_i (@chaintrace) {
	my $fX_i = mapply( $c2f , $X_i );
	push @fchaintrace , $fX_i;
}


# do the symmetric expansion
my $Ts_expand = [ 
	[-1,-1,-1],[-1,-1,0],[-1,-1,1],  [-1,0,-1],[-1,0,0],[-1,0,1],  [-1,1,-1],[-1,1,0],[-1,1,1],
    [0,-1,-1] ,[0,-1,0] ,[0,-1,1] ,  [0,0,-1] ,        ,[0,0,1] ,  [0,1,-1] ,[0,1,0] ,[0,1,1],
	[1,-1,-1] ,[1,-1,0] ,[1,-1,1] ,  [1,0,-1] ,[1,0,0] ,[1,0,1] ,  [1,1,-1] ,[1,1,0] ,[1,1,1] 
];

my %symminterface = ();
foreach my $j_symm (0..($nsymm-1)) {
	## TO DO: early stop
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
				#    another offset within the interaction radius
				#    but not closer than this offset for every (i,j) pair
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


# symmetric interfaces may come in pairs
# find these pairs and place in syminterfaces_paired like [0, 1, 1', 2, 2', ...]
my @syminterfaces = sort { $symminterface{$a} <=> $symminterface{$b} } keys %symminterface;
my @syminterfaces_paired;
my @syminterfaces_unpaired;
push @syminterfaces_unpaired, $syminterfaces[0];   # first the origin
my %paired = (0=>1);

OUTER: foreach my $i (1..$#syminterfaces-1) {
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
	push @syminterfaces_unpaired, $syminterfaces[$i];   # unpaired
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
print "anchor_residue $minRes\n";

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
	my $xyzline = "xyz VRT_".$symmkey;
	$xyzline =~ s/_-(\d)/_n\1/g;

	# X
	my $fX  = mapply( $c2f , [1,0,0] );
	my $sfX  = mapply($Rs->[$j_symm],$fX);
	my $sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	$xyzline = $xyzline." ".$string;

	# Y
	$fX  = mapply( $c2f , [0,1,0] );
	$sfX  = mapply($Rs->[$j_symm],$fX);
	$sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	$xyzline = $xyzline." ".$string;

	# orig
	$fX  = mapply( $c2f , vadd( $CoM, [1,1,1] ) ) ;
	$sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
	$sfX = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
	$sX = mapply( $f2c , $sfX );

	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	$xyzline = $xyzline." ".$string;
	print "$xyzline\n";

	# centers of mass
	# X
	my $xyzline = "xyz VRT_".$symmkey."_base";
	$xyzline =~ s/_-(\d)/_n\1/g;
	my $fX  = mapply( $c2f , [1,0,0] );
	my $sfX  = mapply($Rs->[$j_symm],$fX);
	my $sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	$xyzline = $xyzline." ".$string;

	# Y
	$fX  = mapply( $c2f , [0,1,0] );
	$sfX  = mapply($Rs->[$j_symm],$fX);
	$sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	$xyzline = $xyzline." ".$string;

	# orig
	$fX  = mapply( $c2f , $CoM );
	$sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
	$sfX = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
	$sX = mapply( $f2c , $sfX );

	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
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
if ($cheshire->[0][0] != $cheshire->[0][1] || $cheshire->[1][0] != $cheshire->[1][1] || $cheshire->[2][0] != $cheshire->[2][1] ) {
	print "set_dof JUMP_$syminterfaces_all[0]"."_to_com";
	if ($cheshire->[0][0] != $cheshire->[0][1]) { print " x"; }
	if ($cheshire->[1][0] != $cheshire->[1][1]) { print " y"; }
	if ($cheshire->[2][0] != $cheshire->[2][1]) { print " z"; }
	print "\n";
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
my $chnidx = 0;
my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";
foreach my $symmkey (@syminterfaces_all) {
	my ($j_symm,$shiftX,$shiftY,$shiftZ) = split '_',$symmkey;
	foreach my $line (@filebuf) {
		my $linecopy = $line;

		my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
		my $fX = mapply( $c2f , $X );
		my $sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
		my $sfX1 = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
		my $sX = mapply( $f2c , $sfX1 );

		substr ($linecopy, 30, 8) = sprintf ("%8.3f", $sX->[0]);
		substr ($linecopy, 38, 8) = sprintf ("%8.3f", $sX->[1]);
		substr ($linecopy, 46, 8) = sprintf ("%8.3f", $sX->[2]);
		substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);

		print OUTPDB $linecopy."\n";
	}
	print OUTPDB "TER   \n";
	$chnidx++;
}
close(OUTPDB);


exit 0;

