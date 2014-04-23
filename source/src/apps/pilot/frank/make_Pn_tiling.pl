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
	print STDERR "Wallpaper symmdef file generation. Makes a wallpaper symmdef file from a point symmetric or wallpaper symmetric complex.\n";
	print STDERR "Supported symmetries:\n";
	print STDERR "       p2g  (input C2)                        p4   (input C4)\n";
	print STDERR "       c2m  (input D2)                        p4g  (input C4, with -z)\n";
	print STDERR "       p3   (input C3)                        p4m  (input D4)\n";
	print STDERR "       p31m (input C3 or D3, with -z)         p6   (input C6)\n";
	print STDERR "       p3m1 (input D3)                        p6m  (input D6)\n";
	print STDERR "usage: $0 [options]\n";
	print STDERR "example:   $0 -a A -i B C -x 22.6 -r 12.0 -p mystructure.pdb\n";
	print STDERR "options: \n";
	print STDERR "    -p <string> : Input PDB file\n";
	print STDERR "    -x <real>   : [default 100.0] distance between pt symm groups\n";
	print STDERR "    -a <char>   : [default A] the chain ID of the main chain\n";
	print STDERR "    -i <char>   : [default B] the chain IDs of point-group chain(s) (C2,C3,C4,C6 or D2,D3,D4,D6)\n";
	print STDERR "    -j <char>   : [default none] the chain IDs of a lattice chain\n";
	print STDERR "    -z          : IF -j IS NOT GIVEN, enable alternate wallpaper group from point symmetry\n";
	print STDERR "                :    P31m from C3/D3 or P4g from C4\n";
	print STDERR "    -q          : [default false] don't output PDBs\n";
	exit -1;
}

my $pdbfile;
my $ptgp_sep = 100.0;  # min interaction distance
my $primary_chain = 'A';
my @secondary_chains = ('B');
my $lattice_chain = '';
my $quietMode = 0;
my $mirroredGroupFlag = 0;

## parse options (do this by hand since Getopt does not handle this well)
my $inlinefull = (join ' ',@ARGV)." ";
my @suboptions = split /(-[a-z|A-Z] )/, $inlinefull;
for ( my $i=0; $i<=$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-x " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$ptgp_sep = ( $suboptions[++$i] );
		$ptgp_sep =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-a " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$primary_chain = $suboptions[++$i];
		$primary_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-i " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@secondary_chains = split /[, ]/,$suboptions[++$i];
	} elsif ($suboptions[$i] eq "-j " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$lattice_chain = $suboptions[++$i];
		$lattice_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
		$pdbfile =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] =~ /^-z/ ) {
		$mirroredGroupFlag = 1;
	} elsif ($suboptions[$i] =~ /^-q/ ) {
		$quietMode = 1;
		print STDERR "Not outputting PDBs.\n";
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

### '_' -> ' '
if ($primary_chain eq '_') {
	$primary_chain = ' ';
}
foreach my $i (0..$#secondary_chains) {
	if ($secondary_chains[$i] eq '_') {
		$secondary_chains[$i] = ' ';
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


# recenter the primary chain ...
# ... and add the identity op as root of the NCS op tree
if ( ! defined $chains{ $primary_chain } ) {
	die "Chain '$primary_chain' not in input!\n";
}

my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
my $COM_0 = recenter( $chains{ $primary_chain } );

# find residue closest to CoM of the system
# find the radius of the molecule
my $maxDist2 = 0;
my $minDist2 = 9999999;
my $minRes = 1;
foreach my $i ( 0..scalar( @{ $chains{ $primary_chain } })-1 ) {
	my $dist2 = vnorm2( $chains{ $primary_chain }->[$i] );
	if ($dist2 < $minDist2) {
		$minDist2 = $dist2;
		$minRes = $i+1;
	}
	if ($dist2 > $maxDist2) {
		$maxDist2 = $dist2;
	}
}
my $monomerRadius = sqrt( $maxDist2 );

my $Dsymm=0;
if ($#secondary_chains > 0) { $Dsymm = 1; }
if ($#secondary_chains > 1) { die "Too many secondary chains specified!\n"; }

## make sure 2ary & lattice chains are defined
foreach my $i (0..$#secondary_chains) {
	my @sec_chain_ids = split( ':', $secondary_chains[$i] );
	if ( ! defined $chains{ $sec_chain_ids[0] } ) {
		die "Chain ".$sec_chain_ids[0]." not in input!\n";
	}
	if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $sec_chain_ids[0] } } ) ) {
		print STDERR "ERROR! chains '$primary_chain' and '".$sec_chain_ids[0]."' have different residue counts! (".
		             scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $sec_chain_ids[0] } } ).")\n";
		die "Chain length mismatch!\n";
	}
}
if ( $lattice_chain ne '' && !defined $chains{ $lattice_chain } ) {
	die "Chain $lattice_chain not in input!\n";
}
if ( $lattice_chain ne '' &&
     scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $lattice_chain } } ) ) {
	print STDERR "ERROR! chains '$primary_chain' and '$lattice_chain' have different residue counts! (".
	             scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $lattice_chain } } ).")\n";
	die "Chain length mismatch!\n";
}

my (@Qs, @COMs, @sym_orders, @secondary_chains_filt);

## get symmetric transformations
foreach my $i (0..$#secondary_chains) {
	my @sec_chain_ids = split( ':', $secondary_chains[$i] );
	my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $sec_chain_ids[0] } );
	my $delCOM = vsub ($COM_i, $COM_0);
	my ($X,$Y,$Z,$W)=R2quat($R);
	my $omega = acos(abs($W));
	my $symm_order = int(PI/$omega + 0.5);
	if ($#sec_chain_ids > 0) { $symm_order = $sec_chain_ids[1]; }

	my $Wmult=1;
	if ($W<0) { $Wmult=-1; }
	my $newW = -$Wmult*cos( PI/$symm_order );
	my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );

	push @secondary_chains_filt, $sec_chain_ids[0];
	push @sym_orders, $symm_order;
	push @Qs, [$X*$S , $Y*$S, $Z*$S, $newW];
	push @COMs, $delCOM;
}

# expand D symmetries in proper order
if ($Dsymm && ($sym_orders[ 0 ] == 2 && $sym_orders[ 1 ] != 2) ) {
	my $temp;
	$temp = $Qs[1]; $Qs[1] = $Qs[0]; $Qs[0] = $temp;
	$temp = $COMs[1]; $COMs[1] = $COMs[0]; $COMs[0] = $temp;
	$temp = $sym_orders[1]; $sym_orders[1] = $sym_orders[0]; $sym_orders[0] = $temp;
	$temp = $secondary_chains_filt[1]; $secondary_chains_filt[1] = $secondary_chains_filt[0]; $secondary_chains_filt[0] = $temp;
}

# sanity check
if ($sym_orders[ 0 ] != 2 && $sym_orders[ 0 ] != 3 && $sym_orders[ 0 ] != 4 && $sym_orders[ 0 ] != 6 ) {
	die "Point symmetry of input structure must be of order 2,3,4 or 6 (detected ".$sym_orders[ 0 ].")\n";
}

# outerR
my $outer_sym_order = $sym_orders[ 0 ];
my $symm_group_name;

# P2gg,C2mm
if ($sym_orders[ 0 ] == 2) {
	if ($Dsymm==1) { $symm_group_name="c2m"; }
	else { $symm_group_name="p2g"; }
	$outer_sym_order = 4;
}

# P3, P3m1, and P31m from D3
if ($sym_orders[ 0 ] == 3) {
	if ($Dsymm==1 || $mirroredGroupFlag == 0) { $outer_sym_order = 6; }
	if ($Dsymm==1) { $symm_group_name="p3m1"; }
	elsif ($mirroredGroupFlag == 1) { $symm_group_name="p31m"; }
	else { $symm_group_name="p3"; }
}

if ($sym_orders[ 0 ] == 4) {
	if ($Dsymm==1) { $symm_group_name="p4m"; } 
	elsif ($mirroredGroupFlag == 1) { $symm_group_name="p4g"; }
	else { $symm_group_name="p4"; }
}

if ($sym_orders[ 0 ] == 6) {
	if ($Dsymm==1) { $symm_group_name="p6m"; } else { $symm_group_name="p6"; }
}
print STDERR "Building $symm_group_name lattice!\n";

###
### PROPERTIES
my ($mirroredGroup,$moveInPlane,$latticeRotation,$secondShell) = (0,0,0,0);
if ($symm_group_name eq "c2m" || $symm_group_name eq "p2g") { 
	$moveInPlane = 1;
	$secondShell = 1;
	$latticeRotation = PI/2;
}
if ($symm_group_name eq "p2g" || $symm_group_name eq "p31m" || $symm_group_name eq "p4g") { $mirroredGroup = 1; }

#if ($symm_group_name eq "p4m" || $symm_group_name eq "p3m1" || $symm_group_name eq "p6m") { $mirroredGroup = 1; } #fpd

if ($mirroredGroupFlag == 1 && $mirroredGroup==0) {
	print STDERR "Warning: -z flag not applicable for this point group.  Ignoring!\n";
}
print STDERR "mirroredGroup = $mirroredGroup\n";
print STDERR "moveInPlane = $moveInPlane\n";
print STDERR "latticeRotation = $latticeRotation\n";

# correct transformations
my $COM_pointgp = [0,0,0];
my $COM_Ccomplex = [0,0,0];

##  A->B0
if ($Dsymm == 0) {
	my $newR = quat2R( $Qs[0]->[0], $Qs[0]->[1], $Qs[0]->[2], $Qs[0]->[3] );
	my $err_pos = [0,0,0];
	my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (1..$sym_orders[0]) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$COMs[0] ) );
		$R_i = mmult($newR, $R_i);
	}
	$COMs[0] = vsub( $COMs[0] , [ $err_pos->[0]/$sym_orders[0] , $err_pos->[1]/$sym_orders[0] , $err_pos->[2]/$sym_orders[0] ] );

	$err_pos = [0,0,0];
	$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (1..$sym_orders[0]) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$COMs[0] ) );
		$COM_pointgp = vadd( $COM_pointgp, $err_pos );
		$R_i = mmult($newR, $R_i);
	}
	$COM_pointgp = vscale( 1.0/$sym_orders[0], $COM_pointgp );
	$COM_pointgp = vadd( $COM_0, $COM_pointgp );

	# add "noop" symmops
	push @secondary_chains_filt, '';
	push @sym_orders, 1;
	push @Qs, [0,0,0,1];
	push @COMs, [0,0,0];
} elsif ($Dsymm == 1) {
	my $X = [ $Qs[0]->[0],  $Qs[0]->[1],  $Qs[0]->[2] ];
	my $Y = [ $Qs[1]->[0],  $Qs[1]->[1],  $Qs[1]->[2] ];
	normalize($X); normalize($Y);
	my $Xtgt = vsub( $X , vscale(dot($X,$Y),$Y) );
	my $Ytgt = vsub( $Y , vscale(dot($X,$Y),$X) );
	normalize( $Xtgt ); normalize( $Ytgt );

	my $X0 = [ ($X->[0]+$Xtgt->[0])/2 , ($X->[1]+$Xtgt->[1])/2 , ($X->[2]+$Xtgt->[2])/2 ];
	my $Y0 = [ ($Y->[0]+$Ytgt->[0])/2 , ($Y->[1]+$Ytgt->[1])/2 , ($Y->[2]+$Ytgt->[2])/2 ];
	my $W_x = $Qs[0]->[3];
	my $W_y = $Qs[1]->[3];
	my $S_x = sqrt ( (1-$W_x*$W_x)/vnorm2($X0) );
	my $S_y = sqrt ( (1-$W_y*$W_y)/vnorm2($Y0) );

	$Qs[0] = [ $X0->[0]*$S_x , $X0->[1]*$S_x, $X0->[2]*$S_x, $W_x];
	$Qs[1] = [ $Y0->[0]*$S_y , $Y0->[1]*$S_y, $Y0->[2]*$S_y, $W_y];

	# fix gp 1 transform
	my $err_pos = [0,0,0];
	my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	my $newR = quat2R( $Qs[0]->[0], $Qs[0]->[1], $Qs[0]->[2], $Qs[0]->[3] );
	foreach my $j (1..$sym_orders[0]) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$COMs[0] ) );
		$R_i = mmult($newR, $R_i);
	}
	#print STDERR "   err_trans(1) = ".vnorm( $err_pos )."\n";
	$COMs[0] = vsub( $COMs[0] , [ $err_pos->[0]/$sym_orders[0] , $err_pos->[1]/$sym_orders[0] , $err_pos->[2]/$sym_orders[0] ] );

	# corrected CoM of C complex
	$err_pos = [0,0,0];
	$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (0..$sym_orders[0]-1) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$COMs[0] ) );
		$COM_Ccomplex = vadd( $COM_Ccomplex, $err_pos );
		$R_i = mmult($newR, $R_i);
	}
	$COM_Ccomplex = vscale( 1.0/$sym_orders[0], $COM_Ccomplex );

	# fix gp 2 transform
	$err_pos = [0,0,0];
	$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	$newR = quat2R( $Qs[1]->[0], $Qs[1]->[1], $Qs[1]->[2], $Qs[1]->[3] );
	my $com_secondary = vadd( $COMs[1], mapply( $newR , $COM_Ccomplex ) );

	my $newDelCOM = vsub ( $com_secondary , $COM_Ccomplex );

	foreach my $j (1..$sym_orders[1]) {
		$err_pos = vadd( $err_pos, mapply( $R_i, $newDelCOM ) );
		$R_i = mmult($newR, $R_i);
	}
	my $axis_proj_i =  [ $Qs[0]->[0],  $Qs[0]->[1],  $Qs[0]->[2] ];
	normalize( $axis_proj_i );
	my $del_COM_inplane = vsub( $newDelCOM , vscale(dot($newDelCOM,$axis_proj_i),$axis_proj_i) );
	$err_pos = vscale( $sym_orders[1] , $del_COM_inplane );

	#print STDERR "   err_trans(2) = ".vnorm( $err_pos )."\n";
	$COMs[1] = vsub( $COMs[1] , [ $err_pos->[0]/$sym_orders[1] , $err_pos->[1]/$sym_orders[1] , $err_pos->[2]/$sym_orders[1] ] );

	$err_pos = [0,0,0];
	$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (1..$sym_orders[1]) {
		$err_pos = vadd( $err_pos, mapply( $R_i, $newDelCOM ) );
		$COM_pointgp = vadd( $COM_pointgp, $err_pos );
		$R_i = mmult($newR, $R_i);
	}
	$COM_pointgp = vscale( 1.0/$sym_orders[1], $COM_pointgp );
	$COM_pointgp = vadd( $COM_Ccomplex, $COM_pointgp );
	$COM_pointgp = vadd( $COM_0, $COM_pointgp );
}
my $outerW = cos( PI/$outer_sym_order );
my $S = sqrt ( (1-$outerW*$outerW)/($Qs[0]->[0]*$Qs[0]->[0] + $Qs[0]->[1]*$Qs[0]->[1] + $Qs[0]->[2]*$Qs[0]->[2]) );
my $outerQ = [$Qs[0]->[0]*$S , $Qs[0]->[1]*$S, $Qs[0]->[2]*$S, $outerW];
my $outerR = quat2R( $outerQ->[0], $outerQ->[1], $outerQ->[2], $outerQ->[3] );

## get translation A->C
my $delCOM_lattice;
if ( $lattice_chain ne '' ) {
	my ($R_lattice,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $lattice_chain } );
 	$delCOM_lattice = vsub ($COM_i, $COM_0);
 	if (!is_identity($R_lattice)) {
 		print STDERR "     ".$R_lattice->[0][0].",".$R_lattice->[0][1].",".$R_lattice->[0][2]."\n";
 		print STDERR " R = ".$R_lattice->[1][0].",".$R_lattice->[1][1].",".$R_lattice->[1][2]."\n";
 		print STDERR "     ".$R_lattice->[2][0].",".$R_lattice->[2][1].",".$R_lattice->[2][2]."\n";
 		die "Lattice symmetry ($primary_chain,$lattice_chain) must not have rotation!\n";  ##?? recover from this instead?
 	}

	## correct this to be perpendicular to pointgroup symmaxis
	my $axis_proj_i =  [ $Qs[0]->[0],  $Qs[0]->[1],  $Qs[0]->[2] ];
	normalize( $axis_proj_i );
	$delCOM_lattice = vsub( $delCOM_lattice , vscale(dot($delCOM_lattice,$axis_proj_i),$axis_proj_i) );
 	$ptgp_sep = vnorm( $delCOM_lattice );
} else {
 	$delCOM_lattice = vsub($COM_0, $COM_pointgp);   # point towards the first subunit

	## correct this to be perpendicular to pointgroup symmaxis
	my $axis_proj_i =  [ $Qs[0]->[0],  $Qs[0]->[1],  $Qs[0]->[2] ];
	normalize( $axis_proj_i );
	$delCOM_lattice = vsub( $delCOM_lattice , vscale(dot($delCOM_lattice,$axis_proj_i),$axis_proj_i) );
 	normalize( $delCOM_lattice );

 	$delCOM_lattice = vscale( $ptgp_sep, $delCOM_lattice );
}

###
###  compute the individual symmops
###
my %Rs;
my %Ts;
my %ptgps;
$ptgps{0} = $COM_pointgp; # com of cplx

## 1: expand point group 0
my $R_i = [[1,0,0],[0,1,0],[0,0,1]];
my $T_i = deep_copy( $COM_0 );
my $R0 = quat2R( $Qs[0]->[0], $Qs[0]->[1], $Qs[0]->[2], $Qs[0]->[3] );
my $R1 = quat2R( $Qs[1]->[0], $Qs[1]->[1], $Qs[1]->[2], $Qs[1]->[3] );
foreach my $i (1..$sym_orders[1]) {
	foreach my $j (1..$sym_orders[0]) {
		my $id = "0_".$j."_".$i;
		$Rs{ $id } = $R_i;
		$Ts{ $id } = $T_i;
		$T_i = vadd( $T_i, mapply( $R_i,$COMs[0] ) );
		$R_i = mmult($R_i, $R0);
	}
	$T_i = vadd( $COM_0, $COMs[1] );
	$R_i = deep_copy($R1);
}

# "master" coordinate frame defined by central ring
my ($masterX, $masterY, $masterZ);
$masterZ = [ $Qs[0]->[0],  $Qs[0]->[1],  $Qs[0]->[2] ];
normalize( $masterZ );
$masterY = vsub( $Ts{ "0_2_1" } , $Ts{ "0_1_1" } );
normalize( $masterY );
$masterX = cross( $masterY , $masterZ );
normalize( $masterX );

## 2: expand lattice groups
$R_i = [[1,0,0],[0,1,0],[0,0,1]];
$T_i = vadd( $COM_pointgp , $delCOM_lattice );
foreach my $i (1..$outer_sym_order) {
	$ptgps{$i} = $T_i;
	$R_i = mmult($outerR, $R_i);
	$T_i = vadd( $COM_pointgp, mapply( $R_i,$delCOM_lattice ) );
}


## 3: expand point gps of lattice groups
my $latticeR = [[1,0,0],[0,1,0],[0,0,1]];
my $localLatticeR = [[1,0,0],[0,1,0],[0,0,1]];
if ($mirroredGroup == 1) {
	my $X_l = vsub( $ptgps{1} , $ptgps{0} );
	my $W_l = cos( PI/2 );
	my $S_l = sqrt ( (1-$W_l*$W_l)/vnorm2($X_l) );
	my $latticeR_i = quat2R( $S_l*$X_l->[0], $S_l*$X_l->[1], $S_l*$X_l->[2], $W_l );
	$latticeR = mmult( $latticeR_i, $latticeR );
}

if ($latticeRotation > 0) {
	my $X_l = deep_copy( $masterZ );
	my $W_l = cos( $latticeRotation/2 );
	my $S_l = sqrt ( (1-$W_l*$W_l)/vnorm2($X_l) );
	my $latticeR_i = quat2R( $S_l*$X_l->[0], $S_l*$X_l->[1], $S_l*$X_l->[2], $W_l );
	$latticeR = mmult( $latticeR_i, $latticeR );
	$localLatticeR = quat2R( 0,0,sqrt(1-$W_l*$W_l), $W_l );
}

foreach my $k (1..$outer_sym_order) {
	foreach my $i (1..$sym_orders[1]) {
		foreach my $j (1..$sym_orders[0]) {
			my $id = $k."_".$j."_".$i;
			my $R_ij = $Rs{ "0_".$j."_".$i };
			$Rs{ $id } = mmult($R_ij, $latticeR );
			my $toffset = mapply( $localLatticeR, vsub( $Ts{ "0_".$j."_".$i } , $ptgps{0}) );
			$Ts{ $id } = vadd( $ptgps{$k} , $toffset );
		}
	}
}

## 4: second shell lattice groups for c2m and p2g ... always a translation only
if ($secondShell == 1) {
	$ptgps{12} = vadd( vsub( $ptgps{1}, $ptgps{0} ) , ( $ptgps{2}) );
	$ptgps{23} = vadd( vsub( $ptgps{2}, $ptgps{0} ) , ( $ptgps{3}) );
	$ptgps{34} = vadd( vsub( $ptgps{3}, $ptgps{0} ) , ( $ptgps{4}) );
	$ptgps{41} = vadd( vsub( $ptgps{4}, $ptgps{0} ) , ( $ptgps{1}) );

	foreach my $k (12,23,34,41) {
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
				my $R_ij = 
				$Rs{ $id } = $Rs{ "0_".$j."_".$i };
				$Ts{ $id } = vadd( $ptgps{$k} , vsub( $Ts{ "0_".$j."_".$i } , $ptgps{0}) );
			}
		}
	}

}


##
# find the equation for the energy of the complex
# first find the 1->n interfaces that come in pairs
my @syminterfaces = keys %Rs;
my %energy_counter;
foreach my $interface (@syminterfaces) {
	$energy_counter{ $interface } = 1;
}
# delete self-interface(???)
delete ($energy_counter{ "0_1_1" });
 
## symm file gen
my $symmname = $pdbfile;
$symmname =~ s/\.pdb$//;
$symmname = $symmname."__".$symm_group_name;
print "symmetry_name $symmname\n";
print "E = 2*VRT0_1_1";
foreach my $complex (keys %energy_counter) {
	print " + ".$energy_counter{$complex}."*(VRT0_1_1:VRT".$complex.")";
}
print "\n";
print "anchor_residue COM\n";

######   XYZ
print "virtual_coordinates_start\n";
print "xyz VRT0  ".
			sprintf("%.6f,%.6f,%.6f", 1, 0, 0)."  ".
			sprintf("%.6f,%.6f,%.6f", 0, 1, 0)."  ".
			sprintf("%.6f,%.6f,%.6f", $COM_pointgp->[0]+1, $COM_pointgp->[1], $COM_pointgp->[2])."\n";

# central point group
foreach my $i (1..$sym_orders[1]) {
	foreach my $j (1..$sym_orders[0]) {
		my $id = "0_".$j."_".$i;

		# x points from origin to parent CoM
		my $myX = mapply( $Rs{ $id } , $masterX );
		my $myY = mapply( $Rs{ $id } , $masterY );

		if ($i == 1 && $j == 1) {
			print "xyz VRT0"."_ctrl  ".
						sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
						sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
						sprintf("%.6f,%.6f,%.6f", $COM_pointgp->[0], $COM_pointgp->[1], $COM_pointgp->[2])."\n";
		}
		print "xyz VRT$id  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $COM_pointgp->[0], $COM_pointgp->[1], $COM_pointgp->[2])."\n";
	}
}


my $latticeX = mapply( $localLatticeR, $masterX );
my $latticeY = mapply( $localLatticeR, $masterY );
my $latticeZ = mapply( $localLatticeR, $masterZ );



# redirection to second shell
if ($secondShell == 1) {
	my $myZ = $masterZ;
	foreach my $i (12,23,34,41) {
		# pointing to each point group
		my $parent_com = $COM_pointgp;
	
		my $myX = vsub( $ptgps{0} , $ptgps{$i} );
		normalize( $myX );
		my $myY = cross( $myZ, $myX );
		normalize( $myY );
		print "xyz VRT$i"."_outer  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	
		$parent_com = vscale( 0.5, vadd( $ptgps{0}, $ptgps{$i}) );
		print "xyz VRT$i"."_redir_inner  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	
		$parent_com = $ptgps{$i};
		print "xyz VRT$i"."_redir_outer  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	}

	# connect 12->1,2
	{
		my $parent_com = vscale( 0.5, vadd( $ptgps{0}, $ptgps{12}) );
		my $myX = vsub( $ptgps{1} , $parent_com ); normalize( $myX );
		my $myY = cross( $myZ, $myX ); normalize( $myY );
		print "xyz VRT12_redir_to_VRT1  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";

		print "xyz VRT1_redir  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $ptgps{1}->[0], $ptgps{1}->[1], $ptgps{1}->[2])."\n";

		$myX = vsub( $ptgps{2} , $parent_com ); normalize( $myX );
		$myY = cross( $myZ, $myX ); normalize( $myY );
		print "xyz VRT12_redir_to_VRT2  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";

		print "xyz VRT2_redir  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $ptgps{2}->[0], $ptgps{2}->[1], $ptgps{2}->[2])."\n";
	}

	# connect 34->3,4
	{
		my $parent_com = vscale( 0.5, vadd( $ptgps{0}, $ptgps{34}) );
		my $myX = vsub( $ptgps{3} , $parent_com ); normalize( $myX );
		my $myY = cross( $myZ, $myX ); normalize( $myY );
		print "xyz VRT34_redir_to_VRT3  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";

		print "xyz VRT3_redir  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $ptgps{3}->[0], $ptgps{3}->[1], $ptgps{3}->[2])."\n";


		$myX = vsub( $ptgps{4} , $parent_com ); normalize( $myX );
		$myY = cross( $myZ, $myX ); normalize( $myY );
		print "xyz VRT34_redir_to_VRT4  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";

		print "xyz VRT4_redir  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $ptgps{4}->[0], $ptgps{4}->[1], $ptgps{4}->[2])."\n";
	}

	foreach my $i (1,2,3,4,12,23,34,41) {
		my $parent_com = $ptgps{$i};
		my $myX = mapply( $Rs{ $i."_1_1" } , $masterX );
		my $myY = mapply( $Rs{ $i."_1_1" } , $masterY );
		print "xyz VRT$i"."_ctrl  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	}


} else {
	# redirection layer to outer point groups
	foreach my $i (1..$outer_sym_order) {
		# pointing to each point group
		my $parent_com = $COM_pointgp;
		my $myZ = $masterZ;
	
		if (($moveInPlane==1) && (($i % 2) == 0)) {
			$myZ = [-$masterZ->[0],-$masterZ->[1],-$masterZ->[2]];
		}
	
		my $myX = vsub( $ptgps{0} , $ptgps{$i} );
		normalize( $myX );
		my $myY = cross( $myZ, $myX );
		normalize( $myY );
		print "xyz VRT$i"."_outer  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	
		$parent_com = $ptgps{$i};
		print "xyz VRT$i"."_redir  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	
		my $myX = mapply( $Rs{ $i."_1_1" } , $masterX );
		my $myY = mapply( $Rs{ $i."_1_1" } , $masterY );
		print "xyz VRT$i"."_ctrl  ".
					sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	}
}


##
## finally expand outer point groups
##
foreach my $k (1..$outer_sym_order) {
 	my $parent_com = $ptgps{$k};
 
	foreach my $i (1..$sym_orders[1]) {
		foreach my $j (1..$sym_orders[0]) {
			my $id = $k."_".$j."_".$i;

			# x points from origin to parent CoM
			my $myX = mapply( $Rs{ $id } , $masterX );
			my $myY = mapply( $Rs{ $id } , $masterY );
	
			normalize( $myX );
			$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
			normalize( $myY );

			print "xyz VRT$id  ".
						sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
						sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
						sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
		}
	}
}


if ($secondShell == 1) {
	foreach my $k (12,23,34,41) {
		my $parent_com = $ptgps{$k};
	 
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
	
				# x points from origin to parent CoM
				my $myX = mapply( $Rs{ $id } , $masterX );
				my $myY = mapply( $Rs{ $id } , $masterY );
		
				normalize( $myX );
				$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
				normalize( $myY );
	
				print "xyz VRT$id  ".
							sprintf("%.6f,%.6f,%.6f", $myX->[0], $myX->[1], $myX->[2])."  ".
							sprintf("%.6f,%.6f,%.6f", $myY->[0], $myY->[1], $myY->[2])."  ".
							sprintf("%.6f,%.6f,%.6f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
			}
		}
	}
}


print "virtual_coordinates_stop\n";

######   connect
print "connect_virtual JUMP0 VRT0 VRT0_ctrl\n";
print "connect_virtual JUMP0_1_1 VRT0_ctrl VRT0_1_1\n";

# central point group
foreach my $i (1..$sym_orders[1]) {
	foreach my $j (1..$sym_orders[0]) {
		my $id = "0_".$j."_".$i;
		if ($j != 1 || $i != 1) {
			print "connect_virtual JUMP".$id." VRT0_1_1 VRT$id"."\n";
		}
		print "connect_virtual JUMP".$id."_to_subunit VRT$id"." SUBUNIT"."\n";
	}
}

# outer point groups
if ($secondShell == 1) {
	foreach my $k (12,23,34,41) {
		my $id1 = $k."_1_1";
		print "connect_virtual JUMP$k"."_to_outer  VRT0_ctrl VRT$k"."_outer"."\n";
		print "connect_virtual JUMP$k"."_to_redir_inner VRT$k"."_outer VRT$k"."_redir_inner"."\n";
		print "connect_virtual JUMP$k"."_to_redir_outer VRT$k"."_redir_inner VRT$k"."_redir_outer"."\n";
		print "connect_virtual JUMP$k"."_to_ctrl VRT$k"."_redir_outer VRT$k"."_ctrl"."\n";
		print "connect_virtual JUMP$id1 VRT$k"."_ctrl VRT$id1\n";
	}

	print "connect_virtual JUMP12_to_redir1 VRT12_redir_inner VRT12_redir_to_VRT1"."\n";
	print "connect_virtual JUMP12_to_redir2 VRT12_redir_inner VRT12_redir_to_VRT2"."\n";
	print "connect_virtual JUMP34_to_redir3 VRT34_redir_inner VRT34_redir_to_VRT3"."\n";
	print "connect_virtual JUMP34_to_redir4 VRT34_redir_inner VRT34_redir_to_VRT4"."\n";
	print "connect_virtual JUMP1_to_redir VRT12_redir_to_VRT1 VRT1_redir"."\n";
	print "connect_virtual JUMP2_to_redir VRT12_redir_to_VRT2 VRT2_redir"."\n";
	print "connect_virtual JUMP3_to_redir VRT34_redir_to_VRT3 VRT3_redir"."\n";
	print "connect_virtual JUMP4_to_redir VRT34_redir_to_VRT4 VRT4_redir"."\n";
	print "connect_virtual JUMP1_to_ctrl VRT1_redir VRT1_ctrl"."\n";
	print "connect_virtual JUMP2_to_ctrl VRT2_redir VRT2_ctrl"."\n";
	print "connect_virtual JUMP3_to_ctrl VRT3_redir VRT3_ctrl"."\n";
	print "connect_virtual JUMP4_to_ctrl VRT4_redir VRT4_ctrl"."\n";
	print "connect_virtual JUMP1_1_1 VRT1_ctrl VRT1_1_1\n";
	print "connect_virtual JUMP2_1_1 VRT2_ctrl VRT2_1_1\n";
	print "connect_virtual JUMP3_1_1 VRT3_ctrl VRT3_1_1\n";
	print "connect_virtual JUMP4_1_1 VRT4_ctrl VRT4_1_1\n";


	foreach my $k (1,2,3,4,12,23,34,41) {
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
				if ($j != 1 || $i != 1) {
					print "connect_virtual JUMP".$id." VRT".$k."_1_1 VRT$id"."\n";
				}
				print "connect_virtual JUMP".$id."_to_subunit VRT$id"." SUBUNIT"."\n";
			}
		}
	}

} else {
	foreach my $k (1..$outer_sym_order) {
		my $id1 = $k."_1_1";
		print "connect_virtual JUMP$k"."_to_outer  VRT0_ctrl VRT$k"."_outer"."\n";
		print "connect_virtual JUMP$k"."_to_redir VRT$k"."_outer VRT$k"."_redir"."\n";
		print "connect_virtual JUMP$k"."_to_ctrl VRT$k"."_redir VRT$k"."_ctrl"."\n";
		print "connect_virtual JUMP$id1 VRT$k"."_ctrl VRT$id1\n";
	
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
				if ($j != 1 || $i != 1) {
					print "connect_virtual JUMP".$id." VRT".$k."_1_1 VRT$id"."\n";
				}
				print "connect_virtual JUMP".$id."_to_subunit VRT$id"." SUBUNIT"."\n";
			}
		}
	}
}


## dofs
if ($mirroredGroup == 1) {
	print "set_dof JUMP0_1_1 angle_z(0:360) z[0;-5:5]\n";    # spin
} else {
	print "set_dof JUMP0_1_1 angle_z(0:360)\n";    # spin
}

if ($secondShell == 1) {
	print "set_dof JUMP12_to_redir_inner x($ptgp_sep)\n";
	print "set_dof JUMP23_to_redir_inner x($ptgp_sep)\n";

	# along x
	print "set_jump_group JUMPGROUP0 ";
	print "JUMP12_to_redir_inner JUMP12_to_redir_outer:2 JUMP34_to_redir_inner JUMP34_to_redir_outer:2\n";

	# along y
	print "set_jump_group JUMPGROUP1 ";
	print "JUMP23_to_redir_inner JUMP23_to_redir_outer:2 JUMP41_to_redir_inner JUMP41_to_redir_outer:2 ";
	print "JUMP1_to_redir JUMP2_to_redir JUMP3_to_redir JUMP4_to_redir\n";

} else {
	print "set_dof JUMP1_to_redir x($ptgp_sep)\n";   # lattice spacing

	print "set_jump_group JUMPGROUP1 ";
	foreach my $k (1..$outer_sym_order) {
		print "JUMP$k"."_to_redir ";
	}
	print "\n";
}

print "set_jump_group JUMPGROUP2 ";
foreach my $k (0..$outer_sym_order) {
	my $id1 = $k."_1_1";
	print "JUMP$id1 ";
}
if ($secondShell == 1) {
	foreach my $k (12,23,34,41) {
		my $id1 = $k."_1_1";
		print "JUMP$id1 ";
	}
}
print "\n";

print "set_jump_group JUMPGROUP3 ";
foreach my $k (0..$outer_sym_order) {
	foreach my $i (1..$sym_orders[1]) {
		foreach my $j (1..$sym_orders[0]) {
			my $id = $k."_".$j."_".$i;
			print "JUMP$id"."_to_subunit ";
		}
	}
}
if ($secondShell == 1) {
	foreach my $k (12,23,34,41) {
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
				print "JUMP$id"."_to_subunit ";
			}
		}
	}
}
print "\n";

if ($mirroredGroup == 1) {
	print "slide_type ORDERED_SEQUENTIAL\n";
	if ($moveInPlane == 1) {
		print "slide_order JUMP12_to_redir_inner JUMP23_to_redir_inner\n";
	} else {
		print "slide_order JUMP1_to_redir\n";
	}
} 

if ($quietMode != 1) {
	########################################
	## write output pdb
	########################################
	my $outpdb = $pdbfile;
	my $outmon = $pdbfile;
	my $outmdl = $pdbfile;
	
	if ($outpdb =~ /\.pdb$/) {
		$outpdb =~ s/\.pdb$/_symm.pdb/;
		$outmon =~ s/\.pdb$/_INPUT.pdb/;
	} else {
		$outpdb = $outpdb."_symm.pdb";
		$outmon = $outmon."_INPUT.pdb";
	}
	open (OUTPDB, ">$outpdb");
	open (OUTMON, ">$outmon");
	
	my $chnidx = 0;
	my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";
	foreach my $k (0..$outer_sym_order) {
		foreach my $i (1..$sym_orders[1]) {
			foreach my $j (1..$sym_orders[0]) {
				my $id = $k."_".$j."_".$i;
		
				foreach my $line (@filebuf) {
					my $linecopy = $line;
		
					my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
					my $X_0 = vsub($X,$COM_0);
					my $rX = vadd( mapply($Rs{$id}, $X_0) , $Ts{$id} );
		
					substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
					substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
					substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
					substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);
		
					print OUTPDB $linecopy."\n";
				}
		
				print OUTPDB "TER   \n";
				$chnidx++;
			}
		}
	}
	
	if ($secondShell == 1) {
		foreach my $k (12,23,34,41) {
			foreach my $i (1..$sym_orders[1]) {
				foreach my $j (1..$sym_orders[0]) {
					my $id = $k."_".$j."_".$i;
			
					foreach my $line (@filebuf) {
						my $linecopy = $line;
			
						my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
						my $X_0 = vsub($X,$COM_0);
						my $rX = vadd( mapply($Rs{$id}, $X_0) , $Ts{$id} );
			
						substr ($linecopy, 30, 8) = sprintf ("%8.3f", $rX->[0]);
						substr ($linecopy, 38, 8) = sprintf ("%8.3f", $rX->[1]);
						substr ($linecopy, 46, 8) = sprintf ("%8.3f", $rX->[2]);
						substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);
			
						print OUTPDB $linecopy."\n";
					}
			
					print OUTPDB "TER   \n";
					$chnidx++;
				}
			}
		}
	}

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
