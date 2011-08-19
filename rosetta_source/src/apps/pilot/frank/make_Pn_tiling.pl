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
	print STDERR "example:   $0 -a A -i B C -x 22.6 -r 12.0 -p mystructure.pdb\n";
	print STDERR "options: \n";
	print STDERR "    -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
	print STDERR "    -r <real>   : [default 8.0] the max CA-CA distance between two interacting chains\n";
	print STDERR "    -x <real>   : [default 100.0] distance between pt symm groups\n";
	print STDERR "    -a <char>   : [default A] the chain ID of the main chain\n";
	print STDERR "    -i <char>   : [default B] the chain IDs of a point-group chain\n";
	print STDERR "    -j <char>   : [default none] the chain IDs of a lattice chain\n";
	print STDERR "    -f          : [default false] fast distance checking\n";
	exit -1;
}

my $pdbfile;
my $interact_dist = 8.0;  # min interaction distance
my $ptgp_sep = 100.0;  # min interaction distance
my $primary_chain = 'A';
my $secondary_chain = 'B';
my $lattice_chain = '';
my $fastDistCheck = 0;

## parse options (do this by hand since Getopt does not handle this well)
my $inlinefull = (join ' ',@ARGV)." ";
my @suboptions = split /(-[a-z|A-Z] )/, $inlinefull;
for ( my $i=0; $i<=$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$interact_dist = int( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-x " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$ptgp_sep = ( $suboptions[++$i] );
		$ptgp_sep =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-a " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$primary_chain = $suboptions[++$i];
		$primary_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-i " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$secondary_chain = $suboptions[++$i];
		$secondary_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-j " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$lattice_chain = $suboptions[++$i];
		$lattice_chain =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
		$pdbfile =~ s/\s*(\S+)\s*/$1/;
	} elsif ($suboptions[$i] =~ /^-f/ ) {
		$fastDistCheck = 1;
		print STDERR "Fast distance checking enabled.\n";
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

### '_' -> ' '
if ($primary_chain eq '_') {
	$primary_chain = ' ';
}
if ($secondary_chain eq '_') {
	$secondary_chain = ' ';
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

## make sure 2ary lattice chains are defined
my @sec_chain_ids = split( ':', $secondary_chain );
if ( ! defined $chains{ $sec_chain_ids[0] } ) {
	die "Chain $secondary_chain not in input!\n";
}
if ( $lattice_chain ne '' && !defined $chains{ $lattice_chain } ) {
	die "Chain $lattice_chain not in input!\n";
}

## count # of CA atoms
## TO DO: dont require this; compute optimal superposition (farm this to mammoth?)
if (scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $sec_chain_ids[0] } } ) ) {
	print STDERR "ERROR! chains '$primary_chain' and '$secondary_chain' have different residue counts! (".
	             scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $sec_chain_ids[0] } } ).")\n";
	die "Chain length mismatch!\n";
}
if ( $lattice_chain ne '' &&
     scalar( @{ $chains{ $primary_chain } } ) != scalar( @{ $chains{ $lattice_chain } } ) ) {
	print STDERR "ERROR! chains '$primary_chain' and '$secondary_chain' have different residue counts! (".
	             scalar( @{ $chains{ $primary_chain } } )." vs ".scalar( @{ $chains{ $sec_chain_ids[0] } } ).")\n";
	die "Chain length mismatch!\n";
}


## get superposition A->B
my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $sec_chain_ids[0] } );
my $delCOM = vsub ($COM_i, $COM_0);

my ($X,$Y,$Z,$W)=R2quat($R);
my $Worig = $W;
my $Wmult = 1;
if ($W < 0) { $W = -$W; $Wmult = -1; }
my $omega = acos($W);
my $symm_order = int(PI/$omega + 0.5);

# optionally ... allow input to 'force' a symmetric order
# note that this may result is a system quite far from the input system
if ($#sec_chain_ids > 0) {
	$symm_order = $sec_chain_ids[1];
}
print STDERR "Found ".$symm_order."-fold symmetric complex at chain ".$sec_chain_ids[0]."\n";
if ($symm_order != 3 && $symm_order != 4 && $symm_order != 6) {
	die "Symmetric order must equal 3,4, or 6!\n";
}

my $outer_symm_order = $symm_order;
if ($symm_order == 3) {
	$outer_symm_order = 6;
}


# now make perfectly symmetrical version of superposition
my $newW = -$Wmult *cos( PI/$symm_order );
my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
my $newQ = [$X*$S , $Y*$S, $Z*$S, $newW];
my $newR = quat2R( $newQ->[0], $newQ->[1], $newQ->[2], $newQ->[3] );

my $outerW = -$Wmult *cos( PI/$outer_symm_order );
my $S = sqrt ( (1-$outerW*$outerW)/($X*$X+$Y*$Y+$Z*$Z) );
my $outerQ = [$X*$S , $Y*$S, $Z*$S, $outerW];
my $outerR = quat2R( $outerQ->[0], $outerQ->[1], $outerQ->[2], $outerQ->[3] );


# symmetrize delCOM
my $err_pos = [0,0,0];
my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
foreach my $j (1..$symm_order) {
	$err_pos = vadd( $err_pos, mapply( $R_i,$delCOM ) );
	$R_i = mmult($newR, $R_i);
}
$delCOM = vsub( $delCOM , [ $err_pos->[0]/$symm_order , $err_pos->[1]/$symm_order , $err_pos->[2]/$symm_order ] );

# get the COM of the complex
my $sum_pos = [0,0,0];
my $curr_pos = deep_copy( $COM_0 );
$R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
foreach my $j (1..$symm_order) {
	$curr_pos = vadd( $curr_pos, mapply( $R_i,$delCOM ) );
	$sum_pos = vadd( $sum_pos, $curr_pos );
	$R_i = mmult($newR, $R_i);
}
$sum_pos = vscale( 1/$symm_order , $sum_pos );

## get translation A->G
my $delCOM_lattice;
if ( $lattice_chain ne '' ) {
	my ($R_lattice,$rmsd, $COM_i, $COM_ij) = rms_align( $chains{ $primary_chain } , $chains{ $lattice_chain } );
	$delCOM_lattice = vsub ($COM_i, $COM_0);
	if (!is_identity($R_lattice)) {
		print STDERR "     ".$R->[0][0].",".$R->[0][1].",".$R->[0][2]."\n";
		print STDERR " R = ".$R->[1][0].",".$R->[1][1].",".$R->[1][2]."\n";
		print STDERR "     ".$R->[2][0].",".$R->[2][1].",".$R->[2][2]."\n";
		die "Lattice symmetry must not have rotation!\n";  ##?? recover from this instead?
	}
	$ptgp_sep = vnorm( $delCOM_lattice );
} else {
	# no lattice chain defined; pick a random direction
	$delCOM_lattice = vsub($COM_0, $sum_pos);
	normalize( $delCOM_lattice );
	$delCOM_lattice = vscale( $ptgp_sep, $delCOM_lattice );
}

print STDERR "delCOM_lattice = ".$delCOM_lattice->[0].",".$delCOM_lattice->[1].",".$delCOM_lattice->[2]."\n";


## newR, adj_newDelCOM contain rot,trans
## compute for all subunits
my %Rs;
my %Ts;
my %ptgps;

$ptgps{"0"} = $sum_pos; # com of cplx

## 1: expand point group 0
my $R_i = deep_copy( $R_0 );
my $T_i = deep_copy( $COM_0 );
foreach my $i (1..$symm_order) {
	my $id = "0_".$i;

	$Rs{ $id } = $R_i;
	$Ts{ $id } = $T_i;

	$T_i = vadd( $T_i, mapply( $R_i,$delCOM ) );
	$R_i = mmult($newR, $R_i);
}

## 2: expand lattice groups
$R_i = deep_copy( $R_0 );
$T_i = vadd( $sum_pos , $delCOM_lattice );
foreach my $i (1..$outer_symm_order) {
	$ptgps{$i} = $T_i;

	$R_i = mmult($outerR, $R_i);
	$T_i = vadd( $sum_pos, mapply( $R_i,$delCOM_lattice ) );
}

## 3: expand point gps of lattice groups
foreach my $j (1..$outer_symm_order) {
	foreach my $i (1..$symm_order) {
		my $id  = $j."_".$i;
		$Rs{ $id } = $Rs{ "0_".$i };
		$Ts{ $id } = vadd( $ptgps{$j} , vsub( $Ts{ "0_".$i } , $ptgps{0}) );
	}
}


## dist checks
## first-pass filter throws out monomers very far from the primary
my $counter = 0;
my %excludeinterface = ();
foreach my $j (0..$outer_symm_order) {
	foreach my $i (1..$symm_order) {
		my $id  = $j."_".$i;
		my $delXY = vsub( $COM_0,$Ts{$id} );
		my $dist2XY = vnorm2( $delXY );

		if ( sqrt($dist2XY) > 2*$monomerRadius + $interact_dist ) {
			print STDERR " [$counter] Excluding interface '".$id."'\n";
			$excludeinterface{ $id } = $counter++;
		}
	}
}

my %symminterface = ();
$counter = 0;
if ($fastDistCheck == 1) {
	foreach my $j (0..$outer_symm_order) {
		foreach my $i (1..$symm_order) {
			my $id  = $j."_".$i;
			next if (defined $excludeinterface{ $id });

			# we have a hit! tag NCS copy $id as a symmetic interface
			print STDERR " Adding interface '".$id."'\n";
			$symminterface{ $id } = $counter++;
		}
	}
} else {
	#print "COM = ".$COM_0->[0].",".$COM_0->[1].",".$COM_0->[2]."\n";
	foreach my $X_i ( @{ $chains{ $primary_chain } } ) {
		foreach my $Y_i (  @{ $chains{ $primary_chain } } ) {
			foreach my $j (0..$outer_symm_order) {
				foreach my $i (1..$symm_order) {
					my $id  = $j."_".$i;
					next if (defined $symminterface{ $id });
					next if (defined $excludeinterface{ $id });
	
					#   x_i = R_i * (x_0 - COM_0) + COM_i
					#   The rms function already ofsets x_0 by -COM_0
					my $rX_i = vadd( $X_i , $COM_0 );
					my $rY_j = vadd( mapply($Rs{$id}, $Y_i) , $Ts{$id} );
					my $delXY = vsub( $rY_j,$rX_i );
					my $dist2XY = vnorm2( $delXY );

					if ($dist2XY < $interact_dist*$interact_dist) {
						# we have a hit! tag NCS copy $id as a symmetic interface
						print STDERR " Adding interface '".$id."'\n";
						$symminterface{ $id } = $counter++;
					}
				}
			}
		}
	}
}

##
# find the equation for the energy of the complex
# first find the 1->n interfaces that come in pairs
my @syminterfaces = sort { $symminterface{$a} <=> $symminterface{$b} } keys %symminterface;

# the energy equation ...
my %energy_counter;
foreach my $interface (@syminterfaces) {
	$energy_counter{ $interface } = 1;
}
# delete self-interface(???)
delete ($energy_counter{ $syminterfaces[0] });

####
####  This appears to be calculated incorrectly for certain orientations
####  For now compute every interface energy in Rosetta
####  This will lead to a slowdown in Rosetta but it is not clear how significant
####
# OUTER1: foreach my $i (1..$#syminterfaces) {
# 	next if (!defined  $energy_counter{ $syminterfaces[$i] } );
# 
# 	my $R_i = $Rs{$syminterfaces[$i]};
# 	my $T_i = vsub( $Ts{$syminterfaces[$i]} , $COM_0 );
# 
# 	foreach my $j ($i+1..$#syminterfaces) {
# 		my $R_j = $Rs{$syminterfaces[$j]};
# 		my $T_j = vsub( $Ts{$syminterfaces[$j]} , $COM_0 );
# 
# 		# if transform i is the inverse of transform j we have our (i,j) pair
# 		if ( is_inverse( $R_i,$T_i, $R_j,$T_j ) ) {
# 			$energy_counter{ $syminterfaces[$i] } = 2;
# 			delete ($energy_counter{ $syminterfaces[$j] });
# 			next OUTER1;
# 		}
# 	}
# }
#######################################
#######################################
#######################################


## symm file gen
# symmetry_name c4
# E = 2*VRT2
# anchor_residue 17
my $symmname = $pdbfile;
$symmname =~ s/\.pdb$//;
$symmname = $symmname."_P".$symm_order;
print "symmetry_name $symmname\n";
print "E = 1*VRT".$syminterfaces[0]."_base";
foreach my $complex (sort { $symminterface{$a} <=> $symminterface{$b} } keys %energy_counter) {
	print " + ".$energy_counter{$complex}."*(VRT".$syminterfaces[0]."_base".":VRT".$complex."_base".")";
}
print "\n";
print "anchor_residue $minRes\n";

#######################################
######
######   XYZ
######
#######################################
print "virtual_coordinates_start\n";
print "xyz VRT0  ".
			sprintf("%.4f,%.4f,%.4f", 1, 0, 0)."  ".
			sprintf("%.4f,%.4f,%.4f", 0, 1, 0)."  ".
			sprintf("%.4f,%.4f,%.4f", $sum_pos->[0]+1, $sum_pos->[1], $sum_pos->[2])."\n";

#subunits of first pt_gp
foreach my $i (1..$symm_order) {
	my $id = "0_".$i;
	my $id_sibling = "0_".($i+1);
	if ($i == $symm_order) { $id_sibling = "0_1"; }
	my $parent_com = $sum_pos;

	# x points from origin to parent CoM
	my $myX = vsub( $parent_com , $Ts{ $id } );
	my $myY = vsub( $Ts{ $id_sibling } , $Ts{ $id } );
	#my $myZ = mapply( $Rs{ $id }, [0,0,1]);
	#my $myY = cross( $myZ, $myX );
	# y in pt-gp plane

	normalize( $myX );
	$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
	normalize( $myY );

	print "xyz VRT$id  ".
				sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	if ($i == 1) {
		print "xyz VRT0"."_ctrl  ".
					sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
	}
	print "xyz VRT$id"."_base  ".
				sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $Ts{ $id }->[0], $Ts{ $id }->[1], $Ts{ $id }->[2])."\n";
}


foreach my $i (1..$outer_symm_order) {
	my $id_sibling = ($i+1);
	if ($i == $symm_order) { $id_sibling = 1; }

	# pointing to each point group
	my $parent_com = $sum_pos;
	my $myX = vsub( $ptgps{0} , $ptgps{$i} );
	my $myY = vsub( $ptgps{$id_sibling} , $ptgps{$i} );
	normalize( $myX );
	$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
	normalize( $myY );
	print "xyz VRT$i"."_dir  ".
				sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";

	# redirect each point group
	if ($outer_symm_order > $symm_order) {
		my $id = "0_1";
		$id_sibling = "0_2";
		my $myX = vsub( $parent_com , $Ts{ $id } );
		my $myY = vsub( $Ts{ $id_sibling } , $Ts{ $id } );
		normalize( $myX );
		$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
		normalize( $myY );
	}
	print "xyz VRT$i"."_redir  ".
				sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
}

#coms of each pt_gp
foreach my $j (1..$outer_symm_order) {
	#my $T_global = vsub( $Ts{"0_".$j}, $sum_pos);
	#normalize( $T_global );
	#$T_global = vscale( $ptgp_sep, $T_global );
	my $parent_com = $ptgps{$j}; #vadd( $sum_pos, $T_global );

	foreach my $i (1..$symm_order) {
		my $id = $j."_".$i;
		my $id_sibling = $j."_".($i+1);
		if ($i == $symm_order) { $id_sibling = $j."_1"; }

		# x points from origin to parent CoM
		my $myX = vsub( $parent_com , $Ts{ $id } );
		#my $myZ = mapply( $Rs{ $id }, [0,0,1]);
		#my $myY = cross( $myZ, $myX );
		# y in pt-gp plane
		normalize( $myX );
		my $myY = vsub( $Ts{ $id_sibling } , $Ts{ $id } );
	
		$myY = vsub( $myY , vscale(dot($myY, $myX), $myX) );
		normalize( $myY );

		print "xyz VRT$id  ".
					sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
		if ( $i == 1 ) {
			print "xyz VRT$j"."_ctrl  ".
						sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
						sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
						sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2])."\n";
		}
		print "xyz VRT$id"."_base  ".
					sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $Ts{ $id }->[0], $Ts{ $id }->[1], $Ts{ $id }->[2])."\n";
	}
}
print "virtual_coordinates_stop\n";


#######################################
######
######   connect
######
#######################################
print "connect_virtual JUMP0 VRT0 VRT0_ctrl\n";  # root jump
print "connect_virtual JUMP0_1 VRT0_ctrl VRT0_1\n";  # root jump
foreach my $j (1..$symm_order) {
	my $id = "0_".$j;
	if ($j != 1) {
		print "connect_virtual JUMP".$id." VRT0_1 VRT$id"."\n";
	}
	print "connect_virtual JUMP".$id."_to_com VRT$id VRT$id"."_base"."\n";
	if (exists $symminterface{$id}) {
		print "connect_virtual JUMP".$id."_to_subunit VRT$id"."_base SUBUNIT"."\n";
	}
}

foreach my $j (1..$outer_symm_order) {
	my $id1 = $j."_1";
	#print "connect_virtual JUMP".$j." VRT$id VRT$id1"."\n";
	print "connect_virtual JUMP$j"."_to_dir  VRT0_ctrl VRT$j"."_dir"."\n";
	print "connect_virtual JUMP$j"."_to_redir VRT$j"."_dir VRT$j"."_redir"."\n";
	print "connect_virtual JUMP$j"."_to_ctrl VRT$j"."_redir VRT$j"."_ctrl"."\n";
	print "connect_virtual JUMP$id1 VRT$j"."_ctrl VRT$id1\n";  # root jump
	foreach my $i (1..$symm_order) {
		my $id = $j."_".$i;
		if ($i != 1) {
			print "connect_virtual JUMP".$id." VRT".$j."_1 VRT$id"."\n";
		}
		print "connect_virtual JUMP".$id."_to_com VRT$id VRT$id"."_base"."\n";
		if (exists $symminterface{$id}) {
			print "connect_virtual JUMP".$id."_to_subunit VRT$id"."_base SUBUNIT"."\n";
		}
	}
}

## dofs
print "set_dof JUMP1_to_ctrl x($ptgp_sep)\n";
print "set_dof JUMP0_1 angle_z\n";
print "set_dof JUMP0_1_to_com x\n";
print "set_dof JUMP0_1_to_subunit angle_x angle_y angle_z\n";


## jumpgroup between pointgroups
print "set_jump_group JUMPGROUP1 ";
foreach my $j (1..$outer_symm_order) {
	print "JUMP$j"."_to_ctrl ";
}
print "\n";

print "set_jump_group JUMPGROUP2 ";
foreach my $j (0..$outer_symm_order) {
	my $id1 = $j."_1";
	print "JUMP$id1 ";
}
print "\n";

## jumpgroups to com, subunit
print "set_jump_group JUMPGROUP3 ";
foreach my $j (0..$outer_symm_order) {
	foreach my $i (1..$symm_order) {
		my $id = $j."_".$i;
		print "JUMP$id"."_to_com ";
	}
}
print "\n";
print "set_jump_group JUMPGROUP4 ";
foreach my $j (0..$outer_symm_order) {
	foreach my $i (1..$symm_order) {
		my $id = $j."_".$i;
		if (exists $symminterface{$id}) { print "JUMP$id"."_to_subunit " };
	}
}
print "\n";

########################################
## write output pdb
########################################
my $outpdb = $pdbfile;
my $outmon = $pdbfile;
my $outmdl = $pdbfile;
my $outkin = $pdbfile;

my $suffix = "_model_$primary_chain"."$secondary_chain";
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
foreach my $j (0..$outer_symm_order) {
	foreach my $i (1..$symm_order) {
		my $id  = $j."_".$i;

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

			if (defined $symminterface{ $id }) {
				print OUTMDL $linecopy."\n";
			}
		}

		print OUTPDB "TER   \n";
		if (defined $symminterface{ $id }) {
			print OUTMDL "TER   \n";
		}
		print STDERR "Writing interface ".$id." as chain ".substr ($chains, $chnidx, 1)."\n";
		$chnidx++;
	}
}

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
