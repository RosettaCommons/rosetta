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

## parameters
my $MAX_TRANS_ERR = 5.0; # angstrom
my $MAX_ROT_ERR   = 3.0; # degrees

###############################################################################

if ($#ARGV < 0) {
	print STDERR "usage: $0 [options]\n";
	print STDERR "example:   $0 -a A -r 12.0 -p mystructure.pdb\n";
	print STDERR "NOTE: This is experimental code that tries to automatically determine the symmetry of the\n";
	print STDERR "input system.  Use with caution!\n";
	print STDERR "options: \n";
	print STDERR "    -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
	print STDERR "    -r <real>   : [default 12.0] the max CA-CA distance between two interacting chains\n";
	print STDERR "    -a <char>   : [default A] the chain ID of the main chain\n";
	print STDERR "    -f          : [default false] fast distance checking\n";
	exit -1;
}

my $pdbfile;
my $interact_dist = 12.0;  # min interaction distance
my $primary_chain = 'A';
my @secondary_chains = ();
my $fastDistCheck = 0;

## parse options (do this by hand since Getopt does not handle this well)
my $inlinefull = (join ' ',@ARGV)." ";
my @suboptions = split /(-[a-z|A-Z] )/, $inlinefull;
for ( my $i=0; $i<=$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$interact_dist = int( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-a " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$primary_chain = $suboptions[++$i];
		$primary_chain =~ s/\s*(\S+)\s*/$1/;
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
				push @secondary_chains, $chnid if ($chnid ne $primary_chain);
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

my $NCS_ops = {};
my $R_0 = [ [1,0,0], [0,1,0], [0,0,1] ];  # R_0
my $COM_0 = recenter( $chains{ $primary_chain } );
$NCS_ops->{R} = $R_0;
$NCS_ops->{T} = $COM_0;
$NCS_ops->{PATH} = "";
$NCS_ops->{CHILDREN} = [];


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

# pass 1 -- throw out nonsymm stuff
my @allQs;
my @allCOMs;
my @allCorrCOMs;
my @sym_orders;
my @secondary_chains_filt;

# remember 'best' so far
my $best_score=999;
my $best_order=1;
my @best_chain;
foreach my $sec_chain (@secondary_chains) {
	my $chain1 = deep_copy( $chains{ $primary_chain });
	my $chain2 = deep_copy( $chains{ $sec_chain });
	my ($R,$rmsd, $COM_i, $COM_ij) = rms_align( $chain1,$chain2 );
	my $del_COM = vsub ($COM_i, $COM_0);

	my ($X,$Y,$Z,$W)=R2quat($R);
	my $Worig = $W;
	my $Wmult = 1;
	if ($W < 0) { $W = -$W; $Wmult = -1; }
	my $omega = acos($W);
	my $sym_order = int(PI/$omega + 0.5);

	# now make perfectly symmetrical version of superposition
	my $newW = -$Wmult *cos( PI/$sym_order );
	my $S = sqrt ( (1-$newW*$newW)/($X*$X+$Y*$Y+$Z*$Z) );
	my $newQ = [$X*$S , $Y*$S, $Z*$S, $newW];
	my $newR = quat2R( $newQ->[0], $newQ->[1], $newQ->[2], $newQ->[3] );

	# error (rot)
	my $Werror = 180/PI * abs( $Worig + $newW );

	# error (trans)
	my $err_pos = [0,0,0];
	my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (1..$sym_order) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$del_COM ) );
		$R_i = mmult($newR, $R_i);
	}

	next if (vnorm( $err_pos ) > $MAX_TRANS_ERR);
	next if ($Werror > $MAX_ROT_ERR);

	#print STDERR "$pdbfile: Found ".$sym_order."-fold (err=".$Werror."/".vnorm( $err_pos ).") symmetric complex at chain ".$sec_chain."\n";
	my $adj_newDelCOM = vsub( $del_COM , [ $err_pos->[0]/$sym_order , $err_pos->[1]/$sym_order , $err_pos->[2]/$sym_order ] );
	push @secondary_chains_filt, $sec_chain;
	push @sym_orders, $sym_order;
	push @allQs, $newQ;
	push @allCOMs, $del_COM;
	push @allCorrCOMs, $adj_newDelCOM;

	my $score = $Werror+vnorm( $err_pos );
	if ($sym_order > $best_order || ( $sym_order == $best_order && $best_score > $score) ) {
		$best_order = $sym_order;
		$best_score = $score;
		@best_chain = ($sec_chain);
	}
}

# pass 2 -- try to see if any pairs are Dn-compatible
my @secondary_chains_filt_dn;
my @scores_dn;
my @symm_orders_dn;
foreach my $i (0..$#secondary_chains_filt) {
foreach my $j (0..$#secondary_chains_filt) {
	next if ($i eq $j);
	next if ($sym_orders[$i] != 2);

	# symm agreement
	# the two transformations must be about perpendicular axes
	my $X = [ $allQs[$i]->[0],  $allQs[$i]->[1],  $allQs[$i]->[2] ];
	my $Y = [ $allQs[$j]->[0],  $allQs[$j]->[1],  $allQs[$j]->[2] ];
	normalize($X); normalize($Y);
	my $ang1 = abs( 90-angle_deg($X,$Y) ); # angle between sym axes
	next if ($ang1 > $MAX_ROT_ERR);

	# trans agreement
	my $NCS_opstemp = {};
	my $newR = quat2R( $allQs[$i]->[0], $allQs[$i]->[1], $allQs[$i]->[2], $allQs[$i]->[3] );
	my $adj_newDelCOM = $allCorrCOMs[$i];
	expand_symmops_by_split( $NCS_opstemp, $newR, $adj_newDelCOM, $sym_orders[$i]);

	my $newQ = $allQs[ $j ];
	my $del_COM = $allCOMs[ $j ];
	my $sym_order = $sym_orders[ $j ];
	$newR = quat2R( $newQ->[0], $newQ->[1], $newQ->[2], $newQ->[3] );

	my $com_complex = $NCS_opstemp->{T};
	my $com_secondary = vadd( $del_COM, vadd( $COM_0, mapply( $newR , vsub ( $com_complex , $COM_0 ) ) ) );
	my $newDelCOM = vsub ( $com_secondary , $com_complex );
	my $err_pos = [0,0,0];
	my $R_i = [ [1,0,0], [0,1,0], [0,0,1] ];
	foreach my $j (1..$sym_order) {
		$err_pos = vadd( $err_pos, mapply( $R_i,$newDelCOM ) );
		$R_i = mmult($newR, $R_i);
	}
	my $axis_proj_i =  [ $allQs[1-$i]->[0],  $allQs[1-$i]->[1],  $allQs[1-$i]->[2] ];
	normalize( $axis_proj_i );
	my $del_COM_inplane    = vsub( $newDelCOM , vscale(dot($newDelCOM,$axis_proj_i),$axis_proj_i) );
	$err_pos = vscale( $sym_orders[ $i ] , $del_COM_inplane );
	next if (vnorm( $err_pos ) > $MAX_TRANS_ERR);

	#print STDERR "$pdbfile: Found D".$sym_orders[$j]." (err=".$ang1."/".vnorm( $err_pos ).") symmetric complex at chains ".
	#   $secondary_chains_filt[$i]."/".$secondary_chains_filt[$j]."\n";

	my $score = ($ang1+vnorm( $err_pos ));
	if (2*$sym_orders[$j] > $best_order || ( $sym_order == $best_order && $best_score > $score) ) {
		$best_order = 2*$sym_orders[$j];
		$best_score = $score;
		@best_chain = ($secondary_chains_filt[$i],$secondary_chains_filt[$j]);
	}
}
}

my $symmtype= "C".$best_order;
if ( scalar(@best_chain) == 2 ) {
	$symmtype = "D".($best_order/2);
}

if ($symmtype eq "C1") {
	print STDERR "$pdbfile: No symmetry found.\n";	
	exit;
}
print STDERR "$pdbfile: Found ".$symmtype." symmetric complex!\n";

my $outfile = $pdbfile;
$outfile =~ s/\.pdb/_$symmtype.symm/;

# make symmdef file
@allQs = ();
@allCOMs = ();
@sym_orders = ();

foreach my $sec_chain (@best_chain) {
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
	my $del_COM = vsub ($COM_i, $COM_0);
	push @allCOMs, $del_COM;

	my ($X,$Y,$Z,$W)=R2quat($R);
	my $Worig = $W;
	my $Wmult = 1;
	if ($W < 0) { $W = -$W; $Wmult = -1; }
	my $omega = acos($W);
	my $sym_order = int(PI/$omega + 0.5);

	# optionally ... allow input to 'force' a symmetric order
	# NOTE THAT THIS MAY RESULT IS A SYSTEM QUITE FAR FROM THE INPUT SYSTEM
	if ($#sec_chain_ids > 0) {
		$sym_order = $sec_chain_ids[1];
	}

	push @sym_orders, $sym_order;
	#print STDERR "Found ".$sym_order."-fold (".(PI/$omega).") symmetric complex at chain ".$sec_chain_ids[0]."\n";

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
	my $X = [ $allQs[0]->[0],  $allQs[0]->[1],  $allQs[0]->[2] ];
	my $Y = [ $allQs[1]->[0],  $allQs[1]->[1],  $allQs[1]->[2] ];

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

# 2) Icosahedral symmetries
#   ??? TO DO
#  * 5-fold symm axis must lie in y-z plane
#  * two fold axes on y and z
#  * 3-fold about line x=y=z (?)


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
		my $axis_proj_i =  [ $allQs[1-$i]->[0],  $allQs[1-$i]->[1],  $allQs[1-$i]->[2] ];
		# project $delCOM along this axis
		normalize( $axis_proj_i );
		my $del_COM_inplane    = vsub( $newDelCOM , vscale(dot($newDelCOM,$axis_proj_i),$axis_proj_i) );
		$err_pos = vscale( $sym_orders[ $i ] , $del_COM_inplane );
	}
	#print STDERR "  translation error = ".vnorm( $err_pos )."\n";


	# special case for icosehedral symmetry
	# see above for restrictions

	# get the center of the symmgp
	my $adj_newDelCOM = vsub( $newDelCOM , [ $err_pos->[0]/$sym_order , $err_pos->[1]/$sym_order , $err_pos->[2]/$sym_order ] );

	expand_symmops_by_split( $NCS_ops, $newR, $adj_newDelCOM, $sym_order);
}

my ($nnodes,$nleaves) = tree_size( $NCS_ops );
#print STDERR "Found a total of $nleaves monomers in the symmetric complex.\n";
#print STDERR "Placing $nnodes virtual residues.\n";

#exit -1;

## dist checks
my $symops = tree_traverse( $NCS_ops );
my %symminterface = ();
my $counter = 0;

## first-pass filter throws out monomers very far from the primary
my %excludeinterface = ();
foreach my $symop (@{ $symops }) {
	my $id = $symop->{PATH};
  	my $delXY = vsub( $COM_0,$symop->{T} );
  	my $dist2XY = vnorm2( $delXY );

	if ( sqrt($dist2XY) > 2*$monomerRadius + $interact_dist ) {
		#print STDERR " [$counter] Excluding interface '".$symop->{PATH}."'\n";
		$excludeinterface{ $id } = $counter++;
	}
}

$counter = 0;

if ($fastDistCheck == 1) {
	foreach my $symop (@{ $symops }) {
		my $id = $symop->{PATH};
		next if (defined $excludeinterface{ $id });

		# we have a hit! tag NCS copy $j_symm as a non-symmetic interface
		#print STDERR " Adding interface '".$symop->{PATH}."'\n";
		$symminterface{ $id } = $counter++;
	}
} else {
	#print "COM = ".$COM_0->[0].",".$COM_0->[1].",".$COM_0->[2]."\n";
	foreach my $X_i ( @{ $chains{ $primary_chain } } ) {
		foreach my $Y_i (  @{ $chains{ $primary_chain } } ) {
			foreach my $symop (@{ $symops }) {
				my $id = $symop->{PATH};
				next if (defined $symminterface{ $id });
				next if (defined $excludeinterface{ $id });
	
				#   x_i = R_i * (x_0 - COM_0) + COM_i
				#   The rms function already ofsets x_0 by -COM_0
				my $rX_i = vadd( $X_i , $COM_0 );
				my $rY_j = vadd( mapply($symop->{R}, $Y_i) , $symop->{T} );
				my $delXY = vsub( $rY_j,$rX_i );
				my $dist2XY = vnorm2( $delXY );
	 
				if ($dist2XY < $interact_dist*$interact_dist) {
					# we have a hit! tag NCS copy $j_symm as a non-symmetic interface
					#print STDERR " Adding interface '".$symop->{PATH}."'\n";
					$symminterface{ $id } = $counter++;
				}
			}
		}
	}
}

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
#foreach my $inter_n (keys %energy_counter) {
#	print " $inter_n -> $energy_counter{$inter_n}\n";
#}


# now, for each unique 1->n interface
# sum up the number of equivalent i->j interfaces
#foreach my $symop_i (@{ $symops }) {
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



#######################################
## symm file gen
## write output symm file
# symmetry_name c4
# E = 2*VRT2
# anchor_residue 17
open (OUTSYM, ">$outfile");

my $symmname = $pdbfile;
$symmname =~ s/\.pdb$//;
$symmname = $symmname."_".get_topology( $NCS_ops );
print OUTSYM "symmetry_name $symmname\n";
print OUTSYM "E = ".($nleaves)."*VRT".$syminterfaces[0]."_base";
foreach my $complex (sort { $symminterface{$a} <=> $symminterface{$b} } keys %energy_counter) {
	print OUTSYM " + ".$energy_counter{$complex}."*(VRT".$syminterfaces[0]."_base".":VRT".$complex."_base".")";
}
print OUTSYM "\n";
print OUTSYM "anchor_residue com\n";

# virtual_coordinates_start
# xyz VRT1 -1,0,0 0,1,0 0,0,0
# xyz VRT2 0,-1,0 -1,0,0 0,0,0
# xyz VRT3 1,0,0 0,-1,0 0,0,0
# xyz VRT4 0,1,0 1,0,0 0,0,0
# virtual_coordinates_stop
my $vrts_by_depth = tree_traverse_by_depth( $NCS_ops );
my ($vrt_lines, $connect_lines, $dof_lines, $debug_lines) = fold_tree_from_ncs( $NCS_ops , $vrts_by_depth, \%symminterface );
print OUTSYM "virtual_coordinates_start\n";
foreach my $vrt_line (@{ $vrt_lines }) {
	print OUTSYM $vrt_line."\n";
}
# connect_virtual JUMP1 VRT1 VRT2
# connect_virtual JUMP2 VRT1 VRT3
# connect_virtual JUMP3 VRT1 VRT4
print OUTSYM "virtual_coordinates_stop\n";
foreach my $connect_line (@{ $connect_lines }) {
	print OUTSYM $connect_line."\n";
}
foreach my $dof_line (@{ $dof_lines }) {
	print OUTSYM $dof_line."\n";
}
my $counter = 0;
foreach my $vrt_level ( @{ $vrts_by_depth } ) {
	if ($counter > 1) {
		print OUTSYM "set_jump_group JUMPGROUP$counter";
		foreach my $tag (@{ $vrt_level }) {
			# jump to _first_ child in each group
			if ($tag =~ /_0$/) {
				print OUTSYM " JUMP$tag";
			}
		}
		print OUTSYM "\n";
	}
	$counter++;
}
print OUTSYM "set_jump_group JUMPGROUP$counter";
foreach my $leaf ( @{ $symops } ) {
	my $tag = $leaf->{PATH};
	print OUTSYM " JUMP$tag"."_to_com";
}
print OUTSYM "\n";
$counter++;
print OUTSYM "set_jump_group JUMPGROUP$counter";
foreach my $tag ( keys %symminterface ) {
	print OUTSYM " JUMP$tag"."_to_subunit";
}
print OUTSYM "\n";



#############################################################################
#############################################################################

###
###   SUBROUTINES
###

#############################################################################
#############################################################################

sub angle_deg {
    my ($x, $y) = @_;

    my $x_m = vnorm($x);
    my $y_m = vnorm($y);

	my $xy = ($x->[0]*$y->[0]+$x->[1]*$y->[1]+$x->[2]*$y->[2]);

    return (180*acos ( $xy / ($x_m * $y_m) )/PI);
}

#############
# expand_symmops_by_split( $NCS_ops, $newR, $adj_newDelCOM, $sym_order)
sub expand_symmops_by_split {
	my ( $tree, $newR, $newDelT, $sym_order ) = @_;

	my $newNCSops = {};
	my $COM_0 = [0 , 0 , 0];
	$newNCSops->{R} = [ [1,0,0], [0,1,0], [0,0,1] ];
	$newNCSops->{CHILDREN} = [];

	my $COM_i = [ 0,0,0 ];
	my $R_i   = [ [1,0,0], [0,1,0], [0,0,1] ];
	my $newCOM0 = [0,0,0];

	foreach my $i (0..$sym_order-1) {
		my $newNCSops_i = deep_copy( $NCS_ops );

		# rotate about the center of mass of the subtree
		#    then translate to the new CoM
		apply_transformation ( $newNCSops_i, $R_i, $NCS_ops->{T}, $COM_i, $i );
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
##    get the size of a tree
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
##   traverse the tree
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
##     traverse every node in the tree returning them sorted by depth
sub tree_traverse_by_depth {
	my $tree = shift;
	my $depth = get_depth( $tree );
	my $nodes_by_depth = [  ];
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
##    the meat of the code lives here
sub fold_tree_from_ncs {
	my $tree = shift;
	my $nodes_by_depth = shift;
	my $connected_subunits = shift;

	# root doesnt have parents or siblings
	my $parent_com = vadd( [1,0,0] ,get_com( $tree ) );
	my $sibling_com = vadd( [0,1,0] , get_com( $tree ) );
	my $nsiblings = 1;
	my $vrt_lines = [];
	my $connect_lines = [];
	my $dof_lines = [];
	my $debug_lines = [];

	fold_tree_from_ncs_recursive( $tree , $parent_com, $nsiblings, $sibling_com, $vrt_lines, 
	                              $connect_lines, $dof_lines, $debug_lines,$nodes_by_depth, $connected_subunits );
	return ( $vrt_lines, $connect_lines, $dof_lines , $debug_lines);
}

sub fold_tree_from_ncs_recursive {
	my ($tree,$parent_com,$nsiblings, $sibling_com,
	    $vrt_lines,$connect_lines,$dof_lines,$debug_lines,$nodes_by_depth,$connected_subunits) = @_;

	my $origin = get_com( $tree );
	my $id_string = $tree->{PATH};

	# x points from origin to parent CoM
	my $myX = vsub( $parent_com , $origin );
	my $myZ = mapply( $tree->{R}, [0,0,1]);

	# y is whatever is left
	my $myY = cross( $myZ, $myX );

	normalize( $myX );
	normalize( $myY );

	# 
	my $nchildren = scalar( @{ $tree->{CHILDREN} } );

	## recursive call
	foreach my $child_idx ( 0..$nchildren-1 ) {
		my $child_idstring = $tree->{CHILDREN}->[$child_idx]->{PATH};

		# now recursively call each child
		my $sibling_com = get_com( $child_idx == 0? $tree->{CHILDREN}->[1] : $tree->{CHILDREN}->[0] );
		fold_tree_from_ncs_recursive( $tree->{CHILDREN}->[$child_idx], $origin, $nchildren , $sibling_com, 
									  $vrt_lines, $connect_lines, $dof_lines, $debug_lines, $nodes_by_depth,$connected_subunits);
	}

	# xyz VRT1 -1,0,0 0,1,0 0,0,0
	push @{ $vrt_lines }, "xyz VRT$id_string  ".
				sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
				sprintf("%.4f,%.4f,%.4f", $parent_com->[0], $parent_com->[1], $parent_com->[2]);
	#my $fakePDBline1 = sprintf "ATOM    %3d  C   ORI Z   1     %7.3f %7.3f %7.3f  1.00  0.00\n", 1,
	#					 $parent_com->[0], $parent_com->[1], $parent_com->[2];
	#my $fakePDBline2 = sprintf "ATOM    %3d  O   X   Z   2     %7.3f %7.3f %7.3f  1.00  0.00\n", 2,
	#					 $parent_com->[0]+$myX->[0], $parent_com->[1]+$myX->[1], $parent_com->[2]+$myX->[2];
	#my $fakePDBline3 = sprintf "ATOM    %3d  C   Y   Z   3     %7.3f %7.3f %7.3f  1.00  0.00\n", 3,
	#					 $parent_com->[0]+$myY->[0], $parent_com->[1]+$myY->[1], $parent_com->[2]+$myY->[2];
	#push @{ $debug_lines }, $fakePDBline1;
	#push @{ $debug_lines }, $fakePDBline2;
	#push @{ $debug_lines }, $fakePDBline3;
	my $kineline1 = "{$id_string up} P  ".$parent_com->[0]." , ".$parent_com->[1]." , ".$parent_com->[2]."\n";
	my $kineline2 = "{$id_string down}  ".$origin->[0]." , ".$origin->[1]." , ".$origin->[2]."\n";
	push @{ $debug_lines }, $kineline1;
	push @{ $debug_lines }, $kineline2;

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

		if ($is_controlling >= 1) {  # level 1 jumps don't move, hence '>'
			my $x_dist = vnorm( vsub( $parent_com,$origin) ) ;

			# in icosohedral (and maybe octohedral/tetraherdal) symm we have 
			#     2-fold symm operators where x==0 and must ==0 to keep other symmetries intact
			# for these cases we cannot allow x/angle_x to move
			if ( $x_dist > 1e-3 ) {
				if ($nsiblings == 2) {
					push @{ $dof_lines },     "set_dof JUMP$child0_idstring x($x_dist) angle_x";
					#push @{ $dof_lines },     "set_dof JUMPGROUP$is_controlling x angle_x";
				} elsif ($nsiblings > 2) {
					push @{ $dof_lines },     "set_dof JUMP$child0_idstring x($x_dist)";
					#push @{ $dof_lines },     "set_dof JUMPGROUP$is_controlling x";
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
			#     2-fold symm operators where x==0 and must ==0 to keep other symmetries intact
			# for these cases we cannot allow x/angle_x to move
			if ( $x_dist > 1e-3 ) {
				if ($nsiblings == 2) {
					push @{ $dof_lines },     "set_dof JUMP$id_string"."_to_com x($x_dist) angle_x";
				} elsif ($nsiblings > 2) {
					push @{ $dof_lines },     "set_dof JUMP$id_string"."_to_com x($x_dist)";
				}
			}
		}

		# 2 -- if an interface res, jump to the subunit
		if (defined $connected_subunits->{ $tree->{PATH} }) {
			# if this is "subunit 0"
			if ($connected_subunits->{ $tree->{PATH} } == 0) {
				push @{ $dof_lines },     "set_dof JUMP".$id_string."_to_subunit angle_x angle_y angle_z";
			}

			push @{ $connect_lines }, "connect_virtual JUMP".$id_string."_to_subunit VRT$id_string"."_base SUBUNIT";

		}


		push @{ $vrt_lines }, "xyz VRT$id_string"."_base  ".
					sprintf("%.4f,%.4f,%.4f", $myX->[0], $myX->[1], $myX->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $myY->[0], $myY->[1], $myY->[2])."  ".
					sprintf("%.4f,%.4f,%.4f", $origin->[0], $origin->[1], $origin->[2]);


		#my $fakePDBline1 = sprintf "ATOM    %3d  C   ORI Z   1     %7.3f %7.3f %7.3f  1.00  0.00\n", 1,
		#                     $origin->[0], $origin->[1], $origin->[2];
		#my $fakePDBline2 = sprintf "ATOM    %3d  O   X   Z   2     %7.3f %7.3f %7.3f  1.00  0.00\n", 2,
		#                     $origin->[0]+$myX->[0], $origin->[1]+$myX->[1], $origin->[2]+$myX->[2];
		#my $fakePDBline3 = sprintf "ATOM    %3d  C   Y   Z   3     %7.3f %7.3f %7.3f  1.00  0.00\n", 3,
		#                     $origin->[0]+$myY->[0], $origin->[1]+$myY->[1], $origin->[2]+$myY->[2];
		#push @{ $debug_lines }, $fakePDBline1;
		#push @{ $debug_lines }, $fakePDBline2;
		#push @{ $debug_lines }, $fakePDBline3;
	}
}

##########
## my $x = get_com( $subtree )
##    get the center of mass of a subtree
##    do this by summing the CoM's of all subunits
sub get_com {
	my $tree = shift;
	return $tree->{T};
}

##########
##  Tree iterator
##    depth-first traversal
# sub tree_iterator {
# 	my $tree = shift;
# 
#     return sub {
#         # code to calculate $next_state or $done;
#         return undef if $done;
#         return $current_state = $next_state;   
#     };
# }
