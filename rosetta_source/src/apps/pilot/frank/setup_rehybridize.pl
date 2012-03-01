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
my $DENSGENAPP = "/work/dimaio/rosetta/rosetta_source/bin/pdb_to_map.default.linuxgccrelease".
                 " -database /work/dimaio/rosetta/rosetta_database ".
                 " -edensity:mapreso 10.0 -edensity::grid_spacing 2.0 -in:file:fullatom";
my $THESEUSAPP = "/work/dimaio/bin/theseus  -f ";
my $THESEUSCUT = 4;

## paths
my $NATIVEDIR        = "native/";
my $FRAG3FILE      = "fragments/%s_templatesvall.25.3mers";
my $FRAG9FILE      = "fragments/%s_templatesvall.25.9mers";

my $CSTFILENAMES  = "cluster%d.filt.dist_csts";
my $FLAGFILENAMES = "flags%d";
my $MAPFILENAMES  = "templates%d.mrc";
my $ALNTEMPLATEDIR  = "aligned_templates";
my $UNALNTEMPLATEDIR  = "unaligned_templates";

my $XMLFILENAMES = "hybridize%d.xml";
my $RUNFILENAMES = "run%d.sh";
my $HYBRIDIZEOPTIONS = "batch=2 stage1_increase_cycles=1.0 stage2_increase_cycles=1.0 linmin_only=0";

#####
#####

## structure alignments
#my @RMS_CUTOFFS = (10,5,4,3,2,1.5,1);
#my $CLUSTERCUTOFF = 0.40;
#my $ALIGNCUTOFF   = 0.20;  # to get better superpositions, trade coverage for alignment
my @RMS_CUTOFFS = (10,5,4,3,2,1.5,1);
my $CLUSTERCUTOFF = 0.40;
my $ALIGNCUTOFF   = 0.20;  # to get better superpositions, trade coverage for alignment

## amino acid map
my %one_to_three = ( 
	'G' => 'GLY', 'A' => 'ALA', 'V' => 'VAL', 'L' => 'LEU',
	'I' => 'ILE', 'P' => 'PRO', 'C' => 'CYS', 'M' => 'MET',
	'H' => 'HIS', 'F' => 'PHE', 'Y' => 'TYR', 'W' => 'TRP',
	'N' => 'ASN', 'Q' => 'GLN', 'S' => 'SER', 'T' => 'THR',
	'K' => 'LYS', 'R' => 'ARG', 'D' => 'ASP', 'E' => 'GLU' );

my %three_to_one = (
	'GLY' => 'G', 'ALA' => 'A', 'VAL' => 'V', 'LEU' => 'L', 'ILE' => 'I',
	'PRO' => 'P', 'CYS' => 'C', 'MET' => 'M', 'HIS' => 'H', 'PHE' => 'F',
	'TYR' => 'Y', 'TRP' => 'W', 'ASN' => 'N', 'GLN' => 'Q', 'SER' => 'S',
	'THR' => 'T', 'LYS' => 'K', 'ARG' => 'R', 'ASP' => 'D', 'GLU' => 'E',
	# nonstd
	'5HP' => 'Q', 'ABA' => 'C', 'AGM' => 'R', 'CEA' => 'C', 'CGU' => 'E',
	'CME' => 'C', 'CSB' => 'C', 'CSE' => 'C', 'CSD' => 'C', 'CSO' => 'C',
	'CSP' => 'C', 'CSS' => 'C', 'CSW' => 'C', 'CSX' => 'C', 'CXM' => 'M',
	'CYM' => 'C', 'CYG' => 'C', 'DOH' => 'D', 'FME' => 'M', 'GL3' => 'G',
	'HYP' => 'P', 'KCX' => 'K', 'LLP' => 'K', 'LYZ' => 'K', 'MEN' => 'N',
	'MGN' => 'Q', 'MHS' => 'H', 'MIS' => 'S', 'MLY' => 'K', 'MSE' => 'M',
	'NEP' => 'H', 'OCS' => 'C', 'PCA' => 'Q', 'PTR' => 'Y', 'SAC' => 'S',
	'SEP' => 'S', 'SMC' => 'C', 'STY' => 'Y', 'SVA' => 'S', 'TPO' => 'T',
	'TPQ' => 'Y', 'TRN' => 'W', 'TRO' => 'W', 'YOF' => 'Y'
);

###############################
###############################

## main()
if ($#ARGV < 2) {
	print STDERR "usage: $0 <fasta-file> <silent_file> <working-dir>\n";
	exit -1;
}

my $fastafile  = shift @ARGV;
my $silentfile = shift @ARGV;
my $WORKDIR    = shift @ARGV;

## files inside WORKDIR
$CSTFILENAMES = "$WORKDIR/".$CSTFILENAMES;
$FLAGFILENAMES = "$WORKDIR/".$FLAGFILENAMES;
$XMLFILENAMES  = "$WORKDIR/".$XMLFILENAMES;
$MAPFILENAMES  = "$WORKDIR/".$MAPFILENAMES;
$ALNTEMPLATEDIR  = "$WORKDIR/".$ALNTEMPLATEDIR;
$UNALNTEMPLATEDIR  = "$WORKDIR/".$UNALNTEMPLATEDIR;


## read sequence
open (FASTA, $fastafile) || die "Cannot open $_";
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


## 

###############################
###############################
## STAGE 0 -- make output directories if they don't exist
rmtree($WORKDIR);
mkdir $WORKDIR;
mkdir $UNALNTEMPLATEDIR;

# extract
chdir $UNALNTEMPLATEDIR;
my $cmd = "$EXTRACTAPP -in:file:silent ../../$silentfile -in:file:silent_struct_type binary";
print STDERR $cmd."\n";
system($cmd);
chdir "../..";


###############################
###############################
## STAGE 2 -- clustering + alignment
##
my @THREADED_MDLS=<$UNALNTEMPLATEDIR/*.pdb>;
my %allfrags;
my %bbatoms;
my %resids;
my %fraglens;
my %ss;

# read in all PDBs
my $minres = 9999;
my $maxres = 0;
foreach my $pdb (@THREADED_MDLS) {
	print STDERR $pdb."\n";
	open (PDB, $pdb) || die "Cannot open $pdb";

	$bbatoms{$pdb} = {};
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

		# only care about bb heavy atoms + CB
		next if ($atom ne " CA " && $atom ne " C  " && $atom ne " O  " && $atom ne " N  " && $atom ne " CB ");

		my $x = substr ($line, 30, 8);
		my $y = substr ($line, 38, 8);
		my $z = substr ($line, 46, 8);

		my $id = $atom.$chain.$res;
		$bbatoms{$pdb}->{ $id } = [$x,$y,$z];
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
	foreach ( keys %{ $bbatoms{$THREADED_MDLS[$i]} } ) {
		push @common_atoms, $_ if exists $bbatoms{$THREADED_MDLS[$j]}->{$_};
	}
	my $ncommon_atoms = scalar( @common_atoms );

	if ($ncommon_atoms < 6) {
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
		push @{ $atoms_i }, deep_copy( $bbatoms{$THREADED_MDLS[$i]}->{$_} );
		push @{ $atoms_j }, deep_copy( $bbatoms{$THREADED_MDLS[$j]}->{$_} );
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
			my $x = $bbatoms{$THREADED_MDLS[$i]}->{$_};
			my $y = vadd( mapply( $Rs->[$i][$j] , vsub( $bbatoms{$THREADED_MDLS[$j]}->{$_} , $comis->[$i][$j] ) ), 
			              vadd($comis->[$i][$j] , $comjs->[$i][$j]) );
			my $dist = dist( $x,$y );
			if ($dist <= $RMS_ALIGN) {
				push @new_common_atoms, $_;
				push @{ $atoms_i }, deep_copy( $bbatoms{$THREADED_MDLS[$i]}->{$_} );
				push @{ $atoms_j }, deep_copy( $bbatoms{$THREADED_MDLS[$j]}->{$_} );
			}
		}

		$ncommon_atoms = scalar( @new_common_atoms );

		my $natoms_tot = min( scalar(keys %{ $bbatoms{$THREADED_MDLS[$i]} }),  scalar(keys %{ $bbatoms{$THREADED_MDLS[$j]} }) );

		# penalize very short alignments
		if ($RMS_ALIGN == $RMS_CUTOFFS[2]) {
			$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $natoms_tot;
			print STDERR "overlap(".$THREADED_MDLS[$i].",".$THREADED_MDLS[$j].") over ".@common_atoms."/".$natoms_tot." atoms is ".$overlapscore->[$i][$j]."\n";
		}

		if ($ncommon_atoms < 6) {
		 	$overlapscore->[$i][$j] = 0;
			next;
		}

		last if (($ncommon_atoms < ($ALIGNCUTOFF*$natoms_tot)) && ($RMS_ALIGN<=$RMS_CUTOFFS[2]));

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
print STDERR "SEED $minI (".$THREADED_MDLS[$minI].") with probability $maxprob\n";

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
	print STDERR "ALIGN ".$THREADED_MDLS[$minJ]." (from ".$THREADED_MDLS[$minI].") with overlap ".$overlapscore->[$minI][$minJ]."\n";
	$globalRs->[$minJ] = mmult( $globalRs->[$minI], $Rs->[$minI][$minJ]  );
	$preTs->[$minJ] = deep_copy( $comis->[$minI][$minJ] );

	# check to see if this belongs in a new cluster
	if ($overlapscore->[$minI][$minJ] >= $CLUSTERCUTOFF) {
		$clusterid->[$minJ] = $clusterid->[$minI];
	} else {
		$clusterid->[$minJ] = $nextclusterid;
		$nextclusterid++;
	}

	# post-rotation translation is a bit tricky
	# we need to transform to j's coordinate frame, then apply j's transformation to this result
	my $post_to_j = vadd( $comjs->[$minI][$minJ], $preTs->[$minJ] );
	my $post_j_to_global = vadd( mapply( $globalRs->[$minI] , vsub( $post_to_j , $preTs->[$minI]) ) , 
	                             $postTs->[$minI] );
	$postTs->[$minJ] = $post_j_to_global;

	$aligned[$minJ] = 1;
}

# sort clusters by p(correct)
my %cluster_members;
my @clusterprobs;
foreach my $i (1..$nextclusterid) {
	my @currclust = ();
	$clusterprobs[$i-1] = 0;
	foreach my $j (0..$#THREADED_MDLS) {
		if ($clusterid->[$j] == $i) {
			push @currclust, $j;
			my $prob = 1.0;
			$clusterprobs[$i-1] += $prob;
		}
	}
	$cluster_members{ $clusterprobs[$i-1] } = \@currclust;
}

@clusterprobs = sort {$b <=> $a} @clusterprobs;
my $counter = 0;
my $nclust = 0;
print STDERR "Using $nclust clusters\n";

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

###############################
###############################
## STAGE 3 -- output
mkdir ($WORKDIR);

# aligned templates
mkdir ($ALNTEMPLATEDIR);
foreach my $i (0..$#THREADED_MDLS) {
	next if ($clusterid->[$i] == -1);
	my $clid = $clusterid->[$i];

	my $pdb = $THREADED_MDLS[$i];
	my $nfrags = scalar( @{ $allfrags{$pdb} } );

	my $outpdb = $pdb;
	$outpdb =~ s/.*\///;
	$outpdb =~ s/\.pdb$/_aln.cl$clid.pdb/;
	print STDERR "writing aligned_templates/$outpdb\n";
	open (PDBF, ">$ALNTEMPLATEDIR/$outpdb") || die "Cannot open $_";

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

# thesius
$cmd = "$THESEUSAPP $ALNTEMPLATEDIR/\*.pdb";
system($cmd);
open (THES, 'theseus_variances.txt') || die "Cannot open $_";
my @thes_vars = <THES>;
close(THES);

my %per_res_var;
foreach my $line (@thes_vars) {
	# RES 1           MET      1    77.520451     8.804570    21.785661 CORE
	if ($line =~/RES.*CORE/) {
		my @fields = split ' ', $line;
		$per_res_var{ int($fields[3]) } = $fields[4];
	}
}

#trim
foreach my $i (0..$#THREADED_MDLS) {
	next if ($clusterid->[$i] == -1);
	my $clid = $clusterid->[$i];
	my $pdb = $THREADED_MDLS[$i];
	my $nfrags = scalar( @{ $allfrags{$pdb} } );
	my $outpdb = $pdb;
	$outpdb =~ s/.*\///;
	$outpdb =~ s/\.pdb$/_aln.cl$clid.pdb/;
	print STDERR "trimming aligned_templates/$outpdb\n";
	open (PDBIN, "$ALNTEMPLATEDIR/$outpdb") || die "Cannot open $_";
	my @oldlines = <PDBIN>;
	close PDBIN;

	open (PDBOUT, ">$ALNTEMPLATEDIR/$outpdb") || die "Cannot open $_";
	foreach my $line (@oldlines) {
		my $resid = substr ($line, 22, 4);
		next if ($per_res_var{ int($resid) } >  $THESEUSCUT);
		print PDBOUT $line;
	}
	close PDBOUT;
}


#trim
unlink("theseus_ave.pdb");
unlink("theseus_residuals.txt");
unlink("theseus_sup.pdb");
unlink("theseus_transf2.txt");
unlink("theseus_tree.nxs");
unlink("theseus_variances.txt");


# (d) template density maps
foreach my $i (0..$nclust) {
	my $mapfile = sprintf $MAPFILENAMES , $i;

	my @input_pdbs = ();
	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);

		my $temppdb = $mapfile;
		$temppdb =~ s/\.mrc/.$j.pdb/g;
		open (PDBCAT, ">$temppdb") || die "Cannot open $_";

		my $tmpl_filename = $THREADED_MDLS[ $j ];
		open (TEMPLPDB, "$tmpl_filename") || die "Cannot open $_";
		my @templatepdb_lines = <TEMPLPDB>;
		close TEMPLPDB;

		foreach my $line (@templatepdb_lines) {
			next if ($line !~ /^ATOM/);
			my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
			$newX = vsub( $newX, $preTs->[$j]);
			$newX = mapply( $globalRs->[$j], $newX );
			$newX = vadd( $newX, $postTs->[$j]);
			substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
			print PDBCAT $line;
		}
		close PDBCAT;
		push @input_pdbs, $temppdb;
	}

	my $outmap = sprintf $MAPFILENAMES, $i;
	my $cmd = "$DENSGENAPP -s ".(join ' ',@input_pdbs)." -edensity::mapfile $outmap";
	print STDERR $cmd."\n";
	system($cmd);

	foreach my $file (@input_pdbs) {
		unlink($file);
	}
}

# (e) config/flags file
my $dir = getcwd;
# local hack
$dir =~ s/\/gpfs\/DS3524-1//g;
$dir =~ s/WORK/work/g;
my $dirtag = $dir;
$dirtag =~ s/.*\/([^\/]+$)/$1/g;

foreach my $i (0..$nclust) {
	my $xmlfile = sprintf $XMLFILENAMES, $i;
	open (XML, ">$xmlfile") || die "Cannot open $_";

	print XML  "<dock_design>\n";
	print XML  "    <TASKOPERATIONS>\n";
	print XML  "    </TASKOPERATIONS>\n";
	print XML  "    <SCOREFXNS>\n";
	print XML  "        <fullatom weights=stage3_rlx>\n";
	print XML  "	    </fullatom>\n";
	print XML  "        <stage1 weights=stage1>\n";
	print XML  "        </stage1>\n";
	print XML  "        <stage2 weights=stage2>\n";
	print XML  "        </stage2>\n";
	print XML  "    </SCOREFXNS>\n";
	print XML  "    <FILTERS>\n";
	print XML  "    </FILTERS>\n";
	print XML  "    <MOVERS>\n";
	print XML  "        <Hybridize name=hybridize stage1_scorefxn=stage1 stage2_scorefxn=stage2 fa_scorefxn=fullatom $HYBRIDIZEOPTIONS>\n";
 	my $frag3filepath = sprintf "$dir/$FRAG3FILE", $dirtag;
 	my $frag9filepath = sprintf "$dir/$FRAG9FILE", $dirtag;
	print XML  "            <Fragments 3mers=\"$frag3filepath\" 9mers=\"$frag9filepath\"/>\n";

	my $cstfile = sprintf $CSTFILENAMES, $i;
	$cstfile = $dir."/$cstfile";

 	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);
		my $id = $THREADED_MDLS[$j];
		$id =~ s/.*\///;
		$id =~ s/\.pdb$/_aln.cl$i.pdb/;
		$id = "$dir/$ALNTEMPLATEDIR/$id";

		my $prob = 1;

	 	my $outline = sprintf "pdb=\"$id\" cst_file=\"$cstfile\" weight=$prob";
		print XML  "            <Template $outline/>\n";
	}
	print XML  "        </Hybridize>\n";
	print XML  "    </MOVERS>\n";
	print XML  "    <APPLY_TO_POSE>\n";
	print XML  "    </APPLY_TO_POSE>\n";
	print XML  "    <PROTOCOLS>\n";
	print XML  "        <Add mover=hybridize/>\n";
	print XML  "    </PROTOCOLS>\n";
	print XML  "</dock_design>\n";
	close(XML);

 	# finally write flags
 	my $flagfile = sprintf $FLAGFILENAMES, $i;
 	open (FLAGS, ">$flagfile") || die "Cannot open $_";
 	print FLAGS "-in:file:fasta $dir/$fastafile\n";
 	print FLAGS "-out:prefix $dirtag\n";
 	print FLAGS "-out:file:silent $dir/../$dirtag.cl$i.silent\n";
 	print FLAGS "-in:file:native $dir/$NATIVEDIR/$dirtag.pdb\n";
 	print FLAGS "-parser:protocol $dir/$xmlfile\n";
 	printf FLAGS "-edensity::mapfile $dir/$MAPFILENAMES\n", $i;
}

exit 0;

###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}
