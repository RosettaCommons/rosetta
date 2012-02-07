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
my $PARTIALTHREAD_APP = "/work/dimaio/rosetta/rosetta_source/bin/partial_thread.default.linuxgccrelease".
                        " -database /work/dimaio/rosetta/rosetta_database ".
						"-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer  VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm";
my $DENSGENAPP = "/work/dimaio/rosetta/rosetta_source/bin/pdb_to_map.default.linuxgccrelease".
                 " -database /work/dimaio/rosetta/rosetta_database ".
                 " -edensity:mapreso 5.0 -edensity::grid_spacing 2.0 -in:file:centroid_input";
my $CSTGEN_APP = "/work/dimaio/rosetta/cm_scripts/bin/predict_distances.pl";
my $PDBFILEDIR = "/lab/databases/wwpdb/";


## paths
my $NATIVEDIR        = "native/";
my $TEMPLATEDIR      = "templates/";
my $PARTIALTHREADDIR = "partial_threads/";
my $COORDCSTDIR      = "coordCsts_resOnly/";
my $FRAG3FILE      = "fragments/%s_templatesvall.200.3mers";
my $FRAG9FILE      = "fragments/%s_templatesvall.200.9mers";
my $SSPREDFILE      = "fragments/t000_.psipred_ss2";

my $PCORRFILENAMES  = "p_correct%d.txt";
my $ALNFILENAMES  = "cluster%d.filt";
my $CSTFILENAMES  = "cluster%d.filt.dist_csts";
my $FLAGFILENAMES = "flags%d";
my $MAPFILENAMES  = "templates%d.mrc";
my $ALNTEMPLATEDIR  = "aligned_templates";

my $XMLFILENAMES = "hybridize%d.xml";
my $RUNFILENAMES = "run%d.sh";
my $HYBRIDIZEOPTIONS = "batch=16 stage1_increase_cycles=1.0 stage2_increase_cycles=0.25 linmin_only=1";

#####
#####

## structure alignments
my @RMS_CUTOFFS = (10,5,4,3,2,1.5,1);
my $CLUSTERCUTOFF = 0.40;
my $ALIGNCUTOFF   = 0.20;  # to get better superpositions, trade coverage for alignment

## probability correct
## this could be read from a file
my %template_probs =  (
	201 => 0.15481, 301 => 0.13389, 401 => 0.16318,
	202 => 0.03347, 302 => 0.09623, 402 => 0.03766, 
	203 => 0.02510, 303 => 0.05858, 403 => 0.05021,
	204 => 0.01255, 304 => 0.05439, 404 => 0.02092,
	205 => 0.00837, 305 => 0.02929, 405 => 0.01674,
	206 => 0.00418, 306 => 0.00837, 406 => 0.00418, 
	                307 => 0.00418, 407 => 0.01255,
	                308 => 0.01255, 408 => 0.01255,
	                309 => 0.01255, 409 => 0.01255,
	                310 => 0.00837, 410 => 0.01255,
	0 => 0.00418    # default
);
my $template_M_EST = 0.010; # favor lower-ranked alignments a bit more in sampling (and cst-gen?)

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
	#print STDERR "usage: $0 <fasta-file> <ali-file> <evmap> <working-dir>\n";
	print STDERR "usage: $0 <fasta-file> <ali-file> <working-dir>\n";
	exit -1;
}

my $fastafile = shift @ARGV;
my $alifile   = shift @ARGV;
#my $evmapfile = shift @ARGV;
my $WORKDIR   = shift @ARGV;

## files inside WORKDIR
$PCORRFILENAMES  = "$WORKDIR/".$PCORRFILENAMES;
$ALNFILENAMES  = "$WORKDIR/".$ALNFILENAMES;
$CSTFILENAMES  = "$WORKDIR/".$CSTFILENAMES;
$FLAGFILENAMES = "$WORKDIR/".$FLAGFILENAMES;
$XMLFILENAMES  = "$WORKDIR/".$XMLFILENAMES;
$MAPFILENAMES  = "$WORKDIR/".$MAPFILENAMES;
$ALNTEMPLATEDIR  = "$WORKDIR/".$ALNTEMPLATEDIR;


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


## real ali file
open (ALI, $alifile) || die "Cannot open $_";
my @alilines = <ALI>;
close ALI;
my %alimap = ();
my @currbuff; 
my $currtag = "";
foreach my $line (@alilines) {
	if ($line =~ /^##/) {
		chomp $line;
		if ($currtag ne "") {
			$alimap {$currtag} = dclone(\@currbuff);
		}
		my @fields = split / /, $line;
		$currtag = $fields[$#fields];
		@currbuff = ( $line."\n" );
	} else {
		push @currbuff, $line;
	}
}
if ($currtag ne "") {
	$alimap {$currtag} = dclone(\@currbuff);
}

###############################
###############################
## STAGE 0 -- make output directories if they don't exist
#rmtree($TEMPLATEDIR);
rmtree($PARTIALTHREADDIR);
rmtree($WORKDIR);
mkdir $TEMPLATEDIR;
mkdir $PARTIALTHREADDIR;
mkdir $WORKDIR;


###############################
###############################
## STAGE 1 -- generate partial threads
## (a) grab+sanitize templates
foreach my $tag (keys %alimap) {
	my $template = substr($tag,0,5);
	my $pdbout = $TEMPLATEDIR."/".$template.".pdb";
	if ( ! -f $pdbout ) {
		print STDERR "trying to get $pdbout!\n";
		my $pdbid = substr( $template, 0, 4 );
		my $chain = substr( $template, 4, 1 );
		my $dirid = substr( $template, 1, 2 );
		open (PDB, $PDBFILEDIR."/".$dirid."/".$pdbid.".pdb") || die "Cannot open $_";
		my @pdblines = <PDB>;
		close PDB;
		open (PDBOUT, ">$pdbout") || die "Cannot open $_";
		my $linecount = 0;
		foreach my $line (@pdblines) {
			last if ($line =~ /^ENDMDL/ && $linecount>0);
			next if ($line !~/^ATOM  / && $line !~/^HETATM/);
			my $chainid = substr($line,21,1);
			next if ($chainid ne $chain);
			my $resname = substr($line,17,3);
			next if (!defined $three_to_one{ $resname });
			my $conf = substr($line,16,1);
			next if ($conf ne " " && $conf ne "A");
			# sanitization
			substr($line,0,6) = "ATOM  ";
			substr($line,17,3) = $one_to_three{ $three_to_one{ $resname } };
			# MSE
			if ($resname eq "MSE") {
				my $atomname = substr ($line,12,4);
				if ($atomname eq "SE  ") {
					substr ($line,12,4) = " SD ";
				}
			}
			# all other non-standard ... throw away sidechain
			elsif (substr($line,17,3) ne $resname) {
				my $atomname = substr ($line,12,4);
				next if ($atomname ne " C  " && $atomname ne " CA "  && $atomname ne " O  "  && $atomname ne " N  "  && $atomname ne " CB ");
			}
			print PDBOUT $line;
			$linecount++;
		}
		close PDBOUT
	}
}

## (b) partial threading app
chdir $PARTIALTHREADDIR;
my $cmd = "$PARTIALTHREAD_APP -in::file::fasta ../$fastafile -in::file::alignment ../$alifile -in::file::template_pdb ../$TEMPLATEDIR/*.pdb";
print STDERR $cmd."\n";
system($cmd);
chdir "..";

###############################
###############################
## STAGE 2 -- clustering + alignment
##
my @THREADED_MDLS=<$PARTIALTHREADDIR/*.pdb>;
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
		if ($RMS_ALIGN == $RMS_CUTOFFS[2]) {
			$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $natoms_tot;
			print STDERR "overlap(".$THREADED_MDLS[$i].",".$THREADED_MDLS[$j].") over ".@common_atoms."/".$natoms_tot." atoms is ".$overlapscore->[$i][$j]."\n";
		}

		last if ($ncommon_atoms < 6); #overlap == 0
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
my $maxprob = 0;
my $minI = -1;
foreach my $i (0..$#THREADED_MDLS) {
	my $tag = $THREADED_MDLS[$i]; $tag =~ s/.*_(\d\d\d).*/$1/;
	my $prob = $template_probs{ $tag };
	$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
	if ($prob > $maxprob) {
		$maxprob = $prob;
		$minI = $i;
	}
}
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
			my $tag = $THREADED_MDLS[$j]; $tag =~ s/.*_(\d\d\d).*/$1/;
			my $prob = $template_probs{ $tag };
			$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
			$clusterprobs[$i-1] += $prob;
		}
	}
	$cluster_members{ $clusterprobs[$i-1] } = \@currclust;
}

# find out how many clusters to use
@clusterprobs = sort {$b <=> $a} @clusterprobs;
my $counter = 0;
my $nclust = 0;
foreach my $prob (@clusterprobs) {
	foreach my $member (@{ $cluster_members{ $prob } }) {
		my $filename = $THREADED_MDLS[ $member ];
		$filename =~ s/.*\///;
		if ( $filename =~ /201/ || $filename =~ /301/ || $filename =~ /401/ ) {
			$nclust = $counter;
		}
	}
	$counter++;
}
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

# (a) per-cluster alignment files
#      + normalized p-correct
foreach my $i (0..$nclust) {
	my $clstfile = sprintf $ALNFILENAMES, $i;
	open (CLST, ">$clstfile") || die "Cannot open $_";
	my $pcorrect = sprintf $PCORRFILENAMES, $i;
	open (PCORR, ">$pcorrect") || die "Cannot open $_";

	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);

		my $filename = $THREADED_MDLS[ $j ];
		$filename =~ s/.*\///;
		foreach my $alnline ( @{$alimap {substr($filename,0,9)}} ) {
			print CLST $alnline;
		}

		my $tag = $filename; $tag =~ s/.*_(\d\d\d).*/$1/;
		my $prob = $template_probs{ $tag };
		$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
		$prob = $prob / $clusterprobs[$i];
		print PCORR "$tag $prob\n";
	}
	$counter++;

	close CLST;
	close PCORR;

	# (b) constraint files
	#chdir $TEMPLATEDIR;
	mkdir "temp";
	chdir "temp";
	#~tex/src/cm_scripts/bin/predict_distances.pl $file $pdb/alignment/$pdb.fasta -aln_format grishin -weights_file /work/tex/src/cm_scripts/bin/p-correct.txt
	my $cmd = "$CSTGEN_APP ../$clstfile ../$fastafile -aln_format grishin -weights_file ../$pcorrect"; ## -ev_map_file ../$evmapfile";
	print STDERR $cmd."\n";
	system($cmd);
	chdir "..";
	rmtree("temp");
}

# (c) aligned templates
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

# (d) template density maps
foreach my $i (0..$nclust) {
	my $inputpdbs = "";
 	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);

		my $pdb = $THREADED_MDLS[$i];
		my $nfrags = scalar( @{ $allfrags{$pdb} } );

		my $outpdb = $pdb;
		$outpdb =~ s/.*\///;
		$outpdb =~ s/\.pdb$/_aln.cl$i.pdb/;
		$inputpdbs = $inputpdbs." $ALNTEMPLATEDIR/$outpdb";
	}
	my $outmap = sprintf $MAPFILENAMES, $i;

	my $cmd = "$DENSGENAPP -s $inputpdbs -edensity::mapfile $outmap";
	print STDERR $cmd."\n";
	system($cmd);
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
	print XML  "        <fullatom weights=score12_full>\n";
	print XML  "   	        <Reweight scoretype=cart_bonded weight=0.5/>\n";
	print XML  "   	        <Reweight scoretype=elec_dens_fast weight=2.0/>\n";
	print XML  "	    </fullatom>\n";
	print XML  "        <stage1>\n";
	print XML  "   	        <Reweight scoretype=env weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=pair weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=cbeta weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=cenpack weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=hs_pair weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=ss_pair weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=rsigma weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=sheet weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=vdw weight=0.20/>\n";
	print XML  "   	        <Reweight scoretype=rg weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=rama weight=0.3/>\n";
	print XML  "   	        <Reweight scoretype=atom_pair_constraint weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=elec_dens_fast weight=2.0/>\n";
	print XML  "	       </stage1>\n";
	print XML  "        <stage2>\n";
	print XML  "   	        <Reweight scoretype=hbond_sr_bb weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=hbond_lr_bb weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=rama weight=0.2/>\n";
	print XML  "   	        <Reweight scoretype=omega weight=0.2/>\n";
	print XML  "   	        <Reweight scoretype=rg weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=cen_env_smooth weight=2.0/>\n";
	print XML  "   	        <Reweight scoretype=cen_pair_smooth weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=cbeta_smooth weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=cenpack_smooth weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=vdw weight=1.0/>\n";
	print XML  "   	        <Reweight scoretype=atom_pair_constraint weight=0.5/>\n";
	print XML  "   	        <Reweight scoretype=cart_bonded weight=0.05/>\n";
	print XML  "   	        <Reweight scoretype=elec_dens_fast weight=2.0/>\n";
	print XML  "        </stage2>\n";
	print XML  "    </SCOREFXNS>\n";
	print XML  "    <FILTERS>\n";
	print XML  "    </FILTERS>\n";
	print XML  "    <MOVERS>\n";
	print XML  "        <Hybridize name=hybridize stage1_scorefxn=stage1 stage2_scorefxn=stage2 fa_scorefxn=fullatom $HYBRIDIZEOPTIONS>\n";
 	my $frag3filepath = sprintf "$dir/$FRAG3FILE", $dirtag;
 	my $frag9filepath = sprintf "$dir/$FRAG9FILE", $dirtag;
	print XML  "            <Fragments 3mers=\"$frag3filepath\" 9mers=\"$frag9filepath\"/>\n";

 	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);
		my $id = $THREADED_MDLS[$j];
		$id =~ s/.*\///;
		$id =~ s/\.pdb$/_aln.cl$i.pdb/;
		$id = "$dir/$ALNTEMPLATEDIR/$id";

		my $cstfile = sprintf $CSTFILENAMES, $i;
		$cstfile = $dir."/$cstfile";

		my $tag = $THREADED_MDLS[$j]; $tag =~ s/.*_(\d\d\d).*/$1/;
		my $prob = $template_probs{ $tag };
		$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
		$prob += $template_M_EST;

	 	my $outline = sprintf "pdb=\"$id\" cst_file=\"$cstfile\" weight=$prob";
		print XML  "			            <Template $outline/>\n";
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



# foreach my $i (0..$nclust) {
# 	my $cfgfile = sprintf $CFGFILENAMES, $i;
# 	my $flagfile = sprintf $FLAGFILENAMES, $i;
# 	open (HYB, ">$cfgfile") || die "Cannot open $_";
# 	open (FLAGS, ">$flagfile") || die "Cannot open $_";
# 
# 	foreach my $j (0..$#THREADED_MDLS) {
# 		next if ($clusterid->[$j] != $i);
# 
# 		my $id = $THREADED_MDLS[$j];
# 		$id =~ s/.*\///;
# 		$id =~ s/\.pdb$/_aln.cl$i.pdb/;
# 		$id = "$dir/$ALNTEMPLATEDIR/$id";
# 
# 		my $cstfile = sprintf $CSTFILENAMES, $i;
# 		$cstfile = $dir."/$cstfile";
# 
# 		my $tag = $THREADED_MDLS[$j]; $tag =~ s/.*_(\d\d\d).*/$1/;
# 		my $prob = $template_probs{ $tag };
# 		$prob = $template_probs{ 0 } if !defined $template_probs{ $tag };
# 		$prob += $template_M_EST;
# 
# 	 	my $outline = sprintf "%20s %30s %4d %.5f  ", $id, $cstfile, $i, $prob;
# 
# 		# load coord csts
# 		my $coordcstfile = $COORDCSTDIR."/".substr( $id,0,9 ).".coordCsts";
# 
# 		# stupid renumbering issue
# 		my $coordcstfileAlt = $coordcstfile;
# 		my $newtag = $tag-1;
# 		$coordcstfileAlt =~ s/_$tag/_$newtag/g;
# 
# 		if (-e $coordcstfile) {
# 			open (RSD, $coordcstfile) || die "Cannot open $_";
# 			my @resids = <RSD>;
# 			close RSD;
# 			chomp (@resids);
# 			$outline = $outline.(join ',',@resids);
# 		} elsif (-e $coordcstfileAlt) {
# 			open (RSD, $coordcstfileAlt) || die "Cannot open $_";
# 			my @resids = <RSD>;
# 			close RSD;
# 			chomp (@resids);
# 			$outline = $outline.(join ',',@resids);
# 		}
# 		print HYB $outline."\n";
# 	}
# 
# 	# finally write flags
# 	print FLAGS "-in:file:fasta $dir/$fastafile\n";
# 	print FLAGS "-out:prefix $dirtag\n";
# 	print FLAGS "-out:file:silent $dir/../$dirtag.cl$i.silent\n";
# 	print FLAGS "-in:file:native $dir/$NATIVEDIR/$dirtag.pdb\n";
# 	print FLAGS "-cm::hybridize::template_list $dir/$cfgfile\n";
# 	printf FLAGS "-in::file::frag3 $dir/$FRAG3FILE\n", $dirtag;
# 	printf FLAGS "-in::file::frag9 $dir/$FRAG9FILE\n", $dirtag;
# 	printf FLAGS "-edensity::mapfile $dir/$MAPFILENAMES\n", $i;
# 	print FLAGS "-in:file:psipred_ss2 $SSPREDFILE\n";
# }

exit 0;

###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}
