#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);

use lib (".");
use File::Basename;
use File::Path qw(rmtree);
use File::Copy;
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
my $SYMMDEFFILE_APP = "/work/dimaio/rosetta/rosetta_source/src/apps/pilot/frank/make_NCS.pl -r 999";
my $DDG_APP = "/work/dimaio/rosetta/rosetta_source/bin/ddg.default.linuxgccrelease".
              " -database /work/dimaio/rosetta/rosetta_database -score:weights softrep";
my $DDG_CUTOFF = -10;

my $PDBFILEDIR = "/lab/databases/wwpdb/";
my $BIOPDBFILEDIR = "/lab/databases/biounits/";

## paths
my $NATIVEDIR        = "native/";
my $TEMPLATEDIR      = "templates/";
my $PARTIALTHREADDIR = "partial_threads/";
my $SYMMDIR          = "symmetry/";
my $COORDCSTDIR      = "coordCsts_resOnly/";
my $FRAG3FILE      = "fragments/%s_quicktemplatesvall.25.3mers";
my $FRAG9FILE      = "fragments/%s_quicktemplatesvall.25.9mers";

my $PCORRFILENAMES  = "p_correct%d.txt";
my $ALNFILENAMES  = "cluster%d.filt";
my $CSTFILENAMES  = "cluster%d.filt.dist_csts";
my $FLAGFILENAMES = "flags%d_%s";
my $MAPFILENAMES  = "templates%d.mrc";
my $ALNTEMPLATEDIR  = "aligned_templates";
my $SYMMDEFDIR  = "symmdef";

#my $CFGFILENAMES  = "hybrid%d.config";
my $XMLFILENAMES = "hybridize%d_%s.xml";
my $RUNFILENAMES = "run%d.sh";
my $HYBRIDIZEOPTIONS = "batch=2 stage1_increase_cycles=1.0 stage2_increase_cycles=1.0 linmin_only=0";

#####
#####

## structure alignments
my @RMS_CUTOFFS = (10,5,4,3,2,1.5,1);
my $CLUSTERCUTOFF = 0.40;
my $ALIGNCUTOFF   = 0.20;  # to get better superpositions, trade coverage for alignment

## symmetry
my $MINSYMMPROB = 0.10;

## probability correct
## this could be read from a file
my %template_probs =  (
	101 =>0.1430503, 102 => 0.0770168, 103 => 0.0427618, 104 => 0.0249920, 105 => 0.0157738,
	106 =>0.0109919, 107 => 0.0085113, 108 => 0.0072245, 109 => 0.0065569, 110 => 0.0062107,

	201 =>0.1276504, 202 => 0.0535388, 203 => 0.0243316, 204 => 0.0128212, 205 => 0.0082849,
	206 =>0.0064972, 207 => 0.0057927, 208 => 0.0055150, 209 => 0.0054056, 210 => 0.0053625,

	301 =>0.1161114, 302 => 0.0838580, 303 => 0.0609236, 304 => 0.0446156, 305 => 0.0330195,
	306 =>0.0247739, 307 => 0.0189107, 308 => 0.0147415, 309 => 0.0117769, 310 => 0.0096689,

	401 =>0.1430503, 402 => 0.0770168, 403 => 0.0427618, 404 => 0.0249920, 405 => 0.0157738,
	406 =>0.0109919, 407 => 0.0085113, 408 => 0.0072245, 409 => 0.0065569, 410 => 0.0062107,

	501 =>0.0666667, 502 => 0.0666667, 503 => 0.0666667, 504 => 0.0666667, 505 => 0.0666667,
	506 =>0.0666667, 507 => 0.0666667, 508 => 0.0666667, 509 => 0.0666667, 510 => 0.0666667,

	0 => 0.005    # default
);
my $template_M_EST = 0.01; # favor lower-ranked alignments a bit more in sampling
my $template_M_EST_counts = 1; # when using custom counts to weight, use this M instead

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
$SYMMDEFDIR = "$WORKDIR/".$SYMMDEFDIR;

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
rmtree($SYMMDIR);
rmtree($PARTIALTHREADDIR);
rmtree($WORKDIR);
mkdir $SYMMDIR;
mkdir $TEMPLATEDIR;
mkdir $PARTIALTHREADDIR;
mkdir $WORKDIR;


###############################
###############################
## STAGE 1 -- generate partial threads
## (a) grab+sanitize templates
foreach my $tag (keys %alimap) {
	my ($template,$idtag) = split '_', $tag;
	my $pdbid = substr( $template, 0, 4 );
	my $tgtchain = substr( $template, 4, 1 );
	my $dirid = substr( $template, 1, 2 );

	print STDERR "At  $tag\n";

	my $pdbout1 = $TEMPLATEDIR."/".$template.".pdb";
	unless ( ( -e $pdbout1 ) ) {
		print STDERR "writing $pdbout1!\n";
		open (PDB, $PDBFILEDIR."/".$dirid."/".$pdbid.".pdb") || die "Cannot open $_";
		my @pdblines = <PDB>;
		close PDB;

		open (PDBOUT1, ">$pdbout1") || die "Cannot open $_";
		my $linecount = 0;
		foreach my $line (@pdblines) {
			last if ($line =~ /^ENDMDL/ && $linecount>0);
			next if ($line !~/^ATOM  / && $line !~/^HETATM/);

			my $chainid = substr($line,21,1);
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

			if ($chainid eq $tgtchain) {
				print PDBOUT1 $line;
			}

			$linecount++;
		}
		close PDBOUT1;
	}

	# only consider symmetry of the top 3 templates
	$idtag =~ s/\.pdb//g;
	next if (($idtag % 100) > 3);

	# grab biounits
	my @allbiofiles = <$BIOPDBFILEDIR/$dirid/$pdbid.pdb*>;
	my (@bestbiofile,@bestbiofile_chains);
	my $mostbiochains=1;
	foreach my $biofile (@allbiofiles) {
		# size filter ... 20 MB gzipped???
		# mainly for 3k1qA
		next if ( (-s $biofile) > 20*1024*1024);

		print STDERR "Unzip $biofile\n";
		open( BIOPDB, "-|", "zcat " . $biofile) || die "!";
		my @pdblines = <BIOPDB>;
		close BIOPDB;

		# verify that it contains the target chain
		# remember each chain's fasta
		my %fastas;
		my $mdl = 0;
		foreach my $line (@pdblines) {
			$mdl++ if ($line =~/^ENDMDL/);

			next if ($line !~/^ATOM  / && $line !~/^HETATM/);
			my $resname = substr($line,17,3);
			next if (!defined $three_to_one{ $resname });
			my $conf = substr($line,16,1);
			next if ($conf ne " " && $conf ne "A");

			my $atomname = substr ($line,12,4);
			next if ($atomname ne " CA ");

			my $chainid = substr($line,21,1).$mdl;
			if (!defined ($fastas{$chainid})) {
				$fastas{$chainid} = "";
			}
			$fastas{$chainid} = $fastas{$chainid}.$three_to_one{ $resname };
		}

		# make sure chain exists
		next if !defined( $fastas{$tgtchain.'0'} ); # chain is not part of biounit

		# throw out nonidentical chains
		my @chain_matches = ($tgtchain.'0');
		foreach my $chainid (keys %fastas) {
			next if ($tgtchain eq $chainid);
			my ($tgt_aln, $src_aln, $pdbstart, $maxscore) = alignSW( $fastas{$tgtchain.'0'}, $fastas{$chainid} );
			#my $target = realign( $tgt_aln, $fastas{$tgtchain}, $fastas{$chainid} ); # attempt to remap target seq

			next if ($maxscore <= 0); # <50% aln res
			push @chain_matches, $chainid;
		}

		my $nbiochains = scalar( @chain_matches );
		if ($nbiochains > $mostbiochains) {
			$mostbiochains = $nbiochains;
			@bestbiofile = @pdblines;
			@bestbiofile_chains = @chain_matches;
		}
	}

	next if ($mostbiochains <= 1);
	next if ($mostbiochains > 62); # ...

	my %chainhash;
	my $chainids = "BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()";
	my $currchain = 0;
	foreach my $chainid (@bestbiofile_chains) {
		if ($chainid eq $tgtchain.'0') {
			$chainhash{$chainid} = 'A';
		} else {
			$chainhash{$chainid} = substr($chainids,$currchain,1);
			$currchain++;
		}
	}

	my $pdbout2 = $SYMMDIR."/".$template.".pdb";
	#print STDERR "writing $pdbout2!\n";
	open (PDBOUT2, ">$pdbout2") || die "Cannot open $_";
	my $linecount = 0;
	my $mdl = 0;
	foreach my $line (@bestbiofile) {
		$mdl++ if ($line =~/^ENDMDL/);
		next if ($line !~/^ATOM  / && $line !~/^HETATM/);

		my $chainid = substr($line,21,1);
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

		if (defined $chainhash{$chainid.$mdl}) {
			substr($line,21,1) = $chainhash{$chainid.$mdl};
			print PDBOUT2 $line;
		}


		$linecount++;
	}
	#print STDERR "   ... done\n";
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
my $maxprob = 0;
my $minI = -1;
foreach my $i (0..$#THREADED_MDLS) {
	my $tag = $THREADED_MDLS[$i];
	my $prob = $template_probs{ 0 };
	if ($tag =~ /.*_(\d\d\d).*/) {
		$tag = $1;
		$prob = $template_probs{ $tag } if defined $template_probs{ $tag };
	} elsif ($tag =~ /.*_w(\d+).*/) {
		$prob = $1+$template_M_EST_counts;
	}
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
			my $tag = $THREADED_MDLS[$j];
			my $prob = $template_probs{ 0 };
			if ($tag =~ /.*_(\d\d\d).*/) {
				$tag = $1;
				$prob = $template_probs{ $tag } if defined $template_probs{ $tag };
			} elsif ($tag =~ /.*_w(\d+).*/) {
				$prob = $1+$template_M_EST_counts;
			}
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
		if ( $filename =~ /_101/ || $filename =~ /_201/ || $filename =~ /_301/ || $filename =~ /_401/ ) {
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
## STAGE 2.5 -- symmetry analysis
my %symdeffiles;
foreach my $i (0..$nclust) {
	my @input_pdbs = ();
	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);

		my $tmpl_filename = $THREADED_MDLS[ $j ];
		$tmpl_filename =~ s/$PARTIALTHREADDIR/$SYMMDIR/g;
		my $tmpl_filename_in = $tmpl_filename;
		$tmpl_filename_in =~ s/_w?\d+//g;

		open (TEMPLPDB, "$tmpl_filename_in") || next;
		my @templatepdb_lines = <TEMPLPDB>;
		close TEMPLPDB;

		# get alignment
		my $filename = $THREADED_MDLS[ $j ];
		$filename =~ s/.*\///;
		$filename =~ s/\.pdb//;
		my @target_aln = split(' ',$alimap {$filename}->[3]);
		my @templ_aln = split(' ',$alimap {$filename}->[4]);
		chomp( @target_aln);
		my $target_start = 0;
		my $target_stop  = length( $target_aln[1] )-1;
		while( substr($target_aln[1],$target_start,1) eq '-' ) { $target_start++; }
		while( substr($target_aln[1],$target_stop,1) eq '-' ) { $target_stop--; }
		print STDERR "Trim $tmpl_filename to $target_start,$target_stop\n";

		# trim to common core
		# TODO assumes numbering is right ... 
		my %allreses;
		my %allseqs;
		foreach my $line (@templatepdb_lines) {
			my $chainid = substr($line,21,1);
			my $resid = int( substr($line,22,4) );
			my $resname = $three_to_one{ substr($line,17,3) };
			if (!defined $allreses{$chainid}->{$resid}) {
				$allreses{$chainid}->{$resid} = 1;
				if (!defined $allseqs{$chainid}) {
					$allseqs{$chainid} = $resname;
				} else {
					$allseqs{$chainid} = $allseqs{$chainid}.$resname;
				}
			}
		}
		my %common_reses;
		resloop: foreach my $resid( keys %{ $allreses{'A'} } ) {
			foreach my $chn ( keys %allreses ) {
				next resloop if (!exists $allreses{$chn}->{$resid});
			}
			$common_reses{$resid} = 1;
		}

		# remove unaligned termini + non-common residues
		my @templatepdb_trim1;
		my ($aliptr,$pdbptr)=(-1,-1);
		my ($ali_seq, $pdb_seq);
		my %seenreses;
		my %seenchains;
		my %allcas;
		foreach my $line (@templatepdb_lines) {
			my $chainid = substr($line,21,1);
			my ($tgt_aln, $src_aln);
			if (!defined $seenchains{$chainid}) {
				$seenchains{$chainid} = 1;
				$allcas{$chainid} = [];
				$aliptr=-1;
				$pdbptr=-1;
				%seenreses = ();

				# align 'ali' sequence with 'pdb sequence
				($ali_seq, $pdb_seq) = alignSW($templ_aln[1], $allseqs{$chainid} );
				#print "$ali_seq\n$pdb_seq\n"
			}

			my $resname = $three_to_one{ substr($line,17,3) };
			my $resid = int( substr($line,22,4) );

			# not in chain A
			next if (!defined $allreses{'A'}->{$resid} );

			if (!defined ($seenreses{$resid}) ) {
				$seenreses{$resid} = 1;
				$pdbptr++;

				# extra residues in pdb
				if (substr($ali_seq,$pdbptr,1) eq '-' && substr($pdb_seq,$pdbptr,1) ne '-' ) {
					; #do nothing
				} else {
					$aliptr++;
					# missing density in pdb
					while (substr($pdb_seq,$pdbptr,1) eq '-' && substr($ali_seq,$pdbptr,1) ne '-' ) {
						$pdbptr++;
						$aliptr++;
					}
				}
			}

			# not in some copies
			next if (!defined $common_reses{$resid});

			# termini
			next if (($aliptr < $target_start) || ($aliptr > $target_stop));
			push @templatepdb_trim1, $line;

			# remember ca trace for next step
			if (substr($line,12,4) eq " CA ") {
				push @{ $allcas{$chainid} }, [substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8)];
			}
		}

		# remove non-contacting chains
		my %contacting_chains = ('A' => 1 );
		chainloop: foreach my $chn ( keys %seenchains ) {
			next if ($chn eq 'A');
			foreach my $ca1 (@{ $allcas{'A'} }) {
			foreach my $ca2 (@{ $allcas{$chn} }) {
				if ( dist2( $ca1 , $ca2 ) < 12.0*12.0 ) {
					$contacting_chains{$chn} = 1;
					next chainloop;
				}
			}
			}
		}

		# are we now monomeric?
		next if (scalar( keys %contacting_chains ) == 1);

		open (TEMPLOUT, ">$tmpl_filename") || die "Cannot open $_";
		foreach my $line (@templatepdb_trim1) {
			next if ($line !~ /^ATOM/);
			next if (!defined $contacting_chains{substr($line,21,1)});
			my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
			$newX = vsub( $newX, $preTs->[$j]);
			$newX = mapply( $globalRs->[$j], $newX );
			$newX = vadd( $newX, $postTs->[$j]);
			substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
			substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
			substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
			print TEMPLOUT $line;
		}
		close TEMPLOUT;

		# run "autosymm" code
		my $cmd = "$SYMMDEFFILE_APP -p $tmpl_filename";
		print STDERR $cmd."\n";
		system($cmd);
	}

	# collect symmdef files
	mkdir ($SYMMDEFDIR);
	my @SYMMDEFFILES=<$SYMMDIR/*.symm>;
	
	my %symmap;
	foreach my $symmdef (@SYMMDEFFILES) {
		my $symmstem = $symmdef;
		$symmstem =~ s/.*\/([^\/]+)$/$1/;
		$symmstem =~ s/\.symm//;
		my ($templ, $tag, $symmgp) = split '_', $symmstem;
		if (!defined $symmap{$symmgp}) {
			$symmap{$symmgp} = 0;
		}
		$symmap{$symmgp} += $template_probs{$tag};
	}

	$symdeffiles{$i} = {"C1"=>[]};
	foreach my $symmdef (@SYMMDEFFILES) {
		my $symmstem = $symmdef;
		$symmstem =~ s/.*\/([^\/]+)$/$1/;
		$symmstem =~ s/\.symm//;
		my ($templ, $tag, $symmgp) = split '_', $symmstem;
		if ( $symmap{$symmgp} >= $MINSYMMPROB ) {
			move("$symmdef","$SYMMDEFDIR");
			if (!defined $symdeffiles{$i}->{$symmgp}) { $symdeffiles{$i}->{$symmgp} = []; }
			push @{ $symdeffiles{$i}->{$symmgp} }, $symmstem;
		} else {
			unlink ($symmdef);
		}
	}
}
# delete original alignments
foreach my $j (0..$#THREADED_MDLS) {
	my $tmpl_filename = $THREADED_MDLS[ $j ];
	$tmpl_filename =~ s/$PARTIALTHREADDIR/$SYMMDIR/g;
	$tmpl_filename =~ s/_w?\d+//g;
	unlink($tmpl_filename);
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
		$filename =~ s/\.pdb//;
		foreach my $alnline ( @{$alimap {$filename}} ) {
			print CLST $alnline;
		}

		my $tag = $filename;
		my $prob = $template_probs{ 0 };
		if ($tag =~ /.*_(\d\d\d).*/) {
			$tag = $1;
			$prob = $template_probs{ $tag } if defined $template_probs{ $tag };
		} elsif ($tag =~ /.*_w(\d+).*/) {
			$tag = "w$1";
			$prob = $1+$template_M_EST_counts;
		}
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
	my $mapfile = sprintf $MAPFILENAMES , $i;

	my @input_pdbs = ();
	foreach my $j (0..$#THREADED_MDLS) {
		next if ($clusterid->[$j] != $i);

		my $temppdb = $mapfile;
		$temppdb =~ s/\.mrc/.$j.pdb/g;
		open (PDBCAT, ">$temppdb") || die "Cannot open $_";

		my $tmpl_filename = $THREADED_MDLS[ $j ];
		$tmpl_filename =~ s/$PARTIALTHREADDIR/$TEMPLATEDIR/g;
		$tmpl_filename =~ s/_w?\d+//g;
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
	my @nsymm = keys %{ $symdeffiles{$i} };
	foreach my $symmgp (@nsymm) {
		my $xmlfile = sprintf $XMLFILENAMES, $i, $symmgp;
		open (XML, ">$xmlfile") || die "Cannot open $_";

		my $symmflag = 1;
		if ($symmgp eq "C1") { $symmflag=0; }
	
		print XML  "<dock_design>\n";
		print XML  "    <TASKOPERATIONS>\n";
		print XML  "    </TASKOPERATIONS>\n";
		print XML  "    <SCOREFXNS>\n";
		print XML  "        <fullatom weights=stage3_rlx symmetric=$symmflag>\n";
		print XML  "	    </fullatom>\n";
		print XML  "        <stage1 weights=stage1 symmetric=$symmflag>\n";
		print XML  "        </stage1>\n";
		print XML  "        <stage2 weights=stage2 symmetric=$symmflag>\n";
		print XML  "        </stage2>\n";
		print XML  "    </SCOREFXNS>\n";
		print XML  "    <FILTERS>\n";
		print XML  "    </FILTERS>\n";
		print XML  "    <MOVERS>\n";
		if ($symmflag == 1) {
			my (@allmovers,@allweights);
			foreach my $symmdeffile (@{ $symdeffiles{$i}->{$symmgp} }) {
				my $symmstem = $symmdeffile;
				$symmstem =~ s/\.symm//;
				print XML  "       <SetupForSymmetry name=$symmstem definition=\"$dir/$SYMMDEFDIR/$symmdeffile.symm\"/>\n";
				my ($templ, $tag, $symmgp) = split '_', $symmstem;
				push @allmovers, $symmstem;
				push @allweights, $template_probs{$tag};
			}
			print XML  "       <RandomMover name=setup_symm movers=\"".(join ',',@allmovers)."\" weights=\"".(join ',',@allweights)."\"/>\n";
		}
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
	
			my $tag = $THREADED_MDLS[$j];
			my $prob = $template_probs{ 0 };
			if ($tag =~ /.*_(\d\d\d).*/) {
				$tag = $1;
				$prob = $template_probs{ $tag } if defined $template_probs{ $tag };
				$prob += $template_M_EST;
			} elsif ($tag =~ /.*_w(\d+).*/) {
				$prob = $1+$template_M_EST_counts;
			}

			my $outline = sprintf "pdb=\"$id\" cst_file=\"$cstfile\" weight=$prob";
			print XML  "            <Template $outline/>\n";
		}
		print XML  "        </Hybridize>\n";
		print XML  "    </MOVERS>\n";
		print XML  "    <APPLY_TO_POSE>\n";
		print XML  "    </APPLY_TO_POSE>\n";
		print XML  "    <PROTOCOLS>\n";
		if ($symmflag == 1) {
			print XML  "        <Add mover=setup_symm/>\n";
		}
		print XML  "        <Add mover=hybridize/>\n";
		print XML  "    </PROTOCOLS>\n";
		print XML  "</dock_design>\n";
		close(XML);
	
		# finally write flags
		my $flagfile = sprintf $FLAGFILENAMES, $i, $symmgp;
		open (FLAGS, ">$flagfile") || die "Cannot open $_";
		print FLAGS "-in:file:fasta $dir/$fastafile\n";
		print FLAGS "-out:prefix $dirtag\n";
		print FLAGS "-out:file:silent $dir/../$dirtag.$symmgp.cl$i.silent\n";
		print FLAGS "-in:file:native $dir/$NATIVEDIR/$dirtag.pdb\n";
		print FLAGS "-parser:protocol $dir/$xmlfile\n";
		printf FLAGS "-edensity::mapfile $dir/$MAPFILENAMES\n", $i;
	}
}

exit 0;

###########
###########

sub dist {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}

sub dist2 {
	my ($x, $y) = @_;
	my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
	return ( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}


## transitive alignment realigns hhr sequence to pdb sequence
sub realign {
	my ($templ_aln, $template, $target) = @_;

	my $ptr1=0;
	my $ptr2=0;

	# find gaps in templ_aln
	while (	$ptr1 <= length($templ_aln) && $ptr2 <= length($template) ) {
		if ( substr($template, $ptr2, 1) eq substr($templ_aln, $ptr1, 1) ) {
			$ptr1++;
			$ptr2++;
		} elsif ( substr($templ_aln, $ptr1, 1) eq '-') {
			$target = substr($target, 0, $ptr1).'-'.substr($target, $ptr1);
			$ptr1++;
		} else {
			return ""; # failure
		}
	}
	return $target;
}

## DP sequence alignment
sub alignSW {
	my ($seq1, $seq2) = @_;
	
	# scoring scheme
	my $MATCH    =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP      = -1; # -1 for any gap
	
	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "n";
	for(my $j = 1; $j <= length($seq1); $j++) {
		$matrix[0][$j]{score}   = -$j;
		$matrix[0][$j]{pointer} = "l";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = -$i;
		$matrix[$i][0]{pointer} = "n";
	}

	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = -99999;
	
	for(my $i = 1; $i <= length($seq2); $i++) {
		for(my $j = 1; $j <= length($seq1); $j++) {
			my ($diagonal_score, $left_score, $up_score);
			
			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);       
			if ($letter1 eq $letter2) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			}
			else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}
			
			# calculate gap scores
			if ($letter2 eq '-') {
				$up_score   = $matrix[$i-1][$j]{score} + $MATCH;
			} else {
				$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			}
			if ($letter1 eq '-') {
				$left_score = $matrix[$i][$j-1]{score} + $MATCH;
			} else {
				$left_score = $matrix[$i][$j-1]{score} + $GAP;
			}
			
			#if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
			#	$matrix[$i][$j]{score}   = 0;
			#	$matrix[$i][$j]{pointer} = "n";
			#	next; # terminate this iteration of the loop
			#}
			
			# choose best score
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "d";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "l";
				}
			} else {
				if ($up_score >= $left_score) {
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "u";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "l";
				}
			}
			
			# set maximum score
			#if ($matrix[$i][$j]{score} > $max_score) {
			#	$max_i     = $i;
			#	$max_j     = $j;
			#	$max_score = $matrix[$i][$j]{score};
			#}
		}
	}
	
	# trace-back
	my $align1 = "";
	my $align2 = "";
	

	# find max score
	for(my $i = 1; $i <= length($seq2); $i++) {
		if ($matrix[$i][length($seq1)]{score} > $max_score) {
			$max_i     = $i;
			$max_score = $matrix[$i][length($seq1)]{score};
		}
	}
	my $i = $max_i;
	my $j = length($seq1);

	while (1) {
		last if $matrix[$i][$j]{pointer} eq "n";
		
		if ($matrix[$i][$j]{pointer} eq "d") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			$i--; $j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "l") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= "-";
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "u") {
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
		}   
	}
	
	$align1 = reverse $align1;
	$align2 = reverse $align2;

	return ($align1, $align2, $i, $max_score);
}