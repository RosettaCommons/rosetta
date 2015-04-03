#!/usr/bin/perl
##
##
###############################################################################

use strict;
use Math::Trig;   # inv trig ops
use POSIX qw(ceil floor fmod fabs);
use constant PI    => 4 * atan2(1, 1);

use File::Copy;
use Getopt::Long;
use Net::FTP;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

###############################################################################

Getopt::Long::Configure ('bundling');
my $local;
GetOptions ('l|local=s' => \$local,);

if ($#ARGV < 0) {
	print STDERR "usage: $0 align.hhr\n";
	exit -1;
}

my %three_to_one = ( 
	'GLY' => 'G', 'ALA' => 'A', 'VAL' => 'V', 'LEU' => 'L', 'ILE' => 'I',
	'PRO' => 'P', 'CYS' => 'C', 'MET' => 'M', 'HIS' => 'H', 'PHE' => 'F',
	'TYR' => 'Y', 'TRP' => 'W', 'ASN' => 'N', 'GLN' => 'Q', 'SER' => 'S',
	'THR' => 'T', 'LYS' => 'K', 'ARG' => 'R', 'ASP' => 'D', 'GLU' => 'E',

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

my $hhrfile = shift @ARGV;

## read ali file
open (HHR, $hhrfile) || die "Cannot open $hhrfile.";

my @altpdbs;
my @pdbids;
my @chainids;
my @templates;
my @targets;
my @templatestarts;
my @targetstarts;

##############
## PARSE HHR FILE
##############
while ( my $line = <HHR> ) {
	if (defined $line && $line =~ /^>(\w\w\w\w)_(\w)/) {
		my $pdbid = $1;
		my $chainid = $2;
		my $templseq = '';
		my $targetseq = '';
		my $templstart = -1;
		my $targetstart = -1;
		my $altpdbs="";

		# alt PDBs
		if ($line =~ /PDB: (.*)$/) {
			$altpdbs = $1;
		}

		do {
			my $extraline = <HHR>;
			if ($extraline =~ /^Confidence/) {
				<HHR>;
			}
			<HHR>;
			$line = <HHR>;
			if (defined $line && $line !~/^No/  && $line !~/^Done/) {
				my $targetline = <HHR>;
				if ($targetline =~ /^Q.* (\d+) (\S+) /) {
					if ($targetstart < 0) { $targetstart=$1; }
					$targetseq = $targetseq.$2;
				} else {
					die "(1) Error parsing target line at $targetline\n";
				}
				<HHR>;
				<HHR>;
				<HHR>;
				my $templline = <HHR>;
				if ($templline =~ /^T.* (\d+) (\S+) /) {
					if ($templstart < 0) { $templstart=$1; }
					$templseq = $templseq.$2;
				} else {
					die "(2) Error parsing template line at $templline\n";
				}
				<HHR>;
				<HHR>;
			}
		} while (defined $line && $line !~/^No/ && $line !~/^Done/);

		## sanity check
		if ( length($targetseq) != length($templseq) ) { die "Sequence length mismatch."; }

		# done!  now push onto the arrays
		push @altpdbs, $altpdbs;
		push @pdbids, $pdbid;
		push @chainids, $chainid;
		push @templates, $templseq;
		push @targets, $targetseq;
		push @templatestarts, $templstart;
		push @targetstarts, $targetstart;
	}
}

my $nseqs = $#pdbids+1;
print "Read $nseqs alignments from $hhrfile.\n";

my %seqs_to_skip;

if (defined $local) {
	foreach my $i (0..$nseqs-1) {
		my $pdbid = $pdbids[$i];
		my $folder_id = substr($pdbid,1,2);
		my $copied = copy("$local/$folder_id/$pdbid.pdb","./$pdbid.pdb");

		if (!$copied) {
			# try alts
			my @altpdb_i=split ' ', $altpdbs[$i];
			ALTS: foreach my $altpdb (@altpdb_i) {
				if ($altpdb =~ /(\w\w\w\w)_(\w)/) {
					$copied = copy("$local/$folder_id/$1.pdb","./$1.pdb");
					if ($copied) {
						$pdbids[$i] = $1;
						$chainids[$i] = $2;
						print STDERR "WARNING $pdbid replaced with $1\n";
						last ALTS;
					}
				}
			}
		}

		if (!$copied) {
			print STDERR "ERROR copying $pdbid! skipping...\n";
			$seqs_to_skip{$i} = 1;
		}
	}
} else {
	################
	### connect to PDB ftp
	################
	print "Connecting to ftp.wwpdb.org\n";
	my $ftp = Net::FTP->new("ftp.wwpdb.org", (Passive => 1, Debug => 0))
		  or die "Cannot connect to ftp.wwpdb.org: $@";
	$ftp->login("anonymous",'-anonymous@')
		  or die "Cannot login ", $ftp->message;
	$ftp->binary;

	################
	### fetch PDBS!!
	################
	foreach my $i (0..$nseqs-1) {
		my $pdbid = $pdbids[$i];
		my $folder_id = substr($pdbid,1,2);
	
		#  e.g. ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/00/pdb100d.ent.gz
		#print "cwd /pub/pdb/data/structures/divided/pdb/$folder_id/\n";
		$ftp->cwd("/pub/pdb/data/structures/divided/pdb/$folder_id/") or die "Cannot change working directory ", $ftp->message;
		print "Fetch pdb$pdbid.ent.gz\n";

		my $copied = $ftp->get("pdb$pdbid.ent.gz");
		if (!$copied) {
			# try alts
			my @altpdb_i=split / /, $altpdbs[$i];
			ALTS: foreach my $altpdb (@altpdb_i) {
				if ($altpdb =~ /(\w\w\w\w)_(\w)/) {
					print " --> try pdb$1.ent.gz\n";
					$pdbid = $pdbids[$i];
					my $new_folder_id = substr($1,1,2);
					$ftp->cwd("/pub/pdb/data/structures/divided/pdb/$new_folder_id/") or die "Cannot change working directory ", $ftp->message;
					$copied = $ftp->get("pdb$1.ent.gz");
					if ($copied) {
						$pdbids[$i] = $1;
						$chainids[$i] = $2;
						print STDERR "WARNING $pdbid replaced with $1\n";
						$pdbid = $1;
						last ALTS;
					}
				}
			}
		}
		if (!$copied) {
			print STDERR "ERROR downloading $pdbid! skipping...\n";
			$seqs_to_skip{$i} = 1;
		} else {
			gunzip "pdb$pdbid.ent.gz" => "$pdbid.pdb" or die "gunzip failed: $GunzipError\n";
			unlink( "pdb$pdbid.ent.gz" );
		}
	}

	## disconnect
	$ftp->quit;
}

open (ALIOUT, ">alignments.filt") || die "Cannot open $hhrfile.";
foreach my $i (0..$nseqs-1) {
	next if (defined $seqs_to_skip{$i});

	my $pdbid = $pdbids[$i];
	my $chainid = $chainids[$i];
	my $template = $templates[$i];
	my $target = $targets[$i];
	my $templatestart = $templatestarts[$i];
	my $targetstart = $targetstarts[$i];

	################
	### parse PDB
	################
	my ($pdbseq, $filebuf, $respointers) = parsePDB( "$pdbid.pdb", $chainid );

	################
	### align HHR seq to PDB seq
	################
	my ($templ_aln, $pdb_aln, $pdbstart) = alignSW( $template, $pdbseq );

	## this will happen if either:
	##   (a) the PDB and hhr file's sequence differ dramatically
	##   (b) the PDB contains residues not in the hhr file
	## (b) is okay; (a) is a problem so we should abort
	if ($templ_aln ne $template) {
		# attempt to remap target seq
		($target) = realign( $templ_aln, $template, $target );

		if (length($target) == 0) {
			print "WARNING: possible problem aligning PDB file $pdbid.  Verify sequence.\n";
			print "   ".$templ_aln."\n";
			print "   ".$template."\n";
			next;
		}
	}

	################
	### strip PDB
	################
	open (PDBOUT, ">$pdbid"."$chainid.pdb") || die "Cannot open $hhrfile.";
	my $pdbptr = $pdbstart;
	foreach my $i (0..length($target)-1) {
		# gap
		#if (substr($template, $i, 1) eq "-" || substr($pdb_aln, $i, 1) eq "-" ) { next; }
		# exact match
		if (substr($pdb_aln, $i, 1) eq substr($target, $i, 1)) {
			# take all heavyatoms
			foreach my $j ($respointers->[$pdbptr]..$respointers->[$pdbptr+1]-1) {
				my $line_j = $filebuf->[$j];

				if (defined $three_to_one{ substr($line_j, 17, 3) } && substr($line_j, 13, 1) ne "H" && substr($line_j, 12, 1) ne "H") {
					print PDBOUT $line_j."\n";
				}
			}
			$pdbptr++;
		}
		elsif (substr($pdb_aln, $i, 1) ne "-") {
			foreach my $j ($respointers->[$pdbptr]..$respointers->[$pdbptr+1]-1) {
				my $line_j = $filebuf->[$j];

				if (defined $three_to_one{ substr($line_j, 17, 3) } && substr($line_j, 13, 1) ne "H" && substr($line_j, 12, 1) ne "H") {
					print PDBOUT $line_j."\n";
				}
			}
			$pdbptr++;
		}

		# mismatch
		#elsif ( substr($target, $i, 1) ne "-") {
		#	# take everything to CG
		#	foreach my $j ($respointers->[$pdbptr]..$respointers->[$pdbptr+1]-1) {
		#		my $line_j = $filebuf->[$j];
		#		my $atmid = substr($line_j, 12, 4);
		#		if (defined $three_to_one{ substr($line_j, 17, 3) } &&
		#			($atmid eq " CB " ||  $atmid eq " CA " || $atmid eq " N  " || $atmid eq " O  " 
		#			  || $atmid eq " C  " || $atmid eq " CB " || $atmid eq " CG " || $atmid eq " OG ") ) {
		#			print PDBOUT $line_j."\n";
		#		}
		#	}
		#}
		#if (substr($pdb_aln, $i, 1) ne "-") {
		#	$pdbptr++;
		#}
	}
	close (PDBOUT);

	################
	### write Rosetta ali file
	################
	#  ## 1CRB_ 2qo4b
	#  # hhsearch
	#  scores_from_program: 0 1.00
	print ALIOUT "## 1XXX_ $pdbid"."$chainid"."_".(201+$i)."\n";
	print ALIOUT "# hhsearch\n";
	print ALIOUT "scores_from_program: 0 1.00\n";
	print ALIOUT ($targetstart-1)." $target\n";
	print ALIOUT "0 $pdb_aln\n";
	print ALIOUT "--\n";
}
close (ALIOUT);


#############################################################
#############################################################
## read pdb file
sub parsePDB {
	my $pdbfile = $_[0];
	my $tgt_chain = $_[1];

	my $lastresid=-999;
	my $thisresid=-999;
	my @respointers;
	my $pdbseq = "";
	my @filebuf;

	open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
	while (<PDB>) {
		chomp;
		if (/^ATOM/ || /^HETATM/) {
			my $atom = substr ($_, 12, 4);
			my $chnid = substr ($_, 21, 1);
			my $confid = substr ($_, 16, 1);
	
			# select chain
			next if ($chnid ne $tgt_chain);
	
			# remove alt confs
			next if ($confid ne " " && $confid ne "A");
	
			# map MSE->MET
			my $restype = substr ($_, 17, 3);
			if ($restype eq "MSE") {
				substr ($_, 0, 6) = "ATOM  ";
				if ($atom eq "SE  ") {
					substr ($_, 12, 4) = " SD ";
				}
				substr ($_, 17, 3) = "MET";
			}
	
			push @filebuf, $_;
	
			if ($atom eq " N  ") {
				$lastresid = $thisresid;
				$thisresid = substr ($_, 22, 4);
	
				push @respointers, $#filebuf;
				$pdbseq = $pdbseq.$three_to_one{ $restype };
			}
		}
	}
	push @respointers, $#filebuf+1;
	close (PDB);

	return ($pdbseq, \@filebuf, \@respointers);
}


##############################################################
##############################################################
## transitive alignment realigns hhr sequence to pdb sequence
##   ($target, $is_good) = realign( $templ_aln, $template, $target );
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


##############################################################
##############################################################
## align hhr sequence to pdb sequence
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
		$matrix[0][$j]{score}   = 0;
		$matrix[0][$j]{pointer} = "l";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = 0;
		$matrix[$i][0]{pointer} = "n";
	}
	
	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = 0;
	
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
			if ($matrix[$i][$j]{score} > $max_score) {
				$max_i     = $i;
				$max_j     = $j;
				$max_score = $matrix[$i][$j]{score};
			}
		}
	}
	
	# trace-back
	my $align1 = "";
	my $align2 = "";
	

	# find max score
	for(my $i = 1; $i <= length($seq2); $i++) {
		if ($matrix[$i][length($seq2)]{score} > $max_score) {
			$max_i     = $i;
			$max_score = $matrix[$i][length($seq2)]{score};
		}
	}
	#my $j = $max_j;
	my $i = $max_i;
	#my $i = length($seq2);
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

	return ($align1, $align2, $i);
}
