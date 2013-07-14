#!/usr/bin/perl
##
##
###############################################################################

use strict;
use Math::Trig;   # inv trig ops
use POSIX qw(ceil floor fmod fabs);
use constant PI    => 4 * atan2(1, 1);

###############################################################################

if ($#ARGV < 1) {
        print STDERR "usage: $0 template.pdb template.ali\n";
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

my $pdbfile = shift @ARGV;
my $alifile = shift @ARGV;

## read ali file
open (ALI, $alifile) || die "Cannot open $pdbfile.";
# ignore 1st three lines
<ALI>;
<ALI>;
<ALI>;
my $line = <ALI>;
chomp $line;
my @targetseq = split ' ',$line;
my $line = <ALI>;
chomp $line;
my @templateseq = split ' ',$line;
close (ALI);

## sanity check
if ( length($targetseq[1]) != length($templateseq[1]) ) { die "sequence length mismatch"; }

## read pdb file
my @respointers;
my $pdbseq = "";
my @filebuf;

my $lastresid=-999;
my $thisresid=-999;

open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
while (<PDB>) {
        chomp;
        if (/^ATOM/ || /^HETATM/) {
                my $atom = substr ($_, 12, 4);
                my $confid = substr ($_, 16, 1);

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

                        #if ($lastresid != -999) {
                        #       foreach my $resid ($lastresid..$thisresid-1) {
                        #               $pdbseq = $pdbseq.'-';
                        #               push @respointers, $#filebuf;
                        #       }
                        #}
                        push @respointers, $#filebuf;
                        $pdbseq = $pdbseq.$three_to_one{ $restype };
                }
        }
}
push @respointers, $#filebuf+1;
close (PDB);

##
my $nviolations = 0;
my $MAXVIOLATIONS = 1;
my $pdbptr = $templateseq[0];
foreach my $i (0..length($targetseq[1])-1) {
        if (substr($templateseq[1], $i, 1) eq "-" || substr($pdbseq, $i, 1) eq "-" ) {
                next;
        }

        if ( $pdbptr > length($pdbseq) ||  substr($pdbseq, $pdbptr, 1) ne substr($templateseq[1], $i, 1) ) {
                print STDERR "warning seq mismatch ".$i." ".substr($pdbseq, $pdbptr, 1)." != ".substr($templateseq[1], $i, 1)."\n";
                $nviolations = $nviolations + 1;
        }
        if ( $nviolations > $MAXVIOLATIONS ) { die "sequence mismatch"; }

        if ( $pdbptr > length($pdbseq) ) { next; }

        if ( substr($templateseq[1], $i, 1) eq substr($targetseq[1], $i, 1) ) {
                # take all heavyatoms
                foreach my $j ($respointers[$pdbptr]..$respointers[$pdbptr+1]-1) {
                        my $line_j = $filebuf[$j];

                        if (defined $three_to_one{ substr($line_j, 17, 3) } && substr($line_j, 13, 1) ne "H") {
                                print $line_j."\n";
                        }
                }
        } elsif ( substr($targetseq[1], $i, 1) ne "-") {
                # take everything to CG
                foreach my $j ($respointers[$pdbptr]..$respointers[$pdbptr+1]-1) {
                        my $line_j = $filebuf[$j];
                        my $atmid = substr($line_j, 12, 4);
                        if (defined $three_to_one{ substr($line_j, 17, 3) } &&
                            ($atmid eq " CB " ||  $atmid eq " CA " || $atmid eq " N  " || $atmid eq " O  " 
                              || $atmid eq " C  " || $atmid eq " CB " || $atmid eq " CG " || $atmid eq " OG ") ) {
                                print $line_j."\n";
                        }
                }

        }

        $pdbptr = $pdbptr+1;
}
