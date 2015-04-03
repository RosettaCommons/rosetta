#! /usr/bin/perl -w
# alignblast.pl version 1.1.4 (April 2005)
# Extract a multiple alignment of hits from Blast or PsiBlast output
# Usage:   alignblast.pl [options] blast.out alignment-file
#
# Please report bugs to johannes@soeding.com. Thank you.


use strict;

# Default options
my $E_max=1E-3;
my $cov_thrshd=0;
my $qid_thrshd=0;
my $sc_thrshd=-10;
my $outformat="a2m";
my $append=0;
my $infile;
my $outfile;
my $v=2;

sub fexit()
{
    print("\nExtract a multiple alignment of hits (HSPs) from Blast or PSI-BLAST output\n");
    print("Usage: alignblast.pl blast-file alignment-file [options]\n");
    print("Options for thresholds\n");
    print("  -e   e-value  : maximum e-value (default=$E_max)\n");
    print("  -qid percent  : minimum sequence identity to query in % (default=$qid_thrshd) \n");
    print("                  (seq-id = identities in match columns / HSP residues in match columns)\n");
    print("  -cov coverage : minimum coverage in % (default=$cov_thrshd) \n");
    print("  -s/c value    : minimum score per alignment column in 1/3 bits (default=no cut-off) \n");
    print("Options for output format:\n");
    print("  -fas          : aligned fasta; all residues upper case, all gaps '-'\n");
    print("  -a2m          : like FASTA, but residues aligned to a query residue in upper case,\n");
    print("                  inserts in lower case, deletes '-', gaps aligned to inserts '.' (default)\n");
    print("  -a3m          : like -a2m, but gaps aligned to inserts omitted\n");
    print("  -psi          : format for -B option of PSI-BLAST; inserts relative to query (=first)\n");
    print("                  sequence omitted, capitalization of residues same as query sequence\n");
    print("  -ufas         : unaligned fasta format (without gaps)\n");
    print("Other options:\n");
    print("  -Q   file     : insert A2m/FASTA-formatted query sequence into output alignment;\n");
    print("                  (all query residues will be match states in A2M format)\n");
    print("  -v            : verbose mode (default=off)\n");
    print("  -append       : append output to file (default=overwrite)\n");
    print("  -best         : extract only the best HSP per sequence (default=off)\n");
    print("\n");
    print("Examples: \n");
    print("alignblast.pl 1enh.blast 1enh.fas\n");
    print("alignblast.pl 1enh.blast 1enh.a3m -e 1E-4 -cov 50 -s/c 1 -a3m\n");
    print("\n");
    exit;
}

# Activate autoflushing on STDOUT
$|= 1; 

###################################################################################################
# Process input options

my $ARGC=scalar(@ARGV);
if ($ARGC<2) {&fexit;}

# Variable declarations
my $i;                 # residue index
my $j;                 # residue index
my $k;                 # sequence index
my $options="";
my $line;              # line read in from file
my $query_length=0;    # number of residues in query sequence
my $query_match=0;     # number of upper-case residues (=match states) in query sequence
my $nameline;          # >template_name
my $Evalue;            # e-value of hit
my $score;             # bit score of hit
my $coverage;          # hit-length/$query_length
my $hit_length;        # number of residues in HSP
my $score_col;         # score per column

my $first_res;         # index of first query residue in pairwise alignment
my $last_res;          # index of last  query residue in pairwise alignment
my @query_res;         # query residues from current pairwise alignment 
my @template_res;        # template residues from current pairwise alignment 
my $query_res;         # query residues from current pairwise alignment 
my $template_res;        # template residues from current pairwise alignment 
my $line_number=0;
my $new_hit="";        # new sequence record; is only written if coverage threshold is exceeded
my $nhit=0;            # counts the number of sequences already in alignment
my @hitnames;          # $hitnames[$nhit] is the nameline of the ihit'th hit
my @hitseqs;           # $hitseqs[$nhit] contains the residues of the ihit'th hit
my @match;             # for -q option: $match[$i]=1 if $i'th query residue is capital letter in query, else 0
my $qid;               # $qid is sequence identity with query (for -q option: CONSIDER ONLY MATCH STATES)
my $len;               # $len is sequence number of residues of seq k aligned with a match state in query
my $skip=0;            # skip this template sequence because it might be a synthetic fusion protein
my $best=0;            # extract only the best HSP per sequence
my $query_file="";     # name of A2M/FASTA-formatted file with query sequence
my $query_name;        # name of query sequence
my $queryseq;          # residues of query read in  with -Q option
my @queryseq;


for ($i=0; $i<$ARGC; $i++) {$options.=" $ARGV[$i]";}

# Set E-value thresholds etc. for inclusion in alignment
if ($options=~s/ -e\s+(\S+)//g)    {$E_max=$1;}
if ($options=~s/ -qid\s+(\S+)//g) {$qid_thrshd=$1;}
if ($options=~s/ -cov\s+(\S+)//g) {$cov_thrshd=$1;}
if ($options=~s/ -s\/c\s+(\S+)//g) {$sc_thrshd=$1;}
if ($options=~s/ -best//g)   {$best=1;} 

# Set output format
if ($options=~s/ -psi//g) {$outformat="psi";}
if ($options=~s/ -fas//g) {$outformat="fas";}
if ($options=~s/ -a2m//g) {$outformat="a2m";}
if ($options=~s/ -a3m//g) {$outformat="a3m";}
if ($options=~s/ -ufas//g) {$outformat="ufas";}

# Set input and output file
if ($options=~s/ -i\s+(\S+)//g) {$infile=$1;}
if ($options=~s/ -o\s+(\S+)//g) {$outfile=$1;}
if ($options=~s/ -app//g) {$append=1;}
if ($options=~s/ -Q\s+(\S+)//g)  {$query_file=$1;}

# Verbose mode?
if ($options=~s/ -v\s*(\d)//g) {$v=$1;}
if ($options=~s/ -v//g) {$v=2;}

# Read infile and outfile 
if (!$infile  && $options=~s/^\s*([^- ]\S+)\s*//) {$infile=$1;} 
if (!$outfile && $options=~s/^\s*([^- ]\S+)\s*//) {$outfile=$1;} 

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile || !$outfile) {&fexit();}

# Scan Blast output file for query length (needed for coverage)
open(INFILE,"<$infile") or die ("Error: cannot open $infile: $!\n");
$line_number++;
while ($line=<INFILE>)
{
    if ($line =~ /^\s*\((\d+)\s+letters\)/) {$query_length = $1; last;}
    $line_number++;
}

# Verbose feedback?
if ($v>=3) 
{
    print("\n");
    print("E-value maximum        : $E_max\n");
    print("coverage threshold     : $cov_thrshd\n");
    print("sequence id threshold  : $qid_thrshd\n");
    print("score/column threshold : $sc_thrshd\n");
    print("Blast file             : $infile\n");
    print("Output file            : $outfile\n");
    print("Output format          : $outformat\n");
    print("Query length = $query_length\n");
}
		    
#Include query sequence as first sequence in alignment?
if ($query_file) {
    open(QUERYFILE,"<$query_file") or die ("ERROR: Cannot open $query_file: $!\n");
    while($line=<QUERYFILE>) # Read name line
    {
	if ($line=~/^>(.*)/) 
	{
	    $query_name=$1;
	    last;
	}
    }
    $hitseqs[0]="";
    while($line=<QUERYFILE>)  # Read residues
    {
	if ($line=~/^>/) {last;}
	chomp($line);
	$line=~s/\s+//g;      # remove white space
	$hitseqs[0].=$line;
    }
    close(QUERYFILE);

    # Prepare name line of hit
    if ($outformat eq "psi") {
	$query_name=~/^(\S{1,20})\S*\s*(.*)/;       # delete everything after first block
	$line=sprintf("%s",$1);
	$line=~ tr/ /_/;
	$hitnames[0] = sprintf("%-31.31s ",$line);
    } else {
	$hitnames[0] = sprintf(">%s",$query_name); 
    } 
    $hitseqs[0] =~ tr/-.//d;      # delete all gaps from query
    $queryseq = $hitseqs[0];
    $hitseqs[0] =~ tr/a-z/A-Z/d;  # capitalize hitseq[0] and delete gaps
    $hitseqs[0] =~ tr/Uu/Cc/;
    $nhit=1;

    # Capitalize query
    $queryseq =~ tr/a-z/A-Z/;
    $query_match = ($queryseq=~tr/A-Z/A-Z/);  # count number of match states in query
    
    # Determine match columns as those with upper case residue in query 
    @queryseq=unpack("C*",$queryseq);
    for ($j=0; $j<@queryseq; $j++) {
	if ($queryseq[$j]>=65 && $queryseq[$j]<=90) {$match[$j]=1;} else {$match[$j]=0;}
    }
}


# Search for "Results from round"
# If found, we are looking at PsiBlast output and have to search for the beginning of the last round
$line_number++;
for ($i=1; $i<=9; $i++)
{
    $line=<INFILE>;
    if (!$line) {die("Error in alignblast.pl: blast output file $infile truncated: $!\n");} 
    $line_number++;
    if ($line=~/^Results from round/) 
    {
	# PsiBlast output! Search for line number with last occurence of "Results from round"
	if ($v>=3) {print("PsiBlast output with several rounds detected. Searching for last round...\n");}
	my $last_line=$line_number;
	while ($line=<INFILE>) #scan through PsiBlast-output line by line
	{
	    if ($line=~/^Results from round/) {$last_line=$line_number;}
	    $line_number++;
	}
	# Advance to line with last occurence of "Results from round" 
	close INFILE;
	open INFILE,"<$infile" or die ("cannot open $infile: $!\n");
	for ($j=1; $j<=$last_line; $j++) {<INFILE>;}
	last;
    } 
}


########################################################################################
# Read template namelines and HSPs (Evalue, score etc. and pairwise alignment)

while ($line = <INFILE>) #scan through PsiBlast-output line by line
{
    # New nameline found?
    if ($line=~s/^>//) {
	$line=~s/\s+/ /g;
	chomp($line);
	$nameline=$line;
	while ($line=<INFILE>) {
	    if ($line=~/^\s+Length =/) {last;}
	    chomp($line);
	    $nameline.=$line;
	}
	$nameline=~s/\s+/ /g;
	$nameline=~s/\s+gi\|/   gi\|/g;

	# Is sequence a synthetic fusion protein ? 
	# Remark: excluding fusion proteins from PSI-BLAST profiles may prevent a corruption of the evolving profile
#	if ($nameline=~/(\[synthetic| synthetic|construct|cloning|vector|chimeric|fusion)/i) {$skip=1;} else {$skip=0;}
    }

    # New HSP found? 
    elsif (!$skip && $line=~/^ Score =/) 
    {
	if($best) {$skip=1;} # skip all following hits with same sequence?
	
	# First check whether E-value is small enough
	if($line =~ /^ Score =\s*(\S+)\s*bits\s*\S*\s*Expect =\s*(\S+)/) {
	    $score=$1;
	    $Evalue=$2;
	} else { 
	    print("\nWARNING: wrong format in blast output, line $.. Expecting Score = ... Expect = ..\n$line\n");
	}
	$Evalue=~s/^(e|E)/1$1/;   # Expect = e-123 -> 1e-123
	$Evalue=~tr/,//d; 
	if ($Evalue>$E_max) {$new_hit=""; next;} # reject hit

	# Record sequence identity 
	# (not needed, qid calculated afterwards WITHOUT counting template residues aligned to gaps in query)
	$line=<INFILE>;
	if ($line =~ /^ Identities =\s*\S+\/(\S+)\s+\((\S+)%\)/) {$qid=$2; $line=<INFILE>;}
	
      	# Skip another line and read following line
	$line=<INFILE>;
	
	# Read pairwise alignment 
	$query_res="";
	$template_res="";
	if ($line!~/^Query:\s*(\d+)\s+\S+/) {die("Error 2: wrong format of blast output in $infile, line $.\n");} 
	$line=~/^Query:\s*(\d+)\s+\S+/;
	$first_res=$1;
	while ($line=~/Query:\s*\d+\s+(\S+)\s+(\d+)/) # Cycle in this loop until no new "Query:" lines are found
	{
	    $query_res .= $1;
	    $last_res=$2;
	    $line=<INFILE>;
	    $line=<INFILE>;
	    if ($line!~/^Sbjct:\s*\d+\s+(\S+)/) {die("Error 3: wrong format of blast output in $infile, line $.\n");} 
	    $template_res .= $1;
	    $line=<INFILE>;
	    $line=<INFILE>;
	} # end while(1)
	
	# Check lengths
	if (length($template_res)!=length($query_res)) {
	    print("ERROR: Query and template lines do not have the same length!\n");
	    print("Q: $query_res\n");
	    print("T: $template_res\n");
	    exit 1;
	}
	    
	# Check whether hit has sufficient score per column
	$score_col=$score/length($query_res);
	if ($score_col<$sc_thrshd) {return 0;}            # reject?

	@query_res =unpack("C*",$query_res);
	@template_res=unpack("C*",$template_res);
	
	# Check whether hit has sufficient sequence identity and coverage with query
	if (!$query_file) {
	    $len=0; $qid=0;
	    for ($i=0; $i<scalar(@query_res); $i++) {
		if ($template_res[$i]!=45 && $query_res[$i]!=45) {  # count only non-gap template residues in match columns!
		    $len++;
		    if ($query_res[$i]==$template_res[$i]) {$qid++;}
		}
	    }
	    $coverage = 100*$len/$query_length;
	} else {
	    $len=1; $qid=0; $j=$first_res-1; # if first_res=1 then $j=0
	    for ($i=0; $i<scalar(@query_res); $i++) { 
		if ($query_res[$i]!=45) {
		    if ($template_res[$i]!=45 && $match[$j]) {      # count only non-gap template residues in match columns!
			$len++;
			if ($query_res[$i]==$template_res[$i]) {$qid++;}
		    }
		    $j++;                                         # $j = next position in query
		}
	    }
	    $coverage = 100*$len/$query_match;
	} 
	if ($len==0) {next;}                              # Reject hit? 
	if (100*$qid/$len<$qid_thrshd) {next;}            # Reject hit?
	if ($coverage<$cov_thrshd) {next;}                # Reject hit?
#	print("Q: $query_res\n");
#	print("T: $template_res\n\n");

	if ($v>=3) {printf("nhit=%-2i  qid=%-3i  qlen=%-3i  qid=%-3i%% s/c=%-6.3f\n",$nhit,$qid,$len,100*$qid/$len,$score_col);}
	
	# Record residues  
	$new_hit = "-"x($first_res-1);     	          # Print gaps at beginning of sequence
	if ($outformat eq "psi") {
	    for ($i=0; $i<scalar(@query_res); $i++) {
		if ($query_res[$i]!=45) {                 # residues aligned to gaps are ignored
		    $new_hit .= uc(chr($template_res[$i])); # UPPER case if aligned with a query residue (match state)
		}
	    }
	} else {
	    for ($i=0; $i<scalar(@query_res); $i++) {
		if ($query_res[$i]!=45) {
		    $new_hit .= uc(chr($template_res[$i])); # UPPER case if aligned with a query residue (match state)
		} else {
		    if($template_res[$i]!=45) {
			$new_hit.=lc(chr($template_res[$i])); # lower case if aligned with a gap in the query (insert state)
		    }
		}
	    }
	}
	$new_hit .= "-" x ($query_length-$last_res);      # Print gaps at end of sequence
	$new_hit =~ tr/Uu/Cc/;
	$hitseqs[$nhit] = $new_hit; 
#	printf("%s\n",$new_hit);
	
	# Prepare name line of hit
	if ($outformat eq "psi") {
	    $nameline=~/^(\S{1,20})\S*\s*(.*)/;           # delete everything after first block
	    $line=sprintf("%s:(%i-%i)",$1,$first_res,$last_res);
	    $line=~ tr/ /_/;
	    $hitnames[$nhit] = sprintf("%-31.31s ",$line);
	} else {
	    $nameline=~/^(\S*)\s*(.*)/;                   # delete everything after first block
	    $hitnames[$nhit] = sprintf(">%s:(%i-%i) %s  E=%g s/c=%4.2f id=%.0f%% cov=%.0f%%",
				       $1,$first_res,$last_res,$2,$Evalue,$score_col,100*$qid/$len,$coverage); 
	} 
	
	$nhit++;
    } # end elseif new HSP found

} # end while ($line)
########################################################################################
close INFILE;


# If output format is fasta or a2m we have to insert gaps:
if ($outformat ne "psi")
{
    my @len_ins; # $len_ins[$j] will count the maximum number of inserted residues after match state $j.
    my @inserts; # $inserts[$j] contains the insert (in small case) of sequence $k after the $j'th match state
    my $insert;
    my $ngap;

    # For each match state determine length of LONGEST insert after this match state and store in @len_ins
    for ($k=0; $k<$nhit; $k++) {
	# split into list of single match states and variable-length inserts
	# ([A-Z]|-) is the split pattern. The parenthesis indicate that split patterns are to be included as list elements
	# The '#' symbol is prepended to get rid of a perl bug in split
	$j=0;
 	@inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
#	printf("%3i: %12.12s %s\n",$k,$hitnames[$k],$hitseqs[$k]);
#	printf("Sequence $k: @inserts\n");
	foreach $insert (@inserts) {
	    if( !defined $len_ins[$j] || length($insert)>$len_ins[$j]) {
		$len_ins[$j]=length($insert);
	    }
	    $j++;
#	    printf("$insert|");
	}
#	for (my $i=0; $i<@inserts; $i++) {printf("%s%-2i ",$inserts[$i],$len_ins[$i]);}
#	printf("\n");
    }

    # After each match state insert residues and fill up with gaps to $len_ins[$i] characters
    for ($k=0; $k<$nhit; $k++) {
	# split into list of single match states and variable-length inserts
	@inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
	$j=0;
	
	# append the missing number of gaps after each match state
	foreach $insert (@inserts) {
	    if($outformat eq "fas") {
		for (my $l=length($insert); $l<$len_ins[$j]; $l++) {$insert.="-";}
	    }
	    else {
		for (my $l=length($insert); $l<$len_ins[$j]; $l++) {$insert.=".";}
	    }
	    $j++;
	}
	$hitseqs[$k] = join("",@inserts);
	$hitseqs[$k] =~ tr/\#//d; # remove the '#' symbols inserted at the beginning and end
    }
}

# Remove gaps? Captialize?
if ($outformat eq "ufas") {    
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/a-z.-/A-Z/d;} # Transform to upper case and remove all gaps
} elsif ($outformat eq "fas") {
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/a-z./A-Z-/;}  # Transform to upper case
} elsif ($outformat eq "a3m") {
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/.//d;}        # Remove gaps aligned to inserts
}

# Write sequences into output file
if ($append) {open (OUTFILE, ">>$outfile") or die ("cannot open $outfile:$!\n");}
else {open (OUTFILE, ">$outfile") or die ("cannot open $outfile:$!\n");}
if ($outformat eq "psi") {
    for ($k=0; $k<$nhit; $k++) {
	$hitseqs[$k] =~ tr/./-/;
	printf(OUTFILE "%s %s\n",$hitnames[$k],$hitseqs[$k]);
    }
} else {
    for ($k=0; $k<$nhit; $k++) {
	printf(OUTFILE "%s\n%s\n",$hitnames[$k],$hitseqs[$k]);
    }
}
close OUTFILE;

if ($v>=2) {printf("$nhit sequences extracted from $infile and written to $outfile\n");}

# Return number of hits in one byte (0-255)
if    ($nhit<110)  {exit $nhit;}
elsif ($nhit<1100) {exit 100+int($nhit/10);}
elsif ($nhit<5500) {exit 200+int($nhit/100);}
else                {exit 255;}              
exit(0);

sub System() {
    if ($v>=3) {print("$_[0]\n");} 
    return system($_[0])/256;
}

