#!/usr/bin/perl -W 
###############################################################################
# written by Grant S. Murphy 2008                                             # 
# Script gets RMSD and SCORE info from PDBS created in Rosetta and prints     #
# info in the format used by Matt O'Meara's RMSD_SCORE_PLOTTER in pymol       # 
# ./RMSD_SCORE_WRITER_FOR_PYMOL.pl LIST_ofPDBS.txt                            #
# ListofPDBS.txt needs to include the full path to the pdbs                   #
###############################################################################




MAIN: {

   my $listOfIds = $ARGV[0];
   if (@ARGV < 1) {
      print "This script gets RMSD and SCORE info from pdbs created in Rosetta";
      print "./RMSD_SCORE_WRITER_FOR_PYMOL.pl LIST_ofPDBS.txt";
      print "The list needs to include the full path to the pdbs";
      exit;
   }

   my $returnvalue;
   $returnvalue = open( MY_FILE, "$listOfIds" );
   if (!defined($returnvalue)) {
      print STDERR "Unable to open input file: $!";
      exit;
   }


  use File::Slurp;




   my @PDBS = <MY_FILE>;
   chomp @PDBS;
#   print "@PDBS \n";

   foreach my $PDB (@PDBS) {

  my $text = read_file( $PDB ) ;
  my @lines = read_file( $PDB ) ;
  my $size = @lines;
#  print "SIZE $size \n";
if ($size > 1 ) {
 my  @rms = grep (/rms/,@lines);
chomp @rms;
#print "$rms[0] \n";
@split_rms = split(/:/,$rms[0]);
 my  @score = grep (/score/,@lines);
chomp @score;
#print "$score[0] \n";
@split_score = split(/:/,$score[0]);
#print "$PDB $rms[0] $score[0] \n";
$PDB =~ s/^\s+//;
$split_score[1] =~ s/^\s+//;
$split_rms[1] =~ s/^\s+//;

print "$PDB $split_score[1] $split_rms[1]\n";
}

}
}
