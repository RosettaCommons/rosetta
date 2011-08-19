#!/usr/bin/perl

#script to read in a file containing a data matrix and another file containing requirements for the values in the columns of the matrix
#and then only select those rows of the matrix that meet all the requirements
#
#Florian Richter (floric@u.washington.edu, june 2007 )

use strict;

#function to calculate standard deviation of array of values. IMPORTANT. the average value has to be in the 0 element of the array
sub StDev {
    my @Values = @_;
    my $NumValues = (scalar @Values) - 1;
    #printf STDERR " hiargh %s %s \n",$NumValues, @Values[2];
    my $SDsum = 0;
    for(my $i = 1; $i <= $NumValues; $i++) {
	my $SQelement = $Values[$i] - $Values[0];
	$SDsum = $SDsum + ($SQelement * $SQelement);
    }
    return sqrt($SDsum/$NumValues);
}

#function to calculate stdDev, expects reference to an array with all the values and the average value as input
sub StdDev {

    my $RefValues = shift;
    my $AvValue = shift;

#    printf STDERR "refvalues is %s\n", $RefValues;
    my @Values = @$RefValues;

    my $NumValues = (scalar @Values) - 1;

    #printf STDERR "av value received is %s, number of values in array is %s. \n" , $AvValue, $NumValues;
    my $SDsum = 0;
    for(my $i = 0; $i <= $NumValues; $i++) {
	my $SQelement = $Values[$i] - $AvValue;
	$SDsum = $SDsum + ($SQelement * $SQelement);
    }
    return sqrt($SDsum/$NumValues);

}

#function: Padsp(InpString,len) adds spaces to the end of the input string until the desired length is reached
sub Padsp {
    my $InpString = $_[0];
    my $newlen = $_[1];
    
    my $origlen=length($InpString);
    for (my $i=0; $i<($newlen-$origlen); $i++) {
	$InpString=$InpString." ";
    }
    return $InpString;
}


sub usage {
  printf STDERR "\n";
  printf STDERR "usage: -d <datafile> -c <file with requirements> \n";
  printf STDERR "\n";
  exit 1;
}

#printf STDERR " $#ARGV \n";
if($#ARGV < 2) {&usage(); exit 1}


my $DataFile;
my @DataStrings = ();
my @DataArray = ();
my @ColumnArray = ();

my $RequireFile;
my $NumReqs = 0;
my @ReqStrings = ();
my $OutputSort = 0;
my $OutputColumn = 0;
my $OutOption = 0;
my $ReducedOutput = 0;
my $EZlabels = 1;
my $sort_data = 0;

my $tag_column = 0;

for(my $ii = -1; $ii < $#ARGV;$ii++){
    if ($ARGV[$ii] eq '-d'){$DataFile = $ARGV[$ii+1];}
    if ($ARGV[$ii] eq '-c'){$RequireFile = $ARGV[$ii+1];}
    if ($ARGV[$ii] eq '-output_filename_only'){$OutOption = 1;}
    if ($ARGV[$ii] eq '-short_output'){$ReducedOutput = 1;}
    if ($ARGV[$ii] eq '-easy_labels'){$EZlabels = 1;}
    if ($ARGV[$ii] eq '-column_id_labels'){$EZlabels = 0;}

    if ($ARGV[$ii] eq '-tag_column'){
      if( $ARGV[$ii+1] eq 'last' ) { $tag_column = -1;}
      else { $tag_column = $ARGV[$ii+1]; }
    }
}


#----block to read in the requirements for each column----

open(REQFILE,$RequireFile) || die "Could not open $RequireFile\n";
while(<REQFILE>) {

    chomp;
    my $inline = $_;

    if(substr($inline,0,3) eq 'req') {
	$ReqStrings[$NumReqs] = $inline;
	$NumReqs++;
    }
    if(substr($inline,0,6) eq 'output') {
	my @outarray = split(' ',$inline);
	$OutputSort = $outarray[1];
	$sort_data = 1;
	if($EZlabels == 1) { $OutputColumn = $outarray[2]; }
	else{ $OutputColumn = $outarray[2] - 1; }#remember offset
    }
}
close REQFILE;    

#requirements read in, now read in data

open(DATAFILE,$DataFile) || die "Could not open $DataFile\n"; 
@DataStrings = <DATAFILE>;
close DATAFILE;
my %title_hash = ();
my @title_array = split(' ',$DataStrings[0]);
my @req_columns = ();
for(my $ii = 0; $ii < scalar @title_array; $ii++){
  $title_hash{ $title_array[$ii] } = $ii;
  $req_columns[$ii] = 0;
}

if($EZlabels == 1){ #have to process title line to figure out which column number belongs to which label
  printf STDERR "output label %s ", $OutputColumn;
  $OutputColumn = $title_hash{ $OutputColumn };
  printf STDERR "translates to column %s.\n",$OutputColumn;
  $req_columns[$OutputColumn] = 1;
}
my @low_values;
my @high_values;
my @num_passing_reqs;


my $NumRows = scalar @DataStrings;

my $NumColumns;
my $NumCommentLines = 0;

for(my $ii = 0; $ii < $NumRows; $ii++){

    if( (substr($DataStrings[$ii],0,1) eq '#') || (substr($DataStrings[$ii],0,6) eq 'SCORES' ) || (substr($DataStrings[$ii],0,15) eq '    total_score')) {$NumCommentLines++;}

    else{
	my @CurLineArray = split(' ',$DataStrings[$ii]);

	 #first data line, check how many columns there are and initialize high and low arrays
	if($ii == $NumCommentLines) {
	  $NumColumns = scalar @CurLineArray;
	  if( $tag_column == -1 ) { $tag_column = $NumColumns -1; }
	  @low_values = @CurLineArray;
	  @high_values = @CurLineArray;
	}

	if( $NumColumns != scalar @CurLineArray) {printf STDERR "File $DataFile is corrupted. Rows don't have equal number of columns. First offending entry is line $ii\n"; exit 1;}
	
	$CurLineArray[$NumColumns] = $ii;    #keeping track of what DataArray line belongs to what DataStrings line
	$CurLineArray[$NumColumns+1] = 1;    #indicator of wheter this line fullfills all the requirements
	@{$DataArray[$ii - $NumCommentLines]} = @CurLineArray;

	#buildup column arrays
	for(my $jj = 1; $jj <= $NumColumns; $jj++){
	    $ColumnArray[$jj][$ii-$NumCommentLines] = $CurLineArray[$jj-1];
	    $ColumnArray[0][$jj] = $ColumnArray[0][$jj] + $CurLineArray[$jj-1]; #keep average values in 0th line of column array
	    if( $low_values[$jj-1] > $CurLineArray[$jj-1] ) { $low_values[$jj-1] = $CurLineArray[$jj-1]; }
	    if( $high_values[$jj-1] < $CurLineArray[$jj-1] ){  $high_values[$jj-1] = $CurLineArray[$jj-1]; }
	} 
    }
}

my $NumEntries = $NumRows - $NumCommentLines;

#calculate average values
for(my $jj = 1; $jj <= $NumColumns; $jj++){
    $ColumnArray[0][$jj] = $ColumnArray[0][$jj] / $NumEntries;
    $num_passing_reqs[$jj - 1 ] = $NumEntries;
}

printf STDERR "Data read in, average values for $NumEntries entries with $NumColumns columns determined, now checking for satisfying entries\ncolumn averages: ";
for(my $jj = 1; $jj <= $NumColumns; $jj++){
    printf STDERR "%s ",sprintf("%.2f",$ColumnArray[0][$jj]);
}
printf STDERR "\n";


#---now perform exclusions----

for(my $rr = 0; $rr < $NumReqs; $rr++) {
    my @CurReqArray = split(' ',$ReqStrings[$rr]);
    my $CurCol = -1;
    if($EZlabels == 1){ 
      $CurCol = $title_hash{ $CurReqArray[1] }; 
      if( $CurCol eq "" ){ 
	printf STDERR "\nError: label %s was not found in the data file.\n", $CurReqArray[1];
	exit 1;
      }
      printf STDERR "Label %s translates to column %s, first value is %s \n", $CurReqArray[1],$CurCol, $ColumnArray[$CurCol+1][0]; 
    }
    else{ $CurCol = $CurReqArray[1] - 1;} #remember offset
    $req_columns[$CurCol] = 1;			
    my $CurMode = $CurReqArray[2];
    my $CurSubmode = $CurReqArray[3];
    my $CurValue = $CurReqArray[4];
    printf STDERR "\n$CurCol $CurMode $CurSubmode $CurValue ";

    if($CurMode eq 'value'){    #absolute value mode, select only those entries that have an absolute value higher, lower or equal to a given cutoff
      my $num_passing = $NumEntries;
	if($CurSubmode eq '>'){
	    for(my $ii = 0; $ii < $NumEntries; $ii++){
		if($DataArray[$ii][$CurCol ] <= $CurValue){
		  $DataArray[$ii][$NumColumns+1] = 0;
		  $num_passing = $num_passing - 1;
		}
	    }
	}
	elsif($CurSubmode eq '<'){
	    for(my $ii = 0; $ii < $NumEntries; $ii++){
		if($DataArray[$ii][$CurCol] >= $CurValue){
		  $DataArray[$ii][$NumColumns+1] = 0;
		  $num_passing = $num_passing - 1;
		}
	    }
	
	}
	elsif($CurSubmode eq '=') { #value has to equal something
	   for(my $ii = 0; $ii < $NumEntries; $ii++){
	       #printf STDERR "comparing %s and %s ", $DataArray[$ii][$CurCol], $CurValue;
	       if($DataArray[$ii][$CurCol] != $CurValue){
		 $DataArray[$ii][$NumColumns+1] = 0;
		 $num_passing = $num_passing - 1;
	       }
	       #else{ printf STDERR "accept \n"; }
	    }
       }
      @num_passing_reqs[$CurCol] = $num_passing
    } #if mode eq value

    if($CurMode eq 'fraction'){     #percentage mode, select only those entries which are among the best or worst $CurValue percent entries in a column
	
	my @CurColSort = sort {$a <=> $b} @{$ColumnArray[$CurCol +1]};
	
	if((scalar @CurColSort) != $NumEntries){printf STDERR "Error, not enough entries for column $CurCol.\n"; exit 1;} #sanity check

	$CurColSort[$NumEntries] = $CurColSort[$NumEntries - 1] + 1; #have to put sentinel at end of array that's the biggest number
	my $NumToSelect = sprintf("%.0f",$NumEntries * $CurValue);
	#my $NumToSelect = $NumEntries * $CurValue;
	my $CutOffValue = 0;
	printf STDERR "numto select is $NumToSelect, ";
	
	if($CurSubmode eq '>'){ #highest percentage
	    $CutOffValue = $CurColSort[$NumEntries - $NumToSelect];
	    printf STDERR "cut off is $CutOffValue \n";
	    for(my $ii = 0; $ii < $NumEntries; $ii++){
		if($DataArray[$ii][$CurCol ] < $CutOffValue){$DataArray[$ii][$NumColumns+1] = 0;}
	    }
	} 
	elsif($CurSubmode eq '<'){#lowest percentage
	    $CutOffValue = $CurColSort[$NumToSelect - 1];
	    printf STDERR "cut off is $CutOffValue \n";
	    for(my $ii = 0; $ii < $NumEntries; $ii++){
		if($DataArray[$ii][$CurCol ] > $CutOffValue){$DataArray[$ii][$NumColumns+1] = 0;}
	    }
	} 
	@num_passing_reqs[$CurCol] = $NumToSelect
    }
}
printf STDERR "\n";


#-----exlusions performed, now output
#first StdDevs, 
printf STDERR "The following averages, Standard Deviations, low/high values, and fraction passing  were obtained: \n Field                      Av     StdDev          low        high       num_passing_requirements \n";
for(my $kk = 0; $kk < $NumColumns; $kk++){
    #if($req_columns[$kk] == 1) {
	my @CurStArray = @{$ColumnArray[$kk + 1]};
	#printf STDERR "wtf?!? %s   ", @CurArray[0];
	my $RefCurStArray;
	#printf STDERR " huh? numbers in array is %s \n", scalar @CurArray;
	@{$RefCurStArray} = @CurStArray;
	#printf STDERR "About to call StdDev\n";
	my $CurStdDev = &StdDev( $RefCurStArray, $ColumnArray[0][$kk + 1] );
	printf STDERR "%s    %.2f      %.2f        %.2f       %.2f         %s\n", &Padsp($title_array[$kk], 20),$ColumnArray[0][$kk + 1], $CurStdDev, $low_values[$kk], $high_values[$kk], $num_passing_reqs[$kk]
    #}
}
printf STDERR "\n";


my @SortDataArray = ();

if($sort_data){
    
    if($OutputSort eq 'sortmin'){ printf STDERR "output sorted ascending by column $OutputColumn( %s )\n", $title_array[ $OutputColumn ]; 
				  @SortDataArray = sort {$a->[$OutputColumn ] <=> $b->[$OutputColumn ]} @DataArray; 
				}
    elsif($OutputSort eq 'sortmax'){ printf STDERR "output sorted descending by column $OutputColumn\n"; @SortDataArray = sort {$b->[$OutputColumn ] <=> $a->[$OutputColumn ]} @DataArray; }
    else{@SortDataArray = @DataArray; printf STDERR "output line not understood, not sorting output.\n";}
}
else{@SortDataArray = @DataArray;}

if( ($OutOption != 1) && ($ReducedOutput != 1) && ($NumCommentLines > 0)) {printf STDOUT "%s",$DataStrings[0];} #assuming the first line is the title line
elsif($ReducedOutput == 1){
  my $outstring = $title_array[ $tag_column ]." ";
  for(my $kk = 0; $kk < $NumColumns; $kk++){
    if( ($kk != $tag_column ) && ($req_columns[$kk] == 1 )) { $outstring = $outstring.$title_array[$kk]." "; }
  }
  printf STDOUT "%s\n",$outstring; 
}


for(my $ii = 0; $ii < $NumEntries; $ii++){  
    
    if(($SortDataArray[$ii][$NumColumns + 1]) == 1) {    #means this fulfills all the requirements

	if($OutOption == 1) { printf STDOUT "%s\n",$SortDataArray[$ii][ $tag_column ];} #only print out filenames
	elsif($ReducedOutput == 1){   #only print columns of interest (as specified by reqfile)
	  my $outstring = $SortDataArray[$ii][ $tag_column ]." ";
	  for(my $kk = 0; $kk < $NumColumns; $kk++){
	    if( ($kk != $tag_column ) && ($req_columns[$kk] == 1 )) { $outstring = $outstring.$SortDataArray[$ii][$kk]." "; }
	  }
	  printf STDOUT "%s\n",$outstring;
	}
	else {printf STDOUT "%s",$DataStrings[$SortDataArray[$ii][$NumColumns]];}
    }
}
		
			
