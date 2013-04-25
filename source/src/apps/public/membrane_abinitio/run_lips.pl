#!/usr/bin/perl -w 
###############################################################################
# written by Bjorn Wallner and modified by Vladimir Yarov-Yarovoy             #
# Script takes sequence and transmembrane topology prediction files and       #
# outputs lipophilicity score (*.lips4) for input to rosetta-membrane         # 
# ab initio protocol                                                          #
###############################################################################

if ($#ARGV < 4) {
    print STDERR "usage: $0 <fasta file> <span file> <path to blastpgp> <path to nr database> <path to alignblast.pl script> \n";
    print STDERR "example: $0 BRD4.fasta BRD4.span /work/bjornw/Apps/blast/bin/blastpgp /scratch/shared/genomes/nr  ~bjornw/mini/src/apps/pilot/bjornw/alignblast.pl \n";
    exit -1;
}

$seq=$ARGV[0];
$spanfile=$ARGV[1];
$PSIBLAST=$ARGV[2]; # path to blastpgp
$NR=$ARGV[3]; # path to nr database
$parseblast=$ARGV[4]; #path to alignblast.pl
$sequence=`grep -v '>' $seq`;
$sequence=~s/\n//g;

$pdb_id=substr($seq,0,4);
@helix_starts=();
@helix_ends=();
if(-e $spanfile)
{
    open(SPAN,$spanfile) || die "Cannot open $spanfile\n";
    while(<SPAN>)
    {
	
	if(/\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+/)
	{
	    push(@helix_starts,$1);
	    push(@helix_ends,$2);
	    
	}
    }
    close(SPAN);
}

$len=length($sequence);
@tm=();
$tm="";
$helix_counter=1;
$old_membrane_region=0;
for(my $i=0;$i<$len;$i++)
{
    $pos=$i+1;
    $membrane_region=0;
    $ss="-";
    for(my $i=0;$i<scalar @helix_starts;$i++)
    {
	if($pos>=$helix_starts[$i] && $pos<=$helix_ends[$i])
	{
	    $membrane_region=$helix_counter;
	    $ss="H";
	    last;
	}
    }
    $tm.=$ss;
    push(@tm,$membrane_region);
    $helix_counter++ if($membrane_region == 0 && $old_membrane_region!=0);
    $old_membrane_region=$membrane_region;
}
$helix_counter-- if($old_membrane_region==0);

# psi-blast

$DB=$NR;
$blastout=$seq;
$blastout=~s/fasta/blast/g;

if(!-e $blastout)
{
    $psiblast_command="$PSIBLAST -i $seq -d $DB -j 2 -h 0.001 -b5000 -v5000 -o $blastout";
    print "Running:\n $psiblast_command\n";
    `$psiblast_command`;
}

if(!-e "$blastout.msa")
{
    `$parseblast $blastout $blastout.msa -psi`;
}

@ali=split(//,$sequence);
@ali_matrix=[@ali];
%ali=();
$seq_count=1;
open(MSA,"$blastout.msa");
while(<MSA>)
{
    ($key,$ali)=split(/\s+/);
    @ali=split(//,$ali);
    $ali{$key}=$ali;
    $ali_matrix[$seq_count]=[@ali];
    $seq_count++;
}

close(MSA);

my $lips_high="";
my $lips_low="";
open(OUT,">$pdb_id.lipo");
print OUT "Lipid exposed data: resnum mean-lipo lipophil entropy\n";
open(DATA,">$pdb_id.raw");
$url="http://gila.bioengr.uic.edu/cgi-bin/lips/script.cgi";
for($helix_num=1;$helix_num<=$helix_counter;$helix_num++)
{
    my $sequences="";
    my $first_num=$helix_starts[$helix_num-1];
    for(my $i=0;$i<$seq_count;$i++)
    {
	my $temp_str="";
	for(my $j=0;$j<$len;$j++)
	{
	    if($tm[$j] == $helix_num)
	    {
		$temp_str.=$ali_matrix[$i][$j];
	    }
	}
	if($i<5000 && not($temp_str=~/-/)) 
	{
	    $sequences.="$temp_str\n";
	}
    }
    $data=`curl -s $url -d sequence='$sequences' -d num=$first_num`;
   
    print DATA $data;
   
    my ($lips_high_temp,$lips_low_temp,$raw_data)=parse_lips($data);  
   
    $raw_data=~s/helix_num/$helix_num/g;
    print OUT $raw_data;
    if(length($lips_high)>0)
    {
    	$lips_high.="+$lips_high_temp";
    	$lips_low.="+$lips_low_temp";
    }
    else
    {
    	$lips_high.="+$lips_high_temp";
    	$lips_low.="+$lips_low_temp";
    }
}
close(DATA);
close(OUT);
open(OUT,">".$pdb_id.".pdb.color");
print OUT "select high, resi $lips_high\n";
print OUT "select low, resi $lips_low\n";
close(OUT);

my $outfile_base="$pdb_id";
my $data=`cat $pdb_id.raw`;
my $m=4;
my ($raw_data)=parse_lips3($data,$m);
open(OUT,">$outfile_base.lips$m");
print OUT "Lipid exposed data: resnum mean-lipo lipophil entropy\n";
print OUT $raw_data;
close(OUT);

sub parse_lips
{
    my $data=shift;
    my @data=split(/\n/,$data);
    $surface=0;
    my @resnum=();
    my @lips=(0,0,0,0,0,0,0);
    my $return_str="";
    my $highest_lipo=-99999;
    my $lowest_lipo=99999999;
    my $highest_lip_index="";
    my $lowest_lip_index="";
    foreach my $line(@data)
    {
	if($line=~/SURFACE\s+(\d):/)
	{
	 
	    $surface=$1;
	}
	if($line=~/(\d+)\s+[A-Z]\s+[\d\.\-]+/)
	{
	    push(@{$resnum[$surface]},$1) #+110);
	}
	if($line=~/(\d+)\s+[\d\.\-]+\s+[\d\.\-]+\s+([\d\.\-]+)/)
	{
	    $lips[$1]=$2;
	    if($lips[$1]>$highest_lipo)
	    {
		$highest_lipo=$lips[$1];
		$highest_lipo_index=$1;
	    }
	    if($lips[$1]<$lowest_lipo)
	    {
		$lowest_lipo=$lips[$1];
		$lowest_lipo_index=$1;
	    }
	}
    }

    my $res_high=join('+',@{$resnum[$highest_lipo_index]});
    my $res_low=join('+',@{$resnum[$lowest_lipo_index]});
    
    for(my $i=0;$i<7;$i++)
    {
	my $res=join('+',@{$resnum[$i]});
	$return_str.="helix_num $lips[$i] ( $res )\n";
    }
    
    return ($res_high,$res_low,$return_str);
}

sub parse_lips2
{
    my $data=shift;
    my @data=split(/\n/,$data);
    $surface=0;
    my @resnum=();
    my @lipophilicity=();
    my @entropy=();
    my @lips=(0,0,0,0,0,0,0);
    my $return_str="";
    my $highest_lipo=-99999;
    my $lowest_lipo=99999999;
    my $highest_lip_index="";
    my $lowest_lip_index="";
    foreach my $line(@data)
    {
	if($line=~/SURFACE\s+(\d):/)
	{
	 
	    $surface=$1;
	}
	if($line=~/(\d+)\s+[A-Z]\s+([\d\.\-]+)\s+([\d\.\-]+)/)
	{
	    push(@{$resnum[$surface]},$1);
	    $lipophilicity=2*$2;  #For some reason the webserver multiply the score by 2.
	    $entropy=$3;
	    push(@{$lipophilicity[$surface]},$lipophilicity);
	    push(@{$entropy[$surface]},$entropy);
	}
	if($line=~/(\d+)\s+[\d\.\-]+\s+[\d\.\-]+\s+([\d\.\-]+)/)
	{
	    $lips[$1]=$2;
	    if($lips[$1]>$highest_lipo)
	    {
		$highest_lipo=$lips[$1];
		$highest_lipo_index=$1;
	    }
	    if($lips[$1]<$lowest_lipo)
	    {
		$lowest_lipo=$lips[$1];
		$lowest_lipo_index=$1;
	    }
	}
    }

    my $mean_lips=sprintf("%5.3f",mean(@lips));
    my $std_lips=sprintf("%5.3f",std(@lips));

    my $res_high=join('+',@{$resnum[$highest_lipo_index]});
    for(my $i=0;$i<scalar(@{$resnum[$highest_lipo_index]});$i++)
    {

	$return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$highest_lipo_index]}[$i], $lips[$highest_lipo_index], ${$lipophilicity[$highest_lipo_index]}[$i], ${$entropy[$highest_lipo_index]}[$i]);
    }
    for(my $i=0;$i<scalar(@{$resnum[$lowest_lipo_index]});$i++)
    {

	$return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$lowest_lipo_index]}[$i], -$lips[$lowest_lipo_index], ${$lipophilicity[$lowest_lipo_index]}[$i], ${$entropy[$lowest_lipo_index]}[$i]);
    }
    my $res_low=join('+',@{$resnum[$lowest_lipo_index]});
    return ($res_high,$res_low,$return_str);
}

sub parse_lips3
{
    my $data=shift;
    my $mode=shift;
    my @data=split(/\n/,$data);
    $surface=0;
    my @resnum=();
    my @lipophilicity=();
    my @entropy=();
    my @lips=(0,0,0,0,0,0,0);
    my $return_str="";
    my $highest_lipo=-99999;
    my $lowest_lipo=99999999;
    my $highest_lip_index="";
    my $lowest_lip_index="";
    my @nb_list=();
    @{$nb_list[0]}=(3,4);
    @{$nb_list[1]}=(4,5);
    @{$nb_list[2]}=(5,6);
    @{$nb_list[3]}=(0,6);
    @{$nb_list[4]}=(0,1);
    @{$nb_list[5]}=(1,2);
    @{$nb_list[6]}=(2,3);


    foreach my $line(@data)
    {
	if($line=~/SURFACE\s+(\d):/)
	{
	 
	    $surface=$1;
	}
	if($line=~/(\d+)\s+[A-Z]\s+([\d\.\-]+)\s+([\d\.\-]+)/)
	{
	    push(@{$resnum[$surface]},$1);
	    $lipophilicity=2*$2;  #For some reason the webserver multiply the score by 2.
	    $entropy=$3;
	    push(@{$lipophilicity[$surface]},$lipophilicity);
	    push(@{$entropy[$surface]},$entropy);
	}
	if($line=~/(\d+)\s+[\d\.\-]+\s+[\d\.\-]+\s+([\d\.\-]+)/)
	{
	    $lips[$1]=$2;
	    if($lips[$1]>$highest_lipo)
	    {
		$highest_lipo=$lips[$1];
		$highest_lipo_index=$1;
	    }
	    if($lips[$1]<$lowest_lipo)
	    {
		$lowest_lipo=$lips[$1];
		$lowest_lipo_index=$1;
	    }
	}
    
	if($line=~/<\/body>/)
	{

	    my $mean_lips=sprintf("%5.3f",mean(@lips));
	    my $std_lips=sprintf("%5.3f",std(@lips));
	    
	    my $res_high=join('+',@{$resnum[$highest_lipo_index]});
	    if($mode eq "1")
	    {
		for(my $i=0;$i<scalar(@{$resnum[$highest_lipo_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$highest_lipo_index]}[$i], $lips[$highest_lipo_index], ${$lipophilicity[$highest_lipo_index]}[$i], ${$entropy[$highest_lipo_index]}[$i]);
	    }
	    for(my $i=0;$i<scalar(@{$resnum[$lowest_lipo_index]});$i++)
	    {
		
		$return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$lowest_lipo_index]}[$i], -$lips[$lowest_lipo_index], ${$lipophilicity[$lowest_lipo_index]}[$i], ${$entropy[$lowest_lipo_index]}[$i]);
	    }
	    }
	    if($mode eq "2")
	    {
		for(my $i=0;$i<scalar(@{$resnum[$highest_lipo_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$highest_lipo_index]}[$i]+110, $lips[$highest_lipo_index], ${$lipophilicity[$highest_lipo_index]}[$i], ${$entropy[$highest_lipo_index]}[$i]);
		}
		
		#the other surface should be lowest lipo + not neighbor to highest.
		my $selected_index=0;
		$lowest_lips=9999999999;
		for(my $i=0;$i<scalar @lips;$i++)
		{
		   # print "$i $highest_lipo_index $nb_list[$highest_lipo_index][0] $nb_list[$highest_lipo_index][1]\n";
		    if($i!=$nb_list[$highest_lipo_index][0] &&
		       $i!=$nb_list[$highest_lipo_index][1] &&
		       $lips[$i]<$lowest_lips)
		    {
			$selected_index=$i;
			$lowest_lips=$lips[$i];
		    }
		}
	#	print "$highest_lipo_index $selected_index $lowest_lipo_index\n";
	    
		for(my $i=0;$i<scalar(@{$resnum[$selected_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$selected_index]}[$i], -$lips[$selected_index], ${$lipophilicity[$selected_index]}[$i], ${$entropy[$selected_index]}[$i]);
		}
	    }
	    if($mode eq "3")
	    {
		for(my $i=0;$i<scalar(@{$resnum[$lowest_lipo_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$lowest_lipo_index]}[$i], -$lips[$lowest_lipo_index], ${$lipophilicity[$lowest_lipo_index]}[$i], ${$entropy[$lowest_lipo_index]}[$i]);
		}
		
		#the other surface should be highest lipo + not neighbor to lowest.
		my $selected_index=0;
		$highest_lips=-9999999999;
		for(my $i=0;$i<scalar @lips;$i++)
		{
		   # print "$i $lowest_lipo_index $nb_list[$lowest_lipo_index][0] $nb_list[$lowest_lipo_index][1]\n";
		    if($i!=$nb_list[$lowest_lipo_index][0] &&
		       $i!=$nb_list[$lowest_lipo_index][1] &&
		       $lips[$i]>$highest_lips)
		    {
			$selected_index=$i;
			$highest_lips=$lips[$i];
		    }
		}
	#	print "$lowest_lipo_index $selected_index $highest_lipo_index\n";
	    
		for(my $i=0;$i<scalar(@{$resnum[$selected_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$selected_index]}[$i], $lips[$selected_index], ${$lipophilicity[$selected_index]}[$i], ${$entropy[$selected_index]}[$i]);
		}
	    }
	    if($mode eq "4")
	    {
	
		#exit;
		#the other surface should be highest lipo + not neighbor to lowest and its neighbor.
		my $selected_index=0;
		$highest_lips=-9999999999;
		for(my $i=0;$i<scalar @lips;$i++)
		{
		   # print "$i $lowest_lipo_index $nb_list[$lowest_lipo_index][0] $nb_list[$lowest_lipo_index][1]\n";
		    if($i!=$nb_list[$lowest_lipo_index][0] &&
		       $i!=$nb_list[$lowest_lipo_index][1] &&
		       $nb_list[$i][0]!=$nb_list[$lowest_lipo_index][0] &&
		       $nb_list[$i][1]!=$nb_list[$lowest_lipo_index][1] &&
		       $nb_list[$i][0]!=$nb_list[$lowest_lipo_index][1] &&
		       $nb_list[$i][1]!=$nb_list[$lowest_lipo_index][0] &&
		       $lips[$i]>$highest_lips)
		    {
			$selected_index=$i;
			$highest_lips=$lips[$i];
		    }
		}
		print "$lowest_lipo_index $selected_index $highest_lipo_index\n";
	    

		$n1=$nb_list[$lowest_lipo_index][0];
		$n2=$nb_list[$lowest_lipo_index][1];
		my $normal25=0.5;
		for(my $i=0;$i<scalar(@{$resnum[$lowest_lipo_index]});$i++)
		{
		    
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$lowest_lipo_index]}[$i], -1, ${$lipophilicity[$lowest_lipo_index]}[$i], ${$entropy[$lowest_lipo_index]}[$i]);
		}
		for(my $i=0;$i<scalar(@{$resnum[$n1]});$i++)
		{
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$n1]}[$i], -$normal25, ${$lipophilicity[$n1]}[$i], ${$entropy[$n1]}[$i]);
		}
		for(my $i=0;$i<scalar(@{$resnum[$n2]});$i++)
		{
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$n2]}[$i], -$normal25, ${$lipophilicity[$n2]}[$i], ${$entropy[$n2]}[$i]);
		}

		$ns1=$nb_list[$selected_index][0];
		$ns2=$nb_list[$selected_index][1];
		for(my $i=0;$i<scalar(@{$resnum[$selected_index]});$i++)
		{
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$selected_index]}[$i], 1, ${$lipophilicity[$selected_index]}[$i], ${$entropy[$selected_index]}[$i]);
		}
		for(my $i=0;$i<scalar(@{$resnum[$ns1]});$i++)
		{
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$ns1]}[$i], $normal25, ${$lipophilicity[$ns1]}[$i], ${$entropy[$ns1]}[$i]);
		}
		for(my $i=0;$i<scalar(@{$resnum[$ns2]});$i++)
		{
		    $return_str.=sprintf("%7d %7.3f %7.3f %7.3f\n",${$resnum[$ns2]}[$i], $normal25, ${$lipophilicity[$ns2]}[$i], ${$entropy[$ns2]}[$i]);
		}

	    }
	    

	  #  print "=====\n";
	    @resnum=();
	    @lipophilicity=();
	    @entropy=();
	    @lips=(0,0,0,0,0,0,0);
	    $highest_lipo=-99999;
	    $lowest_lipo=99999999;
	    $highest_lip_index="";
	    $lowest_lip_index="";
	}
    }
    return $return_str;
}

sub std
{
    my @data=@_;
    my $mean=mean(@data);
    my $n=scalar @data;
    if($n!=1)
    {
	my $sum=0;
	foreach my $term(@data)
	{
	    $sum+=($term-$mean)*($term-$mean);
	}
	return sqrt(1/($n-1)*$sum);
    }
    else
    {
	return 0;
    }
}

sub mean
{
    my @data=@_;
    my $sum=0;
    my $number_of_elements=scalar @data;
    foreach my $term(@data)
    {
	$sum+=$term;
    }
    if($number_of_elements==0)
    {
	return 0;
    }
    else
    {
	return $sum/$number_of_elements;
    }
}
