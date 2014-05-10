#!/bin/bash
#Edited by Liza for gnuplot Version 4.2
# Further edited by Emily Koo

echo "Generating secondary structure plot"
ls *.pdb > PDBsInDirectory
DecoyCnt=$( wc -l PDBsInDirectory | cut -f1 -d" " )
current_structure=1

while [ $current_structure -le $DecoyCnt ]
  do
  current_structure_name=$( head -$current_structure PDBsInDirectory | tail -1 )	
  PhiPsi.py $current_structure_name
  current_structure=$[ $current_structure + 1 ]
done
	
ls *.dss > DssInDirectory
DssCnt=$( wc -l DssInDirectory | awk '{print $1}')
FirstDss=$( head -1 DssInDirectory )
LinesInDss=$( wc -l $FirstDss | awk '{print $1}' )
current_line=2

echo '0      0      0     0' >> SecStruct.dat
while [ $current_line -le $LinesInDss ]
  do
  current_dss=1
  Helix=0
  Loop=0
  Extended=0
  while [ $current_dss -le $DssCnt ]
    do
    current_dss_name=$( head -$current_dss DssInDirectory | tail -1 )
    current_res_SS_type=$( head -$current_line $current_dss_name | tail -1 | cut -c 28-36 )
    if [ $current_res_SS_type == "Helix" ]
	then
	Helix=$[ $Helix + 1 ]
    elif [ $current_res_SS_type == "Loop" ]
	then
	Loop=$[ $Loop + 1 ]
    elif [ $current_res_SS_type == "Extended" ]
	then
	Extended=$[ $Extended + 1 ]
    fi
    current_dss=$[ $current_dss + 1 ]
  done
  FractionHelix=$( sh bashcalc.sh $Helix/$DssCnt )
  FractionLoop=$( sh bashcalc.sh $Loop/$DssCnt )
  FractionExtended=$( sh bashcalc.sh $Extended/$DssCnt )
  current_residue=$[ $current_line - 1 ]
  #echo $current_residue
  echo $current_residue '   '  $FractionHelix  '     '    $FractionLoop   '    '  $FractionExtended >> SecStruct.dat
  current_line=$[ $current_line + 1 ]
done

rm DssInDirectory PDBsInDirectory *.dss 

x=$(wc -l SecStruct.dat | cut -f1 -d" ")
X=$(echo $x - 0.5 | bc)
#BoxSize=$(echo 1+"(2*$X/10)" | bc)
BoxSize=6
BoxPosition=$(echo $X/1.2 | bc)
OffSet=$( echo $BoxSize/2 | bc )
TextPosition=$( echo $BoxPosition-$OffSet+0.5 | bc )

# Generate GNU script to make the plot
echo set terminal png size 800, 800 nocrop notransparent enhanced font \"VeraBd,15\" > SecStruct.gnuplot #Emily
echo set output \"SecStruct.png\" >> SecStruct.gnuplot
echo unset key >> SecStruct.gnuplot #Version 4
echo set style data histogram >> SecStruct.gnuplot
echo set style fill solid 1.0 >> SecStruct.gnuplot
echo set style histogram rowstacked >> SecStruct.gnuplot
echo set xrange[0.5:$X] >> SecStruct.gnuplot
echo set xtics out nomirror >> SecStruct.gnuplot #Version 4
echo set mxtics >> SecStruct.gnuplot #Version 4
echo set ytics out nomirror >> SecStruct.gnuplot #Version 4
#echo show xtics >> SecStruct.gnuplot #Version 4
#echo show ytics >> SecStruct.gnuplot #Version 4
echo set border lw 3 >> SecStruct.gnuplot #Version 4
echo set obj 10 rect at $BoxPosition,0.9 size $BoxSize,0.125 >> SecStruct.gnuplot
echo set obj 10 fillstyle solid border -1 front >> SecStruct.gnuplot
echo set label \"Helix\" at $TextPosition,0.925 front textcolor lt 1 >> SecStruct.gnuplot
echo set label \"Turn\" at $TextPosition,0.895 front textcolor lt 2 >> SecStruct.gnuplot
echo set label \"Other\" at $TextPosition,0.865 front textcolor lt 3 >> SecStruct.gnuplot
echo set ylabel \"Fraction of decoys\" font \"VeraBd, 20\" >> SecStruct.gnuplot
echo set xlabel \"Residue\" font \"VeraBd, 20\" >> SecStruct.gnuplot
#echo set title \"Osteocalcin\(Crystal\)\" font \"VeraBd, 20\" >> SecStruct.gnuplot   
echo plot \'SecStruct.dat\' using 2, \'\' using 3, \'\' using 4 >> SecStruct.gnuplot #Version 4

echo "Sec struct plot DONE"
