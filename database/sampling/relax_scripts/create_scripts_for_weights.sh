#!/bin/bash

#It can be tedious to create a bunch of script files when all you really need to do is just change 4 numbers
#This script takes care of all the tedious work for you
#Example
#bash create_scripts_for_weights.sh rosettacon2018 0.079 0.295 0.577 1


#Author: Jack Maguire, August 2018
#Note: There are currently only 2 options, making 4 scripts. This is small enough where it's easiest to just go through each case explicitly. If we ever get more than a few options we should make this more intelligent
#Note: By convention, all extensions are listed alphabetically

if [[ "$#" -ne "5" ]]; then
    echo "Usage: bash create_scripts_for_weights.sh [prefix] [first fa_rep coefficient] [second fa_rep coefficient] [third fa_rep coefficient] [fourth fa_rep coefficient, which is usually 1]"
    echo " "
    echo "Example: bash create_scripts_for_weights.sh rosettacon2018 0.079 0.295 0.577 1"
    exit 1
fi

prefix=$1
weight1=$2
weight2=$3
weight3=$4
weight4=$5


#Normal File
echo "repeat %%nrepeats%%" > $prefix.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.txt
echo "accept_to_best" >> $prefix.txt
echo "endrepeat" >> $prefix.txt



#Dualspace File
echo "switch:torsion" > $prefix.dualspace.txt
echo "repeat %%nrepeats%%" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.dualspace.txt
echo "accept_to_best" >> $prefix.dualspace.txt
echo "endrepeat" >> $prefix.dualspace.txt
echo " " >> $prefix.dualspace.txt
echo "switch:cartesian" >> $prefix.dualspace.txt
echo "repeat 1" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.dualspace.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.dualspace.txt
echo "accept_to_best" >> $prefix.dualspace.txt
echo "endrepeat" >> $prefix.dualspace.txt



#Normal beta_nov16 File
echo "repeat %%nrepeats%%" > $prefix.beta_nov16.txt
echo "reference 0.3     3.1     -2.6     -2.55    4.8     -0.5    0.7      4.5     -1.6     4.0     3.9     -1.7     -2.0     -1.5     -1.0     -2.0    -2.0     4.0     9.0     3.7" >> $prefix.beta_nov16.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.beta_nov16.txt
echo "reference 2.2619  4.8148  -1.6204  -1.6058  2.7602  1.0350  1.3406   2.5006  -0.6895  1.9223  2.3633  -0.3009  -4.2787   0.1077   0.0423  -0.4390 -0.7333  3.2371  4.7077  2.3379" >> $prefix.beta_nov16.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.beta_nov16.txt
echo "reference 2.2619  4.5648  -1.6204  -1.6158  2.5602  1.1350  1.2406   2.3006  -0.7895  1.7223  2.1633  -0.3009  -4.3787   0.1077   0.0423  -0.4390 -0.7333  3.1371  4.4077  2.1379" >> $prefix.beta_nov16.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.beta_nov16.txt
echo "reference 2.2619  4.3148  -1.6204  -1.6358  1.9602  1.4350  0.8406   1.8006  -0.8895  1.3223  1.4633  -0.3009  -4.6787  -0.1077  -0.1423  -0.5390 -0.9333  2.7371  3.7077  1.7379" >> $prefix.beta_nov16.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.beta_nov16.txt
echo "accept_to_best" >> $prefix.beta_nov16.txt
echo "endrepeat" >> $prefix.beta_nov16.txt



#Dualspace beta_nov16 File
echo "switch:torsion" > $prefix.beta_nov16.dualspace.txt
echo "repeat %%nrepeats%%" >> $prefix.beta_nov16.dualspace.txt
echo "reference 0.3     3.1     -2.6     -2.55    4.8     -0.5    0.7      4.5     -1.6     4.0     3.9     -1.7     -2.0     -1.5     -1.0     -2.0    -2.0     4.0     9.0     3.7" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.8148  -1.6204  -1.6058  2.7602  1.0350  1.3406   2.5006  -0.6895  1.9223  2.3633  -0.3009  -4.2787   0.1077   0.0423  -0.4390 -0.7333  3.2371  4.7077  2.3379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.5648  -1.6204  -1.6158  2.5602  1.1350  1.2406   2.3006  -0.7895  1.7223  2.1633  -0.3009  -4.3787   0.1077   0.0423  -0.4390 -0.7333  3.1371  4.4077  2.1379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.3148  -1.6204  -1.6358  1.9602  1.4350  0.8406   1.8006  -0.8895  1.3223  1.4633  -0.3009  -4.6787  -0.1077  -0.1423  -0.5390 -0.9333  2.7371  3.7077  1.7379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.beta_nov16.dualspace.txt
echo "accept_to_best" >> $prefix.beta_nov16.dualspace.txt
echo "endrepeat" >> $prefix.beta_nov16.dualspace.txt

echo " " >> $prefix.beta_nov16.dualspace.txt

echo "switch:cartesian" >> $prefix.beta_nov16.dualspace.txt
echo "repeat 1" >> $prefix.beta_nov16.dualspace.txt
echo "reference 0.3     3.1     -2.6     -2.55    4.8     -0.5    0.7      4.5     -1.6     4.0     3.9     -1.7     -2.0     -1.5     -1.0     -2.0    -2.0     4.0     9.0     3.7" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight1 "0.01     1.0" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.8148  -1.6204  -1.6058  2.7602  1.0350  1.3406   2.5006  -0.6895  1.9223  2.3633  -0.3009  -4.2787   0.1077   0.0423  -0.4390 -0.7333  3.2371  4.7077  2.3379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight2 "0.01     0.5" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.5648  -1.6204  -1.6158  2.5602  1.1350  1.2406   2.3006  -0.7895  1.7223  2.1633  -0.3009  -4.3787   0.1077   0.0423  -0.4390 -0.7333  3.1371  4.4077  2.1379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight3 "0.01     0.0" >> $prefix.beta_nov16.dualspace.txt
echo "reference 2.2619  4.3148  -1.6204  -1.6358  1.9602  1.4350  0.8406   1.8006  -0.8895  1.3223  1.4633  -0.3009  -4.6787  -0.1077  -0.1423  -0.5390 -0.9333  2.7371  3.7077  1.7379" >> $prefix.beta_nov16.dualspace.txt
echo "ramp_repack_min" $weight4 "    0.00001  0.0" >> $prefix.beta_nov16.dualspace.txt
echo "accept_to_best" >> $prefix.beta_nov16.dualspace.txt
echo "endrepeat" >> $prefix.beta_nov16.dualspace.txt

