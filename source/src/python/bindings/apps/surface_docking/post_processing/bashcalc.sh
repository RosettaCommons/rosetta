#!/bin/bash
# The "scale=4" allows for division to be carried out to 4 decimal places. 
# You can edit this to whatever length you wish to carry the division to.
# I suggest putting an alias into your .bashrc file
# Example: alias calc='sh /home/$USER/scripts/bashcalc.sh'
# This allows for your bashcalc.sh script to operate from the shell 
# Example: 
#   crouse@linux:~> calc 3.555+7.999
#   11.554
#   crouse@linux:~>
# Alternatively, you could also make this work system wide if you have
# root access, you could put the script into /usr/bin/calc. Whichever method
# you choose to use, don't forget to "chmod a+x" the script to make it executable.
echo "scale=4; $1" | bc ;exit
