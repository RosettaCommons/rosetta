#!/usr/local/bin/bash
# Script file to execute all post-processing scripts at once
# Emily Koo

# Surface contact map frequency plot (ads only)
nohup PlotSurfaceContactMap.py &

# Secondary structure histogram plot
nohup PlotSecStruct.sh &

# Contact map frequency plot
nohup PlotContactMap.sh &

~/Applications/bin/gnuplot *.gnuplot
