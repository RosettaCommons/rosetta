PyRosetta RosettaSurface Setup

Before running the scripts for the first time, make sure that all the following details are addressed.
 
A) When PyRosetta is installed, the install directory may be different from user to user, so make sure to acticate PyRosetta Python environment before running these scripts.

B) Make sure the python scripts are given the permission to be executed. Run the following command in the directory containing the scripts
        chmod +x *.py
        
C) Gnuplot is the plotting program used to generate the plots, so the program path has to be modified in PostProcessRS.sh for the plots to be generated automatically. The version used is 4.2.
