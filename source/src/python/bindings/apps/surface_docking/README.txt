PyRosetta RosettaSurface Setup

Before running the scripts for the first time, make sure that all the following details are addressed.
 
A) When PyRosetta is installed, the install directory may be different from user to user, so adding the scripts directories to the PATH environment variable in .bashrc file will allow the scripts to be run from any directory, independent from the actual location of the scripts. The following instructions show how the directories are added:

    1) Create/open .bashrc in home directory
        vi ~/.bashrc
        
    2) If file exists, skip to 3. Else, add the following lines to the file:
        # .bashrc

        # Source global definitions
        if [ -f /etc/bashrc ]; then
            . /etc/bashrc
        fi

        export PATH=$PATH
        
    3) Add all directories containing scripts to end of PATH statement, delimited by a colon
    
        export PATH=$PATH:/path/to/scripts:/path/to/more/scripts
        
        where /path/to/scripts and /path/to/more/scripts should be modified to the correct directories.

    4) Make sure that the SetPyRosettaEnvironment.sh is sourced every session. Add 

        source /path/to/PyRosetta/SetPyRosettaEnvironment.sh
                
    5) Save and close file, then source it.
        :wq
        source ~/.bashrc

B) Make sure the python scripts are given the permission to be executed. Run the following command in the directory containing the scripts
        chmod +x *.py
        
C) Gnuplot is the plotting program used to generate the plots, so the program path has to be modified in PostProcessRS.sh for the plots to be generated automatically. The version used is 4.2.