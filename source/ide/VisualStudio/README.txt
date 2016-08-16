(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under license.
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org. Questions about this can be
(c) addressed to University of Washington CoMotion, email: license@uw.edu.

Build Instructions:

You may have to update the VC++ project files for the rosetta libraries.
To do so, run the shell script, rosetta/rosetta_source/Visual Studio/makeFiles.csh

For Boinc Builds do the following first:

1. Get boinc source code.
   svn co http://boinc.berkeley.edu/svn/trunk/boinc
   
2. Place boinc directory in rosetta/rosetta_source/external/

3. For graphics, place glut directory in rosetta/rosetta_source/external/


Open rosetta/rosetta_source/Visual Studio/boinc_mini.sln with Visual Studio. Select appropriate
configuration and then build.


