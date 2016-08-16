#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/window_main/global_variables.py
## @brief  Global variables for PyRosetta Toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#Globals that will change regularly by many functions and classes accross the GUI should go here.


#This gets set to the toolkit directory by pyrosetta_toolkit.py.
current_directory=""


#NOTE:  To use these global variables import the file into your module.  Use 'from window_main import global_variables' as the main toolkit directory is added to PythonPath at launch of the GUI
#These variables are truly global and just by importing the module, the instance of the variable can be accessed and manipulated.
