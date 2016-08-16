#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/tests/Test_Main_GUI.py
## @brief  Test the frames of the Main GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

#Python Imports
import os
import sys

#if sys.platform == "darwin": sys.exit(0)  # skipping this test on Mac OS due to memory error (pointer being freed was not allocated)


from Tkinter import *

#Must have Toolkit added to path!

#Toolkit Imports
from app.pyrosetta_toolkit.pyrosetta_toolkit import MainTracer
from app.pyrosetta_toolkit.pyrosetta_toolkit import main_window


from rosetta import *
rosetta.init()


print "Running main GUI window"
main_window_class = main_window()
main_window_class.TR = MainTracer(main_window_class.output_textbox)
rosetta.basic.Tracer.set_ios_hook(main_window_class.TR, rosetta.basic.Tracer.get_all_channels_string(), False)
   #rosetta.init(extra_options="-mute all")


main_window_class.run(False)


print "Testing main window frames:"
#Here we figure out how to simulate button presses.
main_window_class.quit()
