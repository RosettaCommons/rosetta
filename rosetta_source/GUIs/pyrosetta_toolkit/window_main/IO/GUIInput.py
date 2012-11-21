
#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/IO/GUIInput.py
## @brief  Class responsible for managing input variables of the GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
## @author Steven Combs (steven.combs1@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os.path

#Tkinter Imports
from Tkinter import StringVar
import tkFileDialog

#Toolkit Imports
from window_main import global_variables
from modules.ScoreBase import ScoreBase as score_base

class GUIInput:
    def __init__(self, toolkit):
        self.toolkit = toolkit; #Basically an AP of the toolkit
        self.pdb_path = StringVar(); self.pdb_path.set("0");
        self.PDBLIST = StringVar(); self.PDBLIST.set("")
        
        self.loop_start=StringVar(); #Start of Loop
        self.loop_end=StringVar(); #End of Loop
        self.loop_chain=StringVar(); #Chain of Loop
        self.loop_sequence=StringVar(); #Sequence in Entry
        self.loops_as_strings = []; #Array of Loops: start:end:chain
        self.loops = Loops()

#### POSE INPUT ####

    def load_pose(self):
        
        print self.pdb_path.get()
        pose_from_pdb(self.toolkit.pose, self.pdb_path.get())
        self.toolkit.native_pose.assign(self.toolkit.pose); #Set native pose for RMSD.

        print self.toolkit.pose
        self.toolkit.pymol_class.SendNewPose()
        self.ScoreBaseObject = score_base(self.toolkit.pose, self.toolkit.score_class.score); #Score Base object for controling scoring of loops.
        pdbname = os.path.basename(self.pdb_path.get())
        pdbname = pdbname.split(".")[0]
        self.toolkit.output_class.outname.set(pdbname)
        self.toolkit.DesignDic = dict()
        
    def choose_load_pose(self):
        """
        Loads a File through the tk File Dialog
        """
        f = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title='Pick a file')
        if not f:return
        self.pdb_path.set(f)
        global_variables.current_directory= os.path.dirname(f)
        print global_variables.current_directory
        self.load_pose()
    
    def set_PDBLIST(self):
        f = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory,title='Pick a directory')
        if not f:return
        global_variables.current_directory =os.path.dirname(f)
        print "PDBLIST set"
        self.PDBLIST.set(f)