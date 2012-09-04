#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/input.py
## @brief  Handles input and export of pdbs, xml files, options, etc 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
## @author Steven Combs (steven.combs1@gmail.com)

from rosetta import *
from Tkinter import *
from modules import *
from modules.ScoreBase import ScoreBase as score_base
from modules.tools import input as input_tools
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
import os

class InputFiles():

    def __init__(self, main, toolkit):
        self.toolkit = toolkit
        #self.pwd= self.location()[0]
        self.current_directory = StringVar()
        self.indir = StringVar()
        self.filename = StringVar(); self.filename.set("0");
        self.PDBLIST = StringVar(); self.PDBLIST.set("")
        #self.loops_as_strings = []; #Array of Loops: start:end:chain
        #self.VarStartLoop=StringVar(); #Start of Loop
        #self.VarEndLoop=StringVar(); #End of Loop
        #self.VarLoopChain=StringVar(); #Chain of Loop
        #self.VarLoopSeq=StringVar(); #Sequence in Entry
        #self.loops = Loops()



#### POSE INPUT ####

    def loadp(self):
        
        print self.filename.get()
        pose_from_pdb(self.toolkit.pose, self.filename.get())
        self.toolkit.native_pose.assign(self.toolkit.pose); #Set native pose for RMSD.
        #Note Dec2011 - Giving up on the PList.  Don't really want to try to integrate it....
       #p= Pose()
        
        print self.toolkit.pose
        self.toolkit.PyMOLObject.SendNewPose()
        self.ScoreBaseObject = score_base(self.toolkit.pose, self.toolkit.ScoreObject.score); #Score Base object for controling scoring of loops.

        self.DesignDic = dict(); #Clears the Design Dictionary upon loading of new pose
        #Design1(self.main).clearResData()
        #return self.DesignDic
        pdbname = os.path.basename(self.filename.get())
        pdbname = pdbname.split(".")[0]
        self.toolkit.output_class.outname.set(pdbname)
        self.toolkit.DesignDic = dict()
    def choose_load_pose(self):
        """
        Loads a File through the tk File Dialog
        """
        f = tkFileDialog.askopenfilename(initialdir=self.current_directory.get(), title='Pick a file')
        if not f:
            return
        self.filename.set(f)
        self.current_directory.set(os.path.dirname(f))
        print self.current_directory.get()
        self.loadp()
    
    def set_PDBLIST(self):
        f = tkFileDialog.askopenfilename(initialdir=self.current_directory.get(),title='Pick a directory')
        if not f:
            return
        self.PDBLIST.set(f)