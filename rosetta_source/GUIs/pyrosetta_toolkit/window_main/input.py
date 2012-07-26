#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/input.py
## @brief  Handles the main input component of the GUI 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
from modules import *
from modules.ScoreBase import ScoreBase as score_base
from modules.tools import input as input_tools
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
import os

class Input():
    
    def __init__(self, main, toolkit):
        self.main = main
        self.toolkit = toolkit
        self.current_directory = StringVar()
        self.indir = StringVar()
        self.filename = StringVar(); self.filename.set("0");
        self.PDBLIST = StringVar(); self.PDBLIST.set("")
        self.loops_as_strings = []; #Array of Loops: start:end:chain
        self.VarStartLoop=StringVar(); #Start of Loop
        self.VarEndLoop=StringVar(); #End of Loop
        self.VarLoopChain=StringVar(); #Chain of Loop
        self.VarLoopSeq=StringVar(); #Sequence in Entry
        self.loops = Loops()
        
    def setTk(self):
        self.label_display = Label(self.main, text="PyRosetta Toolkit July 2012", font=("Arial", 14))
        self.label_A=Label(self.main, text="Choose PDBLIST:") #This needs to change to a load file or path widget       
        self.label_B=Label(self.main, text="Input Options", font=("Arial"))     
        self.button_PDBin = Button(self.main, text="Find", command=lambda: self.set_PDBLIST())   
        self.label_BB= Label(self.main, text="Choose PDB:")
        self.button_PDBin2= Button(self.main, text="Find", command=lambda: self.choose_load_pose())
        #self.button_loadD = Button(self.main, text="Load Folder")
        
        ############## LOOPS ##############
        self.label_Loop = Label(self.main, text="Region Selection", font=("Arial"))
        self.AddLoopButton = Button(self.main, text = "Add Region", command=lambda: self.addLoop())
        self.StartLoopLabel=Label(self.main, text="Start of Region:")
        self.EndLoopLabel=Label(self.main, text="End of Region:")
        self.ChainIDLabel=Label(self.main, text="Chain ID:")
        self.loops_listbox = Listbox(self.main)
        self.StartLoopEntry=Entry(self.main, textvariable=self.VarStartLoop)
        self.EndLoopEntry=Entry(self.main, textvariable=self.VarEndLoop)
        self.ChainIDEntry=Entry(self.main, textvariable=self.VarLoopChain)
        ##################################
        
        
    def shoTk(self, r=0, c=0):
        self.label_display.grid(row=r, column=c,columnspan=6, pady=10)
        self.label_A.grid(row=2, column=0, sticky=W)
        self.label_B.grid(row=1, column=0, columnspan=2)
        self.button_PDBin.grid(row=2, column=1, sticky=W+E)
        self.label_BB.grid(row=3, column=0, sticky=W); self.button_PDBin2.grid(row=3, column=1, sticky=W+E)
        #self.button_loadD.grid(row=4, column=1, sticky=W+E); Was to load a folder...
        
        ############# LOOPS #############
        self.label_Loop.grid(row=11, column=0, columnspan=2, pady=15)
        self.loops_listbox.bind("<Double-Button-1>", lambda event: self.remLoop())
        self.StartLoopLabel.grid(row=16, column=0); self.StartLoopEntry.grid(row=16, column=1)
        self.EndLoopLabel.grid(row=17, column=0); self.EndLoopEntry.grid(row=17, column=1)
        self.ChainIDLabel.grid(row=18, column=0); self.ChainIDEntry.grid(row=18, column=1)
        self.loops_listbox.grid(row=20, column=0, rowspan=6, columnspan = 1, sticky = W+E); self.AddLoopButton.grid(row=19, column=1, sticky = W+E)
        ################################
        
#### POSE ####

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
        '''
        Loads a File through the tk File Dialog
        '''
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
#### LOOPS ####

    def addLoop(self):
        looFull = self.VarStartLoop.get()+ ":"+ self.VarEndLoop.get()+":"+self.VarLoopChain.get().upper()
        self.loops_as_strings.append(looFull)
        self.loops_listbox.insert(END, looFull)
        #print self.toolkit.input_class.loops_as_strings
    def remLoop(self):

        self.loops_as_strings.remove(self.loops_listbox.get(self.loops_listbox.curselection()))
        self.loops_listbox.delete(self.loops_listbox.curselection())
        #print self.toolkit.input_class.loops_as_strings
        
    