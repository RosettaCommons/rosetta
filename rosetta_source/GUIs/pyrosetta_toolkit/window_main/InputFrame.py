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

class InputFrame(Frame):
    
    def __init__(self, main, toolkit, **options):
        self.toolkit = toolkit
        #Initialize frame, which becomes master instead of main - Which is why we use self in making the GUI components
        Frame.__init__(self, main, **options)
        self.pwd= self.location()[0]
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
        
        self.create_GUI_objects()
        self.grid_GUI_objects()

        
    def create_GUI_objects(self):
        self.label_pdblist=Label(self, text="Choose PDBLIST:") #This needs to change to a load file or path widget       
        self.main_label=Label(self, text="Input Options", font=("Arial"))     
        self.button_PDBList = Button(self, text="Find", command=lambda: self.set_PDBLIST())   
        self.label_pdb= Label(self, text="Choose PDB:")
        self.button_PDB= Button(self, text="Find", command=lambda: self.choose_load_pose())
        #self.button_loadD = Button(self.main, text="Load Folder")
        
        ############## LOOPS ##############
        self.label_Loop = Label(self, text="Region Selection", font=("Arial"))
        self.AddLoopButton = Button(self, text = "Add Region", command=lambda: self.addLoop())
        self.StartLoopLabel=Label(self, text="Start of Region:")
        self.EndLoopLabel=Label(self, text="End of Region:")
        self.ChainIDLabel=Label(self, text="Chain ID:")
        self.loops_listbox = Listbox(self)
        self.StartLoopEntry=Entry(self, textvariable=self.VarStartLoop)
        self.EndLoopEntry=Entry(self, textvariable=self.VarEndLoop)
        self.ChainIDEntry=Entry(self, textvariable=self.VarLoopChain)
        ##################################
        
        
        ############ SEQUENCE ############
        self.ShoSeqButton = Button(self, text="Show Sequence", command=lambda: self.VarLoopSeq.set(sequence_tools.get_sequence(self.toolkit.pose, self.VarStartLoop.get()+":"+self.VarEndLoop.get()+":"+self.VarLoopChain.get())))
        self.entry_LoopSeq = Entry(self, textvariable=self.VarLoopSeq, justify=CENTER)
        ##################################
        
        
        ############ PYMOL ###############
        self.button_ShowPymol = Button(self, text="Show in PyMol", command=lambda: self.toolkit.PyMOLObject.pymover.apply(self.toolkit.pose)) #Whichever of these is picked, it should hold it in it's memory and be able to do whatever with them?
        self.check_button_Pym = Checkbutton(self, text = "Pymol Observer?", variable=self.toolkit.PyMOLObject.auto_send)
        ##################################
        
           ### Photo ###
      
        DesignPhoto =PhotoImage(file = (self.pwd+ "/media/RosettaLogo.gif"))
        self.Photo = Label(master=self, image=DesignPhoto)
        self.Photo.image = DesignPhoto
      
    def grid_GUI_objects(self):
        self.label_pdblist.grid(row=2, column=0, sticky=W)
        self.main_label.grid(row=1, column=0, columnspan=2)
        self.button_PDBList.grid(row=2, column=1, sticky=W+E)
        self.label_pdb.grid(row=3, column=0, sticky=W); self.button_PDB.grid(row=3, column=1, sticky=W+E)
        
        ############# LOOPS #############
        self.label_Loop.grid(row=11, column=0, columnspan=2, pady=15)
        self.loops_listbox.bind("<Double-Button-1>", lambda event: self.remLoop())
        self.StartLoopLabel.grid(row=16, column=0); self.StartLoopEntry.grid(row=16, column=1)
        self.EndLoopLabel.grid(row=17, column=0); self.EndLoopEntry.grid(row=17, column=1)
        self.ChainIDLabel.grid(row=18, column=0); self.ChainIDEntry.grid(row=18, column=1)
        self.loops_listbox.grid(row=20, column=0, rowspan=6, columnspan = 1, sticky = W+E, padx=3); self.AddLoopButton.grid(row=19, column=1, sticky = W+E)
        ################################
        
        
        ########### SEQUENCE ###########
        self.ShoSeqButton.grid(row=12, column=0)
        self.entry_LoopSeq.grid(row=12, column=1)
        ################################
        
        
        ########### PYMOL ##############
        self.check_button_Pym.grid(row=5, column=0)
        self.button_ShowPymol.grid(row=5, column=1)

        self.Photo.grid(row=20, column=1, rowspan=6, columnspan=1, sticky=W+E, padx=3)
        
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
        
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported/ran from anywhere.
        """
            
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP