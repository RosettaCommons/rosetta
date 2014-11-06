#!/usr/bin/python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/pymol_integration/PyMOL.py
## @brief  Main window and controller of PyMol integration accross the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import glob
import os

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog

#Tookit Imports
#from app.pyrosetta_toolkit.window_modules.scorefunction.ScoreFxnControl import ScoreFxnControl

class AdvancedPyMOL():
    """
    Pymol Visualization Window
    Version 2.0 ONLY
    """

    def __init__(self, pose):
	"""
	This class handles all PyMOL integration with the GUI.  
	"""

	self.auto_send = IntVar()
	self.auto_send.set(False)
        
        self.keep_history=IntVar()
        self.keep_history.set(True)
        
        self.send_energies = IntVar()
        self.send_energies.set(False)
        
	self.pymover = PyMOL_Mover()
	self.pymover.keep_history(True)
	self.observer = PyMOL_Observer()
	self.observer.pymol.keep_history(self.keep_history.get())
	self.observer.pymol.update_energy(self.send_energies.get())
	
        #Tracers for PyMOL Observer.
	self.auto_send.trace_variable('w', self.auto_change_observer)
        self.keep_history.trace_variable('w', self.auto_change_keep_observer_history)
        self.send_energies.trace_variable('w', self.auto_change_send_observer_energies)
        
        self.auto_send_region_colors = IntVar(); self.auto_send_region_colors.set(False)
        self.auto_send_residue_colors = IntVar(); self.auto_send_residue_colors.set(False)
        
	self.pnum = 0
	self.pose = pose
	self.send_label = IntVar();
	self.send_label.set(0)
	
	#Pymol Functions currently in trunk.
	self.pymol_functions = {
	"Send Residue Energy":lambda:self.pymover.send_energy(self.pose, 'total_score', self.send_label.get()),
	"View Hydrogen Bonds": lambda:self.pymover.send_hbonds(self.pose),
	"View Polar Identity": lambda:self.pymover.send_polars(self.pose),
	"View DSSP SS": lambda:self.pymover.send_polars(self.pose),
	"View Foldtree":lambda: self.pymover.send_foldtree(self.pose),
	"View Foldtree Diagram":lambda: self.pymover.view_foldtree_diagram(self.pose)
	#"View MoveMap": lambda: self.pymover.send_movemap(self.pose)
	}
        
    def auto_change_observer(self, name, index, mode):
	varValue = self.auto_send.get()
	if not self.auto_send.get():
	    print "Removed observer"
	    self.observer.remove_observer(self.pose)
	    
	else:
	    print "Added observer"
	    self.observer.add_observer(self.pose)
    
    def auto_change_keep_observer_history(self, name, index, mode):
        varValue = self.keep_history.get()
        if not self.keep_history.get():
            print "Not keeping history"
            self.observer.pymol.keep_history(False)
        else:
            print "Keeping history"
            self.observer.pymol.keep_history(True)
    
    def auto_change_send_observer_energies(self, name, index, mode):
        varValue = self.send_energies.get()
        if not self.send_energies.get():
            print "Not sending energies"
            self.observer.pymol.update_energy(False)
        else:
            print "Sending energies. Slightly slower updating in PyMOL."
            self.observer.pymol.update_energy(True)
            
    def makeWindow(self, row, column, main, ScoreObject):

	self.main = main
	self.main.title("PyMOL Visualization")
	try :
	    print self.pose.pdb_info().name()
	except AttributeError:
	    tkMessageBox.showwarning(message = 'Please Load a Pose...')
	    return
	    
	self.ScoreObject = ScoreObject; self.score = self.ScoreObject.score
	
	#For Komodo Autocomplete only:
	#if 0:
	    #self.main = Tk()
	    #self.ScoreObject = ScoreFxnControl()
	    
	self.viewlabel_ = Label(self.main, text="View Options")
	self.scorelabel_ = Label(self.main, text="Score Terms")
	self.viewList = Listbox(self.main)
        self.view_scroll = Scrollbar(self.main)
        self.viewList.config(yscrollcommand=self.view_scroll.set); self.view_scroll.config(command = self.viewList.yview)
        
	self.scoreList = Listbox(self.main)
        self.score_scroll = Scrollbar(self.main)
        self.scoreList.config(yscrollcommand = self.score_scroll.set); self.score_scroll.config(command = self.scoreList.yview)
        
	self.sendbutton_ = Button(self.main, text = "Send Pose", command = lambda: self.SendPose(self.viewList.get(self.viewList.curselection())))
	self.sendnewbutton_ = Button(self.main, text = "Send New Pose", command = lambda: self.SendNewPose())
	self.autosend_new_checkbutton = Checkbutton(self.main, text = "Send Poses as new Objects", variable = self.auto_send)
	self.label_energies_checkbutton = Checkbutton(self.main, text = "Label Energies", variable = self.send_label)
	
	### Grid ###
        self.label_energies_checkbutton.grid(row=row, column=column, sticky=W)
        self.autosend_new_checkbutton.grid(row = row+1, column=column, sticky=W)
	self.viewlabel_.grid(row = row+2, column=column, padx=5); 
	self.viewList.grid(row = row+3, column=column, padx=5);self.view_scroll.grid(row=row+3, column=column+1, sticky=N+S);
        self.scorelabel_.grid(row = row+4, column = column)
        self.scoreList.grid(row = row+5, column = column);self.score_scroll.grid(row=row+5, column=column+1, sticky=N+S)
        
        
	#self.sendbutton_.grid(row = row+4, column = column+1, sticky=W+E); self.sendnewbutton_.grid(row = row+3, column = column+1, sticky=W+E); 
	
	
	for option in sorted(self.pymol_functions):
	    self.viewList.insert(END, option)

	ZeroTerms, NonZeroTerms = ScoreObject.scoreOption("Breakdown ScoreFxn")
	for option in NonZeroTerms:
	    self.scoreList.insert(END, option)

	self.viewList.bind("<Double-Button-1>", lambda event: self.SendPose(self.viewList.get(self.viewList.curselection())))
	self.scoreList.bind("<Double-Button-1>", lambda event: self.SendEPose(self.scoreList.get(self.scoreList.curselection())))
    
#### Pose Sending #### 
    def SendPose(self, option):
	"""
	Sends pose according to option. Will Fix as callback.
	"""
	print self.score(self.pose)
	
	func = self.pymol_functions[option]
	try:
	    func()
	    return
	except AttributeError:
	    print "Function not working in this version of PyRosetta...."
	    return
	
    def SendEPose(self, option):
	optionSP = option.split(";")
	term = optionSP[0].lstrip()
	print term
	#vars()[term]
	print self.score(self.pose)
	
	try:
	    self.pymover.send_energy(self.pose, term, self.send_label.get())
	except IOError:
	    print "Could not send label for the individual energy function."
	
    def SendNewPose(self):
	"""
	Changes the name of the pose so that the next thing sent is a new object...
	"""
	self.pnum = self.pnum+1
	Newname = os.path.basename(self.pose.pdb_info().name()).split("_")[0]+"_"+repr(self.pnum)
	self.pose.pdb_info().name(Newname)
        try:
            self.pymover.apply(self.pose)
        except PyRosettaException:
            print "Could not send pose to pymol."
            return
	return
    
#### Coloring Functions ####
    def color_residue(self, rosetta_resnum, color="magenta", all_other="blue"):
        """
        Colors residue using PyMol Mover send energies.  I don't know how to *not* color the other residues.
        """
        color_map = dict()
        color_map[rosetta_resnum]=color
        self.pymover.send_colors(self.pose, color_map, all_other)
    
    def color_region(self, resnum_start, resnum_end, color="magenta", all_other="blue"):
        color_map = dict()
        for i in range(resnum_start, resnum_end+1):
            color_map[i]=color
        
        self.pymover.send_colors(self.pose, color_map, all_other)
        
    def color_regions(self, regions, color="magenta", all_other="blue"):
        movemap = regions.get_movemap(self.pose)
        self.pymover.send_movemap(self.pose, movemap, all_other)
