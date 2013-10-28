#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/vicinity_options/vicinity_window.py
## @brief  Vicinity dialog window
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkSimpleDialog
from os import getcwd

pwd = getcwd()
class vicinityrelaxwindow(tkSimpleDialog.Dialog):
    """
    This is the Vaccinity Relax Window to specify options while doing relax.  May be used in other protocols if possible.
    It is set as a tkSimpleDialog, but this can be changed fairly easily to allow more use.
    """
    
    def body(self, main):
        #self.main = Toplevel(main)
        self.main = main
        #self.main.title("Neighbor Options")
        #self.column = column; #Column to begin Grid (Int)
        #self.row = row; #Row to begin Grid (Int)
        #self.pwd = pwd
        row = 0; column= 0
        #Options:
        self.FixTarget = StringVar(); self.FixTarget.set("Open")
        self.FixVaccinity = StringVar(); self.FixVaccinity.set("Open")
        self.FixBoth = StringVar(); self.FixBoth.set("UnSet")
        
        
        #Set Tk
        print pwd
        #Photo
        VacPhoto =PhotoImage(file = (pwd+"/Media/Vaccinity_Smaller.gif"))
        self.Photo = Label(self.main, image=VacPhoto)
        self.Photo.image = VacPhoto
        
        #Button/Labels
        #self.setOptionsbutton_ = Button(self.main, text = "Continue...", command = lambda: self.setOptions())
        self.control = Label(self.main, text = "Control.")
        
        #Fix Options
        self.FixBBlab = Label(self.main, text = "   Loop/Target   ")
        self.FixBBOpt = OptionMenu(self.main, self.FixTarget, "Fix", "Fix BB", "Fix Chi", "Open")
        self.FixChilabel_ = Label(self.main, text = "   Vaccinity   ")
        self.FixChiOpt = OptionMenu(self.main, self.FixVaccinity, "Fix", "Fix BB", "Fix Chi", "Open")
        self.FixBotlabel_ = Label(self.main, text = "Fix BackBone and Rotamers")
        self.FixBothOpt = OptionMenu(self.main, self.FixBoth, "UnSet", "Target", "Vaccinity")
        
        
        #ShoTk

        self.Photo.grid(row =row, column = column+1, rowspan=17, columnspan=17)
        self.FixBBlab.grid(row = row+5, column= column)
        self.FixBBOpt.grid(row = row+6, column = column, sticky="ew")
        
        self.FixChilabel_.grid(row=row+5, column = column+18)
        self.FixChiOpt.grid(row=row+6, column = column+18, sticky="ew")
        
        self.FixBotlabel_.grid(row=row, column = column+8+1)
        self.FixBothOpt.grid(row = row+1, column = column+8+1)
        
        self.control.grid(row=row+18, column=column+8+1)
        #self.setOptionsbutton_.grid(row = row +19, column = column+8+1)
    def apply(self):
        fixloo = self.FixTarget.get()
        fixvac = self.FixVaccinity.get()
        fixboth = self.FixBoth.get()
        self.result = (fixloo, fixvac, fixboth)
        #self.main.destroy()
        
