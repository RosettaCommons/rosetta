
#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/IO/GUIOutput.py
## @brief  Class responsible for managing output variables of the GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Tkinter Imports
from Tkinter import *

#Toolkit Imports
from window_main import global_variables

class GUIOutput:
    def __init__(self, toolkit):
        """
        Clas responsible for managing output variables of the GUI.
        """
        
        self.toolkit = toolkit;
        self.outdir = StringVar(); self.outdir.set(global_variables.current_directory)
        self.outname = StringVar(); self.outname.set("0")
        self.decoys=IntVar(); self.decoys.set(0)
        self.rounds=IntVar(); self.rounds.set(1)
        self.VarLoopSeq = StringVar(); #Sequence of the particular region just queried.
        self.auto_write = IntVar(); self.auto_write.set(0) ; #Auto Save pose after protocol?

        #### Tracers ####
        self.auto_write.trace_variable('w', self.auto_set_decoys)
        self.decoys.trace_variable('w', self.auto_set_write)
        
    def auto_set_decoys(self, name, index, mode):
        """
        Changes decoy number according to auto_write variable
        """
        
        varValue = self.auto_write.get()
        if varValue and not self.decoys.get():
            self.decoys.set(1)
        elif not varValue:
            self.decoys.set(0)
        return
    
    def auto_set_write(self, name, index, mode):
        """
        Changes auto_write according to decoy_number
        """
        try:
            varValue = self.decoys.get()
        except ValueError:
            pass
            return
        if varValue and not self.auto_write.get():
            self.auto_write.set(1)
        elif not varValue:
            self.auto_write.set(0)