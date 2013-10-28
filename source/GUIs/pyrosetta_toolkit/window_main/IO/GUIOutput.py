
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
import tkSimpleDialog

#Toolkit Imports
from window_main import global_variables

class GUIOutput:
    def __init__(self, toolkit):
        """
        Clas responsible for managing output variables of the GUI.
        """
        self.toolkit = toolkit;
        self.processors = IntVar(); self.processors.set(1); # We will only use more if the user sets more to use.
        self.use_boltzmann = IntVar(); self.use_boltzmann.set(False)
        self.recover_low = IntVar(); self.recover_low.set(False)
        self.kT = 1.0; #Main KT used for rounds
        self.outdir = StringVar(); self.outdir.set(global_variables.current_directory)
        self.outname = StringVar();
        self.decoys=IntVar(); self.decoys.set(0)
        self.rounds=IntVar(); self.rounds.set(1)
        self.VarLoopSeq = StringVar(); #Sequence of the particular region just queried.
        self.auto_write = IntVar(); self.auto_write.set(False) ; #Auto Save pose after protocol?
        self.overwrite = IntVar(); self.overwrite.set(False); #Overwrite PDBs?
        self.terminal_output = IntVar(); #Tracer to redirect stdout or not.  Tracer code is in PyRosetta Toolkit.py 0 is textbox, 1 is stdout.

        
        
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
            
    
    
######### Functions that cannot be put in output_tools, as they set a variable within this class. ################
    def set_temperature(self):
        self.kT = tkSimpleDialog.askfloat(title="kT", prompt="Set Temperature", initialvalue=self.kT)
        
        
        
