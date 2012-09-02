
#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/output.py
## @brief  Handles the main output component of the GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
from modules.tools import input as input_tools
from modules.tools import output as output_tools
from modules.tools import sequence as sequence_tools

class OutputFrame(Frame):
    def __init__(self, main, toolkit, **options):
        #Initialize frame, which becomes master instead of main - Which is why we use self in making the GUI components
        Frame.__init__(self, main, **options)
        
        self.outdir = toolkit.input_class.current_directory
        self.outname = StringVar(); self.outname.set("0")
        self.toolkit = toolkit; #Main window and everything in it is global.  Not great coding, but it works fairly well...
        self.decoys=IntVar(); self.decoys.set(0)
        self.VarLoopSeq = StringVar(); #Sequence of the particular region just queried.
        self.auto_write = IntVar(); self.auto_write.set(0) ; #Auto Save pose after protocol?
        
    #### Tracers ####
        self.auto_write.trace_variable('w', self.auto_set_decoys)
        self.decoys.trace_variable('w', self.auto_set_write)
        
        self.create_GUI_objects()
        self.grid_GUI_objects()
        
    def create_GUI_objects(self):
        self.label_output_options = Label(self, text="Output Options", font=("Arial"))
        self.label_output_path = Label(self, text="Output Path:")
        self.label_output_name = Label(self, text="Output Name:")
        self.button_Out=Button(self, text="Choose", command=lambda: self.outdir.set(input_tools.tk_get_directory(self.outdir.get())))
        self.entry_Outname = Entry(self, textvariable = self.outname)
        self.check_button_auto_write = Checkbutton(self, text="Write after protocol?", variable = self.auto_write)
        self.button_dump_pose = Button(self, text="Write Current Pose", command=lambda: output_tools.dumpPDB(self.toolkit.pose, self.outdir.get() + "/" + self.outname.get(), self.toolkit.ScoreObject.score))


        
    def grid_GUI_objects(self):
        self.label_output_options.grid(row=0, column=0, columnspan=2, sticky=W+E)
        self.label_output_path.grid(row=1, column=0, sticky=W)
        self.label_output_name.grid(row=2, column=0, sticky=W)
        self.button_Out.grid(row=1, column=1, sticky=W+E)
        self.entry_Outname.grid(row=2, column=1)
        self.button_dump_pose.grid(row=3, column=1)
        
        self.check_button_auto_write.grid(row=3, column=0)
    #### Tracers ####
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