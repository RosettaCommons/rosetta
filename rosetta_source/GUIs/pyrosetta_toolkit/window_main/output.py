
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

class Outputs():
    def __init__(self, main, toolkit):
        self.main = main
        self.outdir = toolkit.input_class.current_directory
        self.outname = StringVar(); self.outname.set("0")
        self.toolkit = toolkit; #Main window and everything in it is global.  Not great coding, but it works fairly well...
        self.decoys=IntVar(); self.decoys.set(0)
        self.VarLoopSeq = StringVar(); #Sequence of the particular region just queried.
        self.auto_write = IntVar(); self.auto_write.set(0) ; #Auto Save pose after protocol?
        
    #### Tracers ####
        self.auto_write.trace_variable('w', self.auto_set_decoys)
        self.decoys.trace_variable('w', self.auto_set_write)
        
    def setTk(self):
        self.label_output_options = Label(self.main, text="Output Options", font=("Arial"))
        self.label_output_path = Label(self.main, text="Output Path:")
        self.label_output_name = Label(self.main, text="Output Name:")
        self.button_Out=Button(self.main, text="Choose", command=lambda: self.outdir.set(input_tools.tk_get_directory(self.outdir.get())))
        self.entry_Outname = Entry(self.main, textvariable = self.outname)
        self.check_button_auto_write = Checkbutton(self.main, text="Write after protocol?", variable = self.auto_write)
        self.button_dump_pose = Button(self.main, text="Write Current Pose(s)", command=lambda: output_tools.dumpPDB(self.toolkit.pose, self.outdir.get() + "/" + self.outname.get(), self.toolkit.ScoreObject.score))
        #button_GcdrFW = Button(a, text="Write CDR and Frameworks") 
        self.button_ShowPymol = Button(self.main, text="Show in PyMol", command=lambda: self.toolkit.PyMOLObject.pymover.apply(self.toolkit.pose)) #Whichever of these is picked, it should hold it in it's memory and be able to do whatever with them?
        self.check_button_Pym = Checkbutton(self.main, text = "Pymol Observer?", variable=self.toolkit.PyMOLObject.auto_send)
        
        ########### Sequence ##########
        self.ShoSeqButton = Button(self.main, text="Show Sequence", command=lambda: self.VarLoopSeq.set(sequence_tools.get_sequence(self.toolkit.pose, self.toolkit.input_class.VarStartLoop.get()+":"+self.toolkit.input_class.VarEndLoop.get()+":"+self.toolkit.input_class.VarLoopChain.get())))
        self.entry_LoopSeq = Entry(self.main, textvariable=self.VarLoopSeq, justify=CENTER)
    def shoTk(self):
        self.label_output_options.grid(row=1, column=2, columnspan=2, sticky=W+E)
        self.label_output_path.grid(row=2, column=2, sticky=W)
        self.label_output_name.grid(row=3, column=2, sticky=W)
        self.button_Out.grid(row=2, column=3, sticky=W+E)
        self.entry_Outname.grid(row=3, column=3)
        self.button_dump_pose.grid(row=5, column=3)
        #button_GcdrFW.grid(row=5, column=2, columnspan=2, sticky=W+E)
        self.button_ShowPymol.grid(row=5, column=1)
        self.check_button_Pym.grid(row=5, column=0)
        self.check_button_auto_write.grid(row=5, column=2)
        ######### Sequence ###########
        self.ShoSeqButton.grid(row=12, column=0)
        self.entry_LoopSeq.grid(row=12, column=1)
        
    #### Tracers ####
    def auto_set_decoys(self, name, index, mode):
        '''
        Changes decoy number according to auto_write variable
        '''
        
        varValue = self.auto_write.get()
        if varValue and not self.decoys.get():
            self.decoys.set(1)
        elif not varValue:
            self.decoys.set(0)
        return
    
    def auto_set_write(self, name, index, mode):
        '''
        Changes auto_write according to decoy_number
        '''
        try:
            varValue = self.decoys.get()
        except ValueError:
            pass
            return
        if varValue and not self.auto_write.get():
            self.auto_write.set(1)
        elif not varValue:
            self.auto_write.set(0)