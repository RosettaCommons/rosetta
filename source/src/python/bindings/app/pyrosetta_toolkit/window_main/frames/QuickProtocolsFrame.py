#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/frames/quick_min.py
## @brief  Handles the quick protocols section of the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Tkinter Imports
from Tkinter import *
from Tkinter import Frame as TkFrame

#Toolkit Imports
from app.pyrosetta_toolkit.modules.protocols.MinimizationProtocols import MinimizationProtocols
from app.pyrosetta_toolkit.modules.tools import output as output_tools

#from app.pyrosetta_toolkit import main_window
from app.pyrosetta_toolkit.window_main.IO.GUIOutput import GUIOutput

class QuickProtocolsFrame(TkFrame):
    def __init__(self, main, toolkit, output_class, **options):
        self.toolkit = toolkit
        self.output_class = output_class
        TkFrame.__init__(self, main, **options)
        
        self.min_cmd = StringVar()
        self.min_cmd.set("Optimize Rotamers") ; #Default command
        
        self.randomize_loops_bool=StringVar(); #Random Loop check_button_ckbox
        self.set_options_menus()
        
        #Will change when protocols are re-organized.
        self.min_protocols = MinimizationProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.output_class)

        self.create_GUI_objects()
        self.grid_GUI_objects()
        
        #Ignore this. It is for Komodo autocomplete.
        if 0:
            self.main = Tk()
            #self.toolkit = main_window()
            self.output_class = GUIOutput()
            
    def set_options_menus(self):
        
        self.minOPTIONS = {
            "Optimize Rotamers":lambda: self.min_protocols.optimize_rotamers(),
            "Optimize Rotamers (SCWRL)":lambda: self.min_protocols.SCWRL(),
            "Minimize Backbone+Sidechains":lambda: self.min_protocols.minimize(),
            "Minimize Backbone only":lambda:self.min_protocols.minimize(False, False, True),
            "Minimize Sidechains only":lambda:self.min_protocols.minimize(False, False, False, True),
            "FastRelax Backbone+Sidechains":lambda: self.min_protocols.relax(1),
            "FastRelax Backbone only":lambda:self.min_protocols.relax(1, False, True),
            "FastRelax Sidechains only":lambda:self.min_protocols.relax(1, False, False, True),
        }
        
    def create_GUI_objects(self):
        
        #self.label_quick_min = Label(self, text="Quick Protocols", font=("Arial"))
        self.randomize_button=Checkbutton(self, text="Randomize Loops?", variable=self.randomize_loops_bool, command=lambda: self.toolkit.pose.assign(tools.loops.initLoops().RandLoop(self.toolkit.pose, self.toolkit.input_class.loops_as_strings, self.output_class.show_pymol.get())))
        
        self.decoy_entry=Entry(self, textvariable=self.output_class.decoys, justify=CENTER)
        self.rounds_entry=Entry(self, textvariable=self.output_class.rounds, relief=SUNKEN, justify=CENTER)
        self.processors_entry = Entry(self, textvariable=self.output_class.processors, justify=CENTER, relief=SUNKEN)
        #self.rounds_entry.set(1)
        self.general_label = Label(self, text = "General Protocol Options")
        self.decoy_label=Label(self, text="Decoys")
        self.rounds_label=Label(self, text="Rounds")
        self.processors_label = Label(self, text="Processors")
        self.min_options = OptionMenu(self, self.min_cmd, *(sorted(self.minOPTIONS)))
        #self.full_min_options =  OptionMenu(self, self.fullCommand, *(sorted(self.FullMinOPTIONS)))
        self.kickmin_options = Button(self, text = "Run Protocol", command = lambda: self.kick_min_cmd())
        
    def grid_GUI_objects(self):
        """
        Rowspan : 5
        Columnspan : 2
        """
        r = 0; c=0;
        #self.label_quick_min.grid(row=r, column=c, columnspan=2, pady=15)
        self.general_label.grid(row=r+2, column=c, columnspan=2, sticky=W+E)
        self.processors_entry.grid(row=r+3, column=c); self.processors_label.grid(row=r+3, column=c+1, sticky=W)
        self.decoy_entry.grid(row=r+4, column=c); self.decoy_label.grid(row=r+4, column=c+1, sticky=W)
        self.rounds_entry.grid(row=r+5, column=c); self.rounds_label.grid(row=r+5, column=c+1, sticky=W)
        
        #Minimization
        self.min_options.grid(row=r+1, column = c, sticky=W+E); 
        self.kickmin_options.grid(row=r+1, column=c+1, sticky=W+E);

    def kick_min_cmd(self):
        """
        This kicks the minimization of the loop.
        """
        func = self.minOPTIONS[self.min_cmd.get()]
        func()
        return
    
            
    
