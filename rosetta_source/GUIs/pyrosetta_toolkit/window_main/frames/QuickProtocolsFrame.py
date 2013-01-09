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
from modules.protocols.MinimizationProtocols import MinimizationProtocols
from modules.tools import output as output_tools



class QuickProtocolsFrame(TkFrame):
    def __init__(self, main, toolkit, output_class, **options):
        self.toolkit = toolkit
        self.output_class = output_class
        TkFrame.__init__(self, main, **options)
        
        self.min_cmd = StringVar()
        self.min_cmd.set("Optimize Rotamers") ; #Default loop command
        
        self.fast_relax_bool=StringVar(); self.fast_relax_bool.set(1)#Fast Relax check_button_ckbox
        self.randomize_loops_bool=StringVar(); #Random Loop check_button_ckbox
        self.set_options_menus()
        
        #Will change when protocols are re-organized.
        self.min_protocols = MinimizationProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.output_class)

        self.create_GUI_objects()
        self.grid_GUI_objects()
        
    def set_options_menus(self):
        
        self.minOPTIONS = {
            "Optimize Rotamers":lambda: self.min_protocols.optimize_rotamers(),
            "Optimize Rotamers (SCWRL)":lambda: self.min_protocols.SCWRL(),
            "Minimize Structure":lambda: self.min_protocols.minimize(),
            "Minimize BB":lambda:self.min_protocols.minimize(False, False, True),
            "Minimize SC":lambda:self.min_protocols.minimize(False, False, False, True),
            "Relax Structure":lambda: self.min_protocols.relax(self.fast_relax_bool.get()),
            "Relax BB":lambda:self.min_protocols.relax(self.fast_relax_bool.get(), False, True),
            "Relax SC":lambda:self.min_protocols.relax(self.fast_relax_bool.get(), False, False, True),

        }
        
    def create_GUI_objects(self):
        
        #self.label_quick_min = Label(self, text="Quick Protocols", font=("Arial"))
        self.randomize_button=Checkbutton(self, text="Randomize Loops?", variable=self.randomize_loops_bool, command=lambda: self.toolkit.pose.assign(tools.loops.initLoops().RandLoop(self.toolkit.pose, self.toolkit.input_class.loops_as_strings, self.output_class.show_pymol.get())))
        self.fastrelax_button=Checkbutton(self, text="Fast Relax?", variable=self.fast_relax_bool)
        
        self.decoy_entry=Entry(self, textvariable=self.output_class.decoys, justify=CENTER)
        self.rounds_entry=Entry(self, textvariable=self.output_class.rounds, relief=SUNKEN, justify=CENTER)
        #self.rounds_entry.set(1)
        self.decoy_label=Label(self, text="Decoys Desired")
        self.rounds_label=Label(self, text="Rounds Desired")
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
        self.decoy_entry.grid(row=r+2, column=c); self.decoy_label.grid(row=r+2, column=c+1)
        self.rounds_entry.grid(row=r+3, column=c, sticky=W+E); self.rounds_label.grid(row=r+3, column=c+1)
        self.fastrelax_button.grid(row=r+4, column=c, columnspan = 2)
        
        
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
    
            
    
