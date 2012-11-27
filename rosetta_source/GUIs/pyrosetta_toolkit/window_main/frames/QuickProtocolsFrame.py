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
from modules.protocols.loop_minimization import Loop_Min
from modules.protocols.protein_minimization import Protein_Min
from modules.tools import output as output_tools



class QuickProtocolsFrame(TkFrame):
    def __init__(self, main, toolkit, output_class, **options):
        self.toolkit = toolkit
        self.output_class = output_class
        TkFrame.__init__(self, main, **options)
        
        self.loopCommand = StringVar()
        self.loopCommand.set("Relax Loop(s)") ; #Default loop command
        
        self.fullCommand = StringVar()
        self.fullCommand.set("Relax All"); #Default full command
        
        self.fast_relax_bool=StringVar(); self.fast_relax_bool.set(1)#Fast Relax check_button_ckbox
        self.randomize_loops_bool=StringVar(); #Random Loop check_button_ckbox
        self.set_options_menus()
        
        #Will change when protocols are re-organized.
        self.loop_protocols = Loop_Min(self.toolkit.score_class, self.toolkit.pose)
        self.full_protocols = Protein_Min(self.toolkit.score_class, self.toolkit.pose)

        self.create_GUI_objects()
        self.grid_GUI_objects()
        
    def set_options_menus(self):
        
        self.LoopMinOPTIONS = {
            "Relax Loop(s)":lambda: self.loop_protocols.RelaxLoop(self.rounds_entry.get(), self.toolkit.input_class.loops_as_strings, self.fast_relax_bool.get()),
            #"Backrub Loop(s)":self.loop_protocols.LoopBackrubRef(self.rounds_entry.get(), self.toolkit.input_class.loops_as_strings),
            "Minimize Loop(s)":lambda: self.loop_protocols.classicMinLoop(self.rounds_entry.get(), self.toolkit.input_class.loops_as_strings),
            "Optimize Loop Rotamers":lambda: self.loop_protocols.optimizeRotLoop(self.rounds_entry.get(), self.toolkit.input_class.loops_as_strings),
            "Optimize Loop Rotamers (SCWRL)":lambda: self.loop_protocols.SCWRL(self.toolkit.input_class.loops_as_strings, int(self.rounds_entry.get()) )
        }
        
        self.FullMinOPTIONS = {
            "Relax All":lambda: self.full_protocols.Relax(self.rounds_entry.get(), self.fast_relax_bool.get()),
            #"Backrub All"lambda: self.full_protocols.Backrub(self.rounds_entry.get())
            "Minimize All":lambda: self.full_protocols.classicMin(self.rounds_entry.get()),
            "Optimize All Rotamers":lambda: self.full_protocols.optimizeRot(self.rounds_entry.get()),
            "Optimize All Rotamers (SCWRL)":lambda: self.full_protocols.SCWRL(int(self.rounds_entry.get()))
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
        self.loops_min_options = OptionMenu(self, self.loopCommand, *(sorted(self.LoopMinOPTIONS)))
        self.full_min_options =  OptionMenu(self, self.fullCommand, *(sorted(self.FullMinOPTIONS)))
        self.kickloops_min_options = Button(self, text = "Minimize Region(s)", command = lambda: self.kickMinimizationLoop())
        self.kickfull_min_options =  Button(self, text = "Minimize PDB", command = lambda: self.kickMinimizationFull())
        
    def grid_GUI_objects(self):
        """
        Rowspan : 5
        Columnspan : 2
        """
        r = 0; c=0;
        #self.label_quick_min.grid(row=r, column=c, columnspan=2, pady=15)
        self.decoy_entry.grid(row=r+3, column=c); self.decoy_label.grid(row=r+3, column=c+1)
        self.rounds_entry.grid(row=r+4, column=c, sticky=W+E); self.rounds_label.grid(row=r+4, column=c+1)
        self.fastrelax_button.grid(row=r+5, column=c, columnspan = 2)
        
        
        #Minimization
        self.loops_min_options.grid(row=r+1, column = c, sticky=W+E); self.full_min_options.grid(row=r+2, column=c, sticky=W+E)
        self.kickloops_min_options.grid(row=r+1, column=c+1, sticky=W+E); self.kickfull_min_options.grid(row=r+2, column=c+1, sticky=W+E)

    def kickMinimizationLoop(self):
        """
        This kicks the minimization of the loop.
        """
        
        if self.output_class.auto_write.get():
            jd=PyJobDistributor(self.output_class.outdir.get() + "/" + self.output_class.outname.get(), int(self.decoy_entry.get()), self.toolkit.score_class.score);
            jd.native_pose = self.toolkit.native_pose
            for i in range(0, int(self.decoy_entry.get())):
                print "Decoy: "+repr(i)
                func = self.LoopMinOPTIONS[self.loopCommand.get()]
                func()
                jd.output_decoy(self.toolkit.pose)
                
        else:
            func = self.LoopMinOPTIONS[self.loopCommand.get()]
            func()
        return
    
    def kickMinimizationFull(self):
        """
        This kicks the minimization of the protein.
        """
        if self.output_class.auto_write.get():
            jd=PyJobDistributor(self.output_class.outdir.get() + "/" + self.output_class.outname.get(), int(self.decoy_entry.get()), self.toolkit.score_class.score);
            jd.native_pose = self.toolkit.native_pose
            for i in range(0, int(self.decoy_entry.get())):
                print "Decoy: "+repr(i)
                func = self.FullMinOPTIONS[self.fullCommand.get()]
                func()
                jd.output_decoy(self.toolkit.pose)
                
        else:
            func = self.FullMinOPTIONS[self.fullCommand.get()]
            func()
        return
            
    
