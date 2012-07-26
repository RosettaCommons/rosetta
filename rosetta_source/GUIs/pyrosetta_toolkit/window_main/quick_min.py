#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/quick_min.py
## @brief  Handles the quick protocols section of the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from modules.protocols.loop_minimization import Loop_Min
from modules.protocols.protein_minimization import Protein_Min
from modules.tools import output as output_tools
from Tkinter import *
from rosetta import *

class Minimization():
    def __init__(self, main, toolkit):
        self.main = main
        self.toolkit = toolkit
        self.loopCommand = StringVar()
        self.loopCommand.set("Relax Loop(s)") ; #Default loop command
        
        self.fullCommand = StringVar()
        self.fullCommand.set("Relax All"); #Default full command
        
        self.fast_relax_bool=StringVar(); self.fast_relax_bool.set(1)#Fast Relax check_button_ckbox
        self.randomize_loops_bool=StringVar(); #Random Loop check_button_ckbox
        self.set_options_menus()
        
        #Will change when protocols are re-organized.
        self.loop_protocols = Loop_Min(self.toolkit.ScoreObject, self.toolkit.pose)
        self.full_protocols = Protein_Min(self.toolkit.ScoreObject, self.toolkit.pose)


    def set_options_menus(self):
        
        self.LoopMinOPTIONS = {
            "Relax Loop(s)":lambda: self.loop_protocols.RelaxLoop(self.SliHigh.get(), self.toolkit.input_class.loops_as_strings, self.fast_relax_bool.get()),
            #"Backrub Loop(s)":self.loop_protocols.LoopBackrubRef(self.SliHigh.get(), self.toolkit.input_class.loops_as_strings),
            "Minimize Loop(s)":lambda: self.loop_protocols.classicMinLoop(self.SliHigh.get(), self.toolkit.input_class.loops_as_strings),
            "Optimize Loop Rotamers":lambda: self.loop_protocols.optimizeRotLoop(self.SliHigh.get(), self.toolkit.input_class.loops_as_strings),
            "Optimize Loop Rotamers (SCWRL)":lambda: self.loop_protocols.SCWRL(self.toolkit.input_class.loops_as_strings, int(self.SliHigh.get()) )
        }
        
        self.FullMinOPTIONS = {
            "Relax All":lambda: self.full_protocols.Relax(self.SliHigh.get(), self.fast_relax_bool.get()),
            #"Backrub All"lambda: self.full_protocols.Backrub(self.SliHigh.get())
            "Minimize All":lambda: self.full_protocols.classicMin(self.SliHigh.get()),
            "Optimize All Rotamers":lambda: self.full_protocols.optimizeRot(self.SliHigh.get()),
            "Optimize All Rotamers (SCWRL)":lambda: self.full_protocols.SCWRL(int(self.SliHigh.get()))
        }
        
    def setTk(self):
        
        #self.button_Save = Button(self.main, text = "Save Pose")
        self.label_quick_min = Label(self.main, text="Quick Protocols", font=("Arial"))
        self.check_button_Rand=Checkbutton(self.main, text="Randomize Loops?", variable=self.randomize_loops_bool, command=lambda: self.toolkit.pose.assign(tools.loops.initLoops().RandLoop(self.toolkit.pose, self.toolkit.input_class.loops_as_strings, self.toolkit.output_class.show_pymol.get())))
        self.check_button_Relax=Checkbutton(self.main, text="Fast Relax?", variable=self.fast_relax_bool)
        self.button_LoopW=Button(self.main, text="Write New PDB")
        
        ##### CDRS ######
        self.button_L1=Button(self.main, text="CDR L1", command=lambda: setL1()); self.button_L2=Button(self.main, text="CDR L2", command=lambda: setL2()); self.button_L3=Button(self.main, text="CDR L3", command=lambda: setL3())
        self.button_H1=Button(self.main, text="CDR H1", command=lambda: setH1()); self.button_H2=Button(self.main, text="CDR H2", command=lambda: setH2()); self.button_H3=Button(self.main, text="CDR H3", command=lambda: setH3())
        ##### CDRS ######
        
        self.entry_decoy=Entry(self.main, textvariable=self.toolkit.output_class.decoys, justify=CENTER)
        self.SliHigh=Scale(self.main, orient=HORIZONTAL, to=100, resolution=1, from_=1, relief=SUNKEN)
        self.SliHigh.set(1)
        self.label_Dec=Label(self.main, text="Decoys Desired")
        self.label_High=Label(self.main, text="Rounds Desired")
        self.loops_min_options = OptionMenu(self.main, self.loopCommand, *(sorted(self.LoopMinOPTIONS)))
        self.full_min_options =  OptionMenu(self.main, self.fullCommand, *(sorted(self.FullMinOPTIONS)))
        self.kickloops_min_options = Button(self.main, text = "Minimize Region(s)", command = lambda: self.kickMinimizationLoop())
        self.kickfull_min_options =  Button(self.main, text = "Minimize PDB", command = lambda: self.kickMinimizationFull())
        
    def shoTk(self, r=0, c=0):
        '''
        Rowspan : 5
        Columnspan : 2
        '''
        
        #self.button_L1.grid(row=13, column=0, sticky=W+E); self.button_H1.grid(row=13, column=1, sticky=W+E)
        #self.button_L2.grid(row=14, column=0, sticky=W+E); self.button_H2.grid(row=14, column=1, sticky=W+E)
        #self.button_L3.grid(row=15, column=0, sticky=W+E); self.button_H3.grid(row=15, column=1, sticky=W+E)
        self.label_quick_min.grid(row=r, column=c, columnspan=2, pady=15)
        self.entry_decoy.grid(row=r+3, column=c); self.label_Dec.grid(row=r+3, column=c+1)
        self.SliHigh.grid(row=r+4, column=c, sticky=W+E); self.label_High.grid(row=r+4, column=c+1)
        self.check_button_Relax.grid(row=r+5, column=c, columnspan = 2)
        
        
        #Minimization
        self.loops_min_options.grid(row=r+1, column = c, sticky=W+E); self.full_min_options.grid(row=r+2, column=c, sticky=W+E)
        self.kickloops_min_options.grid(row=r+1, column=c+1, sticky=W+E); self.kickfull_min_options.grid(row=r+2, column=c+1, sticky=W+E)

    def kickMinimizationLoop(self):
        '''
        This kicks the minimization of the loop.
        '''
        
        if self.toolkit.output_class.auto_write.get():
            jd=PyJobDistributor(self.toolkit.outdir.get() + "/" + self.toolkit.outname.get(), int(self.entry_decoy.get()), self.toolkit.ScoreObject.score);
            jd.native_pose = self.toolkit.native_pose
            for i in range(0, int(self.entry_decoy.get())):
                print "Decoy: "+repr(i)
                func = self.LoopMinOPTIONS[self.loopCommand.get()]
                func()
                jd.output_decoy(self.toolkit.pose)
                
        else:
            func = self.LoopMinOPTIONS[self.loopCommand.get()]
            func()
        return
    
    def kickMinimizationFull(self):
        '''
        This kicks the minimization of the protein.
        '''
        if self.toolkit.output_class.auto_write.get():
            jd=PyJobDistributor(self.toolkit.outdir.get() + "/" + self.toolkit.outname.get(), int(self.entry_decoy.get()), self.toolkit.ScoreObject.score);
            jd.native_pose = self.toolkit.native_pose
            for i in range(0, int(self.entry_decoy.get())):
                print "Decoy: "+repr(i)
                func = self.FullMinOPTIONS[self.fullCommand.get()]
                func()
                jd.output_decoy(self.toolkit.pose)
                
        else:
            func = self.FullMinOPTIONS[self.fullCommand.get()]
            func()
        return
            
    
