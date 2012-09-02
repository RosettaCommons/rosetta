#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/simple_analysis.py
## @brief  Handles simple analysis component of the GUI 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
from Tkinter import Frame as TkFrame
from rosetta import *
import tkFileDialog
from modules.tools import analysis as analysis_tools
from modules.tools import input as input_tools
from modules.tools import loops as loop_tools


class SimpleAnalysisFrame(TkFrame):
    def __init__(self, main, toolkit, **options):
        TkFrame.__init__(self, main, **options)
        self.toolkit = toolkit
        self.basicOPT = StringVar()
        self.basicOPT.set("Basic")
        self.basicOPTIONS = {
            "Score FA Energies":lambda:self.print_full_energy(),
            "Score Loops":lambda: self.score_loops(),
            #"Analyze Interface":"",
            #"Score Loop Residues"
            "Number of H-bonds":lambda: self.print_hbonds()

        }
        
        self.rmsdOPT = StringVar()
        self.rmsdOPT.set("RMSD")
        self.rmsdOPTIONS = {
            "RMSD: C alpha":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.loops_as_strings, True),
            "RMSD: BB heavy":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.loops_as_strings, False),
            "RMSD: all atom":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.loops_as_strings, False, True),
            "RMSD: - Load Pose":lambda:self.rmsd_load_pose(),
            "Native RMSD: C alpha":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.loops_as_strings, True),
            "Native RMSD: BB heavy":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.loops_as_strings, False),
            "Native RMSD: all atom":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.loops_as_strings, False, True),
        }
        
    #### Set Tracers ####
        self.basicOPT.trace_variable('w', self.basic_option_tracer)
        self.rmsdOPT.trace_variable('w', self.rmsd_option_tracer)
        
        self.create_GUI_objects()
        self.grid_GUI_objects()
        
    def create_GUI_objects(self):
        self.label_widget=Label(self, text="Basic Analysis", font=("Arial"))
        self.option_menu_basic = OptionMenu(self, self.basicOPT, *(sorted(self.basicOPTIONS)))
        self.option_menu_rmsd = OptionMenu(self, self.rmsdOPT, *(sorted(self.rmsdOPTIONS)))
        self.kickBasic = Button(self, text = "Go", command = lambda: self.kickOptions(self.basicOPT.get()))
        #self.button_Chi = Button(self, text="Chi"); self.button_Chi2 = Button(self, text="Write to file")
    
    def grid_GUI_objects(self):
        """
        Columnspan: 2
        Rowspan: 2
        """
        r=0; c=0;
        self.label_widget.grid(row=r, column=c, columnspan=2,)
        self.option_menu_basic.grid(row=r+1, column=c, sticky=W+E); self.option_menu_rmsd.grid(row = r+1, column=c+1, sticky=W+E)
        
        
    #### TRACERS ####
    def basic_option_tracer(self, name, index, mode):
        varValue = self.basicOPT.get()
        func = self.basicOPTIONS[self.basicOPT.get()]
        func()
        return
    
    def rmsd_option_tracer(self, name, index, mode):
        varValue = self.rmsdOPT.get()
        func = self.rmsdOPTIONS[self.rmsdOPT.get()]
        func()
        return
    
    #### FUNCTIONS ####
    
    def rmsd_load_pose(self):
        self.rmsdP = Pose()
        f = tkFileDialog.askopenfilename(initialdir=self.toolkit.current_directory.get(), title='Choose PDB')
        pose_from_pdb(self.rmsdP, f)
        print "Pose to compare loaded..."
        return
      
    def score_loops(self):
        for l in self.toolkit.loops_as_strings:
            print "\n"+l
            residue_array = loop_tools.return_residue_array(self.toolkit.pose, l)
            start_chain = self.toolkit.pose.pdb_info().pose2pdb(residue_array[0])
            start = start_chain.split()[0]; chain = start_chain.split()[1]
            end = self.toolkit.pose.pdb_info().pose2pdb(residue_array[-1]).split()[0]
            e = self.toolkit.input_class.ScoreBaseObject.ret_loop_energy(chain, int(start), int(end))
            print "Total Region energy: %.3f REU"%e
            
            n_e = self.toolkit.input_class.ScoreBaseObject.ret_loop_neighbor_energy(chain, int(start), int(end))
            print "Total Neighbor (ci_2b) Interaction energy: %.3f REU"%n_e
    
    def print_full_energy(self):
        self.toolkit.ScoreObject.score.show(self.toolkit.pose)
        
    def print_hbonds(self):
        n = self.toolkit.input_class.ScoreBaseObject.ret_n_hbonds()
        print "Num Hbonds: "+repr(n)

            

