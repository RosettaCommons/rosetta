#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/frames/simple_analysis.py
## @brief  Handles simple analysis component of the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.scoring import ScoreType
from rosetta.core.scoring.hbonds import *

#Python Imports
import os

#Tkinter Imports
from Tkinter import *
from Tkinter import Frame as TkFrame
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import analysis as analysis_tools
from app.pyrosetta_toolkit.modules.tools import input as input_tools
from app.pyrosetta_toolkit.modules.tools import loops as loop_tools
from app.pyrosetta_toolkit.window_main import global_variables
#from pyrosetta_toolkit import main_window

class SimpleAnalysisFrame(TkFrame):
    def __init__(self, main, toolkit, **options):
        TkFrame.__init__(self, main, **options)
        self.toolkit = toolkit
        self.basicOPT = StringVar()
        self.basicOPT.set("Basic Analysis")
        self.basicOPTIONS = {
            "Score FA Energies":lambda:self.print_full_energy(),
            #"Score Loops":lambda: self.score_loops(),
            "Score Overall":lambda:self.print_energy(),
            #"Score Loop Residues"
            "Print Pose info":lambda:self.print_pose(),
            "Print Option info":lambda:self.print_option_info(),
            "Print Hbond info":lambda:self.print_all_hbonds(),
            "Total Hydrogen Bonds":lambda: self.print_hbonds(),

        }

        self.rmsdOPT = StringVar()
        self.rmsdOPT.set("RMSD Analysis")
        self.rmsdOPTIONS = {
            "RMSD: C alpha":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, True),
            "RMSD: BB heavy":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, False),
            "RMSD: all atom":lambda:analysis_tools.rmsd(self.rmsdP, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, False, True),
            "RMSD: - Load Pose":lambda:self.rmsd_load_pose(),
            "Native RMSD: C alpha":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, True),
            "Native RMSD: BB heavy":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, False),
            "Native RMSD: all atom":lambda:analysis_tools.rmsd(self.toolkit.native_pose, self.toolkit.pose, self.toolkit.input_class.loops_as_strings, False, True),
        }

    #### Set Tracers ####
        self.basicOPT.trace_variable('w', self.basic_option_tracer)
        self.rmsdOPT.trace_variable('w', self.rmsd_option_tracer)

        self.create_GUI_objects()
        self.grid_GUI_objects()

        #Ignore this.  It is for Komodo autocomplete.
        #if 0:
            #self.main = Tk()
            #self.toolkit = main_window()

    def create_GUI_objects(self):
        #self.label_widget=Label(self, text="Basic Analysis", font=("Arial"))
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
        #self.label_widget.grid(row=r, column=c, columnspan=2,)
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

    def rmsd_load_pose(self, infilename = False):
        self.rmsdP = Pose()
        if not infilename:
            infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title='Choose PDB')
        if not infilename:return
        global_variables.current_directory = os.path.dirname(infilename)
        pose_from_file(self.rmsdP, infilename)
        print "Pose to compare loaded..."
        return

    def score_loops(self):
        for l in self.toolkit.input_class.loops_as_strings:
            print "\n"+l
            residue_array = loop_tools.return_residue_array(self.toolkit.pose, l)
            start_chain = self.toolkit.pose.pdb_info().pose2pdb(residue_array[0])
            start = start_chain.split()[0]; chain = start_chain.split()[1]
            end = self.toolkit.pose.pdb_info().pose2pdb(residue_array[-1]).split()[0]
            e = self.toolkit.input_class.regional_score_class.ret_loop_energy(chain, int(start), int(end))
            print "Total Region energy: %.3f REU"%e

            n_e = self.toolkit.input_class.regional_score_class.ret_loop_neighbor_energy(chain, int(start), int(end))
            print "Total Neighbor (ci_2b) Interaction energy: %.3f REU"%n_e

    def print_full_energy(self):
        if self.toolkit.pose.total_residue()==0:return
        """
        if self.toolkit.output_class.terminal_output.get():
            self.toolkit.score_class.score.show(self.toolkit.pose)
        else:
            self.show_energies()
        """
        self.toolkit.score_class.score.show(self.toolkit.pose)
        print "Printed to Terminal"
        #self.toolkit.score_class.score.show(self.toolkit.TR, self.toolkit.pose); No Go.  sys.stdout also does not work.

    def show_energies(self):
        """
        Currently, string inputted into textbox is uneven.  Not using till this is fixed.
        """
        if self.toolkit.pose.total_residue()==0:return
        out = self.toolkit.input_class.regional_score_class.ret_energy_string()
        print out

    def print_energy(self):
        if self.toolkit.pose.total_residue()==0:return

        score = self.toolkit.score_class.score(self.toolkit.pose)
        print "%.3f REU"%score

    def print_pose(self):
        print self.toolkit.pose

    def print_option_info(self):
        self.toolkit.input_class.options_manager.print_current_options()

    def print_hbonds(self):
        if self.toolkit.pose.total_residue()==0:return
        n = self.toolkit.input_class.regional_score_class.ret_n_hbonds()
        print "Detected Hydrogen bonds: "+repr(n)

    def print_all_hbonds(self):
        if self.toolkit.pose.total_residue()==0:return
        set = HBondSet();
        self.toolkit.pose.update_residue_neighbors()
        fill_hbond_set(self.toolkit.pose, True, set);
        set.show(self.toolkit.pose)

