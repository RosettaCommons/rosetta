#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/GraftingProtocols.py
## @brief  Functions for Running GraftMovers and integrated GUI for setup
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
import rosetta.protocols.grafting as graft

#Python Imports
import re
import os

#Tkinter Imports
from Tkinter import *
import tkSimpleDialog
import tkFileDialog

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
from app.pyrosetta_toolkit.window_main import global_variables

class GraftingProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
        pass


class GraftMoverWindow(ProtocolBaseClass):
    """
    Protocol for running the GraftMovers.  Setup is done through a GUI window.
    """

    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
        self.insert_start = StringVar(); self.insert_end = StringVar(); self.insert_chain = StringVar();
        self.scaffold_start = input_class.region_start; self.scaffold_end = input_class.region_end; self.scaffold_chain = input_class.region_chain;
        self.insert_nter_flex = IntVar(); self.insert_cter_flex = IntVar();
        self.scaffold_nter_flex = IntVar(); self.scaffold_cter_flex = IntVar();
        self.scaffold_nter_flex.set(2); self.scaffold_cter_flex.set(2)

        #Grafting Options
        self.cycles = IntVar(); self.cycles.set(200)
        self.graft_type = StringVar()
        self.graft_types = ["Default", "Double Arm", "Double Loop Double Arm", "Double Loop Quad Arm"]

        self.repack_connection = IntVar(); self.repack_connection.set(True)
        self.repack_connection_and_piece = IntVar(); self.repack_connection_and_piece.set(False)

        self.randomize_first = IntVar(); self.randomize_first.set(False)
        self.smooth_centroid = IntVar(); self.smooth_centroid.set(False)


    def setTk(self, main):
        """
        Setup TK objects
        """
        main.title("Grafting Setup")
        self.main = main

        self.button_load_from_pose = Button(self.main, text = "Load Pose", command=lambda:self.load_from_pose())
        self.label_current_pose = Label(self.main, text = "Current Pose")
        self.label_insert = Label(self.main, text = "INSERT")
        self.label_scaffold = Label(self.main, text = "BETWEEN")
        self.label_start = Label(self.main, text = "Start")
        self.label_end = Label(self.main, text = "End")
        self.label_chain = Label(self.main, text = "Chain")
        self.label_nter_flex = Label(self.main, text = "Nter Flex")
        self.label_cter_flex = Label(self.main, text = "Cter Flex")

        #Entries
        self.entry_insert_start = Entry(self.main, textvariable=self.insert_start)
        self.entry_insert_start.config(width=4)
        self.entry_insert_end = Entry(self.main, textvariable=self.insert_end)
        self.entry_insert_end.config(width=4)
        self.entry_insert_chain = Entry(self.main, textvariable = self.insert_chain)
        self.entry_insert_chain.config(width=4)
        self.entry_insert_nter_flex = Entry(self.main, textvariable = self.insert_nter_flex)
        self.entry_insert_nter_flex.config(width=3)
        self.entry_insert_cter_flex = Entry(self.main, textvariable = self.insert_cter_flex)
        self.entry_insert_cter_flex.config(width=3)
        self.entry_scaffold_start = Entry(self.main, textvariable=self.scaffold_start)
        self.entry_scaffold_start.config(width=4)
        self.entry_scaffold_end = Entry(self.main, textvariable = self.scaffold_end)
        self.entry_scaffold_end.config(width=4)
        self.entry_scaffold_chain = Entry(self.main, textvariable = self.scaffold_chain)
        self.entry_scaffold_chain.config(width=4)
        self.entry_scaffold_nter_flex = Entry(self.main, textvariable = self.scaffold_nter_flex)
        self.entry_scaffold_nter_flex.config(width=3)
        self.entry_scaffold_cter_flex = Entry(self.main, textvariable = self.scaffold_cter_flex)
        self.entry_scaffold_cter_flex.config(width = 3)

        self.button_begin = Button(self.main, text = "Run Protocol", command = lambda:self.run_graft())

        #Options
        self.label_cycles = Label(self.main, text = "Cycles")
        self.entry_cycles = Entry(self.main, textvariable = self.cycles)
        self.entry_cycles.config(width=7)
        self.label_graft_method = Label(self.main, text = "Method")
        self.optionmenu_method = OptionMenu(self.main, self.graft_type, *self.graft_types)
        self.optionmenu_method.config(width=7)
        self.checkbutton_randomize = Checkbutton(self.main, text="Randomize flexible residues before protocol?", variable = self.randomize_first)
        self.checkbutton_smooth = Checkbutton(self.main, text="Use smooth centroid e-function for centroid mode?", variable = self.smooth_centroid)
        self.checkbutton_repack_con = Checkbutton(self.main, text = "Repack flexible residues after protocol?", variable = self.repack_connection)
        self.checkbutonn_repack_con_and_insert = Checkbutton(self.main, text = "Repack flexible residues + insert?", variable = self.repack_connection_and_piece)

    def shoTk(self, r, c):
        """
        Show Tk Objects
        """
        self.label_start.grid(row=r, column=c+2); self.label_end.grid(row=r, column=c+3); self.label_chain.grid(row=r, column=c+4); self.label_nter_flex.grid(row=r, column=c+5); self.label_cter_flex.grid(row=r, column=c+6)

        self.button_load_from_pose.grid(row=r+1, column=c, sticky=W+E); self.label_insert.grid(row=r+1, column=c+1, sticky=E);
        self.entry_insert_start.grid(row=r+1, column=c+2); self.entry_insert_end.grid(row=r+1, column=c+3); self.entry_insert_chain.grid(row=r+1, column=c+4); self.entry_insert_nter_flex.grid(row=r+1, column=c+5); self.entry_insert_cter_flex.grid(row=r+1, column=c+6)

        self.label_current_pose.grid(row=r+2, column=c, sticky=W+E); self.label_scaffold.grid(row=r+2, column=c+1, sticky=E)
        self.entry_scaffold_start.grid(row=r+2, column=c+2); self.entry_scaffold_end.grid(row=r+2, column=c+3); self.entry_scaffold_chain.grid(row=r+2, column=c+4); self.entry_scaffold_nter_flex.grid(row=r+2, column=c+5); self.entry_scaffold_cter_flex.grid(row=r+2, column=c+6)

        self.label_cycles.grid(row=r+3, column=c+2, sticky=E); self.entry_cycles.grid(row=r+3, column=c+3, columnspan=2, sticky=E+W),self.checkbutton_randomize.grid(row=r+3, column=c+5, columnspan=5, sticky=W)
        self.label_graft_method.grid(row=r+4, column=c+2, sticky=E); self.optionmenu_method.grid(row=r+4, column=c+3, columnspan=2, sticky=E+W), self.checkbutton_smooth.grid(row=r+4, column=c+5, columnspan=5, sticky=W)


        self.checkbutton_repack_con.grid(row=r+7, column=c+5, columnspan=5, sticky=W)
        self.checkbutonn_repack_con_and_insert.grid(row=r+8, column=c+5, columnspan=5, sticky=W)

        self.button_begin.grid(row=r+9, column=c+1, columnspan=6, sticky=W+E)


    ###Callbacks###
        print "Please Reference: \n Lewis SM, Kuhlman BA. Anchored design of protein-protein interfaces. PLoS One. 2011;6(6):e20872. Epub 2011 Jun 17.\n Gulyani A, Vitriol E, Allen R, Wu J, Gremyachinskiy D, Lewis S, Dewar B, Graves LM, Kay BK, Kuhlman B, Elston T, Hahn KM. A biosensor generated via high-throughput screening quantifies cell edge Src dynamics. Nat Chem Biol. 2011 Jun 12;7(7):437-44. doi: 10.1038/nchembio.585."
        self.graft_type.trace_variable('w', self.graft_method_callback)
        self.graft_type.set(self.graft_types[0])
        self.repack_connection.trace_variable('w', self.repack_connection_callback)
        self.repack_connection_and_piece.trace_variable('w', self.repack_connection_and_piece_callback)


    def graft_method_callback(self, name, index, mode):
        varValue = self.graft_type.get()
        if varValue == "Default":
            print "Using Default method as used in AnchoredPDBCreator application"
            print "This uses 1 CCD arm that includes the insert.  Insert will move in cartesian space to close the loop."
            print "****Nter_loop_start---->Piece----> | Cter_loop_end****"

        elif varValue =="Double Arm":
            print "This uses 2 CCD arms that includes the insert.  Insert will move in cartesion space to close the loop."
            print "****Nter_loop_start---->Piece | <----Nter_loop_end****"

        elif varValue=="Double Loop Double Arm":
            print "Using Two CCD arms and two loops.  Insert is frozen in Cartesian space"
            print " ****Nter_loop_start-----> | Piece | <----Nter_loop_end****"
        elif varValue=="Double Loop Quad Arm":
            print "Using 4 CCD arms - 2 arms for each side of the insert, as in normal loop modeling.  Insert is frozen in Cartesian space"
            print "****Nter_loop_start----> | <---- Piece -----> | <------ Cter_loop_end****"

    def repack_connection_callback(self, name, index, mode):
        if self.repack_connection.get():
            self.repack_connection_and_piece.set(False)
        else:
            self.repack_connection_and_piece.set(True)

    def repack_connection_and_piece_callback(self, name, index, mode):
        if self.repack_connection_and_piece.get():
            self.repack_connection.set(False)
        else:
            self.repack_connection.set(True)

    ###Functions###
    def load_from_pose(self):
        filename = tkFileDialog.askopenfilename(title="PDB with insert OR insert PDB", initialdir=global_variables.current_directory)
        if not filename:return
        global_variables.current_directory = os.path.dirname(filename)
        self.from_pose = self.input_class.return_loaded_pose(filename)
        print self.from_pose
        print "Donor pose loaded..."

    def run_graft(self):
        anchor_start = self.pose.pdb_info().pdb2pose(self.scaffold_chain.get(), int(self.scaffold_start.get()))
        anchor_end = self.pose.pdb_info().pdb2pose(self.scaffold_chain.get(), int(self.scaffold_end.get()))

        graftmover = graft.AnchoredGraftMover(anchor_start, anchor_end)

        insert_start = self.from_pose.pdb_info().pdb2pose(self.insert_chain.get(), int(self.insert_start.get()))
        insert_end = self.from_pose.pdb_info().pdb2pose(self.insert_chain.get(), int(self.insert_end.get()))

        if re.search("Double Loop", self.graft_type.get()):
            nter_overhang = tkSimpleDialog.askinteger(title='N-terminal Overhang', prompt='Please enter the number of residues N-terminal from insert_start to use for superposition', initialvalue=3)
            if not nter_overhang:return
            cter_overhang = tkSimpleDialog.askinteger(title='C-terminal Overhang', prompt='Please enter the number of residues C-terminal from insert_end to use for superposition', initialvalue=3)
            insert = graft.return_region(self.from_pose, insert_start-nter_overhang, insert_end+cter_overhang)
            if not cter_overhang:return
            print "Insert with overhang:"
            print insert
            graftmover.set_piece(insert, nter_overhang, cter_overhang)
            graftmover.superimpose_overhangs_heavy(self.pose, False, False)
        else:
            insert = graft.return_region(self.from_pose, insert_start, insert_end)
            print "Insert: "
            print insert
            graftmover.set_piece(insert, 0, 0)

        if self.graft_type.get()=="Double Arm":
            graftmover.set_use_single_loop_double_CCD_arms(True)
        elif self.graft_type.get()=="Double Loop Double Arm":
            graftmover.set_use_double_loop_double_CCD_arms(True)
        elif self.graft_type.get()=="Double Loop Quad Arm":
            graftmover.set_use_double_loop_quad_CCD_arms(True)
        else:
            print "Using default graft type"

        if self.randomize_first.get():
            graftmover.set_test_control_mode(True)

        graftmover.set_cycles(self.cycles.get())
        graftmover.set_scaffold_flexibility(self.scaffold_nter_flex.get(), self.scaffold_cter_flex.get())
        graftmover.set_insert_flexibility(self.insert_nter_flex.get(), self.insert_cter_flex.get())

        new_graftmover = GraftMover(graftmover)
        new_graftmover.set_repack_connection(self.repack_connection.get(), self.score_class.score)
        new_graftmover.set_repack_connection_and_insert(self.repack_connection_and_piece.get(), self.score_class.score)

        self.run_protocol(new_graftmover)
        self.main.destroy()

class GraftMover:
    """
    Small class to repack during apply method - Passed to run_mover in base class.  Workaround till/if we update the C++ class
    """

    def __init__(self, grafter):
        self.grafter = grafter
        self.repack_connection = False
        self.repack_connection_and_insert = False

    def set_repack_connection(self, boolean, fa_scorefunction):
        self.repack_connection = boolean
        self.fa_scorefxn = fa_scorefunction

    def set_repack_connection_and_insert(self, boolean, fa_scorefunction):
        self.repack_connection_and_insert=boolean
        self.fa_scorefxn = fa_scorefunction

    def apply(self, pose):
        self.grafter.apply(pose)
        if self.repack_connection:
            self.grafter.repack_connection_and_residues_in_movemap(pose, self.fa_scorefxn)
        elif self.repack_connection_and_insert:
            self.grafter.repack_connection_and_residues_in_movemap_and_piece(pose, self.fa_scorefxn)


