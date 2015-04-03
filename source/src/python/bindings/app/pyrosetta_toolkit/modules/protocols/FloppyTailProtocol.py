#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/FloppyTailProtocol.py
## @brief  Class for setting up and running the FloppyTailMover protocol.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.floppy_tail import FloppyTailMover
from rosetta.basic.options import set_integer_option
from rosetta.basic.options import set_boolean_option
from rosetta.basic.options import set_file_vector_option
from rosetta.basic.options import set_string_vector_option
from rosetta.basic.options import set_string_option

#Python Imports
import os
import time

#Tkinter Imports
from Tkinter import *
try:
    import ttk
except ImportError:
    pass

import tkMessageBox
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.Region import *
from ProtocolBaseClass import ProtocolBaseClass
from app.pyrosetta_toolkit.window_main import global_variables

class FloppyTailProtocol(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)

    def setup_protocol(self, main):
        """
        Runs a UI for setting options.
        """

        print "Please reference Kleiger G, Saha A, Lewis S, Kuhlman B, Deshaies RJ. Rapid E2-E3 assembly and disassembly enable processive ubiquitylation of cullin-RING ubiquitin ligase substrates. Cell. 2009 Nov 25;139(5):957-68. PubMed PMID: 19945379."
        print "Please see the rosetta user manual for full documentation of this protocol before use. Other options can be set using the options system."


        if not self.input_class.regions.get_num_regions():
            print "No regions set.  Returning"
            return
        result = tkMessageBox.askokcancel(title="Continue?", message="Production runs require ~30,000 decoys at ~ 1 decoy/hr/processor. Continue?")
        if not result: return

        self.main = Toplevel(main)
        self.main.title("FloppyTail Setup")


        #Set Tk:
        self.set_constraints = IntVar(); self.set_constraints.set(True)
        self.repack_only = IntVar(); self.repack_only.set(True)
        self.c_root = IntVar(); self.c_root.set(False)
        self.linear = IntVar(); self.linear.set(True)
        self.perturb_cycles = IntVar(); self.perturb_cycles.set(5000)
        self.refine_cycles = IntVar(); self.refine_cycles.set(3000)
        self.refine_repack_cycles = IntVar(); self.refine_repack_cycles.set(30)

        self.resfile = StringVar();
        self.frag3 = StringVar();

        self.check_constraints = Checkbutton(self.main, text="centroid_constraints", variable=self.set_constraints)
        self.label_constraints = Label(self.main, text = "Use any set constraints in centroid step")
        self.check_repack = Checkbutton(self.main, text = "-packing:repack_only", variable=self.repack_only)
        self.label_repack = Label(self.main, text = "Do not allow protein sequences to mutate arbitrarily")

        self.check_c_root = Checkbutton(self.main, text = "-FloppyTail:C_root", variable=self.c_root)
        self.label_c_root = Label(self.main, text = "Reroot the fold_tree to the C-terminus. Use if you have an N-Terminal tail")
        self.check_linear = Checkbutton(self.main, text = "-FloppyTail:force_linear_fold_tree", variable = self.linear)
        self.label_linear = Label(self.main, text = "Force a linear fold tree.Use C_root and reordering the chains in your input PDB to ensure correct kinematics")

        self.entry_perturb = Entry(self.main, textvariable=self.perturb_cycles, justify=CENTER)
        self.entry_refine = Entry(self.main, textvariable=self.refine_cycles, justify=CENTER)
        self.entry_refine_repack = Entry(self.main, textvariable=self.refine_repack_cycles, justify=CENTER)
        self.label_perturb = Label(self.main, text = "Perturbation phase cycles (-FloppyTail:perturb_cycles)")
        self.label_refine = Label(self.main, text = "Refinement phase cycles (-FloppyTail:refine_cycles)")
        self.label_refine_repack = Label(self.main, text = "Refinement phase cycles repack every __ cycles (-FloppyTail:refine_repack_cycles)")

        #Python less then 2.7 does not have ttk it looks like
        try:
            self.sep = ttk.Separator(self.main, orient=HORIZONTAL)
        except NameError:
            pass

        self.label_resfile = Label(self.main, text = "Optional Resfile")
        self.label_fragment = Label(self.main, text = "Optional Fragfile (Len3)")
        #self.entry_resfile = Entry(self.main, textvariable=self.resfile)
        #self.entry_frag3 = Entry(self.main, textvariable = self.frag3)
        self.button_resfile = Button(self.main, text="Open", command = lambda:self.get_and_set_resfile())
        self.button_fragment = Button(self.main, text = "Open", command = lambda:self.get_and_set_fragfile())

        self.button_go = Button(self.main, text = "Run Protocol", command = lambda:self.setup_mover_and_run())


        #ShoTk
        self.check_constraints.grid(row=0, column=0, sticky=W); self.label_constraints.grid(row=0, column=1, sticky=W)
        self.check_repack.grid(row=1, column=0, sticky=W); self.label_repack.grid(row=1, column=1, sticky=W)
        self.check_c_root.grid(row=2, column=0, sticky=W); self.label_c_root.grid(row=2, column=1, sticky=W)
        self.check_linear.grid(row=3, column=0, sticky=W); self.label_linear.grid(row=3, column=1, sticky=W)
        self.entry_perturb.grid(row=4, column=0); self.label_perturb.grid(row=4, column=1, sticky=W)
        self.entry_refine.grid(row=5, column=0); self.label_refine.grid(row=5, column=1, sticky=W)
        self.entry_refine_repack.grid(row=6, column=0); self.label_refine_repack.grid(row=6, column=1, sticky=W)
        try:
            self.sep.grid(row=7, column=0, sticky=W+E, pady=4)
        except AttributeError:
            pass

        self.button_resfile.grid(row=8, column=0, sticky=W+E); self.label_resfile.grid(row=8, column=1, sticky=W)
        self.button_fragment.grid(row=9, column=0, sticky=W+E); self.label_fragment.grid(row=9, column=1, sticky=W)

        self.button_go.grid(row=10, column=0, columnspan = 2, sticky=W+E)

    def setup_mover_and_run(self):
        #Option System setters
        if self.linear.get():
            set_boolean_option('FloppyTail:force_linear_fold_tree', True)
        if self.c_root.get():
            set_boolean_option('FloppyTail:C_root', True)
        if (self.set_constraints.get() and self.input_class.constraint_file_path.get()):
            set_string_vector_option('constraints:cst_file', rosetta.Vector1(self.input_class.constraint_file_path.get()))
        if self.resfile.get():
            set_file_vector_option('packing:resfile', rosetta.Vector1(self.resfile.get()))
        if self.frag3.get():
            set_string_option('in:file:frag3', self.frag3.get())
        set_boolean_option('packing:repack_only', self.repack_only.get())
        set_integer_option('FloppyTail:perturb_cycles', self.perturb_cycles.get())
        set_integer_option('FloppyTail:refine_cycles', self.refine_cycles.get())
        set_integer_option('FloppyTail:refine_repack_cycles', self.refine_repack_cycles.get())

        mover = FloppyTailMover()
        if self.pose.total_residue()==0:
            print "Please load a pose."
            return

        if not self.input_class.regions.get_num_regions():
            print "Floppytail requires regions to be set"
            return


        movemap = self.input_class.regions.get_movemap(self.pose)
        movemap.show()
        mover.set_movemap(movemap)
        mover.set_fa_scorefxn(self.score_class.score)
        #Cen Scorefunction - not yet.


        self.run_protocol(mover)
        self.main.destroy()

    def get_and_set_resfile(self):
        self.resfile.set(tkFileDialog.askopenfilename(title="Resfile", initialdir = global_variables.current_directory))
        if not self.resfile.get():
            return
        global_variables.current_directory = os.path.dirname(self.resfile.get())
        #Here we need to set a vector1 string.

    def get_and_set_fragfile(self):
        self.frag3.set(tkFileDialog.askopenfilename(title="Frag3", initialdir = global_variables.current_directory))
        if not self.frag3.get():
            return
        global_variables.current_directory = os.path.dirname(self.frag3.get())
