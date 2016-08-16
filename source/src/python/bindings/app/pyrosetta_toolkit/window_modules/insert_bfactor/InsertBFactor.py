#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/InsertBFactor.py
## @brief  Simple window for inserting arbitrary data from a file to BFactor column of PDBs by residue or atom
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#Python Imports
import os
import re
import copy

#Tkinter Imports
from Tkinter import *
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.window_main import global_variables
from app.pyrosetta_toolkit.modules.PythonPDB import PythonPDB

class InsertBFactor:
    """
    A simple window for inserting data from CSVs, tab deliminated, space delimated files into the B-factor column of a PDB for residues or atoms.
    """

    def __init__(self, input_class):
        self.filename = input_class.pdb_path.get()

        if not self.filename:
            print "Please load a pose."
            return

        self.PythonPDB = PythonPDB(self.filename)
        self.input_class = input_class

        self.resnum_column = StringVar(); self.resnum_column.set("1")
        self.chain_column = StringVar(); self.chain_column.set("2")
        self.data_column = StringVar(); self.data_column.set("3")
        self.atomname_column = StringVar(); self.atomname_column.set("NA")

        self.delim_type = StringVar(); self.delim_type.set(" ")


    def show_window(self, main):
        if not self.filename:
            print "Please load a pose."
            return

        self.main = Toplevel(main)
        self.main.title("Insert data into B Factor")

        self.resnum_label = Label(self.main, text="Resnum\nColumn")
        self.resnum_entry = Entry(self.main, textvariable=self.resnum_column, justify=CENTER)
        self.chain_label = Label(self.main, text ="Chain\nColumn")
        self.chain_entry = Entry(self.main, textvariable=self.chain_column, justify=CENTER)
        self.data_label = Label(self.main, text = "Data\nColumn")
        self.data_entry = Entry(self.main, textvariable=self.data_column, justify=CENTER)
        self.atomname_label = Label(self.main, text = "Atom Name Column\n(optional)")
        self.atomname_entry = Entry(self.main, textvariable=self.atomname_column, justify=CENTER)

        self.space_delim = Radiobutton(self.main, text = "space", variable = self.delim_type, value=" ")
        self.tab_delim = Radiobutton(self.main, text = "tab", variable = self.delim_type, value="\t")
        self.comma_delim = Radiobutton(self.main, text = "comma", variable = self.delim_type, value = ",")
        self.delim_label = Label(self.main, text = "Deliminator: ")

        self.insert_button = Button(self.main, text="Insert", command = lambda:self.insert())
        self.insert_load_button = Button(self.main, text = "Insert and Load PDB", command = lambda: self.insert_and_reload())

        self.resnum_label.grid(row=0, column=0); self.chain_label.grid(row=0, column=1)
        self.data_label.grid(row=0, column=2); self.atomname_label.grid(row=0, column=3)

        self.resnum_entry.grid(row=1, column=0); self.chain_entry.grid(row=1, column=1)
        self.data_entry.grid(row=1, column=2); self.atomname_entry.grid(row=1, column=3)

        self.delim_label.grid(row=2, column=0, sticky=W)
        self.space_delim.grid(row=2, column=1, sticky=W)
        self.tab_delim.grid(row=2, column=2, sticky=W)
        self.comma_delim.grid(row=2, column=3, sticky=W)

        self.insert_button.grid(row=3, column=0, columnspan=2, sticky=W+E)
        self.insert_load_button.grid(row=3, column=2, columnspan=2, sticky=W+E)

### Functions ###

    def insert(self):
        filename = tkFileDialog.askopenfilename(title="Data file", initialdir=global_variables.current_directory)
        if not filename:return
        global_variables.current_directory = os.path.dirname(filename)

        if self.atomname_column.get()=="NA":
            self.PythonPDB.read_file_and_replace_b_factors(self.delim_type.get(), filename, \
                int(self.resnum_column.get()), int(self.chain_column.get()), int(self.data_column.get()))
        else:
            self.PythonPDB.read_file_and_replace_b_factors(self.delim_type.get(), filename, \
                int(self.resnum_column.get()), int(self.chain_column.get()), int(self.data_column.get()), int(self.atomname_column.get()))

        outname = tkFileDialog.asksaveasfilename(title="Output PDB Name", initialdir=global_variables.current_directory)
        if not outname:return
        if not re.search("pdb", outname):
            outname = outname+".pdb"

        self.PythonPDB.save_PDB(outname)
        print "PDB bfactor insertion complete.  Model saved."
        self.main.destroy()

    def insert_and_reload(self):
        filename = tkFileDialog.askopenfilename(title="Data file", initialdir=global_variables.current_directory)
        if not filename:return
        global_variables.current_directory = os.path.dirname(filename)

        if self.atomname_column.get()=="NA":
            self.PythonPDB.read_file_and_replace_b_factors(self.delim_type.get(), filename, \
                int(self.resnum_column.get()), int(self.chain_column.get()), int(self.data_column.get()))
        else:
            self.PythonPDB.read_file_and_replace_b_factors(self.delim_type.get(), filename, \
                int(self.resnum_column.get()), int(self.chain_column.get()), int(self.data_column.get()), int(self.atomname_column.get()))

        outname = tkFileDialog.asksaveasfilename(title="Output PDB Name", initialdir=global_variables.current_directory)
        if not outname:return
        if not re.search("pdb", outname):
            outname = outname+".pdb"

        self.PythonPDB.save_PDB(outname)
        print "PDB bfactor insertion complete.  Model saved."
        self.input_class.load_pose(outname)
        self.main.destroy()

