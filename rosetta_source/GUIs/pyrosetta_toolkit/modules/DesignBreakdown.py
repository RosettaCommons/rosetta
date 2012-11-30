#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/DesignBreakdown.py
## @brief  Class for analyzing design results from a fasta of sequences
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#Tkinter Imports
from Tkinter import *

#Python Imports
import re

#Toolkit Imports
from prettytable import prettytable
from PDB import *

class DesignBreakdown:
    """
    This class functions in organizing breakdowns for design.  Perhaps you can do this using features.  Probably.  Anyhow, this will work for the GUI.
    It can output a text file, a database, and run an R script to graph the results.
    It also has some GUI functions, which is why it's here.
    """
    
    def __init__(self, fasta_path):
        self.total_sequences = 0
        self.sequences = []
        self.frequency_map = dict()
        self.load_sequences(self, fasta_path)
        self.region = ""
        # Native PDB?
        
    def load_sequences(self, fasta_path):
        """
        Opens Fasta file and appends a list of sequences.
        """
        
        FILE = open(fasta_path, 'r')
        for line in FILE:
            if re.search(">", line):
                continue
            line = line.strip()
            self.sequences.append(line)
            
    def calculate_frequencies(self):
        pass
    
    def output_prettytable(self):
        pass
    
    def output_database(self):
        pass
    
    def run_R_script(self):
        pass
    
    def run_GUI_dialog(self):
        """
        Runs a series of GUI dialogs to get information from the user
        """
        pass