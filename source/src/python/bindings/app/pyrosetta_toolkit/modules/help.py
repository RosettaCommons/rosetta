#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/help.py
## @brief  Functions for the help menu item.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Tkinter Imports
import tkMessageBox
import tkSimpleDialog

#### Collection of help functions. #####

def about():
    """
    This will be a menu pop up
    """
    message = "GUI Originally Created by Jared Adolf-Bryfogle, March 2011-2013.\n"
    message = message + "Roland Dunbrack / FCCC / Drexel College of Medicine.\n"
    message = message + "PyRosetta created by Jeffrey J. Gray, Sergey Lyskov, and the PyRosetta Team at JHU.\n"
    message = message + "All Rosetta code is a community effort by the extraordinary RosettaCommons Team."
    tkMessageBox.showinfo(title="About", message=message)


def citation():
    pass

def region_selection():
    message = "Residue numbers are in PDB numbering.  \n \
*To specify only a chain, add a chain with no residues \n \
*To specify only a tail: only specify end for Nter tail, only specify start for Cter tail \n \
*To specify a single residue, use the same PDB number for start and end."
    tkMessageBox.showinfo(title="Region Selection", message=message)

def mutSA():
    message = "Accessible Surface Area: Square Angstroms, (Chothia 1976)"
    tkMessageBox.showinfo(title="Design Help", message=message)
def mutRM():
    message = "Relative Mutability: Arbitrary units (ALA=100).  Probability Residue will mutate in given time. (Dayhoff 1978)"
    tkMessageBox.showinfo(title="Design Help", message=message)
def mutSP():
    message = "Surface Probability: P(x). Liklihood that 5% or more of the surface area will be exposed to solution. (Chothia 1976)"
    tkMessageBox.showinfo(title="Design Help", message=message)

def show_license():

    message = "(c) Copyright Rosetta Commons Member Institutions. \n \
(c) This file is part of the Rosetta software suite and is made available under license.\n \
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.\n \
(c) For more information, see http://www.rosettacommons.org. Questions about this can be \n \
(c) addressed to University of Washington CoMotion, email: license@uw.edu."
    tkMessageBox.showinfo(title="License", message=message)

def print_glossary():
    """
    Added by Steven Combs.
    """

    message = "*all-atom = in the case of sampling, synonymous with fine movements and often including side chain information; also referred to as high-resolution \n \
*benchmark = another word for a test of a method, scoring function, algorithm, etc. by comparing results from the method to accepted methods/models \n \
*binary file = a file in machine-readable language that can be executed to do something in silico \n \
*BioPython = a set of tools for biological computing written and compatible with Python http://biopython.org/wiki/Biopython \n \
*build = to compile the source code so it can be used as a program \n \
*centroid = in Rosetta centroid mode, side chains are represented as unified spheres centered at the residue?s center of mass\n \
*cluster center = the geometric center of a cluster, or group, of models \n \
*clustering = in this case, grouping models with similar structure together \n \
*comparative model = a protein model where the primary sequence from one protein (target) is placed, or threaded, onto the three dimensional coordinates of a protein of known structure (template)language (binary)"
    tkMessageBox.showinfo(title="Rosetta Glossary from Nature Protocols Paper. Written by Stephanie Hirst DeLuca", message=message)


def scwrl():
    message = "To use SCWRL with the GUI, please copy compiled/installed scwrl directory (Scwrl4) pyrosetta_toolit/scwrl/[my_platform]"
    tkMessageBox.showinfo(title="SCWRL Info", message=message)


