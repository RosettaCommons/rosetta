#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/input.py
## @brief  general input functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os
import re

#Tkinter Imports
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox

#Toolkit Imports
import interfaces
from app.pyrosetta_toolkit.window_main import global_variables

def tk_get_directory():
    """
    Loads a Folder through the tk File Dialog.  Uses and Sets current directory.
    """

    dir_name = tkFileDialog.askdirectory(initialdir=global_variables.current_directory, title='Pick a directory')

    if not dir_name:
        return None
    print dir_name
    global_variables.current_directory = dir_name
    return dir_name

def tk_get_file():
    """
    Loads a File through the tk File Dialog. Uses and Sets current directory.
    """
    filename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title='Pick a file')
    if not filename:
        return
    filename = self.fixFilename(filename)
    global_variables.current_directory = os.path.dirname(filename)
    print filename
    return filename

def fixFilename(file):
    """
    Rosetta cannot use spaces in the filename...
    """

    x=file.split()
    if len(x) > 1:
        newFile = x[0]
        for i in range (1, len(x)):
            newFile = newFile + "\ " + x[i]
        return newFile
    else:
        return file


def get_pdb_list_from_dir(dir):
    """
    Returns a list of PDB files within a directory.
    check_button_ck if it returns the full path..
    """

    directorylist =os.listdir(dir)
    directorylist = list(directorylist) #May already be a list, but I forget.
    for Files in directorylist:
        if re.search(".pdb", Files):
            pass
        else:
            directorylist.remove(Files)
    return directorylist

def add_constraints_to_pose_and_scorefunction(pose, score, default_weight = 1.0, constraint_types = False, constraint_file=False):
    """
    Adds constraint from file to pose and score.  Sets all constraint_types to 1.0.
    Can pass an array of constraint_types.
    """
    if pose.total_residue()==0:
        print "Please load a pose."
        return ""
    if not constraint_types:
        constraint_types = [atom_pair_constraint, angle_constraint, dihedral_constraint, coordinate_constraint, constant_constraint]

    if not constraint_file:
        constraint_file = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title = "Open Constraint File")
    if not constraint_file:return
    global_variables.current_directory = os.path.dirname(constraint_file)
    print "Setting constraints to pose and scorefunction at default weight of 1.0 if not already set.  "
    setup = ConstraintSetMover()
    setup.constraint_file(constraint_file)
    setup.apply(pose)

    for constraint in constraint_types:
        if score.get_weight(constraint)==0:
            score.set_weight(constraint, default_weight)
    return constraint_file

def get_residuetypeset_from_path_array(param_path_array, loaded_path_array):
    """
    Returns ResidueTypeSet from an array of paths.
    """

    #So that there are no duplicates.
    params = dict()
    for p in param_path_array:
        #Only add to unique dictionary if path is not in loaded_path_array
        try:
            ind = loaded_path_array.index(p)
        except ValueError:
            params[p]=0
            loaded_path_array.append(p)
    if not params: return
    params_paths = utility.vector1_string()
    params_paths.extend(list(params.keys()))
    residuetypeset = generate_nonstandard_residue_set(params_paths)
    print "Nonstandard residue type set loaded."
    return residuetypeset, loaded_path_array
