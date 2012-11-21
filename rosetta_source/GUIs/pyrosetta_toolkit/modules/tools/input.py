#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
from window_main import global_variables

def tk_get_directory():
    """
    Loads a Folder through the tk File Dialog
    """

    dir_name = tkFileDialog.askdirectory(initialdir=global_variables.current_directory, title='Pick a directory')
    
    if not dir_name:
        return      
    print dir_name
    global_variables.current_directory = dir_name
    return dir_name

def tk_get_file():
    """
    Loads a File through the tk File Dialog
    """
    filename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title='Pick a file')
    if not filename:
        return
    filename = self.fixFilename(filename)
    global_variables.current_directory = os.path.dirname(filename)
    print filename
    return filename

def load_pdb(file):
    """
    Loads a PDB into a pose
    """
    print file
    p = Pose()
    pose_from_pdb(p, file)
    return p

def load_loop_file(pose, filename):
    """
    Loads the loop.  Doesn't work.
    """
    INFILE = open(filename, 'r')
    loops_as_strings = []
    for line in INFILE:
        lineSP=line.split()
        start_pdb = pose.pdb_info().pose2pdb(int(lineSP[1]))
        end_pdb   = pose.pdb_info().pose2pdb(int(lineSP[2]))
        string = repr(start_pdb[0])+":"+repr(end_pdb[0])+":"+start_pdb[1]
        loops_as_strings.append(string)
    INFILE.close()
    return loops_as_strings

def load_vicinity(p, lisloop, vaccinity):
    """
    Takes LisLoop, Vaccinity Settings, and returns the vacDic.
    See vaccinity window module for more.
    """
       
        #Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
    t = time.time()
    tempdir = pwd + "/temp2/" + "_"+repr(t)
    os.mkdir(tempdir)
    p.dump_pdb(tempdir+"/temp.pdb"); temp = tempdir + "/temp.pdb"
        
        #Gets Residues that are in Vaccinity using python. ([res:chain]=atomic contact #)
    vacDic = tools.interfaces.around().getVaccinity(p, lisloop, temp, vaccinity)
    
    rmtree(tempdir)
    #os.mkdir(pwd+"/temp")
    return vacDic


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
    
