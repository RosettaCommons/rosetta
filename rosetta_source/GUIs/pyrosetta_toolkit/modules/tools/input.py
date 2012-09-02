#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/input.py
## @brief  general input functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox
import interfaces
import os
import re


#rosetta.init()

pwd = os.getcwd()

"""
PyRosetta Loading controls with TK
"""


def tk_get_directory(open):
    """
    Loads a Folder through the tk File Dialog
    """
    dir = os.getcwd()
    if open==0:
        dirname1 = tkFileDialog.askdirectory(initialdir=dir,title='Pick a directory')
    else:
        dirname1 = tkFileDialog.askdirectory(initialdir=open, title='Pick a directory')
    if not dirname1:
        return
    x=dirname1.split()       
    print dirname1
    return dirname1

def tk_get_file():
    """
    Loads a File through the tk File Dialog
    """
    dir = os.getcwd()
    filename1 = tkFileDialog.askopenfilename(initialdir=dir, title='Pick a file')
    if not filenmame1:
        return
    filename1 = self.fixFilename(filename1)
    print filename1
    return filename1

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
    
