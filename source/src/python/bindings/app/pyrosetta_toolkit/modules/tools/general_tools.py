#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/general_tools.py
## @brief  general functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os
import math
import re
from sys import platform
from shutil import rmtree
import time

#Tkinter Imports
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox

#Toolkit Imports
from app.pyrosetta_toolkit.modules.Region import Region
from app.pyrosetta_toolkit.modules.Region import Regions

pwd = os.getcwd()


def renameandsave(inFolder, outFolder, outName, contains):
    """
    Renames all files in a particular directory recursively from 1 - N. Useful due to some apps not being JD2 compatible!
    """
    
    filenum = 1
    if contains=="all":
        contains = ""
    #print outFolder
    AllFiles = []
    for root, dirs, files in os.walk(inFolder, topdown=True):
        #print "Root" + root
        for f in files:
            if re.search(contains, f):
                print "File_"+repr(filenum)+"_"+f
                p = os.path.join(root, f)
                AllFiles.append(p)

                
                                  
    for f in AllFiles:
        os.system("cp "+f+" "+outFolder)
        print f
        fSP = os.path.split(f)
        fileSP = fSP[1].split(".")
        newName = outName+"_"+repr(filenum)+"."+fileSP[1]
        print newName
        newPath = os.path.join(outFolder, newName)
        oldPath = os.path.join(outFolder, fSP[1])
        os.system("mv "+oldPath+" "+newPath)
        filenum+=1
    print "Files Copied.."
    return

    
def getDist(p, res1, res2, atom1, atom2):
    """
    Gets distance between atom one and two of two residues.
    """
    xyz1 = p.residue(res1).xyz(atom1)
    xyz2 = p.residue(res2).xyz(atom2)
    return getDistGen(xyz1, xyz2)

def getDistGen( xyz1, xyz2):
    """
    Gets distance bt two coord vectors(list)
    """
    
    #xyz1 is a list with (x, y, z)
    d = math.sqrt(pow(xyz1[0]-xyz2[0], 2)+pow(xyz1[1]-xyz2[1], 2)+pow(xyz1[2]-xyz2[2], 2))
    return d

def getOS():
    """
    Get OS of the particular platform the toolkit is being run on.
    """
    
    plat = sys.platform
    if re.search("darwin", plat):
        return "Mac"
    elif re.search("linux", plat):
        return "Linux"
    elif re.search("win", plat):
        return "Windows"
    else:
        print "Platform Not Found"
        return "error"
    
def loop_string_to_region(loop_string):
    """
    Loop string (start:end:chain) conversion to newer Region class.
    """
         
    start = loop_string.split(":")[0]; end = loop_string.split(":")[1]; chain = loop_string.split(":")[2]
    #Chain
         
    if (start == "" and end==""):
        region = Region(chain.upper(), None, None)
    #Nter
    elif start=="":
        region = Region(chain.upper(), None, int(end))
    #Cter
    elif end=="":
        region = Region(chain.upper(), int(start), None)
    #Loop
    else: 
        region = Region(chain.upper(), int(start), int(end))
    
    return region

def loops_as_strings_to_regions(loops_as_strings):
    """
    Loops as strings representation of regions to newer Regions class.
    """
    reg = Regions()
    for loop_string in loops_as_strings:
        reg.add_region(loop_string_to_region(loop_string))
    
    return reg
    
