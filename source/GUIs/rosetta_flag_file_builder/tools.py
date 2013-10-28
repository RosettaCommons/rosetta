#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/tools.py
## @brief  Functions used by the rosetta gui window
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from Tkinter import *
import tkFileDialog
import os

def setdefaultdir(defaultdir, pwddir):
    """
    Sets the default directory as pwd, saves it.
    """
    directory = tkFileDialog.askdirectory(initialdir = defaultdir)
    SETTINGSFILE = open(pwddir+"/ROSETTATOOLS_SETTINGS.txt", 'w')
    SETTINGSFILE.write("DEFAULTDIR=="+directory)
    SETTINGSFILE.close()

def loaddefaultdir(dir):
    """
    Parses the Setting file for defaultdir
    """
    if os.path.exists(dir+"/ROSETTATOOLS_SETTINGS.txt"):
        SETTINGSFILE = open(dir+"/ROSETTATOOLS_SETTINGS.txt", 'r')
        for line in SETTINGSFILE:
            lineSP = line.split("==")
            if lineSP[0]=="DEFAULTDIR":
                dir = lineSP[1]
                return dir
    else:
        return dir
    
"""    
def loadConfiguration(self):

    Load a rosetta cmd config file.
    
    FILE = tkFileDialog.askopenfile(initialdir = self.defaultdir)
    config = FILE.read()
    #Parse config, take out database and app type. Set curselection to app type.
    configSP = config.split()
    apppath = configSP[0]; configSP.pop(0)
    app = os.path.split(apppath)[1].split('.')[0]
    print app
    for stuff in configSP:
        if re.search("database", stuff):
            ind = configSP.index(stuff)
            configSP.pop(ind)
            configSP.pop(ind)
            break
    config = ' '.join(configSP)
"""
