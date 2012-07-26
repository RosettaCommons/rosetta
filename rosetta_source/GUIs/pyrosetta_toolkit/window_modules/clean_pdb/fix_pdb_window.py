#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/clean_pdb/fix_pdb_window.py
## @brief  Simple GUI for cleaning pdbs for rosetta.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkFileDialog
from modules.tools import pdbs as pdb_tools
import os
import re
class FixPDB:
     '''
     Quick and dirty cleaning of PDB.
     '''
   
     def runfixPDBWindow(self, m, r, c):
        '''
        This First asks for a list of PDB paths.  Then it can do two things, depending on the choices:  Remove HETATM, fix Alternate Residues/Occupancies.
        Gives a GUI, because, well, why not? Should be rewritten a bit and expanded...
        m = Main Window
        r = Row
        c = Column
        '''
        self.defaultdir = os.getcwd(); #Needs to get this from toolkit...
        fixWindow = Toplevel(m)
        fixWindow.title("Setup PDBs for Rosetta")
        self.PDBlistPath = StringVar()
        self.hetVar= IntVar(); self.altVar = IntVar(); self.occVar=IntVar()
        self.hetVar.set(True); self.altVar.set(True); self.occVar.set(True)
        #hetVar, altVar, occVar = 1
        #fixDictionary = dict()
        self.pathEntry = Entry(fixWindow, textvariable=self.PDBlistPath)
        self.pathbutton_ = Button(fixWindow, text = "PDBLIST or PDB", command = lambda: self.PDBlistPath.set(self.getfile()))
        self.remHet = Checkbutton(fixWindow, variable=self.hetVar, text = "Remove HETATM")
        self.remAlt = Checkbutton(fixWindow, variable = self.altVar, text = "Fix Alternative Residues")
        self.chaOcc = Checkbutton(fixWindow, variable = self.occVar, text = "Change Occupancies to 1")
        self.gobutton_ = Button(fixWindow, text = "GO", command = lambda: self.runFixPDB())
                            
        self.pathEntry.grid(row = r, column=c)
        self.pathbutton_.grid(row=r, column = c+1)
        self.remHet.grid(row=r, column=c+2, sticky=W)
        self.remHet.config(state=DISABLED)
        self.remAlt.grid(row=r+1, column=c+2, sticky=W)
        self.chaOcc.grid(row=r+2, column=c+2, sticky=W)
        
        self.gobutton_.grid(row=r+3, column=c, columnspan = 3, sticky=W+E)
        #fixDictionary: [het, alt, occ]=num
     def runFixPDB(self):
        '''
        Runs the fixpdbprotocol using options set in fixDictionary.
        '''
        print self.PDBlistPath.get()
        filename = self.PDBlistPath.get()
        outdirectory = tkFileDialog.askdirectory(title = "Output DIR", initialdir = self.defaultdir)
        if re.search(".pdb", filename):
            self.fixonepdb(filename, outdirectory)
        elif re.search(".txt", filename):
            print "Fixing All PDBs in list.."
            LIST = open(filename, 'r')
            for pdbpath in LIST:
                pdbpath = pdbpath.rstrip()
                print pdbpath
                self.fixonepdb(pdbpath, outdirectory)
            LIST.close()
            print "PDBs Fixed..."
        else:
            print "Please choose a .pdb or .txt file..."
            return
    
     def fixonepdb(self, filename, outdir):
        pdbDic = dict()
        pdbDic = pdb_tools.pdbTools().parsePDB(filename, pdbDic)
        if self.hetVar.get()==1:
            print "Removing HETATM Residues"
        if self.altVar.get()==1:
            pdbDic = pdb_tools.pdbTools().RemAlt(pdbDic)
        if self.occVar.get()==1:
            pdbDic = pdb_tools.pdbTools().chaOcc(pdbDic)
        
        pdbname = os.path.basename(filename)
        pdb_tools.pdbTools().savePDB(pdbDic, outdir+"/"+pdbname)
        print "File Saved..."
     def getfile(self):
        '''
        simply gets a filename...wish I could do this within the command syntax of tkinter...
        '''
        filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
        return filepath