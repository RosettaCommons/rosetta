#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/IO/SessionIO.py
## @brief  Class responsible for managing saving and loading a GUI session.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import shutil

#Tkinter Imports
import Tkinter
import tkSimpleDialog
import tkFileDialog

#Rosetta Imports
from rosetta import *

#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.modules.tools import input as input_tools
from app.pyrosetta_toolkit.window_main import global_variables

class SessionIO:
    """
    This class is responsible for saving and loading a GUI session.
    If you have variables you need to save and load, enter them here.
    This class is not very pretty, but it's not too bad.
    All input files are copied over to the directory so that you can save the session and use it on a different cpu.  Could zip it at the end.
    
    """
    def __init__(self, toolkit):
        self.toolkit = toolkit
        self.indir=""
        self.outdir=""
        
        
        #### GUI Tkinter Variables required to save and load should go here. ####
        self.stringvar_settings = {
        #Output Options
            "outname":self.toolkit.output_class.outname,
            "outdir":self.toolkit.output_class.outdir,
        #Input Options
            "pdb_path":self.toolkit.input_class.pdb_path,
            "region_start":self.toolkit.input_class.region_start,
            "region_end":self.toolkit.input_class.region_end,
            "region_chain":self.toolkit.input_class.region_chain,
            "region_sequence":self.toolkit.input_class.region_sequence,
        #PyMOL Options
            "scorefunction":self.toolkit.score_class.ScoreType,
            "patch":self.toolkit.score_class.ScorePatch,
        }
        
        self.intvar_settings = {
        #Output Options
            "processors":self.toolkit.output_class.processors,
            "decoys":self.toolkit.output_class.decoys,
            "kt":self.toolkit.output_class.decoys,
            "overwrite":self.toolkit.output_class.overwrite,
            "use_boltzmann":self.toolkit.output_class.use_boltzmann,
            "recover_low":self.toolkit.output_class.recover_low,
            "rounds":self.toolkit.output_class.rounds,
            "autosave":self.toolkit.output_class.auto_write,
            "terminal_output":self.toolkit.output_class.terminal_output,
        #PyMOL Options
            "auto_send":self.toolkit.pymol_class.auto_send,
            "keep_history":self.toolkit.pymol_class.keep_history,
            "send_energies":self.toolkit.pymol_class.send_energies,
            "auto_send_region_colors":self.toolkit.pymol_class.auto_send_region_colors,
            "auto_send_residue_colors":self.toolkit.pymol_class.auto_send_residue_colors,
            "send_label":self.toolkit.pymol_class.send_label,
        }
        #########################################################################
        
        
    def save_session(self):
        """
        Saves the current session to a directory.
        That directory has individual settings, and any specific files that are tied to the session (loops, resfile, PDB).
        """
        outdir = self.toolkit.toolkit_home+"/SESSIONS"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        name = tkSimpleDialog.askstring(title="Session Name", prompt="Please enter a name for this session.")
        if not name:return
        outdir = outdir+"/"+name
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        output_start = self.toolkit.output_class.terminal_output.get();
        self.toolkit.output_class.terminal_output.set(1); #Redirect to stdout.
        
        #Saving Main settings
        OUTFILE = open(outdir+"/settings.txt", 'w')
        for option in sorted(self.stringvar_settings):
            if self.stringvar_settings[option].get():
                OUTFILE.write(option+" "+self.stringvar_settings[option].get()+"\n")
            else:
                OUTFILE.write(option+" "+"NA"+"\n")
        for option in sorted(self.intvar_settings):
            OUTFILE.write(option+" "+str(self.intvar_settings[option].get())+"\n")
        OUTFILE.close()


        
        #PDB Save
        if self.toolkit.pose.total_residue()!=0:
            output_tools.dumpPDB(self.toolkit.pose, self.toolkit.native_pose, outdir+"/pose", self.toolkit.score_class.score, True)
            output_tools.dumpPDB(self.toolkit.native_pose, self.toolkit.native_pose, outdir+"/native_pose", self.toolkit.score_class.score, True)
        
        #Loop Save
        output_tools.save_loop_file(self.toolkit.pose, self.toolkit.input_class.regions, False, outdir+"/temp.loop")
        
        #Params Save
        self.save_params()
        
        #PDBList
        if self.toolkit.input_class.PDBLIST.get():
            shutil.copy2(self.toolkit.input_class.PDBLIST.get(), outdir+"/PDBLIST.txt")
            
        #ScoredPDBList
        if self.toolkit.menu_class.score_analyzer.filepath:
            shutil.copy2(self.toolkit.menu_class.score_analyzer.filepath, outdir+"ScoredPDBLIST.txt")
        
        #Constraints
        self.save_constraints()
        
        #Scoring Options
        self.toolkit.score_class.saveAS(outdir+"/current_scorefunction.wts")
        
        self.toolkit.output_class.terminal_output.set(output_start); #Redirect to stdout. 
    
        print "Session saved to: "+outdir
        
    def load_session(self):
        """
        Loads the current session from a directory.
        """
        
        self.indir = tkFileDialog.askdirectory(title="Open session directory", initialdir=self.toolkit.toolkit_home+"/SESSIONS")
        if not self.indir:return
        
        if self.indir == self.toolkit.toolkit_home+"/SESSIONS":
            print "Directory not chosen.  Please try again."
            return
        
        #Order is very important here, or everything would be in a map
        self.load_params()
                
        #PDBList
        if os.path.exists(self.indir+"/PDBLIST.txt"):
            self.toolkit.input_class.PDBLIST.set(self.indir+"/PDBLIST.txt")
        
        #ScoredPDBList
        if os.path.exists(self.indir+"/ScoredPDBLIST.txt"):
            self.toolkit.menu_class.score_analyzer.set_filepath(self.indir+"/ScoredPDBLIST.txt")
        
        #Load Poses
        if os.path.exists(self.indir+"/pose_1.pdb"):
            self.toolkit.input_class.load_pose(self.indir+"/pose_1.pdb")
            pose_from_file(self.toolkit.native_pose, self.indir+"/native_pose_1.pdb")
        
        #Load Score +Settings
        OUTFILE = open(self.indir+"/settings.txt", 'r')
        for line in OUTFILE:
            line = line.strip()
            lineSP = line.split()
            if lineSP[0]=="dir":
                global_variables.current_directory = lineSP[1]
                continue
            if self.stringvar_settings.has_key(lineSP[0]):
                if lineSP[1]=="NA":
                    self.stringvar_settings[lineSP[0]].set("")
                else:
                    self.stringvar_settings[lineSP[0]].set(lineSP[1])
            elif self.intvar_settings.has_key(lineSP[0]):
                self.intvar_settings[lineSP[0]].set(int(lineSP[1]))
            else:
                continue
            
        self.toolkit.score_class.set_scorefunction(self.stringvar_settings["scorefunction"].get(), self.stringvar_settings["patch"].get())
        
        #Constraints (Load score first)
        self.load_constraints()
        
        #Load Loops
        if os.path.exists(self.indir+"/loops.loop"):
            self.toolkit.input_class.load_loop(self.indir+"/loops.loop")
        
        
    #### Functions ####   
    def save_params(self):
        if not self.toolkit.input_class.param_paths: return
        
        param_path = outdir+"/PARAMS"
        if not os.path.exists(param_path):
            os.mkdir(param_path)
        
        for path in self.toolkit.input_class.param_paths:
            shutil.copy2(path, param_path)
    
    def save_constraints(self):
        if self.toolkit.input_class.constraint_file_paths:
            if not os.path.exists(outdir+"/CONSTRAINTS"):
                os.mkdir(outdir+"/CONSTRAINTS")
            for path in self.toolkit.input_class.constraint_file_paths:
                shutil.copy2(path, outdir+"/CONSTRAINTS")
    
    def load_params(self):
        if not os.path.exists(self.indir+"/PARAMS"):return
        FILES = os.listdir(self.indir+"/PARAMS")
        if len(FILES)<1:return
        
        FILEOUT = open('w', self.indir+"/PARAMS/paramPATHS.txt")
        for name in FILES:
            if name == "paramPATHS.txt":continue
            
            path = self.indir+"/PARAMS/"+name
            FILEOUT.write(path+"\n")
        FILEOUT.close()
        self.toolkit.input_class.load_param_list(self.indir+"/PARAMS/paramPATHS.txt")
    
    def load_constraints(self):
        if not os.path.exists(self.indir+"/CONSTRAINTS"):return
        FILES = os.listdir(self.indir+"/CONSTRAINTS")
        if len(FILES)<1:return
        for name in FILES:
            path = self.indir+"/CONSTRAINTS/"+name
            input_tools.add_constraints_to_pose_and_scorefunction(self.toolkit.pose, self.toolkit.score_class.score, 1.0, False, path)
    
        
