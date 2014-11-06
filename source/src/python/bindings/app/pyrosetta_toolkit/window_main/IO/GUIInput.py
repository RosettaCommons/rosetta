
#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/IO/GUIInput.py
## @brief  Class responsible for managing input variables of the GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
## @author Steven Combs (steven.combs1@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import urllib2
import os.path
import re

#Tkinter Imports
from Tkinter import *
from Tkinter import StringVar
import tkFileDialog
import tkSimpleDialog

#Toolkit Imports
from app.pyrosetta_toolkit.window_main import global_variables
from app.pyrosetta_toolkit.modules.RegionalScoring import RegionalScoring
from app.pyrosetta_toolkit.modules.Region import Region
from app.pyrosetta_toolkit.modules.Region import Regions
from app.pyrosetta_toolkit.modules.tools import input as input_tools

from app.pyrosetta_toolkit.window_modules.clean_pdb.FixPDBWindow import FixPDBWindow
from app.pyrosetta_toolkit.window_modules.options_system.OptionSystemManager import OptionSystemManager
#from app.pyrosetta_toolkit import main_window

class GUIInput:
    def __init__(self, toolkit):
        self.toolkit = toolkit; #Basically an AP of the toolkit
        self.pose = self.toolkit.pose
        self.pdb_path = StringVar(); self.pdb_path.set("");
        self.PDBLIST = StringVar(); self.PDBLIST.set("")
        
        self.region_start=StringVar(); #Start of Region
        self.region_end=StringVar(); #End of Region
        self.region_chain=StringVar(); #Chain of Region
        self.region_sequence=StringVar(); #Sequence in Entry
        self.loops_as_strings = []; #Array of Regions: start:end:chain
        self.regions = Regions(); #This will replace loops_as_strings
        self.loops = Loops()
        
        #These are set if a user selects a residue from the sequence
        self.residue_string = StringVar(); #String that is displayed when individual reside is selected
        self.residue_resnum = StringVar();
        self.residue_rosetta_resnum = StringVar();
        self.residue_chain = StringVar()
        
        self.constraint_file_paths = []; #Path to constraint file if loaded through GUI (Not options system).
        self.param_pathlist_file = ""; #Path to a file which lists paths to all params to use.  One on each line
        self.param_paths = []; #Array of parameter paths.
        self.loaded_paths = [];  #Since ResidueTypeSet is a singleton, with horrible exception handling, WE need to keep track of it.
        self.nonstandard_ResidueTypeSet = ""; #This is set through the ncaa window or loading a param path file.  
        
        self.options_manager= OptionSystemManager(); #This is due to Protocols needing Rosetta to be reinitialized without loosing already set options -  to set the seed up before multiprocessing runs!
        
        self.options_manager= OptionSystemManager(); #This is due to Protocols needing Rosetta to be reinitialized without loosing already set options -  to set the seed up before multiprocessing runs!
        
        self.pdb_url = "http://www.rcsb.org/pdb/files"
        #if 0: self.toolkit = main_window()
        
        
        
        
        
        
        
        
######### Functions that cannot be put in input_tools or do not belong in the main frame, as they set a variable within this class. ################

    def choose_load_pose(self, message="Load Pose"):
        """
        Loads a Pose through the tk File Dialog
        """
        infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title=message)
        if not infilename:return
        
        global_variables.current_directory= os.path.dirname(infilename)
        print global_variables.current_directory
        self.load_pose(infilename)
        
    def fetch_pdb(self):
        """
        Fetches the PDB, opens FixPDBWindow to allow the user to clean the PDB before trying to load it.
        """
        #Create default directory
        outpath = self.toolkit.toolkit_home+"/PDBs"
        if not os.path.exists(outpath):os.mkdir(outpath)
        global_variables.current_directory = outpath
        
        #Ask for PDB id.
        pdbID = tkSimpleDialog.askstring(title="Fetch pdb", prompt="Please enter PDB ID.")
        if not pdbID: return
        
        #Open and Write the PDB
        FILE = urllib2.urlopen(self.pdb_url+'/'+pdbID.lower()+'.pdb')
        OUTFILE = open(outpath+'/'+pdbID.upper()+'.pdb', 'w')
        for line in FILE:
            OUTFILE.write(line)
        OUTFILE.close()
        
        fetched_pdb = outpath+'/'+pdbID.upper()+'.pdb'
        print "PDB saved to pyrosetta_toolkit/PDBs"
        cleaner = FixPDBWindow(self, self.toolkit.score_class, self.toolkit.pose, fetched_pdb)
        cleaner.runfixPDBWindow(self.toolkit._tk_, 0, 0)
    
    def select_pose_then_launch_fixpdb(self):
        """
        This way of loading a pose asks the user to open a PDB, then launches the fixPDBWindow as per Labonte's suggestion.
        """
        
        infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory, title="Select PDB file")
        if not infilename:return
        
        global_variables.current_directory= os.path.dirname(infilename)
        print global_variables.current_directory
        cleaner = FixPDBWindow(self, self.toolkit.score_class, self.toolkit.pose, infilename)
        cleaner.cleaned_pdb_path.set(infilename)
        cleaner.runfixPDBWindow(self.toolkit._tk_, 0, 0)
        cleaner.enable_load()
        
    def load_pose(self, path):
        """
        Load a pose into the toolkit.pose variable.  Setup nessessary variables/etc for objects and window objects of the toolkit.  Can have NCAA that have been enabled.
        Please use this when loading the main pose into the toolkit.
        """
        if not os.path.exists(path):
            print "PDB path does not exist.  Cannot load pose."
            return
        
        self.pdb_path.set(path)
        print self.pdb_path.get()
        
        #Turn off PyMOL Observer if on - It will try to update on new pose.
        observer_on = self.toolkit.pymol_class.auto_send.get()
        if observer_on:
            self.toolkit.pymol_class.auto_send.set(False)
            
        #Load Pose
        if self.nonstandard_ResidueTypeSet:
            self.toolkit.pose.assign(pose_from_pdb(self.nonstandard_ResidueTypeSet, path))
        else:
            pose_from_pdb(self.toolkit.pose, self.pdb_path.get())
        self.toolkit.native_pose.assign(self.toolkit.pose); #Set native pose for RMSD.

        print self.toolkit.pose
        
        #Reinitialize PyMOL
        self.toolkit.pymol_class.SendNewPose()
        if observer_on:
            self.toolkit.pymol_class.auto_send.set(True)
            
        self.regional_score_class = RegionalScoring(self.toolkit.pose, self.toolkit.score_class.score);
        
        #Reinitialize Output
        pdbname = os.path.basename(self.pdb_path.get())
        pdbname = pdbname.split(".")[0]
        self.toolkit.output_class.outname.set(pdbname)
        self.toolkit.output_class.outdir.set(os.path.dirname(self.pdb_path.get()))
        
        

        
        #Reinitialize Sequence + Regions
        self.reinit_regions_on_new_pose()
        self.region_sequence.set(self.toolkit.pose.sequence())
        self.residue_rosetta_resnum.set("")
        
        #Reinitialize Design
        self.toolkit.DesignDic = dict()
        
    def return_loaded_pose(self, path):
        """
        Load and return a pose.  Can have NCAA that have been enabled.
        """
        p = Pose()
        if self.nonstandard_ResidueTypeSet:
            p.assign(pose_from_pdb(self.nonstandard_ResidueTypeSet, path))
        else:
            pose_from_pdb(p, self.pdb_path.get())
        return p
    
    def set_PDBLIST(self):
        infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory,title='Open PDBLIST')
        if not infilename:return
        global_variables.current_directory =os.path.dirname(infilename)
        print "PDBLIST set"
        self.PDBLIST.set(infilename)
    
    def load_param_list(self, infilename=False):
        """
        Loads paths from param path files into an array.  Creates a residue type set from the params.
        """
        if not infilename:
            infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory,title='Open param pathList file')
        if not infilename: return
        self.param_pathlist_file = infilename
        global_variables.current_directory =os.path.dirname(infilename)
        FILE = open(infilename, 'r')
        for line in FILE:
            if re.search("#", line): continue
            line = line.strip()
            self.param_paths.append(line)
        
        self.nonstandard_ResidueTypeSet, self.loaded_paths = input_tools.get_residuetypeset_from_path_array(self.param_paths, self.loaded_paths)
        FILE.close()
        
    def load_loop(self, infilename=False):
        """
        Returns a loops_as_strings array to be used by the InputFrame
        """
        
        if not self.toolkit.pose.total_residue():
            print "\nNumbering conversion requires a pose to be loaded...please load a pose.."
            return
        
        if not infilename:
            infilename = tkFileDialog.askopenfilename(initialdir=global_variables.current_directory,title='Open Loop file')
        if not infilename: return
        global_variables.current_directory =os.path.dirname(infilename)
        FILE = open(infilename, 'r')
        loops_as_strings = []
        for line in FILE:
            print line
            line = line.strip()
            lineSP = line.split()
            start = lineSP[1]
            end = lineSP[2]
            chain_start = self.toolkit.pose.pdb_info().chain(int(start))
            chain_end = self.toolkit.pose.pdb_info().chain(int(end))
            if chain_start != chain_end:
                print "Invalid loop for GUI.  Start and end residues are on different chains.\n"
                return
            loop_string = repr(self.toolkit.pose.pdb_info().number(int(start)))+":"+repr(self.toolkit.pose.pdb_info().number(int(end)))+":"+chain_start
            loops_as_strings.append(loop_string)
        FILE.close()
        return loops_as_strings

    def reinit_regions_on_new_pose(self):
        if not self.regions:return #If regions not yet set
        for region in self.regions:
            if not region.region_exists(self.pose):
                loop_string = region.get_region_string()
                self.loops_as_strings.remove(loop_string)
                self.regions.remove_region(loop_string)
            else:
                loop_string = region.get_region_string()
                print loop_string +" found in new Pose"
        
        #Here is where we have to actually interact with the frame.
        
        self.toolkit.input_frame.loops_listbox.delete(0, END)
        
        for loop_string in self.loops_as_strings:
            self.toolkit.input_frame.loops_listbox.insert(END, loop_string)
###Region Setting + Residue Setting ####
    def set_residue_of_interest(self, resnum, chain, rosetta_resnum):
        """
        Sets current individual residue information.
        All Strings.
        """
        self.residue_resnum.set(resnum)
        self.residue_chain.set(chain)
        self.residue_rosetta_resnum.set(rosetta_resnum)
           
