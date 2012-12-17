#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/for_design.py
## @brief  design protocols.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.basic.options import get_string_option

#from rosetta.protocols.forge.remodel import *

#Python Imports
import time

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
from window_main import global_variables

class DesignProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
        
    def packDesign(self, resFile=False):
        """
        Follows fixbb.cc to allow most user defined options within the GUI.
        Limitations: No symmetry.  Annealers should work fine through the options system.
        """
        
        if not resfile:
            resfile = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Resfile..")
            if not resfile:return
        
        task = TaskFactory.create_packer_task(self.pose)
        parse_resfile(self.pose, task, resFile)
        s_mover = SequenceMover()
        result = tkMessageBox.askyesno(title='-min_pack', message="Pack and minimize sidechains simultaneously?")
        if result:
            design_mover = MinPackMover(self.score_class.score, task)
            s_mover.add_mover(design_mover)
        else:
            design_mover = PackRotamersMover(self.score_class.score, task)
            s_mover.add_mover(design_mover)
        
        result = tkMessageBox.askyesno(title='-minimize_sidecahins', message="Do minimization of side chains after rotamer packing (slower)")
        
        if result:
            mm = MoveMap(); #Empty movemap
            min_mover = MinMover(mm, self.score_class.score, get_string_option('run:min_type'), 0.01, True)
            task_min_mover = TaskAwareMinMover(min_mover, task)
            s_mover.add_mover(task_min_mover)
            
        
        print "Please cite the many references included for Fixed Backbone Design in the Rosetta Manual."
        time.sleep(5)
        self.run_protocol(s_mover)
    
    def remodel(self, blueprint=False):
        """
        Ya.  Not yet unfortunately.
        """
        pass
    
        if not blueprint:
            blueprint = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Blueprint..")
            if not blueprint:return
            
    def pack_residue(self, res, chain):
        """
        Packs an individual residue.  Scores before and after.
        """
        
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        task = TaskFactory.create_packer_task(self.pose)
        task.restrict_to_repacking()
        task.temporarily_fix_everything()
        #task = self._get_set_pack_neighbors(res, task)
        task.temporarily_set_pack_residue(res, True)
        pack_mover = PackRotamersMover(self.score_class.score, task)
        print self.score_class.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_class.score(self.pose)
    
    def design_residue(self, res, chain):
        """
        Designs a residue by task not restricting to repacking.
        """
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        task = TaskFactory.create_packer_task(self.pose)
        task.temporarily_fix_everything()
        #task = self._get_set_pack_neighbors(res, task)
        task.temporarily_set_pack_residue(res, True)
        pack_mover = PackRotamersMover(self.score_class.score, task)
        print self.score_class.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_class.score(self.pose)
        
    def mutateRes(self, res, chain, new_res): 
        new_res = new_res.split(":")
        new_res = new_res[2]
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        #radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=0.0)
        #print radius
        print self.score_class.score(self.pose)
        print "Mutating to " + new_res
        #self.pose = mutate_residue(self.pose, res, new_res, radius, self.score_class.score); TypeError
        self.pose = mutate_residue(self.pose, res, new_res)
        print self.score_class.score(self.pose)
        print "Mutagenesis Complete."
        return self.pose
    
    def _get_set_pack_neighbors(self, res, task):
        """
        Asks for packing radius.  Sets task to prevent repacking for all but within that radius.
        Original Author: Evan H. Baugh, Johns Hopkins University. from mutants.py.
        Somehow not working with task
        """
        radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=0.0)
        
        center = self.pose.residue( res ).nbr_atom_xyz()
        for i in range( 1 , self.pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
            if not i == res or center.distance_squared(self.pose.residue( i ).nbr_atom_xyz() ) > radius**2:
                task.nonconst_residue_task( i ).prevent_repacking()
                
        return task

