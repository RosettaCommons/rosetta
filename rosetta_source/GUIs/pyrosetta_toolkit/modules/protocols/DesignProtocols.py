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
#from rosetta.protocols.forge.remodel import *

#Tkinter Imports
from Tkinter import *
import tkFileDialog

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
from window_main import global_variables

class DesignProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
        
    def packDesign(self, resFile=False):
        
        if not resfile:
            resfile = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Resfile..")
            if not resfile:return
            
        task = TaskFactory.create_packer_task(self.pose)
        task.read_resfile(resFile)
        design_mover = PackRotamersMover(self.score_class.score, task)
        self.run_protocol(design_mover)
    
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
        task.temporarily_fix_everything()
        task.temporarily_set_pack_residue(res, True)
        pack_mover = PackRotamersMover(self.score_class.score, task)
        print self.score_class.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_class.score(self.pose)
        
    def mutateRes(self, res, chain, new_res): 
        new_res = new_res.split(":")
        new_res = new_res[2]
        res = self.pose.pdb_info().pdb2pose(chain, int(res))

        print self.score_class.score(self.pose)
        print "Mutating to " + new_res
        mutate_residue(self.pose, res, new_res)
        print self.score_class.score(self.pose)
        print "Mutagenesis Complete."
        return self.pose