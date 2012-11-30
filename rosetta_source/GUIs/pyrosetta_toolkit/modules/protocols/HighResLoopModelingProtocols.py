#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_modeling_high.py
## @brief  High res loop modeling protocols. 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.loops.loop_mover.refine import *

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
import modules.tools.loops as loop_tools

class HighResLoopModelingProtocols(ProtocolBaseClass):
    """
    Simple class for handling HighRes loop modeling...
    """
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
         
    def default_CCD(self):
        """
        Uses LoopMover_Refine_CCD as High Res Mover
        """
        movemap=MoveMap()
        ft = self.pose.fold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        #ft.simple_tree(self.posetotal_residue())
        print len(self.input_class.loops_as_strings)
        ft, movemap, loops_object=loop_tools.InitializeLoops(self.pose, self.input_class.loops_as_strings, ft, movemap)
        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.pose.fold_tree(ft)        
        loop_refine = LoopMover_Refine_CCD(loops_object)
        loop_refine.set_scorefxn(self.score_class.score)
        self.run_protocol(loop_refine)
        ft.assign(ft_o)
        self.pose.fold_tree(ft)
          
    def default_KIC(self):
        """
        Uses LoopMover_Refine_KIC as High Res Mover
        """

        movemap=MoveMap()
        ft = self.pose.fold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        #ft.simple_tree(self.posetotal_residue())
        print len(self.input_class.loops_as_strings)
        ft, movemap, loops_object=loop_tools.InitializeLoops(self.pose, self.input_class.loops_as_strings, ft, movemap)

        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.pose.fold_tree(ft)
        loop_refine = LoopMover_Refine_KIC(loops_object)
        loop_refine.set_scorefxn(self.score_class.score)
        self.run_protocol(loop_refine)
        ft.assign(ft_o)
        self.pose.fold_tree(ft)
