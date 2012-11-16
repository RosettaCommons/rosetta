#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_modeling_low.py
## @brief  Low res loop modeling protocols
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
import tools.loops as loop_tools
class lowRes_Loop_Modeling:
     def __init__(self, score_object, pose):
        self.score_object = score_object
        if not self.score_object.centroid_score:
            self.score_object.centroid_score = create_score_function('cen_std')
        
#default_CCD        
     def default_CCD(self, rounds, loops_as_string_array, ob, fragset, fraglength):
         """
         LoopMover_Perturb_CCD for Low Resolution Loop Modeling
         REQUIRES FRAGSET!!
         """
         
         movemap=MoveMap()
         ft = self.pose.fold_tree(); ft_o = FoldTree()
         ft_o.assign(ft)
         ft.clear()
         #ft.simple_tree(self.pose.total_residue())
         print len(loops_as_string_array)
         ft, movemap, loops_object=loop_tools.InitializeLoops(p, loops_as_string_array, ft, movemap)
         print "Fold Tree Correct? " + repr(ft.check_fold_tree())
         self.pose.fold_tree(ft)
         Frag = ConstantLengthFragSet(fraglength, fragset)
         x = LoopMover_Perturb_CCD(loops_object, self.score_object.centroid_score, Frag)
         for i in range(1, rounds+1):
             print self.score_object.centroid_score(self.pose); x.apply(self.pose); print self.score_object.centroid_score(self.pose)
         ft.assign(ft_o)
         self.pose.fold_tree(ft)

         
#default_KIC            
     def default_KIC(self, p, rounds, loops_as_string_array, ob):
         """
         LoopMover_Perturb_KIC for Low Resolution Loop Modeling
         """
         movemap=MoveMap()
         ft = self.pose.fold_tree(); ft_o = FoldTree()
         ft_o.assign(ft)
         ft.clear()
         #ft.simple_tree(self.pose.total_residue())
         print len(loops_as_string_array)
         ft, movemap, loops_object=loop_tools.InitializeLoops(p, loops_as_string_array, ft, movemap)
         print "Fold Tree Correct? " + repr(ft.check_fold_tree())
         self.pose.fold_tree(ft)
         kic_mover= LoopMover_Perturb_KIC(loops_object, self.score_object.centroid_score)
         for i in range(1, rounds+1):
             print self.score_object.centroid_score(self.pose); kic_mover.apply(self.pose); print self.score_object.centroid_score(self.pose)
         
         ft.assign(ft_o)
         self.pose.fold_tree(ft)
