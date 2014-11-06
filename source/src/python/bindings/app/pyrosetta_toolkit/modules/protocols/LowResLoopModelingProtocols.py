#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_modeling_low.py
## @brief  Low res loop modeling protocols
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.loops.loop_mover.perturb import *

#Tkinter Imports
from Tkinter import *
import tkMessageBox
import tkSimpleDialog
import tkFileDialog

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
import app.pyrosetta_toolkit.modules.tools.loops as loop_tools
from app.pyrosetta_toolkit.window_main import global_variables

class LowResLoopModelingProtocols(ProtocolBaseClass):
     def __init__(self, pose, score_class, input_class, output_class):
          ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
              
     def default_CCD(self, fragset=False, fraglength=False):
          """
          LoopMover_Perturb_CCD for Low Resolution Loop Modeling
          """
          print "Please cite: Wang C, Bradley P, Baker D (2007). Protein-protein docking with backbone flexibility. J. Mol. Biol. 373, 503."
          print "As well as: Canutescu A, Dunbrack R., Jr (2003) Cyclic coordinate descent: A robotics algorithm for protein loop closure. Protein Sci. 12, 963."
          print "Additional options can be set using the options system.  Symmetry is not supported at this time."
          
          if self.score_class.ScoreType.get()!="cen_std" != self.score_class.ScorePatch.get() != "score4L":
               result = tkMessageBox.askyesnocancel(title="Set loop modeling scorefunction?", message="Standard loop modeling scorefunction not set (cen_std/score4L patch). Temporarily set scorefunction?")
               if result:
                    start_scorefxn = ScoreFunction()
                    start_scorefxn.assign(self.score_class.score)
                    self.score_class.score.assign(create_score_function_ws_patch("cen_std", "score4L"))
               elif result == None:return
          
          if not fragset:
               fragset = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Fragset File")
               if not fragset:return
               
          if not fraglength:
               fraglength = tkSimpleDialog.askinteger(title="Fragment Length", initialvalue=3)
               
          extend = tkMessageBox.askyesno(title="Extended", message="Start with idealized bond lenghths, angles, and discard original phi/psi?") 
          movemap=MoveMap()
          ft = self.pose.fold_tree(); ft_o = FoldTree()
          ft_o.assign(ft)
          ft.clear()
          #ft.simple_tree(self.pose.total_residue())
          print len(loops_as_string_array)
          ft, movemap, loops_object=loop_tools.InitializeLoops(p, loops_as_string_array, ft, movemap)
          #print "Fold Tree Correct? " + repr(ft.check_fold_tree())
          self.pose.fold_tree(ft)
          Frag = ConstantLengthFragSet(fraglength, fragset)
          
          switch = SwitchResidueTypeSetMover("centroid")
          recover = ReturnSidechainMover(self.pose)
          switch.apply(self.pose)
          
          #Setup extended loop to throw away initial loop structure
          if extend:
               for loop in loops_object:
                    loop.set_extended(True)
          x = LoopMover_Perturb_CCD(loops_object, self.score_class.score, Frag)
          self.run_protocol(x)
          recover.apply(self.pose)
          ft.assign(ft_o)
          if result:
               self.score_class.score.assign(start_scorefxn)
          self.pose.fold_tree(ft)

                 
     def default_KIC(self):
          """
          LoopMover_Perturb_KIC for Low Resolution Loop Modeling
          """
          print "Please cite: Mandell DJ, Coutsias EA, Kortemme T. (2009). Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling. Nature Methods 6(8):551-2."
          print "Additional options can be set using the options system.  Symmetry is not supported at this time."
          
          if self.score_class.ScoreType.get()!="cen_std" != self.score_class.ScorePatch.get() != "score4L":
               result = tkMessageBox.askyesnocancel(title="Set loop modeling scorefunction?", message="Standard loop modeling scorefunction not set (cen_std/score4L patch). Temporarily set scorefunction?")
               if result:
                    start_scorefxn = ScoreFunction()
                    start_scorefxn.assign(self.score_class.score)
                    self.score_class.score.assign(create_score_function_ws_patch("cen_std", "score4L"))
               elif result == None:return
               
          extend = tkMessageBox.askyesno(title="Extended", message="Start with idealized bond lenghths, angles, and discard original phi/psi?")
          movemap=MoveMap()
          ft = self.pose.fold_tree(); ft_o = FoldTree()
          ft_o.assign(ft)
          ft.clear()
          #ft.simple_tree(self.pose.total_residue())
          ft, movemap, loops_object=loop_tools.InitializeLoops(self.pose, self.input_class.loops_as_strings, ft, movemap)
          print "Fold Tree Correct? " + repr(ft.check_fold_tree())
          self.pose.fold_tree(ft)
          
          #Setup extended loop to throw away initial loop structure
          if extend:
               for loop in loops_object:
                    loop.set_extended(True)
                    
          kic_mover= LoopMover_Perturb_KIC(loops_object, self.score_class.score)
          
          switch = SwitchResidueTypeSetMover("centroid")
          recover = ReturnSidechainMover(self.pose)
          switch.apply(self.pose)
          self.run_protocol(kic_mover)
          recover.apply(self.pose)
          ft.assign(ft_o)
          self.pose.fold_tree(ft)
          
          if result:
               self.score_class.score.assign(start_scorefxn)
