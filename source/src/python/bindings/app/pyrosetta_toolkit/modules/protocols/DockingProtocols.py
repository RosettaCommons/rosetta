#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/docking.py
## @brief  Docking protocols frontend
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.toolbox.task_operations import *

#Tkinter Imports
from Tkinter import *
import tkMessageBox
import tkSimpleDialog

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass

class DockingProtocols(ProtocolBaseClass):
    """
    Basic frontend for interacting with rosetta docking protocols
    Originally, was going to exlclude Low_Res because of the high number of decoys; but with multiprocessing incorporated and better computers every year, might as well.
    """
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)

    def low_res_dock(self):
        """
        Low res docking protocol
        This will have its own UI shortly.
        """
        print "This protocol is not set to run global docking.  Please see the rosettacommons documentation for global docking runs."

        if self.score_class.ScoreType.get()!= "interchain_cen":
            result = tkMessageBox.askyesnocancel(title="Set Docking scorefunction?", message="Docking requires the appropriate scorefunction (interchain_cen).  Temporarily set docking scorefunction?")
            if result:
                start_scorefxn = ScoreFunction()
                start_scorefxn.assign(self.score_class.score)
                low_res_scorefxn = create_score_function("interchain_cen")
            elif result==None:
                return
            else:
                low_res_scorefxn = start_scorefxn


        to_dock = tkSimpleDialog.askstring(title = "Input", prompt="Please supply two partners (ex - A_B or A_BC ).  The first will be held rigid.")
        if not to_dock: return
        to_dock = to_dock.upper()

        repack_interface = tkMessageBox.askyesno(message="Repack interface post-dock (talaris2013)?")
        low_res_mover = DockingLowRes(low_res_scorefxn, 1)
        low_res_wrapper = LowResWrapper(to_dock, repack_interface, low_res_mover, low_res_scorefxn)

        self.run_protocol(low_res_wrapper)

        if result:
            self.score_class.score.assign(start_scorefxn)


    def high_res_dock(self):
        """
        Uses DockMCMProtocol for high-res docking.  Checks scorefunction for docking patch.
        """
        if self.score_class.ScoreType.get()!="docking":
            result = tkMessageBox.askyesnocancel(title="Set docking scorefunction?", message="Standard docking scorefunction not set. Temporarily set scorefunction?")
            if result:
                start_scorefxn = ScoreFunction()
                start_scorefxn.assign(self.score_class.score)
                self.score_class.score.assign(create_score_function_ws_patch("docking", "docking_min"))
            elif result == None:return

        to_dock = tkSimpleDialog.askstring(title = "Input", prompt="Please supply two partners (ex - A_B or A_BC ).  The first will be held rigid.")
        to_dock = to_dock.upper()
        if not to_dock: return

        dock_mover = DockMCMProtocol(Vector1([1]), self.score_class.score, create_score_function("talaris2013"))

        original_ft = self.pose.fold_tree()
        setup_foldtree(self.pose, to_dock, Vector1([1]))
        self.run_protocol(dock_mover)
        if result:
            self.score_class.score.assign(start_scorefxn)
        self.pose.fold_tree(original_ft)

class LowResWrapper:
    """
    Wrapper to return sc + repack interface rotamers after dock.
    """

    def __init__(self, to_dock_string, repack_interface, low_res_mover, low_res_scorefxn):

        self.to_dock_string = to_dock_string
        self.repack_interface = repack_interface
        self.low_res_scorefxn = low_res_scorefxn

        self.low_res_mover = low_res_mover
        self.switch = SwitchResidueTypeSetMover("centroid")

        self.tf = TaskFactory()
        self.initialize_taskfactory()

    def initialize_taskfactory(self):
        self.tf.push_back(InitializeFromCommandline())
        self.tf.push_back(RestrictToRepacking())
        self.tf.push_back(RestrictToInterface(1, 8.0))
        print "Task factory initialized"

    def apply(self, pose):
        original_ft = pose.fold_tree()
        setup_foldtree(pose, self.to_dock_string, Vector1([1]))


        recover = ReturnSidechainMover(pose)
        self.switch.apply(pose)
        original_score = self.low_res_scorefxn(pose)
        self.low_res_mover.apply(pose)
        print "LowResScore Start: "+ repr(original_score)
        print "LowResScore End: "+ repr(self.low_res_scorefxn(pose))
        recover.apply(pose)

        if self.repack_interface:
            pack_mover=PackRotamersMover(create_score_function("talaris2013"))
            pack_mover.task_factory(self.tf)
            pack_mover.apply(pose)

        pose.fold_tree(original_ft)



