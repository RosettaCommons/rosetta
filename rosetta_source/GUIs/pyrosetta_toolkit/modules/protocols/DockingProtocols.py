#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/docking.py
## @brief  Docking protocols frontend
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

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
        """
        result = tkMessageBox.askokcancel(title="Continue?", message="Docking requires 10^5-10^6 decoys and the appropriate scorefunction (interchain_cen).  Continue?")
        if not result: return
        to_dock = tkSimpleDialog.askstring(title = "Input", prompt="Please supply two partners (ex - A_B or A_BC ).  The first will be held rigid.")
        if not to_dock: return
        to_dock = to_dock.upper()
        setup_foldtree(self.pose, to_dock, Vector1([1]))
        switch = SwitchResidueTypeSetMover("centroid")
        recover = ReturnSidechainMover(self.pose)
        switch.apply(self.pose)
        low_res_mover = DockingLowRes(self.score_class.score, 1)
        self.run_protocol(low_res_mover)
    
    def high_res_dock(self):
        """
        Uses DockMCMProtocol for high-res docking.  Checks scorefunction for docking patch.
        """
        
        to_dock = tkSimpleDialog.askstring(title = "Input", prompt="Please supply two partners (ex - A_B or A_BC ).  The first will be held rigid.")
        to_dock = to_dock.upper()
        if not to_dock: return
        if self.score_class.ScorePatch.get()!="docking":
            result = tkMessageBox.askokcancel(title="Continue?", message="Docking patch not set in scorefunction (use standard as main weights). Continue?")
            if not result:return
        dock_mover = DockMCMProtocol()
        dock_mover.set_scorefxn(self.score_class.score)
        dock_mover.set_partners(to_dock)
        self.run_protocol(dock_mover)

        