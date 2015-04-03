#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/AnalysisProtocols
## @brief  Analysis functions that output decoys - Basically VIP Mover for now.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.vip import *
from rosetta.basic.options import get_integer_option

#Python Imports
import time

#Tkinter Imports
from Tkinter import *
import tkMessageBox
import tkSimpleDialog

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass


class AnalysisProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)


    def analyze_vip(self):
        """
        Uses VIP mover to get Mutational information. Takes a while. Other options can be set through the Options System.
        """
        vip_mover = VIP_Mover()
        
        
        
        print "This code uses the RosettaHoles approach to identify problematic buried cavities, and suggests a set of mutations that are predicted to improve stability as determined by improvements to the RosettaHoles and Rosetta fullatom energy scores."
        print "NOTE: For full options, please us the options system."
        print "Please see Borgo, B., Havranek, J.J. (2012), 'Automated selection of stabilizing mutations in designed and natural proteins', Proc. Natl. Acad. Sci. USA, v.109(5) pp.1494-99."
        message = "This is going to take some time....and will modify the current pose if not outputting >1 decoys"
        tkMessageBox.showinfo(title = 'note', message = message)
        if not (tkMessageBox.askyesno(message="Continue?")): return
        cycles = tkSimpleDialog.askinteger(title = "Cycles", prompt="Please enter max cycles (0 indicates that the protocol will run until favorable mutations are found)", initialvalue=get_integer_option('cp:ncycles'))
        
        VIP = VIP_Wrapper(vip_mover, self.score_class.score, cycles)
        self.run_protocol(VIP)
            
class VIP_Wrapper:
    """
    Basically a wrapper of the main protocol into a mover-like object to pass to ProtocolBaseClass
    """
    def __init__(self, vip_mover, scorefxn, cycles):
        self.vip_mover = vip_mover
        self.scorefxn = scorefxn
        self.cycles = cycles
        
    def apply(self, p):
        # Rewritten in python From VIP.cc so the same behavior is met.  Just an interface through PyRosetta to the application code.
        vip_mover.set_initial_pose(p)
        old_energy = scorefxn(p)
        not_finished=True
        improved = False
        i=1
        out = Pose()
        while (not_finished):
    
            self.vip_mover.apply()
            out = self.vip_mover.get_final_pose()
            print "Comparing new energy " + repr(self.scorefxn(out)) + " with old energy " + repr(old_energy)
            if (old_energy>self.scorefxn(out)):
                improved = True
            else:
                improved = False
                
            if(improved):
                for j in range(1, p.total_residue()+1):
                    if( out.residue(j).name() != p.residue(j).name() ):
                        position = out.pdb_info().number(j)
                        pdb_chain = out.pdb_info().chain(j)
                        print "Accepting mutation at position "+repr(pdb_position)+" chain "+pdb_chain +" from " +p.residue(j).name() +" to " +out.residue(j).name()
                old_energy = self.scorefxn(out)
            
            i+=1
            if self.cycles==0:not_finished=improved
            else:not_finished=(i<=self.cycles)
        
        p.assign(out)
