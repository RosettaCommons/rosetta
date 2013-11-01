#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/analysis.py
## @brief  Functions for analysis in the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.protocols.analysis import *
from rosetta.protocols.vip import *

#Tkinter Imports
import tkMessageBox
import tkSimpleDialog
import tkFileDialog

#Python Imports
import time
import math
import os.path

#Toolkit Imports

import loops as loop_tools
from window_main import global_variables

def print_all_phi_psi(pose):
    tot=pose.total_residue()
    print "Residue - Phi - Psi"
    for i in range(1, tot+1):
        print repr(pose.pdb_info().number(i))+"    "+repr(pose.phi(i))+"  "+repr(pose.psi(i))

def print_all_residue_energies(pose):
    for i in range(0, pose.total_residue()):
        pose.energies().show(i)

def rmsd(native, p, loops_as_strings, ca_only = False, all_atom=False):
    """
    Prints + Returns RMSD for Full Protein, as well as any loops in loops_as_strings.
    """
    rms = ""
    if ca_only:
        rms = CA_rmsd(native, p)
        print "C-Alpha:"
        print "%.3f RMSD"%rms
    elif all_atom:
        rms = all_atom_rmsd(native, p)
        print "All Atom:"
        print "%.3f RMSD"%rms
    else:
        rms = bb_rmsd(native, p)
        print "Backbone:"
        print "%.3f RMSD"%rms
        
    #if start !=0 and end!=0:
    #    start = p.pdb_info().pdb2pose(chain, int(start)); end = p.pdb_info().pdb2pose(chain, int(end))
    #    cut = start+((end-start)/2); loo=Loop(start,end, cut)
    #    loops = Loops()
    #    loops.add_loop(loo)
        
    #    lrms = loop_rmsd(p, native, loops, False)
    #    print "Loop RMSD:" + str(lrms)
    
    if all_atom:bb_only = False
    else: bb_only=True
    loop_rmsd_map = dict(); #[loop_string]:[lrms] and ['total']:[total average loop rms]
    if loops_as_strings:
        all_rosetta_loops = Loops()
        for loop_string in loops_as_strings:
            rosetta_loop = loop_tools.return_rosetta_Loop(p, loop_string)
            single_loops = Loops()
            
            single_loops.add_loop(rosetta_loop)
            all_rosetta_loops.add_loop(rosetta_loop)
            
            lrms = loop_rmsd(p, native, single_loops, ca_only, bb_only)
            print "Region "+loop_string
            print "%.3f RMSD"%lrms
            loop_rmsd_map[loop_string]=lrms
        lrms = loop_rmsd(p, native, all_rosetta_loops, ca_only, bb_only)
        loop_rmsd_map["total"]=lrms
        print "Regions Average: "
        print "%.3f RMSD"%lrms
        
    return rms, loop_rmsd_map

def readFASC(fileName):
    File = open(fileName, 'r')
    fascData = dict()
    for line in File:
        lineSplit = line.split()
        if lineSplit[0]!="pdb":
            x = 3
            for i in range(3, (len(lineSplit)/2)+2):
                if not fascData.has_key(lineSplit[1]):
                    fascData[lineSplit[1]]=dict()
                fascData[lineSplit[1]][lineSplit[x-1]]=lineSplit[x]
                x = x+2
    return fascData
        
#### Analyze Movers ####

def analyze_packing(p):
    pack_mover = PackStatMover()
    print "\nSee the paper on RosettaHoles to find out more about this statistic (Protein Sci. 2009 Jan;18(1):229-39.)"
    pack_mover.apply(p)

def analyze_interface(p, scorefxn, toolkit = None):
    print "Analyzing Interface.  "
    print "\nNo references are directly associated with this protocol. It was used with Documentation for AnchoredDesign application (see that app's documentation) and CAPRI21 interface descrimination. (Steven Lewis)"
    print "The Mover will seperate chains defined in the interface to calculate energy differences.  Repacking is recommended."
    chains = ""
    if (p.conformation().num_chains()==2):
        pass
    else:
        chains = tkSimpleDialog.askstring(title = "-fixedchains", prompt = "Please input chains to keep fixed  - seperated by a space")
        chains.upper()
        if (chains):
            chains = chains.split()
        else:return
        
    pack_together = tkMessageBox.askyesno(message="rePack before separation")
    pack_separated = tkMessageBox.askyesno(message="rePack after separation (Recommended)")

    if (p.conformation().num_chains()==2):
        interface_mover = InterfaceAnalyzerMover(1, True, scorefxn, False, pack_together, pack_separated, False)
        interface_mover.apply(p)
    
    else:
        chain_ids = []
        #Get ChainIDs
        for chain in chains:
            for i in range(1, p.total_residue()+1):
                if (p.pdb_info().chain( i ) == chain):
                    chain_ids.append( p.chain(i))
                    break
        
        #Pass in the set.
        print "Fixedchains: "+repr(chains)
        interface_mover = InterfaceAnalyzerMover(rosetta.Set(chain_ids), True, scorefxn, False, pack_together, pack_separated, False )
        interface_mover.use_tracer(True)
        interface_mover.apply(p)
        print "Interface Analyzer complete."
        
def analyze_loops(p, loops_as_strings):
    """
    Uses AnalyzeLoopMover to print loop information.
    To accurately use this, add cutpoint variants using the FullControl Window.
    """
    
    if not loops_as_strings: return
    loops_object = loop_tools.InitializeLoops(p, loops_as_strings)
    loop_mover = LoopAnalyzerMover(loops_object, True)
    loop_mover.apply(p)

def analyze_vip(p, scorefxn):
    """
    Uses VIP mover to get Mutational information.
    Should be threaded, which will be added soon.
    """
    vip_mover = VIP_Mover()
    vip_mover.set_initial_pose(p)
    old_energy = scorefxn(p)
    print "\nThis is going to take some time...."
    print "This code uses the RosettaHoles approach to identify problematic buried cavities, and suggests a set of mutations that are predicted to improve stability as determined by improvements to the RosettaHoles and Rosetta fullatom energy scores."
    print "NOTE: For full options, please see the Rosetta application."
    print "Please see Borgo, B., Havranek, J.J. (2012), 'Automated selection of stabilizing mutations in designed and natural proteins', Proc. Natl. Acad. Sci. USA, v.109(5) pp.1494-99."
    time.sleep(6)
    if (tkMessageBox.askyesno(message="Continue?")):
        pass
    else:
        return
    cycles = tkSimpleDialog.askinteger(title = "Cycles", prompt="Please enter max cycles", initialvalue=rosetta.basic.options.get_integer_option('cp:ncycles'))
    
    # Rewritten in python From VIP.cc so the same behavior is met.  Just an interface through PyRosetta to the application code.
    not_finished=True
    improved = False
    i=1
    while (not_finished):

        vip_mover.apply()
        out = vip_mover.get_final_pose()
        print "Comparing new energy " + repr(scorefxn(out)) + " with old energy " + repr(old_energy)
        if (old_energy>scorefxn(out)):
            improved = True
        else:
            improved = False
            
        if(improved):
            for j in range(1, p.total_residue()+1):
                if( out.residue(j).name() != p.residue(j).name() ):
                    position = out.pdb_info().number(j)
                    pdb_chain = out.pdb_info().chain(j)
                    print "Accepting mutation at position "+repr(pdb_position)+" chain "+pdb_chain +" from " +p.residue(j).name() +" to " +out.residue(j).name()
            old_energy = scorefxn(out)
        
#### Rotamers ####
    """
    Used in full control window to get energy and probability of the rotamer.
    """
    
def return_energy(p, res, chain=0, type = fa_dun):
    """
    Returns the energy of the rotamer/scoretype by talaris2013.
    """
    
    if chain !=0:
        res = p.pdb_info().pdb2pose(chain, int(res))
    score = create_score_function("talaris2013")
    emap = core.scoring.EMapVector()
    score.eval_ci_1b(p.residue(res), p, emap)
    e = emap[type]
    return e
    
def return_probability(p, res, chain=0):
    """
    Returns probability of rotamer based on -ln(p)=E -> p = e^-E
    """
    
    e = return_energy(p, res, chain)
    p = math.exp(-e)
    return p
    

                
