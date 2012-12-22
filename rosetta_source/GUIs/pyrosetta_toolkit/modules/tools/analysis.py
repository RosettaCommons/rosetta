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

#Tkinter Imports
import tkMessageBox
import tkSimpleDialog
import tkFileDialog

#Python Imports
import time
import math
import os.path

#Toolkit Imports
from rosetta.protocols.analysis import *
import loops as loop_tools
from window_main import global_variables

def RetPhiPsi(p):
    tot=p.total_residue()
    print "Residue - Phi - Psi"
    for i in range(1, tot+1):
        print repr(p.pdb_info().number(i))+":    "+repr(p.phi(i))+"  "+repr(p.psi(i))

def RetFAEnergyAll(p):
    for i in range(0, p.total_residue()):
        p.energies().show(i)

def rmsd(native, p, loops_as_strings, ca_only = False, all_atom=False):
    """
    Prints + Returns RMSD for Full Protein, as well as any loops in loops_as_strings.
    """
    rms = ""
    if ca_only:
        rms = CA_rmsd(native, p)
        print "\nCA RMSD %.3f"%rms
    elif all_atom:
        rms = all_atom_rmsd(native, p)
        print "\nAll Atom RMSD %.3f"%rms
    else:
        rms = bb_rmsd(native, p)
        print "\nBB RMSD %.3f"%rms
        
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
            print "\n"+loop_string+" RMSD %.3f"%lrms
            loop_rmsd_map[loop_string]=lrms
        lrms = loop_rmsd(p, native, all_rosetta_loops, ca_only, bb_only)
        loop_rmsd_map["total"]=lrms
        print "\nALL Loop RMSD:%.3f"%lrms
        
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

def analyze_interface(p, scorefxn):
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
        interface_mover = InterfaceAnalyzerMover(1, True, scorefxn, False, pack_together, pack_separated)
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
        interface_mover = InterfaceAnalyzerMover(rosetta.Set(chain_ids), True, scorefxn, False, pack_together, pack_separated)
        interface_mover.apply(p)
        
def analyze_loops(p, loops_as_strings):
    """
    Uses AnalyzeLoopMover to print loop information.
    To accurately use this, add cutpoint variants using the FullControl Window.
    """
    
    if not loops_as_strings: return
    loops_object = loop_tools.InitializeLoops(p, loops_as_strings)
    loop_mover = LoopAnalyzerMover(loops_object, True)
    loop_mover.apply(p)
        
#### Rotamers ####
    """
    Used in full control window to get energy and probability of the rotamer.
    """
    
def return_energy(p, res, chain=0, type = fa_dun):
    """
    Returns the energy of the rotamer/scoretype.
    """
    
    if chain !=0:
        res = p.pdb_info().pdb2pose(chain, int(res))
    score = create_score_function_ws_patch('standard', 'score12')
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
    

                
