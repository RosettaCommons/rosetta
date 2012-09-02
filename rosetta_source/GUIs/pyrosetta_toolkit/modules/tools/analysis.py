#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/analysis.py
## @brief  Functions for analysis in the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from rosetta import *
import math
import loops as loop_tools


def RetPhiPsi(p):
    tot=p.total_residue()
    print "Residue - Phi - Psi"
    for i in range(1, tot+1):
        print repr(p.pdb_info().number(i))+":    "+repr(p.phi(i))+"  "+repr(p.psi(i))

def RetFAEnergyAll(p):
    for i in range(0, p.total_residue()):
        p.energies().show(i)
        
def rmsd(rmsdP, p, loops_as_strings, ca_only = False, all_atom=False):
    """
    Returns RMSD for Full Protein, as well as any loops in loops_as_strings.
    """
    
    if ca_only:
        print "\nCA RMSD %.3f"%CA_rmsd(rmsdP, p)
    elif all_atom:
        print "\nAll Atom RMSD %.3f"%all_atom_rmsd(rmsdP, p)
    else:
        print "\nBB RMSD %.3f"%bb_rmsd(rmsdP, p)
        
    #if start !=0 and end!=0:
    #    start = p.pdb_info().pdb2pose(chain, int(start)); end = p.pdb_info().pdb2pose(chain, int(end))
    #    cut = start+((end-start)/2); loo=Loop(start,end, cut)
    #    loops = Loops()
    #    loops.add_loop(loo)
        
    #    lrms = loop_rmsd(p, rmsdP, loops, False)
    #    print "Loop RMSD:" + str(lrms)
    
    if all_atom:bb_only = False
    else: bb_only=True
    if loops_as_strings:
        all_rosetta_loops = Loops()
        for loop_string in loops_as_strings:
            rosetta_loop = loop_tools.return_rosetta_Loop(p, loop_string)
            single_loops = Loops()
            
            single_loops.add_loop(rosetta_loop)
            all_rosetta_loops.add_loop(rosetta_loop)
            
            lrms = loop_rmsd(p, rmsdP, single_loops, ca_only, bb_only)
            print "\n"+loop_string+" RMSD %.3f"%lrms
        lrms = loop_rmsd(p, rmsdP, all_rosetta_loops, ca_only, bb_only)
        print "\nALL Loop RMSD:%.3f"%lrms
    return

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
    Returns probability of rotamer
    """
    
    e = return_energy(p, res, chain)
    p = math.exp(-e)
    return p
    

                
