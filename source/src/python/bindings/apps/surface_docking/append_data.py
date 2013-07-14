# Appends data to end of PDB output file
# Emily Koo

# Import python modules
import os
import math
import sys

# Import Rosetta modules
from rosetta import *

#========================= Data ========================================
def append_scores(file, pose, scorefxn):
    # Appends energy scores to end of PDB
    prot_start = pose.num_jump() + 1
    prot_end = pose.total_residue()

    pose_file = open(file,'a')

    res = 1 - prot_start

    # Full atom scores
    if pose.is_fullatom():
        list = [fa_atr,fa_rep,fa_sol,fa_intra_rep,fa_pair,fa_dun, hbond_sr_bb,\
                hbond_lr_bb,hbond_bb_sc,hbond_sc,fa_elec,p_aa_pp,ref, pro_close,\
                atom_pair_constraint, dihedral_constraint]
    else:
        # Centroid scores
        list = [vdw, pair, env, cbeta]

    # Write total energy
    pose_file.write("\nTotal score: \n" )
    stream = rosetta.utility.OStringStream()
    scorefxn(pose)
    scorefxn.show(stream, pose)
    pose_file.write(stream.str())

    # Write all energy scores for each residue in protein
    pose_file.write("\nEnergies by residue: ")
    for i in range(prot_start, prot_end+1):
        # Print res as actual protein residue (for adsorbed) without surface res
        res += i
        pose_file.write("\nResidue no.: "+str(res)+"\nResidue type: "+\
        str(pose.residue(i).name()))

        for  energy in list:
            score = '{0:.3f}'.format(pose.energies().residue_total_energies(i)[energy])
            pose_file.write('{0:20}'.format("\n"+str(energy)) + "\t"+str(score).rjust(10))

        res -= i

        pose_file.write("\n")
    pose_file.close()

def append_hbonds(file, pose, scorefxn):
# Appends hbond scores to end of PDB

    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    pose.update_residue_neighbors()
    rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbond_set)
    stream = rosetta.utility.OStringStream()
    hbond_set.show(pose, True, stream)

    pose_file = open(file, 'a')
    pose_file.write("\nbegin protein intra-molecular hydrogen bonds")
    pose_file.write("\ndonor atom acceptor atom energy\n")
    pose_file.write(stream.str())
    pose_file.write("end protein intra-molecular hydrogen bonds\n")
    pose_file.close()

def append_constraints(file, pose, cst, scorefxn):
    if cst:

        # Appends constraint data to end of PDB
        prot_start = pose.num_jump() + 1
        prot_end = pose.total_residue()

        pose_file = open(file,'a')
        scorefxn(pose)
                    
        print "Constraints energy:"
        scorefxn.show(pose)

        pose_file.write("\nConstraint Energies")
        pose_file.write("\nAtom pair energy: " + str(pose.energies().total_energies()[atom_pair_constraint]))
        pose_file.write("\nDihedral energy: " + str(pose.energies().total_energies()[dihedral_constraint]))

        """
        # Booleans for headers
        m = 0
        n = 0
        lowest_dist = []
        prev_res1 = 0
        amb = 0
        
        # {{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
        for cst_file in cst:
            for line in open(cst_file.strip(), 'r'):
                # if line is not empty
                if line.strip() is not "":
                    data = line.split()
                    # Dihedrals ==============================================================================================
                    if data[0] == 'Dihedral':
                        if m == 0:
                            pose_file.write("\nType\t\t\tRes\tActual\tConstr\tError\t\tScore")
                            m = 1
                            
                        res1 = int(data[4])
                        # Normalize the value such that it is the actual protein res number
                        res1 +=  - prot_start + 1

                        # Phi/Psi
                        if data[1] == 'C':
                            type = "Phi"
                            angle = '{0:.3f}'.format(math.radians(pose.phi(int(data[4]))))
                        else:
                            type = "Psi"
                            angle = '{0:.3f}'.format(math.radians(pose.psi(int(data[4]))))

                        score = '{0:.3f}'.format(pose.energies().residue_total_energies(res1)[dihedral_constraint])
                        constr = '{0:.3f}'.format(float(data[10]) + float(data[12]))
                        error = data[12]
                        pose_file.write("\n"+str(data[0])+"\t"+type+"\t"+str(res1)+\
                                        "\t"+str(angle)+ "\t"+str(constr)+"\t"+str(error)+"\t\t"+\
                                        str(score))

                    # Atom pair ================================================================================================
                    elif data[0] == 'AtomPair' or data[1] == 'AtomPair':
                        if n == 0:
                             pose_file.write("\nType\t\t\tAtom1\tRes1\tAtom2\tRes2\tActual\tConstr\tError\tScore")
                             n = 1
                             
                        if data[1] == 'AtomPair':
                            data = data[1:]
                            amb = 1
                            type = "AmbiConstr"

                        else:
                            amb = 0
                            type = "AtomPair"

                        res1 = int(data[2])
                        res2 = int(data[4])
                        # Normalize the value such that it is the actual protein res number
                        res1 += - prot_start + 1
                        if type == "AtomPair" or res2 > prot_start -1:
                            res2 += - prot_start + 1

                        if type == "AmbiConstr" and prev_res1 == 0:
                            prev_res1 = res1
                                
                        # Calculate distance between the two atoms
                        start = pose.residue(res1).xyz(str(data[1]))
                        end = pose.residue(res2).xyz(str(data[3]))
                        dist = '{0:.3f}'.format((end - start).norm)
                        score = '{0:.3f}'.format(pose.energies().residue_total_energies(res1)[atom_pair_constraint])
                        constr = '{0:.3f}'.format(float(data[6]) + float(data[8]))
                        error = '{0:.3f}'.format(float(data[8]))

                        # For ambiguous constraints
                        if amb == 1:
                        # Only output pair with lowest distance
                        # Keep going if res in curr line is same as res in prev line (investigating same ambiguous constr)
                            if res1 is prev_res1:
                                lowest_dist.append([dist, str(data[1]), str(res1), str(data[3]), str(res2), str(constr), str(error), str(score)])
                                
                            else:
     
                                low = min(lowest_dist)
                                pose_file.write("\n"+str(type)+"\t\t"+low[1]+"\t\t"+low[2]\
                                            +"\t\t"+low[3]+"\t\t"+low[4]+"\t\t"+low[0]+"\t"+low[5]\
                                            +"\t"+low[6]+"\t"+low[7])
                                # clear list
                                lowest_dist = []

                            prev_res1 = res1
                            
                        else:
                            pose_file.write("\n"+str(type)+"\t\t"+str(data[1])+"\t\t"+str(res1)\
                                            +"\t\t"+str(data[3])+"\t\t"+str(res2)+"\t\t"+str(dist)+"\t"+str(constr)\
                                            +"\t"+str(error)+"\t"+str(score))

        # }}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

        # Append last amb constraint (no more change in res so it doesn't go into else statement)
        low = min(lowest_dist)
        pose_file.write("\n"+str(type)+"\t\t"+low[1]+"\t\t"+low[2]\
                    +"\t\t"+low[3]+"\t\t"+low[4]+"\t\t"+low[0]+"\t"+low[5]\
                    +"\t"+low[6]+"\t"+low[7])
        
        """
        pose_file.close()
                
def append_rmsd(file, pose, native):
    pose_file = open(file, 'r')
    
    for line in pose_file:
        if line[:22] == " Total weighted score:":
            score = line.split()[3]
            break
    pose_file = open(file, 'a')

    if file.startswith("Ads"):
	rmsd = 0
    else:
	rmsd = CA_rmsd(native, pose)

    pose_file.write("\nRMSD and Total Weighted Score: \n"+str(rmsd)+"\t" + str(score))
    pose_file.close()

